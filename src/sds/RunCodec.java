package sds;

import java.util.ArrayList;
import java.util.List;

/**
 * 论文 2.1.3 节 "Burrows-Wheeler transform" + 2.1.6 节 "GBWT index"。
 *
 * 每个节点的 BWT 子串（BWT_v）用游程编码存储：出边数少时用单字节基础编码，
 * 游程长度超过阈值时溢出到变长续码。这里同时给出编码和解码，互为逆操作。
 */
public class RunCodec {

    // ---------------- ByteCode：LEB128 风格变长整数 ----------------
    public static long readByteCode(byte[] a, int[] pos) {
        int i = pos[0];
        long res = a[i] & 0x7FL;
        int shift = 0;
        while ((a[i] & 0x80) != 0) {
            i++;
            shift += 7;
            res += ((long) (a[i] & 0x7F)) << shift;
        }
        i++;
        pos[0] = i;
        return res;
    }

    public static void writeByteCode(SdsWriter w, long value) {
        while (value >= 0x80) {
            w.rawByte((byte) ((value & 0x7F) | 0x80));
            value >>>= 7;
        }
        w.rawByte((byte) (value & 0x7F));
    }

    // ---------------- Run：单字节基础编码 + 续码 ----------------

    public final long sigma;        // 出边数，即当前记录的字母表大小
    public final long runContinues; // sigma<255 时 = 256/sigma，否则恒走双变长码

    public RunCodec(long sigma) {
        this.sigma = sigma;
        long maxCode = 255;
        this.runContinues = (sigma < maxCode) ? (maxCode + 1) / sigma : 0;
    }

    /** 解码一个游程，返回 {outrank, runLength}，并推进 pos。 */
    public long[] read(byte[] a, int[] pos) {
        if (runContinues == 0) {
            long value = readByteCode(a, pos);
            long len = readByteCode(a, pos) + 1;
            return new long[]{value, len};
        }
        int code = a[pos[0]] & 0xFF;
        pos[0]++;
        long value = code % sigma;
        long len = code / sigma + 1;
        if (len >= runContinues) {
            len += readByteCode(a, pos);
        }
        return new long[]{value, len};
    }

    /** 编码一个游程（value=outrank, len=游程长度，len>=1）。跟 read() 互为逆操作。 */
    public void write(SdsWriter w, long value, long len) {
        if (runContinues == 0) {
            writeByteCode(w, value);
            writeByteCode(w, len - 1);
            return;
        }
        if (len < runContinues) {
            w.rawByte((byte) (value + sigma * (len - 1)));
        } else {
            w.rawByte((byte) (value + sigma * (runContinues - 1)));
            writeByteCode(w, len - runContinues);
        }
    }

    // ---------------- 单节点记录的完整结构 ----------------

    /** 一个节点的记录：出边表 + 游程列表。这是"人类可读"的中间表示，方便你自己的图往这边转换。 */
    public static class NodeRecord {
        /** 每个元素是 {目标 node id, 目标记录内的 rank offset}，顺序即出边的 rank。 */
        public List<long[]> outgoing = new ArrayList<>();
        /** 每个元素是 {outrank, runLength}，按顺序拼起来就是这个节点完整的 BWT 子串。 */
        public List<long[]> runs = new ArrayList<>();
    }

    /** 把一个 NodeRecord 编码进 data 字节流（sigma + 出边表 + 游程体）。 */
    public static void encodeRecord(SdsWriter w, NodeRecord nr) {
        writeByteCode(w, nr.outgoing.size());
        long prev = 0;
        for (long[] edge : nr.outgoing) {
            writeByteCode(w, edge[0] - prev); // 出边目标 node id 增量编码
            prev = edge[0];
            writeByteCode(w, edge[1]);
        }
        long sigma = nr.outgoing.size();
        if (sigma > 0) {
            RunCodec codec = new RunCodec(sigma);
            for (long[] run : nr.runs) codec.write(w, run[0], run[1]);
        }
    }

    /** 从 data[start,limit) 解出一个 NodeRecord。 */
    public static NodeRecord decodeRecord(byte[] data, int start, int limit) {
        int[] pos = {start};
        long sigma = readByteCode(data, pos);
        NodeRecord nr = new NodeRecord();
        long prev = 0;
        for (int i = 0; i < sigma; i++) {
            long to = readByteCode(data, pos) + prev;
            prev = to;
            long rankOffset = readByteCode(data, pos);
            nr.outgoing.add(new long[]{to, rankOffset});
        }
        if (sigma > 0) {
            RunCodec codec = new RunCodec(sigma);
            while (pos[0] < limit) {
                nr.runs.add(codec.read(data, pos));
            }
        }
        return nr;
    }

    /**
     * 零分配地消费一条记录：只做解码和边界检查，不构建 NodeRecord。
     * 返回消费后的结束位置；记录不合法时返回 -1（sigma 异常、游程越界、
     * 或未恰好消费到 limit）。批量校验记录边界时用这个，避免产生大量小对象。
     */
    public static int skipRecord(byte[] data, int start, int limit) {
        int[] pos = {start};
        long sigma = readByteCode(data, pos);
        if (sigma < 0 || sigma > limit - start) return -1; // 出边数不可能超过记录字节数
        for (int i = 0; i < sigma; i++) {
            readByteCode(data, pos); // 目标增量
            readByteCode(data, pos); // rank offset
        }
        if (sigma > 0) {
            long runContinues = (sigma < 255) ? 256 / sigma : 0;
            while (pos[0] < limit) {
                skipRun(data, pos, sigma, runContinues);
            }
        }
        return pos[0] == limit ? pos[0] : -1;
    }

    /** 静态游程消费：与 {@link #read} 相同的步进逻辑，但不创建 RunCodec 实例。 */
    private static void skipRun(byte[] a, int[] pos, long sigma, long runContinues) {
        if (runContinues == 0) {
            readByteCode(a, pos);
            readByteCode(a, pos);
            return;
        }
        int code = a[pos[0]] & 0xFF;
        pos[0]++;
        long len = code / sigma + 1;
        if (len >= runContinues) {
            readByteCode(a, pos);
        }
    }
}
