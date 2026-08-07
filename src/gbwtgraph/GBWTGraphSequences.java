package gbwtgraph;

import sds.SdsCodec;
import sds.SdsReader;
import sds.SdsWriter;
import sds.StringCodec;

import java.io.IOException;

/**
 * 论文 2.3.4 节："While the in-memory data structure ... stores the labels in both
 * orientations for faster access, serializing the reverse orientation is clearly
 * unnecessary." 所以文件里只存正向序列，下标 = nodeId - 1（节点从 1 开始编号）；
 * 具体压缩方式由 GBWTGraphHeader 决定：
 *   version 3（如 y.giraffe.gbz）：普通的字母表压缩 StringArray；
 *   version 4+ 且 FLAG_ZSTD 置位：zstd 压缩整块字节流
 *   （论文里叫 compress_even，因为在内存模型里正向节点用偶数下标）。
 */
public class GBWTGraphSequences implements SdsCodec {
    /** 下标 = nodeId - 1 的正向序列。 */
    public String[] forward = new String[0];

    /** 打包存储形式（大文件用）：设置后优先于 forward 使用，forward 留空以省内存。 */
    public PackedSequences packed = null;

    /** 文件里的存储格式，decode 前由 GBWTGraph 按 header 设置；encode 时按同一格式写回。 */
    public boolean zstd = false;

    public String getForward(long nodeId) {
        return packed != null ? packed.decode((int) (nodeId - 1)) : forward[(int) (nodeId - 1)];
    }

    /** 序列条数（两种存储形式通用）。 */
    public long count() {
        return packed != null ? packed.count : forward.length;
    }

    public String getSequence(long nodeId, boolean isForward) {
        String fwd = getForward(nodeId);
        return isForward ? fwd : reverseComplement(fwd);
    }

    public static String reverseComplement(String s) {
        StringBuilder sb = new StringBuilder(s.length());
        for (int i = s.length() - 1; i >= 0; i--) sb.append(complement(s.charAt(i)));
        return sb.toString();
    }

    private static char complement(char c) {
        switch (c) {
            case 'A': return 'T'; case 'T': return 'A';
            case 'C': return 'G'; case 'G': return 'C';
            case 'a': return 't'; case 't': return 'a';
            case 'c': return 'g'; case 'g': return 'c';
            default: return c;
        }
    }

    @Override
    public void decode(SdsReader r) throws IOException {
        forward = zstd ? StringCodec.decodeEvenCompressed(r) : StringCodec.decodeStringArray(r);
    }

    @Override
    public void encode(SdsWriter w) throws IOException {
        String[] seqs = forward;
        if (packed != null) {
            // 从打包存储还原出字符串数组再编码
            // （几亿条序列的整体写回应等流式 SdsWriter，目前只适合中小规模文件）
            seqs = new String[packed.count];
            for (int i = 0; i < seqs.length; i++) seqs[i] = packed.decode(i);
        }
        if (zstd) StringCodec.encodeEvenCompressed(w, seqs);
        else StringCodec.encodeStringArray(w, seqs);
    }
}
