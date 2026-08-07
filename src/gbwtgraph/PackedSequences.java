package gbwtgraph;

import edu.sysu.pmglab.container.list.IntDList;
import sds.DNACodec;
import sds.SdsPrimitives.IntVector;
import sds.SdsPrimitives.SparseVector;
import sds.SdsReader;

import java.io.IOException;

/**
 * 大批量节点序列的紧凑内存存储。
 *
 * 每条序列独立用 {@link DNACodec} 打包（含"满/不满"标志位的 int 流），
 * 全部序列的 int 连续存放在 16M-int 的分块数组里（绕开单数组 2^31 上限），
 * 再用 starts/lengths 两个索引数组定位每条序列。
 *
 * 内存对比（hprc-v2.1，2.13 亿条序列）：String[] 需要十几 GB 且 GC 压力大；
 * 这里只有 打包数据(约 1.3 GB) + starts(8 B/条) + lengths(4 B/条)。
 *
 * 访问全部是零分配接口：charAt/byteAt/codeAt 单点取值，
 * decodeInto/fillBytes/copyWords 把结果写进调用方复用的数组。
 */
public class PackedSequences {
    private static final int BLOCK_SHIFT = 24;              // 16M ints = 64 MiB
    private static final int BLOCK_INTS = 1 << BLOCK_SHIFT;
    private static final long BLOCK_MASK = BLOCK_INTS - 1L;

    public final DNACodec.Format format;
    public final int count;
    public final long totalBases;

    private final int[][] blocks;
    private final long[] starts;   // 每条序列在全局 int 流里的起始下标
    private final int[] lengths;

    private PackedSequences(DNACodec.Format format, int count, long totalBases,
                            int[][] blocks, long[] starts, int[] lengths) {
        this.format = format;
        this.count = count;
        this.totalBases = totalBases;
        this.blocks = blocks;
        this.starts = starts;
        this.lengths = lengths;
    }

    // ---------------- 访问 ----------------

    public int length(int i) { return lengths[i]; }

    private int wordAt(long g) {
        return blocks[(int) (g >>> BLOCK_SHIFT)][(int) (g & BLOCK_MASK)];
    }

    /** 第 i 条序列第 pos 个碱基的编码值（A=0,C=1,G=2,T=3,N=4）。 */
    public int codeAt(int i, int pos) {
        long g = starts[i] + pos / format.basesPerInt;
        return (wordAt(g) >>> ((pos % format.basesPerInt) * format.bits)) & format.mask;
    }

    /** 第 i 条序列第 pos 个碱基解析为字符（如 'A'）。 */
    public char charAt(int i, int pos) {
        return DNACodec.charOfCode(codeAt(i, pos));
    }

    /** 第 i 条序列第 pos 个碱基解析为字符对应的 byte 值（如 (byte)'A'）。 */
    public byte byteAt(int i, int pos) {
        return DNACodec.byteOfCode(codeAt(i, pos));
    }

    /** 零分配解码第 i 条序列到 out[outOff..]（调用方复用 out，长度须 >= length(i)）。 */
    public void decodeInto(int i, char[] out, int outOff) {
        int len = lengths[i];
        long g0 = starts[i];
        int bpi = format.basesPerInt, bits = format.bits, mask = format.mask;
        for (int p = 0; p < len; p++) {
            long g = g0 + p / bpi;
            out[outOff + p] = DNACodec.charOfCode(
                    (wordAt(g) >>> ((p % bpi) * bits)) & mask);
        }
    }

    /** 同 {@link #decodeInto}，输出 byte（碱基字符 ASCII 值）。 */
    public void fillBytes(int i, byte[] out, int outOff) {
        int len = lengths[i];
        long g0 = starts[i];
        int bpi = format.basesPerInt, bits = format.bits, mask = format.mask;
        for (int p = 0; p < len; p++) {
            long g = g0 + p / bpi;
            out[outOff + p] = DNACodec.byteOfCode(
                    (wordAt(g) >>> ((p % bpi) * bits)) & mask);
        }
    }

    /** 便捷整串解码（新建 String；热路径请用 decodeInto 复用数组）。 */
    public String decode(int i) {
        char[] out = new char[lengths[i]];
        decodeInto(i, out, 0);
        return new String(out);
    }

    /** 把第 i 条序列的打包 int（含标志位）拷到 out[outOff..]，供下游直接用 DNACodec 处理。 */
    public void copyWords(int i, int[] out, int outOff) {
        int nWords = DNACodec.wordsFor(lengths[i], format);
        long g = starts[i];
        for (int k = 0; k < nWords; k++, g++) {
            out[outOff + k] = wordAt(g);
        }
    }

    public void copyWords(int i, IntDList out, int outOff) {
        int nWords = DNACodec.wordsFor(lengths[i], format);
        long g = starts[i];
        for (int k = 0; k < nWords; k++, g++) {
            out.add(wordAt(g));
        }
    }
    // ---------------- 构建 ----------------

    /**
     * 从文件的 StringArray 结构直接构建（GBWTGraph 节点序列的读取路径）。
     * 全程不产生 String/中间 char 数组：文件里的位压缩字符流逐碱基转写成
     * DNACodec 的带标志位 int 流。alphabet 必须是 ACGTN 的子集。
     */
    public static PackedSequences fromStringArray(SdsReader r) throws IOException {
        SparseVector index = SparseVector.decode(r);
        long alphabetSize = r.usize();
        byte[] compToChar = r.rawBytes(alphabetSize);
        r.align8();
        IntVector chars = IntVector.decode(r);

        // 文件字母表 -> DNACodec 编码；同时定格式（含 N 则 THREE_BIT）
        int[] fileToCode = new int[(int) alphabetSize];
        boolean hasN = false;
        for (int i = 0; i < alphabetSize; i++) {
            fileToCode[i] = DNACodec.codeOf((char) (compToChar[i] & 0xFF));
            if (fileToCode[i] > 3) hasN = true;
        }
        DNACodec.Format format = hasN ? DNACodec.Format.THREE_BIT : DNACodec.Format.TWO_BIT;

        int n = (int) index.size();
        long[] starts = new long[n];
        int[] lengths = new int[n];
        long totalBases = chars.length;

        // 先算每条序列长度和打包后的起始 int 下标（前缀和），算完即可释放 index
        long acc = 0;
        long[] startsOfStrings = index.values;
        for (int i = 0; i < n; i++) {
            long start = startsOfStrings[i];
            long end = (i + 1 < n) ? startsOfStrings[i + 1] : totalBases;
            int len = (int) (end - start);
            lengths[i] = len;
            starts[i] = acc;
            acc += DNACodec.wordsFor(len, format);
        }
        index = null;               // 大文件下这 1.x GB 的 long[] 让 GC 尽早回收
        startsOfStrings = null;

        int nBlocks = (int) ((acc + BLOCK_INTS - 1) >>> BLOCK_SHIFT);
        int[][] blocks = new int[nBlocks][];
        for (int b = 0; b < nBlocks; b++) {
            blocks[b] = new int[(int) Math.min(BLOCK_INTS, acc - ((long) b << BLOCK_SHIFT))];
        }

        // 逐序列转写（不产生中间对象）
        int bpi = format.basesPerInt, bits = format.bits;
        long charPos = 0;
        for (int i = 0; i < n; i++) {
            int len = lengths[i];
            long g = starts[i];
            int slot = 0;
            int word = 0;
            for (int p = 0; p < len; p++, charPos++) {
                int code = fileToCode[(int) chars.get(charPos)];
                word |= code << (slot * bits);
                slot++;
                if (slot == bpi) {
                    blocks[(int) (g >>> BLOCK_SHIFT)][(int) (g & BLOCK_MASK)] = word | DNACodec.FLAG_FULL;
                    g++;
                    slot = 0;
                    word = 0;
                }
            }
            if (slot > 0) {
                blocks[(int) (g >>> BLOCK_SHIFT)][(int) (g & BLOCK_MASK)] = word; // 不满：标志位为 0
            }
        }

        return new PackedSequences(format, n, totalBases, blocks, starts, lengths);
    }

    /** 便捷构建（小数据量/测试用）：把若干字符串打包成一个存储。 */
    public static PackedSequences of(String[] seqs, DNACodec.Format format) {
        int n = seqs.length;
        long[] starts = new long[n];
        int[] lengths = new int[n];
        long acc = 0, totalBases = 0;
        for (int i = 0; i < n; i++) {
            lengths[i] = seqs[i].length();
            starts[i] = acc;
            acc += DNACodec.wordsFor(lengths[i], format);
            totalBases += lengths[i];
        }
        int nBlocks = (int) ((acc + BLOCK_INTS - 1) >>> BLOCK_SHIFT);
        int[][] blocks = new int[nBlocks][];
        for (int b = 0; b < nBlocks; b++) {
            blocks[b] = new int[(int) Math.min(BLOCK_INTS, acc - ((long) b << BLOCK_SHIFT))];
        }
        for (int i = 0; i < n; i++) {
            int[] packed = DNACodec.encode(seqs[i], format);
            long g = starts[i];
            for (int k = 0; k < packed.length; k++, g++) {
                blocks[(int) (g >>> BLOCK_SHIFT)][(int) (g & BLOCK_MASK)] = packed[k];
            }
        }
        return new PackedSequences(format, n, totalBases, blocks, starts, lengths);
    }
}
