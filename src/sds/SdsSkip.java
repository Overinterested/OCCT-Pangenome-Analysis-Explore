package sds;

/**
 * 结构化跳过：每种 simple-sds 结构的字节长度都能从它的长度字段算出来，
 * 所以"跳过"只需要读几个 u64 头字段，数据部分直接 skipBytes（流式下是 seek），
 * 不需要把数据实例化到内存里。这是 GBZStreamer 各 skipX 方法的基础。
 */
public final class SdsSkip {
    private SdsSkip() {}

    /** IntVector：length, width, bitSize, wordCount, words[wordCount]。 */
    public static void intVector(SdsReader r) {
        r.usize(); // length
        r.usize(); // width
        r.usize(); // bitSize
        long words = r.usize();
        r.skipBytes(words * 8);
    }

    /** 普通位图：ones, bits, wordCount, words[wordCount]，再跟 rank/select_1/select_0 三个 option。 */
    public static void bitVector(SdsReader r) {
        r.usize(); // ones
        r.usize(); // bits
        long words = r.usize();
        r.skipBytes(words * 8);
        for (int i = 0; i < 3; i++) {
            long elements = r.usize();
            r.skipBytes(elements * 8);
        }
    }

    /** SparseVector（Elias-Fano）：universe + 高位位图 + 低位 IntVector。 */
    public static void sparseVector(SdsReader r) {
        r.usize(); // universe
        bitVector(r);
        intVector(r);
    }

    /** StringArray：起始位置 SparseVector + 字母表（长度前缀 + padding）+ 字符 IntVector。 */
    public static void stringArray(SdsReader r) {
        sparseVector(r);
        long alphabetSize = r.usize();
        r.skipBytes(alphabetSize);
        r.align8();
        intVector(r);
    }

    /** Dictionary（sample/contig 名字典）：StringArray + 字典序下标 IntVector。 */
    public static void dictionary(SdsReader r) {
        stringArray(r);
        intVector(r);
    }

    /** Vec&lt;u8&gt;：长度前缀 + 原始字节 + padding。 */
    public static void byteVector(SdsReader r) {
        long len = r.usize();
        r.skipBytes(len);
        r.align8();
    }
}
