package sds;

/**
 * DNA 序列的 int 位压缩存储。
 *
 * 每个 int 的最高位（bit 31，首位）是标志位：1 = 该 int 存满，0 = 没存满
 * （只有序列末尾的最后一个 int 可能没存满）。其余 31 位从低位起依次存放碱基：
 *   第 0 个碱基在 bits [0..w)，第 1 个在 bits [w..2w)，……（w 为位宽）
 *
 * 两种格式：
 *   TWO_BIT   —— 只含 A/C/G/T，2 bit 一个碱基，每个 int 存 15 个（30 bit + 标志位）；
 *   THREE_BIT —— 可能含 N，3 bit 一个碱基，每个 int 存 10 个（30 bit + 标志位）。
 *
 * 碱基编码：A=0, C=1, G=2, T=3, N=4（N 只在 THREE_BIT 里合法）。
 *
 * 解析侧全部提供零分配的快速函数：charAt/byteAt/codeAt 单点取值；
 * 批量解码用 decodeInto/fillBytes 把结果写进调用方提供的复用数组，
 * 避免每条序列都新建 char[]/String（200K 样本规模下这点很重要）。
 */
public final class DNACodec {
    private DNACodec() {}

    public enum Format {
        TWO_BIT(2, 15),
        THREE_BIT(3, 10);

        /** 每个碱基占的 bit 数。 */
        public final int bits;
        /** 每个 int 最多存放的碱基数（31 个可用位 / bits）。 */
        public final int basesPerInt;
        /** 单个碱基的掩码。 */
        public final int mask;

        Format(int bits, int basesPerInt) {
            this.bits = bits;
            this.basesPerInt = basesPerInt;
            this.mask = (1 << bits) - 1;
        }
    }

    /** 标志位：int 首位（bit 31），1 = 存满，0 = 未存满。 */
    public static final int FLAG_FULL = 0x80000000;

    private static final char[] CODE_TO_CHAR = {'A', 'C', 'G', 'T', 'N'};
    private static final byte[] CODE_TO_BYTE = {'A', 'C', 'G', 'T', 'N'};

    /** ASCII -> 编码的查表（大小写都接受），非法字符为 -1。 */
    private static final byte[] CHAR_TO_CODE = new byte[128];
    static {
        java.util.Arrays.fill(CHAR_TO_CODE, (byte) -1);
        CHAR_TO_CODE['A'] = CHAR_TO_CODE['a'] = 0;
        CHAR_TO_CODE['C'] = CHAR_TO_CODE['c'] = 1;
        CHAR_TO_CODE['G'] = CHAR_TO_CODE['g'] = 2;
        CHAR_TO_CODE['T'] = CHAR_TO_CODE['t'] = 3;
        CHAR_TO_CODE['N'] = CHAR_TO_CODE['n'] = 4;
    }

    // ---------------- 编码 ----------------

    /** 存储 length 个碱基需要的 int 数。 */
    public static int wordsFor(int length, Format f) {
        return (length + f.basesPerInt - 1) / f.basesPerInt;
    }

    public static boolean isFull(int packedWord) {
        return (packedWord & FLAG_FULL) != 0;
    }

    /** 碱基字符 -> 编码（大小写均可），非法字符抛 IllegalArgumentException。 */
    public static int codeOf(char c) {
        if (c >= 128 || CHAR_TO_CODE[c] < 0) {
            throw new IllegalArgumentException("非法 DNA 碱基字符: " + c);
        }
        return CHAR_TO_CODE[c];
    }

    /** 编码值 -> 字符（0..4 -> 'A','C','G','T','N'）。 */
    public static char charOfCode(int code) {
        return CODE_TO_CHAR[code];
    }

    /** 编码值 -> 字符对应的 byte 值。 */
    public static byte byteOfCode(int code) {
        return CODE_TO_BYTE[code];
    }

    public static int[] encode(CharSequence seq, Format f) {
        int n = seq.length();
        int[] out = new int[wordsFor(n, f)];
        encodeRange(seq, 0, n, f, out, 0);
        // 标志位：前 n/basesPerInt 个满 int 置 1；末尾 int 视是否正好存满而定
        markFullWords(out, n, f);
        return out;
    }

    /** 同 {@link #encode(CharSequence, Format)}，输入为 ASCII 字节序列。 */
    public static int[] encode(byte[] seq, Format f) {
        int n = seq.length;
        int[] out = new int[wordsFor(n, f)];
        for (int i = 0; i < n; i++) {
            int code = codeOf((char) (seq[i] & 0xFF));
            checkCode(code, f, seq[i]);
            out[i / f.basesPerInt] |= code << ((i % f.basesPerInt) * f.bits);
        }
        markFullWords(out, n, f);
        return out;
    }

    /**
     * 把 seq[off, off+len) 追加编码到 out 里以 baseIndex 起始的位置上（调用方负责
     * out 足够大）。用于把多条序列连续打包进同一个数组，不产生中间对象。
     */
    public static void encodeRange(CharSequence seq, int off, int len, Format f, int[] out, int baseIndex) {
        for (int i = 0; i < len; i++) {
            char c = seq.charAt(off + i);
            int code = codeOf(c);
            checkCode(code, f, (byte) c);
            int g = baseIndex + i;
            out[g / f.basesPerInt] |= code << ((g % f.basesPerInt) * f.bits);
        }
    }

    private static void markFullWords(int[] out, int n, Format f) {
        if (out.length == 0) return;
        int fullWords = n / f.basesPerInt;
        for (int i = 0; i < fullWords; i++) out[i] |= FLAG_FULL;
        if (n % f.basesPerInt == 0) out[out.length - 1] |= FLAG_FULL; // 正好存满
    }

    private static void checkCode(int code, Format f, byte raw) {
        if (f == Format.TWO_BIT && code > 3) {
            throw new IllegalArgumentException("TWO_BIT 格式只接受 A/C/G/T，遇到: " + (char) raw);
        }
    }

    // ---------------- 快速单点解析（零分配） ----------------

    /** 第 i 个碱基的编码值（A=0,C=1,G=2,T=3,N=4）。 */
    public static int codeAt(int[] packed, int i, Format f) {
        return (packed[i / f.basesPerInt] >>> ((i % f.basesPerInt) * f.bits)) & f.mask;
    }

    /** 第 i 个碱基解析为字符（如 'A'）。 */
    public static char charAt(int[] packed, int i, Format f) {
        return CODE_TO_CHAR[codeAt(packed, i, f)];
    }

    /** 第 i 个碱基解析为字符对应的 byte 值（如 (byte)'A'）。 */
    public static byte byteAt(int[] packed, int i, Format f) {
        return CODE_TO_BYTE[codeAt(packed, i, f)];
    }

    // ---------------- 批量解析（复用调用方数组，零分配） ----------------

    /** 解码 length 个碱基到 out[outOff..]，调用方负责复用 out。 */
    public static void decodeInto(int[] packed, int length, Format f, char[] out, int outOff) {
        for (int i = 0; i < length; i++) {
            out[outOff + i] = CODE_TO_CHAR[(packed[i / f.basesPerInt] >>> ((i % f.basesPerInt) * f.bits)) & f.mask];
        }
    }

    /** 同 {@link #decodeInto}，输出 byte（碱基字符的 ASCII 值）。 */
    public static void fillBytes(int[] packed, int length, Format f, byte[] out, int outOff) {
        for (int i = 0; i < length; i++) {
            out[outOff + i] = CODE_TO_BYTE[(packed[i / f.basesPerInt] >>> ((i % f.basesPerInt) * f.bits)) & f.mask];
        }
    }

    /** 便捷整串解码（会新建 String；热路径请用 decodeInto 复用数组）。 */
    public static String decode(int[] packed, int length, Format f) {
        char[] out = new char[length];
        decodeInto(packed, length, f, out, 0);
        return new String(out);
    }
}
