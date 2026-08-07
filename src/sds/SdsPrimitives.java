package sds;

/**
 * 论文 2.3.2 节 "Building blocks" 的三层地基：
 *   IntVector（位压缩定长整数数组）
 *   RawBitVector（普通位图，用作 SparseVector 的 high 部分）
 *   SparseVector（Elias-Fano 稀疏位图，论文 2.1.2/2.2 节）
 *
 * 这几个是"值对象"，不是可变的领域对象，所以用 静态 decode()/build() 工厂 + 实例 encode()
 * 的模式，而不是 SdsCodec 那种"先 new 空对象再 decode 填充"的模式。
 */
public class SdsPrimitives {

    // ---------------- IntVector ----------------

    public static class IntVector {
        public long length;
        public int width; // 0..64（0 是合法的边界情况：每个元素占 0 位，即所有值都隐含为 0）
        public long[] words;

        public static IntVector decode(SdsReader r) {
            long length = r.usize();
            long width = r.usize();
            long bitSize = r.usize();
            long wordCount = r.usize();
            if (length * width != bitSize) throw new IllegalStateException("IntVector: bitSize 字段不一致");
            long expected = (bitSize + 63) / 64;
            if (wordCount != expected) throw new IllegalStateException("IntVector: wordCount 字段不一致");
            IntVector v = new IntVector();
            v.length = length;
            v.width = (int) width;
            v.words = r.rawWords(wordCount);
            return v;
        }

        public void encode(SdsWriter w) {
            long bitSize = length * width;
            long wordCount = (bitSize + 63) / 64;
            w.usize(length);
            w.usize(width);
            w.usize(bitSize);
            w.usize(wordCount);
            for (long word : words) w.u64(word);
        }

        public long get(long i) {
            if (width == 0) return 0;
            long bitOff = i * width;
            int w0 = (int) (bitOff >>> 6);
            int b0 = (int) (bitOff & 63);
            long lo = words[w0] >>> b0;
            if (b0 + width > 64) lo |= words[w0 + 1] << (64 - b0);
            if (width < 64) lo &= (1L << width) - 1;
            return lo;
        }

        public void set(long i, long value) {
            if (width == 0) return;
            long bitOff = i * width;
            int w0 = (int) (bitOff >>> 6);
            int b0 = (int) (bitOff & 63);
            long mask = (width == 64) ? -1L : ((1L << width) - 1);
            value &= mask;
            words[w0] &= ~(mask << b0);
            words[w0] |= (value << b0);
            if (b0 + width > 64) {
                int rem = b0 + width - 64;
                long maskHi = (1L << rem) - 1;
                words[w0 + 1] &= ~maskHi;
                words[w0 + 1] |= (value >>> (64 - b0));
            }
        }

        /** 从一组值构建一个位宽固定为 width 的 IntVector（用于写文件时打包数据）。 */
        public static IntVector pack(long[] values, int width) {
            IntVector v = new IntVector();
            v.length = values.length;
            v.width = width;
            long bitSize = (long) values.length * width;
            long wordCount = (bitSize + 63) / 64;
            v.words = new long[(int) wordCount];
            for (int i = 0; i < values.length; i++) v.set(i, values[i]);
            return v;
        }

        /** 表示 maxValue 所需的最小位宽（至少 1 位）。 */
        public static int widthFor(long maxValue) {
            if (maxValue <= 0) return 1;
            return Math.max(64 - Long.numberOfLeadingZeros(maxValue), 1);
        }
    }

    // ---------------- RawBitVector ----------------

    public static class RawBitVector {
        public long ones;
        public long length;
        public long[] words;

        public static RawBitVector decode(SdsReader r) {
            long ones = r.usize();
            long bitLength = r.usize();
            long wordCount = r.usize();
            long expected = (bitLength + 63) / 64;
            if (wordCount != expected) throw new IllegalStateException("BitVector: wordCount 字段不一致");
            RawBitVector v = new RawBitVector();
            v.ones = ones;
            v.length = bitLength;
            v.words = r.rawWords(wordCount);
            skipEmptyOption(r); // rank support（从不写进文件，加载后按需重建）
            skipEmptyOption(r); // select_1 support
            skipEmptyOption(r); // select_0 support
            return v;
        }

        public void encode(SdsWriter w) {
            long wordCount = (length + 63) / 64;
            w.usize(ones);
            w.usize(length);
            w.usize(wordCount);
            for (long word : words) w.u64(word);
            w.usize(0); // 3 个恒为空的 option
            w.usize(0);
            w.usize(0);
        }

        public boolean get(long i) {
            return ((words[(int) (i >>> 6)] >>> (i & 63)) & 1L) != 0;
        }

        private long[] cumPopcount;
        private void ensureIndex() {
            if (cumPopcount != null) return;
            cumPopcount = new long[words.length];
            long acc = 0;
            for (int i = 0; i < words.length; i++) {
                acc += Long.bitCount(words[i]);
                cumPopcount[i] = acc;
            }
        }

        /** 第 rank 个置位比特的位置（rank 从 1 开始）。 */
        public long selectOne(long rank) {
            ensureIndex();
            int lo = 0, hi = cumPopcount.length - 1, word = hi;
            while (lo <= hi) {
                int mid = (lo + hi) >>> 1;
                if (cumPopcount[mid] >= rank) { word = mid; hi = mid - 1; }
                else lo = mid + 1;
            }
            long before = word == 0 ? 0 : cumPopcount[word - 1];
            long remaining = rank - before;
            long w = words[word];
            for (int b = 0; b < 64; b++) {
                if (((w >>> b) & 1L) != 0) {
                    remaining--;
                    if (remaining == 0) return ((long) word << 6) + b;
                }
            }
            throw new IllegalStateException("selectOne: rank 越界");
        }

        /** 从一组置位位置构建一个长度为 length 的位图。 */
        public static RawBitVector fromPositions(long length, long[] positions) {
            RawBitVector v = new RawBitVector();
            v.length = length;
            v.ones = positions.length;
            v.words = new long[(int) ((length + 63) / 64)];
            for (long p : positions) v.words[(int) (p >>> 6)] |= (1L << (p & 63));
            return v;
        }

        static void skipEmptyOption(SdsReader r) {
            long elements = r.usize();
            if (elements > 0) r.skipBytes(elements * 8);
        }
    }

    // ---------------- SparseVector（Elias-Fano） ----------------

    public static class SparseVector {
        public long universe;
        public long[] values; // 升序排列的整数集合

        public long size() { return values.length; }

        public static SparseVector decode(SdsReader r) {
            long universe = r.usize();
            RawBitVector high = RawBitVector.decode(r);
            IntVector low = IntVector.decode(r);
            int w = low.width;
            long m = low.length;
            long[] values = new long[(int) m];
            for (long i = 0; i < m; i++) {
                long highPart = high.selectOne(i + 1) - i;
                values[(int) i] = low.get(i) | (highPart << w);
            }
            SparseVector v = new SparseVector();
            v.universe = universe;
            v.values = values;
            return v;
        }

        public void encode(SdsWriter w) {
            w.usize(universe);
            if (values.length == 0) {
                // C++ 端写出的空稀疏向量：high 是全空的位图，low 是 width=64 的空 int vector
                RawBitVector.fromPositions(0, new long[0]).encode(w);
                IntVector.pack(new long[0], 64).encode(w);
                return;
            }
            int width = chooseWidth();
            int m = values.length;
            long[] positions = new long[m];
            long[] lowValues = new long[m];
            long mask = (width == 64) ? -1L : ((1L << width) - 1);
            for (int i = 0; i < m; i++) {
                long v = values[i];
                positions[i] = (v >>> width) + i;
                lowValues[i] = v & mask;
            }
            // high 位图长度 = 桶数 ceil(universe / 2^width) + 值个数（Elias-Fano 一元码部分），
            // 不是"最后一个位置 + 1"——C++ 端按整个 universe 预留桶。
            long buckets = (width >= 64) ? 1 : ((universe + (1L << width) - 1) >>> width);
            long highLength = buckets + m;
            RawBitVector.fromPositions(highLength, positions).encode(w);
            IntVector.pack(lowValues, width).encode(w);
        }

        /**
         * 低位宽度 = max(1, floor(log2(universe / size)))，下限是 1 位——
         * simple-sds 对很密的集合也至少留 1 个低位比特。
         */
        private int chooseWidth() {
            long avgGap = Math.max(universe / values.length, 1);
            int w = 1;
            while (w < 62 && (1L << (w + 1)) <= avgGap) w++;
            return w;
        }

        /**
         * 直接从一组已排序的值构建（构建时机通常是"我知道每个记录的起始偏移了，打个索引"）。
         * universe 取 C++ 端同样的约定：最后一个值 + 1（空集为 0），
         * 而不是整个底层数据的总长度——这样写出来的字节才能和 vg/gbwt 完全一致。
         */
        public static SparseVector build(long[] sortedValues) {
            SparseVector v = new SparseVector();
            v.universe = (sortedValues.length == 0) ? 0 : sortedValues[sortedValues.length - 1] + 1;
            v.values = sortedValues;
            return v;
        }
    }
}
