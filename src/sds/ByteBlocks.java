package sds;

/**
 * 分块字节存储：64 MiB 一个 block 的 byte[][]，突破单个 byte[] 最多 2^31-1 字节的限制。
 * BWT records 在大图上超过 2 GiB（hprc-v2.1 实测 3.16 GiB），单个 byte[] 放不下。
 * 逻辑上它就是一段 length 字节的连续空间，按 (block, offset) 寻址。
 */
public class ByteBlocks {
    public static final int BLOCK_SHIFT = 26;              // 64 MiB
    public static final int BLOCK_SIZE = 1 << BLOCK_SHIFT;
    public static final long BLOCK_MASK = BLOCK_SIZE - 1L;

    public final long length;
    public final byte[][] blocks;

    private ByteBlocks(long length, byte[][] blocks) {
        this.length = length;
        this.blocks = blocks;
    }

    public static ByteBlocks allocate(long length) {
        if (length < 0) throw new IllegalArgumentException("length < 0: " + length);
        int n = (int) ((length + BLOCK_SIZE - 1) >>> BLOCK_SHIFT);
        byte[][] blocks = new byte[n][];
        for (int i = 0; i < n; i++) {
            long remain = length - ((long) i << BLOCK_SHIFT);
            blocks[i] = new byte[(int) Math.min(remain, BLOCK_SIZE)];
        }
        return new ByteBlocks(length, blocks);
    }

    /** 包装一个已有的 byte[]（小于一个 block 的数据不产生拷贝）。 */
    public static ByteBlocks of(byte[] data) {
        if (data.length <= BLOCK_SIZE) {
            return new ByteBlocks(data.length, new byte[][]{data});
        }
        ByteBlocks out = allocate(data.length);
        int copied = 0;
        for (byte[] block : out.blocks) {
            System.arraycopy(data, copied, block, 0, block.length);
            copied += block.length;
        }
        return out;
    }

    public byte get(long i) {
        return blocks[(int) (i >>> BLOCK_SHIFT)][(int) (i & BLOCK_MASK)];
    }

    public void set(long i, byte v) {
        blocks[(int) (i >>> BLOCK_SHIFT)][(int) (i & BLOCK_MASK)] = v;
    }

    /** 取出 [start, end) 区间（用于解码单条 BWT 记录这种小段数据）。 */
    public byte[] copyRange(long start, long end) {
        int len = (int) (end - start);
        byte[] out = new byte[len];
        int copied = 0;
        while (copied < len) {
            long p = start + copied;
            int bi = (int) (p >>> BLOCK_SHIFT);
            int off = (int) (p & BLOCK_MASK);
            int n = Math.min(len - copied, blocks[bi].length - off);
            System.arraycopy(blocks[bi], off, out, copied, n);
            copied += n;
        }
        return out;
    }
}
