package sds;

import java.io.IOException;

/**
 * GBZ 论文 2.3.1 节："A file is an array of elements: unsigned little-endian 64-bit integers."
 *
 * 读取端抽象。所有 decode() 只依赖这里声明的几个原语，因此底层可以
 * 整文件载入内存（{@link SdsBufferReader}），也可以按模块分块流式读取
 * （{@link SdsStreamReader}），上层的解析代码完全不用变。
 */
public abstract class SdsReader {
    /** JVM 单个数组的安全上限（超过就分配失败，也是 NegativeArraySizeException 的来源）。 */
    public static final long MAX_ARRAY = Integer.MAX_VALUE - 8;

    public abstract long u64();

    public long usize() { return u64(); }

    public abstract boolean hasRemaining();

    /** 已消费的绝对字节数（8 字节对齐和模块定位都靠它）。 */
    public abstract long position();

    /** 剩余可消费字节数；无法确切知道时返回 Long.MAX_VALUE（即不做限制）。 */
    public long bytesRemaining() { return Long.MAX_VALUE; }

    public abstract void skipBytes(long n);

    /** 跳到 8 字节对齐（字符串类向量末尾常有 0~7 字节 padding）。 */
    public void align8() {
        int rem = (int) (position() % 8);
        if (rem != 0) skipBytes(8 - rem);
    }

    public abstract byte[] rawBytes(long nBytes);

    /**
     * 分块读取大块原始字节（BWT records 这类可能超过 2 GiB 的块必须走这里，
     * rawBytes 受 JVM 单数组上限约束，超过 {@link #MAX_ARRAY} 会直接抛错）。
     */
    public abstract ByteBlocks rawBlocks(long nBytes);

    public long[] rawWords(long nWords) {
        checkArraySize(nWords, "rawWords");
        long[] out = new long[(int) nWords];
        for (int i = 0; i < out.length; i++) out[i] = u64();
        return out;
    }

    /** 数组分配前统一检查：小于 0 或超过 JVM 上限时给出明确错误，而不是 NegativeArraySizeException。 */
    public static void checkArraySize(long n, String what) {
        if (n < 0 || n > MAX_ARRAY) {
            throw new IllegalStateException(what + ": 需要分配 " + n
                    + " 个元素，超过 JVM 单数组上限 " + MAX_ARRAY
                    + "（文件损坏/解析错位会读出这种荒谬长度；真实的超大块请改用分块读取 rawBlocks）");
        }
    }

    /** 整文件载入内存（适合小文件；大文件请用 openStream）。 */
    public static SdsReader fromFile(String path) throws IOException {
        return SdsBufferReader.fromFile(path);
    }

    public static SdsReader fromBytes(byte[] data) {
        return SdsBufferReader.fromBytes(data);
    }

    /** 分块流式读取：内存里只保留一个读缓冲块，模块一个一个解析。 */
    public static SdsStreamReader openStream(String path) throws IOException {
        return new SdsStreamReader(path);
    }
}
