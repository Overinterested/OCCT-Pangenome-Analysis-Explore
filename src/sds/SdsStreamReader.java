package sds;

import edu.sysu.pmglab.bytecode.ByteStream;
import edu.sysu.pmglab.io.reader.ChannelReaderStream;

import java.io.EOFException;
import java.io.IOException;
import java.io.UncheckedIOException;

/**
 * 分块流式 SdsReader：底层是 ChannelReaderStream，上层用一个 ByteStream 做读缓冲块，
 * 读取模式就是示例里的 "fill cache -> 消费 -> clear -> 再 fill" 循环，
 * 只是把它封装进了 u64()/rawBytes()/skipBytes() 这些原语里。
 *
 * 两个实现要点：
 *   1) 用 ChannelReaderStream 而不是 ReaderStream：ReaderStream 内部还有一层私有缓冲，
 *      seek 时无法精确换算逻辑位置；ChannelReaderStream 每次 read 调用之间不留存字节，
 *      tell()/seek() 就是精确的文件位置，大块跳过可以直接 seek，不必逐字节读完。
 *   2) rawBytes(n) 在剩余需求量大于缓冲块容量时直接读进目标数组，不在 cache 里倒第二手，
 *      这样 BWT records 这种大块数据也能流式进内存，峰值内存就是目标数组本身。
 */
public class SdsStreamReader extends SdsReader implements AutoCloseable {
    private static final int BUFFER_SIZE = 1 << 16; // 64 KB

    private final ChannelReaderStream channel;
    private final ByteStream cache;
    private long position;
    private boolean eof;

    public SdsStreamReader(String path) throws IOException {
        this.cache = new ByteStream(BUFFER_SIZE);
        this.channel = new ChannelReaderStream(path);
    }

    // ---------------- 缓冲块管理 ----------------

    /** cache 已全部消费时调用：清空并重新填满（或读到文件尾）。 */
    private void refill() {
        cache.clear();
        byte[] b = cache.bytes();
        int total = 0;
        try {
            while (total < b.length) {
                int got = channel.read(b, total, b.length - total);
                if (got < 0) break;
                total += got;
            }
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
        if (total == 0) eof = true;
        cache.wSeek(total);
    }

    /** 保证 cache 里至少有 n 个未消费字节（要求 n <= 缓冲块容量）。 */
    private void ensure(int n) {
        while (cache.rRemaining() < n) {
            if (eof) throw eof();
            // 把没消费完的数据 compact 到开头，再继续读
            byte[] b = cache.bytes();
            int rem = cache.rRemaining();
            System.arraycopy(b, cache.offset() + cache.rTell(), b, 0, rem);
            cache.rSeek(0);
            cache.wSeek(rem);
            try {
                while (rem < b.length) {
                    int got = channel.read(b, rem, b.length - rem);
                    if (got < 0) { eof = true; break; }
                    rem += got;
                }
            } catch (IOException e) {
                throw new UncheckedIOException(e);
            }
            cache.wSeek(rem);
        }
    }

    private UncheckedIOException eof() {
        return new UncheckedIOException(new EOFException("文件在偏移 " + position + " 处意外结束"));
    }

    // ---------------- SdsReader 原语 ----------------

    @Override
    // TODO: optime
    public long u64() {
        ensure(8);
        byte[] b = cache.bytes();
        int p = cache.offset() + cache.rTell();
        long v = (b[p] & 0xFFL)
                | ((b[p + 1] & 0xFFL) << 8)
                | ((b[p + 2] & 0xFFL) << 16)
                | ((b[p + 3] & 0xFFL) << 24)
                | ((b[p + 4] & 0xFFL) << 32)
                | ((b[p + 5] & 0xFFL) << 40)
                | ((b[p + 6] & 0xFFL) << 48)
                | ((b[p + 7] & 0xFFL) << 56);
        cache.rSkip(8);
        position += 8;
        return v;
    }

    @Override
    public boolean hasRemaining() {
        if (cache.rRemaining() == 0 && !eof) refill();
        return cache.rRemaining() > 0;
    }

    @Override
    public long position() {
        return position;
    }

    /** 文件中尚未消费的字节数（按逻辑位置算，含 cache 里未读走的部分）。 */
    @Override
    public long bytesRemaining() {
        try {
            return channel.length() - position;
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }

    @Override
    public void skipBytes(long n) {
        if (n < 0) {
            throw new IllegalStateException("skipBytes: 负长度 " + n
                    + "（element 计数溢出或文件损坏），当前偏移 " + position);
        }
        long rest = n;
        int take = (int) Math.min(rest, cache.rRemaining());
        cache.rSkip(take);
        rest -= take;
        if (rest > 0) {
            // cache 已空，此时 channel 的物理位置 == 逻辑位置，直接 seek
            try {
                long target = channel.tell() + rest;
                if (target > channel.length()) {
                    throw new UncheckedIOException(new EOFException(
                            "skipBytes(" + n + ") 越过文件末尾（偏移 " + position
                                    + "，文件总长 " + channel.length() + "）——文件损坏或上游解析错位"));
                }
                channel.seek(channel.tell() + rest);
            } catch (IOException e) {
                throw new UncheckedIOException(e);
            }
        }
        position += n;
    }

    @Override
    public byte[] rawBytes(long nBytes) {
        checkArraySize(nBytes, "rawBytes");
        byte[] out = new byte[(int) nBytes];
        int copied = 0;
        while (copied < out.length) {
            if (cache.rRemaining() == 0) {
                int remaining = out.length - copied;
                if (remaining >= cache.bytes().length) {
                    // 大块数据直接读进目标数组，绕过 cache
                    try {
                        while (remaining > 0) {
                            int got = channel.read(out, copied, remaining);
                            if (got < 0) { eof = true; throw eof(); }
                            copied += got;
                            remaining -= got;
                        }
                    } catch (IOException e) {
                        throw new UncheckedIOException(e);
                    }
                } else {
                    refill();
                    if (cache.rRemaining() == 0) throw eof();
                }
            } else {
                copied += cache.read(out, copied, out.length - copied);
            }
        }
        position += nBytes;
        return out;
    }

    @Override
    public ByteBlocks rawBlocks(long nBytes) {
        ByteBlocks out = ByteBlocks.allocate(nBytes);
        for (byte[] block : out.blocks) {
            int copied = 0;
            while (copied < block.length) {
                if (cache.rRemaining() == 0) {
                    refill();
                    if (cache.rRemaining() == 0) throw eof();
                } else {
                    copied += cache.read(block, copied, block.length - copied);
                }
            }
        }
        position += nBytes;
        return out;
    }


    /**
     * 直接跳到指定逻辑位置：丢弃当前缓冲，seek 底层 channel，重置 position。
     * 用于 GBWTGraph header 不存在时回退到 sequences 起始位置。
     */
    public void seekToPosition(long targetPosition) throws IOException {
        if (targetPosition < 0 || targetPosition > channel.length()) {
            throw new IllegalArgumentException("seekToPosition: 越界 " + targetPosition
                    + "（文件长度 " + channel.length() + "）");
        }
        channel.seek(targetPosition);
        position = targetPosition;
        cache.clear();
        eof = false;
    }

    @Override
    public void close() throws IOException {
        channel.close();
        cache.close();
    }
}
