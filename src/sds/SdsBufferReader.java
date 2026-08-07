package sds;

import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.channels.FileChannel;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;

/** 整文件载入内存的 SdsReader（原 SdsReader 的实现，行为不变）。 */
public class SdsBufferReader extends SdsReader {
    private final ByteBuffer buf;

    private SdsBufferReader(ByteBuffer bb) { this.buf = bb; }

    public static SdsBufferReader fromFile(String path) throws IOException {
        try (FileChannel ch = FileChannel.open(Paths.get(path), StandardOpenOption.READ)) {
            ByteBuffer bb = ByteBuffer.allocate((int) ch.size());
            while (bb.hasRemaining()) { if (ch.read(bb) < 0) break; }
            bb.flip();
            bb.order(ByteOrder.LITTLE_ENDIAN);
            return new SdsBufferReader(bb);
        }
    }

    public static SdsBufferReader fromBytes(byte[] data) {
        return new SdsBufferReader(ByteBuffer.wrap(data).order(ByteOrder.LITTLE_ENDIAN));
    }

    @Override
    public long u64() { return buf.getLong(); }

    @Override
    public boolean hasRemaining() { return buf.hasRemaining(); }

    @Override
    public long position() { return buf.position(); }

    @Override
    public long bytesRemaining() { return buf.remaining(); }

    @Override
    public void skipBytes(long n) { buf.position((int) (buf.position() + n)); }

    @Override
    public byte[] rawBytes(long nBytes) {
        checkArraySize(nBytes, "rawBytes");
        byte[] out = new byte[(int) nBytes];
        buf.get(out);
        return out;
    }

    @Override
    public ByteBlocks rawBlocks(long nBytes) {
        ByteBlocks out = ByteBlocks.allocate(nBytes);
        for (byte[] block : out.blocks) {
            buf.get(block);
        }
        return out;
    }
}
