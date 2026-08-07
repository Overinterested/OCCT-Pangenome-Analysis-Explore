package sds;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;

/**
 * 与 {@link SdsReader} 对称的写入器。设计成一个可增长的字节缓冲区，
 * u64/align8/writeByteVector 这些是在"原始字节"之上叠加的语义层——
 * 只要保证每个 encode() 方法最终产出的总长度是 8 的倍数（靠 align8 保证），
 * 拼起来的整个文件依然满足 simple-sds 的 element 对齐要求。
 */
public class SdsWriter {
    private final ByteArrayOutputStream buf = new ByteArrayOutputStream();

    public void u64(long v) {
        for (int i = 0; i < 8; i++) {
            buf.write((int) ((v >>> (8 * i)) & 0xFF));
        }
    }

    public void usize(long v) { u64(v); }

    public void rawByte(byte b) { buf.write(b); }

    public void rawBytes(byte[] data) { buf.writeBytes(data); }

    /** 补 0 到 8 字节对齐。 */
    public void align8() {
        int rem = buf.size() % 8;
        if (rem != 0) {
            for (int i = 0; i < 8 - rem; i++) buf.write(0);
        }
    }

    /** Vec&lt;u8&gt;：length(u64) + 原始字节 + padding。 */
    public void writeByteVector(byte[] data) {
        usize(data.length);
        rawBytes(data);
        align8();
    }

    public int size() { return buf.size(); }

    public byte[] toByteArray() {
        if (buf.size() > SdsReader.MAX_ARRAY) {
            throw new IllegalStateException("输出超过 " + SdsReader.MAX_ARRAY
                    + " 字节：当前 SdsWriter 在内存里整段缓冲，写超大 GBZ 需要流式 writer（后续工作）");
        }
        return buf.toByteArray();
    }

    public void writeToFile(String path) throws IOException {
        Files.write(Paths.get(path), toByteArray());
    }

    /** 把一个 SdsCodec 编码进一段独立的字节数组，主要用于 Option&lt;T&gt; 场景（先测量长度再写头）。 */
    public static byte[] capture(SdsCodec obj) throws IOException {
        SdsWriter w = new SdsWriter();
        obj.encode(w);
        return w.toByteArray();
    }
}
