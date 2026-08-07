package gbwt;

import sds.SdsCodec;
import sds.SdsReader;
import sds.SdsWriter;

import java.io.IOException;

/** 论文 2.3.3 节：GBWT 序列化以 header 开头，定长 6 个 u64 word（48 字节）。 */
public class GBWTHeader implements SdsCodec {
    public static final int TAG = 0x6B376B37;
    public static final long FLAG_BIDIRECTIONAL = 0x0001;
    public static final long FLAG_METADATA = 0x0002;
    public static final long FLAG_SIMPLE_SDS = 0x0004;

    public int tag = TAG;
    public int version = 5; // simple_sds 格式对应的 TAGS_VERSION
    public long sequences;
    public long size;
    public long offset;
    public long alphabetSize;
    public long flags = FLAG_SIMPLE_SDS;

    public boolean isSimpleSds() { return (flags & FLAG_SIMPLE_SDS) != 0; }
    public boolean hasMetadata() { return (flags & FLAG_METADATA) != 0; }
    public boolean isBidirectional() { return (flags & FLAG_BIDIRECTIONAL) != 0; }

    @Override
    public void decode(SdsReader r) throws IOException {
        long w0 = r.u64();
        tag = (int) (w0 & 0xFFFFFFFFL);
        version = (int) (w0 >>> 32);
        sequences = r.u64();
        size = r.u64();
        offset = r.u64();
        alphabetSize = r.u64();
        flags = r.u64();
        if (tag != TAG) {
            throw new IllegalStateException(String.format("GBWTHeader: tag 不匹配 (0x%08X)", tag));
        }
    }

    @Override
    public void encode(SdsWriter w) throws IOException {
        long w0 = (tag & 0xFFFFFFFFL) | ((long) version << 32);
        w.u64(w0);
        w.u64(sequences);
        w.u64(size);
        w.u64(offset);
        w.u64(alphabetSize);
        w.u64(flags);
    }
}
