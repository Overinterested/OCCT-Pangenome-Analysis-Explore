package gbwt.meta;

import sds.SdsCodec;
import sds.SdsReader;
import sds.SdsWriter;

import java.io.IOException;

/** Metadata 的定长头部，6 个 u64 word：tag+version 打包一个 word，再跟 3 个计数 + flags。 */
public class MetaHeader implements SdsCodec {
    public static final int TAG = 0x6B375E7A;

    public int tag = TAG;
    public int version = 2;
    public long sampleCount;
    public long haplotypeCount;
    public long contigCount;
    public long flags;

    @Override
    public void decode(SdsReader r) throws IOException {
        long w0 = r.u64();
        tag = (int) (w0 & 0xFFFFFFFFL);
        version = (int) (w0 >>> 32);
        sampleCount = r.u64();
        haplotypeCount = r.u64();
        contigCount = r.u64();
        flags = r.u64();
        if (tag != TAG) {
            throw new IllegalStateException(String.format("MetaHeader: tag 不匹配 (0x%08X)", tag));
        }
    }

    @Override
    public void encode(SdsWriter w) throws IOException {
        long w0 = (tag & 0xFFFFFFFFL) | ((long) version << 32);
        w.u64(w0);
        w.u64(sampleCount);
        w.u64(haplotypeCount);
        w.u64(contigCount);
        w.u64(flags);
    }
}
