import sds.SdsCodec;
import sds.SdsReader;
import sds.SdsWriter;

import java.io.IOException;

/** 论文 2.3.5 节：GBZ 容器自己的 header，2 个 u64 word（tag+version 打包一个 word，再跟 flags）。 */
public class GBZHeader implements SdsCodec {
    public static final int TAG = 0x205A4247; // "GBZ "

    public int tag = TAG;
    public int version = 1;
    public long flags = 0;

    @Override
    public void decode(SdsReader r) throws IOException {
        long w0 = r.u64();
        tag = (int) (w0 & 0xFFFFFFFFL);
        version = (int) (w0 >>> 32);
        flags = r.u64();
        if (tag != TAG) {
            throw new IllegalStateException(String.format("GBZHeader: tag 不匹配 (0x%08X)，不是有效的 .gbz 文件", tag));
        }
    }

    @Override
    public void encode(SdsWriter w) throws IOException {
        long w0 = (tag & 0xFFFFFFFFL) | ((long) version << 32);
        w.u64(w0);
        w.u64(flags);
    }
}
