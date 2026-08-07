package gbwtgraph;

import sds.SdsCodec;
import sds.SdsReader;
import sds.SdsWriter;

import java.io.IOException;

/**
 * 论文 2.3.4 节：GBWTGraph 的 header，3 个 u64 word（tag+version 打包一个 word，再跟 nodes、flags）。
 *
 * 注意 flags 的 0x0002 这一位含义随版本变化：
 *   version 3：只是 simple_sds 标记，节点序列恒为普通的字母表压缩 StringArray；
 *   version 4+：作为 FLAG_ZSTD，置位时节点序列是 zstd 压缩的 even 字符串数组。
 * 实测依据：y.giraffe.gbz 是 version=3、flags=0x2 的文件，其序列区逐字节符合
 * 普通 StringArray 布局（字母表 ACGT + 位压缩字符），并不是 zstd。
 */
public class GBWTGraphHeader implements SdsCodec {
    public static final int TAG = 0x6B3764AF;
    public static final long FLAG_TRANSLATION = 0x0001;
    public static final long FLAG_SIMPLE_SDS = 0x0002;
    public static final long FLAG_ZSTD = 0x0002; // version >= 4 时的含义
    public static final int VERSION_ZSTD_SEQUENCES = 3;

    public int tag = TAG;
    public int version = 4;
    public long nodes;
    public long flags = FLAG_SIMPLE_SDS;

    public boolean hasTranslation() { return (flags & FLAG_TRANSLATION) != 0; }

    /** 节点序列是否用 zstd 压缩（version 3 及更早的文件一律是普通 StringArray）。 */
    public boolean usesZstdSequences() {
        return version > VERSION_ZSTD_SEQUENCES && (flags & FLAG_ZSTD) != 0;
    }

    @Override
    public void decode(SdsReader r) throws IOException {
        long w0 = r.u64();
        tag = (int) (w0 & 0xFFFFFFFFL);
        version = (int) (w0 >>> 32);
        nodes = r.u64();
        flags = r.u64();
        if (tag != TAG) {
            throw new IllegalStateException(String.format("GBWTGraphHeader: tag 不匹配 (0x%08X)", tag));
        }
    }

    @Override
    public void encode(SdsWriter w) throws IOException {
        long w0 = (tag & 0xFFFFFFFFL) | ((long) version << 32);
        w.u64(w0);
        w.u64(nodes);
        w.u64(flags);
    }
}
