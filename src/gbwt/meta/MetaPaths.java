package gbwt.meta;

import sds.SdsCodec;
import sds.SdsReader;
import sds.SdsWriter;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * 论文 2.1.6 节："sample identifier, contig identifier, haplotype/phase identifier, fragment identifier"。
 *
 * 注意：simple_sds 序列化时每个 PathName 是 4 个 uint32（共 16 字节），
 * 不是 4 个 u64。C++ 端用 32 位字段的 ShortPathName 落盘
 * （gbwt metadata.cpp），所以两个 u32 打包在一个 8 字节 element 里。
 * 这决定了 metadata 区后续所有字段的对齐位置，写成 u64 会导致
 * sample/contig 字典整体错位 3*16 字节。
 */
public class MetaPaths implements SdsCodec {

    public static class PathName {
        public long sample, contig, phase, count;

        public PathName() {}

        public PathName(long sample, long contig, long phase, long count) {
            this.sample = sample; this.contig = contig; this.phase = phase; this.count = count;
        }
    }

    /** 32 位全 1 是 gbwt 的"无此值"哨兵（如 _gbwt_ref 参考路径的 phase）。 */
    public static final long NO_VALUE = 0xFFFFFFFFL;

    public List<PathName> paths = new ArrayList<>();

    @Override
    public void decode(SdsReader r) throws IOException {
        long count = r.usize();
        paths = new ArrayList<>();
        for (int i = 0; i < count; i++) {
            long w0 = r.u64(); // sample(低32位) | contig(高32位)
            long w1 = r.u64(); // phase (低32位) | count (高32位)
            paths.add(new PathName(
                    w0 & 0xFFFFFFFFL, w0 >>> 32,
                    w1 & 0xFFFFFFFFL, w1 >>> 32));
        }
    }

    @Override
    public void encode(SdsWriter w) throws IOException {
        w.usize(paths.size());
        for (PathName p : paths) {
            w.u64((p.sample & 0xFFFFFFFFL) | ((p.contig & 0xFFFFFFFFL) << 32));
            w.u64((p.phase & 0xFFFFFFFFL) | ((p.count & 0xFFFFFFFFL) << 32));
        }
    }
}
