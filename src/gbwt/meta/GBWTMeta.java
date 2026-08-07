package gbwt.meta;

import sds.SdsCodec;
import sds.SdsReader;
import sds.SdsWriter;

import java.io.IOException;

/**
 * 论文 2.3.3 节："Metadata starts with a header, followed by a vector of path names.
 * Sample and contig names are serialized as dictionaries."
 * 顺序固定：header -> path_names -> sample_names -> contig_names。
 *
 * 这是 Option&lt;Metadata&gt;：GBWT.header.hasMetadata()==false 时，present 应该是 false，
 * decode/encode 自己处理开头的长度前缀。
 */
public class GBWTMeta implements SdsCodec {
    public boolean present = false;

    public MetaHeader metaHeader = new MetaHeader();
    public MetaPaths metaPaths = new MetaPaths();
    public MetaSamples metaSamples = new MetaSamples();
    public MetaCotigs metaCotigs = new MetaCotigs();

    @Override
    public void decode(SdsReader r) throws IOException {
        long elements = r.usize();
        if (elements == 0) {
            present = false;
            return;
        }
        if (elements < 0 || elements > r.bytesRemaining() / 8) {
            throw new IllegalStateException("Metadata: element 计数 " + elements
                    + " 不合法（剩余可消费字节 " + r.bytesRemaining() + "）——文件损坏或上游解析错位");
        }
        present = true;
        metaHeader.decode(r);
        metaPaths.decode(r);
        metaSamples.decode(r);
        metaCotigs.decode(r);
    }

    @Override
    public void encode(SdsWriter w) throws IOException {
        if (!present) {
            w.usize(0);
            return;
        }
        SdsWriter sub = new SdsWriter();
        metaHeader.encode(sub);
        metaPaths.encode(sub);
        metaSamples.encode(sub);
        metaCotigs.encode(sub);
        byte[] bytes = sub.toByteArray();
        w.usize(bytes.length / 8L);
        w.rawBytes(bytes);
    }
}
