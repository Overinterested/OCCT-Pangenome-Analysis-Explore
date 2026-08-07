package gbwt;

import gbwt.bwt.BWT;
import gbwt.meta.GBWTMeta;
import sds.SdsCodec;
import sds.SdsReader;
import sds.SdsWriter;

import java.io.IOException;

/**
 * 论文 2.3.3 节 "GBWT"：
 *   header -> tags -> bwt(压缩 BWT) -> Option&lt;DASamples&gt; -> Option&lt;Metadata&gt;
 * 顺序固定，改了顺序文件就读不回来了。
 */
public class GBWT implements SdsCodec {
    public GBWTHeader header = new GBWTHeader();
    public GBWTTags tags = new GBWTTags();
    public BWT bwt = new BWT();
    public GBWTDASamples samples = new GBWTDASamples();
    public GBWTMeta meta = new GBWTMeta();

    @Override
    public void decode(SdsReader r) throws IOException {
        header.decode(r);
        if (!header.isSimpleSds()) {
            throw new IllegalStateException("这是旧版 SDSL 格式（非 simple_sds），本项目只覆盖现代 simple_sds 格式。");
        }

        tags = new GBWTTags();
        tags.decode(r);

        bwt = new BWT();
        bwt.decode(r);

        samples = new GBWTDASamples();
        samples.decode(r);

        meta = new GBWTMeta();
        meta.decode(r);

        if (meta.present != header.hasMetadata()) {
            throw new IllegalStateException("header.hasMetadata() 和实际读到的 metadata 是否存在不一致，文件可能损坏");
        }
    }

    @Override
    public void encode(SdsWriter w) throws IOException {
        // 写之前保持 header 的 flags 跟实际内容一致，避免手滑忘了同步
        header.flags = (header.flags & ~GBWTHeader.FLAG_METADATA)
                | (meta.present ? GBWTHeader.FLAG_METADATA : 0)
                | GBWTHeader.FLAG_SIMPLE_SDS;

        header.encode(w);
        tags.encode(w);
        bwt.encode(w);
        samples.encode(w);
        meta.encode(w);
    }
}
