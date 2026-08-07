package gbwtgraph.translation;

import sds.SdsCodec;
import sds.SdsReader;
import sds.SdsWriter;

import java.io.IOException;

/**
 * 论文 2.3.4 节 "node-to-segment translation"：
 * 顺序固定是先 segments 后 mapping。这两块在文件里总是存在（不是 Option 包裹的），
 * 没有做过 translation 时（GBWTGraphHeader.hasTranslation()==false）它们就是空的。
 */
public class Translation implements SdsCodec {
    public Segments segments = new Segments();
    public Mapping mapping = new Mapping();

    /** 节点所在的 segment 名字。没有 translation 时约定 segment 名就是 node id 本身。 */
    public String segmentNameForNode(long nodeId, boolean translationEnabled) {
        if (!translationEnabled) return String.valueOf(nodeId);
        return segments.names[mapping.segmentIndexForNode(nodeId)];
    }

    @Override
    public void decode(SdsReader r) throws IOException {
        segments = new Segments();
        segments.decode(r);
        mapping = new Mapping();
        mapping.decode(r);
    }

    @Override
    public void encode(SdsWriter w) throws IOException {
        segments.encode(w);
        mapping.encode(w);
    }
}
