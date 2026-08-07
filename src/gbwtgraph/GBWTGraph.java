package gbwtgraph;

import gbwtgraph.translation.Translation;
import sds.SdsCodec;
import sds.SdsReader;
import sds.SdsWriter;

import java.io.IOException;

/**
 * 论文 2.3.4 节 "GBWTGraph"："it only needs to store a header and node labels" ——
 * 拓扑（谁连谁）全部来自 GBWT 那份压缩 BWT，GBWTGraph 自己只多存序列和 segment 转换关系。
 * 顺序固定：header -> sequences(zstd) -> segments -> node_to_segment（这两个在 Translation 里）。
 */
public class GBWTGraph implements SdsCodec {
    public GBWTGraphHeader header = new GBWTGraphHeader();
    public GBWTGraphSequences sequences = new GBWTGraphSequences();
    public Translation translation = new Translation();

    public long nodeCount() { return header.nodes; }

    public String getSequence(long nodeId, boolean forward) {
        return sequences.getSequence(nodeId, forward);
    }

    public String segmentNameForNode(long nodeId) {
        return translation.segmentNameForNode(nodeId, header.hasTranslation());
    }

    @Override
    public void decode(SdsReader r) throws IOException {
        header.decode(r);
        sequences = new GBWTGraphSequences();
        sequences.zstd = header.usesZstdSequences();
        sequences.decode(r);
        translation = new Translation();
        translation.decode(r);
    }

    @Override
    public void encode(SdsWriter w) throws IOException {
        header.nodes = sequences.count();
        header.encode(w);
        sequences.encode(w);
        translation.encode(w);
    }
}
