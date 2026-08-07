package gbwtgraph.translation;

import sds.SdsCodec;
import sds.SdsReader;
import sds.SdsWriter;
import sds.SdsPrimitives.SparseVector;

import java.io.IOException;
import java.util.Arrays;

/**
 * 论文 2.3.4 节："a sparse bitvector B mapping node ranges to segments ...
 * Segment Si is the concatenation of nodes v in [B.select(i,1) ... B.select(i+1,1))"。
 * 也就是说 nodeToSegment.values[i] 是第 i 个 segment 对应的第一个 node id。
 */
public class Mapping implements SdsCodec {
    public SparseVector nodeToSegment;

    /** 给定一个 node id，找它属于第几个 segment（下标对应 Segments.names）。 */
    public int segmentIndexForNode(long nodeId) {
        int idx = Arrays.binarySearch(nodeToSegment.values, nodeId);
        if (idx < 0) idx = -idx - 2; // 找最后一个 <= nodeId 的起始位置
        return idx;
    }

    @Override
    public void decode(SdsReader r) throws IOException {
        nodeToSegment = SparseVector.decode(r);
    }

    @Override
    public void encode(SdsWriter w) throws IOException {
        nodeToSegment.encode(w);
    }
}
