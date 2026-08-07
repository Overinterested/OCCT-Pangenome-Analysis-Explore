package occt;

import sds.DNACodec;
import sds.DNACodec.Format;

/**
 * A single node row in the OCCT (Ordered Compressed Chain Tree) format.
 *
 * Each record stores a node's sequence, outgoing edges, genomic context,
 * structural flags (block boundaries, cycles), and optional path-distribution data.
 */
public class OCCTRecord {

    // ---- Flag bit definitions (16-bit) ----
    public static final short FLAG_BLOCK_START      = (short) (1 << 0);
    public static final short FLAG_BLOCK_END        = (short) (1 << 1);
    public static final short FLAG_CYCLE_TRUNCATED  = (short) (1 << 2);
    public static final short FLAG_REFERENCE        = (short) (1 << 3);
    public static final short FLAG_SENTINEL         = (short) (1 << 4);
    public static final short FLAG_HAS_RANK_OFFSET  = (short) (1 << 5);
    // bits 6-7: orientation_in_reference (0=forward, 1=reverse, 2=both, 3=none)
    public static final short ORIENT_REF_SHIFT = 6;

    // ---- Core fields ----
    public long nodeId;          // 1-based node identifier
    public int  seqLen;          // sequence length in bases
    public int[] sequence;       // packed DNA (DNACodec int[] with FLAG_FULL markers), forward strand only

    // ---- Edge fields ----
    public short edgeCount;
    public long[] edgeTargets;      // [edgeCount] absolute target node IDs
    public byte[] edgeOrientation;  // [edgeCount] packed: (from_end << 1) | to_end

    // ---- Genomic context ----
    public int  contigId;        // index into global contig dictionary, -1 if unknown
    public long position;        // approximate position on contig (0-based bp), -1 if unknown

    // ---- Structural flags ----
    public short flags;

    // ---- Path distribution (optional, corresponds to GBWT BWT substrings) ----
    public int[] rankOffsetArray; // cumulative path counts per edge, null if absent

    // ---- Transient fields (not serialized, used during conversion) ----
    public transient int inDegree;
    public transient int outDegree;

    // ---- Edge orientation constants ----
    /** Standard forward: from right (1) to left (0). */
    public static final byte EDGE_LR = (byte) 0b10;
    /** Target traversed in reverse: from right (1) to right (1). */
    public static final byte EDGE_RR = (byte) 0b11;
    /** Current traversed in reverse: from left (0) to left (0). */
    public static final byte EDGE_LL = (byte) 0b00;
    /** Both reverse: from left (0) to right (1). */
    public static final byte EDGE_LR_REV = (byte) 0b01;

    public OCCTRecord() {}

    /** Minimal constructor for converter use. */
    public OCCTRecord(long nodeId, int[] sequence, int seqLen) {
        this.nodeId = nodeId;
        this.sequence = sequence;
        this.seqLen = seqLen;
        this.contigId = -1;
        this.position = -1;
    }

    public void setEdges(long[] targets, byte[] orientations) {
        this.edgeTargets = targets;
        this.edgeOrientation = orientations;
        this.edgeCount = (short) targets.length;
    }

    public void setEdge(int idx, long target, byte orientation) {
        if (edgeTargets == null || idx >= edgeTargets.length) {
            int newLen = Math.max(idx + 1, edgeTargets == null ? 4 : edgeTargets.length * 2);
            long[] newTargets = new long[newLen];
            byte[] newOrs = new byte[newLen];
            if (edgeTargets != null) {
                System.arraycopy(edgeTargets, 0, newTargets, 0, edgeTargets.length);
                System.arraycopy(edgeOrientation, 0, newOrs, 0, edgeOrientation.length);
            }
            edgeTargets = newTargets;
            edgeOrientation = newOrs;
        }
        edgeTargets[idx] = target;
        edgeOrientation[idx] = orientation;
        edgeCount = (short) Math.max(edgeCount, idx + 1);
    }

    public boolean isBlockStart()     { return (flags & FLAG_BLOCK_START) != 0; }
    public boolean isBlockEnd()       { return (flags & FLAG_BLOCK_END) != 0; }
    public boolean isCycleTruncated() { return (flags & FLAG_CYCLE_TRUNCATED) != 0; }

    public void setBlockStart()  { flags |= FLAG_BLOCK_START; }
    public void setBlockEnd()    { flags |= FLAG_BLOCK_END; }
    public void setCycleTruncated() { flags |= FLAG_CYCLE_TRUNCATED; }

    public int fromEnd(int edgeIdx) { return (edgeOrientation[edgeIdx] >> 1) & 1; }
    public int toEnd(int edgeIdx)   { return edgeOrientation[edgeIdx] & 1; }

    /** Pack orientation. */
    public static byte orientation(int fromEnd, int toEnd) {
        return (byte) ((fromEnd << 1) | toEnd);
    }

    public void setRankOffsetArray(int[] cumulative) {
        this.rankOffsetArray = cumulative;
        if (cumulative != null && cumulative.length > 0) {
            flags |= FLAG_HAS_RANK_OFFSET;
        }
    }

    /** Number of paths taking edge i. */
    public int pathCountForEdge(int edgeIdx) {
        if (rankOffsetArray == null || edgeIdx + 1 >= rankOffsetArray.length) return 0;
        return rankOffsetArray[edgeIdx + 1] - rankOffsetArray[edgeIdx];
    }

    /** Total paths through this node. */
    public int totalPaths() {
        if (rankOffsetArray == null || rankOffsetArray.length == 0) return 0;
        return rankOffsetArray[rankOffsetArray.length - 1];
    }
}
