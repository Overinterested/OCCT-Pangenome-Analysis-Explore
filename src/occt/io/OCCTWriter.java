package occt.io;

import edu.sysu.pmglab.bytecode.ByteStream;
import edu.sysu.pmglab.bytecode.Bytes;
import edu.sysu.pmglab.container.list.IntDList;
import edu.sysu.pmglab.container.list.LongDList;
import edu.sysu.pmglab.ecc.ECCWriter;
import edu.sysu.pmglab.ecc.record.BoxRecord;
import edu.sysu.pmglab.ecc.type.FieldType;

import java.io.IOException;

/**
 * Writes OCCT records using ECC columnar writer.
 * <p>
 * Field layout:
 * 0: node_id_delta     (int64)
 * 1: seq_len            (int32)
 * 2: edge_count         (int16)
 * 3: flags              (int16)
 * 4: sequence           (int32List)  packed DNA words
 * 5: edges              (int64List)  delta-encoded target node IDs
 * 6: edge_orientations  (bytecode)   raw orientation bytes
 * 7: contig_id          (int32)
 * 8: position           (int64)
 * 9: rank_offsets       (int32List)  cumulative path counts
 */
public class OCCTWriter implements AutoCloseable {

    private final ECCWriter writer;
    private final BoxRecord record;
    private long prevNodeId;
    private long recordsWritten;
    private final LongDList deltas = LongDList.init();
    private final ByteStream ors = new ByteStream();

    public OCCTWriter(String path) throws IOException {
        this.writer = ECCWriter.setOutput(path)
                .addField("node_id_delta", FieldType.int64)
                .addField("seq_len", FieldType.int32)
                .addField("edge_count", FieldType.int16)
                .addField("flags", FieldType.int16)
                .addField("sequence", FieldType.int32List)
                .addField("edges", FieldType.int64List)
                .addField("edge_orientations", FieldType.bytecode)
                .addField("contig_id", FieldType.int32)
                .addField("position", FieldType.int64)
                .addField("rank_offsets", FieldType.int32List)
                .instance();
        this.record = writer.getRecord();
        this.prevNodeId = 0;
        this.recordsWritten = 0;
    }

    /**
     * Write a single node record.
     *
     * @param nodeId      absolute node ID
     * @param seqLen      sequence length in bases
     * @param edgeCnt     number of outgoing edges
     * @param flags       bitfield (block start/end, cycle, etc.)
     * @param packedSeq   packed DNA words (int[] from DNACodec), may be null
     * @param edgeTargets absolute target node IDs
     * @param edgeOrs     orientation bytes per edge
     * @param contigId    contig index, -1 if unknown
     * @param position    contig position, -1 if unknown
     * @param rankOffsets cumulative path counts per edge, null if absent
     */
    public void writeRecord(long nodeId, int seqLen, short edgeCnt, short flags,
                            IntDList packedSeq, long[] edgeTargets, byte[] edgeOrs,
                            int contigId, long position, int[] rankOffsets) throws IOException {
        record.clear();

        // Field 0: delta node_id
        record.set(0, nodeId - prevNodeId);
        prevNodeId = nodeId;

        // Field 1: seq_len
        record.set(1, seqLen);

        // Field 2: edge_count
        record.set(2, (int) edgeCnt);

        // Field 3: flags
        record.set(3, (int) flags);

        // Field 4: sequence (packed DNA as int32List)
        if (packedSeq != null && !packedSeq.isEmpty()) {
            record.set(4, packedSeq);
        } else {
            record.set(4, IntDList.EMPTY());
        }

        // Field 5: edges (delta-encoded as int64List)
        if (edgeCnt > 0 && edgeTargets != null) {
            long prev = nodeId;
            for (int i = 0; i < edgeCnt; i++) {
                deltas.add(edgeTargets[i] - prev);
                prev = edgeTargets[i];
            }
            record.set(5, deltas);
        } else {
            record.set(5, LongDList.EMPTY());
        }

        // Field 6: edge_orientations (raw bytes)
        if (edgeCnt > 0 && edgeOrs != null) {
            for (int i = 0; i < edgeCnt; i++) {
                ors.write(edgeOrs[i]);
            }
            record.set(6, ors.asBytes());
        } else {
            record.set(6, Bytes.EMPTY);
        }

        // Field 7: contig_id
        record.set(7, contigId);

        // Field 8: position
        record.set(8, position);

        // Field 9: rank_offsets (int32List)
        if (rankOffsets != null && rankOffsets.length > 0) {
            record.set(9, IntDList.wrap(rankOffsets));
        } else {
            record.set(9, IntDList.EMPTY());
        }

        writer.write(record);
        record.clear();
        deltas.clear();
        ors.clear();
        if (packedSeq != null) packedSeq.clear();
        recordsWritten++;
    }

    public long recordsWritten() {
        return recordsWritten;
    }

    @Override
    public void close() throws IOException {
        writer.close();
    }
}
