package occt.io;

import edu.sysu.pmglab.ecc.ECCReader;
import edu.sysu.pmglab.ecc.record.BoxRecord;
import edu.sysu.pmglab.io.livefile.LiveFile;
import occt.OCCTRecord;

import java.io.IOException;
/**
 * Writes OCCT records using ECC columnar writer.
 *
 * Field layout:
 *   0: node_id_delta     (int64)
 *   1: seq_len            (int32)
 *   2: edge_count         (int16)
 *   3: flags              (int16)
 *   4: sequence           (int32List)  packed DNA words
 *   5: edges              (int64List)  delta-encoded target node IDs
 *   6: edge_orientations  (bytecode)   raw orientation bytes
 *   7: contig_id          (int32)
 *   8: position           (int64)
 *   9: rank_offsets       (int32List)  cumulative path counts
 */
/**
 * Reads OCCT files written by OCCTWriter.
 * Field layout matches OCCTWriter exactly.
 */
public class OCCTReader implements AutoCloseable {

    private final ECCReader reader;
    private final BoxRecord record;
    private long prevNodeId;
    private long recordsRead;

    public OCCTReader(String path) throws IOException {
        this.reader = new ECCReader(LiveFile.of(path));
        this.record = reader.getRecord();
        this.prevNodeId = 0;
        this.recordsRead = 0;
    }

    public OCCTRecord read() throws IOException {
        if (!reader.read(record)) return null;
        OCCTRecord rec = new OCCTRecord();

        long delta = record.get(0);
        rec.nodeId = prevNodeId + delta;
        prevNodeId = rec.nodeId;

        rec.seqLen = record.get(1);
        rec.edgeCount = record.get(2);
        rec.flags = record.get(3);
        rec.contigId = record.get(7);
        rec.position = record.get(8);

        recordsRead++;
        return rec;
    }

    public long recordsRead() { return recordsRead; }

    @Override
    public void close() throws IOException { reader.close(); }
}
