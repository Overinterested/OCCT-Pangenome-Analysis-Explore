package occt;

/**
 * OCCT file header.
 *
 * Stored at the beginning of every .occt file. Contains format identification,
 * version, node count, and global metadata sufficient to initialize the reader.
 */
public class OCCTHeader {
    /** Magic tag: "OCCT" in little-endian ASCII = 0x5443434F */
    public static final int TAG = 0x5443434F;

    public int  tag     = TAG;
    public int  version = 1;
    public long nodeCount;       // total number of node records in the file
    public long blockCount;      // total number of traversal blocks
    public long totalBases;      // sum of all sequence lengths
    public long totalEdges;      // sum of all outgoing edge counts
    public int  formatVersion;   // DNACodec format: 0=TWO_BIT, 1=THREE_BIT

    /** Contig dictionary: contigNames[i] is the name for contigId=i. */
    public String[] contigNames;

    public OCCTHeader() {}

    public boolean hasContigNames() {
        return contigNames != null && contigNames.length > 0;
    }

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        sb.append("OCCT v").append(version).append(", ");
        sb.append("nodes=").append(String.format("%,d", nodeCount));
        sb.append(", blocks=").append(String.format("%,d", blockCount));
        sb.append(", bases=").append(String.format("%,d", totalBases));
        sb.append(", edges=").append(String.format("%,d", totalEdges));
        return sb.toString();
    }
}
