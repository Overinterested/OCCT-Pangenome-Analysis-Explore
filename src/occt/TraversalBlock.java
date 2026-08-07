package occt;

import java.util.ArrayList;
import java.util.List;

/**
 * A traversal block: a maximal subgraph interval in the chain-tree DFS order,
 * bounded by degree-change nodes.
 *
 * Linear blocks correspond to unitigs (all internal nodes have degree 1).
 * Branching blocks correspond to relaxed superbubbles (snarl-like regions).
 *
 * The block is the unit of I/O: nodes within a block are loaded/written together.
 */
public class TraversalBlock {
    public int blockId;                  // sequential block ID
    public long startNodeId;             // first node in this block
    public long endNodeId;               // last node in this block (may equal startNodeId for single-node blocks)
    public boolean isLinear;             // true if all internal nodes have inDegree=outDegree=1
    public boolean isCyclic;             // true if this block forms a cycle

    /** Ordered list of OCCTRecords in this block. */
    public final List<OCCTRecord> nodes;

    /** Optional: child blocks for nested structures. */
    public final List<TraversalBlock> children;

    public TraversalBlock() {
        this.nodes = new ArrayList<>();
        this.children = new ArrayList<>();
        this.blockId = -1;
    }

    public TraversalBlock(int blockId) {
        this();
        this.blockId = blockId;
    }

    public void addNode(OCCTRecord node) {
        if (nodes.isEmpty()) {
            startNodeId = node.nodeId;
        }
        nodes.add(node);
        endNodeId = node.nodeId;
    }

    public int size() { return nodes.size(); }
    public boolean isEmpty() { return nodes.isEmpty(); }

    /** Compute block flags for each node based on position within the block. */
    public void finalizeBlock() {
        if (nodes.isEmpty()) return;
        nodes.get(0).setBlockStart();
        nodes.get(nodes.size() - 1).setBlockEnd();

        // Determine if linear
        isLinear = true;
        for (int i = 0; i < nodes.size(); i++) {
            OCCTRecord n = nodes.get(i);
            if (n.inDegree > 1 || n.outDegree > 1) {
                isLinear = false;
            }
        }
    }

    /** Flatten this block and all children into a single ordered list. */
    public List<OCCTRecord> flatten() {
        List<OCCTRecord> result = new ArrayList<>();
        flattenInto(result);
        return result;
    }

    private void flattenInto(List<OCCTRecord> out) {
        out.addAll(nodes);
        for (TraversalBlock child : children) {
            child.flattenInto(out);
        }
    }

    @Override
    public String toString() {
        return String.format("Block[%d] nodes=%d children=%d %s%s [%d..%d]",
                blockId, nodes.size(), children.size(),
                isLinear ? "linear" : "branching",
                isCyclic ? " cyclic" : "",
                startNodeId, endNodeId);
    }
}
