import gbwt.GBWT;
import gbwt.GBWTHeader;
import gbwt.bwt.BWT;
import gbwt.meta.GBWTMeta;
import gbwt.meta.MetaCotigs;
import gbwt.meta.MetaHeader;
import gbwt.meta.MetaPaths;
import gbwt.meta.MetaSamples;
import gbwtgraph.GBWTGraph;
import gbwtgraph.GBWTGraphHeader;
import gbwtgraph.GBWTGraphSequences;
import gbwtgraph.translation.Segments;
import gbwtgraph.translation.Translation;
import sds.RunCodec;
import sds.SdsWriter;

import java.io.*;
import java.nio.charset.StandardCharsets;
import java.util.*;
import java.util.zip.GZIPInputStream;

/**
 * GFA → GBZ file converter.
 *
 * <p>Converts a pangenome GFA file (from tools like minigraph-cactus, pggb)
 * into the GBZ binary format used by vg giraffe for haplotype-aware read
 * mapping.  The converter handles three dimensions:
 *
 * <ol>
 *   <li><b>Topology</b> — GFA S-lines become GBWT nodes; L-lines become
 *       outgoing edges in each node's BWT record.</li>
 *   <li><b>Paths</b> — Links are grouped by {@code SR:i:N} tag.
 *       Each SR value traces one haplotype walk; these walks are encoded
 *       as BWT runs (one run per walk per node).</li>
 *   <li><b>Sequences</b> — Segment sequences are packed into
 *       GBWTGraphSequences (2-bit per base for ACGT, 3-bit if N is present).</li>
 * </ol>
 *
 * <p>Usage:
 * <pre>{@code
 *   GFA2GBZ.convert("input.gfa", "output.gbz");
 * }</pre>
 */
public class GFA2GBZ {

    // ---- Internal data structures -------------------------------------

    /** One GFA segment, parsed from an S-line. */
    private static class GfaSegment {
        String id;
        String sequence;
        final Map<String, String> tags = new HashMap<>();
        int nodeId = -1; // assigned GBWT node ID (0-based internal)
    }

    /** One GFA link, parsed from an L-line. */
    private static class GfaLink {
        String fromId;
        boolean fromForward;
        String toId;
        boolean toForward;
        String overlap;
        final Map<String, String> tags = new HashMap<>();
    }

    // ---- Public API ---------------------------------------------------

    /**
     * Convert a GFA file to GBZ format.
     *
     * @param gfaPath  path to the input GFA file (.gfa or .gfa.gz)
     * @param gbzPath  path for the output GBZ file
     */
    public static void convert(String gfaPath, String gbzPath) throws IOException {
        long t0 = System.currentTimeMillis();

        // ---- Phase 1: parse GFA ----
        List<GfaSegment> segments = new ArrayList<>();
        List<GfaLink> links = new ArrayList<>();
        parseGFA(gfaPath, segments, links);
        System.out.println("[GFA→GBZ] Parsed " + segments.size()
                + " segments, " + links.size() + " links");

        // ---- Phase 2: assign node IDs ----
        Map<String, Integer> segIdToNode = assignNodeIds(segments);

        // ---- Phase 3: build path index (SR → {segmentId → LinkTarget}) ----
        Map<Integer, Map<String, LinkTarget>> srIndex = buildSRIndex(links);
        int[] sortedRanks = srIndex.keySet().stream().sorted().mapToInt(i -> i).toArray();

        // ---- Phase 4: trace paths and collect per-node path info ----
        // nodePaths[seg.nodeId] = list of (outrank, pathIndex) for paths through this node
        int totalPaths = sortedRanks.length;
        List<List<int[]>> nodePaths = new ArrayList<>(segments.size());
        for (int i = 0; i < segments.size(); i++) nodePaths.add(new ArrayList<>());
        long totalPathBp = 0;

        for (int pi = 0; pi < sortedRanks.length; pi++) {
            int sr = sortedRanks[pi];
            Map<String, LinkTarget> srLinks = srIndex.get(sr);
            // Find start node
            Set<String> targets = new HashSet<>();
            for (LinkTarget lt : srLinks.values()) targets.add(lt.toId);
            String start = null;
            for (String from : srLinks.keySet()) {
                if (!targets.contains(from)) { start = from; break; }
            }
            if (start == null) continue;

            String cur = start;
            Set<String> visited = new HashSet<>();
            while (cur != null && visited.add(cur)) {
                Integer nid = segIdToNode.get(cur);
                if (nid == null) break;
                GfaSegment seg = segments.get(nid);
                // Determine outrank for this path at this node
                // Build the out-edge list for this node and find rank of SR-specific edge
                List<String> orderedTargets = getOutTargets(links, cur);
                LinkTarget lt = srLinks.get(cur);
                int outrank = lt != null ? orderedTargets.indexOf(lt.toId) : 0;
                if (outrank < 0) outrank = 0;
                nodePaths.get(nid).add(new int[]{outrank, pi});
                totalPathBp += seg.sequence.length();
                cur = lt != null ? lt.toId : null;
            }
        }

        // ---- Phase 5: build BWT NodeRecords ----
        List<RunCodec.NodeRecord> nodeRecords = new ArrayList<>(segments.size());
        for (int i = 0; i < segments.size(); i++) {
            GfaSegment seg = segments.get(i);
            RunCodec.NodeRecord rec = new RunCodec.NodeRecord();
            // Outgoing edges from this node
            List<String> targetList = getOutTargets(links, seg.id);
            for (String tgt : targetList) {
                Integer tgtNode = segIdToNode.get(tgt);
                if (tgtNode != null) {
                    rec.outgoing.add(new long[]{tgtNode, 0});
                }
            }
            // BWT runs: one run per path through this node, sorted by path index
            List<int[]> pathInfo = nodePaths.get(i);
            pathInfo.sort(Comparator.comparingInt(a -> a[1])); // sort by path index
            if (!pathInfo.isEmpty()) {
                int curRank = pathInfo.get(0)[0];
                int runLen = 1;
                for (int j = 1; j < pathInfo.size(); j++) {
                    int rank = pathInfo.get(j)[0];
                    if (rank == curRank) {
                        runLen++;
                    } else {
                        rec.runs.add(new long[]{curRank, runLen});
                        curRank = rank;
                        runLen = 1;
                    }
                }
                rec.runs.add(new long[]{curRank, runLen});
            }
            nodeRecords.add(rec);
        }

        System.out.println("[GFA→GBZ] " + totalPaths + " paths, "
                + totalPathBp + " total path bp");

        // ---- Phase 6: build sequences ----
        String[] sequences = new String[segments.size()];
        for (int i = 0; i < segments.size(); i++) {
            sequences[i] = segments.get(i).sequence;
        }

        // ---- Phase 7: assemble GBZ ----
        GBZFile gbz = new GBZFile();

        // GBZHeader
        gbz.header = new GBZHeader();
        gbz.header.tag = GBZHeader.TAG;
        gbz.header.version = 1;
        gbz.header.flags = 0;

        // GBZTags - store source info
        gbz.tags = new GBZTags();
        gbz.tags.entries.put("source", "GFA2GBZ");

        // GBWT
        gbz.gbwt = new GBWT();
        gbz.gbwt.header = new GBWTHeader();
        gbz.gbwt.header.tag = GBWTHeader.TAG;
        gbz.gbwt.header.version = 5;
        gbz.gbwt.header.offset = 1; // nodes start from 1
        gbz.gbwt.header.sequences = totalPaths; // number of paths × 1 (single-stranded for haplotype walking)
        gbz.gbwt.header.alphabetSize = nodeRecords.size() + 1; // nodes + sentinel
        gbz.gbwt.header.flags = GBWTHeader.FLAG_SIMPLE_SDS | GBWTHeader.FLAG_BIDIRECTIONAL;

        // BWT
        gbz.gbwt.bwt = BWT.fromNodeRecords(nodeRecords);
        gbz.gbwt.header.size = gbz.gbwt.bwt.records.data.length;

        // Metadata (store path -> SR mapping)
        gbz.gbwt.meta = new GBWTMeta();
        gbz.gbwt.meta.present = true;
        gbz.gbwt.meta.metaHeader = new MetaHeader();
        gbz.gbwt.meta.metaHeader.tag = MetaHeader.TAG;
        gbz.gbwt.meta.metaHeader.version = 2;
        gbz.gbwt.meta.metaHeader.sampleCount = 1;
        gbz.gbwt.meta.metaHeader.haplotypeCount = totalPaths;
        gbz.gbwt.meta.metaHeader.contigCount = 1;
        gbz.gbwt.meta.metaHeader.flags = 0;

        gbz.gbwt.meta.metaSamples = new MetaSamples();
        gbz.gbwt.meta.metaSamples.names = new String[]{"GFA_SAMPLE"};

        gbz.gbwt.meta.metaCotigs = new MetaCotigs();
        // Use SN tag from first segment as contig name
        String contigName = "chr1";
        if (!segments.isEmpty()) {
            String sn = segments.get(0).tags.get("SN");
            if (sn != null && !sn.isEmpty()) contigName = sn;
        }
        gbz.gbwt.meta.metaCotigs.names = new String[]{contigName};

        gbz.gbwt.meta.metaPaths = new MetaPaths();
        for (int i = 0; i < totalPaths; i++) {
            MetaPaths.PathName pn = new MetaPaths.PathName();
            pn.sample = 0;
            pn.contig = 0;
            pn.phase = sortedRanks[i];
            pn.count = 0;
            gbz.gbwt.meta.metaPaths.paths.add(pn);
        }

        gbz.gbwt.header.flags |= GBWTHeader.FLAG_METADATA;

        // GBWTGraph
        gbz.graph = new GBWTGraph();
        gbz.graph.header = new GBWTGraphHeader();
        gbz.graph.header.tag = GBWTGraphHeader.TAG;
        gbz.graph.header.version = 3;
        gbz.graph.header.nodes = sequences.length;
        gbz.graph.header.flags = GBWTGraphHeader.FLAG_SIMPLE_SDS;

        gbz.graph.sequences = new GBWTGraphSequences();
        gbz.graph.sequences.forward = sequences;
        gbz.graph.sequences.zstd = false;

        // Translation: segment name → node mapping
        gbz.graph.translation = new Translation();
        gbz.graph.translation.segments = new Segments();
        gbz.graph.translation.segments.names = new String[segments.size()];
        for (int i = 0; i < segments.size(); i++) {
            gbz.graph.translation.segments.names[i] = segments.get(i).id;
        }

        // ---- Phase 8: write ----
        System.out.println("[GFA→GBZ] Writing " + gbzPath + " ...");
        gbz.parseTo(gbzPath);
        long t1 = System.currentTimeMillis();
        System.out.printf("[GFA→GBZ] Done in %.1f s%n", (t1 - t0) / 1000.0);
    }

    // ---- GFA parser ---------------------------------------------------

    private static final int DEFAULT_BUF = 1 << 20;

    private static void parseGFA(String path, List<GfaSegment> segments,
                                  List<GfaLink> links) throws IOException {
        try (BufferedReader br = openReader(path)) {
            String line;
            while ((line = br.readLine()) != null) {
                if (line.isEmpty()) continue;
                char type = line.charAt(0);
                switch (type) {
                    case 'S' -> {
                        GfaSegment s = new GfaSegment();
                        int p1 = line.indexOf('\t', 2);
                        int p2 = line.indexOf('\t', p1 + 1);
                        s.id = line.substring(2, p1);
                        s.sequence = line.substring(p1 + 1, p2);
                        parseTags(line, p2 + 1, s.tags);
                        segments.add(s);
                    }
                    case 'L' -> {
                        GfaLink l = new GfaLink();
                        int p1 = line.indexOf('\t', 2);
                        int p2 = line.indexOf('\t', p1 + 1);
                        int p3 = line.indexOf('\t', p2 + 1);
                        int p4 = line.indexOf('\t', p3 + 1);
                        int p5 = line.indexOf('\t', p4 + 1);
                        l.fromId = line.substring(2, p1);
                        l.fromForward = line.charAt(p1 + 1) == '+';
                        l.toId = line.substring(p2 + 1, p3);
                        l.toForward = line.charAt(p3 + 1) == '+';
                        if (p5 > 0) {
                            l.overlap = line.substring(p4 + 1, p5);
                            parseTags(line, p5 + 1, l.tags);
                        } else {
                            l.overlap = line.substring(p4 + 1);
                        }
                        links.add(l);
                    }
                }
            }
        }
    }

    private static void parseTags(String line, int start, Map<String, String> tags) {
        int len = line.length();
        int pos = start;
        while (pos < len) {
            int next = line.indexOf('\t', pos);
            if (next < 0) next = len;
            String field = line.substring(pos, next);
            int c1 = field.indexOf(':');
            int c2 = field.indexOf(':', c1 + 1);
            if (c1 >= 0 && c2 >= 0) {
                tags.put(field.substring(0, c1), field.substring(c2 + 1));
            }
            pos = next + 1;
        }
    }

    private static BufferedReader openReader(String path) throws IOException {
        InputStream in = new BufferedInputStream(new FileInputStream(path), DEFAULT_BUF);
        if (path.endsWith(".gz")) {
            in = new GZIPInputStream(in, DEFAULT_BUF);
        }
        return new BufferedReader(new InputStreamReader(in, StandardCharsets.UTF_8), DEFAULT_BUF);
    }

    // ---- Helpers ------------------------------------------------------

    private static Map<String, Integer> assignNodeIds(List<GfaSegment> segments) {
        Map<String, Integer> map = new LinkedHashMap<>();
        for (int i = 0; i < segments.size(); i++) {
            map.put(segments.get(i).id, i);
            segments.get(i).nodeId = i;
        }
        return map;
    }

    private static Map<Integer, Map<String, LinkTarget>> buildSRIndex(List<GfaLink> links) {
        Map<Integer, Map<String, LinkTarget>> index = new HashMap<>();
        for (GfaLink l : links) {
            String srStr = l.tags.get("SR");
            int sr = srStr != null ? Integer.parseInt(srStr) : 0;
            index.computeIfAbsent(sr, k -> new HashMap<>())
                 .put(l.fromId, new LinkTarget(l.toId, l.toForward));
        }
        return index;
    }

    private static List<String> getOutTargets(List<GfaLink> links, String segmentId) {
        List<String> targets = new ArrayList<>();
        for (GfaLink l : links) {
            if (l.fromId.equals(segmentId) && !targets.contains(l.toId)) {
                targets.add(l.toId);
            }
        }
        return targets;
    }

    // ---- Inner types --------------------------------------------------

    private static class LinkTarget {
        final String toId;
        final boolean toForward;
        LinkTarget(String toId, boolean toForward) {
            this.toId = toId; this.toForward = toForward;
        }
    }

    // ---- Entry point (for standalone use) -----------------------------

    public static void main(String[] args) throws IOException {
        if (args.length < 2) {
            System.err.println("Usage: java GFA2GBZ <input.gfa> <output.gbz>");
            System.exit(1);
        }
        convert(args[0], args[1]);
    }
}
