import edu.sysu.pmglab.container.list.IntDList;
import gbwt.GBWT;
import gbwt.bwt.BWT;
import gbwtgraph.PackedSequences;
import gbwtgraph.GBWTGraphSequences;
import occt.io.OCCTWriter;
import sds.DNACodec;
import sds.RunCodec;

import java.io.*;
import java.util.*;
import java.util.concurrent.*;
import java.util.concurrent.atomic.AtomicLong;

public class GBZ2OCCT {

    private final String gbzPath, outputDir;
    private BWT bwt;
    private PackedSequences sequences;
    private long nodeCount;

    private BitSet inDegGt1, inDegZero;

    static class ContigRange {
        final String name;
        final long startNode, endNode;
        ContigRange(String n, long s, long e) { name = n; startNode = s; endNode = e; }
    }

    private List<ContigRange> contigs;
    private final AtomicLong globalRecords = new AtomicLong();
    private final AtomicLong globalBlocks = new AtomicLong();
    private final AtomicLong globalBases = new AtomicLong();
    private final AtomicLong globalEdges = new AtomicLong();
    private final int numThreads;

    public GBZ2OCCT(String gbzPath, String outputDir, int numThreads) {
        this.gbzPath = gbzPath; this.outputDir = outputDir; this.numThreads = numThreads;
    }

    public void convert() throws IOException {
        log("[OCCT] Loading GBZ: " + gbzPath);
        loadGBZ();
        log("[OCCT] Detecting contig boundaries from sentinel...");
        detectContigBoundaries();
        log("[OCCT] Computing in-degree map...");
        computeInDegGt1();
        int threads = numThreads > 0 ? numThreads : Runtime.getRuntime().availableProcessors();
        log("[OCCT] Parallel mode: " + contigs.size() + " ranges, " + threads + " threads");
        convertParallel(threads);
        printStats();
        log("[OCCT] Done. Output: " + outputDir);
    }

    private void loadGBZ() throws IOException {
        try (GBZStreamer s = new GBZStreamer(gbzPath)) {
            s.readGBZHeader(); s.skipGBZTags();
            GBWT gbwt = s.readGBWT();
            this.bwt = gbwt.bwt;
            gbwt = null; System.gc();
            log(String.format("  BWT: %,d records, mem=%d MB", bwt.recordCount(), usedMemMB()));
            s.readGraphHeader();
            GBWTGraphSequences seq = s.readSequencesPacked();
            this.sequences = seq.packed; this.nodeCount = seq.packed.count;
            seq = null; System.gc();
            log(String.format("  seq: %,d nodes, %,d bases, mem=%d MB", nodeCount, sequences.totalBases, usedMemMB()));
            s.skipTranslation();
        }
    }

    private void detectContigBoundaries() {
        long[] sr = new long[2];
        bwt.index.rangeInto(0, bwt.records.data.length, sr);
        if (sr[1] <= sr[0]) {
            contigs = new ArrayList<>();
            contigs.add(new ContigRange("full", 1, nodeCount));
            log("  sentinel empty, using single range");
            return;
        }
        int sl = (int)(sr[1] - sr[0]);
        byte[] buf = new byte[sl];
        copyFromBlocks(bwt.records.data, sr[0], buf, sl);
        int[] sp = {0};
        long sigma = RunCodec.readByteCode(buf, sp);
        long[] endpoints = new long[(int) sigma];
        long prev = 0;
        for (int e = 0; e < sigma; e++) {
            long to = RunCodec.readByteCode(buf, sp) + prev;
            prev = to; RunCodec.readByteCode(buf, sp);
            endpoints[e] = to >> 1;
        }
        Arrays.sort(endpoints);
        log(String.format("  sentinel children: %,d endpoints", endpoints.length));
        long maxGap = Math.max(1000, nodeCount / 500);
        contigs = new ArrayList<>();
        long rangeStart = 1;
        int contigIdx = 0;
        for (int i = 0; i < endpoints.length; i++) {
            if (i > 0 && endpoints[i] - endpoints[i-1] > maxGap) {
                long contigEnd = endpoints[i-1];
                if (contigEnd >= rangeStart) {
                    contigs.add(new ContigRange("ctg_" + String.format("%05d", contigIdx), rangeStart, contigEnd));
                    contigIdx++;
                }
                rangeStart = contigEnd + 1;
            }
        }
        if (rangeStart <= nodeCount) {
            contigs.add(new ContigRange("ctg_" + String.format("%05d", contigIdx), rangeStart, nodeCount));
        }
        log(String.format("  detected %,d contig ranges (gap threshold=%,d)", contigs.size(), maxGap));
    }

    private void computeInDegGt1() {
        inDegGt1 = new BitSet((int) nodeCount + 1);
        inDegZero = new BitSet((int) nodeCount + 1);
        inDegZero.set(1, (int) nodeCount + 1);
        byte[] inDegCount = new byte[(int) nodeCount + 1];
        long[] range = new long[2];
        byte[] buf = new byte[1 << 18];
        int totalRecords = (int) bwt.recordCount();
        int ri = Math.max(1, totalRecords / 20);
        for (int i = 0; i < totalRecords; i++) {
            bwt.index.rangeInto(i, bwt.records.data.length, range);
            if (range[1] <= range[0]) continue;
            int len = (int)(range[1] - range[0]);
            if (buf.length < len) buf = new byte[Math.max(len, buf.length * 2)];
            copyFromBlocks(bwt.records.data, range[0], buf, len);
            int[] pos = {0}; long sigma = RunCodec.readByteCode(buf, pos);
            long prev = 0;
            for (int e = 0; e < sigma; e++) {
                long to = RunCodec.readByteCode(buf, pos) + prev;
                prev = to; RunCodec.readByteCode(buf, pos);
                long tid = to >> 1;
                if (tid > 0 && tid < inDegCount.length) {
                    int t = (int)tid; inDegZero.clear(t);
                    byte c = inDegCount[t];
                    if (c < 2) { inDegCount[t] = (byte)(c+1); if (c==1) inDegGt1.set(t); }
                }
            }
            if (i > 0 && i % ri == 0)
                log(String.format("  in-deg: %,d / %,d  mem=%d MB", i, totalRecords, usedMemMB()));
        }
        System.gc();
        log(String.format("  inDeg>1=%,d  inDeg=0=%,d  mem=%d MB", inDegGt1.cardinality(), inDegZero.cardinality(), usedMemMB()));
    }

    private void convertParallel(int threads) throws IOException {
        new File(outputDir).mkdirs();
        ExecutorService pool = Executors.newFixedThreadPool(threads);
        List<Future<?>> futures = new ArrayList<>();
        contigs.sort((a, b) -> Long.compare(b.endNode - b.startNode, a.endNode - a.startNode));
        for (ContigRange cr : contigs) {
            String outPath = outputDir + "/" + sanitizeFilename(cr.name) + ".occt";
            futures.add(pool.submit(() -> {
                try { dfsContig(cr.startNode, cr.endNode, outPath, cr.name); }
                catch (IOException e) { throw new RuntimeException(e); }
            }));
        }
        pool.shutdown();
        for (Future<?> f : futures) { try { f.get(); } catch (Exception e) { log("[ERROR] " + e); } }
    }

    private void dfsContig(long minNode, long maxNode, String occtPath, String contigName) throws IOException {
        byte[] buf = new byte[1 << 18];
        BitSet visited = new BitSet((int) nodeCount + 1);
        BitSet inPath  = new BitSet((int) nodeCount + 1);
        long[] stack = new long[1 << 20]; int sp = 0;
        long blockCount = 0, bases = 0, edgeTotal = 0, records = 0;
        int nodesInBlock = 0; boolean blockStarted = false;
        OCCTWriter writer = new OCCTWriter(occtPath);
        long[] et = new long[256]; byte[] eo = new byte[256]; int[] ro = new int[256];
        long[] range = new long[2]; IntDList seq = IntDList.init();
        for (int cur = inDegZero.nextSetBit((int) minNode); cur > 0 && cur <= maxNode; cur = inDegZero.nextSetBit(cur + 1))
            if (!visited.get(cur)) { if (sp >= stack.length) stack = Arrays.copyOf(stack, stack.length * 2); stack[sp++] = pack(cur, true); }
        while (true) {
            if (sp == 0) {
                if (nodesInBlock > 0) { blockCount++; nodesInBlock = 0; blockStarted = false; }
                int next = visited.nextClearBit((int) minNode);
                if (next <= 0 || next > maxNode) break;
                if (sp >= stack.length) stack = Arrays.copyOf(stack, stack.length * 2);
                stack[sp++] = pack(next, true);
            }
            long packed = stack[--sp];
            int node = (int)(packed & 0x7FFFFFFFL);
            boolean isBlockStart = (packed & 0x80000000L) != 0;
            if (inPath.get(node) || visited.get(node)) continue;
            visited.set(node); inPath.set(node);
            bwt.index.rangeInto(node * 2, bwt.records.data.length, range);
            int rlen = (int)(range[1] - range[0]);
            if (rlen <= 0) { inPath.clear(node); continue; }
            if (buf.length < rlen) buf = new byte[Math.max(rlen, buf.length * 2)];
            copyFromBlocks(bwt.records.data, range[0], buf, rlen);
            int[] rp = {0}; long sigma = RunCodec.readByteCode(buf, rp);
            long prevTarget = 0;
            if (sigma > et.length) { et = Arrays.copyOf(et, (int)sigma*2); eo = Arrays.copyOf(eo, (int)sigma*2); }
            for (int e = 0; e < sigma; e++) {
                long to = RunCodec.readByteCode(buf, rp) + prevTarget;
                prevTarget = to; RunCodec.readByteCode(buf, rp);
                et[e] = to >> 1; eo[e] = (byte)((to & 1) == 0 ? 2 : 3);
            }
            short edgeCnt = (short)sigma;
            int rankOffLen = 0;
            if (sigma > 0) { if (ro.length < sigma + 2) ro = new int[(int)sigma + 2]; ro[0] = 0; int cum = 0, ri = 0; RunCodec codec = new RunCodec(sigma); while (rp[0] < rlen && ri < ro.length - 1) { long[] run = codec.read(buf, rp); cum += (int)run[1]; ri++; if (ri < ro.length) ro[ri] = cum; } rankOffLen = ri + 1; }
            int seqLen = 0; if (sequences != null && node >= 1 && node <= sequences.count) seqLen = sequences.length(node - 1);
            short flags = 0; if (!blockStarted || isBlockStart) { if (nodesInBlock > 0) { blockCount++; nodesInBlock = 0; } flags |= 1; blockStarted = true; }
            nodesInBlock++; boolean endsBlock = (edgeCnt != 1); if (!endsBlock) { int child = (int)et[0]; endsBlock = (child < minNode || child > maxNode || visited.get(child)); } if (endsBlock) flags |= 2;
            if (seqLen > 0) sequences.copyWords(node - 1, seq, 0);
            int[] rOff = rankOffLen > 0 ? Arrays.copyOf(ro, rankOffLen) : null;
            writer.writeRecord(node, seqLen, edgeCnt, flags, seq, et, eo, -1, -1, rOff);
            bases += seqLen; edgeTotal += edgeCnt; records++;
            if (records % 100000 == 0) log(String.format("  [%s] %,d nodes, %,d blocks, mem=%d MB", contigName, records, blockCount, usedMemMB()));
            if (edgeCnt == 1) { int child = (int)et[0]; if (child >= minNode && child <= maxNode && !visited.get(child)) { if (sp >= stack.length) stack = Arrays.copyOf(stack, stack.length * 2); stack[sp++] = pack(child, endsBlock); } }
            else if (edgeCnt > 1) { for (int e = (int)edgeCnt-1; e >= 0; e--) { int child = (int)et[e]; if (child >= minNode && child <= maxNode && !visited.get(child)) { if (sp >= stack.length) stack = Arrays.copyOf(stack, stack.length * 2); stack[sp++] = pack(child, true); } } }
            inPath.clear(node); if (endsBlock && nodesInBlock > 0) { blockCount++; nodesInBlock = 0; blockStarted = false; }
        }
        if (nodesInBlock > 0) blockCount++;
        writer.close();
        globalRecords.addAndGet(records); globalBlocks.addAndGet(blockCount); globalBases.addAndGet(bases); globalEdges.addAndGet(edgeTotal);
        log(String.format("  [%s] %,d nodes, %,d blocks", contigName, records, blockCount));
    }

    private static long pack(int node, boolean blockStart) { long v = node; if (blockStart) v |= 0x80000000L; return v; }
    private static void copyFromBlocks(sds.ByteBlocks d, long start, byte[] buf, int len) { int c = 0; while (c < len) { long p = start + c; int bi = (int)(p >>> sds.ByteBlocks.BLOCK_SHIFT); int off = (int)(p & sds.ByteBlocks.BLOCK_MASK); int n = Math.min(len - c, d.blocks[bi].length - off); System.arraycopy(d.blocks[bi], off, buf, c, n); c += n; } }
    private static String sanitizeFilename(String name) { return name.replaceAll("[^a-zA-Z0-9._-]", "_"); }
    private void printStats() { log(""); log("=== OCCT ==="); log(String.format("  nodes: %,d  records: %,d  blocks: %,d", nodeCount, globalRecords.get(), globalBlocks.get())); log(String.format("  bases: %,d  edges: %,d  mem: %d MB", globalBases.get(), globalEdges.get(), usedMemMB())); }
    private static long usedMemMB() { Runtime rt = Runtime.getRuntime(); return (rt.totalMemory() - rt.freeMemory()) >>> 20; }
    private static void log(String msg) { String ts = String.format("%tT", new java.util.Date()); System.out.println(ts + " " + msg); System.out.flush(); }

    public static void main(String[] args) throws IOException {
        String gbzPath = "/Users/wenjiepeng/Downloads/hprc-v1.1-mc-grch38.gbz";
        String outputDir = "/Users/wenjiepeng/Downloads/occt_output";
        String logPath = "/Users/wenjiepeng/Desktop/tmp/tmp.txt";
        if (args.length >= 1) gbzPath = args[0];
        if (args.length >= 2) outputDir = args[1];
        if (args.length >= 3) logPath = args[2];
        int threads = 0;
        if (args.length >= 4) threads = Integer.parseInt(args[3]);
        System.setOut(new PrintStream(new FileOutputStream(logPath)));
        System.setErr(System.out);
        long start = System.currentTimeMillis();
        try { new GBZ2OCCT(gbzPath, outputDir, threads).convert(); } catch (Exception e) { e.printStackTrace(); }
        log(String.format("%nTotal: %.1fs", (System.currentTimeMillis()-start)/1000.0));
    }
}
