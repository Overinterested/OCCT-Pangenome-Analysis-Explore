import gbwt.GBWT;
import gbwt.GBWTHeader;
import gbwt.bwt.BWT;
import gbwt.bwt.BWTIndex;
import gbwt.bwt.BWTRecords;
import gbwt.meta.GBWTMeta;
import gbwt.meta.MetaPaths;
import gbwtgraph.GBWTGraphSequences;
import gbwtgraph.PackedSequences;
import gbwtgraph.translation.Mapping;
import gbwtgraph.translation.Segments;
import gbwtgraph.translation.Translation;
import sds.DNACodec;
import sds.RunCodec;

import java.io.IOException;

/**
 * GBZ 核心类使用示例，以 hprc-v2.1-mc-grch38.gbz 为例。
 * 采用 GBZStreamer 分模块流式读取。
 *
 * 运行方式：java -Xmx16g Usage
 *
 * @author Wenjie Peng
 */
public class Usage {
    static String path = "/Users/wenjiepeng/Downloads/hprc-v2.1-mc-grch38.gbz";

    public static void main(String[] args) throws IOException {
        if (args.length > 0) path = args[0];

        try (GBZStreamer s = new GBZStreamer(path)) {
            GBZHeader gbzHdr = s.readGBZHeader();
            s.readGBZTags();
            System.out.printf("GBZ version=%d%n", gbzHdr.version);

            // --- GBWT（readGBWT 内置 try-rewind 自动处理 DA/meta） ---
            GBWT gbwt = s.readGBWT();
            GBWTHeader hdr = gbwt.header;

            System.out.printf("节点总数: %,d%n", hdr.alphabetSize - hdr.offset);
            System.out.printf("原始输入序列: %,d 条%n", hdr.sequences);

            // ============================================================
            // 1. BWT / BWTIndex / BWTRecords —— 图的拓扑
            // ============================================================
            System.out.println("\n=== 1. BWT 压缩索引（图的拓扑） ===");

            BWT bwt = gbwt.bwt;
            BWTIndex index = bwt.index;
            BWTRecords records = bwt.records;
            System.out.printf("节点记录数: %,d%n", index.recordCount());
            System.out.printf("压缩数据: %,d 字节 (%.2f GB)%n",
                    records.data.length, records.data.length / (1024.0 * 1024 * 1024));

            // 每条 BWT 记录 = 出边表(delta-编码 node id + rank offset) + 游程体(outrank+len)
            System.out.println("\n--- node 100 ---");
            RunCodec.NodeRecord rec = bwt.getRecord(100);
            System.out.printf("  出边数 σ=%d  游程数=%d%n", rec.outgoing.size(), rec.runs.size());
            if (!rec.outgoing.isEmpty()) {
                long[] e = rec.outgoing.get(0);
                System.out.printf("  第1边: →node %,d rank_offset=%d%n", e[0], e[1]);
            }
            if (!rec.runs.isEmpty()) {
                long[] r = rec.runs.get(0);
                System.out.printf("  第1游程: outrank=%d len=%d%n", r[0], r[1]);
            }

            // 前 20 个节点出边数（node 0 = sentinel，σ 很大因为所有路径终点都汇集于此）
            System.out.println("\n--- 前 20 个节点出边数 ---");
            for (int i = 0; i < 20; i++) {
                RunCodec.NodeRecord nr = bwt.getRecord(i);
                System.out.printf("  node %d: σ=%d%s%n", i, nr.outgoing.size(),
                        i == 0 ? "  ← sentinel（所有路径终点）" : "");
            }

            // ============================================================
            // 2. Metadata
            // ============================================================
            System.out.println("\n=== 2. Metadata（样本/contig/路径） ===");

            GBWTMeta meta = gbwt.meta;
            if (meta.present) {
                System.out.printf("样本数: %d  contig数: %d  路径数: %d%n",
                        meta.metaSamples.names.length,
                        meta.metaCotigs.names.length,
                        meta.metaPaths.paths.size());

                System.out.println("\n--- 前 5 个样本 ---");
                for (int i = 0; i < Math.min(5, meta.metaSamples.names.length); i++)
                    System.out.printf("  [%d] %s%n", i, meta.metaSamples.get(i));

                System.out.println("\n--- 前 5 条路径 ---");
                for (int i = 0; i < Math.min(5, meta.metaPaths.paths.size()); i++) {
                    MetaPaths.PathName p = meta.metaPaths.paths.get(i);
                    String sn = p.sample < meta.metaSamples.names.length
                            ? meta.metaSamples.get(p.sample) : "?";
                    String cn = p.contig < meta.metaCotigs.names.length
                            ? meta.metaCotigs.get(p.contig) : "?";
                    System.out.printf("  [%d] sample=%s contig=%s phase=%d frag=%d%n",
                            i, sn, cn, p.phase, p.count);
                }

                // 统计每个样本的路径数
                System.out.println("\n--- 前 5 个样本的路径数 ---");
                long[] pc = new long[(int) meta.metaSamples.names.length];
                for (MetaPaths.PathName pn : meta.metaPaths.paths) {
                    if (pn.sample < pc.length) pc[(int) pn.sample]++;
                }
                for (int i = 0, shown = 0; i < pc.length && shown < 5; i++) {
                    if (pc[i] > 0) {
                        System.out.printf("  %s: %,d 条路径%n", meta.metaSamples.get(i), pc[i]);
                        shown++;
                    }
                }
            }
            gbwt = null;
            System.gc();

            // ============================================================
            // 3. GBWTGraph header
            // ============================================================
            s.readGraphHeader();
            long nodeCount = s.graphHeader.nodes;
            System.out.printf("%n--- GBWTGraph: %,d 个存在节点 (graph.nodes) ---%n", nodeCount);
            System.out.println("说明：BWT 有 425M 个节点（含空节点），但 " +
                    "graph.nodes 只统计被至少一条路径穿过的节点（即局部字母表 σ>0 的节点）");

            // ============================================================
            // 4. PackedSequences / GBWTGraphSequences —— 节点序列
            // ============================================================
            System.out.println("\n=== 3. PackedSequences（节点 DNA 序列） ===");

            GBWTGraphSequences seq = s.readSequencesPacked();
            PackedSequences packed = seq.packed;
            System.out.printf("序列数: %,d  format=%s  总碱基: %,d%n",
                    packed.count, packed.format, packed.totalBases);

            System.out.println("\n--- 前 10 条 ---");
            for (int i = 0; i < 10; i++) {
                String dna = packed.decode(i);
                System.out.printf("  node %d (len=%d): %s%n", i, packed.length(i),
                        dna.length() > 60 ? dna.substring(0, 57) + "..." : dna);
            }

            // 零分配单碱基查询
            System.out.println("\n--- DNACodec 零分配查询（node 5 前 5 bp）---");
            int[] buf = new int[DNACodec.wordsFor(packed.length(5), packed.format)];
            packed.copyWords(5, buf, 0);
            for (int p = 0; p < Math.min(5, packed.length(5)); p++)
                System.out.printf("  pos %d: %c%n", p, DNACodec.charAt(buf, p, packed.format));

            // 正向/反向互补
            String fwd = packed.decode(0);
            String rev = GBWTGraphSequences.reverseComplement(fwd);
            System.out.printf("%n  正向: %s%n", fwd.length() > 50 ? fwd.substring(0, 47) + "..." : fwd);
            System.out.printf("  反向: %s%n", rev.length() > 50 ? rev.substring(0, 47) + "..." : rev);

            seq = null;
            System.gc();

            // ============================================================
            // 5. Translation: Segments / Mapping
            // ============================================================
            System.out.println("\n=== 4. Translation（node → segment 映射） ===");

            Translation trans = s.readTranslation();
            Segments segs = trans.segments;
            Mapping map = trans.mapping;

            System.out.printf("segment 数: %,d  mapping 边界: %,d%n",
                    segs.names.length, map.nodeToSegment.size());

            if (segs.names.length > 0) {
                System.out.println("\n--- 前 5 个 segment ---");
                for (int i = 0; i < Math.min(5, segs.names.length); i++)
                    System.out.printf("  [%d] %s%n", i, segs.names[i]);

                System.out.println("\n--- segment → 节点区间 ---");
                long[] sv = map.nodeToSegment.values;
                for (int i = 0; i < Math.min(5, sv.length); i++) {
                    long st = sv[i], ed = (i + 1 < sv.length) ? sv[i + 1] : nodeCount + 1;
                    System.out.printf("  %s: nodes [%,d, %,d) 长度 %,d%n",
                            segs.names[i], st, ed, ed - st);
                }

                System.out.println("\n--- node → segment 查询 ---");
                for (long nid : new long[]{1, 100, 1_000_000, nodeCount}) {
                    int si = map.segmentIndexForNode(nid);
                    System.out.printf("  node %,d → seg[%d] %s%n", nid, si,
                            si < segs.names.length ? segs.names[si] : "?");
                }
            } else {
                System.out.println("此文件没有 segment 映射（translation 为空）");
            }
        }
    }
}
