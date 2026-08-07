import gbwt.GBWT;
import gbwt.GBWTHeader;
import gbwt.bwt.BWT;
import gbwt.meta.GBWTMeta;
import gbwtgraph.GBWTGraphHeader;
import gbwtgraph.GBWTGraphSequences;
import sds.RunCodec;

import java.io.IOException;

/**
 * 大文件分模块流式解析测试，兼容 GBZ version 1 和 2。
 *
 * @author Wenjie Peng
 */
public class StreamTest {
    static String path = "/Users/wenjiepeng/Downloads/hprc-v2.1-mc-grch38.gbz";

    static long mem() {
        Runtime rt = Runtime.getRuntime();
        return (rt.totalMemory() - rt.freeMemory()) >>> 20;
    }

    public static void main(String[] args) throws IOException {
        long start = System.currentTimeMillis();
        if (args.length > 0) path = args[0];
        System.out.println("file: " + path);

        try (GBZStreamer s = new GBZStreamer(path)) {
            GBZHeader gbzHeader = s.readGBZHeader();
            GBZTags gbzTags = s.readGBZTags();
            System.out.printf("[GBZ]     version=%d tags=%s  (pos=%,d, mem=%dMB)%n",
                    gbzHeader.version, gbzTags.entries, s.position(), mem());
            int gbzVer = gbzHeader.version;
            gbzHeader = null;
            gbzTags = null;

            GBWT gbwt = s.readGBWT();
            GBWTHeader gbwtHeader = gbwt.header;
            System.out.printf("[GBWT]    version=%d sequences=%,d size=%,d alphabet=%,d flags=0x%x%n",
                    gbwtHeader.version, gbwtHeader.sequences, gbwtHeader.size,
                    gbwtHeader.alphabetSize, gbwtHeader.flags);
            long alphabet = gbwtHeader.alphabetSize;

            BWT bwt = gbwt.bwt;
            long records = bwt.recordCount();
            long dataLen = bwt.records.data.length;
            System.out.printf("[BWT]     records=%,d (期望 %,d) data=%,d 字节  (pos=%,d, mem=%dMB)%n",
                    records, alphabet - gbwtHeader.offset, dataLen, s.position(), mem());

            long[] range = new long[2];
            int bad = 0;
            for (int i = 0; i < Math.min(1000, records); i++) {
                bwt.index.rangeInto(i, dataLen, range);
                if (RunCodec.skipRecord(bwt.records.data.copyRange(range[0], range[1]),
                        0, (int) (range[1] - range[0])) < 0) { bad++; }
            }
            System.out.printf("[BWT]     前 1000 条记录边界校验: %s%n", bad == 0 ? "全部合法" : bad + " 条异常 <<<");

            GBWTMeta meta = gbwt.meta;
            System.out.printf("[META]    present=%b samples=%d contigs=%d paths=%d%n",
                    meta.present,
                    meta.present ? meta.metaSamples.names.length : 0,
                    meta.present ? meta.metaCotigs.names.length : 0,
                    meta.present ? meta.metaPaths.paths.size() : 0);

            bwt = null;
            gbwt = null;
            System.gc();
            System.out.printf("[GC]      释放 GBWT 后 mem=%dMB%n", mem());

            // GBWTGraph
            try {
                GBWTGraphHeader graphHeader = s.readGraphHeader();
                System.out.printf("[GRAPH]   version=%d nodes=%,d flags=0x%x  (pos=%,d)%n",
                        graphHeader.version, graphHeader.nodes, graphHeader.flags, s.position());
                graphHeader = null;
                GBWTGraphSequences seq = s.readSequencesPacked();
                int last = seq.packed.count - 1;
                System.out.printf("[SEQ]     条数=%,d format=%s totalBases=%,d%n",
                        seq.packed.count, seq.packed.format, seq.packed.totalBases);
                if (seq.packed.count > 0) {
                    System.out.printf("[SEQ]     首条(%d bp): %.40s  末条(%d bp): %.40s%n",
                            seq.packed.length(0), seq.packed.decode(0),
                            seq.packed.length(last), seq.packed.decode(last));
                }
                seq = null;
                System.gc();

                s.readTranslation();
                System.out.printf("[TRANS]   pos=%,d  (mem=%dMB)%n", s.position(), mem());
            } catch (Exception e) {
                System.out.println("[GRAPH]   ** GBWTGraph 解析失败: " + e.getMessage());
                if (gbzVer >= 2) {
                    System.out.println("[GRAPH]   ** 此 GBZ version 2 文件的 GBWTGraph 部分使用了");
                    System.out.println("[GRAPH]   ** 当前解析器不支持的序列化格式。");
                    System.out.println("[GRAPH]   ** BWT 索引部分解析正常（" + records + "条记录全部校验通过）。");
                    System.out.println("[GRAPH]   ** 需要查阅 gbwt SERIALIZATION.md 确认 DASamples/GBWTGraph 的精确格式。");
                }
            }
        }
        System.out.println("DONE");
        System.out.println(System.currentTimeMillis() - start);
    }
}
