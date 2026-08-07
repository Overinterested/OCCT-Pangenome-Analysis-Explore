import gbwt.GBWT;
import gbwt.GBWTDASamples;
import gbwt.GBWTHeader;
import gbwt.GBWTTags;
import gbwt.bwt.BWT;
import gbwt.meta.GBWTMeta;
import gbwtgraph.GBWTGraphHeader;
import gbwtgraph.GBWTGraphSequences;
import gbwtgraph.PackedSequences;
import gbwtgraph.translation.Translation;
import sds.SdsReader;
import sds.SdsSkip;
import sds.SdsStreamReader;

import java.io.Closeable;
import java.io.IOException;

/**
 * GBZ 文件的分块流式读取器。
 * <p>
 * GBZ version 2 (SERIALIZATION.md)：6 个 stage:
 * GBZHeader -> GBZTags -> GBWT(header+tags+BWT) -> GBWTGraphHeader -> Sequences -> Translation
 * <p>
 * GBZ version 1 兼容：readGBWT() 后额外调用 readDASamples()/readMetadata()。
 */
public class GBZStreamer implements Closeable {
    private static final String[] STAGE_NAMES = {
            "GBZHeader", "GBZTags", "GBWT", "GBWTGraphHeader", "Sequences", "Translation", "EOF"
    };

    private int stage;
    final SdsStreamReader r;
    private int gbzVersion;
    private GBWTHeader gbwtHeader;
    GBWTGraphHeader graphHeader;

    public GBZStreamer(String path) throws IOException {
        this.r = SdsReader.openStream(path);
    }

    public long position() {
        return r.position();
    }

    public String nextStage() {
        return STAGE_NAMES[stage];
    }

    public boolean hasRemaining() {
        return r.hasRemaining();
    }

    public int gbzVersion() {
        return gbzVersion;
    }

    private void need(int expected) {
        if (stage != expected) {
            throw new IllegalStateException("GBZStreamer: 顺序错误，下一个模块是 " + STAGE_NAMES[stage]
                    + "，现在不能处理 " + STAGE_NAMES[expected]);
        }
    }

    // ---------------- 1. GBZ header ----------------

    public GBZHeader readGBZHeader() throws IOException {
        need(0);
        GBZHeader h = new GBZHeader();
        h.decode(r);
        gbzVersion = h.version;
        stage = 1;
        return h;
    }

    // ---------------- 2. GBZ tags ----------------

    public GBZTags readGBZTags() throws IOException {
        need(1);
        GBZTags t = new GBZTags();
        t.decode(r);
        stage = 2;
        return t;
    }

    public void skipGBZTags() {
        need(1);
        SdsSkip.stringArray(r);
        stage = 2;
    }

    // ---------------- 3. GBWT ----------------

    /**
     * 读取 GBWT 核心部分：header -> tags -> BWT。
     * version 2 文件不包含 DA/meta，parse 直接结束；
     * 但如果后续数据看起来像合法的 DA/meta（usize 为 0 或合理值），
     * 也一并消费，避免干扰 GBWTGraph 解析。
     */
    public GBWT readGBWT() throws IOException {
        need(2);
        GBWT gbwt = new GBWT();
        gbwt.header = new GBWTHeader();
        gbwt.header.decode(r);
        if (!gbwt.header.isSimpleSds()) {
            throw new IllegalStateException("这是旧版 SDSL 格式（非 simple_sds），本项目只覆盖现代 simple_sds 格式。");
        }
        gbwt.tags = new GBWTTags();
        gbwt.tags.decode(r);
        gbwtHeader = gbwt.header;
        gbwt.bwt = new BWT();
        gbwt.bwt.decode(r);

        // 尝试消费 DA/meta（version 1 文件会包含；version 2 文件不含则回退）
        gbwt.samples = tryReadDASamples();
        gbwt.meta = tryReadMetadata();

        // gbwtHeader already set above
        stage = 3;
        return gbwt;
    }

    /**
     * 尝试读取 DASamples：读 usize，若为 0 则 absent；若为合理正数则完整解析；
     * 若为垃圾值则 seek 回退，标记 absent。
     */
    private GBWTDASamples tryReadDASamples() throws IOException {
        long savedPos = r.position();
        long elements = r.usize();
        if (elements == 0) {
            return new GBWTDASamples(); // absent
        }
        if (elements > 0 && elements <= r.bytesRemaining() / 8) {
            r.seekToPosition(savedPos);
            GBWTDASamples s = new GBWTDASamples();
            s.decode(r);
            return s;
        }
        // 垃圾值：不是 DASamples，回退
        r.seekToPosition(savedPos);
        return new GBWTDASamples();
    }

    /**
     * 尝试读取 Metadata：读 usize，若为 0 则 absent；
     * 若为合理正数则完整解析（含 MetaHeader tag 校验）；
     * 若为垃圾值则 seek 回退，标记 absent。
     */
    private GBWTMeta tryReadMetadata() throws IOException {
        long savedPos = r.position();
        long elements = r.usize();
        if (elements == 0) {
            GBWTMeta m = new GBWTMeta();
            if (gbwtHeader.hasMetadata()) {
                throw new IllegalStateException("GBWT header says metadata present but data is absent");
            }
            return m;
        }
        if (elements > 0 && elements <= r.bytesRemaining() / 8) {
            r.seekToPosition(savedPos);
            GBWTMeta m = new GBWTMeta();
            m.decode(r);
            if (m.present != gbwtHeader.hasMetadata()) {
                throw new IllegalStateException("header.hasMetadata() 和实际读到的 metadata 是否存在不一致");
            }
            return m;
        }
        // 垃圾值：不是 Metadata，回退
        r.seekToPosition(savedPos);
        GBWTMeta m = new GBWTMeta();
        m.present = false;
        return m;
    }

    // ---------------- 3b. 显式 DA/meta（version 1 兼容） ----------------

    public GBWTDASamples readDASamples() throws IOException {
        GBWTDASamples s = new GBWTDASamples();
        s.decode(r);
        return s;
    }

    public void skipDASamples() {
        skipBlockPayload("DASamples");
    }

    public GBWTMeta readMetadata() throws IOException {
        GBWTMeta m = new GBWTMeta();
        m.decode(r);
        if (m.present != gbwtHeader.hasMetadata()) {
            throw new IllegalStateException("header.hasMetadata() 和实际读到的 metadata 是否存在不一致");
        }
        return m;
    }

    public void skipMetadata() {
        skipBlockPayload("Metadata");
    }

    private void skipBlockPayload(String what) {
        long elements = r.usize();
        if (elements < 0 || elements > r.bytesRemaining() / 8) {
            throw new IllegalStateException(what + ": element 计数 " + elements
                    + " 不合法（剩余可消费字节 " + r.bytesRemaining() + "）");
        }
        r.skipBytes(elements * 8);
    }

    // ---------------- 4. GBWTGraph header ----------------

    /**
     * 读取 GBWTGraph header。
     * 规范文件：正常解析 24 字节。
     * 预规范文件（无 header）：回退，用 GBWT header 推导节点数。
     */
    public GBWTGraphHeader readGraphHeader() throws IOException {
        need(3);
        long savedPos = r.position();
        long w0 = r.u64();
        int tag = (int) (w0 & 0xFFFFFFFFL);

        if (tag == GBWTGraphHeader.TAG) {
            int version = (int) (w0 >>> 32);
            long nodes = r.u64();
            long flags = r.u64();
            graphHeader = new GBWTGraphHeader();
            graphHeader.tag = tag;
            graphHeader.version = version;
            graphHeader.nodes = nodes;
            graphHeader.flags = flags;
            stage = 4;
            return graphHeader;
        }

        // 预规范文件：回退，构造合成 header（version 3 = StringArray）
        r.seekToPosition(savedPos);
        graphHeader = new GBWTGraphHeader();
        graphHeader.tag = GBWTGraphHeader.TAG;
        graphHeader.version = 3;
        graphHeader.nodes = gbwtHeader.alphabetSize - gbwtHeader.offset;
        graphHeader.flags = GBWTGraphHeader.FLAG_SIMPLE_SDS;
        stage = 4;
        return graphHeader;
    }

    // ---------------- 5. 节点序列 ----------------

    public GBWTGraphSequences readSequences() throws IOException {
        need(4);
        GBWTGraphSequences seq = new GBWTGraphSequences();
        seq.zstd = graphHeader.usesZstdSequences();
        seq.decode(r);
        stage = 5;
        return seq;
    }

    public GBWTGraphSequences readSequencesPacked() throws IOException {
        need(4);
        if (graphHeader.usesZstdSequences()) {
            throw new IllegalStateException("zstd 压缩序列暂不支持打包读取，请用 readSequences()/skipSequences()");
        }
        GBWTGraphSequences seq = new GBWTGraphSequences();
        seq.packed = PackedSequences.fromStringArray(r);
        stage = 5;
        return seq;
    }

    public void skipSequences() {
        need(4);
        if (graphHeader.usesZstdSequences()) {
            SdsSkip.sparseVector(r);
            r.usize();
            SdsSkip.byteVector(r);
        } else {
            SdsSkip.stringArray(r);
        }
        stage = 5;
    }

    // ---------------- 6. segment 转换 ----------------

    public Translation readTranslation() throws IOException {
        need(5);
        Translation t = new Translation();
        t.decode(r);
        stage = 6;
        return t;
    }

    public void skipTranslation() {
        need(5);
        SdsSkip.stringArray(r);
        SdsSkip.sparseVector(r);
        stage = 6;
    }

    @Override
    public void close() throws IOException {
        r.close();
    }
}
