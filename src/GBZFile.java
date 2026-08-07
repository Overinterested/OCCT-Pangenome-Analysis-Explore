import edu.sysu.pmglab.io.reader.ChannelReaderStream;
import edu.sysu.pmglab.io.reader.ReaderStream;
import gbwt.GBWT;
import gbwtgraph.GBWTGraph;
import sds.SdsCodec;
import sds.SdsReader;
import sds.SdsWriter;

import java.io.IOException;

/**
 * 论文 2.3.5 节 "GBZ"。
 * GBZ version 2：GBWT 区域只含 header + tags + BWT。
 * GBZ version 1：GBWT 区域之后还有 DA/meta（readGBWT 内置 try-rewind 自动处理）。
 */
public class GBZFile implements SdsCodec {
    public GBZHeader header = new GBZHeader();
    public GBZTags tags = new GBZTags();
    public GBWT gbwt = new GBWT();
    public GBWTGraph graph = new GBWTGraph();

    public static GBZFile parse(String path) throws IOException {
        return parse(path, GBZFilter.LOAD_ALL);
    }

    public static GBZFile parse(String path, GBZFilter filter) throws IOException {
        GBZFile file = new GBZFile();
        try (GBZStreamer reader = new GBZStreamer(path)) {
            file.header = reader.readGBZHeader();
            if (filter.gbzTags) file.tags = reader.readGBZTags(); else reader.skipGBZTags();

            // readGBWT 内置 try-rewind 自动处理 version 1 的 DA/meta
            file.gbwt = reader.readGBWT();

            file.graph.header = reader.readGraphHeader();
            if (filter.graphSequences) file.graph.sequences = reader.readSequences();
            else reader.skipSequences();
            if (filter.graphTranslation) file.graph.translation = reader.readTranslation();
            else reader.skipTranslation();
        }
        return file;
    }

    public void parseTo(String path) throws IOException {
        SdsWriter writer = new SdsWriter();
        this.encode(writer);
        writer.writeToFile(path);
    }

    @Override
    public void decode(SdsReader r) throws IOException {
        header = new GBZHeader();
        header.decode(r);
        tags = new GBZTags();
        tags.decode(r);
        gbwt = new GBWT();
        gbwt.decode(r);
        graph = new GBWTGraph();
        graph.decode(r);
    }

    @Override
    public void encode(SdsWriter writer) throws IOException {
        header.encode(writer);
        tags.encode(writer);
        gbwt.encode(writer);
        graph.encode(writer);
    }
}
