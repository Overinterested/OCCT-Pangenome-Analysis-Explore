package gbwt.meta;

import sds.SdsCodec;
import sds.SdsReader;
import sds.SdsWriter;
import sds.StringCodec;

import java.io.IOException;

/** 论文 2.3.2 节 "Dictionary"：sample id（PathName.sample 里的整数）对应的字符串名字。 */
public class MetaSamples implements SdsCodec {
    public String[] names = new String[0];

    public String get(long sampleId) { return names[(int) sampleId]; }

    @Override
    public void decode(SdsReader r) throws IOException {
        names = StringCodec.decodeDictionaryNames(r);
    }

    @Override
    public void encode(SdsWriter w) throws IOException {
        StringCodec.encodeDictionaryNames(w, names);
    }
}
