package gbwtgraph.translation;

import edu.sysu.pmglab.container.indexable.DynamicIndexableMap;
import edu.sysu.pmglab.container.indexable.IndexableSet;
import sds.SdsCodec;
import sds.SdsReader;
import sds.SdsWriter;
import sds.StringCodec;

import java.io.IOException;

/** 论文 2.3.4 节："a string array storing segment names S0,...,Sm-1"。没有 translation 时这里是空数组。 */
public class Segments implements SdsCodec {
    public String[] names = new String[0];

    @Override
    public void decode(SdsReader r) throws IOException {
        names = StringCodec.decodeStringArray(r);
    }

    @Override
    public void encode(SdsWriter w) throws IOException {
        StringCodec.encodeStringArray(w, names);
    }
}
