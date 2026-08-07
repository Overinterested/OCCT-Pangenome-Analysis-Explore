package gbwt;

import sds.SdsCodec;
import sds.SdsReader;
import sds.SdsWriter;
import sds.StringCodec;

import java.io.IOException;
import java.util.LinkedHashMap;
import java.util.Map;

/** 论文 2.3.2 节 "Tags"：key/value 交替存成一个 StringArray。key 按惯例全小写、不重复。 */
public class GBWTTags implements SdsCodec {
    public Map<String, String> entries = new LinkedHashMap<>();

    @Override
    public void decode(SdsReader r) throws IOException {
        String[] flat = StringCodec.decodeStringArray(r);
        entries = new LinkedHashMap<>();
        for (int i = 0; i + 1 < flat.length; i += 2) {
            entries.put(flat[i], flat[i + 1]);
        }
    }

    @Override
    public void encode(SdsWriter w) throws IOException {
        String[] flat = new String[entries.size() * 2];
        int i = 0;
        for (Map.Entry<String, String> e : entries.entrySet()) {
            flat[i++] = e.getKey().toLowerCase();
            flat[i++] = e.getValue();
        }
        StringCodec.encodeStringArray(w, flat);
    }
}
