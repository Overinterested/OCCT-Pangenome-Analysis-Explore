package gbwt.bwt;

import sds.SdsCodec;
import sds.SdsReader;
import sds.SdsWriter;
import sds.SdsPrimitives.SparseVector;

import java.io.IOException;

/** RecordArray 的 index 部分：sparse bitvector，universe=数据总字节数，值=每条记录的起始偏移。 */
public class BWTIndex implements SdsCodec {
    public SparseVector recordStarts;

    public long recordCount() { return recordStarts.size(); }

    /** 第 i 条记录在 data 里的 [start, end) 字节范围，需要传入 data 总长度来确定最后一条记录的结尾。 */
    public long[] range(int i, long dataLength) {
        long start = recordStarts.values[i];
        long end = (i + 1 < recordStarts.values.length) ? recordStarts.values[i + 1] : dataLength;
        return new long[]{start, end};
    }

    /** 零分配版本：把范围写进调用方复用的 out[2]（批量遍历记录时用）。 */
    public void rangeInto(int i, long dataLength, long[] out) {
        out[0] = recordStarts.values[i];
        out[1] = (i + 1 < recordStarts.values.length) ? recordStarts.values[i + 1] : dataLength;
    }

    @Override
    public void decode(SdsReader r) throws IOException {
        recordStarts = SparseVector.decode(r);
    }

    @Override
    public void encode(SdsWriter w) throws IOException {
        recordStarts.encode(w);
    }
}
