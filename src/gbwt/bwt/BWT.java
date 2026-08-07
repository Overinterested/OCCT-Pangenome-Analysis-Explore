package gbwt.bwt;

import sds.ByteBlocks;
import sds.RunCodec;
import sds.SdsCodec;
import sds.SdsReader;
import sds.SdsWriter;
import sds.SdsPrimitives.SparseVector;

import java.io.IOException;
import java.util.List;

/**
 * 论文 2.1.6 节 "GBWT index" + 2.3.3 节：完整的压缩 BWT（对应 gbwt::RecordArray）。
 * 布局顺序是先 index 后 data（BWTIndex 需要先知道每条记录从哪开始，
 * 读的时候才能对 BWTRecords 里的大字节块做切片）。
 */
public class BWT implements SdsCodec {
    public BWTIndex index;
    public BWTRecords records;

    public long recordCount() { return index.recordCount(); }

    /** 取第 i 个节点的记录（0-based，对应内部 node id 从某个基准开始的顺序）。 */
    public RunCodec.NodeRecord getRecord(int i) {
        long[] range = index.range(i, records.data.length);
        return records.decodeRecord(range[0], range[1]);
    }

    /** 从一组"人类可读"的节点记录构建 BWT，方便把你自己的图结构转换成 GBZ 需要的格式。 */
    public static BWT fromNodeRecords(List<RunCodec.NodeRecord> nodes) {
        SdsWriter dataWriter = new SdsWriter();
        long[] starts = new long[nodes.size()];
        for (int i = 0; i < nodes.size(); i++) {
            starts[i] = dataWriter.size();
            RunCodec.encodeRecord(dataWriter, nodes.get(i));
        }
        byte[] data = dataWriter.toByteArray();

        BWT bwt = new BWT();
        bwt.records = new BWTRecords();
        bwt.records.data = ByteBlocks.of(data);
        bwt.index = new BWTIndex();
        bwt.index.recordStarts = SparseVector.build(starts);
        return bwt;
    }

    @Override
    public void decode(SdsReader r) throws IOException {
        index = new BWTIndex();
        index.decode(r);
        records = new BWTRecords();
        records.decode(r);
    }

    @Override
    public void encode(SdsWriter w) throws IOException {
        index.encode(w);
        records.encode(w);
    }
}
