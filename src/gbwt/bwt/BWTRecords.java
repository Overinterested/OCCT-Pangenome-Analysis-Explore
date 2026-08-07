package gbwt.bwt;

import sds.ByteBlocks;
import sds.RunCodec;
import sds.SdsCodec;
import sds.SdsReader;
import sds.SdsWriter;

import java.io.IOException;

/**
 * RecordArray 的 data 部分：Vec&lt;u8&gt;，所有节点的压缩记录首尾相连。
 * 大图上这块会超过 2 GiB（hprc-v2.1 实测 3.16 GiB），单个 byte[] 放不下，
 * 所以用分块的 {@link ByteBlocks} 存储。
 */
public class BWTRecords implements SdsCodec {
    public ByteBlocks data = ByteBlocks.of(new byte[0]);

    /**
     * 解出某条记录（配合 BWTIndex 给出的字节范围）。
     * 每条记录都要新建一个小 byte[] 和 NodeRecord；批量扫描请用
     * {@link #validateRecord} 或 {@link RunCodec#skipRecord}，零分配。
     */
    public RunCodec.NodeRecord decodeRecord(long start, long end) {
        byte[] slice = data.copyRange(start, end);
        return RunCodec.decodeRecord(slice, 0, slice.length);
    }

    /** 零分配校验：第 start..end 字节是否正好构成一条合法记录（消耗恰好到 end）。 */
    public boolean validateRecord(long start, long end) {
        byte[] slice = data.copyRange(start, end);
        return RunCodec.skipRecord(slice, 0, slice.length) == slice.length;
    }

    @Override
    public void decode(SdsReader r) throws IOException {
        long len = r.usize();
        data = r.rawBlocks(len);
        r.align8();
    }

    @Override
    public void encode(SdsWriter w) throws IOException {
        w.usize(data.length);
        for (byte[] block : data.blocks) {
            w.rawBytes(block);
        }
        w.align8();
    }
}
