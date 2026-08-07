package gbwt;

import sds.SdsCodec;
import sds.SdsReader;
import sds.SdsWriter;
import sds.SdsPrimitives.IntVector;
import sds.SdsPrimitives.RawBitVector;
import sds.SdsPrimitives.SparseVector;

import java.io.IOException;

/**
 * 论文原话："serialized as an optional structure in an unspecified format" ——
 * 官方没有把它标准化，具体布局取决于写文件的实现版本。下面这套是当前 C++ 实现
 * （gbwt::DASamples）实际写出来的布局。如果你只关心"序列/拓扑/样本"，
 * 可以完全不用管这个类的内部字段，present=false 时它就是空的。
 *
 * 这是 Option&lt;DASamples&gt;：decode/encode 自己处理开头那个"多少个 element"的
 * 长度前缀，0 表示不存在。
 */
public class GBWTDASamples implements SdsCodec {
    public boolean present = false;

    public RawBitVector sampledRecords;
    public SparseVector bwtRanges;
    public SparseVector sampledOffsets;
    public IntVector array;

    @Override
    public void decode(SdsReader r) throws IOException {
        long elements = r.usize();
        if (elements == 0) {
            present = false;
            return;
        }
        if (elements < 0 || elements > r.bytesRemaining() / 8) {
            throw new IllegalStateException("DASamples: element 计数 " + elements
                    + " 不合法（剩余可消费字节 " + r.bytesRemaining() + "）——文件损坏或上游解析错位");
        }
        present = true;
        sampledRecords = RawBitVector.decode(r);
        bwtRanges = SparseVector.decode(r);
        sampledOffsets = SparseVector.decode(r);
        array = IntVector.decode(r);
    }

    @Override
    public void encode(SdsWriter w) throws IOException {
        if (!present) {
            w.usize(0);
            return;
        }
        SdsWriter sub = new SdsWriter();
        sampledRecords.encode(sub);
        bwtRanges.encode(sub);
        sampledOffsets.encode(sub);
        array.encode(sub);
        byte[] bytes = sub.toByteArray();
        w.usize(bytes.length / 8L);
        w.rawBytes(bytes);
    }
}
