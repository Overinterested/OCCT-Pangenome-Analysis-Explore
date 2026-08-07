import edu.sysu.pmglab.bytecode.ByteStream;
import edu.sysu.pmglab.container.list.IntDList;
import edu.sysu.pmglab.ecc.ECCWriter;
import edu.sysu.pmglab.ecc.compressor.zstd.ZSTDCompressor;
import edu.sysu.pmglab.ecc.compressor.zstd.ZSTDDecompressor;
import edu.sysu.pmglab.ecc.record.BoxRecord;
import edu.sysu.pmglab.ecc.type.FieldType;
import edu.sysu.pmglab.io.reader.ChannelReaderStream;
import edu.sysu.pmglab.io.reader.ReaderStream;

import java.io.IOException;

/**
 * @author Wenjie Peng
 * @create 2026-07-31 17:09
 * @description
 */
public class Test {
    public static void main(String[] args) throws IOException {
        GBZFile file = GBZFile.parse("/Users/wenjiepeng/Downloads/y.giraffe.gbz");
        file.parseTo("/Users/wenjiepeng/Downloads/y.giraffe_1.gbz");

        // 一套自定义的文件按行读取的操作示例
        // 注意两点：1) 构造参数换成真实路径；2) ByteStream 要指定容量——
        // 无参构造的内部数组长度为 0，reader.read() 会一直读到 0 个字节，循环永不结束。
        ByteStream cache = new ByteStream(1 << 16);
        ReaderStream reader = new ReaderStream(new ChannelReaderStream("/Users/wenjiepeng/Downloads/y.giraffe.gbz"));
        while(reader.read(cache.bytes())!=-1){
            byte[] item = cache.bytes();
            // 操作
            cache.clear();
        }
        cache.close();
        reader.close();
    }
}
