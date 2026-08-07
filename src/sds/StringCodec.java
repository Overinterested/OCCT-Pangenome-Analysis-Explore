package sds;

// TODO: 换成你项目里 ZSTDCompressor / ZSTDDecompressor 实际所在的包，例如：
// import your.pkg.zstd.ZSTDCompressor;
// import your.pkg.zstd.ZSTDDecompressor;
// 下面假设它们各有一个方法：
//   byte[] compress(byte[] input)
//   byte[] decompress(byte[] input, int originalSize)
// 如果你项目里的方法名/签名不一样，把 encodeEvenCompressed / decodeEvenCompressed 里对应两行改掉就行。

import com.github.luben.zstd.ZstdDecompressCtx;
import edu.sysu.pmglab.ecc.compressor.zstd.ZSTDCompressor;
import edu.sysu.pmglab.ecc.compressor.zstd.ZSTDDecompressor;
import sds.SdsPrimitives.IntVector;
import sds.SdsPrimitives.SparseVector;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * 论文 2.2 节 "Improved sparse bitvectors" + 2.3.2 节 "String array"。
 * <p>
 * 普通版本：字符串拼成一条大字符串，sparse bitvector 记录每个字符串起点，
 * 字符本身做字母表压缩（比如 DNA 只有 ACGTN，重映射成 0..4 再按最小位宽打包）。
 * even 压缩版本：GBWTGraph 存节点序列专用，只存正向（偶数下标）字符串，
 * 不做字母表压缩，改用 zstd 压缩整块字节流（论文 2.3.4 节）。
 */
public class StringCodec {

    // ---------------- 普通 StringArray（字母表压缩） ----------------

    public static String[] decodeStringArray(SdsReader r) {
        SparseVector index = SparseVector.decode(r);

        long alphabetSize = r.usize();
        byte[] compToChar = r.rawBytes(alphabetSize);
        r.align8();

        IntVector compressed = IntVector.decode(r);
        long totalChars = compressed.length;
        SdsReader.checkArraySize(totalChars, "StringArray 字符数据（大图的节点序列请用 GBZStreamer.readSequencesPacked() 打包读取）");
        byte[] chars = new byte[(int) totalChars];
        for (long i = 0; i < totalChars; i++) {
            chars[(int) i] = compToChar[(int) compressed.get(i)];
        }

        int n = (int) index.size();
        String[] result = new String[n];
        for (int i = 0; i < n; i++) {
            long start = index.values[i];
            long end = (i + 1 < n) ? index.values[i + 1] : totalChars;
            result[i] = new String(chars, (int) start, (int) (end - start), StandardCharsets.UTF_8);
        }
        return result;
    }

    public static void encodeStringArray(SdsWriter w, String[] items) {
        StringBuilder concat = new StringBuilder();
        long[] starts = new long[items.length];
        for (int i = 0; i < items.length; i++) {
            starts[i] = concat.length();
            concat.append(items[i]);
        }
        byte[] allChars = concat.toString().getBytes(StandardCharsets.UTF_8);

        // 只收集实际出现过的字节，得到字母表（comp_to_char）
        boolean[] seen = new boolean[256];
        for (byte b : allChars) seen[b & 0xFF] = true;
        int[] charToComp = new int[256];
        Arrays.fill(charToComp, -1);
        List<Byte> alphabet = new ArrayList<>();
        for (int c = 0; c < 256; c++) {
            if (seen[c]) {
                charToComp[c] = alphabet.size();
                alphabet.add((byte) c);
            }
        }
        byte[] compToChar = new byte[alphabet.size()];
        // 将Byte转位byte
        for (int i = 0; i < compToChar.length; i++) compToChar[i] = alphabet.get(i);

        int width = IntVector.widthFor(Math.max(alphabet.size() - 1, 0));
        long[] codes = new long[allChars.length];
        for (int i = 0; i < allChars.length; i++) codes[i] = charToComp[allChars[i] & 0xFF];

        SparseVector.build(starts).encode(w);
        w.usize(compToChar.length);
        w.rawBytes(compToChar);
        w.align8();
        IntVector.pack(codes, width).encode(w);
    }

    // ---------------- Dictionary（StringArray + 按字典序排列的下标） ----------------

    public static String[] decodeDictionaryNames(SdsReader r) {
        String[] names = decodeStringArray(r);
        IntVector.decode(r); // sortedIds：这里用不到具体查找功能，读掉即可，不然后面字段会错位
        return names;
    }

    public static void encodeDictionaryNames(SdsWriter w, String[] names) {
        encodeStringArray(w, names);
        Integer[] order = new Integer[names.length];
        for (int i = 0; i < names.length; i++) order[i] = i;
        Arrays.sort(order, (a, b) -> names[a].compareTo(names[b]));
        long[] sortedIds = new long[names.length];
        for (int i = 0; i < names.length; i++) sortedIds[i] = order[i];
        int width = IntVector.widthFor(Math.max(names.length - 1, 0));
        IntVector.pack(sortedIds, width).encode(w);
    }

    // ---------------- even 压缩版本（zstd，GBWTGraph 节点序列专用） ----------------

    public static String[] decodeEvenCompressed(SdsReader r) throws IOException {
        SparseVector index = SparseVector.decode(r);
        long totalLen = r.usize();
        byte[] compressedBytes = readByteVector(r);
        byte[] decompressed = new ZstdDecompressCtx().decompress(compressedBytes, (int) totalLen);
        if (decompressed.length != totalLen) {
            throw new IllegalStateException("zstd 解压后长度与文件声明不一致");
        }

        int n = (int) index.size();
        String[] result = new String[n];
        for (int i = 0; i < n; i++) {
            long start = index.values[i];
            long end = (i + 1 < n) ? index.values[i + 1] : totalLen;
            result[i] = new String(decompressed, (int) start, (int) (end - start), StandardCharsets.UTF_8);
        }
        return result;
    }

    public static void encodeEvenCompressed(SdsWriter w, String[] forwardSequences) throws IOException {
        StringBuilder concat = new StringBuilder();
        long[] starts = new long[forwardSequences.length];
        for (int i = 0; i < forwardSequences.length; i++) {
            starts[i] = concat.length();
            concat.append(forwardSequences[i]);
        }
        byte[] allBytes = concat.toString().getBytes(StandardCharsets.UTF_8);

        SparseVector.build(starts).encode(w);
        w.usize(allBytes.length);
        ZSTDCompressor compressor = new ZSTDCompressor();
        byte[] compressed = compressor.compress(allBytes).bytes();
        w.writeByteVector(compressed);
    }

    // ---------------- 小工具 ----------------

    private static byte[] readByteVector(SdsReader r) {
        long len = r.usize();
        byte[] data = r.rawBytes(len);
        r.align8();
        return data;
    }
}
