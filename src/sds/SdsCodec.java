package sds;

import java.io.IOException;

/**
 * 全项目统一的 encode/decode 接口，对应论文里 C++ 那套
 * simple_sds_serialize(out) / simple_sds_load(in) 的写法：
 *   decode(r) —— 从字节流里把字段填进 this
 *   encode(w) —— 把 this 当前的字段写成字节流
 * 两者互为逆操作，这也是为什么整个项目能做"读一遍再原样写回去"这种回归测试。
 */
public interface SdsCodec {
    void decode(SdsReader r) throws IOException;
    void encode(SdsWriter w) throws IOException;
}
