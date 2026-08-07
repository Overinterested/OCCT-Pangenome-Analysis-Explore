# GBZ 文件格式逐字节解析文档

> 本文档基于项目源码逐字段分析，描述 `simple_sds` 序列化格式的 GBZ 文件（`.gbz`），对应论文 2.3.5 节。文件底层是 **unsigned little-endian 64-bit integers** 构成的元素数组（`u64`），所有解析均从文件偏移 0 开始顺序消费。

---

## 基础约定

- **u64**: 8 字节 little-endian 无符号整数。文件里所有整数字段（包括长度前缀）都是 u64。
- **usize**: 语义上等同于 u64，用作长度/计数前缀。
- **对齐**: 字节向量（`Vec<u8>`）末尾有 0~7 字节的 padding 补齐到 8 字节边界。
- **字节序**: 全部 little-endian，由 [SdsWriter.java](src/sds/SdsWriter.java:16-20) 和 [SdsBufferReader.java](src/sds/SdsBufferReader.java:12) 保证。
- **无全局索引**: 文件没有 TOC 或中央目录。解析器读完一个模块，usize 长度字段自然决定了下一个模块的起始偏移。

---

## 模块总览

GBZ 文件由 10 个模块按固定顺序串联而成：

| 序号 | 模块 | 定长/变长 | 源码位置 |
|------|------|-----------|----------|
| 1 | GBZ Header | **定长 16 字节** | [GBZHeader.java](src/GBZHeader.java) |
| 2 | GBZ Tags | StringArray 变长 | [GBZTags.java](src/GBZTags.java) |
| 3 | GBWT Header | **定长 48 字节** | [GBWTHeader.java](src/gbwt/GBWTHeader.java) |
| 4 | GBWT Tags | StringArray 变长 | [GBWTTags.java](src/gbwt/GBWTTags.java) |
| 5 | BWT | SparseVector + ByteVector 变长 | [BWT.java](src/gbwt/bwt/BWT.java) |
| 6 | DASamples | usize 前缀 (0=无) + 内部模块 | [GBWTDASamples.java](src/gbwt/GBWTDASamples.java) |
| 7 | Metadata | usize 前缀 (0=无) + 内部模块 | [GBWTMeta.java](src/gbwt/meta/GBWTMeta.java) |
| 8 | GBWTGraph Header | **定长 24 字节** | [GBWTGraphHeader.java](src/gbwtgraph/GBWTGraphHeader.java) |
| 9 | Sequences | StringArray 或 zstd 压缩变长 | [GBWTGraphSequences.java](src/gbwtgraph/GBWTGraphSequences.java) |
| 10 | Translation | StringArray + SparseVector 变长 | [Translation.java](src/gbwtgraph/translation/Translation.java) |

解析总调度见 [GBZStreamer.java](src/GBZStreamer.java:36-45)，模块顺序由 `STAGE_NAMES` 数组定义，必须严格按序消费，乱序抛异常。

---

## 一、GBZ Header（定长 16 字节）

来源：[GBZHeader.java](src/GBZHeader.java)，[GBZStreamer.java:53-59](src/GBZStreamer.java)

**定长 16 字节，无长度前缀。** 这是整个文件的起点，从偏移 0 开始。

| 文件偏移 | 字段 | 字节数 | 解析方式 | 含义 |
|----------|------|--------|----------|------|
| 0 (byte 0-7) | `w0` | 8 | `r.u64()` | 低 32 位 = tag 魔数，高 32 位 = version |
| 8 (byte 8-15) | `flags` | 8 | `r.u64()` | 预留标志位，当前恒为 0 |

**解析步骤:**

1. 读 8 字节 u64 → `w0`。拆分：`tag = (int)(w0 & 0xFFFFFFFFL)`, `version = (int)(w0 >>> 32)`
2. 校验 `tag == 0x205A4247`（小端字节序下即 ASCII "GBZ "），不匹配则抛出 `IllegalStateException`
3. 读 8 字节 u64 → `flags`

**示例（y.giraffe.gbz 典型值）:**
```
offset 0:  47 42 5A 20 01 00 00 00   → tag=0x205A4247, version=1
offset 8:  00 00 00 00 00 00 00 00   → flags=0
```

---

## 二、GBZ Tags（StringArray 变长）

来源：[GBZTags.java](src/GBZTags.java)，[GBZStreamer.java:63-72](src/GBZStreamer.java)

key-value 标签对，存储为 **StringArray**。结构：

```
SparseVector(index)  → 每条 tag 字符串在拼接大串中的起始位置
usize(alphabetSize)  → 字母表大小
byte[alphabetSize]   → 压缩编码→原始字符映射表 (comp_to_char)
padding to 8-byte
IntVector(chars)     → 位压缩的字符编码序列
```

### StringArray 子结构详解

对应代码：[StringCodec.decodeStringArray](src/sds/StringCodec.java:31-47)

#### 步骤 2.1: 读 SparseVector (index)

SparseVector 由 3 部分组成，对应 [SdsPrimitives.SparseVector.decode](src/sds/SdsPrimitives.java:118-135):

| 读取操作 | 字段 | 字节数 | 含义 |
|----------|------|--------|------|
| `r.usize()` | `universe` | 8 | Elias-Fano 的全集大小 |
| 读 RawBitVector | `high` | 变长 | Elias-Fano 高位部分 |
| 读 IntVector | `low` | 变长 | Elias-Fano 低位部分 |

**RawBitVector 子结构**（[SdsPrimitives.RawBitVector.decode](src/sds/SdsPrimitives.java:60-71)）:

| 读取操作 | 字段 | 字节数 | 含义 |
|----------|------|--------|------|
| `r.usize()` | `ones` | 8 | 置位比特数 |
| `r.usize()` | `bitLength` | 8 | 位图总长度（bit 数） |
| `r.usize()` | `wordCount` | 8 | ceil(bitLength/64) |
| `r.rawWords(wordCount)` | `words` | wordCount×8 | 位图数据 |
| `r.usize()` + skip | 3 个空 option | 3×8 | rank/sel1/sel0（运行时重建） |

**IntVector 子结构**（[SdsPrimitives.IntVector.decode](src/sds/SdsPrimitives.java:24-35)）:

| 读取操作 | 字段 | 字节数 | 含义 |
|----------|------|--------|------|
| `r.usize()` | `length` | 8 | 元素个数 |
| `r.usize()` | `width` | 8 | 每个元素的位宽 (0..64) |
| `r.usize()` | `bitSize` | 8 | 必须 == length×width |
| `r.usize()` | `wordCount` | 8 | ceil(bitSize/64) |
| `r.rawWords(wordCount)` | `words` | wordCount×8 | 打包的位数据 |

**SparseVector 重建值**（从 high 和 low 算出）:
```
对于 i = 0..m-1:
  highPart = high.selectOne(i+1) - i
  values[i] = low.get(i) | (highPart << low.width)
```
其中 `low.width = max(1, floor(log2(universe / m)))`，`m = low.length`。

至此得到 `index.values[]`——一个升序 long 数组，表示每条字符串在总拼接串中的起始偏移。

#### 步骤 2.2: 读字母表

| 读取操作 | 字段 | 字节数 | 含义 |
|----------|------|--------|------|
| `r.usize()` | `alphabetSize` | 8 | 字母表大小 |
| `r.rawBytes(alphabetSize)` | `compToChar` | alphabetSize | 索引→原始字符的映射 |
| `r.align8()` | padding | 0~7 | 补齐到 8 字节 |

#### 步骤 2.3: 读压缩字符数据 (IntVector)

同步骤 2.1 中的 IntVector 结构。得到 `chars`，每个元素是压缩编码（0..alphabetSize-1）。

#### 步骤 2.4: 重建字符串

```
n = index.size()
for i = 0..n-1:
  start = index.values[i]
  end = (i+1 < n) ? index.values[i+1] : chars.length
  s[i] = new String(chars[start..end), 每个 char 查 compToChar 表还原)
```

#### 步骤 2.5: 组装 key-value

标签是 key-value 交替: `flat[0]=key0, flat[1]=val0, flat[2]=key1, flat[3]=val1, ...`
偶数条目的 StringArray 意味着无标签。

**跳过策略**（[SdsSkip.stringArray](src/sds/SdsSkip.java:31-36)）:
读 SparseVector(universe→RawBitVector→IntVector)、读 alphabetSize 并 skip、align8、再读并 skip IntVector 即可。

---

## 三、GBWT Header（定长 48 字节）

来源：[GBWTHeader.java](src/gbwt/GBWTHeader.java)，[GBZStreamer.java:76-82](src/GBZStreamer.java)

**定长 48 字节，无长度前缀。** 这个区域不是通过 usize 长度字段来界定边界的——它就是固定 6 个 u64。

| 文件偏移 | 读取操作 | 字段 | 含义 |
|----------|----------|------|------|
| +0  | `r.u64()` | `w0` | 低 32 位: tag=0x6B376B37, 高 32 位: version=5 |
| +8  | `r.u64()` | `sequences` | 输入序列总数 |
| +16 | `r.u64()` | `size` | BWT 总大小（字符数） |
| +24 | `r.u64()` | `offset` | 序列起始偏移 |
| +32 | `r.u64()` | `alphabetSize` | 字母表大小 |
| +40 | `r.u64()` | `flags` | 标志位 |

**flags 关键位:**
```
FLAG_BIDIRECTIONAL  = 0x0001  → 双向索引
FLAG_METADATA       = 0x0002  → 后面 Metadata 区域有数据
FLAG_SIMPLE_SDS     = 0x0004  → simple_sds 格式（本项目只支持此种）
```

解析完后校验 `isSimpleSds()`（即 `flags & 0x0004 != 0`），不满足说明是旧版 SDSL 格式，直接抛异常（[GBZStreamer.java:80](src/GBZStreamer.java)）。

---

## 四、GBWT Tags（StringArray 变长）

来源：[GBWTTags.java](src/gbwt/GBWTTags.java)，[GBZStreamer.java:86-95](src/GBZStreamer.java)

与 **区域二 (GBZ Tags)** 结构完全相同：一个 StringArray，key-value 交替。解析/跳过方式一样。不再重复。

---

## 五、BWT（SparseVector + ByteVector）

来源：[BWT.java](src/gbwt/bwt/BWT.java)，[GBZStreamer.java:99-110](src/GBZStreamer.java)

BWT 由两段组成：**BWTIndex (SparseVector)** + **BWTRecords (ByteVector)**，顺序固定（[BWT.decode](src/gbwt/bwt/BWT.java:50-53)）。

### 5.1 BWTIndex: SparseVector

解析方式同 2.1 节。得到 `recordStarts.values[]`：第 i 个值 = 第 i 条节点记录在 BWTRecords.data 中的起始字节偏移。

- `recordCount() = recordStarts.size()` = 节点数量
- 第 i 条记录的字节范围: `[recordStarts.values[i], recordStarts.values[i+1])`，最后一条到 data 末尾
- 代码见 [BWTIndex.range](src/gbwt/bwt/BWTIndex.java:23-26)

### 5.2 BWTRecords: ByteVector

| 读取操作 | 字段 | 字节数 | 含义 |
|----------|------|--------|------|
| `r.usize()` | `len` | 8 | 原始数据总字节数 |
| `r.rawBlocks(len)` 或 `r.rawBytes(len)` | `data` | len | 所有节点记录首尾相连的压缩数据 |
| `r.align8()` | padding | 0~7 | 补齐到 8 字节 |

**大文件注意**: BWT records 可能超过 2 GiB（hprc-v2.1 实测 3.16 GiB），使用 `rawBlocks` 分块读取（[ByteBlocks](src/sds/ByteBlocks.java)），而不是受 JVM 单数组上限约束的 `rawBytes`。

### 5.3 单条 BWT 记录的 RunCodec 解码

拿到 `[start, end)` 字节范围后，解码一条节点记录（[RunCodec.decodeRecord](src/sds/RunCodec.java:92-103)）：

#### 步骤 5.3.1: 读出边数 sigma (ByteCode 变长整数)

从 `data[start]` 开始读 LEB128 风格的变长整数（[RunCodec.readByteCode](src/sds/RunCodec.java:18-25)）：

```
byte b = data[pos]
result = b & 0x7F
shift = 0
while (b & 0x80):
  pos++
  shift += 7
  result += ((data[pos] & 0x7F) << shift)
pos++
```

例如:
- `0x05` → 值为 5，消费 1 字节
- `0x80 0x01` → 值为 128，消费 2 字节
- `0xE5 0x8E 0x26` → 值为 624869，消费 3 字节

#### 步骤 5.3.2: 读出边表

连续读 sigma 次，每次读 2 个 ByteCode:
- 第 1 个: 目标 node id **增量** (delta 编码)。当前边目标 = 前一个目标 + 增量（第一个增量的前值是 0）。
- 第 2 个: rank offset（目标节点记录内的游程序号偏移）

```java
long prev = 0;
for (int i = 0; i < sigma; i++) {
    long to = readByteCode(data, pos) + prev;  // 增量解码
    prev = to;
    long rankOffset = readByteCode(data, pos);
    outgoing.add(new long[]{to, rankOffset});
}
```

#### 步骤 5.3.3: 读游程体 (RunCodec)

如果 sigma == 0: 无游程，记录到此结束。

如果 sigma > 0: 剩余字节全是游程序列。每个游程编码为 `{outrank, runLength}`，用 `RunCodec` 解码。

**RunCodec 编码规则**（[RunCodec.read](src/sds/RunCodec.java:41-54)）：

情况 A: `sigma < 255`，`runContinues = 256 / sigma`

| 条件 | 编码 | 消费字节 |
|------|------|----------|
| `len < runContinues` | 1 字节: `code = value + sigma * (len-1)` | 1 |
| `len >= runContinues` | 1 字节: `code = value + sigma * (runContinues-1)` + ByteCode(len - runContinues) | 1 + ByteCode 字节数 |

解码: `code = data[pos] & 0xFF`，然后:
- `value = code % sigma`
- `len = code / sigma + 1`
- 若 `len >= runContinues`，再读一个 ByteCode 加到 len 上

情况 B: `sigma >= 255`，`runContinues = 0`

全走双变长码: `ByteCode(value)` + `ByteCode(len-1)`。每个游程消费 2 个 ByteCode 的字节数。

游程按顺序拼接即该节点完整的 BWT 子串。

**跳过单条记录（零分配）**: [RunCodec.skipRecord](src/sds/RunCodec.java:105-118)，只消费字节不构建对象，校验是否恰好消费到 limit 位置。

---

## 六、DASamples（Option，usize 前缀）

来源：[GBWTDASamples.java](src/gbwt/GBWTDASamples.java)，[GBZStreamer.java:114-122](src/GBZStreamer.java)

**首先读 8 字节 usize 作为"是否存在"标志：**

| 读取操作 | 字段 | 含义 |
|----------|------|------|
| `r.usize()` | `elements` | 0 = 不存在；>0 = 数据区 element 数 |

### 情况 1: elements == 0

DASamples 不存在，`present = false`。该区域结束，共消费 8 字节。

### 情况 2: elements > 0

DASamples 存在，`present = true`。后续为一个字节块，总字节数 = `elements * 8`（[GBWTDASamples.decode](src/gbwt/GBWTDASamples.java:31-36)）。在这个字节块内，按以下顺序解析：

| 子模块 | 结构 | 含义 |
|--------|------|------|
| `sampledRecords` | RawBitVector | 被采样的记录位图 |
| `bwtRanges` | SparseVector | BWT 区间映射 |
| `sampledOffsets` | SparseVector | 采样偏移 |
| `array` | IntVector | 采样数据数组 |

RawBitVector 结构见 2.1 节，SparseVector 结构见 2.1 节，IntVector 结构见 2.1 节。

**注意**: DASamples 在论文中描述为 "serialized as an optional structure in an unspecified format"，实际布局取决于写文件的 C++ 实现版本。如果只关心序列/拓扑/样本信息，可以直接跳过。

**跳过策略**（[GBZStreamer.skipDASamples](src/GBZStreamer.java:124-128)）:
```
elements = r.usize()
if (elements > 0) r.skipBytes(elements * 8)
```

**合法性检查**: `elements < 0 || elements > r.bytesRemaining() / 8` 时说明文件损坏，直接抛异常。

---

## 七、Metadata（Option，usize 前缀）

来源：[GBWTMeta.java](src/gbwt/meta/GBWTMeta.java)，[GBZStreamer.java:132-142](src/GBZStreamer.java)

**首先读 8 字节 usize 作为"是否存在"标志：**

| 读取操作 | 字段 | 含义 |
|----------|------|------|
| `r.usize()` | `elements` | 0 = 不存在；>0 = 数据区 element 数 |

### 情况 1: elements == 0

Metadata 不存在，`present = false`。与 GBWT header 的 `hasMetadata()` 必须一致。共消费 8 字节。

### 情况 2: elements > 0

Metadata 存在，`present = true`。后续 `elements * 8` 字节为一个独立的子数据区。在这个数据区内按顺序解析：

### 7.1 MetaHeader（定长 48 字节 = 6 × u64）

来源：[MetaHeader.java](src/gbwt/meta/MetaHeader.java)

与 GBWT Header 类似，6 个 u64 无前缀：

| 偏移 | 字段 | 含义 |
|------|------|------|
| +0  | `w0` | 低 32 位: tag=0x6B375E7A, 高 32 位: version=2 |
| +8  | `sampleCount` | 样本数 |
| +16 | `haplotypeCount` | 单倍型数 |
| +24 | `contigCount` | contig 数 |
| +32 | `flags` | 标志位 |

tag 校验同前。

### 7.2 MetaPaths（变长）

来源：[MetaPaths.java](src/gbwt/meta/MetaPaths.java)

| 读取操作 | 字段 | 字节数 | 含义 |
|----------|------|--------|------|
| `r.usize()` | `count` | 8 | 路径条数 |
| `count × r.u64()` | `paths` | count×16 | 每个 PathName 2 个 u64 |

**关键细节**: C++ 端用 `ShortPathName` 落盘（4 个 uint32），所以 Java 端 2 个 u32 打包成一个 u64 读取:

```
w0 = r.u64()  →  sample = w0 & 0xFFFFFFFFL          (低 32 位)
              →  contig = w0 >>> 32                  (高 32 位)
w1 = r.u64()  →  phase  = w1 & 0xFFFFFFFFL          (低 32 位)
              →  count  = w1 >>> 32                  (高 32 位)
```

**哨兵值**: `NO_VALUE = 0xFFFFFFFFL` 表示"无此值"（如参考路径 `_gbwt_ref` 的 phase）。

### 7.3 MetaSamples（Dictionary）

来源：[MetaSamples.java](src/gbwt/meta/MetaSamples.java)

Dictionary = StringArray + 字典序下标 IntVector（[StringCodec.decodeDictionaryNames](src/sds/StringCodec.java:70-73)）：

| 读取操作 | 字段 | 含义 |
|----------|------|------|
| 读 StringArray | `names` | sample id → 样本名字符串 |
| 读 IntVector | `sortedIds` | 按字典序排好的下标（解码时可忽略但必须读掉） |

### 7.4 MetaCotigs（Dictionary）

来源：[MetaCotigs.java](src/gbwt/meta/MetaCotigs.java)

结构同 7.3 MetaSamples：StringArray(contig 名) + IntVector(sortedIds)。

### Metadata 整体编码方式

Metadata 作为 Option，encode 时先把内部全部模块写进一个临时 `SdsWriter`，算出总字节数，然后写 `usize(totalBytes/8)` + `rawBytes(临时数据)`（[GBWTMeta.encode](src/gbwt/meta/GBWTMeta.java:33-40)）。所以 `elements` 字段实际表示内部子数据区占用了多少个 8 字节 element。

**跳过策略**（[GBZStreamer.skipMetadata](src/GBZStreamer.java:138-141)）:
```
elements = r.usize()
r.skipBytes(elements * 8)
```

---

## 八、GBWTGraph Header（定长 24 字节）

来源：[GBWTGraphHeader.java](src/gbwtgraph/GBWTGraphHeader.java)，[GBZStreamer.java:146-152](src/GBZStreamer.java)

**定长 24 字节，无长度前缀。** 3 个 u64:

| 文件偏移 | 读取操作 | 字段 | 含义 |
|----------|----------|------|------|
| +0  | `r.u64()` | `w0` | 低 32 位: tag=0x6B3764AF, 高 32 位: version |
| +8  | `r.u64()` | `nodes` | 节点总数 |
| +16 | `r.u64()` | `flags` | 标志位 |

**flags 关键位:**
```
FLAG_TRANSLATION = 0x0001  → 区域十 Translation 有数据
FLAG_SIMPLE_SDS   = 0x0002  → simple_sds 标记
FLAG_ZSTD         = 0x0002  → version >= 4 时: 节点序列用 zstd 压缩
```

**version 决定了序列区的解码方式**（[GBWTGraphHeader.usesZstdSequences](src/gbwtgraph/GBWTGraphHeader.java:27-28)）:
- `version < 4`: 序列区是普通 StringArray（字母表压缩）
- `version >= 4 && flags & 0x0002 != 0`: 序列区是 zstd 压缩

---

## 九、Sequences（变长）

来源：[GBWTGraphSequences.java](src/gbwtgraph/GBWTGraphSequences.java)，[GBZStreamer.java:156-159](src/GBZStreamer.java)

根据 `GBWTGraphHeader.version` 决定解析路径。

### 情况 A: 普通 StringArray（version ≤ 3）

解析方式同区域二 (GBZ Tags) 中的 StringArray。字母表通常是 ACGT 或 ACGTN。

解析结果: `forward[i]` = 第 i+1 个节点的正向 DNA 序列（节点从 1 开始编号）。

反向序列不存储，按需计算 `reverseComplement(forward[i])`（[GBWTGraphSequences.reverseComplement](src/gbwtgraph/GBWTGraphSequences.java:33-38)）。

**大文件优化**: `readSequencesPacked()` 不走 StringArray 全量解析，而是边读边转写成 `PackedSequences`（DNACodec 带标志位的 int 流），避免生成几千万个 String 对象（[GBZStreamer.readSequencesPacked](src/GBZStreamer.java:165-173)）。

### 情况 B: zstd 压缩（version ≥ 4 且 FLAG_ZSTD 置位）

来源：[StringCodec.decodeEvenCompressed](src/sds/StringCodec.java:79-91)

| 读取操作 | 字段 | 含义 |
|----------|------|------|
| 读 SparseVector | `index` | 每条序列在解压后大串中的起始位置 |
| `r.usize()` | `totalLen` | 解压后所有序列拼接的总字节数 |
| 读 ByteVector | `compressed` | zstd 压缩的原始字节 |

解析流程:
1. 读 SparseVector (同 2.1 节)
2. 读 usize = 解压后总长度
3. 读 usize = 压缩后长度，再读对应字节（ByteVector 格式），然后 align8
4. zstd 解压 → 得到 totalLen 字节的大字符串，所有正向序列首尾相连
5. 用 index.values[] 定位每条序列的边界

---

## 十、Translation（变长）

来源：[Translation.java](src/gbwtgraph/translation/Translation.java)，[GBZStreamer.java:181-188](src/GBZStreamer.java)

两段：[Segments](src/gbwtgraph/translation/Segments.java) + [Mapping](src/gbwtgraph/translation/Mapping.java)，顺序固定。

### 10.1 Segments（StringArray）

StringArray，同 2.1 节解析。`segments.names[i]` = 第 i 个 segment 名字符串（如 "chr1", "chr2"）。

### 10.2 Mapping（SparseVector）

SparseVector，同 2.1 节解析。`nodeToSegment.values[i]` = 第 i 个 segment 对应的第一个 node id。

**查询节点所属 segment**（[Mapping.segmentIndexForNode](src/gbwtgraph/translation/Mapping.java:19-22)）:
```
idx = binarySearch(nodeToSegment.values, nodeId)
// 找最后一个 <= nodeId 的起始位置
if (idx < 0) idx = -idx - 2
return idx  // segments.names[idx] 即该节点所属 segment 名
```

**无 translation 时**: `GBWTGraphHeader.hasTranslation() == false`，Segments 是空 StringArray，Mapping 是空 SparseVector。此时 segment 名约定为 node id 本身的字符串。

**跳过策略**（[GBZStreamer.skipTranslation](src/GBZStreamer.java:190-193)）:
```
SdsSkip.stringArray(r);   // Segments
SdsSkip.sparseVector(r);  // Mapping
```

---

## 附录 A: RunCodec 完整编码表

来源：[RunCodec.java](src/sds/RunCodec.java)

### ByteCode (LEB128 风格变长整数)

每个 ByteCode 消费的字节数由最高位决定:

| 值范围 | 编码字节数 | 示例 |
|--------|------------|------|
| 0..127 | 1 | `0x05` → 5 |
| 128..16383 | 2 | `0x80 0x01` → 128 |
| 16384..2097151 | 3 | `0xE5 0x8E 0x26` → 624869 |
| ... | ... | ... |

### 单节点记录结构

```
[ByteCode: sigma]                           -- 出边数
[for i=0..sigma-1:
    ByteCode: delta_node_id                 -- 目标节点增量编码
    ByteCode: rank_offset]                  -- rank 偏移
[for each run:
    RunCodec 编码的 {outrank, runLength}]    -- BWT 子串游程
```

### RunCodec 游程编码（sigma < 255 时）

设 `sigma = 7`, `runContinues = 256/7 = 36`:

| runLength | outrank | 编码字节 | 说明 |
|-----------|---------|----------|------|
| 1 | 3 | `0x03` | code = 3 + 7×0 = 3 |
| 1 | 5 | `0x05` | code = 5 + 7×0 = 5 |
| 35 | 2 | `0xF7` | code = 2 + 7×34 = 240 = 0xF0, 但 len=35：code = 2 + 7×(35-1) = 240 = 0xF0 |
| 36 | 1 | `0xFA 0x00` | 基础 code = 1 + 7×35 = 246 = 0xF6, 续码 = 36-36 = 0 |
| 100 | 4 | `0xFC 0x40` | 基础 code = 4 + 7×35 = 249 = 0xF9, 续码 = 100-36 = 64 |

---

## 附录 B: DNA 碱基编码表

来源：[DNACodec.java](src/sds/DNACodec.java)

| 碱基 | 编码值 | TWO_BIT 位宽 | THREE_BIT 位宽 |
|------|--------|-------------|----------------|
| A | 0 | 2 bit | 3 bit |
| C | 1 | 2 bit | 3 bit |
| G | 2 | 2 bit | 3 bit |
| T | 3 | 2 bit | 3 bit |
| N | 4 | 不允许 | 3 bit |

每个 int 数据结构: `bit31=标志位(1=存满) | bits[30..0]=碱基序列`。
TWO_BIT 每 int 存 15 个碱基，THREE_BIT 每 int 存 10 个。

---

## 附录 C: SDS 基础原语总表

| 原语 | 定长/变长 | 头部字段 | 数据字段 |
|------|-----------|----------|----------|
| u64/usize | 定长 8 字节 | 无 | 8 字节 LE |
| IntVector | 变长 | usize×4 | u64[wordCount] |
| RawBitVector | 变长 | usize×3 | u64[wordCount] + usize×3(空) |
| SparseVector | 变长 | usize(universe) | RawBitVector + IntVector |
| StringArray | 变长 | SparseVector + usize + bytes | IntVector |
| ByteVector | 变长 | usize(len) | bytes[len] + padding |
| Dictionary | 变长 | StringArray | IntVector |

---

## 附录 D: 魔术编码汇总

| 区域 | Tag (16进制) | ASCII | 源码 |
|------|-------------|-------|------|
| GBZ Header | `0x205A4247` | `GBZ ` | [GBZHeader.java:11](src/GBZHeader.java) |
| GBWT Header | `0x6B376B37` | `k7k7` | [GBWTHeader.java:12](src/gbwt/GBWTHeader.java) |
| MetaHeader | `0x6B375E7A` | `k7^z` | [MetaHeader.java:13](src/gbwt/meta/MetaHeader.java) |
| GBWTGraph Header | `0x6B3764AF` | `k7d¯` | [GBWTGraphHeader.java:26](src/gbwtgraph/GBWTGraphHeader.java) |

所有 tag 均存储在小端字节序的低 32 位，高 32 位是 version。

---

## 附录 E: 实际文件示例（y.giraffe.gbz）

一个真实的小型 GBZ 文件的字节布局示意：

```
offset      内容                    区域
------      ----                    ----
0x000000    47 42 5A 20 01 00...   GBZ Header (16 bytes)
0x000010    00 00 00 00 00 00...   GBZ Tags (StringArray)
...
            6B 37 6B 37 05 00...   GBWT Header (48 bytes)
            00 00 00 00 00 00...   GBWT Tags (StringArray)
            ...                    BWT (SparseVector + ByteVector)
            00 00 00 00 00 00...   DASamples (usize=0, 8 bytes)
            00 00 00 00 00 00...   Metadata (usize=0, 8 bytes)
            6B 37 64 AF 03 00...   GBWTGraph Header (24 bytes)
            ...                    Sequences (StringArray)
            ...                    Translation (StringArray + SparseVector)
```

可以用以下命令查看文件头部的 tag:
```
xxd -l 16 y.giraffe.gbz        # GBZ header
xxd -s <offset> -l 48 y.giraffe.gbz  # GBWT header
```
