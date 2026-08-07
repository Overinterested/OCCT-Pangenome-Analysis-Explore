# 	GBWT 与 GBWTGraph 详解

> 基于 hprc-v2.1-mc-grch38.gbz 真实数据，结合 [SERIALIZATION.md](https://github.com/jltsiren/gbwtgraph/blob/main/SERIALIZATION.md) 规范和本项目 Java 源码。

---

## 一、GBWT：图的拓扑引擎

### 1.1 GBWT 存的是什么

GBWT 不存碱基序列，它只存**图的拓扑结构**——谁连谁、怎么连。

可以把 GBWT 理解成一张图的压缩邻接表。图中的顶点叫**节点（node）**，每个节点是一个整数 ID。GBWT 为每个节点存一条**记录（record）**，记录里有两样东西：

1. **出边表**：从这个节点出发可以去哪些其他节点
2. **BWT 子串**：所有穿过该节点的路径在这里被压缩成游程

### 1.2 什么是"节点"

在 pangenome 的 de Bruijn 图里，一个节点 = 基因组中一段**在所有样本里完全相同的 DNA 序列片段**。如果某个位置有变异（SNP/indel），就会分叉成多个节点。

**实例**：假设 3 个人的基因组在同一个位置分别是：

```
样本A: ...AAG | C T A| GCA...
样本B: ...AAG | C T G| GCA...
样本C: ...AAG | C - -| GCA...   ← 缺失 TA
```

这会产生 5 个节点（用 `|` 分隔）：

```
node 1: AAG（三个样本共用）
node 2: C（三个样本共用）
node 3: TA（只有 A 有）     ─┬─ 分叉区域
node 4: TG（只有 B 有）     ─┤
node 5: GCA（三个样本又合并）─┘
```

每个节点在 GBWT 里有一条记录，记录描述了"从该节点可以去哪里"和"哪些路径穿过了它"。

### 1.3 如何存储

GBWT 由三部分组成：

#### 1.3.1 GBWTHeader（48 字节，定长）

| 字段 | 含义 | 示例值 |
|------|------|--------|
| tag | 魔数 0x6B376B37 | — |
| version | 格式版本 | 5 |
| sequences | 输入序列总数 | 106,300（= 53,150 条路径 × 2 方向） |
| size | BWT 总字符数 | 78,199,541,546 |
| offset | 节点起始编号 | 1（节点从 1 开始） |
| alphabetSize | 字母表大小（= 节点总数 + offset） | 425,853,422 |
| flags | 标志位 | 0x7（双向 + 有 metadata + simple_sds） |

#### 1.3.2 BWTIndex（SparseVector，Elias-Fano 编码）

存每条记录在 data 中的**起始字节偏移**。425M 条记录 → 425M 个 long 值，Elias-Fano 编码约 230 MB。

获取第 i 条记录：
```java
BWTIndex index = bwt.index;
long[] range = new long[2];
index.rangeInto(i, records.data.length, range);
// range[0] = 起始字节, range[1] = 结束字节（或 data.length）
```

#### 1.3.3 BWTRecords（ByteVector，3.16 GB）

所有 425M 条记录的首尾相连的**压缩字节流**。每条记录编码：

```
[ByteCode: sigma]                          ← 出边数
[sigma 对 × (ByteCode: delta_node, ByteCode: rank_offset)]  ← 出边表
[RunCodec 游程序列]                         ← BWT 子串
```

解码记录：
```java
byte[] slice = records.data.copyRange(range[0], range[1]);
RunCodec.NodeRecord record = RunCodec.decodeRecord(slice, 0, slice.length);
// record.outgoing: 出边列表 [{目标node, rank_offset}, ...]
// record.runs: 游程列表 [{outrank, length}, ...]
```

### 1.4 出边表和游程体详解

以 hprc-v2.1-mc-grch38.gbz 的真实节点为例：

```
node 0 (sentinel): σ=20502, 游程不展开
  含义：所有 53,150 条路径的终点都汇集到 sentinel（类似字符串末尾的 $）
  20502 = 各路径的不同"终止模式"数

node 1-19: 大部分 σ=1
  含义：这些节点在基因组中非常保守，只有一种"出去的路"
  游程体：所有路径在这里的 BWT 字符都是同一个值
```

一个更典型的例子——假设某个多态位点的节点：

```
出边表: [{→node 12345, rank=0}, {→node 67890, rank=5}]
游程体: [{outrank=0, len=150}, {outrank=1, len=83}]

含义：
  - 该节点有 2 种走出去的方式
  - 150 条路径选第 0 条边（→node 12345）
  - 83 条路径选第 1 条边（→node 67890）
  - 总共 233 条路径穿过该节点（= 样本数）
```

### 1.5 GBWT 的用途：序列到图的比对

**场景**：给一条 query 序列 `ACGTAC`，找到它在 pangenome 图中的位置。

**步骤（概念层面）**：

```
1. 找种子节点：用 k-mer 索引找到 query 开头的 k-mer（如 ACGTA）所在的节点集合

2. 正向扩展：从种子节点出发，通过出边表找到下一层节点
   - 对每个候选节点，用 GBWTGraph 查节点序列（见第二部分）
   - 如果某节点的序列前缀匹配 query 的下一段，继续扩展
   - 不匹配的剪枝

3. 反向扩展：同理，从种子节点反向遍历（利用双向 BWT）

4. 定位：最终找到的节点 ID 通过 Translation（见 2.4 节）映射回
   染色体坐标，就知道 query 在哪个样本的哪条染色体的哪个位置
```

**实现代码草稿**（正向扩展一步）：

```java
// 给定当前节点 curNode，向后走一步到达的节点集合
RunCodec.NodeRecord rec = bwt.getRecord(curNode);
for (long[] edge : rec.outgoing) {
    long nextNode = edge[0];     // 目标节点
    long rankOffset = edge[1];   // 目标节点记录中的 BWT rank 偏移
    // 走到 nextNode，其 BWT 子串中第 rankOffset 个字符对应 curNode
    String nextSeq = seq.getSequence(nextNode, true);
    // 用 nextSeq 和 query 的下一段比对...
}
```

---

## 二、GBWTGraph：节点的序列标签

### 2.1 GBWTGraph 存的是什么

GBWTGraph 存的是**每个节点的 DNA 序列**和**节点到染色体坐标的映射**。

GBWT 只存拓扑（谁连谁），不存序列。要拿到某个节点的碱基串，需要查 GBWTGraph。

### 2.2 GBWTGraphSequences / PackedSequences

存所有节点（包括空节点）的正向 DNA 序列。下标 = nodeId - 1（offset=1 时）。

**两种存储形态：**

| 形态 | 适用场景 | 内存 |
|------|---------|------|
| `String[] forward` | 小文件 | 每字符 2 字节（char） |
| `PackedSequences packed` | 大文件（本项目默认） | 每碱基 2-3 bit |

**PackedSequences 内部结构：**

- `format`: `TWO_BIT`（只含 ACGT）或 `THREE_BIT`（含 N）
- `starts[i]`: 第 i 条序列在全局 int 流中的起始下标
- `lengths[i]`: 第 i 条序列的碱基数
- `blocks[][]`: 16M-int 分块存储，突破 JVM 2^31 数组上限

**获取序列：**

```java
// 方式1：整串解码（新建 String，适合少量查询）
String seq = packed.decode(nodeId);

// 方式2：零分配逐碱基（热路径，如 k-mer 遍历）
int[] buf = new int[DNACodec.wordsFor(packed.length(nodeId), packed.format)];
packed.copyWords(nodeId, buf, 0);
char base = DNACodec.charAt(buf, pos, packed.format);
```

**反向互补：**

GBWTGraph 文件里只存正向序列。反向互补按需计算：

```java
String rev = GBWTGraphSequences.reverseComplement("ACGT");  // → "ACGT"
// A→T, C→G, G→C, T→A
```

### 2.3 真实数据示例

hprc-v2.1-mc-grch38.gbz：

```
序列数: 212,926,710（含空节点占位）  ← 比 graph.nodes(140M) 多，因为空节点也有占位
format: THREE_BIT  总碱基: 3,417,672,786

前10条全是 NNNN...  ← sentinel 及 padding 节点
实际有意义的序列从某个非 N 节点开始（如 chr1 的第一个 variant 位点）
```

空节点的序列是空字符串——它虽然在 BWT 里有记录（σ=0），但没有任何路径穿过它，所以 GBWTGraph 给它一个占位。

### 2.4 Translation：节点 → 染色体坐标

GBWT 的节点 ID 是整数（1, 2, 3, ...），但 GFA 格式需要字符串 segment 名（如 `chr1`, `chr2`）。Translation 做这个映射。

**Segments**：segment 名字符串数组
```
segments[0] = ""
segments[1] = ""  
...（很多空 segment，对应未放置的 scaffold）
segments[k] = "chr1"
segments[k+1] = "chr2"
...
```

**Mapping**：Elias-Fano 稀疏位图，存第 i 个 segment 对应的**起始 node ID**。

```
mapping.values[0] = 1     ← segment[0] 对应 node 1 开始
mapping.values[1] = 1000  ← segment[1] 对应 node 1000 开始
...
```

查询 node 12345 属于哪个 segment：
```java
int segIdx = mapping.segmentIndexForNode(12345);
String name = segments.names[segIdx];  // → "chr1"
```

**注意**：hprc-v2.1-mc-grch38.gbz 的 Translation 为空（0 segment）。这是因为节点 ID 本身就是按 contig 顺序排列的，不需要额外映射。hg38 参考基因组的节点 ID 可以直接按区间推导染色体坐标。

---

## 三、综合示例：从 query 序列到基因组坐标

### 场景

你有一条 query 序列 `GCAATC...`（来自某个测序 read），想知道它在 pangenome 图的哪个位置，以及它存在于哪些样本中。

### 步骤

#### 1. 加载文件

```java
GBZFile file = GBZFile.parse("hprc-v2.1-mc-grch38.gbz");
BWT bwt = file.gbwt.bwt;
PackedSequences seq = file.graph.sequences.packed;
```

#### 2. 找种子节点（k-mer 查找，简化版）

假设你已有 k-mer 索引。query 的前 k-mer `GCAA` 在 node 5000000：

```java
long curNode = 5_000_000;
int pos = 0;  // query 中已匹配的位置
```

#### 3. 沿 BWT 边扩展

```java
while (pos < query.length()) {
    String nodeSeq = seq.decode((int) curNode - 1);  // 节点序列（offset=1）
    int matchLen = matchPrefixLength(query, pos, nodeSeq);
    pos += matchLen;

    if (pos >= query.length()) break;  // query 全部匹配

    // 找下一层节点：query 的下一段应该在出边的某个目标节点里
    RunCodec.NodeRecord rec = bwt.getRecord((int) curNode);
    boolean found = false;
    for (long[] edge : rec.outgoing) {
        long nextNode = edge[0];
        String nextSeq = seq.decode((int) nextNode - 1);
        if (nextSeq.length() > 0 && nextSeq.charAt(0) == query.charAt(pos)) {
            curNode = nextNode;
            found = true;
            break;
        }
    }
    if (!found) break;  // 图中没有匹配
}
```

#### 4. 从 BWT 游程中看哪些样本经过该节点

```java
RunCodec.NodeRecord rec = bwt.getRecord((int) curNode);
// rec.runs 中的游程长度之和 = 穿过该节点的路径总数
long totalPaths = 0;
for (long[] run : rec.runs) totalPaths += run[1];
System.out.println("该节点被 " + totalPaths + " 条路径穿过");

// 如果 totalPaths == 样本数 × 2（双向），说明所有样本在该位置一致
// 如果 totalPaths < 样本数 × 2，说明该节点是某个样本特有的 variant
```

#### 5. 查样本名

```java
// 找到某条路径对应的样本
MetaPaths.PathName path = file.gbwt.meta.metaPaths.paths.get(somePathId);
String sample = file.gbwt.meta.metaSamples.get((int) path.sample);
String contig = file.gbwt.meta.metaCotigs.get((int) path.contig);
System.out.println("存在于样本 " + sample + " 的 " + contig);
```

---

## 四、Metadata：结构化元信息

### 4.1 每个字段的含义

| 字段 | 含义 | 示例 |
|------|------|------|
| `sampleCount` | 样本数 | 233 |
| `haplotypeCount` | 单倍型总数（每人每个 contig 1 条） | 53,150 |
| `contigCount` | 不同 contig 名总数 | 18,089 |

### 4.2 PathName 的含义

每条 PathName = `{sample, contig, phase, fragment}`：

| 字段 | 含义 | 示例 |
|------|------|------|
| sample | 样本索引（查 sample 名字典） | 0 → CHM13 |
| contig | contig 索引（查 contig 名字典） | — → chr1 |
| phase | 单倍型编号（0/1） | 0 = 父源, 1 = 母源 |
| fragment | 序列起始偏移（bp） | 12668 |

一条 path = "样本 X 在 contig Y 上的第 Z 号单倍型，从第 fragment bp 开始"。

多个 fragment 覆盖同一样本同一 contig 是因为基因组被分成多段（分块 build GBWT 的结果）。

---

## 附录：关键类速查

| 类 | 一句话作用 |
|---|---|
| `BWT` | 图拓扑：出边 + BWT 子串 |
| `BWTIndex` | 存每条记录在 data 中的起始偏移 |
| `BWTRecords` | 所有记录首尾相连的压缩字节 |
| `RunCodec.NodeRecord` | 单节点解码结果：出边表 + 游程列表 |
| `GBWTGraphSequences` | 节点 DNA 序列（String[] 或 PackedSequences） |
| `PackedSequences` | 2/3 bit 打包的 DNA 序列，零分配访问 |
| `DNACodec` | 从打包 int 中逐碱基解码 |
| `GBWTMeta` | 样本名、contig 名、路径信息 |
| `Segments` | segment 名字符串数组 |
| `Mapping` | node → segment 的区间映射 |
