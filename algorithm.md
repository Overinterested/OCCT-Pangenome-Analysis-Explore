# GBZ2OCCT 算法详解

## 1. 总体架构

```
                        GBZ 文件 (5.5 GB)
                             │
              ┌──────────────┼──────────────┐
              ▼              ▼              ▼
         loadGBZ()    computeInDegGt1()   dfsContig() × N
         (BWT + seq)   (出入度 BitSet)    (并行 per-contig)
                             │              │
                             └──────┬───────┘
                                    ▼
                              OCCT 文件组
```

算法分为三个阶段，其中 Phase 3 支持多线程并行（按 contig 或人造区间分片）。

---

## 2. Phase 1：`loadGBZ()` — 加载 GBZ

### 2.1 加载内容

| 模块 | 用途 | 大小 (HPRC v2.1) |
|------|------|-------------------|
| BWT (`bwt`) | 节点出边表 + BWT 子串（haplotype 游程） | 3.4 GB |
| PackedSequences (`sequences`) | 节点 DNA 序列（2/3-bit 压缩） | 1.4 GB |
| Translation (`translation`) | 节点→contig 映射（可选） | 0~数 MB |

### 2.2 Contig 分组策略

```java
if (hasTranslation && segments <= 50000) {
    buildContigRanges();     // 从 Translation 提取 contig 边界
} else {
    buildArtificialRanges(); // 按节点 ID 均匀分片
}
```

**有 Translation：** `Mapping.nodeToSegment.values` 是升序的 segment 起始节点 ID 数组。第 i 个 segment 的节点范围是 `[values[i], values[i+1] - 1]`。

**无 Translation（如 HPRC v2.1）：** 将 `[1, nodeCount]` 均分为 `numThreads` 个区间。每个区间内节点 ID 连续，拓扑可能跨区。

### 2.3 内存优化

- `gbwt = null`：加载完 BWT 后立即释放 GBWT 包装对象
- `seq = null`：加载完 PackedSequences 后释放中间对象
- `System.gc()`：显式触发 GC 回收中间分配

---

## 3. Phase 2：`computeInDegGt1()` — 计算出入度

### 3.1 算法

```
输入: BWT 记录 (每个节点一条正向 + 一条反向)
输出: inDegGt1 (入度 > 1 的节点集合), inDegZero (入度 = 0 的节点集合)

for each BWT record i in [0, bwt.recordCount()):
    sigma = readByteCode(record)    // 出边数
    for each edge (target, rank_offset):
        target_node = target >> 1   // node_type → node_id
        inDegCount[target_node]++   // 入度 +1
        inDegZero.clear(target_node) // 该节点有入边，移出零入度集合
        if inDegCount[target_node] == 2:
            inDegGt1.set(target_node) // 入度 > 1，标记
```

### 3.2 关键细节

**`to >> 1` 的必要性：** GBWT 的边存储的是 `node_type = (node_id << 1) | orientation`。低 1 位是方向（0=正向, 1=反向），高位是 node_id。`>> 1` 提取纯 node_id。这是之前导致覆盖率从 0.06% 跳到 100% 的 bug 修复。

**扫描所有记录：** 遍历 `i = 0..totalRecords-1`（包含正向和反向记录），因为反向记录中的边也会贡献目标节点的入度。

**`inDegCount` 只计数到 2：** 对于 block 边界判定，只需区分入度 =0、=1、>1 三种情况。`byte[]` 节省内存（213 MB vs short[] 的 426 MB）。

### 3.3 输出含义

- `inDegGt1`：入度 > 1 的节点 → merge point，多个路径在此汇聚
- `inDegZero`：入度 = 0 的节点 → path start，无前驱节点

---

## 4. Phase 3：`dfsContig()` — 核心 DFS + 贪心 Block 分解

### 4.1 数据结构

| 变量 | 类型 | 大小 | 用途 |
|------|------|------|------|
| `visited` | BitSet | 27 MB (213M bits) | 已访问节点 |
| `inPath` | BitSet | 27 MB | 当前 DFS 路径上的节点（环检测） |
| `stack` | long[] | 8 MB 初始，按需扩容 | DFS 栈，每元素编码 (node_id, isBlockStart) |
| `recordBuf` | byte[] | 256 KB 初始，按需扩容 | BWT 记录解析缓冲区（线程局部） |
| `et`, `eo`, `ro` | long[], byte[], int[] | ~2 KB 初始 | 当前节点的出边表、方向、rank_offset |

### 4.2 栈元素编码

```
long stackEntry = node_id | (isBlockStart ? 0x80000000L : 0)
                   └──31 bits──┘   └── 1 bit ──┘
```

每个 8 字节的 `long` 同时存储节点 ID（低 31 位）和 block 起始标记（最高位）。相比 `Stack<Node>` 对象方式节省 ~10× 内存。

### 4.3 起始节点选择（三步）

```
Step 1: Sentinel children（哨兵子节点）
  解析 BWT 记录 0（哨兵），其出边指向所有 haplotype 路径的终点。
  筛选在 [minNode, maxNode] 范围内且未访问的节点入栈。

Step 2: In-degree-0 节点
  遍历 inDegZero 中 [minNode, maxNode] 范围内的所有置位。
  这些是没有入边的路径起点，直接入栈。

Step 3: 回退扫描（nextClearBit）
  当 Step 1 和 Step 2 耗尽后，如果仍有未访问节点，
  通过 visited.nextClearBit(minNode) 找到第一个未访问节点入栈。
  这覆盖了孤立连通分量和入度 >0 但未被 Step 1 覆盖的节点。
```

### 4.4 主循环——节点处理

```
while (true):
    if 栈空:
        结算当前 block，回到起始节点选择 Step 3
        if 无未访问节点: break

    node = stack.pop()
    解码: node_id = packed & 0x7FFFFFFF
          isBlockStart = packed 最高位

    if node 在 inPath 或 visited: continue  // 环截断 或 已处理

    visited.set(node_id)
    inPath.set(node_id)

    // === 解析 BWT 正向记录 ===
    record = bwt.index.rangeInto(node_id * 2)  // 正向记录
    sigma = ByteCode(record)                    // 出边数
    for each edge:
        to = ByteCode(record) + prevTarget      // delta 解码
        rank_offset = ByteCode(record)
        edgeTargets[e] = to >> 1                // node_type → node_id
        edgeOrientations[e] = (to & 1) == 0 ? 2 : 3  // 0→LR, 1→RR

    // === 解析 rank_offset 数组 ===
    runs = RunCodec.decode(record)  // 游程 {(outrank, length), ...}
    rankOffsets = prefixSum(runs)   // [0, cum1, cum2, ...]

    // === Block 边界判定 ===
    flags = 0
    if (!blockStarted || isBlockStart):
        flags |= BLOCK_START; nodesInBlock = 0; blockStarted = true
    nodesInBlock++

    endsBlock = (sigma != 1)  // 分支点(>1) 或 汇点(0)
    if (!endsBlock && sigma == 1):
        child = edgeTargets[0]
        if child 不在 [minNode, maxNode] 或 child 已访问:
            endsBlock = true  // merge point
    if endsBlock: flags |= BLOCK_END

    // === 写入 OCCT 记录 ===
    writer.writeRecord(node_id, seqLen, sigma, flags,
                       packedSeq, edgeTargets, edgeOrientations,
                       contig_id=-1, position=-1, rankOffset)

    // === 推送子节点 ===
    if sigma == 1:
        push(child, isBlockStart=endsBlock)  // 延续或新 block
    elif sigma > 1:
        for child in reverse(edgeTargets):   // 逆序保证 DFS 前序
            push(child, isBlockStart=true)   // 分支总是新 block

    inPath.clear(node_id)
    if endsBlock: 结算当前 block
```

### 4.5 贪心 Block 分解算法

核心思想：**线性链（unitig）聚合为连续 block，分支点和 merge point 作为 block 边界。**

```
判定规则:
  ┌─────────────────────────────────────────────────────────┐
  │ sigma = 0 (汇点)                                        │
  │   → 当前节点结束 block（孤点 block）                     │
  │                                                         │
  │ sigma = 1 (单出边)                                      │
  │   child 未访问 且 在区间内                               │
  │     → 线性延续，子节点同 block                          │
  │   child 已访问 或 越界                                  │
  │     → merge point，当前节点结束 block                   │
  │                                                         │
  │ sigma > 1 (分支点)                                      │
  │   → 当前节点结束 block                                  │
  │   → 每个分支子节点开始新的 block                        │
  └─────────────────────────────────────────────────────────┘
```

与 vg 的 snarl decomposition 对比：

| | snarl decomposition | 贪心出入度 |
|---|---|---|
| 算法 | 仙人掌图 → 桥森林 | DFS + 度数判定 |
| 时间复杂度 | O(N + M)，常数大 | O(N + M)，常数小 |
| 边界精度 | 精确（数学保证） | 近似 |
| 嵌套变异 | 正确识别 | 需递归处理 |
| 环处理 | 桥森林自然处理 | inPath 截断 |

**为什么使用 DFS visited 而非预计算的 inDegGt1 做 merge 判定？**

预计算的 `inDegGt1` 包含所有方向的入度（正向 + 反向），会导致过多节点被判为 merge point（几乎所有单边节点都终止 block）。使用运行时 `visited.get(child)` 仅当子节点在当前 DFS 的已完成路径中才触发 block 终止——这更精确地反映正向 traversal 中的 merge 语义。

### 4.6 环处理

```
if (inPath.get(node)):
    continue  // 节点已在当前 DFS 路径上 → 环！跳过不处理
```

当子节点已存在于当前路径时，DFS 处于一个环中。此时不推送子节点（截断），当前节点仍正常写入。环的"回边"以正常出边形式存储在 OCCT 记录中，只是 DFS 不跟踪它。

### 4.7 跨区边处理（含并行分片时）

```
if (child < minNode || child > maxNode):
    // 子节点不在当前分片范围内
    // 边信息完整保留在 OCCT 记录中
    // 但 DFS 不跟踪，block 在此结束
```

跨区边在 OCCT 记录中完整保留（`edgeTargets` 包含越界目标），但 DFS 将其视为 merge point 来结束当前 block。这是为了实现并行化所做的妥协——后续需要合并步骤来修复被切断的 chain。

---

## 5. 并行化设计

### 5.1 `convertParallel()`

```java
ExecutorService pool = Executors.newFixedThreadPool(threads);
contigs.sort((a, b) -> Long.compare(b.size, a.size));  // 大的先处理
for (ContigRange cr : contigs):
    pool.submit(() -> dfsContig(cr.start, cr.end, cr.name + ".occt", cr.name));
pool.shutdown();
for (Future f : futures): f.get();  // 等待全部完成
```

### 5.2 线程安全分析

| 资源 | 共享方式 | 安全 |
|------|---------|------|
| `bwt` | 只读 | 安全 |
| `sequences` | 只读 | 安全 |
| `inDegGt1`, `inDegZero` | 只读（Phase 2 完成后不变） | 安全 |
| `recordBuf` | 线程局部变量 | 安全（修复前有竞态 bug） |
| `visited`, `inPath`, `stack` | 每个线程独立 | 安全 |
| `OCCTWriter` | 每个线程独立 | 安全 |
| `globalRecords` 等 | `AtomicLong` | 安全 |

### 5.3 人造分片的局限性

无 Translation 时按节点 ID 均匀分片，可能恰好切在环或长链中间。这会导致：
- **链断裂**：跨区线性链被拆为两个 block
- **环丢失**：跨区环在两个分片中各为半环，无一个检测到 cycle
- **边数据完整**：跨区边在源节点的 OCCT 记录中保留

修复需要后续的 `mergeOcctFiles()` 步骤（当前为 stub）。

---

## 6. 复杂度分析

| 阶段 | 时间 | 空间 |
|------|------|------|
| Phase 1: 加载 GBZ | O(文件大小) I/O | O(BWT + Seq) ≈ 5 GB |
| Phase 2: 出入度 | O(R × avg_sigma) ≈ 425M × 2 | O(N) BitSet ≈ 54 MB |
| Phase 3: DFS (单线程) | O(N + E + Runs) ≈ 213M + 5M + 78B | O(N) BitSet + O(D) stack |
| Phase 3: DFS (并行) | O((N+E+Runs)/threads) | 同上 × threads |

其中 N = 节点数，E = 总边数，Runs = 总游程数，D = DFS 最大深度。

### 6.1 每节点处理开销分解（实测估计）

| 操作 | 耗时 (μs) | 占比 |
|------|----------|------|
| ECC writer.write() | ~350 | 50% |
| BWT 解析（copyFromBlocks + RunCodec） | ~150 | 21% |
| 序列拷贝（sequences.copyWords） | ~80 | 11% |
| Block 判定 + BitSet | ~50 | 7% |
| rank_offset 计算 | ~40 | 6% |
| 栈操作 | ~30 | 4% |
