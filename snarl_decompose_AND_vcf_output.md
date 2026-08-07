## vg 的 Snarl 分解与 GBZ 到 VCF 的转换过程

### 1. 引言

vg (variant graph) 使用 **snarl** 作为变异位点的基本抽象单元。Snarl 分解将一个复杂的序列图（bi-directed sequence graph）转换成具有树形嵌套结构（snarl tree）和线性链接结构（chain）的层次化表示，这套结构与 VCF 中的变异位点（locus）和单倍型（haplotype）具有自然的对应关系。GBZ 文件 = GBWT（Graph Burrows-Wheeler Transform，单倍型索引）+ GBWTGraph（节点序列和边）+ tags，存储了完整的图拓扑和单倍型穿行信息。`vg deconstruct` 命令正是利用 snarl 分解从 GBZ 中恢复 VCF。

#### 1.1 关键术语

| 术语 | 含义 |
|------|------|
| **Snarl** | 图上的一个变异位点，由起始边界节点和终止边界节点定义，内部包含替代等位基因 |
| **Chain** | 一系列首尾相连的 snarl，共享边界节点 |
| **Ultrabubble** | 一种特殊的 snarl，满足无环、起始可达终止、无 turnaround 等性质 |
| **Unary snarl** | 起始节点与终止节点相同的 snarl（如倒位形成的自环） |
| **Cactus graph** | 仙人掌图：任意两个简单环至多共享一个顶点的图 |
| **Bridge forest** | 桥森林：合并所有环边后的图，形成森林结构 |
| **3-edge-connected component** | 3-边连通分量：删除任意两条边后仍连通的最大子图 |
| **NetGraph** | 将子 snarl/chain 视为不透明节点的抽象图，用于分析 snarl 内部连通性 |

---

### 2. Snarl 分解算法（`IntegratedSnarlFinder`）

vg 目前使用 `IntegratedSnarlFinder`（位于 [`src/integrated_snarl_finder.cpp`](/Users/wenjiepeng/Downloads/vg-3e1b3096828c5f758138d502a10f03ea57ad3e32/src/integrated_snarl_finder.cpp)）作为主要的 snarl 发现器。这是一种**纯 C++ 的内部实现**，无需依赖外部的 pinchesAndCacti C 库（旧版的 `CactusSnarlFinder` 使用它）。算法整体流程如下：

```
输入: HandleGraph G = (V, E)
输出: SnarlManager (包含 snarl tree 和 chain 结构)

Step 1: 计算邻接分量 (Adjacency Components)
Step 2: 计算 3-边连通分量，合并为 Cactus Graph
Step 3: 在 Cactus Graph 中找到所有简单环
Step 4: 合并环边，得到 Bridge Forest
Step 5: 在 Bridge Forest 中找最长路径，定根
Step 6: 遍历分解结构，生成 Snarl 和 Chain
```

#### 2.1 Step 1: 邻接分量 (Adjacency Components)

**数学定义**：对于双向序列图 $G = (V, E)$，定义图的每条边 $(u, v)$ 连接端点 $u$ 的右侧和 $v$ 的左侧（考虑链方向）。通过 Union-Find 结构，将每条边两侧的 handle 合并到同一个邻接分量中（[`integrated_snarl_finder.cpp:165-193`](/Users/wenjiepeng/Downloads/vg-3e1b3096828c5f758138d502a10f03ea57ad3e32/src/integrated_snarl_finder.cpp:165)）：

```
for each edge e = (u, d) in E:
    into_a = u  (inward-facing handle into u)
    into_b = flip(d)  (inward-facing handle into d)
    union_find.merge(into_a, into_b)
```

这等价于将线性路径上的节点收缩为一个超节点：如果节点 $x$ 的出边唯一且节点 $y$ 的入边唯一，且它们之间的边是唯一连接，则它们被划入同一分量。

**MergedAdjacencyGraph** 类实现了一个通用的分量合并图，每个分量由其 "head" handle 标识。关键操作包括：
- `merge(a, b)`: 合并两个分量
- `find(h)`: 查找 handle 所属分量的 head
- `for_each_member(head, f)`: 遍历分量内所有成员 handle

#### 2.2 Step 2: 3-边连通分量与 Cactus 图

**理论基础**：仙人掌图（cactus graph）是任意两个简单环最多共享一个顶点的图。序列图中，3-边连通分量（删除任意两条边后仍连通的部分）的收缩可以将一般图转化为仙人掌图。

vg 使用的 3-边连通分量算法来自 **Norouzi & Tsin (2014)** "A simple 3-edge connected component algorithm revisited"，实现在 [`src/algorithms/three_edge_connected_components.cpp`](/Users/wenjiepeng/Downloads/vg-3e1b3096828c5f758138d502a10f03ea57ad3e32/src/algorithms/three_edge_connected_components.cpp)中。

**算法核心——Absorb-Eject 操作**：

定义节点 $u$ 在边 $e = (u, v)$ 上的 **absorb** 操作：

1. 将 $v$ 的所有邻边（除 $e$ 外）转移给 $u$
2. 删除边 $e$
3. 如果 $v$ 的度变为 0 或 1，则将其移除（或保持在单独的 3-边连通分量中）

在深度优先搜索中执行 absorb-eject，通过维护节点的"有效度"（effective degree）来追踪图的拓扑变化，而无需实际修改图结构（[`integrated_snarl_finder.cpp:1130-1160`](/Users/wenjiepeng/Downloads/vg-3e1b3096828c5f758138d502a10f03ea57ad3e32/src/integrated_snarl_finder.cpp:1130)）：

```
algorithms::three_edge_connected_component_merges<handle_t>(
    /* emit all adjacency component heads as nodes */
    [&](emit_node) { cactus.for_each_head(emit_node); },
    /* emit edges: follow each member handle to its connected head */
    [&](node, emit_edge) {
        cactus.for_each_member(node, [&](other) {
            emit_edge(cactus.find(graph->flip(other)));
        });
    },
    /* record merges (deferred to avoid invalidating the algorithm) */
    [&](a, b) { merge_list.emplace_back(a, b); }
);
// Execute deferred merges
for (auto& ab : merge_list) { cactus.merge(ab.first, ab.second); }
```

算法完成后，`cactus` 中的每个分量代表仙人掌图中的一个顶点（即原图中的一个 3-边连通分量）。

#### 2.3 Step 3: 在仙人掌图中找环

使用 DFS 遍历仙人掌图来发现所有简单环（[`integrated_snarl_finder.cpp:272` `cycles_in_cactus()`](/Users/wenjiepeng/Downloads/vg-3e1b3096828c5f758138d502a10f03ea57ad3e32/src/integrated_snarl_finder.cpp:272)）：

```python
# 伪代码
for each component head (unvisited):
    DFS from component head:
        when a back edge is encountered:
            # back edge closes a cycle
            cycle = walk_stack_from_target_frame_to_current_frame()
            cycle_length = sum(edge_lengths on cycle)
            for each consecutive pair of edges on cycle:
                next_edge[edge_i] = edge_{i+1}
            if cycle_length > longest_cycle_length[component]:
                longest_cycle[component] = (cycle_length, first_edge)
```

**数据结构**：
- `next_along_cycle`: `unordered_map<handle_t, handle_t>`，映射环上每条入边到其顺时针方向的下一入边（仅在一个方向上存储）
- `longest_cycles`: `vector<pair<size_t, handle_t>>`，每个连通分量的最长环的（长度，环上某入边）

由于仙人掌图的性质，任意两个环最多共享一个顶点，这意味着 DFS 找到的每个 back edge 对应的环不会与其他环重叠边，因此一次遍历即可找到所有环。

#### 2.4 Step 4: 桥森林 (Bridge Forest)

复制一份仙人掌图作为桥森林（[`integrated_snarl_finder.cpp:1180`](/Users/wenjiepeng/Downloads/vg-3e1b3096828c5f758138d502a10f03ea57ad3e32/src/integrated_snarl_finder.cpp:1180)）：

```
MergedAdjacencyGraph forest(cactus);  // deep copy
for each cycle edge e:
    forest.merge(e, next_along_cycle[e]);  // 合并环上所有边
```

将所有环边合并后，仙人掌图变为一个森林（每个树对应原图的一个桥连通分量）。

**最长路径计算**（[`integrated_snarl_finder.cpp:547` `longest_paths_in_forest()`](/Users/wenjiepeng/Downloads/vg-3e1b3096828c5f758138d502a10f03ea57ad3e32/src/integrated_snarl_finder.cpp:547)）：

对于每个树，通过两次 BFS/DFS 找到树的直径（最长的 leaf-to-leaf 路径）：

1. 从任意节点出发，找到最远节点 $u$
2. 从 $u$ 出发，找到最远节点 $v$
3. 路径 $u \to v$ 即树的直径

同时计算 `towards_deepest_leaf`：对于森林中的每个分量 head，指向其最深的叶节点方向的边。这用于后续将"悬挂"在最长路径上的分支正确地定位到 chain 中。

**边长度**：每条边的权重 = $\text{节点序列长度} + \text{extra\_weight}(\text{node\_id})$（[`integrated_snarl_finder.cpp:158`](/Users/wenjiepeng/Downloads/vg-3e1b3096828c5f758138d502a10f03ea57ad3e32/src/integrated_snarl_finder.cpp:158)）。

#### 2.5 Step 5: 定根 (Rooting)

将最长环和最长路径按长度排序，选择全局最长的结构作为 snarl 树的根（[`integrated_snarl_finder.cpp:1228-1237`](/Users/wenjiepeng/Downloads/vg-3e1b3096828c5f758138d502a10f03ea57ad3e32/src/integrated_snarl_finder.cpp:1228)）：

```python
sort(longest_cycles)   # ascending by length
sort(longest_paths)    # ascending by length

while longest_cycles or longest_paths:
    if longest_path wins (or no cycle left):
        root_snarl = longest_path.pop()
        emit root chain along the path
    else:
        # longest cycle wins
        root_snarl = longest_cycle.pop()
        emit cyclic chain
```

**对参考基因组的定向**：`vg deconstruct` 在调用 `IntegratedSnarlFinder` 时，通过 `extra_node_weight` 参数给每个参考路径的首尾节点加上极大的额外权重（$10^{10}$，见 [`deconstruct_main.cpp:301-307`](/Users/wenjiepeng/Downloads/vg-3e1b3096828c5f758138d502a10f03ea57ad3e32/src/subcommand/deconstruct_main.cpp:301)），这保证了含有参考路径端点节点所在的路径/环在长度比较中"胜出"，从而 snarl 树总是沿参考基因组定向。

#### 2.6 Step 6: 遍历分解结构

`traverse_computed_decomposition()`（[`integrated_snarl_finder.cpp:1283`](/Users/wenjiepeng/Downloads/vg-3e1b3096828c5f758138d502a10f03ea57ad3e32/src/integrated_snarl_finder.cpp:1283)）执行一个双层的栈式遍历：

- **外层栈**：snarl 栈（以最外层虚拟 snarl 为根，无边界节点）
- **内层**：chain 栈（snarl 内的 chain 序列）

遍历发射事件：
```
begin_chain(handle)   -- 进入一个 chain
  begin_snarl(handle) -- 进入 chain 中的第一个 snarl
    ... 子 chain 和子 snarl 的递归遍历 ...
  end_snarl(handle)   -- 离开 snarl
  begin_snarl(handle) -- 进入 chain 中的下一个 snarl
  ...
end_chain(handle)     -- 离开 chain
```

对于环链（cyclic chain），`begin_chain` 和 `end_chain` 使用相同的 handle（环在物理上无头无尾，但遍历中需要一个逻辑起点）。

#### 2.7 Snarl 分类

对每个 snarl 进行类型判定（[`cactus_snarl_finder.cpp:195-260`](/Users/wenjiepeng/Downloads/vg-3e1b3096828c5f758138d502a10f03ea57ad3e32/src/cactus_snarl_finder.cpp:195)），通过构建 `NetGraph` 分析三种连通性：

- $R_{ss}$: **start-self-reachable** — 从 start 出发能否绕回 start 本身（反向）
- $R_{ee}$: **end-self-reachable** — 从 end（反向）出发能否到达 end（正向）
- $R_{se}$: **start-end-reachable** — 从 start 出发能否到达 end

分类逻辑：
```
if start.node_id == end.node_id:
    type = UNARY
elif not R_se:
    type = UNCLASSIFIED
elif R_ss or R_ee:
    type = UNCLASSIFIED  (存在 directed cycle)
elif all children are ULTRABUBBLE and net graph is directed acyclic:
    type = ULTRABUBBLE
else:
    type = UNCLASSIFIED
```

其中 **directed acyclic net graph（DAG）** 的判定通过 `handlealgs::is_directed_acyclic()` 在忽略子 snarl 内部连通的 NetGraph 上完成（仅考虑直接连接）。

#### 2.8 SnarlManager 的数据结构

`SnarlManager`（[`src/snarls.hpp`](/Users/wenjiepeng/Downloads/vg-3e1b3096828c5f758138d502a10f03ea57ad3e32/src/snarls.hpp)）维护以下核心索引：

| 索引 | 类型 | 用途 |
|------|------|------|
| `snarls` | `deque<SnarlRecord>` | 主存储，指针永不失效 |
| `children_of(snarl)` | `vector<const Snarl*>` | 子 snarl 列表 |
| `parent_of(snarl)` | `const Snarl*` | 父 snarl 指针 |
| `chain_of(snarl)` | `const Chain*` | 所属 chain 指针 |
| `snarl_into(node_id, direction)` | `const Snarl*` | 节点遍历方向到 snarl 的映射 |
| `roots` | `vector<const Snarl*>` | 顶层 snarl 集合 |

`SnarlRecord` 是一个巧妙的内存布局：它将 `Snarl`（protobuf 消息）和附加元数据（children, parent, chain, chain_index）连续存放，通过 `reinterpret_cast` 在 `Snarl*` 和 `SnarlRecord*` 之间转换（称为 "record" 惯用法）。

---

### 3. GBZ 到 VCF 的转换过程

#### 3.1 GBZ 文件结构

GBZ 是一个复合文件格式（[`src/gbzgraph.hpp`](/Users/wenjiepeng/Downloads/vg-3e1b3096828c5f758138d502a10f03ea57ad3e32/src/gbzgraph.hpp)）：

```
GBZ = GBWT (单倍型索引) + GBWTGraph (图拓扑 + 节点序列) + Tags (图名等元数据)
```

- **GBWT**（Graph BWT）：存储所有单倍型穿行（haplotype threads），每个 thread 是一系列 `node_type`（node id + orientation）的序列。支持高效的 `find(node)` 和 `extend(state, node)` 操作。
- **GBWTGraph**：从 GBWT 中提取的隐式图，实现 `PathHandleGraph` 接口，提供节点序列访问、路径遍历、步（step）的位置查询等功能。
- **Tags**：存储参考样本名等元数据。

`GBZGraph` 通过 `bdsg::PathHandleGraphProxy<gbwtgraph::GBWTGraph>` 将这一切包装成统一的 PathHandleGraph。

#### 3.2 `vg deconstruct` 命令流程

主入口在 [`src/subcommand/deconstruct_main.cpp`](/Users/wenjiepeng/Downloads/vg-3e1b3096828c5f758138d502a10f03ea57ad3e32/src/subcommand/deconstruct_main.cpp)：

```
输入: GBZ 文件 + 参考路径名 (可选)
输出: VCF (stdout)

Step 1: 加载 GBZ → GBZGraph (PathHandleGraph + GBWT)
Step 2: 用 IntegratedSnarlFinder 计算 snarls (或从文件加载)
Step 3: 扫描路径元数据，识别参考路径和样本单倍型
Step 4: 对每个 snarl，提取穿行 → 映射到等位基因 → 计算基因型 → 输出 VCF 记录
Step 5: 排序并写出最终 VCF
```

#### 3.3 Step 2: Snarl 计算的偏置

`vg deconstruct` 在构造 `IntegratedSnarlFinder` 时传递 `extra_node_weight`（[`deconstruct_main.cpp:301-307`](/Users/wenjiepeng/Downloads/vg-3e1b3096828c5f758138d502a10f03ea57ad3e32/src/subcommand/deconstruct_main.cpp:301)）：

```cpp
for each reference path:
    extra_node_weight[path_first_node] += 10^10;
    extra_node_weight[path_last_node]  += 10^10;
```

这确保了 snarl 分解以参考路径首尾节点所在的桥森林路径/环为根，使得 snarl 的边界节点与参考基因组的坐标对齐。

#### 3.4 Step 3: 样本和参考路径识别

`Deconstructor::deconstruct()`（[`deconstructor.cpp`](/Users/wenjiepeng/Downloads/vg-3e1b3096828c5f758138d502a10f03ea57ad3e32/src/deconstructor.cpp)）扫描图中所有路径：

- **参考路径**：由 `-p` 或 `-P` 指定的路径。默认所有 `REFERENCE` 和 `GENERIC` 路径均为参考。
- **样本单倍型**：`HAPLOTYPE` 路径（来自 GBWT threads 或嵌入路径）。每个样本的倍性由 phase 范围推断。

#### 3.5 Step 4: 对每个 Snarl 的 VCF 记录生成

核心函数 `deconstruct_site()`（[`deconstructor.cpp:636`](/Users/wenjiepeng/Downloads/vg-3e1b3096828c5f758138d502a10f03ea57ad3e32/src/deconstructor.cpp:636)）：

##### 3.5.1 获取穿行 (Traversals)

`get_traversals()` 合并两类穿行：

1. **路径穿行**（`PathTraversalFinder`）：从嵌入路径中提取穿过 snarl 的子路径。每个穿行携带：
   - `Traversal`：handle 序列
   - `path_name`：来源路径名
   - `step_handle` 对：穿行在参考路径上的起止步，用于确定坐标

2. **GBWT 穿行**（`GBWTTraversalFinder`，[`traversal_finder.cpp:3441`](/Users/wenjiepeng/Downloads/vg-3e1b3096828c5f758138d502a10f03ea57ad3e32/src/traversal_finder.cpp:3441)）：
   通过 `list_haplotypes()`（[`haplotype_extracter.cpp:183`](/Users/wenjiepeng/Downloads/vg-3e1b3096828c5f758138d502a10f03ea57ad3e32/src/haplotype_extracter.cpp:183)）执行 BFS：

```
从 snarl_start 出发，使用 GBWT 进行广度优先搜索：
    search_intermediates = [(thread=[start_node], state=gbwt.find(start_node))]
    while search_intermediates:
        pop last thread + state
        for each outgoing edge from thread's last node:
            new_state = gbwt.extend(state, next_node)
            if new_state is not empty:    # GBWT 中存在此延伸
                if next_node == snarl_end:
                    emit result
                else:
                    push to search_intermediates
```

这种搜索只沿 GBWT 中实际存在的单倍型穿行扩展，因此效率远高于穷举搜索。结果去重后，每个穿行携带对应 GBWT path identifiers。

##### 3.5.2 参考穿行的选择

对于有多个参考穿行的 snarl（因参考路径在该位点形成环），通过 **Jaccard 系数** 将替代等位基因分配至正确的参考位置：

$$
J(A, B) = \frac{|A \cap B|}{|A \cup B|}
$$

其中 $A$、$B$ 是参考路径和替代路径在 snarl 两侧窗口内的节点 ID 集合（窗口大小默认为 10000 bp，`-c` 参数可调）。

##### 3.5.3 等位基因字符串提取

对于每个穿行 $\vec{t} = (h_1, h_2, \dots, h_k)$（跳过首尾边界节点 $h_1, h_k$），等位基因字符串为：

$$
\text{allele\_string} = \text{concat}\left(\text{sequence}(h_2), \text{sequence}(h_3), \dots, \text{sequence}(h_{k-1})\right)
$$

参考穿行对应 `REF` 等位基因（索引 0），其他对应 `ALT` 等位基因（索引 1, 2, ...）。重复的序列映射到同一等位基因。

##### 3.5.4 穿行聚类

（可选，通过 `-L` 启用）使用 handle Jaccard 系数对穿行进行聚类：

$$
J(\vec{t}_i, \vec{t}_j) = \frac{|\text{handles}(\vec{t}_i) \cap \text{handles}(\vec{t}_j)|}{|\text{handles}(\vec{t}_i) \cup \text{handles}(\vec{t}_j)|}
$$

相似度 $\ge$ 阈值的穿行合并到同一簇，同一簇的所有穿行共享一个等位基因。`cluster_min_allele_len` 参数可限制此聚类仅作用于 SV 位点。

##### 3.5.5 基因型推断

对于每个样本，`choose_traversals()` 将所有穿行映射到等位基因后，按以下规则选择基因型：

1. 统计各等位基因在该样本穿行中的出现频率
2. 按频率降序排列，优先选择 ALT 等位基因（ref allele 优先被替代）
3. 根据样本倍性（从 GBWT phase 范围推断）选择 top-k 穿行
4. 若无冲突（同一 phase 的多个穿行映射到不同等位基因），输出 `allele1|allele2|...`
5. 若有冲突，输出 `.|.` 并标记 `CONFLICT` INFO 字段

##### 3.5.6 VCF 记录坐标

snarl 的 VCF 坐标由参考穿行在参考路径上的位置确定：

```
v.position = reference_path_position_of_snarl_boundary_node + reference_path_subrange_offset + 1
```

对于 indel，`v.position` 额外左移一位（VCF 规范要求 indel 记录包含前一个碱基）。

##### 3.5.7 Star Allele 支持

当使用 `-a -R`（嵌套模式 + star allele）时，采用自顶向下的处理（[`deconstruct_graph_top_down()`](/Users/wenjiepeng/Downloads/vg-3e1b3096828c5f758138d502a10f03ea57ad3e32/src/deconstructor.cpp:1350)）：
- 父 snarl 处理完后，将样本在父 snarl 中的所有单倍型传递给子 snarl
- 对于在父 snarl 中存在但在当前子 snarl 中无穿行的样本单倍型，生成 `*` 等位基因（表示"父位点的变异跳过了子位点"）
- 这对应 VCF 中的 overlapping deletion 场景

---

### 4. 两者的内在联系

Snarl 分解与 GBZ→VCF 转换的关系可以概括为一个映射管道：

```
GBZ (图 + 单倍型)
    │
    ▼
IntegratedSnarlFinder  ─── SnarlManager (变异位点树)
    │                              │
    │                     snarl boundaries → reference coordinates
    │                     snarl traversals → alleles
    │                     GBWT threads → genotypes
    ▼
Deconstructor ─── VCF 输出
```

**详细的对应关系**：

| 图论概念 | 数学/算法对应 | VCF 对应 |
|----------|-------------|----------|
| Snarl $(s, t)$ | 仙人掌图中连接桥边路径或环边对的子图 | 一个 VCF 记录（变异位点） |
| Snarl 边界节点 | 桥森林中连接两个分量的唯一边 | 参考位置坐标 `POS` |
| Snarl 内部穿行 | GBWT BFS 从 `s` 到 `t` 的所有单倍型路径 | `REF` / `ALT` 等位基因 |
| Chain $(S_1, S_2, \dots, S_k)$ | 桥森林中的路径/环，连接一系列 snarl | 线性排列的变异位点 |
| 嵌套 Snarl | 仙人掌图中环与环、环与桥的层叠关系 | `LV`（level）、`PS`（parent snarl ID）标签 |
| GBWT thread | 图上的完整单倍型穿行 | 样本基因型 GT 字段 |
| 额外权重 `extra_node_weight` | 给参考路径端点的偏置权重 $10^{10}$ | 确保 snarl 树沿参考方向定向 |

**关键洞察**：

1. **Snarl 分解确定"哪里是变异位点"**：仙人掌图的数学结构保证了 snarl 分解的良定性——每个 3-边连通分量是一个 snarl 的"原子"内容，环决定 chain 的结构，桥决定嵌套层序。

2. **GBWT 提供"谁有什么变异"**：GBWT 存储了所有单倍型穿行，`GBWTTraversalFinder` 只是在每个 snarl 内对 GBWT 做局部 BFS，提取过该位点的那部分穿行。

3. **参考基因组定向至关重要**：`extra_node_weight` 机制将 VCF 的参考坐标系（contig + position）与图的拓扑结构（桥森林中的最长路径）耦合起来，使得 snarl 边界在参考路径上有明确的位置。

4. **复合格式 GBZ 使这个流程自包含**：GBZ = 图 + 单倍型，因此 `vg deconstruct` 可以仅从一个 GBZ 文件出发，完成从图拓扑分析到 VCF 输出的全流程，无需额外的 GBWT 或参考基因组文件。

---

### 5. 代码结构索引

| 文件 | 内容 |
|------|------|
| [`src/integrated_snarl_finder.hpp`](/Users/wenjiepeng/Downloads/vg-3e1b3096828c5f758138d502a10f03ea57ad3e32/src/integrated_snarl_finder.hpp) | IntegratedSnarlFinder 接口 |
| [`src/integrated_snarl_finder.cpp`](/Users/wenjiepeng/Downloads/vg-3e1b3096828c5f758138d502a10f03ea57ad3e32/src/integrated_snarl_finder.cpp) | Snarl 发现主算法（MergedAdjacencyGraph、环检测、桥森林、遍历） |
| [`src/algorithms/three_edge_connected_components.cpp`](/Users/wenjiepeng/Downloads/vg-3e1b3096828c5f758138d502a10f03ea57ad3e32/src/algorithms/three_edge_connected_components.cpp) | Tsin 3-边连通分量算法 |
| [`src/snarls.hpp`](/Users/wenjiepeng/Downloads/vg-3e1b3096828c5f758138d502a10f03ea57ad3e32/src/snarls.hpp) | Snarl, Chain, SnarlManager, NetGraph 数据结构 |
| [`src/cactus_snarl_finder.cpp`](/Users/wenjiepeng/Downloads/vg-3e1b3096828c5f758138d502a10f03ea57ad3e32/src/cactus_snarl_finder.cpp) | 基于外部 C 库的 CactusSnarlFinder（含 snarl 分类逻辑） |
| [`src/deconstructor.hpp`](/Users/wenjiepeng/Downloads/vg-3e1b3096828c5f758138d502a10f03ea57ad3e32/src/deconstructor.hpp) | Deconstructor 接口 |
| [`src/deconstructor.cpp`](/Users/wenjiepeng/Downloads/vg-3e1b3096828c5f758138d502a10f03ea57ad3e32/src/deconstructor.cpp) | 穿行提取、等位基因生成、基因型推断 |
| [`src/subcommand/deconstruct_main.cpp`](/Users/wenjiepeng/Downloads/vg-3e1b3096828c5f758138d502a10f03ea57ad3e32/src/subcommand/deconstruct_main.cpp) | `vg deconstruct` 命令行入口 |
| [`src/traversal_finder.hpp`](/Users/wenjiepeng/Downloads/vg-3e1b3096828c5f758138d502a10f03ea57ad3e32/src/traversal_finder.hpp) | 各种 TraversalFinder 接口定义 |
| [`src/traversal_finder.cpp`](/Users/wenjiepeng/Downloads/vg-3e1b3096828c5f758138d502a10f03ea57ad3e32/src/traversal_finder.cpp) | GBWTTraversalFinder 实现（GBWT BFS 穿行提取） |
| [`src/haplotype_extracter.cpp`](/Users/wenjiepeng/Downloads/vg-3e1b3096828c5f758138d502a10f03ea57ad3e32/src/haplotype_extracter.cpp) | `list_haplotypes()` 函数 |
| [`src/gbzgraph.hpp`](/Users/wenjiepeng/Downloads/vg-3e1b3096828c5f758138d502a10f03ea57ad3e32/src/gbzgraph.hpp) | GBZGraph 包装类 |

---

### 6. 关键数学定义总结

**定义 1 (Snarl)**：在有向双向图 $G$ 中，snarl 是一个三元组 $(s, t, G'_{st})$，其中 $s, t$ 是边界节点，$G'_{st}$ 是 $s$ 和 $t$ 之间的内部子图，满足 $s$ 可达 $t$（对于 ultrabubble）。

**定义 2 (Chain)**：Chain 是 snarl 的有序序列 $C = (S_1, S_2, \dots, S_k)$，其中 $S_i$ 的终止边界与 $S_{i+1}$ 的起始边界共享同一节点。每个 snarl 可带有方向（forward/backward）。

**定义 3 (仙人掌图)**：图 $C$ 是仙人掌图当且仅当 $C$ 的任意两个简单环至多共享一个顶点。等价定义：$C$ 的每个块（block）要么是一条边，要么是一个简单环。

**定义 4 (3-边连通分量)**：图 $G$ 的一个极大子图 $H$，满足对于 $H$ 中任意两个顶点 $u, v$，在 $H$ 中至少有 3 条边不相交的 $u$-$v$ 路径。

**定理 (Tsin 2007 / Norouzi-Tsin 2014)**：任意无向图可以通过合并其 3-边连通分量转化为唯一的仙人掌图。该转化可在 $O(n + m)$ 时间内完成。

**引理 (Snarl 分解的良定性)**：图 $G$ 的仙人掌图中，每个环对应一个 chain，每个桥对应一个嵌套层。桥森林的最长路径/环选择给出了 snarl 树的定根方案。

**命题 (GBWT 局部 BFS 的正确性)**：给定 snarl $(s, t)$ 和 GBWT 索引，从 $s$ 出发的 `gbwt.extend()` BFS 恰好产生该 snarl 内所有 GBWT 中存在的单倍型穿行。（证明：GBWT 的 `extend` 操作等价于 FM-index 上的后缀数组查询，BFS 保证了完备性。）

### 7. 总结

vg 的 snarl 分解通过**仙人掌图 → 桥森林 → 最长路径定根**的算法流水线，将复杂的序列图转化为层次化的嵌套位点表示。GBZ→VCF 转换在此基础上，利用 GBWT 的 BFS 在 snarl 边界间提取单倍型穿行，将其映射为等位基因和基因型。两类算法的核心连接点在于：**snarl 定义了"哪里"是变异位点（用图拓扑），GBWT 提供了"谁"有什么变异（用压缩的单倍型索引），而参考基因组偏置保证了输出 VCF 的坐标体系与标准参考对齐**。







