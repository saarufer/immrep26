# 从 TCR 序列与 peptide–HLA 到识别谱，再到功能激活

**一个机制驱动、可证伪、可复现的 TCR–pMHC 研究计划（2026-10-26，v1.0）**

> **核心路线**  
> `TCR sequence + peptide–HLA → Recognition profile → Functional activation`
>
> 研究目标不是再为 IMMREP26 追加一个相似度特征，而是建立两级可验证的生物学描述：**哪些 peptide 能被某个 TCR 识别？** 以及 **这种识别在给定呈递与细胞条件下，什么时候足以引发功能？**

**项目定位**：机制研究 + 下一代预测模型的可行性验证；可以服务 IRIS 等后续项目。**IMMREP26 当前提交 `density3_struct03.csv` 保持冻结，不在此项目中重新搜索公开表位的融合权重、提交新预测或推算私榜 AUC。**

**证据说明**：本计划的公开资源、字段和论文链接截至 2026-10-26 已对照作者/数据仓库说明；**尚未实际下载并分析原始 SEEKER/PRP CSV**。具体列名、过滤阈值、实验重叠度、可建模样本量应在 Gate 0 以文件审计结果为准。文中的样本量、时间、模型和阈值建议属于**预注册的分析方案**，不是已获得的实验结果。

---

## 0. 核心科学假设及判别标准

### H0｜为什么“序列距离”不是“功能距离”

相似 TCR 可以识别不同 peptide；相距较远的 TCR 也可能共享识别谱。序列密度可反映重排、抗原选择、重复发现和研究来源的组合，不保证刻画单次分子识别的物理机制。

**关键区别**：

1. **Sequence similarity**：两条 TCR 的 CDR/V 序列有多近。
2. **Recognition similarity**：在**相同 HLA 与可比 peptide 面板**上，它们的结合/识别读出有多相似。
3. **Functional similarity**：在**相同细胞、表达量、pMHC 呈递、剂量、时间**条件下，它们的激活曲线有多相似。

前两者不等于第三者；后两者也不应通过同一个未经校准的阈值直接转换。

### H1｜识别谱的表示能否跨 TCR 泛化？

- **待检验**：使用 peptide 识别谱学习到的表征，在**未见过的 TCR 家族**上，比单纯 TCR 序列相似度更能解释识别差异。
- **反证条件**：所谓增益只存在随机行划分/同一 TCR 内部、消失于家族留出，或依赖跨面板缺失值被错误编码为零。
- **终点**：跨 TCR（严格家族留出）的连续结合/富集预测、功能谱距离的可复现性；先报告置信区间，不假定一定改善。

### H2｜识别到激活存在条件依赖的非线性映射吗？

- **待检验**：结合/筛选信号与 NFAT/ERK/细胞功能之间可能存在阈值、饱和、时间与剂量依赖；变化也可能是激活细胞比例而非单细胞强度。
- **反证条件**：控制 pMHC 呈递、TCR 表达、细胞背景后，阈值模型相对简约连续模型没有独立预测优势。
- **注意**：公开数据若只有二元结果或选中的 hit，**不能识别** Hill slope、EC50、动力学校对步数或严格的“全或无”机制。

### H3｜跨 peptide 的突变容忍谱是否比通用替换矩阵更接近识别生物学？

- **待检验**：相同 peptide 变体面板上的位点效应、互作效应及可替代氨基酸群，是否解释“序列相似但识别不同”的实例。
- **反证条件**：性能只来自 peptide anchor 对 HLA 稳定性的影响，不能分离 TCR 特异性接触信号；在独立 peptide 家族上失效。

### H4｜模型能否表达“我不知道”？

对远离训练 TCR / peptide / HLA 的输入应扩大不确定性，拒绝过度自信，而不是通过虚构阴性提高表面 AUC。

---

## 1. 数据层：三种不同实验读出，禁止混标签

| 数据源 | 已核对的具体内容 | 真实读出 | 主要用途 | 关键限制 |
|---|---|---|---|---|
| **PRP / Wang et al., Nature Biotechnology 2026** | Dryad 约 **183.88 MB**，2 个 naïve library + 21 个 TCR 选择文件；包括 **16 个患者 TCR 及 5 个工程变体**，HLA-B*27:05 | Yeast-display peptide 选择前/后各轮 read count，**偏向结合富集** | 学习并比较 TCR 的 peptide recognition profile | 不等于 T 细胞激活；P2/P8 与 P2/P9 固定锚点面板不同，工程变体不能跨亲本泄漏 |
| **SEEKER / Li et al., bioRxiv + Dryad 2026** | 约 **27.95 MB**，TCR59/TCR63 的 HLA-B*27:05 筛选含 R0–R2 与部分 flow；四个参考 TCR 对同一 9-mer 文库的 top-10,000；另含患者 PacBio 联合筛选与 101-TCR 参考表 | Jurkat-derived NFAT reporter 筛选富集；**部分有单/高拷贝 flow 激活比例** | 与功能层关联；验证 peptide 变体响应及特异性 | 四受体 top-10,000 是**截断后的正向候选集**，未上榜≠阴性；筛选富集不等于直接测得的结合常数 |
| **IMMREP26 历史复验集** | 2,636 αβ TCR、8 表位；5 Public 的明确标签，3 Private 无标签；252 Public ambiguous | 特定实验协议下的**二元功能复验** | 分析外部一致性与任务差异 | 非连续效力测量；历史标签与 SEEKER/PRP 的测量体系不同，不能直接合并监督 |
| **Hi-TPH / Zhu et al., Nature Computational Science 2026（可选）** | 分层 TCR–pHLA 数据及不同 nonbinder 构造策略 | 异构的结合/非结合证据 | 外部基线与负样本审计 | **构造的 nonbinder ≠ 实测功能阴性**；官方声明 CC BY-NC 4.0，科研外使用需单独审查 |

### 来源与固定版本

- PRP 论文：<https://doi.org/10.1038/s41587-026-03128-x>
- PRP 数据：<https://doi.org/10.5061/dryad.h70rxwdzk>
- PRP 代码：<https://github.com/akds/PRP-TCR-Specificity>；作者仓库包含部分 activation 数据/模型配置，**实际可链接的 TCR–peptide 对应关系须先核实**。
- SEEKER 论文：<https://doi.org/10.64898/2026.09.25.753432>（预印本）。
- SEEKER 数据：<https://doi.org/10.5061/dryad.79cnp5jct>；附带 README 数据字典。
- SEEKER 代码：<https://github.com/ldi2009/SEEKER-2026>；数据页标注相应提交 `3f4d822f541669f5e46eb56c182563ef320bdee0`。
- Hi-TPH：<https://doi.org/10.1038/s43588-026-01073-y>；<https://github.com/Jiadong001/Hi-TpH>。

**不可跨越的证据界限**：

- PRP 主要测量**体外结合/选择富集**，SEEKER 测量**功能 reporter 选择富集**。不能为了拼出两级监督而把 SEEKER 富集当作 binding ground truth。
- SEEKER 的四个参考 TCR 平行筛选可比较共享文库中的**候选排名**，但只发布 top-10,000 时，不能把未出现的 peptide 列当零或负例；比较集应保留审查到的检测/覆盖信息。
- HLA-B*27:05 的 PRP 结果不能直接推及 IMMREP26 的 A*01:01/A*02:01 或 IRIS 的 A*02:01/A*11:01。跨 HLA 模型属于后续阶段，必须有跨 HLA 实测监督。
- 同一研究的工程 TCR（如 PRP 19.2 亲本与 C1–C5）不能当作互相独立的测试 TCR；患者与序列家族也应成组划分。

---

## 2. 统一的机器可读数据协议（Gate 0）

### 2.1 目录与文件元信息

建立 `raw_manifest.tsv`，逐文件记录 `source`, `doi`, `source_url`, `version_date`, `filename`, `sha256`, `license`, `size_bytes`, `n_rows`, `n_unique_tcr`, `n_unique_peptides`, `hla`, `assay`, `rounds`, `notes`。不覆盖原始数据。

### 2.2 统一长表 `recognition_observation.parquet`

每一行 = **一个受体 × peptide × 实验条件 × 重复 × 轮次/时间点** 的实测或经处理的一个读出。

```text
observation_id
source_dataset           # PRP, SEEKER_B27, SEEKER_4TCR, SEEKER_PACBIO, OTHER
study_id                 # 原始研究/实验批次，不等于 publication DOI
assay_type               # yeast_binding_enrichment / NFAT_screen / flow_activation / ...
tcr_id, parent_tcr_id, tcr_family_id
tra_cdr3, trb_cdr3, trav, traj, trbv, trbj
peptide_aa, peptide_length, hla_allele
library_id, library_anchor_scheme, peptide_present_in_input
condition_id, replicate_id, selection_round, stimulation_time
raw_count, library_total_count, input_count
readout_value, readout_unit, readout_uncertainty
label_status             # measured / measured_below_detection / not_observed / not_tested
expression_control, pmhc_presentation_control, assay_metadata
source_file, source_row_id
```

**规则：** `not_observed`（筛选未检出）、`not_tested`（没有实验）和 `measured_below_detection`（实验测量且低于检测限）分别编码。任何下游处理都不能把三者转换成同一个 0。

### 2.3 数据 QC

- 精确标准化 HLA 命名、氨基酸序列合法性、CDR3 配对/唯一性；记录多来源复用的序列。
- PRP：核对两个 anchor library 的 peptide 位点约束；仅在共同可观测范围内做跨 TCR 对比。
- SEEKER：核对 R0 是否真的存在、候选是否受 top-k 截断、flow 是否只验证高排名 hits、单/高拷贝读出是否可比。
- 分别审计 exact peptide overlap、near-neighbor overlap 和 TCR overlap；汇总不同库之间有多少可用于 paired binding–activation calibration 的 **同 TCR / 同 peptide / 同 HLA** 样本。
- 整理来源协议中阴性、空白/无关 peptide、TCR-null 等真实实验对照（仅在原始来源确有该项时使用）。
- 导出 `quality_report.md`、`data_dictionary.md`、`overlap_matrix.csv`；不得因缺失交叉信息伪造映射。

### Gate 0：是否有足够的成对读出？

**可继续：** 有可比较的同受体同肽数据及独立重复/对照，或有明确的 activation validation 子集。  
**只做识别谱、不拟合 binding→activation：** 两个数据源完全没有可对齐的成对读出，或 paired 数据严重选择偏差、无法估计不确定性。  
**必须停止任何跨实验数值统一：** 只有 top-k 和缺失标签，却需要假设其余均为阴性才能推进。

---

## 3. Stage I — 构建定量 Recognition Profile

### 3.1 什么是 profile？

对于相同 HLA、相同 peptide panel 的受体 `t`，定义：

\[
R_t = (r_{t,p_1}, \ldots, r_{t,p_M}), \qquad M_{t,p}\in\{0,1\}
\]

- `r` 是所选读出的**经校正实测富集或模型预测**；不能无条件解释为亲和力。
- `M` 是**真实可比较**的观测掩码；缺失不等于零。
- 附带 `assay`, `hla`, `panel`, `uncertainty`。不同 assay 不直接合成一个向量。
- 如果采用绝对值刻度，需要对输入库频率/测序深度做校正；否则使用每受体内部排名或标准化效果并报告丢失的跨 TCR 幅度信息。

### 3.2 PRP 分析：优先从原始 count 而非最终 0/1 标签入手

**独立分开分析 P2/P8 和 P2/P9 两种锚点库。** 对有足够输入 counts 的 peptide，计算例如：

\[
r_{t,p}^{(r)} = \log_2\left(\frac{c_{t,p,r}+\alpha}{N_{t,r}+\alpha K_r}\right)
-\log_2\left(\frac{c_{p,0}+\alpha}{N_0+\alpha K_0}\right)
\]

其中 `c` 是 count，`N` 为对应样本总计数，`K` 为该轮候选序列数，`α` 为预先固定的伪计数。**这是候选的输入校正富集量，不是论文既定效应估计量**；必须与论文推荐的数据预处理和标签定义并列复现。

注意不同 selection round 会引入饱和、瓶颈和 PCR 放大偏差；优先检查 R1→R4 的一致性，不直接把最终轮次 abundance 解释为结合强度。零输入且只在后轮出现的 peptide 应单独标识，而不是计算无限富集。

**基础输出**：

1. 每条 TCR 的 peptide 覆盖、富集值分布、round trajectory、重复一致性/不确定性。
2. 每条 TCR 的位点特异 amino-acid preference（PWM）、pairwise position interaction、anchor vs putative TCR-facing positions 的对照。
3. 同库 TCR 的共同 peptide 观测数量及可计算功能谱相似度的有效样本量。
4. 对亲本和工程变体单独做 paired perturbation 分析（**不能算独立 TCR 数量**）。

### 3.3 SEEKER 分析：将 NFAT 筛选与 flow 验证分开

- `TCR59/TCR63`：优先使用有 R0 的筛选记录，计算 R1/R0、R2/R1，保留输入 coverage；对具有 `add_flow_single_copy` / `add_flow_high_copy` 的条目，分别计算与 flow activated fraction 的关系。
- `4TCR_screen`：共享 NNK 9-mer 库的四受体 top-10,000，先做 top-k **交集/排名稳定性**、peptide motif 对照和位点偏好；**不可**由未入榜构造 negative label。
- `PacBio`：αβ 配对 × peptide 轮次富集可以检验共同识别谱，但没有明确 R0 的项目不按输入校正强度解释；先审计有多少唯一 clonotype、来源关联与交叉覆盖。
- 每个 assay 的 readout 单独保留；只在 flow validation 条件可比时比较绝对激活比例。

### 3.4 Functional distance vs sequence distance：第一个关键判别实验

同一面板中选择共同观测质量合格的 peptide 集合 \(\mathcal P_{ij}\)，比较：

\[
d_{\mathrm{profile}}(i,j)=1-\operatorname{corr}_{p\in\mathcal P_{ij}}(r_{i,p},r_{j,p})
\]

并行使用 Spearman、余弦、Jensen–Shannon（**仅当已合法转换成归一化、非负分布**）以及加权距离。比较对象为 `tcrdist3`、SCEPTR、CDR3β-only、V-gene-conditional 及对照 random similarity。

**重要**：如果两条 TCR 的共同 peptide 很少，或其 peptide 主要是 selection hits，距离估计会偏；设预先确定的最小 overlap / 质量门槛，并报告支持矩阵与置信区间。近亲工程变体组和独立患者组分别报告。

**不允许的 shortcut**：把 PRP 论文中已用 test labels 微调得到的模型输出重新当成完全独立的生物学标签。

### 3.5 Gate 1

- 观察矩阵在共同 panel 中具有足够有效重叠，且谱距离对预处理、读数阈值不高度敏感；
- 能找到明确的**序列近但识别远**、**序列远但识别近**的高质量实例，最好有独立 round / flow 支持；
- 若高质量重叠不足，则只做每受体识别 motif 和扰动分析，不强行发表“跨 TCR 功能空间”。

---

## 4. Stage II — 从 sequence + peptide–HLA 预测 Recognition Profile

这是对 Stage I 的 **out-of-sample 预测**，不是把实测 profile 直接当成输入泄漏给测试 TCR。

### 4.1 任务形式

\[
\hat r_{t,p}=f_\theta(\mathrm{TRA},\mathrm{TRB},\mathrm{V/J},p,\mathrm{HLA},\mathrm{assay\ context})
\]

**适用域必须声明**：在 PRP 样本里 HLA 基本恒为 B*27:05，**无法单靠该数据训练或检验 HLA 泛化**。SEEKER 也不提供足够的任意 HLA 负样本来克服此限制。跨 HLA 只能在外部已测实验域验证。

### 4.2 三层模型，按照复杂度递进

| 层级 | 模型 | 学习的性质 | 作用 |
|---|---|---|---|
| M0 简单基线 | 每受体 PWM；peptide position additive model；实验背景频率；anchor-only HLA propensity | 位点独立 / 基础呈递 | 判断昂贵模型是否必要 |
| M1 交互 | pairwise amino-acid interactions、低秩 peptide feature / TCR embedding 双线性模型、稀疏正则化 | 识别位点间相互作用 | 检验非加性是否真的存在 |
| M2 条件模型（可选） | PRP 作者模型复现；ESM2/蛋白语言模型特征 + cross interaction，冻结大部分参数并限制容量 | 非线性、跨 TCR 共享 | **只有 M1 严格失败且数据足够才启用** |

HLA 呈递要作为独立基线或协变量控制，而非与 TCR-specific response 混成同一个不可解释的分数。

### 4.3 预注册的数据划分与泄漏阻断

至少报告四种严格程度：

1. **Within-TCR held-out peptide**：同一 TCR 下，留下未用于训练的 peptide；测试学习变体规律，但不能声称泛化到新 TCR。
2. **Peptide family/scaffold holdout**：按 Hamming/序列簇或明确的原始家族分组留出，防止同一近邻变体泄漏；保留 anchor 约束。
3. **Leave-TCR-family-out**：把亲本/工程变体与共享 CDR3/V 基因近邻家族绑定，整组留出；**这才是新 TCR 泛化的关键终点**。
4. **Joint TCR-family × peptide-family holdout**：最难的跨两维外推；若缺少足够独立受体，就标记为 exploratory/underpowered，而不是强报性能。

划分规则**在训练和任何模型选择前冻结**；同一 peptide 在不同实验批次的重复，不能同时进入 train 与 test。对工程衍生 TCR 实施亲本家族强制绑定。

### 4.4 性能与成功标准

- 若有可靠连续 readout：主要报告 Spearman、相对未校准基线的增益、误差与 calibration；按独立 TCR/家族汇总，不按数百万 peptide 行虚增样本量。
- 对经实验定义的阳/阴：AUROC + AUPRC，并附 prevalence、覆盖、阴性构造规则。
- 独立对照：`presentation-only`、`TCR-only`、`shuffled cognate pair`、位置频率匹配 decoy、固定预测器。
- 初始成功门槛（预注册建议）：在 **家族级留出** 的主终点，模型相对最强 M0/M1 基线有稳健正增益，区间/重复划分支持，且没有显著恶化 out-of-domain calibration；不预设任意绝对 AUROC 门槛。
- 如果从 PRP 的 16 个患者受体无法得到足够独立家族进行验证，必须诚实改为 proof-of-concept，不宣称可泛化任意 αβ TCR。

---

## 5. Stage III — Recognition → Functional activation：条件非线性映射

### 5.1 先明确测量条件

理论上功能反应依赖：

\[
A_{t,p,c,\tau}=g(r_{t,p},\,\mathrm{pMHC\ density},\,\mathrm{TCR\ expression},\,c,\,\tau,\,\mathrm{cell\ state})
\]

`c` 代表刺激剂量，`τ` 代表时间。`r` 可以是结合富集代理，但**不能假设其与 KD、koff 存在线性或单调的一一对应**。pMHC 表达/稳定性、共刺激与细胞状态也是潜在主要影响因素。

**目前公开数据识别边界**：

- SEEKER 的 NFAT 富集与流式 activated fraction 属于**同一功能通路中的筛选与验证**，可以检验从 screening readout 到 flow readout 的可迁移性，但本身不足以证明 **物理结合 → 细胞激活**。
- PRP 的 yeast binding enrichment 与作者仓库提供的激活验证资料，**只有存在可核实的同 TCR/peptide/HLA 对应**，才允许训练真正的两级映射。
- 若没有多剂量、多时间点资料，就**不能**声称估计了 EC50、Hill slope、kinetic proofreading 参数或数字式/模拟式单细胞机制。

### 5.2 两类可立即完成的计算任务

**F1：SEEKER 筛选到 flow 的校准**

- 在 TCR59/63 有真实 flow 的候选中，比较 R1/R0、R2/R1、round trajectory 对 `add_flow_single_copy` 的相关性与校准。
- 单拷贝和高拷贝条件**分开拟合并直接对照**，以显示抗原表达量可能改变功能标签。
- 由于 flow 被挑选测量，记录其条件抽样机制；**仅在测量子集上评估，不推断未测 peptide 的总体激活率**。
- 用简约的单调非参数校准（isotonic）和 logistic 作为基线；避免过度拟合 4 参数 Hill 曲线。

**F2：PRP 与外部 activation 的可链接配对**

- 检查作者仓库 `data/activation_binary.csv` 等实测激活资料（作者 README 提及），建立严格的 `tcr_id + peptide + HLA + assay` join。
- 只对真正 matched 条目比较 binding-enrichment 模型、presentation-only 模型、binding+presentation 的结果。
- 若只有 binary activation，报告校准与 ranking；阈值可作为预测函数的统计形状，**不解释为单细胞确定的生物学开关**。
- 将可链接样本数、每个 TCR 的阳/阴数、覆盖率、验证样本选择流程列成审计表；不足则停止 F2。

### 5.3 未来 wetlab：真正检验“全或无”的最低可行设计

这是**研究建议，不代表目前已有实验或经验证的最优样本量**。应与实验合作方确定 assay、细胞平台、重复单元和预实验功效。

**优先 pilot**：

- 选 8–12 条完整 αβ TCR（包括高密度功能阳性、阳性 singleton、近序列 discordant pair、功能阴性；尽量同 HLA 限制）。
- 每条 TCR 挑选 8–12 条 peptide：原始 cognate、少量不同位点的单氨基酸替换、presentation-matched decoy、确定的阴性对照；所有选择在看到结果前记录。
- 对每个受体–peptide 测 **至少 6 个预先规划的抗原剂量**（覆盖基线、过渡和饱和区；实际浓度由先导实验决定），至少 2–3 个独立生物学重复；不要用技术孔替代生物学重复。
- 用统一的 TCR-null 背景/受体表达条件，记录表面 TCR、HLA/pMHC 呈递、细胞活性与批次；必要时对关键 pair 独立测结合或 multimer、SPR/BLI（各自有适用限制）。
- 最好同时采集单细胞 **activated fraction** 与 **activated-cell conditional intensity**；NFAT/ERK（不同时间点）、CD69/细胞因子等分开解释。
- 对弱激动剂、部分激动剂、低呈递、表达差异给出独立控制；特别检查 anchor 突变是不是主要改变 HLA 稳定性。

**建议分析模型**（只有得到剂量曲线后）：

\[
P_{\mathrm{act}}(c)=b+\frac{E_{\max}-b}{1+(EC_{50}/c)^h}
\]

与简约单调模型、含随机效应的层级 logistic 比较；数据不支持则不强行自由估计 `h`。分别对 activated fraction 和 conditional intensity 建模，检验“增加刺激主要改变响应细胞比例”是否成立。

**判别框架**：

| 观测 | 支持的解释 | 不能推出 |
|---|---|---|
| 孤立阳性 TCR 与公共家族具有相近的激活曲线 | Sequence singleton 未必功能弱 | 它们在真实人群中没有相似受体 |
| 结合可检测但功能接近零 | 结合与信号之间存在差异 | 一定由特定 proofreading 步数造成 |
| 对 peptide 变体响应出现陡峭变化 | 存在识别/呈递/信号的非线性 | 单凭此证明严格的全或无细胞开关 |
| activated fraction 变而激活细胞强度近似不变 | 与数字式细胞响应相容 | 所有时间点和通路都完全数字化 |
| anchor 变体明显降低呈递而非 TCR 结合 | peptide presentation 是主要限制 | TCR 结合界面无需建模 |

### Gate 3

- **可以提出 binding→activation 假说**：同 pair 具备独立 binding + activation，并已控制表达与呈递。
- **只能提出 screen→flow 关联**：SEEKER 有 reporter 筛选/flow，但无独立结合数据。
- **只能讨论模型**：没有剂量与时间数据时，不能验证完整动力学校对或 Hill 参数。

---

## 6. Stage IV — 将两级模型联接，且允许输出“不确定/拒绝判断”

### 6.1 结构化输出（研究系统，不是临床结论）

```text
input:
  tcr: full paired alpha/beta + V/J + sequence QC
  target: peptide + HLA allele
  assay_context: optional dose, pMHC presentation, cell line, time
prediction:
  recognition_profile: ranked peptide affinities/enrichment proxies + coverage mask
  binding_readout_type: yeast_selection_proxy / measured_binding / unknown
  presentation_score: independent peptide-HLA control when valid
  activation_probability: only for validated assay context
  uncertainty: in_domain / extrapolation / abstain + reason
provenance:
  model_hash, training_split_hash, source_versions, calibration_domain
```

**禁止**给无 assay context 的跨 HLA 结果附上精确 activation probability；可以仅提供 recognition ranking 和模型域外警告。

### 6.2 三个独立评价层

1. **Recognition**：held-out peptide 与 held-out TCR-family 的读数/识别谱预测。
2. **Activation**：真实 flow 或独立功能 assay 下的校准与排序，按受体/实验单元 bootstrap。
3. **Transfer**：跨平台（yeast → cell）、跨 HLA、跨来源；**只有真实对应数据时才检验**。缺乏对应集不能把移植失败或成功编造为实验结论。

### 6.3 不确定性与 abstention

- 标记未知 HLA、远离训练序列家族、peptide scaffold OOD、输入 panel 非覆盖区域。
- 校准最好采用不与训练共享 TCR family / 原始研究的验证集。
- 报告 selective risk–coverage 曲线，而不是为了高 AUC 选择性删除难样本却不说明覆盖率。

---

## 7. 预注册的负对照和审计（必须先于模型结论）

| 负对照 | 排除的虚假机制 |
|---|---|
| 打乱 TCR–peptide 对应但保留各自边际组成 | 模型是否只靠 TCR/peptide 的典型性 |
| 同 HLA、同 anchor/呈递倾向的 matched decoy | HLA anchor/presentation 被误认成 TCR specificity |
| PRP naïve library 的 peptide count | 高 read count 只是合成/展示频率 |
| 仅 TCR、仅 peptide、仅 HLA-presentation 模型 | 是否真的学习了两者条件交互 |
| 重复数据、同一序列、同一亲本工程 TCR 强制同 split | exact/near duplicate 数据泄漏 |
| 同一 TCR 的跨轮次与跨实验扰动分析 | 时间/筛选 bottleneck 引发的虚假富集 |
| 真实未测量 = missing，非 0 | 截断 top-k 导致的假阴性和虚高 AUC |
| 同研究/同患者全体留出 | 来源、供者和实验批次的假泛化 |
| 固定模型，只做预先定义的评估 | 避免五个 Public 表位上反复选择配置 |

**统计单元**：PRP TCR（及亲本家族）、SEEKER 独立 TCR/原始实验批次、wetlab 独立生物学重复。数百万 peptide 行不是数百万个独立的 TCR 验证单元。

**报告**：每个估计给有效独立单元数、95% 置信区间或 bootstrap 区间（由合适聚类单元重抽样）、预先定义的主要终点与多重检验处理。若无独立单元支持 CI，则报告描述性结果，不伪装成精确推断。

---

## 8. 计算实施（本机为主，Nextflow DSL2；不使用 Snakemake）

### 8.1 推荐拓扑

```text
workflow/recognition_profile.nf
  FETCH_MANIFEST (手动确认许可及文件版本后)
  AUDIT_RAW_DATA
  NORMALIZE_PRP
  NORMALIZE_SEEKER
  AUDIT_PANEL_OVERLAP
  PROFILE_QC
  PROFILE_DISTANCE
  PREDEFINE_SPLITS
  FIT_PWM_BASELINE
  FIT_INTERACTION_BASELINE
  (OPTIONAL) FIT_PRP_PLM
  SCORE_HELDOUT
  AUDIT_FLOW_JOIN
  CALIBRATE_SCREEN_TO_FLOW
  (OPTIONAL) LINK_BINDING_TO_ACTIVATION
  SUMMARIZE_AND_FREEZE
```

**项目资源假设**：本机 WSL2 16 核、约 94 GB RAM，具有 RTX 5090；常规 count 处理与统计分析可优先用 CPU、Parquet/Polars；GPU 只在确认数据量和独立样本数足够时用于预训练模型推断或低秩微调。PRP 作者使用 ESM2-650M 的方案作为外部可复现基线，不默认首先重新全量训练它。

### 8.2 建议目录

```text
recognition_function/
  README.md
  nextflow.config
  main.nf
  conf/{local,gpu}.config
  modules/
  scripts/
    00_fetch_audit.py
    01_normalize_prp.py
    02_normalize_seeker.py
    03_panel_overlap.py
    04_profile_qc.py
    05_profile_distance.py
    06_make_splits.py
    07_fit_pwm.py
    08_fit_interactions.py
    09_eval_recognition.py
    10_flow_calibration.py
    11_activation_linkability.py
    12_eval_activation.py
  metadata/{sources.tsv,manifest.tsv,analysis_plan.md}
  tests/
  data/raw/                 # immutable, ignore in git
  data/derived/
  results/{figures,tables,models,reports}/
```

**这些是待实现的脚本名，而不是已经存在的可执行文件。** 每个 Nextflow process 都须写出输入文件、参数、版本和输出；整个结果目录由 pipeline manifest 与 Git commit 标记。

### 8.3 必须通过的自动化检查

```text
[ ] 原始文件 SHA256 固定；每个下载 URL 与许可证已记录
[ ] PRP anchors P2/P8 与 P2/P9 不被当作同一 peptide panel
[ ] 工程 TCR 和亲本绝不跨 train/valid/test
[ ] 所有跨 profile 距离仅用共同且可观察的 peptide
[ ] 不存在 selected-only SEEKER top-k 缺失值→负样本
[ ] 正负样本的检测流程与 assay_type 清楚区分
[ ] train-test exact TCR pair / peptide scaffold / parent-family leakage 0（按所定义 split）
[ ] 全部 flow 数据有可追踪原始行；无数据的 pair 不补值
[ ] 各主指标按正确独立单元重复评估
[ ] IMMREP26 冻结预测器与提交文件 SHA256 完全不变
```

---

## 9. 工作包、时间与停止门槛

以下是**相对工作阶段**，不是保证完成日期；重点是每阶段可独立交付、遇到不可识别即停止。

| 阶段 | 建议工作量 | 交付物 | 进入下一阶段所需的证据 |
|---|---|---|---|
| **P0 数据审计** | 1–2 个工作日 | 数据字典、重叠矩阵、来源/许可/负样本分类 | 数据确实可对齐；所有 missing/selected 状态明确 |
| **P1 描述性谱** | 2–4 天 | PRP 位点效应、SEEKER 筛选/flow 关联、共同 panel 矩阵 | 距离与 readout 在预处理变动下稳定 |
| **P2 功能距离** | 2–4 天 | sequence vs recognition distance 的反例与不确定性 | 有独立 TCR 及共同 peptide 覆盖，反例有实测支持 |
| **P3 预测基线** | 3–6 天 | M0/M1 严格 TCR-family/peptide-family 留出比较 | 明确优于基础模型且不靠泄漏 |
| **P4 功能解码** | 2–5 天 | flow 校准；binding→activation linkability 判定 | 有同 pair、不同读出、足够样本支持；否则停止映射训练 |
| **P5 wetlab 设计** | 依合作方条件 | 预注册 TCR/peptide 选择、剂量、重复与质控方案 | 实验平台、HLA/表达对照及伦理规范可行 |

**推荐首周只做 P0–P2**。若 P0 发现 PRP 与真实 activation 的配对覆盖不足，不要“想办法把两套数据拼成完整标签”；这正是该计划最重要的可证伪 Gate。

---

## 10. 预期科学产出与否证价值

### 成功路径

1. 数据驱动地复现 **sequence space ≠ recognition space** 的明确案例（两个方向均需要）；
2. 证明其超出 Pgen、HLA 呈递、实验读数偏差和近邻重复的解释；
3. 在独立的同 pair 功能实验上证明识别分数与激活既相关又不完全等价；
4. 在多剂量/单细胞实验中估计何时出现阈值或数字化，而不是根据二元标签假定全或无。

### 有价值的失败路径

- 无共同 peptide panel：研究 **within-TCR mutational landscapes**，不夸大跨受体距离。
- PRP 可预测 binding 但预测 activation 失败：说明 binding–function gap 真实且需要新的呈递/信号数据。
- 全部增益来自 HLA anchor：聚焦 pMHC presentation，而非宣称发现 TCR contact specificity。
- 只有少量患者来源或同亲本工程变体：只作机制个案，不宣称普适 αβ 泛化。
- 所有严格 OOD 验证失败：保留结果、明确技术与数据外推边界；不通过随机 split 营造成功。

---

## 11. 给代码 Agent 的执行指令（可直接复制）

> 目标：在不修改 IMMREP26 冻结提交的前提下，建立 `TCR sequence + peptide–HLA → recognition profile → functional activation` 的可复现实验框架。  
> 第一个任务仅实现 Gate 0、Stage I 描述性分析以及功能距离对照；**不得立即训练多任务大模型**。  
> 
> 1. 读取并严格对照 PRP Dryad、SEEKER Dryad 官方 README 与相应论文/代码仓库；下载与版本锁定需记录 SHA256、来源 URL、许可证和数据字典，不得猜测 CSV 列名。  
> 2. 构建保留 assay、HLA、peptide panel、轮次、重复、来源、`measured/not_observed/not_tested` 状态的统一长表；输出观察矩阵及可比较掩码。  
> 3. 分离 PRP yeast binding enrichment、SEEKER Jurkat-NFAT screen enrichment、flow activation；四参考受体 top-10,000 缺失条目绝不自动标负。  
> 4. 首先验证官方可复现的基础 QC、归一化与数值，再计算同面板 PRP recognition distances。将工程变体与亲本作为同一 TCR family。  
> 5. 做 sequence-distance 与 recognition-distance 的两种反例分析，使用共同观测 peptide、可信 interval、必要的背景对照。  
> 6. 审计 SEEKER flow 和 PRP 作者提供的 activation 数据：仅允许同 TCR、同 peptide、同 HLA 且条件清楚的 pair 进入 binding→activation 校准；匹配不足则标记 `NOT_IDENTIFIABLE`。  
> 7. 后续预测器从 PWM/线性模型开始；仅在严格 family/scaffold holdout 下产生增量时引入大型 pLM，不允许随机打散同一亲本或近邻变体。  
> 8. 交付 `analysis_plan.md`、`quality_report.md`、`overlap_matrix.csv`、`profile_distance.parquet`、`counterexamples.tsv`、`flow_linkability.md`、Nextflow DSL2 可复现流程、测试和哈希报告。  
> 9. 不更改 IMMREP26 模型、分数、提交文件和源码冻结清单，不得以 Public 标签选择本项目外推参数。  
> 10. 每个科学结论附原始证据、统计独立单元、控制项、替代解释与停止条件；**无法识别的参数必须明确标注，不能补充臆测数值。**

---

## 12. 参考文献与直接入口

1. Wang N. et al. (2026). *Deep peptide recognition profiling decodes TCR specificity and enables disease-associated antigen discovery*. **Nature Biotechnology**. <https://doi.org/10.1038/s41587-026-03128-x>. 数据 <https://doi.org/10.5061/dryad.h70rxwdzk>；代码 <https://github.com/akds/PRP-TCR-Specificity>.
2. Li S. et al. (2026). *SEEKER: A genome-scale library-on-library screening platform for deciphering T cell recognition of antigen*. **bioRxiv** 预印本. <https://doi.org/10.64898/2026.09.25.753432>. 数据 <https://doi.org/10.5061/dryad.79cnp5jct>；代码 <https://github.com/ldi2009/SEEKER-2026>.
3. Zhu X. et al. (2026). *Machine learning dataset and benchmark for accurate T cell receptor–pHLA binding prediction*. **Nature Computational Science**. <https://doi.org/10.1038/s43588-026-01073-y>. <https://github.com/Jiadong001/Hi-TpH>.
4. Morita S., Groves J. T. (2025). *Parallel reactions on a single T cell receptor offer a robust kinetic proofreading mechanism*. **PNAS**. <https://doi.org/10.1073/pnas.2514057122>.
5. *Antigen perception in T cells by long-term Erk and NFAT signaling dynamics* (2023/2024). <https://pubmed.ncbi.nlm.nih.gov/38113261/>. 有助于区分短时数字化和长期连续信号。

---

## 一句话决策

**第一阶段只回答：我们能否从真实数据中构建可靠、可比较的 TCR peptide recognition profiles，并证明它们提供了超出序列邻近性的可复现信息？**  
**第二阶段才回答：这种识别信息在什么呈递、剂量、表达与时间条件下转化为 T 细胞激活？**  
**如果缺少同 TCR/peptide/HLA 的双层实测数据，就停止两级联合拟合，而不是把筛选富集冒充结合亲和力。**
