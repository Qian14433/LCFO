# LCFO: Life-Cycle Caterpillar Fungus Optimizer

**生命周期冬虫夏草优化算法** — 对原始 CFO（Yang et al., 2026）的诊断驱动生命周期重构版本。

*A diagnosis-driven life-cycle re-formulation of the Caterpillar Fungus Optimizer (CFO).
A mechanism-level diagnosis shows that CFO's high-dimensional degradation stems from four
structural defects, including the generation-wise repeated parasitism phase that degenerates
into replicating the current best solution. LCFO abolishes the repeated parasitism behavior,
redefines infection as a one-shot life-cycle transition, and organizes the post-infection
search around five components: multi-host differential propagation, partial colonization,
successful-colonization memory, a fixed-capacity dormancy archive, and seasonal
active-population regulation.*

---

## 文件清单 / Files

| 路径 | 说明 |
|---|---|
| `code/LCFO.m` | 本文提出的算法 / The proposed algorithm |
| `code/CFO.m` | 原始 CFO 基线 / Original CFO baseline |
| `code/initialization.m` | 种群初始化（共用）/ Shared initializer |
| `code/main_demo.m` | 演示脚本 / Demo script |
| `data/CEC2017_D{30,50,100}_results.mat` | 36 种算法 × 29 函数的完整运行结果（runs 29×30、curves 29×30×100、elapsed） |
| `data/CEC2017_D{30,50,100}_final_errors.csv` | 逐次运行最终误差（长表格式） |
| `results/表S1–S3_*.csv` | LCFO 与 CEC 冠军算法的逐函数 Mean/Std/Wilcoxon 统计（论文补充材料表S1–S3） |

说明：函数编号 F1, F3–F30（F2 按 CEC2017 官方惯例剔除）。

## 使用方法 / Usage

```matlab
fobj = @(x) sum(x.^2);          % 目标函数（最小化）
dim = 30; lb = -100; ub = 100;

% LCFO：FE 预算驱动，初始种群固定 100
[Best_score, Best_pos, curve] = LCFO(100, 30000, lb, ub, dim, fobj);

% 原始 CFO：迭代驱动
[Best_score, Best_pos, curve] = CFO(30, 499, lb, ub, dim, fobj);
```

## 主要结果 / Benchmark highlights

CEC2017（29 函数，MaxFEs = 30000，30 次独立运行）：

- 对 CEC 竞赛获胜算法（JADE、SHADE、L-SHADE、jSO）的 Friedman 平均秩：
  **1.24（D=30）/ 1.17（D=50）/ 1.34（D=100）**，均列第一；
- 在 13 种经典算法 + 18 种新近算法 + 4 种 CEC 获胜算法共 36 种算法的大规模比较中，
  三维平均秩 **1.44**，排名第一；
- CEC2020-RW 24 个真实世界约束问题：可行性优先规则下平均秩 **2.29**，优于全部对比算法。

## 环境 / Requirements

- MATLAB R2020a+（开发环境 R2025a），无工具箱依赖
- CEC2017 函数请从竞赛官网获取（`cec17_func` MEX + `input_data`，版权原因未收录）

## 引用 / Citation

```bibtex
% 论文发表后更新 / To be updated upon publication
@article{LCFO2026,
  title  = {From repeated parasitism to life-cycle search: a mechanism-guided
            caterpillar fungus optimizer for high-dimensional optimization},
  author = {Anonymous},
  year   = {2026},
  note   = {Manuscript under review}
}
```

## 参考文献

- Yang et al. A novel bio-inspired caterpillar fungus (*Ophiocordyceps sinensis*) optimizer
  for SOFC parameter identification via GRNN. *Renewable Energy* 256 (2026) 123995.
- Tanabe & Fukunaga. Success-history based parameter adaptation for Differential Evolution
  (SHADE). IEEE CEC 2013: 71-78.
