# 文件命名规范化更新

## 文件重命名对照表

### 输入数据文件 (Input Files)
| 原文件名 | 新文件名 | 说明 |
|---------|---------|------|
| `all_pathway_genes.csv` | `input_gene_pathway_data.csv` | 基础基因路径数据 |
| `pathway_gene_summary.csv` | `input_cd3_pathway_genes.csv` | CD3路径基因汇总 |
| `bcr_pathway_summary.csv` | `input_bcell_pathway_genes.csv` | B细胞受体路径汇总 |
| `macrophage_pathway_gene_summary.csv` | `input_macrophage_pathway_genes.csv` | 巨噬细胞路径基因汇总 |
| `gene_protein_attention_scores.csv` | `input_attention_scores.csv` | 基因-蛋白质注意力得分数据 |
| `merged_edges.csv` | `input_gene_peak_network.csv` | 基因-Peak网络连接数据 |

### 输出数据文件 (Output Files)
| 原文件名 | 新文件名 | 说明 |
|---------|---------|------|
| `genes_with_pathway_info.csv` | `output_genes_with_pathways.csv` | 包含路径信息的完整基因数据 |
| `analysis_results.csv` | `output_analysis_results.csv` | 分析结果数据 |
| `pathway_comparison.png` | `output_pathway_comparison.png` | 路径比较可视化图像 |
| `pathway_comparison.pdf` | `output_pathway_comparison.pdf` | 路径比较可视化PDF |

### 保持不变的文件
- `rna_data0320.h5ad` - RNA单细胞数据
- `protein_data0320.h5ad` - 蛋白质单细胞数据
- `gene_peak_attention_matrix.npz` - 注意力矩阵文件

## 代码修改说明

### 0412Visualize.py 修改内容
1. **增强文件查找功能**：
   - 修改 `load_data_files()` 方法，使用 `_find_file()` 智能查找文件
   - 支持新旧文件名的自动匹配，确保向后兼容

2. **文件路径更新**：
   - 主基因数据：优先使用 `input_gene_pathway_data.csv`
   - CD3路径：优先使用 `input_cd3_pathway_genes.csv`
   - B细胞路径：优先使用 `input_bcell_pathway_genes.csv`
   - 巨噬细胞路径：优先使用 `input_macrophage_pathway_genes.csv`
   - 注意力得分：优先使用 `input_attention_scores.csv`

3. **输出文件名规范化**：
   - 输出文件改为 `output_genes_with_pathways.csv`

### analysis_vis.py 修改内容
1. **输入文件路径更新**：
   - 基因边缘数据：优先查找 `output_genes_with_pathways.csv`
   - 网络数据：优先使用 `input_gene_peak_network.csv`

2. **输出文件名规范化**：
   - 分析结果：改为 `output_analysis_results.csv`
   - 可视化图像：改为 `output_pathway_comparison.png/pdf`

## 命名规范原则

### 1. 前缀规范
- `input_`: 输入数据文件
- `output_`: 程序生成的输出文件

### 2. 名称结构
- 使用下划线分隔单词
- 采用描述性命名，明确表达文件内容
- 避免在文件名中包含参数信息

### 3. 向后兼容性
- 所有修改都包含fallback机制
- 优先使用新文件名，如果不存在则自动使用旧文件名
- 确保代码在新旧文件环境下都能正常运行

## 使用说明

### 运行测试
```bash
# 测试基因路径处理
python 0412Visualize.py
# 预期输出: Processing complete. Updated file saved as './output_genes_with_pathways.csv'

# 测试注意力分析和可视化
python analysis_vis.py
# 预期输出: Analysis completed successfully!
#          Visualization generated: output_pathway_comparison.png
```

### 文件依赖关系
1. `0412Visualize.py` 处理原始数据，生成 `output_genes_with_pathways.csv`
2. `analysis_vis.py` 使用该文件进行注意力分析和可视化

## 优势特性

1. **自动适配**：代码会自动查找可用的输入文件
2. **智能选择**：优先使用规范化的新文件名
3. **向后兼容**：支持旧文件名，确保代码健壮性
4. **清晰命名**：文件名直观表达内容和用途
5. **统一规范**：采用一致的命名约定 