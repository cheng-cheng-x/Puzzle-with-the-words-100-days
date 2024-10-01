import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.manifold import TSNE

# 读取和处理 fpkm 数据
fpkm = pd.read_csv("/home/gongfengcz/RNAseq/data/脓毒症/fpkm_LPS.txt", sep="\t", index_col=0)
fpkm_t = fpkm.transpose()
organs = [name.split("_")[0] for name in fpkm_t.index]
fpkm_t = pd.DataFrame(fpkm_t)
fpkm_t['Organ'] = organs
fpkm_data = fpkm_t.drop(columns=['Organ'])

# 运行 t-SNE 降维
tsne = TSNE(n_components=2, perplexity=30, n_iter=1000, random_state=913223)
tsne_result = tsne.fit_transform(fpkm_data)
tsne_data = pd.DataFrame(tsne_result, columns=["Dim1", "Dim2"])
tsne_data['Sample'] = fpkm_t.index
tsne_data['Organ'] = fpkm_t['Organ']

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import chi2
from matplotlib.patches import Ellipse  # 导入 Ellipse

# 假设 Sample 列中的值是器官名 
tsne_data['Organ'] = tsne_data['Sample'].str.split('_').str[0]  
print(tsne_data.head()) 
# 可视化 t-SNE 结果
plt.figure(figsize=(20, 20))
colors = plt.cm.tab20
x_min, x_max, y_min, y_max = np.inf, -np.inf, np.inf, -np.inf

for i, organ in enumerate(tsne_data['Organ'].unique()):
    subset = tsne_data[tsne_data['Organ'] == organ]
    print(f"Organ: {organ}, Number of samples: {len(subset)}")
    plt.scatter(subset['Dim1'], subset['Dim2'], s=50, color=colors(i), alpha=0.8)
    plt.text(subset['Dim1'].mean(), subset['Dim2'].mean(), organ, fontsize=14, fontweight='bold', ha='center', va='center', color='black', family='serif')

    # 计算和绘制置信椭圆
    mean = subset[['Dim1', 'Dim2']].mean().values
    cov = np.cov(subset[['Dim1', 'Dim2']], rowvar=False)
    chi2_val = chi2.ppf(0.99, df=2)
    eigenvalues, eigenvectors = np.linalg.eig(cov)

    angle = np.arctan2(eigenvectors[1, 0], eigenvectors[0, 0])
    ellipse_width = np.sqrt(eigenvalues[0] * chi2_val)
    ellipse_height = np.sqrt(eigenvalues[1] * chi2_val)
    
    ellipse = Ellipse(mean, width=ellipse_width * 2, height=ellipse_height * 2, angle=np.degrees(angle),
                      color=colors(i), alpha=0.2, linewidth=2, fill=True)
    plt.gca().add_artist(ellipse)

    x_min = min(x_min, mean[0] - ellipse_width)
    x_max = max(x_max, mean[0] + ellipse_width)
    y_min = min(y_min, mean[1] - ellipse_height)
    y_max = max(y_max, mean[1] + ellipse_height)

# 设置坐标轴范围，添加标题和图例
plt.xlim(x_min - 0.5, x_max + 0.5)
plt.ylim(y_min - 0.5, y_max + 0.5)
plt.title("t-SNE of RNA-seq Data by Organ", fontsize=20)
plt.xlabel("t-SNE 1", fontsize=16)
plt.ylabel("t-SNE 2", fontsize=16)

legend_labels = {organ: colors(i) for i, organ in enumerate(tsne_data['Organ'].unique())}
handles = [plt.Line2D([0], [0], marker='o', color='w', label=organ, 
                       markerfacecolor=color, markersize=10) for organ, color in legend_labels.items()]
plt.legend(handles=handles, title="Organ", bbox_to_anchor=(1.05, 1), loc='upper left')

plt.grid(False)
plt.show()

plt.savefig("/home/mahc/data/nongdu/剔除样本/RNAseq/figure/tSNE/tSNE_30_plot.png", dpi=300, bbox_inches='tight')

#交集基因
import pandas as pd

#线粒体基因集
mito_genes = pd.read_excel("/home/gongfengcz/scRNA-heart-mitochodria/data/Mouse.MitoCarta3.0.xls")
mito_genes['Genes'] = mito_genes['Genes'].str.split(', ')  # 分割基因列表
mito_genes_list = mito_genes.explode('Genes')['Genes'].tolist()  # 展开并转换为列表

# 免疫基因集
immune_genes = pd.read_csv("/home/hutb/RNAseq/data/inflammatory.csv")
immune_genes_list = immune_genes['gene_symbol'].tolist()  # 假设免疫基因集列名为 'gene_symbol'


fpkm = pd.read_csv("/home/gongfengcz/RNAseq/data/脓毒症/fpkm_LPS.txt", sep="\t", index_col=0)
mito_fpkm_intersection = fpkm.index.intersection(mito_genes_list)
immune_fpkm_intersection = fpkm.index.intersection(immune_genes_list)

mito_fpkm_data = fpkm.loc[mito_fpkm_intersection]
immune_fpkm_data = fpkm.loc[immune_fpkm_intersection]

print("线粒体基因与FPKM的交集基因数量：", len(mito_fpkm_intersection))
print("免疫基因与FPKM的交集基因数量：", len(immune_fpkm_intersection))

print("线粒体基因与FPKM的交集基因表达数据：\n", mito_fpkm_data)
print("免疫基因与FPKM的交集基因表达数据：\n", immune_fpkm_data)

#线粒体基因交集tSNE图
import pandas as pd
from sklearn.manifold import TSNE
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import chi2
from matplotlib.patches import Ellipse  # 导入 Ellipse

# 读取和处理线粒体基因交集的 fpkm 数据
mito_fpkm_data_t = mito_fpkm_data.transpose()
organs = [name.split("_")[0] for name in mito_fpkm_data_t.index]
mito_fpkm_data_t = pd.DataFrame(mito_fpkm_data_t)
mito_fpkm_data_t['Organ'] = organs
mito_fpkm_data_clean = mito_fpkm_data_t.drop(columns=['Organ'])

# 运行 t-SNE 降维
tsne = TSNE(n_components=2, perplexity=15, n_iter=500, random_state=913223)
tsne_result = tsne.fit_transform(mito_fpkm_data_clean)
tsne_data = pd.DataFrame(tsne_result, columns=["Dim1", "Dim2"])
tsne_data['Sample'] = mito_fpkm_data_t.index
tsne_data['Organ'] = mito_fpkm_data_t['Organ']

# 假设 Sample 列中的值是器官名
tsne_data['Organ'] = tsne_data['Sample'].str.split('_').str[0]
print(tsne_data.head()) 

# 可视化 t-SNE 结果
plt.figure(figsize=(20, 20))
colors = plt.cm.tab20
x_min, x_max, y_min, y_max = np.inf, -np.inf, np.inf, -np.inf

for i, organ in enumerate(tsne_data['Organ'].unique()):
    subset = tsne_data[tsne_data['Organ'] == organ]
    print(f"Organ: {organ}, Number of samples: {len(subset)}")
    plt.scatter(subset['Dim1'], subset['Dim2'], s=50, color=colors(i), alpha=0.8)
    plt.text(subset['Dim1'].mean(), subset['Dim2'].mean(), organ, fontsize=14, fontweight='bold', ha='center', va='center', color='black', family='serif')

    # 计算和绘制置信椭圆
    mean = subset[['Dim1', 'Dim2']].mean().values
    cov = np.cov(subset[['Dim1', 'Dim2']], rowvar=False)
    chi2_val = chi2.ppf(0.99, df=2)
    eigenvalues, eigenvectors = np.linalg.eig(cov)

    angle = np.arctan2(eigenvectors[1, 0], eigenvectors[0, 0])
    ellipse_width = np.sqrt(eigenvalues[0] * chi2_val)
    ellipse_height = np.sqrt(eigenvalues[1] * chi2_val)
    
    ellipse = Ellipse(mean, width=ellipse_width * 2, height=ellipse_height * 2, angle=np.degrees(angle),
                      color=colors(i), alpha=0.2, linewidth=2, fill=True)
    plt.gca().add_artist(ellipse)

    x_min = min(x_min, mean[0] - ellipse_width)
    x_max = max(x_max, mean[0] + ellipse_width)
    y_min = min(y_min, mean[1] - ellipse_height)
    y_max = max(y_max, mean[1] + ellipse_height)

# 设置坐标轴范围，添加标题和图例
plt.xlim(x_min - 0.5, x_max + 1)
plt.ylim(y_min - 0.5, y_max + 0.5)
plt.title("t-SNE of Mitochondrial Gene Expression by Organ", fontsize=20)
plt.xlabel("t-SNE 1", fontsize=16)
plt.ylabel("t-SNE 2", fontsize=16)

legend_labels = {organ: colors(i) for i, organ in enumerate(tsne_data['Organ'].unique())}
handles = [plt.Line2D([0], [0], marker='o', color='w', label=organ, 
                       markerfacecolor=color, markersize=10) for organ, color in legend_labels.items()]
plt.legend(handles=handles, title="Organ", bbox_to_anchor=(1.05, 1), loc='upper left')

plt.grid(False)
plt.show()

# 保存图片
plt.savefig("/home/mahc/data/nongdu/剔除样本/RNAseq/figure/tSNE/tSNE_15_plot_mitochondrial.png", dpi=300, bbox_inches='tight')


#炎症基因交集tSNE图
import pandas as pd
from sklearn.manifold import TSNE
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import chi2
from matplotlib.patches import Ellipse  # 导入 Ellipse

# 读取和处理炎症基因交集的 fpkm 数据
immune_fpkm_data_t = immune_fpkm_data.transpose()
organs = [name.split("_")[0] for name in immune_fpkm_data_t.index]
immune_fpkm_data_t = pd.DataFrame(immune_fpkm_data_t)
immune_fpkm_data_t['Organ'] = organs
immune_fpkm_data_clean = immune_fpkm_data_t.drop(columns=['Organ'])

# 运行 t-SNE 降维
tsne = TSNE(n_components=2, perplexity=22, n_iter=1000, random_state=913223)
tsne_result = tsne.fit_transform(immune_fpkm_data_clean)
tsne_data = pd.DataFrame(tsne_result, columns=["Dim1", "Dim2"])
tsne_data['Sample'] = immune_fpkm_data_t.index
tsne_data['Organ'] = immune_fpkm_data_t['Organ']

# 假设 Sample 列中的值是器官名
tsne_data['Organ'] = tsne_data['Sample'].str.split('_').str[0]
print(tsne_data.head()) 

# 可视化 t-SNE 结果
plt.figure(figsize=(20, 20))
colors = plt.cm.tab20
x_min, x_max, y_min, y_max = np.inf, -np.inf, np.inf, -np.inf

for i, organ in enumerate(tsne_data['Organ'].unique()):
    subset = tsne_data[tsne_data['Organ'] == organ]
    print(f"Organ: {organ}, Number of samples: {len(subset)}")
    plt.scatter(subset['Dim1'], subset['Dim2'], s=50, color=colors(i), alpha=0.8)
    plt.text(subset['Dim1'].mean(), subset['Dim2'].mean(), organ, fontsize=14, fontweight='bold', ha='center', va='center', color='black', family='serif')

    # 计算和绘制置信椭圆
    mean = subset[['Dim1', 'Dim2']].mean().values
    cov = np.cov(subset[['Dim1', 'Dim2']], rowvar=False)
    chi2_val = chi2.ppf(0.99, df=2)
    eigenvalues, eigenvectors = np.linalg.eig(cov)

    angle = np.arctan2(eigenvectors[1, 0], eigenvectors[0, 0])
    ellipse_width = np.sqrt(eigenvalues[0] * chi2_val)
    ellipse_height = np.sqrt(eigenvalues[1] * chi2_val)
    
    ellipse = Ellipse(mean, width=ellipse_width * 2, height=ellipse_height * 2, angle=np.degrees(angle),
                      color=colors(i), alpha=0.2, linewidth=2, fill=True)
    plt.gca().add_artist(ellipse)

    x_min = min(x_min, mean[0] - ellipse_width)
    x_max = max(x_max, mean[0] + ellipse_width)
    y_min = min(y_min, mean[1] - ellipse_height)
    y_max = max(y_max, mean[1] + ellipse_height)

# 设置坐标轴范围，添加标题和图例
plt.xlim(x_min - 1, x_max + 1)
plt.ylim(y_min - 1, y_max + 1)
plt.title("t-SNE of Inflammatory Gene Expression by Organ", fontsize=20)
plt.xlabel("t-SNE 1", fontsize=16)
plt.ylabel("t-SNE 2", fontsize=16)

legend_labels = {organ: colors(i) for i, organ in enumerate(tsne_data['Organ'].unique())}
handles = [plt.Line2D([0], [0], marker='o', color='w', label=organ, 
                       markerfacecolor=color, markersize=10) for organ, color in legend_labels.items()]
plt.legend(handles=handles, title="Organ", bbox_to_anchor=(1.05, 1), loc='upper left')

plt.grid(False)
plt.show()

# 保存图片
plt.savefig("/home/mahc/data/nongdu/剔除样本/RNAseq/figure/tSNE/tSNE_22_plot_inflammatory.png", dpi=300, bbox_inches='tight')
