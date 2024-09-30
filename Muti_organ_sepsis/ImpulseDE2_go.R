library(clusterProfiler)
library(org.Mm.eg.db)
library(reshape2)
library(ggplot2)
library(cowplot)
library(readxl)

# 创建一个扩展的函数，查询对应器官的模式、基因数量和基因名，并赋值给指定变量
query_gene_patterns <- function(organ, pattern = NULL, assign_var = NULL) {
  
  if (!organ %in% names(all_gene_patterns)) {
    stop(paste("Organ", organ, "not found in the data."))
  }
  
  
  organ_patterns <- all_gene_patterns[[organ]]
  
  
  if (!is.null(pattern)) {
    if (pattern %in% names(organ_patterns)) {
      pattern_data <- organ_patterns[[pattern]]
      gene_names <- pattern_data$genes
      
      if (!is.null(assign_var)) {
        assign(assign_var, gene_names, envir = .GlobalEnv)
      }
      
      return(list(Pattern = pattern, Count = pattern_data$count, Genes = gene_names))
    } else {
      return(paste("Pattern", pattern, "not found for organ", organ))
    }
  }
  
  pattern_summary <- sapply(organ_patterns, function(x) x$count)
  return(data.frame(Pattern = names(organ_patterns), Count = pattern_summary))
}

load("~/data/nongdu/剔除样本/RNAseq/figure/差异基因IDE2/filtered_all_gene_patterns.RData")
# 定义一个函数，对每个器官执行GO富集分析并保存图像
perform_go_analysis <- function(organ, pattern, output_dir) {
  # 查询基因模式
  query_gene_patterns(organ, pattern = pattern, assign_var = "genes_negative_pattern")
  
  # 执行GO富集分析
  ego <- enrichGO(
    gene = genes_negative_pattern,
    OrgDb = org.Mm.eg.db,
    ont = "BP",  # 生物过程 (Biological Process)
    pAdjustMethod = 'BH',
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    keyType = 'SYMBOL'
  )
  
  # 保存barplot
  barplot(ego, showCategory = 20)
  barplot_filename <- paste0(output_dir, "/barplot_", organ, ".png")
  ggsave(
    filename = barplot_filename,
    width = 15, height = 15,
    dpi = 500
  )
  
  # 生成dotplot
  dotplot_obj <- dotplot(ego, showCategory = 20)
  
  # 保存dotplot
  dotplot_filename <- paste0(output_dir, "/dotplot_", organ, ".png")
  ggsave(
    filename = dotplot_filename,
    plot = dotplot_obj,
    width = 15, height = 15,
    dpi = 500
  )
}
















# 读取Mouse.MitoCarta3.0的基因列表
dbs <- read_excel("/home/gongfengcz/scRNA-heart-mitochodria/data/Mouse.MitoCarta3.0.xls")
# 根据需要，可能需要对dbs数据框进行适当的调整
gene_df <- dbs %>%
  mutate(gene = strsplit(Genes, ", ")) %>%
  unnest(cols = gene) # 确保这个步骤按照数据结构来调整



# 筛选差异基因
counts <- read.table(
  file = "/home/gongfengcz/RNAseq/data/脓毒症/GSE224127_1_counts_LPS_timecourse.txt",
  row.names = 1, header = TRUE)
filtered_counts <- counts[rownames(counts) %in% gene_df$gene, ]
counts <-filtered_counts 
sample_info <- colnames(counts)
sample_time <- sapply(strsplit(sample_info, "_"), function(x) paste(x[1], x[2], sep = "_"))

aggregated_counts <- list()


for (time in unique(sample_time)) {
  matching_samples <- grep(paste0("^", time), sample_info, value = TRUE)
  aggregated_counts[[time]] <- rowSums(counts[, matching_samples])
}

aggregated_counts_matrix <- do.call(cbind, aggregated_counts)
colnames(aggregated_counts_matrix) <- unique(sample_time)


pbmc_samples <- grep("^PBMC", colnames(aggregated_counts_matrix), value = TRUE)
pbmc_counts_matrix <- aggregated_counts_matrix[, pbmc_samples]


pbmc_genes <- read.csv("~/data/nongdu/剔除样本/RNAseq/figure/差异基因IDE2/ImpulseDE2_PBMC_results.csv")

# 设置筛选阈值
p_value_threshold <- 0.05
padj_threshold <- 0.05


significant_genes <- pbmc_genes[pbmc_genes$p < p_value_threshold & pbmc_genes$padj < padj_threshold, ]


diff_gene_names <- significant_genes$Gene  



ego <- enrichGO(
  gene = diff_gene_names,
  OrgDb = org.Mm.eg.db,
  ont = "BP",  # 生物过程 (Biological Process)
  pAdjustMethod = 'BH',
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  keyType = 'SYMBOL'  )



query_gene_patterns("HE", pattern = "+--+-+",assign_var = "genes_negative_pattern")

ego <- enrichGO(
  gene = genes_negative_pattern,
  OrgDb = org.Mm.eg.db,
  ont = "BP",  # 生物过程 (Biological Process)
  pAdjustMethod = 'BH',
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  keyType = 'SYMBOL' 
)




# 定义所有器官的列表和模式
organs <- c("BM", "TH", "SP", "iLN", "PBMC", "BR", "HE", "LU", "LI", "KI", "SI", "CO", "SK")
pattern <- "+--+-+"
output_dir <- "~/data/nongdu/剔除样本/RNAseq/figure/差异基因IDE2/GO/Mito/+--+-+"

# 为每个器官执行分析
for (organ in organs) {
  perform_go_analysis(organ, pattern, output_dir)
}
















# 显示前20个富集的GO分类
barplot(ego, showCategory = 20)

dotplot(ego, showCategory = 20)

# 显示富集结果
query_gene_patterns("HE", pattern = "+--+-+",assign_var = "genes_negative_pattern")

ego <- enrichGO(
  gene = genes_negative_pattern,
  OrgDb = org.Mm.eg.db,
  ont = "BP",  # 生物过程 (Biological Process)
  pAdjustMethod = 'BH',
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  keyType = 'SYMBOL' 
)
barplot(ego, showCategory = 20)
dotplot_obj <- dotplot(ego, showCategory = 20)
ggsave(
  filename = "~/data/nongdu/剔除样本/RNAseq/figure/差异基因IDE2/GO/Mito/+--+-+/dotplot_HE.png", 
  plot = dotplot_obj,                      
  width = 15, height = 15,                   
  dpi = 500                                
)

cnetplot(ego,circular=FALSE,colorEdge = TRUE)
dotplot_obj <- cnetplot(ego, circular = FALSE, color.params = list(edge = TRUE))

ggsave(
  filename = "~/data/nongdu/剔除样本/RNAseq/figure/差异基因IDE2/GO/Mito/PBMC_cnetplot.png", 
  plot = dotplot_obj,                      
  width = 30, height = 30,                   
  dpi = 300                                
)

plotGOgraph_obj<-plotGOgraph(ego)
# 打开图像设备
png("~/data/nongdu/剔除样本/RNAseq/figure/差异基因IDE2/GO/Mito/plotGOgraph.png",width = 2000, height = 1600, res = 300)

# 生成图形
plotGOgraph(ego)

# 关闭图像设备，保存文件
dev.off()


plotGOgraph(ego)

library(enrichplot)
emapplot(ego)

heatplot(ego)

gseaplot(ego)

ego_sim <- pairwise_termsim(ego)
emapplot(ego_sim, showCategory = 50, color = "p.adjust", layout = "kk")

# 将 emapplot 的结果保存为对象
emap_obj <- emapplot(ego_sim, showCategory = 50, color = "p.adjust", layout = "kk")

# 使用 ggsave 保存图像
ggsave(
  filename = "~/data/nongdu/剔除样本/RNAseq/figure/差异基因IDE2/GO/Mito/emapplot_GO_network.png",  # 替换为你的实际保存路径
  plot = emap_obj,  # 使用 emapplot 的对象
  width = 20, height = 20,  # 设置图片尺寸
  dpi = 300  # 设置图片分辨率
)



cnet_obj <- cnetplot(ego, showCategory = 10, circular = TRUE, colorEdge = TRUE, 
                     label_category = TRUE, label_gene = TRUE)


# 使用 ggsave 保存图像
ggsave(
  filename = "~/data/nongdu/剔除样本/RNAseq/figure/差异基因IDE2/GO/Mito/cnetplot_GO_chord.png",  # 替换为你的实际保存路径
  plot = cnet_obj,  # 使用 cnetplot 的对象
  width = 15, height = 15,  # 设置图片尺寸
  dpi = 300  # 设置图片分辨率
)

library(GOplot)
# 提取富集结果中的 GO 数据
go_data <- ego@result
go_data <- go_data[, c("ID", "Description", "p.adjust")]

# 过滤显著富集的 GO 条目（可根据需要设置阈值）
go_data <- go_data[go_data$p.adjust < 0.05, ]