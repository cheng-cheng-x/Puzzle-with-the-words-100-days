
if (!requireNamespace("ImpulseDE2", quietly = TRUE)) {
  BiocManager::install("ImpulseDE2")
}
install.packages("ImpulseDE2")
library(ImpulseDE2)

# 读取数据
counts <- read.table(
  file = "/home/gongfengcz/RNAseq/data/脓毒症/GSE224127_1_counts_LPS_timecourse.txt",
  row.names = 1,
  header = TRUE  
)

# 将表达量小于 10 的基因设置为 0
#counts[counts < 10] <- 0
####样本信息-----
metadata <- t(counts)
metadata <- as.data.frame(metadata)
metadata$tissue <- sapply(strsplit(rownames(metadata), split = "_"), "[", 1)
metadata$time <- sapply(strsplit(rownames(metadata), split = "_"), "[", 2)
metadata$tissue<-as.factor(metadata$tissue)
metadata$sample<-rownames(metadata)
metadata <- as.data.frame(metadata)
metadata$tissue <- as.factor(metadata$tissue)
metadata$time <- as.factor(metadata$time)

metadata$group <- paste0(metadata$tissue,"_",metadata$time)
metadata$sample <- paste0(rownames(metadata),".gencode.vM19")
rownames(metadata) <- metadata$sample


organs <- c("BM", "TH", "SP", "iLN", "PBMC", "BR", "HE", "LU", "LI", "KI", "SI", "CO", "SK")


for (organ in organs) {
  
  
  organ_columns <- grep(paste0("^", organ), colnames(counts), value = TRUE)
  counts_organ <- counts[, organ_columns]
  

  Sample <- colnames(counts_organ)
  

  Time <- gsub(paste0(organ, "_(d[0-9.]+).*"), "\\1", Sample)
  
  Condition <- rep("case", length(Sample)) 
  Batch <- rep("B_NULL", length(Sample))
  
  dfAnnotation <- data.frame(Sample, Condition, Time, Batch, row.names = Sample)
  
  counts_organ_matrix <- as.matrix(counts_organ)

  dfAnnotation$Time <- as.numeric(gsub("d", "", dfAnnotation$Time))
  
  # 运行 ImpulseDE2 分析
  objectImpulseDE2 <- runImpulseDE2(
    matCountData    = counts_organ_matrix,
    dfAnnotation    = dfAnnotation,
    boolCaseCtrl    = FALSE,
    vecConfounders  = NULL,
    scaNProc        = 55
  )
  

  resultTable <- objectImpulseDE2$dfImpulseDE2Results
  

  cleaned_results <- resultTable[!is.na(resultTable$p), ]
  significant_genes <- cleaned_results[cleaned_results$p < 0.05, ]
  
  # 保存结果
  output_path <- paste0("~/data/nongdu/剔除样本/RNAseq/figure/差异基因IDE2/ImpulseDE2_", organ, "_results.csv")
  write.csv(resultTable, output_path, row.names = FALSE)
  

  cat("Finished processing", organ, "\n")
}



# 筛选 PBMC 的列
pbmc_columns <- grep("^PBMC", colnames(counts), value = TRUE)
counts_pbmc <- counts[, pbmc_columns]


Sample <- colnames(counts_pbmc)

Time <- gsub("PBMC_(d[0-9.]+).*", "\\1", Sample)

Condition <- rep("case", length(Sample))  
Batch <- rep("B_NULL", length(Sample))

dfAnnotation <- data.frame(Sample, Condition, Time, Batch, row.names = Sample)

# 转换为矩阵
counts_pbmc_matrix <- as.matrix(counts_pbmc)



dfAnnotation$Time <- as.numeric(gsub("d", "", dfAnnotation$Time))



objectImpulseDE2 <- runImpulseDE2(
  matCountData    = counts_pbmc_matrix ,
  dfAnnotation    = dfAnnotation,
  boolCaseCtrl    = FALSE,
  vecConfounders  = NULL,
  scaNProc        = 55
)

# 查看结果
resultTable <- objectImpulseDE2$dfImpulseDE2Results
head(resultTable)
# 筛选显著差异表达的基因
cleaned_results <- objectImpulseDE2$dfImpulseDE2Results[!is.na(objectImpulseDE2$dfImpulseDE2Results$p), ]
significant_genes <- cleaned_results[cleaned_results$p < 0.05, ]


write.csv(resultTable, "~/data/nongdu/剔除样本/RNAseq/figure/差异基因IDE2/ImpulseDE2_LI_results.csv", row.names = FALSE)








# 筛选差异基因
counts <- read.table(
  file = "/home/gongfengcz/RNAseq/data/脓毒症/GSE224127_1_counts_LPS_timecourse.txt",
  row.names = 1, header = TRUE)

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


pbmc_diff_counts <- pbmc_counts_matrix[rownames(pbmc_counts_matrix) %in% diff_gene_names, ]


time_points <- gsub("PBMC_", "", pbmc_samples)
sorted_time_points <- time_points[order(as.numeric(sub("d", "", time_points)))]


pbmc_counts_matrix <- pbmc_counts_matrix[, paste0("PBMC_", sorted_time_points)]


gene_patterns <- list()


for (i in 1:nrow(pbmc_diff_counts)) {
  gene <- pbmc_diff_counts[i, ]
  pattern <- ""
  
  for (j in 2:length(sorted_time_points)) {
    if (gene[j] > gene[j-1]) {
      pattern <- paste0(pattern, "+")
    } else if (gene[j] < gene[j-1]) {
      pattern <- paste0(pattern, "-")
    } else {
      pattern <- paste0(pattern, "0")
    }
  }
  

  if (!pattern %in% names(gene_patterns)) {
    gene_patterns[[pattern]] <- list(genes = c(), count = 0)
  }
  gene_patterns[[pattern]]$genes <- c(gene_patterns[[pattern]]$genes, rownames(pbmc_diff_counts)[i])
  gene_patterns[[pattern]]$count <- gene_patterns[[pattern]]$count + 1
}


heatmap_matrix <- do.call(rbind, lapply(names(gene_patterns), function(pattern) {
  sapply(strsplit(pattern, "")[[1]], function(char) {
    if (char == "+") {
      return(1)  # 红色
    } else if (char == "-") {
      return(-1)  # 蓝色
    } else {
      return(0)  # 灰色
    }
  })
}))


time_points <- c("d0.25", "d0.5", "d1", "d2", "d3", "d5")
heatmap_matrix <- matrix(heatmap_matrix, ncol = length(time_points))


pattern_names <- names(gene_patterns)
heatmap_data <- as.data.frame(heatmap_matrix)
colnames(heatmap_data) <- time_points 
heatmap_data$Pattern <- pattern_names 


heatmap_data_long <- melt(heatmap_data, id.vars = "Pattern")  # 转换为长格式
colnames(heatmap_data_long) <- c("Pattern", "Time", "Value")

heatmap_data_long$Value <- factor(heatmap_data_long$Value, levels = c(1, 0, -1), labels = c("red", "gray", "blue"))


pattern_df <- data.frame(
  Pattern = pattern_names,
  Count = sapply(gene_patterns, function(x) x$count)
)


ordered_patterns <- pattern_df$Pattern[order(-pattern_df$Count)]
heatmap_data_long$Pattern <- factor(heatmap_data_long$Pattern, levels = ordered_patterns)

bar_plot <- ggplot(pattern_df, aes(x = reorder(Pattern, -Count), y = Count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_minimal() +
  labs(title = "Gene Expression Patterns in PBMC", x = "Expression Pattern", y = "Number of Genes") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
        plot.margin = unit(c(1,1,-0.5,1), "lines"))  # 调整 margin 使得图形对齐


heatmap_plot <- ggplot(heatmap_data_long, aes(x = Pattern, y = Time)) +
  geom_tile(aes(fill = Value), color = "white") +
  scale_fill_manual(values = c("red" = "red", "blue" = "blue", "gray" = "gray")) +  # 自定义颜色
  theme_minimal() +
  labs(y = "Time Points", x = "Expression Pattern") +
  theme(panel.grid = element_blank(), legend.position = "none", 
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        plot.margin = unit(c(-0.5,1,1,1), "lines"))  # 调整 margin 使得图形对齐
library(cowplot)  # 用于精确对齐

aligned_plots <- plot_grid(bar_plot, heatmap_plot, ncol = 1, align = "v", rel_heights = c(2, 1))


ggsave("~/data/nongdu/剔除样本/RNAseq/figure/差异基因IDE2/figure/PBMC_combined_plot.png", plot = aligned_plots, width = 20, height = 6, dpi = 300)







# 创建柱状图
pattern_df <- data.frame(
  Pattern = names(gene_patterns),
  Count = sapply(gene_patterns, function(x) x$count)
)


plot <- ggplot(pattern_df, aes(x = Pattern, y = Count)) +
  geom_bar(stat = "identity", fill = "steelblue") +  
  theme_minimal() +  
  labs(title = "Gene Expression Patterns in PBMC",  
       x = "Expression Pattern", 
       y = "Number of Genes") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

print(plot)


pattern_df <- pattern_df[order(-pattern_df$Count), ]

plot <- ggplot(pattern_df, aes(x = reorder(Pattern, -Count), y = Count)) + 
  geom_bar(stat = "identity", fill = "steelblue") +  
  theme_minimal() +  
  labs(title = "Gene Expression Patterns in PBMC",  
       x = "Expression Pattern", 
       y = "Number of Genes") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 


print(plot)

# 保存图表到文件
output_plot_path <- "~/data/nongdu/剔除样本/RNAseq/figure/差异基因IDE2/figure/PBMC_gene_patterns_plot.png"
ggsave(output_plot_path, plot = plot, width = 20, height = 6)



# 去除含有 '0' 的模式
pattern_df <- pattern_df[!grepl("0", pattern_df$Pattern), ]


pattern_df <- pattern_df[order(-pattern_df$Count), ]


plot <- ggplot(pattern_df, aes(x = reorder(Pattern, -Count), y = Count)) +  # 按数量降序排序
  geom_bar(stat = "identity", fill = "steelblue") +  # 创建柱状图，填充色为蓝色
  theme_minimal() +  # 使用简洁主题
  labs(title = "Gene Expression Patterns in PBMC (No Zeros)",  # 设置标题和轴标签
       x = "Expression Pattern", 
       y = "Number of Genes") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # X 轴文字旋转45度以避免重叠


print(plot)

# 保存图表到文件
output_plot_path <- "~/data/nongdu/剔除样本/RNAseq/figure/差异基因IDE2/PBMC_gene_patterns_no_zeros_plot.png"
ggsave(output_plot_path, plot = plot, width = 8, height = 6)



library(ggplot2)
library(reshape2)
library(gridExtra)

heatmap_matrix <- do.call(rbind, lapply(names(gene_patterns), function(pattern) {
  sapply(strsplit(pattern, "")[[1]], function(char) {
    if (char == "+") {
      return(1)  # 红色
    } else if (char == "-") {
      return(-1)  # 蓝色
    } else {
      return(0)  # 灰色
    }
  })
}))


time_points <- c("d0.25", "d0.5", "d1", "d2", "d3", "d5")
heatmap_matrix <- matrix(heatmap_matrix, ncol = length(time_points))

pattern_names <- names(gene_patterns)
heatmap_data <- as.data.frame(heatmap_matrix)
colnames(heatmap_data) <- time_points  # 每一列对应一个时间点
heatmap_data$Pattern <- pattern_names  # 添加模式名称

heatmap_data_long <- melt(heatmap_data, id.vars = "Pattern")  # 转换为长格式
colnames(heatmap_data_long) <- c("Pattern", "Time", "Value")


heatmap_data_long$Value <- factor(heatmap_data_long$Value, levels = c(1, 0, -1), labels = c("red", "gray", "blue"))

pattern_df <- data.frame(
  Pattern = pattern_names,
  Count = sapply(gene_patterns, function(x) x$count)
)

ordered_patterns <- pattern_df$Pattern[order(-pattern_df$Count)]
heatmap_data_long$Pattern <- factor(heatmap_data_long$Pattern, levels = ordered_patterns)


bar_plot <- ggplot(pattern_df, aes(x = reorder(Pattern, -Count), y = Count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_minimal() +
  labs(title = "Gene Expression Patterns in PBMC", x = "Expression Pattern", y = "Number of Genes") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
        plot.margin = unit(c(1,1,-0.5,1), "lines"))  # 调整 margin 使得图形对齐


heatmap_plot <- ggplot(heatmap_data_long, aes(x = Pattern, y = Time)) +
  geom_tile(aes(fill = Value), color = "white") +
  scale_fill_manual(values = c("red" = "red", "blue" = "blue", "gray" = "gray")) +  # 自定义颜色
  theme_minimal() +
  labs(y = "Time Points", x = "Expression Pattern") +
  theme(panel.grid = element_blank(), legend.position = "none", 
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        plot.margin = unit(c(-0.5,1,1,1), "lines"))  # 调整 margin 使得图形对齐
library(cowplot)  # 用于精确对齐

aligned_plots <- plot_grid(bar_plot, heatmap_plot, ncol = 1, align = "v", rel_heights = c(2, 1))


aligned_plots


ggsave("~/data/nongdu/剔除样本/RNAseq/figure/差异基因IDE2/figure/PBMC_combined_plot.png", plot = aligned_plots, width = 20, height = 6, dpi = 300)



library(reshape2)
library(ggplot2)
library(cowplot)
library(readxl)

# 读取Mouse.MitoCarta3.0的基因列表
dbs <- read_excel("/home/gongfengcz/scRNA-heart-mitochodria/data/Mouse.MitoCarta3.0.xls")
# 根据需要，可能需要对dbs数据框进行适当的调整
gene_df <- dbs %>%
  mutate(gene = strsplit(Genes, ", ")) %>%
  unnest(cols = gene) # 确保这个步骤按照数据结构来调整



counts <- read.table(
  file = "/home/gongfengcz/RNAseq/data/脓毒症/GSE224127_1_counts_LPS_timecourse.txt",
  row.names = 1, header = TRUE
)
#交集
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

organs <- c("BM", "TH", "SP", "iLN", "PBMC", "BR", "HE", "LU", "LI", "KI", "SI", "CO", "SK")
p_value_threshold <- 0.05
padj_threshold <- 0.05

for (organ in organs) {
  organ_samples <- grep(paste0("^", organ), colnames(aggregated_counts_matrix), value = TRUE)
  organ_counts_matrix <- aggregated_counts_matrix[, organ_samples]
  
  organ_genes <- read.csv(paste0("~/data/nongdu/剔除样本/RNAseq/figure/差异基因IDE2/ImpulseDE2_", organ, "_results.csv"))
  
  significant_genes <- organ_genes[organ_genes$p < p_value_threshold & organ_genes$padj < padj_threshold, ]
  
  diff_gene_names <- significant_genes$Gene
  
  organ_diff_counts <- organ_counts_matrix[rownames(organ_counts_matrix) %in% diff_gene_names, ]
  
  time_points <- gsub(paste0(organ, "_"), "", organ_samples)
  sorted_time_points <- time_points[order(as.numeric(sub("d", "", time_points)))]
  
  organ_counts_matrix <- organ_counts_matrix[, paste0(organ, "_", sorted_time_points)]
  
  gene_patterns <- list()
  
  for (i in 1:nrow(organ_diff_counts)) {
    gene <- organ_diff_counts[i, ]
    pattern <- ""
    
    for (j in 2:length(sorted_time_points)) {
      if (gene[j] > gene[j-1]) {
        pattern <- paste0(pattern, "+")
      } else if (gene[j] < gene[j-1]) {
        pattern <- paste0(pattern, "-")
      } else {
        pattern <- paste0(pattern, "0")
      }
    }
    
    if (!pattern %in% names(gene_patterns)) {
      gene_patterns[[pattern]] <- list(genes = c(), count = 0)
    }
    gene_patterns[[pattern]]$genes <- c(gene_patterns[[pattern]]$genes, rownames(organ_diff_counts)[i])
    gene_patterns[[pattern]]$count <- gene_patterns[[pattern]]$count + 1
  }
  
  heatmap_matrix <- do.call(rbind, lapply(names(gene_patterns), function(pattern) {
    sapply(strsplit(pattern, "")[[1]], function(char) {
      if (char == "+") {
        return(1)
      } else if (char == "-") {
        return(-1)
      } else {
        return(0)
      }
    })
  }))
  
  time_points <- c("d0.25", "d0.5", "d1", "d2", "d3", "d5")
  heatmap_matrix <- matrix(heatmap_matrix, ncol = length(time_points))
  
  pattern_names <- names(gene_patterns)
  heatmap_data <- as.data.frame(heatmap_matrix)
  colnames(heatmap_data) <- time_points
  heatmap_data$Pattern <- pattern_names
  
  heatmap_data_long <- melt(heatmap_data, id.vars = "Pattern")
  colnames(heatmap_data_long) <- c("Pattern", "Time", "Value")
  
  heatmap_data_long$Value <- factor(heatmap_data_long$Value, levels = c(1, 0, -1), labels = c("red", "gray", "blue"))
  
  pattern_df <- data.frame(
    Pattern = pattern_names,
    Count = sapply(gene_patterns, function(x) x$count)
  )
  
  ordered_patterns <- pattern_df$Pattern[order(-pattern_df$Count)]
  heatmap_data_long$Pattern <- factor(heatmap_data_long$Pattern, levels = ordered_patterns)
  
  bar_plot <- ggplot(pattern_df, aes(x = reorder(Pattern, -Count), y = Count)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    theme_minimal() +
    labs(title = paste0("Gene Expression Patterns in ", organ), x = "Expression Pattern", y = "Number of Genes") +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
          plot.margin = unit(c(1,1,-0.5,1), "lines"))
  
  heatmap_plot <- ggplot(heatmap_data_long, aes(x = Pattern, y = Time)) +
    geom_tile(aes(fill = Value), color = "white") +
    scale_fill_manual(values = c("red" = "red", "blue" = "blue", "gray" = "gray")) +
    theme_minimal() +
    labs(y = "Time Points", x = "Expression Pattern") +
    theme(panel.grid = element_blank(), legend.position = "none", 
          axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          plot.margin = unit(c(-0.5,1,1,1), "lines"))
  
  aligned_plots <- plot_grid(bar_plot, heatmap_plot, ncol = 1, align = "v", rel_heights = c(2, 1))
  
  output_path <- paste0("~/data/nongdu/剔除样本/RNAseq/figure/差异基因IDE2/figure/filtered_", organ, "_combined_plot.png")
  ggsave(output_path, plot = aligned_plots, width = 10, height = 6, dpi = 300)
}







counts <- read.table(
  file = "/home/gongfengcz/RNAseq/data/脓毒症/GSE224127_1_counts_LPS_timecourse.txt",
  row.names = 1, header = TRUE
)

#交集
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

organs <- c("BM", "TH", "SP", "iLN", "PBMC", "BR", "HE", "LU", "LI", "KI", "SI", "CO", "SK")
p_value_threshold <- 0.05
padj_threshold <- 0.05

# 创建一个空列表用于保存所有器官的基因模式
all_gene_patterns <- list()

for (organ in organs) {
  organ_samples <- grep(paste0("^", organ), colnames(aggregated_counts_matrix), value = TRUE)
  organ_counts_matrix <- aggregated_counts_matrix[, organ_samples]
  
  organ_genes <- read.csv(paste0("~/data/nongdu/剔除样本/RNAseq/figure/差异基因IDE2/ImpulseDE2_", organ, "_results.csv"))
  
  significant_genes <- organ_genes[organ_genes$p < p_value_threshold & organ_genes$padj < padj_threshold, ]
  
  diff_gene_names <- significant_genes$Gene
  
  organ_diff_counts <- organ_counts_matrix[rownames(organ_counts_matrix) %in% diff_gene_names, ]
  
  time_points <- gsub(paste0(organ, "_"), "", organ_samples)
  sorted_time_points <- time_points[order(as.numeric(sub("d", "", time_points)))]
  
  organ_counts_matrix <- organ_counts_matrix[, paste0(organ, "_", sorted_time_points)]
  
  # 计算每个基因的表达模式
  for (i in 1:nrow(organ_diff_counts)) {
    gene <- organ_diff_counts[i, ]
    pattern <- ""
    
    for (j in 2:length(sorted_time_points)) {
      if (gene[j] > gene[j-1]) {
        pattern <- paste0(pattern, "+")
      } else if (gene[j] < gene[j-1]) {
        pattern <- paste0(pattern, "-")
      } else {
        pattern <- paste0(pattern, "0")
      }
    }
    
    # 将模式加入到全局基因模式列表中
    if (!pattern %in% names(all_gene_patterns)) {
      all_gene_patterns[[pattern]] <- list(genes = c(), count = 0)
    }
    all_gene_patterns[[pattern]]$genes <- c(all_gene_patterns[[pattern]]$genes, rownames(organ_diff_counts)[i])
    all_gene_patterns[[pattern]]$count <- all_gene_patterns[[pattern]]$count + 1
  }
}

# 生成热图矩阵
heatmap_matrix <- do.call(rbind, lapply(names(all_gene_patterns), function(pattern) {
  sapply(strsplit(pattern, "")[[1]], function(char) {
    if (char == "+") {
      return(1)
    } else if (char == "-") {
      return(-1)
    } else {
      return(0)
    }
  })
}))

time_points <- c("d0.25", "d0.5", "d1", "d2", "d3", "d5")
heatmap_matrix <- matrix(heatmap_matrix, ncol = length(time_points))

pattern_names <- names(all_gene_patterns)
heatmap_data <- as.data.frame(heatmap_matrix)
colnames(heatmap_data) <- time_points
heatmap_data$Pattern <- pattern_names

heatmap_data_long <- melt(heatmap_data, id.vars = "Pattern")
colnames(heatmap_data_long) <- c("Pattern", "Time", "Value")

# 重新编码数值，1 对应红色，-1 对应蓝色，0 对应灰色
heatmap_data_long$Value <- factor(heatmap_data_long$Value, levels = c(1, 0, -1), labels = c("red", "gray", "blue"))

# 创建柱状图，显示每个模式的基因数量
pattern_df <- data.frame(
  Pattern = pattern_names,
  Count = sapply(all_gene_patterns, function(x) x$count)
)

# 确保热图的 Pattern 顺序和柱状图一致
ordered_patterns <- pattern_df$Pattern[order(-pattern_df$Count)]
heatmap_data_long$Pattern <- factor(heatmap_data_long$Pattern, levels = ordered_patterns)

# 创建柱状图
bar_plot <- ggplot(pattern_df, aes(x = reorder(Pattern, -Count), y = Count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_minimal() +
  labs(title = "Gene Expression Patterns Across Organs", x = "Expression Pattern", y = "Number of Genes") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
        plot.margin = unit(c(1,1,-0.5,1), "lines"))

# 创建热图
heatmap_plot <- ggplot(heatmap_data_long, aes(x = Pattern, y = Time)) +
  geom_tile(aes(fill = Value), color = "white") +
  scale_fill_manual(values = c("red" = "red", "blue" = "blue", "gray" = "gray")) +
  theme_minimal() +
  labs(y = "Time Points", x = "Expression Pattern") +
  theme(panel.grid = element_blank(), legend.position = "none", 
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        plot.margin = unit(c(-0.5,1,1,1), "lines"))

library(cowplot)


aligned_plots <- plot_grid(bar_plot, heatmap_plot, ncol = 1, align = "v", rel_heights = c(2, 1))


ggsave("~/data/nongdu/剔除样本/RNAseq/figure/差异基因IDE2/figure/filtered_All_organs_combined_plot.png", plot = aligned_plots, width = 20, height = 6, dpi = 300)








counts <- read.table(
  file = "/home/gongfengcz/RNAseq/data/脓毒症/GSE224127_1_counts_LPS_timecourse.txt",
  row.names = 1, header = TRUE
)

sample_info <- colnames(counts)
sample_time <- sapply(strsplit(sample_info, "_"), function(x) paste(x[1], x[2], sep = "_"))

aggregated_counts <- list()

for (time in unique(sample_time)) {
  matching_samples <- grep(paste0("^", time), sample_info, value = TRUE)
  aggregated_counts[[time]] <- rowSums(counts[, matching_samples])
}

aggregated_counts_matrix <- do.call(cbind, aggregated_counts)
colnames(aggregated_counts_matrix) <- unique(sample_time)

organs <- c("BM", "TH", "SP", "iLN", "PBMC", "BR", "HE", "LU", "LI", "KI", "SI", "CO", "SK")
p_value_threshold <- 0.05
padj_threshold <- 0.05

# 初始化全局 gene_patterns 列表
all_gene_patterns <- list()

for (organ in organs) {
  organ_samples <- grep(paste0("^", organ), colnames(aggregated_counts_matrix), value = TRUE)
  organ_counts_matrix <- aggregated_counts_matrix[, organ_samples]
  
  organ_genes <- read.csv(paste0("~/data/nongdu/剔除样本/RNAseq/figure/差异基因IDE2/ImpulseDE2_", organ, "_results.csv"))
  
  significant_genes <- organ_genes[organ_genes$p < p_value_threshold & organ_genes$padj < padj_threshold, ]
  
  diff_gene_names <- significant_genes$Gene
  
  organ_diff_counts <- organ_counts_matrix[rownames(organ_counts_matrix) %in% diff_gene_names, ]
  
  time_points <- gsub(paste0(organ, "_"), "", organ_samples)
  sorted_time_points <- time_points[order(as.numeric(sub("d", "", time_points)))]
  
  organ_counts_matrix <- organ_counts_matrix[, paste0(organ, "_", sorted_time_points)]
  
  # 保存每个器官的 gene_patterns
  gene_patterns <- list()
  
  for (i in 1:nrow(organ_diff_counts)) {
    gene <- organ_diff_counts[i, ]
    pattern <- ""
    
    for (j in 2:length(sorted_time_points)) {
      if (gene[j] > gene[j-1]) {
        pattern <- paste0(pattern, "+")
      } else if (gene[j] < gene[j-1]) {
        pattern <- paste0(pattern, "-")
      } else {
        pattern <- paste0(pattern, "0")
      }
    }
    
    if (!pattern %in% names(gene_patterns)) {
      gene_patterns[[pattern]] <- list(genes = c(), count = 0)
    }
    gene_patterns[[pattern]]$genes <- c(gene_patterns[[pattern]]$genes, rownames(organ_diff_counts)[i])
    gene_patterns[[pattern]]$count <- gene_patterns[[pattern]]$count + 1
  }
  
  # 将每个器官的基因模式存入 all_gene_patterns
  all_gene_patterns[[organ]] <- gene_patterns
}

# 现在 all_gene_patterns 中包含每个器官的基因表达模式
save(all_gene_patterns, file = "~/data/nongdu/剔除样本/RNAseq/figure/差异基因IDE2/filtered_all_gene_patterns.RData")





load("~/data/nongdu/剔除样本/RNAseq/figure/差异基因IDE2/filtered_all_gene_patterns.RData")

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

query_gene_patterns("PBMC")


query_gene_patterns("PBMC", pattern = "++-")




query_gene_patterns("PBMC", pattern = "+--+-+",assign_var = "genes_negative_pattern")




#提取模式为 -00000 的基因名称
genes_negative_pattern <- gene_patterns[["-00000"]]$genes

selected_genes_diff_counts <- pbmc_diff_counts[rownames(pbmc_diff_counts) %in% genes_negative_pattern, ]




# 绘制密度图
pbmc_values <- as.vector(as.matrix(pbmc_diff_counts))

value_df <- data.frame(Value = pbmc_values)

value_count <- table(pbmc_values)


print(value_count)
# 检查 pbmc_values 的基础统计信息
summary(pbmc_values)

density_plot <- ggplot(value_df, aes(x = Value)) +
  geom_density(fill = "lightgreen", color = "black", alpha = 0.6) +  # 绘制密度图
  theme_minimal() +  # 使用简洁主题
  labs(title = "Density Plot of PBMC Diff Counts", 
       x = "Counts", 
       y = "Density")


print(density_plot)
output_density_plot_path <- "~/data/nongdu/剔除样本/RNAseq/figure/差异基因IDE2/figur/PBMC_diff_counts_density_plot.png"
ggsave(output_density_plot_path, plot = density_plot, width = 8, height = 6)




