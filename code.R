#对照组数量
control <- 3
#实验组数量
treat <- 3
#二选一
input <- "Expression profiling by array"
# input <- "Expression profiling by high throughput sequencing"
# 获取当前工作目录
current_dir <- getwd()
# 将工作目录设置为当前目录（尽管这没有实际改变，但这是您要求的操作）
setwd(current_dir)
# 验证工作目录
print(getwd())
#输入文件准备----
library(readxl)
library(dplyr)
gpl <- read_excel("GPL.xlsx")
rawdata <- read_excel("GSE14561.xlsx")
mergedata <- merge(gpl, rawdata, by = "ID")
mergedata <- mergedata[, -1]
colnames(mergedata)[1] <- "gene_name"
#不同行取平均值
mergedata <- mergedata %>%
  group_by(gene_name) %>%
  summarise(across(everything(), mean))
write.table(mergedata, "1.rawdata.xls", quote = FALSE,
            sep = "\t", row.names = FALSE)
# 2. 判断是否需要进行ID转换（小鼠与人类的转换、人类的各种ID之间的转换）----
change <- read.table("1.rawdata.xls", header = TRUE, sep = "\t")
#判断第一列元素中有小写字母的占所有元素的比例，达到50%以上可以认为是小鼠
# 提取第一列
first_column <- change[[1]]
# 判断每个元素是否包含小写字母
contains_lowercase <- grepl("[a-z]", first_column)
# 计算包含小写字母的元素所占的比例
lowercase_ratio <- sum(contains_lowercase) / length(first_column)
# 判断是否达到50%以上
if (lowercase_ratio > 0.5) {
  print("小鼠数据")
  #ID转化小鼠的gene_symbol转换为人类的gene_symbol
  library(clusterProfiler)
  library(dplyr)
  library(org.Mm.eg.db) # 小鼠数据库
  library(biomaRt)
  human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl",
                   host = "https://dec2021.archive.ensembl.org/")
  mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl",
                   host = "https://dec2021.archive.ensembl.org/")
  hsa2mus_all <- getLDS(attributes = c("mgi_symbol"),
                        filters = "mgi_symbol",
                        values = change$gene_name,
                        mart = mouse,
                        attributesL = c("hgnc_symbol"),
                        martL = human,
                        uniqueRows = TRUE)
  merged_data <- merge(hsa2mus_all, change,
                       by.x = "MGI.symbol", by.y = "gene_name")
  merged_data <- merged_data[, -1]
  colnames(merged_data)[1] <- "gene_name"
  #不同行取平均值
  merged_data <- merged_data %>%
    group_by(gene_name) %>%
    summarise(across(everything(), mean))
  write.table(merged_data, "1.rawdata.xls", quote = FALSE,
              sep = "\t", row.names = FALSE)
  print("小鼠gene_symbol已经转换为人类gene_symbol")
} else if (all(grepl("^\\d+$", first_column))) {
  #都是数字
  print("是 ENTREZ_ID 数据")
  library("org.Hs.eg.db")
  library(dplyr)
  # 假设第一列是Entrez ID
  entrez_ids <- first_column
  # 使用biomaRt进行ID转换
  human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  # 进行Entrez ID到基因符号的转换
  converted_ids <- getBM(attributes = c("entrezgene_id", "hgnc_symbol"),
                         filters = "entrezgene_id",
                         values = entrez_ids,
                         mart = human,
                         uniqueRows = TRUE)
  merged_data <- merge(converted_ids, change,
                       by.x = "entrezgene_id", by.y = "gene_name")
  merged_data <- merged_data[, -1]
  colnames(merged_data)[1] <- "gene_name"
  #不同行取平均值
  merged_data <- merged_data %>%
    group_by(gene_name) %>%
    summarise(across(everything(), mean))
  write.table(merged_data, "1.rawdata.xls", quote = FALSE,
              sep = "\t", row.names = FALSE)
} else if (all(grepl("^ENS", first_column))) {
  #以ENS开头
  print("是 ENSEMBL_ID 数据")
  library(openxlsx)
  library(clusterProfiler)
  library(org.Hs.eg.db) # 人类数据库
  # 提取基因的 ENSEMBL ID，去除版本号，并转换为向量
  gene_ens_id <- sapply(change$gene_name, function(x) gsub("\\..*", "", x))
  # 转换 ENSEMBL ID 到基因符号
  gene_symbol <- bitr(geneID = gene_ens_id,
                      fromType = "ENSEMBL",
                      toType = c("SYMBOL"),
                      OrgDb = org.Hs.eg.db)
  merged_data <- merge(gene_symbol, change,
                       by.x = "ENSEMBL", by.y = "gene_name")
  merged_data <- merged_data[, -1]
  colnames(merged_data)[1] <- "gene_name"
  #不同行取平均值
  merged_data <- merged_data %>%
    group_by(gene_name) %>%
    summarise(across(everything(), mean))
  write.table(merged_data, "1.rawdata.xls", quote = FALSE,
              sep = "\t", row.names = FALSE)
} else {
  print("已经是gene_symbol数据")
}
# 1. 差异分析：对于芯片数据和高通量测序数据分不同的执行代码----
if (input == "Expression profiling by array") {
  print("数据为芯片数据")
  library(gplots)
  library(limma)
  library(impute)
  raw_data <- read.table("1.rawdata.xls", sep = "\t", header = TRUE)
  colnames(raw_data)[1] <- "gene_name"
  #不同行取平均值
  raw_data <- raw_data %>%
    group_by(gene_name) %>%
    summarise(across(everything(), mean))
  raw_data <- as.matrix(raw_data)
  rownames(raw_data) <- raw_data[, 1]
  exp <- raw_data[, 2:ncol(raw_data)]
  dimnames <- list(rownames(exp), colnames(exp))
  exp <- matrix(as.numeric(as.matrix(exp)), nrow = nrow(exp),
                dimnames = dimnames)
  # impute missing expression data
  mat <- impute.knn(exp)
  raw_data <- mat$data
  raw_data <- avereps(raw_data)
  # differential
  sample_class <- c(rep("nom", control), rep("lch", treat))
  design <- model.matrix(~0 + factor(sample_class))
  colnames(design) <- c("nom", "lch")
  fit <- lmFit(raw_data, design)
  contrast_matrix <- makeContrasts(nom - lch, levels = design)
  fit2 <- contrasts.fit(fit, contrast_matrix)
  fit2 <- eBayes(fit2)
  all_diff <- topTable(fit2, adjust = "fdr", number = 200000)
  # 将行名转换为第一列
  all_diff <- cbind(gene_name = rownames(all_diff), all_diff)
  rownames(all_diff) <- NULL
  # 修改列名
  colnames(all_diff)[1] <- "gene_name"
  # 保存数据框到文件
  write.table(all_diff, file = "2.diffALL.xls", sep = "\t",
              quote = FALSE, row.names = FALSE)
  # 设置调整后的p值
  adjust_p <- 0.05
  # 定义一个向量，包含我们想要的log_fold_change值
  log_fold_changes <- c(0, 0.58, 1, 2)
  # 遍历每个log_fold_change值
  for (i in seq_along(log_fold_changes)) {
    log_fold_change <- log_fold_changes[i]
    # 筛选差异显著的基因
    diff_sig <- all_diff[with(all_diff, (abs(logFC) > log_fold_change &
                                           P.Value < adjust_p)), ]
    # 保存差异显著的基因到文件
    write.table(diff_sig, file = paste0("3.", format(log_fold_change,
                                                     digits = 2), "_diff.xls"),
                sep = "\t", quote = FALSE, row.names = FALSE)
    # 保存上调基因到文件
    diff_up <- all_diff[with(all_diff, (logFC > log_fold_change &
                                          P.Value < adjust_p)), ]
    write.table(diff_up, file = paste0("4.", format(log_fold_change,
                                                    digits = 2), "_up.xls"),
                sep = "\t", quote = FALSE, row.names = FALSE)
    # 保存下调基因到文件
    diff_down <- all_diff[with(all_diff, (logFC < -log_fold_change &
                                            P.Value < adjust_p)), ]
    write.table(diff_down, file = paste0("5.", format(log_fold_change,
                                                      digits = 2), "_down.xls"),
                sep = "\t", quote = FALSE, row.names = FALSE)
  }
} else if (input == "Expression profiling by high throughput sequencing") {
  print("数据为高通量数据")
  library(DESeq2)
  library(openxlsx)
  library(limma)
  library(readxl)
  ####差异分析
  rawdata <- read.table("1.rawdata.xls", header = TRUE, sep = "\t")
  colnames(rawdata)[1] <- "gene_name"
  rawdata <- rawdata %>%
    group_by(gene_name) %>%
    summarise(across(everything(), mean))
  data <- rawdata[, -1]
  rawdata <- avereps(data, ID = rawdata$gene_name) # 自定义
  rawdata <- rawdata[rowMeans(rawdata) > 1, ] # 自定义
  diffcount <- rawdata[, seq_len(ncol(rawdata))]
  ncontrol <- control
  ntreat <- treat
  sample_table <- data.frame(condition = factor(c(rep("C", ncontrol),
                                                  rep("T", ntreat))))
  ##可以依据实际情况修改各个条件的名称及数量
  #构建dds矩阵
  dds <- DESeqDataSetFromMatrix(countData = round(diffcount),
                                colData = sample_table, design = ~condition)
  #对原始dds进行normalize
  dds <- DESeq(dds)
  res <- results(dds)# 将结果用results()函数来获取，赋值给res变量
  #保存全部差异结果
  all_diff <- as.data.frame(res)
  all_diff$gene_id <- rownames(all_diff)
  all_diff <- all_diff[, colnames(all_diff)[c(7, 1:6)]]
  all_diff <- na.omit(all_diff)
  # 重命名列名
  colnames(all_diff)[colnames(all_diff) == "log2FoldChange"] <- "logFC"
  colnames(all_diff)[colnames(all_diff) == "pvalue"] <- "P.Value"
  write.table(all_diff, file = "2.diffALL.xls",
              sep = "\t", row.names = FALSE, quote = FALSE)
  # 设置调整后的p值
  adjust_p <- 0.05
  # 定义一个向量，包含我们想要的log_fold_change值
  log_fold_changes <- c(0, 0.58, 1, 2)
  # 遍历每个log_fold_change值
  for (i in seq_along(log_fold_changes)) {
    log_fold_change <- log_fold_changes[i]
    # 筛选差异显著的基因
    diff_sig <- all_diff[with(all_diff, (abs(logFC) > log_fold_change &
                                           P.Value < adjust_p)), ]
    # 保存差异显著的基因到文件
    write.table(diff_sig, file = paste0("3.", format(log_fold_change,
                                                     digits = 2), "_diff.xls"),
                sep = "\t", quote = FALSE, row.names = FALSE)
    # 保存上调基因到文件
    diff_up <- all_diff[with(all_diff, (logFC > log_fold_change &
                                          P.Value < adjust_p)), ]
    write.table(diff_up, file = paste0("4.", format(log_fold_change,
                                                    digits = 2), "_up.xls"),
                sep = "\t", quote = FALSE, row.names = FALSE)
    # 保存下调基因到文件
    diff_down <- all_diff[with(all_diff, (logFC < -log_fold_change &
                                            P.Value < adjust_p)), ]
    write.table(diff_down, file = paste0("5.", format(log_fold_change,
                                                      digits = 2), "_down.xls"),
                sep = "\t", quote = FALSE, row.names = FALSE)
  }
}
#火山图和热图
library(ggplot2)
library(gplots)
# 假设 all_diff 是你的原始差异表达数据框，raw_data 是你的原始表达矩阵
# 设置调整后的p值
adjust_p <- 0.05
# 定义log_fold_change的值
log_fold_changes <- c(0, 0.58, 1, 2)
# 确保all_diff中没有NA值
all_diff[is.na(all_diff)] <- 0
# 遍历每个log_fold_change值
for (i in seq_along(log_fold_changes)) {
  log_fold_change <- log_fold_changes[i]
  # 火山图绘制
  # 为当前log_fold_change添加change列
  all_diff$change <- with(all_diff, ifelse(
    P.Value < adjust_p & abs(logFC) >= abs(log_fold_change),
    ifelse(logFC > log_fold_change, "Up", "Down"),
    "Stable"
  ))
  # 设置火山图文件名
  volcano_filename <- paste0("volcano_", log_fold_change, ".tiff")
  # 绘制火山图
  p <- ggplot(all_diff, aes(x = logFC, y = -log10(P.Value), colour = change)) +
    geom_point(alpha = 0.4, size = 3.5) +
    scale_color_manual(values = c("#546de5", "#d2dae2", "#ff4757")) +
    geom_vline(xintercept = c(-log_fold_change, log_fold_change),
               lty = 4, col = "black", lwd = 0.8) +
    geom_hline(yintercept = -log10(adjust_p), lty = 4,
               col = "black", lwd = 0.8) +
    labs(x = "log2(fold change)", y = "-log10 (pvalue)") +
    theme_bw() +
    theme(legend.position = "right", legend.title = element_blank())
  # 保存火山图为TIFF文件
  tiff(file = volcano_filename, width = 10, height = 8,
       units = "in", res = 300, compression = "lzw")
  print(p)
  dev.off()
  # 清除change列
  all_diff$change <- NULL
  # 筛选差异表达的基因
  diff_sig <- all_diff[with(all_diff, (P.Value < adjust_p &
                                         abs(logFC) >= abs(log_fold_change))), ]
  # 如果diff_sig不为空，则绘制热图
  if (nrow(diff_sig) > 0) {
    # 确保diff_sig中的基因在raw_data的行名中
    common_genes <- intersect(diff_sig$gene_name, rownames(raw_data))
    if (length(common_genes) != nrow(diff_sig)) {
      warning("Some genes in diff_sig are not present in raw_data.
              Using common genes for heatmap.")
      diff_sig <- diff_sig[common_genes, ]
    }
    # 提取表达数据并进行对数变换
    heatmap_data <- raw_data[common_genes, ]
    hm_exp <- log10(heatmap_data + 0.001)
    hm_mat <- as.matrix(hm_exp)
    # 设置热图文件名
    heatmap_filename <- paste0("heatmap_", log_fold_change, ".tiff")
    # 绘制热图
    tiff(file = heatmap_filename, width = 10, height = 10, units = "in",
         res = 300, compression = "lzw")
    par(oma = c(10, 3, 3, 7)) # 设置图形边距
    heatmap.2(hm_mat, col = bluered(75), trace = "none")
    dev.off()
  } else {
    message("当前log_fold_change值为", log_fold_change, "时，没有差异表达的基因，跳过热图绘制。")
  }
}
# 7. 进行差异基因富集分析----
library(org.Hs.eg.db)
library(clusterProfiler)
library(topGO)
library(DOSE)
library(ggplot2)
library(readxl)
library(ggalluvial)
library(dplyr)
library(stringr)
# 定义要分析的文件名前缀和后缀
file_prefix <- "3."
file_suffix <- "_diff.xls"
files <- c("0", "0.58", "1", "2") # 文件名中间的部分
# 自定义保存TIFF文件的函数
save_tiff <- function(plot, filename, width = 25,
                      height = 25, units = "cm", res = 300) {
  tiff(filename, width = width, height = height, units = units,
       res = res, compression = "lzw")
  print(plot)
  dev.off()
}
# 循环处理每个文件
for (file_mid in files) {
  file_name <- paste0(file_prefix, file_mid, file_suffix)
  # 读取差异表达基因列表
  rt <- read.table(file_name, header = TRUE, sep = "\t")
  gene_names <- rt$gene_name
  # 将基因名转换为ENTREZID
  genes <- bitr(gene_names, fromType = "SYMBOL", toType = "ENTREZID",
                OrgDb = org.Hs.eg.db)
  # 处理1:多映射问题，只保留每个SYMBOL的第一个ENTREZID
  genes <- genes[!duplicated(genes$SYMBOL), ]
  gene_entrezids <- genes$ENTREZID
  # GO富集分析
  bp <- enrichGO(gene_entrezids, OrgDb = org.Hs.eg.db,
                 ont = "BP", pvalueCutoff = 0.05,
                 pAdjustMethod = "bonferroni", qvalueCutoff = 0.2,
                 readable = TRUE)
  cc <- enrichGO(gene_entrezids, OrgDb = org.Hs.eg.db,
                 ont = "CC", pvalueCutoff = 0.05,
                 pAdjustMethod = "bonferroni", qvalueCutoff = 0.2,
                 readable = TRUE)
  mf <- enrichGO(gene_entrezids, OrgDb = org.Hs.eg.db,
                 ont = "MF", pvalueCutoff = 0.05,
                 pAdjustMethod = "bonferroni", qvalueCutoff = 0.2,
                 readable = TRUE)
  all <- enrichGO(gene_entrezids, OrgDb = org.Hs.eg.db,
                  ont = "ALL", pvalueCutoff = 0.05,
                  pAdjustMethod = "bonferroni", qvalueCutoff = 0.2,
                  readable = TRUE)
  filename <- paste0("6.DEGs_GO_", file_mid, ".xls")
  # 将数据框写入文件
  write.table(all, file = filename, row.names = FALSE, sep = "\t",
              quote = FALSE, fileEncoding = "UTF-8")
  # 读取并绘制GO富集结果
  go_file <- paste0("6.DEGs_GO_", file_mid, ".xls")
  if (file.exists(go_file)) {
    go <- tryCatch(
      {
        read.delim(go_file, sep = "\t", header = TRUE, check.names = FALSE)
      },
      error = function(e) {
        message("Error reading file: ", go_file)
        return(NULL)
      }
    )
    if (!is.null(go)) {
      go_sort <- go[order(go$ONTOLOGY, -as.numeric(go$Count)), ]
      # 获取各个分类的前10个条目
      top_m <- head(go_sort[go_sort$ONTOLOGY == "MF", ],
                    min(10, nrow(go_sort[go_sort$ONTOLOGY == "MF", ])))
      top_c <- head(go_sort[go_sort$ONTOLOGY == "CC", ],
                    min(10, nrow(go_sort[go_sort$ONTOLOGY == "CC", ])))
      top_b <- head(go_sort[go_sort$ONTOLOGY == "BP", ],
                    min(10, nrow(go_sort[go_sort$ONTOLOGY == "BP", ])))
      slimgo <- rbind(top_b, top_c, top_m)
      slimgo$Description <- factor(slimgo$Description,
                                   levels = slimgo$Description)
      p <- ggplot(data = slimgo, mapping = aes(x = Description,
                                               y = Count, fill = ONTOLOGY))
      p + geom_bar(stat = "identity")
      pp <- p + geom_bar(stat = "identity") + coord_flip()
      p <- p + geom_bar(stat = "identity") + coord_flip() +
        scale_x_discrete(limits = rev(levels(slimgo$Description)))
      pp
      tiff_filename <- paste0("Enrichment_DEGs_GO_", file_mid, ".tiff")
      save_tiff(pp, tiff_filename)
    }
  } else {
    message("File does not exist: ", go_file)
  }
  # KEGG富集分析
  kk <- enrichKEGG(gene = gene_entrezids, organism = "hsa",
                   pvalueCutoff = 0.05, qvalueCutoff = 0.2)
  data <- setReadable(kk, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  # 保存KEGG富集结果
  write.table(data, paste0("7.DEGs_KEGG_", file_mid, ".xls"), sep = "\t",
              quote = FALSE, row.names = FALSE)
  # 创建KEGG富集图
  barplot_kk <- barplot(kk, drop = TRUE, showCategory = 20, title = "pathway")
  tiff_filename_kegg <- paste0("Enrichment_DEGs_KEGG_", file_mid, ".tiff")
  save_tiff(barplot_kk, tiff_filename_kegg)
  # 桑基图：基因——通路——二级代谢通路
  file_path <- paste0("7.DEGs_KEGG_", file_mid, ".xls")
  data <- read.table(file_path, header = TRUE, sep = "\t")
  # 筛选出pvalue小于0.05的行，并按照Count值降序排列
  filtered_data <- data %>%
    filter(pvalue < 0.05) %>%
    arrange(desc(Count))
  # 选出ID、subcategory、category三列
  data <- filtered_data[, c("ID", "subcategory", "category")]
  data <- data[seq_len(nrow(data)), ]
  # 定义函数来拆分超过20个字符的标签
  split_long_labels <- function(label, max_length = 20) {
    if (nchar(label) > max_length) {
      label <- str_wrap(label, width = max_length, exdent = 0)
    }
    return(label)
  }
  # 应用标签拆分函数
  data$ID <- sapply(data$ID, split_long_labels)
  data$subcategory <- sapply(data$subcategory, split_long_labels)
  data$category <- sapply(data$category, split_long_labels)
  # 转换数据为lodes形式
  cor_lodes <- to_lodes_form(data, axes = seq_len(ncol(data)), id = "Cohort")
  # 获取stratum的唯一值数量
  num_strata <- length(unique(cor_lodes$stratum))
  # 定义足够多的颜色
  mycol <- colorRampPalette(c(
    "#029149", "#E0367A", "#D8D155", "#D20A13", "#91612D",
    "#FFD121", "#088247", "#11AA4D", "#58CDD9", "#5D90BA", "#7CC767"
  ))(num_strata)
  # 绘制桑基图
  plot_sankey <- ggplot(cor_lodes, aes(x = x, stratum = stratum,
                                       alluvium = Cohort,
                                       fill = stratum, label = stratum)) +
    scale_x_discrete(expand = c(0.1, 0.1)) +
    geom_flow(width = 0.4, aes.flow = "forward") +
    geom_stratum(alpha = 0.9, width = 0.4) +
    scale_fill_manual(values = mycol) +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)),
              size = 5, color = "black") +
    xlab("") +
    ylab("") +
    theme_bw() +
    theme(
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      axis.text.y = element_blank(),
      axis.text.x = element_text(size = 26),
      panel.grid = element_blank(),
      panel.border = element_blank()
    ) +
    ggtitle("") +
    guides(fill = "none")
  # 使用save_tiff函数保存文件
  sankey_tiff_filename <- paste0("sankey_ggalluvial_", file_mid, ".tiff")
  save_tiff(plot_sankey, sankey_tiff_filename)
}
# 8. 进行GSEA基因集富集分析----
#GSEA基因集富集分析
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(enrichplot)
library(DOSE)
library(gridExtra)
# 读取数据
data <- read.table("2.diffALL.xls", header = TRUE, sep = "\t")
colnames(data)[1] <- "gene_name"
# 添加 ENTREZID 列并按 log2FoldChange 排序
ids_symbol_and_entrezid <- bitr(geneID = data$gene_name,
                                fromType = "SYMBOL",
                                toType = "ENTREZID",
                                OrgDb = org.Hs.eg.db)
data$ENTREZID <- ids_symbol_and_entrezid[match(data$gene_name,
                                               ids_symbol_and_entrezid$SYMBOL),
                                         2]
data <- na.omit(data)
# 处理重复的 log2FoldChange 值:logFC扰动
set.seed(123)
data$logFC <- jitter(data$logFC)
data <- data[order(data$logFC, decreasing = TRUE), ]
# 构建包含 ENTREZID 名称的 log2FoldChange 向量
gene_list <- data$logFC
names(gene_list) <- data$ENTREZID
# 使用 gseKEGG 进行 KEGG 通路富集分析
gsea_result <- gseKEGG(geneList = gene_list,
                       organism = "hsa",
                       minGSSize = 10,
                       maxGSSize = 10000,
                       pvalueCutoff = 0.5,
                       verbose = FALSE,
                       eps = 0)
# 将 GSEA 结果转换为数据框
gsea_result_df <- as.data.frame(gsea_result)
# 过滤 P 值小于 0.05 的结果
filtered_results <- gsea_result_df[gsea_result_df$pvalue < 0.05, ]
write.table(filtered_results, "8.GSEA_KEGG_result.xls",
            quote = FALSE, sep = "\t")
# 按 NES 值降序排列正的 NES 值并选择前五个
positive_nes <- filtered_results[filtered_results$NES > 0, ]
positive_nes <- positive_nes[order(positive_nes$NES, decreasing = TRUE), ]
top5_positive_nes <- head(positive_nes, 5)
# 按 NES 值升序排列负的 NES 值并选择前五个
negative_nes <- filtered_results[filtered_results$NES < 0, ]
negative_nes <- negative_nes[order(negative_nes$NES, decreasing = FALSE), ]
top5_negative_nes <- head(negative_nes, 5)
# 检查是否找到感兴趣的通路
if (nrow(top5_positive_nes) == 0 && nrow(top5_negative_nes) == 0) {
  stop("没有找到符合条件的通路")
}
# 循环遍历前五个正的 NES 通路
for (i in 1:5) {
  # 可视化每个通路
  plot <- gseaplot2(gsea_result,
                    geneSetID = top5_positive_nes$ID[i],
                    color = "red",
                    base_size = 11,
                    rel_heights = c(1.5, 0.5, 1),
                    subplots = 1:3,
                    pvalue_table = TRUE,
                    ES_geom = "line")
  # 生成文件名
  file_name <- paste0("GSEA_positive_NES_pathway_", i, ".tiff")
  # 保存为 TIFF 格式文件
  ggsave(file_name, plot = plot, device = "tiff",
         width = 10, height = 8, units = "in", dpi = 300)
}
# 循环遍历前五个正的 NES 通路
for (i in 1:5) {
  # 可视化每个通路
  plot <- gseaplot2(gsea_result,
                    geneSetID = top5_negative_nes$ID[i],
                    color = "blue",
                    base_size = 11,
                    rel_heights = c(1.5, 0.5, 1),
                    subplots = 1:3,
                    pvalue_table = TRUE,
                    ES_geom = "line")
  # 生成文件名
  file_name <- paste0("GSEA_negative_NES_pathway_", i, ".tiff")
  # 保存为 TIFF 格式文件
  ggsave(file_name, plot = plot, device = "tiff",
         width = 10, height = 8, units = "in", dpi = 300)
}
# GSEA结果做基因棒棒糖图
#棒棒糖图：反应显著的通路的正负向调节
library(ggplot2)
library(dplyr)
library(tidyr)
library(cols4all)
library(stringr)
# 读取 GSEA 数据文件
dat <- read.table("8.GSEA_KEGG_result.xls", header = TRUE, sep = "\t")
# 根据NES确定通路的方向
dat <- dat %>%
  mutate(Direction = ifelse(NES > 0, "Up", "Down"))
# 分别挑选上调和下调通路的前10个，按Count排序
top_up <- dat %>%
  filter(Direction == "Up") %>%
  arrange(desc(setSize)) %>%
  head(10)
top_down <- dat %>%
  filter(Direction == "Down") %>%
  arrange(desc(setSize)) %>%
  head(10)
# 合并上调和下调的前10个通路
top_paths <- bind_rows(top_up, top_down)
# 为了使得 up 向右延伸，down 向左延伸，我们将 NES 进行调整
top_paths <- top_paths %>%
  mutate(adjusted_NES = ifelse(Direction == "Up", NES, NES))
# 将长标签分成两行
top_paths <- top_paths %>%
  mutate(Description = str_wrap(Description, width = 30)) # 这里假设宽度为30，根据实际情况调整
# 创建图形
p1 <- ggplot(top_paths, aes(x = reorder(Description,
                                        adjusted_NES), y = adjusted_NES)) +
  coord_flip() +
  geom_col(aes(fill = Direction), width = 0.05) + # 修改 width 参数调整柱子的宽度
  geom_point(aes(size = setSize, color = p.adjust)) + # 添加散点/气泡，颜色表示p值的大小
  scale_size_continuous(range = c(4, 10)) +
  scale_color_gradient2(midpoint = 0.05,
                        low = "#16888b", mid = "#c4d3d2", high = "green",
                        limits = c(0, 0.05), breaks = seq(0, 0.05, 0.01),
                        labels = seq(0, 0.05, 0.01)) + # 手动设置颜色渐变
  scale_fill_manual(values = c("Up" = "#e3c4a5",
                               "Down" = "#97a2dc")) + # 手动设置填充颜色
  geom_hline(yintercept = 0, color = "white",
             size = 0.5, lty = "dashed") +
  geom_vline(xintercept = 0, color = "white",
             size = 0.5, lty = "dashed") + # 添加垂直线分割象限
  labs(x = NULL, y = "NES",
       title = "KEGG Pathway Enrichment",
       fill = "Direction", color = "p.adjust") +
  scale_y_continuous(expand = expansion(add = c(0.1, 0.1)),
                     limits = c(-4, 4),
                     breaks = seq(-4, 4, 1), #根据调整后NES值修改X轴刻度
                     labels = c(-4, -3, -2, -1, 0, 1, 2, 3, 4))
# 添加主题调整
p1 <- p1 + theme_bw() +
  theme(axis_text_y <- element_text(size = 15,
                                    hjust = 1, angle = 20)) + # 调整 y 轴标签大小和位置
  scale_x_discrete(labels = function(x) str_wrap(x, width = 30)) + # 将 x 轴标签换行
  theme(axis.text.x = element_text(size = 8,
                                   angle = 0,
                                   vjust = 1,
                                   hjust = 1)) + # 调整 x 轴刻度标签
  guides(size = guide_legend(title = "setSize"),
         fill = guide_legend(title = "Direction"))
# 保存图形为 TIFF 格式
ggsave(filename = "BBT_GSEA_result.tiff", plot = p1,
       device = "tiff", dpi = 300, width = 7, height = 10)
# 3. 进行TPM转换----
#####看每一列和是不是1×10的六次方，不是的话就转换
library(GEOquery)
library(tidyr)
library(dplyr)
library(limma)
library(openxlsx)
library(WGCNA)
library(biomaRt)
# 需要注意使用TPM格式数据做，而不是count格式
# 查看基因组参数
mart <- useMart("ensembl")
# 以人类为例 获取基因组信息
bmart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                 dataset = "hsapiens_gene_ensembl",
                 verbose = TRUE)
file_path <- "1.rawdata.xls"
data <- read.table(file_path, header = TRUE, sep = "\t")
# 检查数据列和是否为1*10^6  # nolint
column_sums <- colSums(data[, -1])
if (all(abs(column_sums - 1e6) < 1e-2)) {
  cat("数据列和已经为1*10^6，不需要重新计算TPM\n")
  tpm <- read.table("1.rawdata.xls", header = TRUE, sep = "\t")
  write.table(tpm, "1.TPM_rawdata.xls", quote = FALSE,
              sep = "\t", row.names = FALSE)
} else {
  cat("数据列和不是1*10^6,重新计算TPM\n")
  temp_data <- data[2:ncol(data)]
  # 使用apply函数，获取重复基因的均值
  temp_list <- apply(temp_data, 2, function(x) {
                                                aggregate(
                                                          x ~ data$gene_name,
                                                          data = temp_data,
                                                          mean)})
  # 对结果进行merge
  temp_merge <- c()
  for (i in seq_along(temp_list)) {
    temp_merge <- cbind(temp_merge, temp_list[[i]]$x)
  }
  temp_df <- data.frame(temp_merge)
  colnames(temp_df) <- names(temp_list)
  temp_df$gene_symbol <- temp_list[[1]][[1]]
  # 获取无重复geneSymbol的表达矩阵
  rt <- temp_df
  exp <- rt[, c(ncol(rt), seq_len(ncol(rt)) - 1)]
  exp_matrix <- exp
  feature_ids <- exp_matrix[, 1]
  attributes <- c("ensembl_gene_id", "hgnc_symbol",
                  "chromosome_name", "start_position",
                  "end_position")
  filters <- "hgnc_symbol"
  feature_info <- biomaRt::getBM(attributes = attributes,
                                 filters = filters, values = feature_ids,
                                 mart = bmart)
  special_string <- "CHR_"
  feature_info_1 <- feature_info[!grepl(special_string,
                                        feature_info$chromosome_name), ]
  feature_info_1 <- feature_info_1[!duplicated(feature_info_1$hgnc_symbol), ]
  # 计算基因的有效长度eff_length
  eff_length <- abs(feature_info_1$end_position - feature_info_1$start_position)
  eff_length <- as.data.frame(eff_length)
  rownames(eff_length) <- feature_info_1$hgnc_symbol
  # 将eff_length的行名变成第一列
  eff_length <- data.frame(hgnc_symbol = rownames(eff_length),
                           eff_length, row.names = NULL)
  eff_length <- eff_length[, c(2, 1)]
  colnames(eff_length) <- c("eff_length", "gene_id")
  rownames(eff_length) <- eff_length$gene_id
  rownames(eff_length) <- do.call(rbind,
                                  strsplit(eff_length$gene_id, "\\."))[, 1]
  feature_ids <- exp_matrix[, 1]
  if (!all(feature_ids %in% rownames(eff_length))) {
    tbl <- table(feature_ids %in% rownames(eff_length))
    msg1 <- sprintf("%i gene is shared, %i gene is specified",
                    tbl[[1]], tbl[[2]])
    warning(msg1)
  }
  if (!identical(feature_ids, rownames(eff_length))) {
    msg2 <- sprintf("Given GTF file only contain %i gene,
                    but expression matrix has %i gene",
                    nrow(eff_length), nrow(exp_matrix))
    warning(msg2)
  }
  # 修剪表达矩阵和有效基因长度
  exp_matrix <- exp_matrix[feature_ids %in% rownames(eff_length), ]
  rownames(exp_matrix) <- exp_matrix[, 1]
  exp_matrix <- exp_matrix[, -1]
  mm <- match(rownames(exp_matrix), rownames(eff_length))
  eff_length <- eff_length[mm, ]
  if (identical(rownames(eff_length), rownames(exp_matrix))) {
    print("GTF和表达矩阵现在有相同的基因,并且基因顺序相同")
  }
  # 计算TPM
  x <- exp_matrix / eff_length$eff_length
  exp_matrix_tpm <- t(t(x) / colSums(x)) * 1e6
  # 检查每列的总和是否为1
  colSums(exp_matrix_tpm)
  exp_matrix_tpm_1 <- as.data.frame(exp_matrix_tpm)
  exp_matrix_tpm_1$gene_symbol <- rownames(exp_matrix_tpm_1)
  exp_matrix_tpm_1 <- exp_matrix_tpm_1[, c(ncol(exp_matrix_tpm_1),
                                           seq_len(ncol(exp_matrix_tpm_1)) - 1)]
  data <- exp_matrix_tpm_1
  temp_data <- data[2:ncol(data)]
  temp_list <- apply(temp_data, 2, function(x) {
                                                aggregate(x ~ data$gene_symbol,
                                                          data = temp_data,
                                                          mean)})
  temp_merge <- c()
  for (i in seq_along(temp_list)) {
    temp_merge <- cbind(temp_merge, temp_list[[i]]$x)
  }
  temp_df <- data.frame(temp_merge)
  colnames(temp_df) <- names(temp_list)
  temp_df$gene_symbol <- temp_list[[1]][[1]]
  exp_data_all <- temp_df
  exp_data_all <- exp_data_all[, c(ncol(exp_data_all),
                                   seq_len(ncol(exp_data_all)) - 1)]
  rownames(exp_data_all) <- exp_data_all[, 1]
  exp_data_all <- exp_data_all[, -1]
  exp_data_all_1 <- exp_data_all[rowSums(exp_data_all) > 10, ]
  # 将exp_data_all_1的行名变成第一列
  exp_data_all_1 <- data.frame(hgnc_symbol = rownames(exp_data_all_1),
                               exp_data_all_1, row.names = NULL)
  colnames(exp_data_all_1)[1] <- "gene_name"
  write.table(exp_data_all_1, "1.TPM_rawdata.xls", quote = FALSE,
              sep = "\t", row.names = FALSE)
}
# 10. 进行WGCNA分析----
######判断软阈值是否在合理范围
######不在合理范围内就把之前的文件都删除，包括转成功的TPM格式文件
library("GEOquery")
library(tidyr)
library(dplyr)
library(limma)
library(openxlsx)
library("WGCNA")
file_path <- "1.TPM_rawdata.xls"
exp <- read.table(file_path, header = TRUE, sep = "\t")
rownames(exp) <- exp[, 1]
exp <- exp[, -1]
exp <- as.data.frame(t(exp))
sample_tree <- hclust(dist(exp), method = "average")
# 设置TIFF文件参数
tiff(file = "WGCNA_1_sample_clustering.tiff", width = 12, height = 8,
     units = "in", res = 300, compression = "lzw")
# 设置图形参数
par(cex = 0.9)
par(mar = c(0, 4, 2, 0))
# 绘制样本聚类图
plot(sample_tree, main = "Sample clustering to detect outliers",
     sub = "", xlab = "", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
# 关闭TIFF设备
dev.off()
#cutreeStatic模块检测在分层树状图使用恒定高度树切割。只有其大小至少为minSize的分支被保留。
clust <- cutreeStatic(sample_tree, cutHeight = 6000,
                      minSize = 5)
table(clust)
powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))
#多软阈值函数的无尺度拓扑分析。目的是为了帮助用户选择一个合适的软阈值权值用于网络构建。
sft <- pickSoftThreshold(exp, powerVector = powers, verbose = 5)
# 假设sft是WGCNA pickSoftThreshold函数的输出
# 首先，找到SFT.R.sq列中的最大值和对应的power
max_r_squared <- max(sft$fitIndices[, "SFT.R.sq"])
max_power <- sft$fitIndices[which.max(sft$fitIndices[, "SFT.R.sq"]), "Power"]
# 根据提供的规则选择软阈值
if (max_r_squared > 0.9) {
  # 如果存在大于0.9的值，找到第一个大于0.9的值对应的power
  threshold_index <- which(sft$fitIndices[, "SFT.R.sq"] > 0.9)[1]
  selected_power <- sft$fitIndices[threshold_index, "Power"]
} else if (max_r_squared > 0.8 && max_r_squared <= 0.9) {
  # 如果没有大于0.9的值，但是有大于0.8的值，选择最大的纵坐标对应的横坐标
  selected_power <- max_power
} else {
  # 如果所有值都小于0.8，输出无法进行WGCNA
  selected_power <- NA
  stop("无法进行WGCNA,因为SFT.R.sq的最大值小于0.8")
}
# 输出选择的软阈值
if (!is.na(selected_power)) {
  print(paste("选择的软阈值为:", selected_power))
}
cut <- selected_power
# 设置TIFF文件参数
tiff(file = "WGCNA_2_soft_threshold.tiff", width = 10, height = 8,
     units = "in", res = 300, compression = "lzw")
# 绘制第一个图表：Scale-free topology fit index
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit,signed R^2", type = "n",
     main = paste("Scale independence"))
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = 1.5, col = "red")
# 添加参考线
abline(h = 0.9, col = "red")
dev.off()
# 绘制第二个图表：Mean connectivity
tiff(file = "WGCNA_3_Mean_connectivity.tiff", width = 10, height = 8,
     units = "in", res = 300, compression = "lzw")
plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[, 1], sft$fitIndices[, 5],
     labels = powers, cex = 1.5, col = "red")
# 关闭TIFF设备
dev.off()
dat_expr <- as.data.frame(exp)
net <- blockwiseModules(dat_expr, power = cut,
                        TOMType = "unsigned", minModuleSize = 30,
                        reassignThreshold = 0, mergeCutHeight = 0.25,
                        numericLabels = TRUE, pamRespectsDendro = FALSE,
                        saveTOMs = TRUE,
                        saveTOMFileBase = "MedBioInfoCloud_TOM",
                        verbose = 3)
merged_colors <- labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
# 设置TIFF文件参数
tiff(file = "WGCNA_4_clustering_tree.tiff", width = 10, height = 8,
     units = "in", res = 300, compression = "lzw")
# 绘制聚类树图
plotDendroAndColors(net$dendrograms[[1]], merged_colors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
# 关闭TIFF设备
dev.off()
module_labels <- net$colors
module_colors <- labels2colors(net$colors)
table(module_labels)
mes <- net$MEs
gene_tree <- net$dendrograms[[1]]
save(mes, module_labels, module_colors, gene_tree,
     file = "MedBioInfoCloud-auto.RData")
# 定义基因和样本的数量
n_genes <- ncol(dat_expr)
n_samples <- nrow(dat_expr)
# 读取数据
clinical <- read.table("1.TPM_rawdata.xls", header = TRUE, sep = "\t")
# 提取第二列到最后一列的列名
colname <- colnames(clinical)[2:ncol(clinical)]
colname <- data.frame(sampleID = colname)
# 获取行数
nrow_colname <- nrow(colname)
# 创建 control 列，填充 control 个 1 和 nrow(colname) - control 个 0
colname$control <- rep(0, nrow_colname)
colname$control[1:control] <- 1
# 创建 treat 列，填充 nrow(colname) - treat 个 0 和 treat 个 1
colname$treat <- rep(0, nrow_colname)
colname$treat[(nrow_colname - treat + 1):nrow_colname] <- 1
py <- colname
# ___________________________________
#|___sampleID__|__control__|__treat__|
#|___GSMXXXXX__|_____1_____|____0____|
#|___GSMXXXXX__|_____1_____|____0____|
#|___GSMXXXXX__|_____1_____|____0____|
#|___GSMXXXXX__|_____0_____|____1____|
#|___GSMXXXXX__|_____0_____|____1____|
#|___GSMXXXXX__|_____0_____|____1____|
rownames(py) <- py[, 1]
py <- py[, -1]
# Recalculate MEs with color labels
mes0 <- moduleEigengenes(dat_expr, module_colors)$eigengenes
mes <- orderMEs(mes0)
module_trait_cor <- cor(mes, py, use = "p")
module_trait_cor <- data.frame(module_trait_cor)
module_trait_cor <- module_trait_cor[module_trait_cor$treat > 0.5 |
                                       module_trait_cor$control > 0.5, ]
model <- rownames(module_trait_cor)
module_trait_cor <- as.matrix(module_trait_cor)
module_trait_pvalue <- corPvalueStudent(module_trait_cor, n_samples)
mes <- mes[, model]
text_matrix <- paste(signif(module_trait_cor, 2), "\n(",
                     signif(module_trait_pvalue, 0.05), ")", sep = "")
dim(text_matrix) <- dim(module_trait_cor)
# 设置TIFF文件参数
tiff(file = "WGCNA_5_relative_heatmap.tiff", width = 10, height = 8,
     units = "in", res = 300, compression = "lzw")
# 设置图形参数
par(mar = c(6, 8.5, 3, 3))
# 绘制热图
labeledHeatmap(Matrix = module_trait_cor,
               xLabels = names(py),
               yLabels = names(mes),
               ySymbols = names(mes),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = text_matrix,
               setStdMargins = FALSE,
               cex.text = 1,
               zlim = c(-1, 1),
               main = paste("Module-trait relationships"))
# 关闭TIFF设备
dev.off()
interested_module <- py$treat
interested_module <- as.data.frame(interested_module)
names(interested_module) <- "InterestedModule"
min_pvalue_rowname <- rownames(
    module_trait_pvalue)[ # nolint
        which.min( # nolint
            as.matrix( # nolint
                module_trait_pvalue))]
cleaned_rowname <- sub("^ME", "", min_pvalue_rowname)
#module
module <- cleaned_rowname
mod_names <- substring(names(mes), 3)
gene_module_membership <- as.data.frame(cor(dat_expr, mes, use = "p"))
mm_pvalue <- as.data.frame(corPvalueStudent(as.matrix(gene_module_membership),
                                            n_samples))
names(gene_module_membership) <- paste("MM", mod_names, sep = "")
names(mm_pvalue) <- paste("p.MM", mod_names, sep = "")
gene_traitsignificance <- as.data.frame(cor(dat_expr,
                                            interested_module, use = "p"))
gs_pvalue <- as.data.frame(corPvalueStudent(as.matrix(gene_traitsignificance),
                                            n_samples))
names(gene_traitsignificance) <- paste("GS", names(interested_module), sep = "")
names(gs_pvalue) <- paste("p.GS", names(interested_module), sep = "")
dat_expr <- as.data.frame(dat_expr)
names(dat_expr)
names(dat_expr)[module_colors == module]
table(module_colors == module)
column <- match(module, mod_names)
module_colors <- merged_colors
module_genes <- module_colors == module
all_genes <- t(dat_expr[module_colors == module])
all_genes <- all_genes[, -c(seq_len(ncol(all_genes)))]
all_genes <- data.frame(hgnc_symbol = rownames(all_genes),
                        all_genes, row.names = NULL)
colnames(all_genes)[1] <- "gene_name"
filename <- paste0("9.WGCNA_result_genes_", cleaned_rowname, ".xls")
# 将数据框写入文件
write.table(all_genes, file = filename, row.names = FALSE,
            quote = FALSE, fileEncoding = "UTF-8")
# 设置TIFF文件参数
tiff(file = "WGCNA_plots.tiff", width = 15, height = 10,
     units = "in", res = 300, compression = "lzw")
# 设置图形参数
par(mfrow = c(1, 1))
# 绘制散点图
verboseScatterplot(abs(gene_module_membership[module_genes, column]),
                   abs(gene_traitsignificance[module_genes, 1]),
                   xlab = "x",
                   ylab = "y", main = "title",
                   cex.main = 1.2, cex.axis = 1.2, cex.lab = 1.2,
                   col = "blue",
                   abline = FALSE,
                   lmFnc = lm,
                   pch = 1)
# 添加参考线
abline(h = 0.8, v = 0.8, col = "red", lty = 1)
# 关闭TIFF设备
dev.off()
# 12. 将WGCNA模块基因和差异基因取交集----
#WGCNA和差异基因取交集，核心基因
# 读取文件并重命名列名
module_gene <- all_genes
colnames(module_gene) <- "gene_name"
# 数据预处理：去除前后空格，统一大小写
module_gene$gene_name <- toupper(trimws(module_gene$gene_name))
# 定义文件名前缀和后缀
file_prefix <- "3."
file_suffix <- "_diff.xls"
files <- c("0", "0.58", "1", "2") # 文件名中间的部分
# 循环处理每个差异基因文件
for (file_mid in files) {
  diff_file_name <- paste0(file_prefix, file_mid, file_suffix)
  # 读取差异基因文件
  diff <- read.table(diff_file_name, header = TRUE, sep = "\t")
  diff <- diff[, 1]
  diff <- as.data.frame(diff)
  colnames(diff) <- "gene_name"
  # 数据预处理：去除前后空格，统一大小写
  diff$gene_name <- toupper(trimws(diff$gene_name))
  # 取交集
  intersection_column <- intersect(module_gene$gene_name, diff$gene_name)
  intersection_column_df <- as.data.frame(intersection_column)
  colnames(intersection_column_df) <- "gene_name"
  # 查看交集的前几行
  print(head(intersection_column_df))
  # 构建输出文件名
  output_file_name <- paste0("10.WGCNA_diff_intersection_", file_mid, ".xls")
  # 保存交集结果
  write.table(intersection_column_df, output_file_name, sep = "\t",
              quote = FALSE, row.names = FALSE)
}
#如何通过R语言得到cytoHubba
library(igraph)
library(STRINGdb)
library(Cairo)
# 定义保存数据框的函数
save_dataframe <- function(filename, df) {
  write.table(df, file = filename, sep = "\t", quote = FALSE, row.names = FALSE)
}
# 定义一个函数来计算中心性指标并保存核心基因
calculate_core_genes <- function(file_mid) {
  # 确保所有图形设备都被关闭
  graphics.off()
  # 读取基因文件
  gene_file_name <- paste0("10.WGCNA_diff_intersection_", file_mid, ".xls")
  gene_names_df <- read.table(gene_file_name, header = TRUE, sep = "\t")
  gene_names <- gene_names_df$gene_name
  # 初始化 STRINGdb 对象
  string_db <- STRINGdb$new(version = "11.0", species = 9606,
                            score_threshold = 400, input_directory = "")
  # 映射基因名到 STRING ID
  mapped_genes <- string_db$map(data.frame(gene = gene_names),
                                "gene", removeUnmappedRows = TRUE)
  # 检查未成功映射的基因
  unmapped_genes <- setdiff(gene_names, mapped_genes$gene)
  if (length(unmapped_genes) > 0) {
    warning("未能映射的基因数: ", length(unmapped_genes))
    print("未能映射的基因: ")
    print(unmapped_genes)
  }
  # 获取相互作用信息
  interactions <- string_db$get_interactions(mapped_genes$STRING_id)
  # 将 STRING ID 映射回基因名
  string_to_gene <- setNames(mapped_genes$gene, mapped_genes$STRING_id)
  interactions$from_gene_name <- string_to_gene[as.character(interactions$from)]
  interactions$to_gene_name <- string_to_gene[as.character(interactions$to)]
  interaction_genes_df <- data.frame(from = interactions$from_gene_name,
                                     to = interactions$to_gene_name)
  # 创建 igraph 对象
  g <- graph_from_data_frame(d = interaction_genes_df, directed = FALSE)
  # 计算中心性指标
  degree_centrality <- degree(g)
  mnc_centrality <- calculate_mnc(g)
  dmnc_centrality <- calculate_dmnc(g)
  # 创建中心性指标数据框
  centrality_df <- data.frame(
    gene = names(degree_centrality),
    Degree = degree_centrality,
    MNC = mnc_centrality,
    DMNC = dmnc_centrality
  )
  # 移除含有 NA 值的行
  centrality_df <- na.omit(centrality_df)
  # 保存所有中心性指标的基因数据框
  output_file_name <- paste0("11.centrality_", file_mid, ".xls")
  save_dataframe(output_file_name, centrality_df)
  # 取前100个或全部基因
  top_genes <- function(df, metric) {
    df <- df[order(-df[[metric]]), ]
    if (nrow(df) > 100) {
      return(df$gene[1:100])
    } else {
      return(df$gene)
    }
  }
  top_degree_genes <- top_genes(centrality_df, "Degree")
  top_mnc_genes <- top_genes(centrality_df, "MNC")
  top_dmnc_genes <- top_genes(centrality_df, "DMNC")
  # 创建一个包含所有中心性指标前100基因的列表
  gene_lists <- list(
    Degree = top_degree_genes,
    MNC = top_mnc_genes,
    DMNC = top_dmnc_genes
  )
  # 取所有中心性指标的交集
  core_genes <- Reduce(intersect, gene_lists)
  # 保存核心基因
  core_genes_df <- data.frame(gene = core_genes)
  core_output_file_name <- paste0("11.core_genes_", file_mid, ".xls")
  save_dataframe(core_output_file_name, core_genes_df)
}
# 定义计算MNC的函数
calculate_mnc <- function(graph) {
  mnc <- sapply(V(graph), function(v) {
    neighbors <- neighbors(graph, v)
    subgraph <- induced_subgraph(graph, neighbors)
    components <- components(subgraph)
    max(components$csize)
  })
  return(mnc)
}
# 定义计算DMNC的函数
calculate_dmnc <- function(graph) {
  dmnc <- sapply(V(graph), function(v) {
    neighbors <- neighbors(graph, v)
    subgraph <- induced_subgraph(graph, neighbors)
    components <- components(subgraph)
    max_component <- which.max(components$csize)
    max_subgraph <- induced_subgraph(subgraph,
                                     V(subgraph)[components$membership ==
                                                   max_component])
    ecount(max_subgraph) / (vcount(max_subgraph) *
                              (vcount(max_subgraph) - 1) / 2)
  })
  return(dmnc)
}
# 定义文件名前缀和后缀
files <- c("0", "0.58", "1", "2") # 文件名中间的部分
# 循环处理每个文件
for (file_mid in files) {
  calculate_core_genes(file_mid)
}
# 14. 核心基因富集分析----
#进行富集分析
library(org.Hs.eg.db)
library(clusterProfiler)
library(topGO)
library(DOSE)
library(ggplot2)
library(readxl)
library(ggalluvial)
library(dplyr)
library(stringr)
# 定义要分析的文件名前缀和后缀
file_prefix <- "11.core_genes_"
file_suffix <- ".xls"
files <- c("0", "0.58", "1", "2") # 文件名中间的部分
# 自定义保存TIFF文件的函数
save_tiff <- function(plot, filename, width = 25,
                      height = 25, units = "cm", res = 300) {
  tiff(filename, width = width, height = height, units = units,
       res = res, compression = "lzw")
  print(plot)
  dev.off()
}
# 循环处理每个文件
for (file_mid in files) {
  tryCatch({
    file_name <- paste0(file_prefix, file_mid, file_suffix)
    # 读取差异表达基因列表
    rt <- read.table(file_name, header = TRUE, sep = "\t")
    colnames(rt)[1] <- "gene_name"
    gene_names <- rt$gene_name
    # 将基因名转换为ENTREZID
    genes <- bitr(gene_names, fromType = "SYMBOL", toType = "ENTREZID",
                  OrgDb = org.Hs.eg.db)
    # 处理1:多映射问题，只保留每个SYMBOL的第一个ENTREZID
    genes <- genes[!duplicated(genes$SYMBOL), ]
    gene_entrezids <- genes$ENTREZID
    # GO富集分析
    bp <- enrichGO(gene_entrezids, OrgDb = org.Hs.eg.db,
                   ont = "BP", pvalueCutoff = 0.05,
                   pAdjustMethod = "bonferroni", qvalueCutoff = 0.2,
                   readable = TRUE)
    cc <- enrichGO(gene_entrezids, OrgDb = org.Hs.eg.db,
                   ont = "CC", pvalueCutoff = 0.05,
                   pAdjustMethod = "bonferroni", qvalueCutoff = 0.2,
                   readable = TRUE)
    mf <- enrichGO(gene_entrezids, OrgDb = org.Hs.eg.db,
                   ont = "MF", pvalueCutoff = 0.05,
                   pAdjustMethod = "bonferroni", qvalueCutoff = 0.2,
                   readable = TRUE)
    all <- enrichGO(gene_entrezids, OrgDb = org.Hs.eg.db,
                    ont = "ALL", pvalueCutoff = 0.05,
                    pAdjustMethod = "bonferroni", qvalueCutoff = 0.2,
                    readable = TRUE)
    filename <- paste0("12.core_genes_GO_", file_mid, ".xls")
    # 将数据框写入文件
    write.table(all, file = filename, row.names = FALSE, sep = "\t",
                quote = FALSE, fileEncoding = "UTF-8")
    # 读取并绘制GO富集结果
    go_file <- paste0("12.core_genes_GO_", file_mid, ".xls")
    if (file.exists(go_file)) {
      go <- tryCatch(
        {
          read.delim(go_file, sep = "\t", header = TRUE, check.names = FALSE)
        },
        error = function(e) {
          message("Error reading file: ", go_file)
          return(NULL)
        }
      )
      if (!is.null(go)) {
        go_sort <- go[order(go$ONTOLOGY, -as.numeric(go$Count)), ]
        # 获取各个分类的前10个条目
        top_m <- head(go_sort[go_sort$ONTOLOGY == "MF", ],
                      min(10, nrow(go_sort[go_sort$ONTOLOGY == "MF", ])))
        top_c <- head(go_sort[go_sort$ONTOLOGY == "CC", ],
                      min(10, nrow(go_sort[go_sort$ONTOLOGY == "CC", ])))
        top_b <- head(go_sort[go_sort$ONTOLOGY == "BP", ],
                      min(10, nrow(go_sort[go_sort$ONTOLOGY == "BP", ])))
        slimgo <- rbind(top_b, top_c, top_m)
        slimgo$Description <- factor(slimgo$Description,
                                     levels = slimgo$Description)
        p <- ggplot(data = slimgo, mapping = aes(x = Description,
                                                 y = Count, fill = ONTOLOGY)) +
          geom_bar(stat = "identity") + coord_flip() +
          scale_x_discrete(limits = rev(levels(slimgo$Description)))
        tiff_filename <- paste0("Enrichment_core_genes_GO_", file_mid, ".tiff")
        save_tiff(p, tiff_filename)
      }
    } else {
      message("File does not exist: ", go_file)
    }
    # KEGG富集分析
    kk <- enrichKEGG(gene = gene_entrezids, organism = "hsa",
                     pvalueCutoff = 0.05, qvalueCutoff = 0.2)
    data <- setReadable(kk, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
    # 保存KEGG富集结果
    write.table(data, paste0("13.core_genes_KEGG_", file_mid, ".xls"),
                sep = "\t",
                quote = FALSE,
                row.names = FALSE)
    # 创建KEGG富集图
    barplot_kk <- barplot(kk, drop = TRUE, showCategory = 20, title = "pathway")
    tiff_filename_kegg <- paste0("Enrichment_core_genes_KEGG_",
                                 file_mid, ".tiff")
    save_tiff(barplot_kk, tiff_filename_kegg)
  }, error = function(e) {
    message("Error processing file: ", file_mid, " - ", e$message)
  })
}
# 15. 免疫浸润----
library(e1071)
library(preprocessCore)
library(parallel)
library(tidyverse)
library(ggsci)
library(tidyr)
library(ggpubr)
# 核心算法
core_alg <- function(X, y) {
  svn_itor <- 3
  res <- function(i) {
    nus <- ifelse(i == 1, 0.25, ifelse(i == 2, 0.5, 0.75))
    model <- e1071::svm(X, y, type = "nu-regression",
                        kernel = "linear", nu = nus, scale = FALSE)
    model
  }
  out <- if (Sys.info()["sysname"] == "Windows") {
    parallel::mclapply(1:svn_itor, res, mc.cores = 1)
  } else {
    parallel::mclapply(1:svn_itor, res, mc.cores = svn_itor)
  }
  nusvm <- numeric(svn_itor)
  corrv <- numeric(svn_itor)
  for (t in 1:svn_itor) {
    weights <- t(out[[t]]$coefs) %*% out[[t]]$SV
    weights[weights < 0] <- 0
    w <- weights / sum(weights)
    u <- sweep(X, 2, w, "*")
    k <- rowSums(u)
    nusvm[t] <- sqrt(mean((k - y)^2))
    corrv[t] <- cor(k, y)
  }
  mn <- which.min(nusvm)
  model <- out[[mn]]
  q <- t(model$coefs) %*% model$SV
  q[q < 0] <- 0
  w <- q / sum(q)
  list("w" = w, "mix_rmse" = nusvm[mn], "mix_r" = corrv[mn])
}
# 置换函数
do_perm <- function(perm, X, Y) {
  Ylist <- as.list(data.matrix(Y))
  dist <- numeric()
  for (itor in 1:perm) {
    yr <- as.numeric(Ylist[sample(length(Ylist), nrow(X))])
    yr <- (yr - mean(yr)) / sd(yr)
    result <- core_alg(X, yr)
    dist <- c(dist, result$mix_r)
  }
  list("dist" = dist)
}
# 主函数
CIBERSORT <- function(sig_matrix, mixture_file, perm = 0, QN = TRUE) {
  X <- read.table(sig_matrix, header = TRUE, sep = "\t",
                  row.names = 1, check.names = FALSE)
  Y <- read.table(mixture_file, header = TRUE, sep = "\t",
                  row.names = 1, check.names = FALSE)
  X <- data.matrix(X)
  Y <- data.matrix(Y)
  X <- X[order(rownames(X)), ]
  Y <- Y[order(rownames(Y)), ]
  if (max(Y) < 50) Y <- 2^Y
  if (QN) {
    tmpc <- colnames(Y)
    tmpr <- rownames(Y)
    Y <- preprocessCore::normalize.quantiles(Y)
    colnames(Y) <- tmpc
    rownames(Y) <- tmpr
  }
  common_genes <- intersect(rownames(X), rownames(Y))
  X <- X[common_genes, ]
  Y <- Y[common_genes, ]
  X <- (X - mean(X)) / sd(as.vector(X))
  if (perm > 0) nulldist <- sort(do_perm(perm, X, Y)$dist)
  header <- c("Mixture", colnames(X), "P-value", "Correlation", "RMSE")
  output <- matrix()
  mixtures <- ncol(Y)
  for (itor in 1:mixtures) {
    y <- Y[, itor]
    y <- (y - mean(y)) / sd(y)
    result <- core_alg(X, y)
    w <- result$w
    mix_r <- result$mix_r
    mix_rmse <- result$mix_rmse
    pval <- if (perm > 0) 1 - (which.min(abs(nulldist - mix_r)) /
                                 length(nulldist)) else 9999
    out <- c(colnames(Y)[itor], w, pval, mix_r, mix_rmse)
    output <- if (itor == 1) out else rbind(output, out)
  }
  write.table(rbind(header, output), file = "17.CIBERSORT-Results.txt",
              sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  obj <- rbind(header, output)[-1, -1]
  obj <- matrix(as.numeric(unlist(obj)), nrow = nrow(obj))
  rownames(obj) <- colnames(Y)
  colnames(obj) <- c(colnames(X), "P-value", "Correlation", "RMSE")
  obj
}
data <- read.table("1.rawdata.xls", sep = "\t",
                   check.names = FALSE, header = TRUE)
temp_data <- data[2:ncol(data)]
#使用apply函数，获取重复基因的均值
temp_list <- apply(temp_data, 2,
  function(x) aggregate(x ~ data$gene_name, data = temp_data, mean)
)
#对结果进行merge
temp_merge <- c()
for (i in seq_along(temp_list)){
  temp_merge <- cbind(temp_merge, temp_list[[i]]$x)
}
temp_df <- data.frame(temp_merge)
colnames(temp_df) <- names(temp_list)
temp_df$gene_symbol <- temp_list[[1]][[1]]
#获取无重复geneSymbol的表达矩阵
rt <- temp_df
length(unique(rt$gene_symbol))
rt <- rt[, c(ncol(rt), seq_len(ncol(rt)) - 1)]
write.table(rt, "1.rawdata.txt", quote = FALSE, sep = "\t", row.names = FALSE)
#构建group.txt文件
# 读取数据
group <- read.table("1.rawdata.xls", header = TRUE, sep = "\t")
# 提取第二列到最后一列的列名
rt <- colnames(group)[2:ncol(group)]
rt <- data.frame(sampleID = rt)
# 获取行数
nrow_col_name <- nrow(rt)
rt$group <- c(rep("con", control),
              rep("tre", treat),
              rep(0, nrow_col_name - control - treat))
write.table(rt, "16.group.txt", quote = FALSE,
            sep = "\t", row.names = FALSE)
# 加载必要的包
library(readxl)
# 定义文件下载的URL和目标文件名
url <- "https://static-content.springer.com/esm/art%3A10.1038%2Fnmeth.3337/MediaObjects/41592_2015_BFnmeth3337_MOESM207_ESM.xls" # nolint
destfile <- "14.Supplementary_table_1.xls"
# 下载文件
download.file(url, destfile, mode = "wb")
# 读取Excel文件的第二个工作表
matrix_data <- read_excel(destfile, sheet = "SuppTable1_GEP_Matrix")
matrix_data <- matrix_data[-c(1:13), -c(1:2)]
# 由于图片显示了一些列名，我们需要重命名这些列
colnames(matrix_data) <- c("gene_name",
                           "B cells naive",
                           "B cells memory",
                           "Plasma cells",
                           "T cells CD8",
                           "T cells CD4 naive",
                           "T cells CD4 memory resting",
                           "T cells CD4 memory activated",
                           "T cells follicular helper",
                           "T cells regulatory (Tregs)",
                           "T cells gamma delta",
                           "NK cells resting",
                           "NK cells activated",
                           "Monocytes",
                           "Macrophages M0",
                           "Macrophages M1",
                           "Macrophages M2",
                           "Dendritic cells resting",
                           "Dendritic cells activated",
                           "Mast cells resting",
                           "Mast cells activated",
                           "Eosinophils",
                           "Neutrophils")
# 定义输出文件路径
output_file <- "15.LM22.txt"
# 将数据写入制表符分割的txt文件，粘贴方式为纯文本粘贴
write.table(matrix_data, file = output_file, sep = "\t",
            row.names = FALSE, col.names = TRUE, quote = FALSE)
if (input == "Expression profiling by array") {
  result1 <- CIBERSORT("15.LM22.txt", "1.rawdata.txt",
                       perm = 1000, QN = TRUE)
} else {
  result1 <- CIBERSORT("15.LM22.txt", "1.rawdata.txt",
                       perm = 1000, QN = FALSE)
}
#perm：表示置换次数，数字越大运行时间越长，一般文章都设置为1000；
#QN：如果为芯片数据这里设为“T”；如果为测序数据设为“F”
#最后会在工作目录下自动生成结果文件“CIBERSORT-Results.txt”
#LM22.txt获取方法：
#https://www.nature.com/articles/nmeth.3337#MOESM207 下载Supplementry table 1，
#然后只将其中矩阵部分(sheet2)保留并保存为制表符分割的txt
#要求两个都是txt格式文件，可以整理下Cibersort_1.R进行修改
b <- read.table("16.group.txt", sep = "\t",
                row.names = 1, check.names = FALSE, header = TRUE)
# _______________________
#|__sample_ID__|__group__|
#|__GSM364419__|___tre___|
#|__GSM364420__|___tre___|
#|__GSM364421__|___tre___|
#|__GSM364422__|___con___|
#|__GSM364423__|___con___|
#|__GSM364424__|___con___|
a <- read.table("17.CIBERSORT-Results.txt", sep = "\t",
                row.names = 1, check.names = FALSE, header = TRUE)
a <- a[, 1:22] #去除浸润分析结果文件后三列
class(b$group)
a$group <- b$group #将分组信息添加到a后面
a <- a %>% rownames_to_column("sample") #把行名变成一列，并加上“sample”
b <- gather(a, key = CIBERSORT, value = Proportion, -c(group, sample))
# 创建一个TIFF文件来保存图像
tiff("CIBERSORT.tiff", width = 10, height = 7, units = "in", res = 300)
# 绘制箱线图
ggboxplot(b, x = "CIBERSORT", y = "Proportion",
          fill = "group", palette = "lancet") +
  stat_compare_means(aes(group = group),
                     method = "wilcox.test",
                     label = "p.signif",
                     symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                        symbols = c("***", "**", "*", "ns"))) +
  theme(text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1))
# 关闭TIFF设备
dev.off()
