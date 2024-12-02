library(DESeq2) 
library(dplyr)
library(clusterProfiler)
library(tibble)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(enrichplot)
library(tidyverse)
library(patchwork)



# GO & KEGG barchart ------------------------------------------------------


Sys.setenv(language = "en")
options(stringsAsFactors = F)

setwd("~/Documents/去掉A4转录分析")

cds <- read.csv(file = "file/AB_new.csv",row.names = 1) # 表达矩阵
colData <- read.csv(file = "file/mymetaAB_new.csv",row.names = 1) # 分组信息
annoData <- read.csv("gene_expression.addAnno.csv") |> 
  dplyr::select(Gene.id, Gene.Name)

# 表达矩阵和分组信息匹配
if(!all(rownames(colData) == colnames(cds))){
  cds <- cds[colnames(cds),]
}

# 分组
colData$dex <- factor(colData$dex,levels = c("control","KO"))


# 创建DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = cds,
                              colData = colData,
                              design = ~ dex)
# tidy = T)  # 在没有行名的情况下可以使用

# 过滤
# counts(dds) 提取原始表达矩阵
keep <- rowSums(counts(dds)) > 1  # 行和大于1的Ensembl Gene ID
dds <- dds[keep,]

# 分析
dds <- DESeq(dds)
resultsNames(dds) # 查看dds的结果名，便于指定想要提取的列
degDeSeq2 <- results(dds,name = "dex_KO_vs_control") %>% 
  as.data.frame() %>% 
  na.omit()

# 通过读取的数据将ensgene转换成symbol
degDeSeq2$Gene.id <- rownames(degDeSeq2)
degSymbol <- inner_join(degDeSeq2,annoData,by = "Gene.id") %>%  # 通过ensgene合并
  dplyr::select(- c(Gene.id)) %>% 
  distinct(Gene.Name,.keep_all = T) %>% 
  column_to_rownames(var = "Gene.Name")

# 设置差异基因筛选条件
pValue <- 0.05

# 上调
pUp <- degSymbol$pvalue < pValue 
upLogFC <- degSymbol$log2FoldChange > 0
up <- pUp & upLogFC


# 下调
downLogFC <- degSymbol$log2FoldChange < 0
down <- pUp & downLogFC

# 合并上调和下调
degUp <- rownames(degSymbol)[up]
degDown <- rownames(degSymbol)[down]
deg <- c(degUp,degDown)


save(degSymbol,degUp,degDown,deg,file = "deg.Rdata")

#GO上调
goMFup <- enrichGO(degUp,
                   keyType = "SYMBOL",    # degEntrezIdl类型
                   OrgDb = "org.Hs.eg.db",  # 人的数据库
                   ont = "MF") 

goBPup <- enrichGO(degUp,
                   keyType = "SYMBOL",
                   OrgDb = "org.Hs.eg.db",
                   ont = "BP")

goCCup <- enrichGO(degUp,
                   keyType = "SYMBOL",
                   OrgDb = "org.Hs.eg.db",
                   ont = "CC") # 细胞组分

goAllup <- enrichGO(degUp,
                    keyType = "SYMBOL",
                    OrgDb = "org.Hs.eg.db",
                    ont = "ALL")

#GO MF上调
goMFup <- as.data.frame(goMFup)
goMFup$negLogP <- -log10(goMFup$pvalue)

goMFup = goMFup |> 
  arrange(desc(negLogP))

write_csv(goMFup, "富集A vs B/GO/file/goMFup.csv")

goMFup10 <- goMFup  %>%
  slice_head(n = 10)

write_csv(goMFup10, "富集A vs B/GO/file/goMFup10.csv")

GOMFUP = ggplot(goMFup10, aes(x = reorder(Description, negLogP), y = negLogP)) +
  geom_col(fill = "#DC143C") +  # 使用统一的颜色
  coord_flip() +  # 翻转坐标轴，使条形图水平显示
  labs(x = "Biological Process", y = "-log10(P value)", title = "Upregulated Molecular Function in GO") +
  theme_minimal() +
  theme(legend.position = "none")

#GO BP上调
goBPup <- as.data.frame(goBPup)
goBPup$negLogP <- -log10(goBPup$pvalue)

goBPup = goBPup |> 
  arrange(desc(negLogP))

write_csv(goBPup, "富集A vs B/GO/file/goBPup.csv")

goBPup10 <- goBPup %>%
  slice_head(n = 10)

write_csv(goBPup10, "富集A vs B/GO/file/goBPup10.csv")

GOBPUP = ggplot(goBPup10, aes(x = reorder(Description, negLogP), y = negLogP)) +
  geom_col(fill = "#DC143C") +  # 使用统一的颜色
  coord_flip() +  # 翻转坐标轴，使条形图水平显示
  labs(x = "Biological Process", y = "-log10(P value)", title = "Upregulated Biological Processes in GO") +
  theme_minimal() +
  theme(legend.position = "none")

#GO CC上调
goCCup <- as.data.frame(goCCup)
goCCup$negLogP <- -log10(goCCup$pvalue)

goCCup = goCCup |> 
  arrange(desc(negLogP))

write_csv(goCCup, "富集A vs B/GO/file/goCCup.csv")

goCCup10 <- goCCup  %>%
  slice_head(n = 10)

write_csv(goCCup10, "富集A vs B/GO/file/goCCup10.csv")

GOCCUP = ggplot(goCCup10, aes(x = reorder(Description, negLogP), y = negLogP)) +
  geom_col(fill = "#DC143C") +  # 使用统一的颜色
  coord_flip() +  # 翻转坐标轴，使条形图水平显示
  labs(x = "Biological Process", y = "-log10(P value)", title = "Upregulated Cellular Component in GO") +
  theme_minimal() +
  theme(legend.position = "none")

#GO ALL上调
goAllup <- as.data.frame(goAllup)
goAllup$negLogP <- -log10(goAllup$pvalue)

goAllup = goAllup |> 
  arrange(desc(negLogP))

write_csv(goAllup, "富集A vs B/GO/file/goAllup.csv")

goAllup10 <- goAllup %>%
  slice_head(n = 10)

write_csv(goAllup10, "富集A vs B/GO/file/goAllup10.csv")


GOALLUP = ggplot(goAllup10, aes(x = reorder(Description, negLogP), y = negLogP)) +
  geom_col(fill = "#DC143C") +  # 使用统一的颜色
  coord_flip() +  # 翻转坐标轴，使条形图水平显示
  labs(x = "Biological Process", y = "-log10(P value)", title = "Up regulation in GO") +
  theme_minimal() +
  theme(legend.position = "none")


#GO下调
goMFdown <- enrichGO(degDown,
                     keyType = "SYMBOL",    # degEntrezIdl类型
                     OrgDb = "org.Hs.eg.db",  # 人的数据库
                     ont = "MF") 

goBPdown <- enrichGO(degDown,
                     keyType = "SYMBOL",
                     OrgDb = "org.Hs.eg.db",
                     ont = "BP")

goCCdown <- enrichGO(degDown,
                     keyType = "SYMBOL",
                     OrgDb = "org.Hs.eg.db",
                     ont = "CC") # 细胞组分

goAlldown <- enrichGO(degDown,
                      keyType = "SYMBOL",
                      OrgDb = "org.Hs.eg.db",
                      ont = "ALL")

#GO MF下调
goMFdown <- as.data.frame(goMFdown)
goMFdown$negLogP <- -log10(goMFdown$pvalue)

goMFdown10 <- goMFdown %>%
  arrange(desc(negLogP)) %>%
  slice_head(n = 10)

GOMFDOWN = ggplot(goMFdown10, aes(x = reorder(Description, negLogP), y = negLogP)) +
  geom_col(fill = "#1E90FF") +  # 使用统一的颜色
  coord_flip() +  # 翻转坐标轴，使条形图水平显示
  labs(x = "Biological Process", y = "-log10(P value)", title = "Downregulated Molecular Function in GO") +
  theme_minimal() +
  theme(legend.position = "none")

#GO BP下调
goBPdown <- as.data.frame(goBPdown)
goBPdown$negLogP <- -log10(goBPdown$pvalue)

goBPdown10 <- goBPdown %>%
  arrange(desc(negLogP)) %>%
  slice_head(n = 10)

GOBPDOWN = ggplot(goBPdown10, aes(x = reorder(Description, negLogP), y = negLogP)) +
  geom_col(fill = "#1E90FF") +  # 使用统一的颜色
  coord_flip() +  # 翻转坐标轴，使条形图水平显示
  labs(x = "Biological Process", y = "-log10(P value)", title = "Downregulated Biological Processes in GO") +
  theme_minimal() +
  theme(legend.position = "none")

#GO CC下调
goCCdown <- as.data.frame(goCCdown)
goCCdown$negLogP <- -log10(goCCdown$pvalue)

goCCdown10 <- goCCdown %>%
  arrange(desc(negLogP)) %>%
  slice_head(n = 10)

GOCCDOWN = ggplot(goCCdown10, aes(x = reorder(Description, negLogP), y = negLogP)) +
  geom_col(fill = "#1E90FF") +  # 使用统一的颜色
  coord_flip() +  # 翻转坐标轴，使条形图水平显示
  labs(x = "Biological Process", y = "-log10(P value)", title = "Downregulated Cellular Component in GO") +
  theme_minimal() +
  theme(legend.position = "none")

#GO ALL下调
goAlldown <- as.data.frame(goAlldown)
goAlldown$negLogP <- -log10(goAlldown$pvalue)

goAlldown10 <- goAlldown %>%
  arrange(desc(negLogP)) %>%
  slice_head(n = 10)

GOALLDOWN = ggplot(goAlldown10, aes(x = reorder(Description, negLogP), y = negLogP)) +
  geom_col(fill = "#1E90FF") +  # 使用统一的颜色
  coord_flip() +  # 翻转坐标轴，使条形图水平显示
  labs(x = "Biological Process", y = "-log10(P value)", title = "Down regulation in GO") +
  theme_minimal() +
  theme(legend.position = "none")

#GO
goMF <- enrichGO(deg,
                 keyType = "SYMBOL",    # degEntrezIdl类型
                 OrgDb = "org.Hs.eg.db",  # 鼠的数据库
                 ont = "MF") 

goBP <- enrichGO(deg,
                 keyType = "SYMBOL",
                 OrgDb = "org.Hs.eg.db",
                 ont = "BP")

goCC <- enrichGO(deg,
                 keyType = "SYMBOL",
                 OrgDb = "org.Hs.eg.db",
                 ont = "CC") # 细胞组分

goAll <- enrichGO(deg,
                  keyType = "SYMBOL",
                  OrgDb = "org.Hs.eg.db",
                  ont = "ALL")

#GO MF
goMF <- as.data.frame(goMF)
goMF$negLogP <- -log10(goMF$pvalue)

goMF10 <- goMF %>%
  arrange(desc(negLogP)) %>%
  slice_head(n = 10)

GOMF = ggplot(goMF10, aes(x = reorder(Description, negLogP), y = negLogP)) +
  geom_col(fill = "#32CD32") +  # 使用统一的颜色
  coord_flip() +  # 翻转坐标轴，使条形图水平显示
  labs(y = "-log10(P value)", title = "MF") +
  theme_minimal() +
  theme(legend.position = "none")

#GO BP
goBP <- as.data.frame(goBP)
goBP$negLogP <- -log10(goBP$pvalue)

goBP10 <- goBP %>%
  arrange(desc(negLogP)) %>%
  slice_head(n = 10)

GOBP = ggplot(goBP10, aes(x = reorder(Description, negLogP), y = negLogP)) +
  geom_col(fill = "#32CD32") +  # 使用统一的颜色
  coord_flip() +  # 翻转坐标轴，使条形图水平显示
  labs( y = "-log10(P value)", title = "BP") +
  theme_minimal() +
  theme(legend.position = "none")

#GO CC
goCC <- as.data.frame(goCC)
goCC$negLogP <- -log10(goCC$pvalue)

goCC10 <- goCC %>%
  arrange(desc(negLogP)) %>%
  slice_head(n = 10)

GOCC = ggplot(goCC10, aes(x = reorder(Description, negLogP), y = negLogP)) +
  geom_col(fill = "#32CD32") +  # 使用统一的颜色
  coord_flip() +  # 翻转坐标轴，使条形图水平显示
  labs( y = "-log10(P value)", title = "CC") +
  theme_minimal() +
  theme(legend.position = "none")

#GO ALL
goAll <- as.data.frame(goAll)
goAll$negLogP <- -log10(goAll$pvalue)

goAll10 <- goAll %>%
  arrange(desc(negLogP)) %>%
  slice_head(n = 10)

GOALL = ggplot(goAll10, aes(x = reorder(Description, negLogP), y = negLogP)) +
  geom_col(fill = "#32CD32") +  # 使用统一的颜色
  coord_flip() +  # 翻转坐标轴，使条形图水平显示
  labs( y = "-log10(P value)", title = "ALL") +
  theme_minimal() +
  theme(legend.position = "none")

plot_layout_GOUP <- GOALLUP/GOMFUP/GOBPUP/GOCCUP
plot_layout_GODOWN <- GOALLDOWN/GOMFDOWN/GOBPDOWN/GOCCDOWN
plot_layout_GOALL = GOALL/GOMF/GOBP/GOCC

ggsave("去掉A4转录分析/富集/GO/GOUPMF.png", GOMFUP,dpi = 500, width = 8, height = 6, bg = "white")
ggsave("去掉A4转录分析/富集/GO/GOUPBP.png", GOBPUP,dpi = 500, width = 8, height = 6, bg = "white")
ggsave("去掉A4转录分析/富集/GO/GOUPCC.png", GOCCUP,dpi = 500, width = 8, height = 6, bg = "white")
ggsave("去掉A4转录分析/富集/GO/GOUPALL.png", GOALLUP,dpi = 500, width = 8, height = 6, bg = "white")

ggsave("去掉A4转录分析/富集/GO/GODOWNMF.png", GOMFDOWN,dpi = 500, width = 10, height = 6, bg = "white")
ggsave("去掉A4转录分析/富集/GO/GODOWNBP.png", GOBPDOWN,dpi = 500, width = 8, height = 6, bg = "white")
ggsave("去掉A4转录分析/富集/GO/GODOWNCC.png", GOCCDOWN,dpi = 500, width = 8, height = 6, bg = "white")
ggsave("去掉A4转录分析/富集/GO/GODOWNALL.png", GOALLDOWN,dpi = 500, width = 8, height = 6, bg = "white")

ggsave("去掉A4转录分析/富集/GO/GOALLMF.png", GOMF,dpi = 500, width = 8, height = 6, bg = "white")
ggsave("去掉A4转录分析/富集/GO/GOALLBP.png", GOBP,dpi = 500, width = 8, height = 6, bg = "white")
ggsave("去掉A4转录分析/富集/GO/GOALLCC.png", GOCC,dpi = 500, width = 8, height = 6, bg = "white")
ggsave("去掉A4转录分析/富集/GO/GOALLALL.png", GOALL,dpi = 500, width = 8, height = 6, bg = "white")



ggsave("去掉A4转录分析/富集/GO/GOUP.png", plot_layout_GOUP,dpi = 500, width = 8, height = 8)
ggsave("去掉A4转录分析/富集/GO/GODOWN.png", plot_layout_GODOWN,dpi = 500, width = 10, height = 8)
ggsave("去掉A4转录分析/富集/GO/GOALL.png", plot_layout_GOALL, dpi = 500, width = 8, height = 8)


#KEGG

#ALL
degId <- mapIds(org.Hs.eg.db,
                keys = deg,  # 需要转化的基因
                column = "ENTREZID", # 转后的类型
                keytype = "SYMBOL")  # 转换前类型

degKEGG <- enrichKEGG(degId,
                      organism = "hsa")

degKEGG@result <- transform(degKEGG@result, negLogPvalue = -log10(pvalue))
top10_degKEGG <- degKEGG@result %>%
  arrange(desc(negLogPvalue)) %>%
  slice_head(n = 10)

KEGGALL = ggplot(top10_degKEGG, aes(x = reorder(Description, negLogPvalue), y = negLogPvalue)) +
  geom_col(fill = "#32CD32") +
  coord_flip() +  # 水平条形图
  labs(x = "KEGG Pathway", 
       y = "-log10(P value)", 
       title = "Top 10 Significant KEGG Pathways") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#UP
degIdup <- mapIds(org.Hs.eg.db,
                  keys = degUp,  # 需要转化的基因
                  column = "ENTREZID", # 转后的类型
                  keytype = "SYMBOL")  # 转换前类型

degKEGGup <- enrichKEGG(degIdup,
                        organism = "hsa")

degKEGGup@result <- transform(degKEGGup@result, negLogPvalue = -log10(pvalue))
top10_degKEGGup <- degKEGGup@result %>%
  arrange(desc(negLogPvalue)) %>%
  slice_head(n = 10)

KEGGUP = ggplot(top10_degKEGGup, aes(x = reorder(Description, negLogPvalue), y = negLogPvalue)) +
  geom_col(fill = "#DC143C") +
  coord_flip() +  # 水平条形图
  labs(x = "KEGG Pathway", 
       y = "-log10(P value)", 
       title = "Up Regulation in KEGG Pathways") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#Down
degIddown <- mapIds(org.Hs.eg.db,
                    keys = degDown,  # 需要转化的基因
                    column = "ENTREZID", # 转后的类型
                    keytype = "SYMBOL")  # 转换前类型

degKEGGdown <- enrichKEGG(degIddown,
                          organism = "hsa")

degKEGGdown@result <- transform(degKEGGdown@result, negLogPvalue = -log10(pvalue))
top10_degKEGGdown <- degKEGGdown@result %>%
  arrange(desc(negLogPvalue)) %>%
  slice_head(n = 10)

KEGGDOWN = ggplot(top10_degKEGGdown, aes(x = reorder(Description, negLogPvalue), y = negLogPvalue)) +
  geom_col(fill = "#1E90FF") +
  coord_flip() +  # 水平条形图
  labs(x = "KEGG Pathway", 
       y = "-log10(P value)", 
       title = "Down Regulation in KEGG Pathways") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("去掉A4转录分析/富集/KEGG/KEGGUP.png", KEGGUP,dpi = 500, width = 8, height = 6, bg = "white")
ggsave("去掉A4转录分析/富集/KEGG/KEGGDOWN.png", KEGGDOWN,dpi = 500, width = 8, height = 6, bg = "white")
ggsave("去掉A4转录分析/富集/KEGG/KEGGALL.png", KEGGALL, dpi = 500, width = 8, height = 6, bg = "white")
plot_layout_KEGG <- KEGGALL/KEGGUP/KEGGDOWN
ggsave("去掉A4转录分析/富集/KEGG/combined_KEGG.png", plot_layout_KEGG, dpi = 500, width = 8, height = 6)

# heatmap ---------------------------------------------------------------

#版本1
library(tidyverse)

rm(list = ls())
setwd("~/Documents/去掉A4转录分析")
gene_name1 = read_csv("gene_tpm_20240704.csv") 
gene_name2 = read_csv("关联图所有基因整理  QLH 去重20240703.csv") |> 
  distinct(`Gene Name`)

data_old3 = gene_name1 |> 
  filter(GeneName %in% gene_name2$`Gene Name`) |> 
  select(-`Gene id`) |> 
  rename(`Gene name` = `GeneName`)

setwd("~/Documents/以往测序的补充数据")



data_old3 <- data_old3 %>%
  left_join(gene_name2, by = "Gene name") %>%
  arrange(match(`Gene name`, gene_name2$`Gene name`)) %>%
  select(names(data_old3))

data_old3 <- as.data.frame(data_old3)

row.names(data_old3) <- data_old3$`Gene name`

# 移除原始的 gene_name 列
data_old3 <- data_old3[,-1]  # 假设 gene_name 是第一列

data_old3 = data_old3 |> 
  mutate("A"       = (`A3`+ `A6`)/2,
         "B"       = (`B3`+ `B4`+ `B6`)/3,
         "C"       = (`C1`+`C2`+`C3`)/3) |> 
  select(A, B, C)

# Create the heatmap
# Load libraries
library(ComplexHeatmap)
library(circlize)

z_score_normalized_data <- t(scale(t(data_old3)))

# Define color mapping
color_mapping <- colorRamp2(c(-1.5, 0, 1.5), c("blue", "white", "red"))

# Define top annotation if needed (commented out here because it seems incorrect)
# top_annotation <- HeatmapAnnotation(df = data.frame(Scale = c(-2, 2)),
#                                     col = list(Scale = c("-2" = "blue", "2" = "red")))

# Create heatmap
heatmap <- Heatmap(z_score_normalized_data,
                   name = "Expression",
                   col = color_mapping,
                   show_row_names = TRUE,
                   show_column_names = TRUE,
                   cluster_rows = TRUE,         
                   cluster_columns = TRUE, 
                   clustering_distance_rows = "euclidean",
                   clustering_distance_columns = "euclidean",
                   clustering_method_rows = "complete",
                   clustering_method_columns = "complete",
                   row_names_side = "right",
                   column_names_side = "bottom",
                   row_dend_side = "left")

png("heatmap20240704-1.png", width = 2609, height = 3548, res = 300)
draw(heatmap, heatmap_legend_side = "right")
dev.off()


#版本2
library(tidyverse)

rm(list = ls())
setwd("~/Documents/去掉A4转录分析")
gene_name1 = read_csv("gene_tpm_20240704.csv") 
gene_name2 = read_csv("关联图所有基因整理  QLH 去重20240704.csv") |> 
  distinct(`Gene Name`)

data_old3 = gene_name1 |> 
  filter(GeneName %in% gene_name2$`Gene Name`) |> 
  select(-`Gene id`) |> 
  rename(`Gene name` = `GeneName`)



data_old3 <- as.data.frame(data_old3)

row.names(data_old3) <- data_old3$`Gene name`

# 移除原始的 gene_name 列
data_old3 <- data_old3[,-1]  # 假设 gene_name 是第一列

data_old3 = data_old3 |> 
  mutate("A"       = (`A3`+ `A6`)/2,
         "B"       = (`B3`+ `B4`+ `B6`)/3,
         "C"       = (`C1`+`C2`+`C3`)/3) |> 
  select(A, B, C)

# Create the heatmap
# Load libraries
library(ComplexHeatmap)
library(circlize)

z_score_normalized_data <- t(scale(t(data_old3)))

# Define color mapping
color_mapping <- colorRamp2(c(-1.5, 0, 1.5), c("blue", "white", "red"))

# Define top annotation if needed (commented out here because it seems incorrect)
# top_annotation <- HeatmapAnnotation(df = data.frame(Scale = c(-2, 2)),
#                                     col = list(Scale = c("-2" = "blue", "2" = "red")))

# Create heatmap
heatmap <- Heatmap(z_score_normalized_data,
                   name = "Expression",
                   col = color_mapping,
                   show_row_names = TRUE,
                   show_column_names = TRUE,
                   cluster_rows = TRUE,         
                   cluster_columns = TRUE, 
                   clustering_distance_rows = "euclidean",
                   clustering_distance_columns = "euclidean",
                   clustering_method_rows = "complete",
                   clustering_method_columns = "complete",
                   row_names_side = "right",
                   column_names_side = "bottom",
                   row_dend_side = "left")

png("heatmap20240704-2.png", width = 2609, height = 3548, res = 300)
draw(heatmap, heatmap_legend_side = "right")
dev.off()



# Venn diagram-----------------------------------------

library(VennDiagram)
#所有差异基因
geneABnewp = geneABnew |> 
  filter(pValue < 0.05)
geneACnewp = geneACnew |> 
  filter(pValue < 0.05)
geneBCnewp = geneBCnew |> 
  filter(pValue < 0.05)


geneABnewp_id = geneABnewp$`Gene id`
geneACnewp_id = geneACnewp$`Gene id`
geneBCnewp_id = geneBCnewp$`Gene id`


ABn12 <- length(intersect(geneABnewp_id, geneACnewp_id))
ABn23 <- length(intersect(geneACnewp_id, geneBCnewp_id))
ABn13 <- length(intersect(geneABnewp_id, geneBCnewp_id))
ABn123 <- length(intersect(intersect(geneABnewp_id, geneACnewp_id), geneBCnewp_id))

#提取余博要求的基因list2275
ABn12_id  = intersect(geneABnewp_id, geneACnewp_id)
ABn123_id = intersect(intersect(geneABnewp_id, geneACnewp_id), geneBCnewp_id)
result_id = setdiff(ABn12_id, ABn123_id)
result_id_name = read_csv("去掉A4转录分析/gene_expression.addAnno.csv") |> 
  select(`Gene id`, `Gene Name`) |> 
  filter(`Gene id` %in% result_id)
write_csv(result_id_name, "去掉A4转录分析/韦恩图/所有差异基因韦恩图余博要求基因列表_20240703_1ed.csv")

#提取余博要求的基因list643

result_id_name_643 = read_csv("去掉A4转录分析/gene_expression.addAnno.csv") |> 
  select(`Gene id`, `Gene Name`) |> 
  filter(`Gene id` %in% ABn123_id)
write_csv(result_id_name_643, "去掉A4转录分析/韦恩图/所有差异基因韦恩图余博要求基因列表643_20240703_1ed.csv")

#提取余博要求的基因list522
ABn23_id  = intersect(geneABnewp_id, geneBCnewp_id)
ABn123_id = intersect(intersect(geneABnewp_id, geneACnewp_id), geneBCnewp_id)
result_id_522 = setdiff(ABn23_id, ABn123_id)
result_id_name_522 = read_csv("去掉A4转录分析/gene_expression.addAnno.csv") |> 
  select(`Gene id`, `Gene Name`) |> 
  filter(`Gene id` %in% result_id_522)
write_csv(result_id_name_522, "去掉A4转录分析/韦恩图/所有差异基因韦恩图余博要求基因列表522_20240703_1ed.csv")


#提取AB，AC共有的，AB，BC共有的，然后取并集
ABn12_id  = intersect(geneABnewp_id, geneACnewp_id)
ABn23_id  = intersect(geneABnewp_id, geneBCnewp_id)
union_id  = union_vector <- union(ABn12_id, ABn23_id)

venn.plot <- draw.triple.venn(
  area1 = length(geneABnewp_id),
  area2 = length(geneACnewp_id),
  area3 = length(geneBCnewp_id),
  n12 = ABn12,
  n23 = ABn23,
  n13 = ABn13,
  n123 = ABn123,
  category = c("Control vs STAT5A-KO", "Control vs UCMSC", "STAT5A-KO vs UCMSC"),
  col = c("#CA0000","#0087FF", "#00BA1D"),
  lwd = 4,
  cat.col = c("#CA0000","#0087FF", "#00BA1D"), 
  cat.cex = 2,
  cex = 2.5,
  fontface = "bold" )

pdf("去掉A4转录分析/韦恩图/显著差异基因.pdf", width=12, height=12)
grid.draw(venn.plot)
dev.off()
grid.draw(venn.plot)



#分化
geneAB12 = read_csv("去掉A4转录分析/新/newD A B.csv")
geneAB34 = read_csv("去掉A4转录分析/新/newD A C.csv")
geneAB56 = read_csv("去掉A4转录分析/新/newD B C.csv")

geneAB12_id = geneAB12$`Gene id`
geneAB34_id = geneAB34$`Gene id`
geneAB56_id = geneAB56$`Gene id`


ABn12 <- length(intersect(geneAB12_id, geneAB34_id))
ABn23 <- length(intersect(geneAB34_id, geneAB56_id))
ABn13 <- length(intersect(geneAB12_id, geneAB56_id))
ABn123 <- length(intersect(intersect(geneAB12_id, geneAB34_id), geneAB56_id))


#提取余博要求的基因list 分化
ABn12_id  = intersect(geneAB12_id, geneAB34_id)
ABn123_id = intersect(intersect(geneAB12_id, geneAB34_id), geneAB56_id)
result_id = setdiff(ABn12_id, ABn123_id)
result_id_name = read_csv("去掉A4转录分析/gene_expression.addAnno.csv") |> 
  select(`Gene id`, `Gene Name`) |> 
  filter(`Gene id` %in% result_id)
write_csv(result_id_name, "去掉A4转录分析/韦恩图/分化差异基因韦恩图余博要求基因列表_20240703_1ed.csv")

venn.plot <- draw.triple.venn(
  area1 = length(geneAB12),
  area2 = length(geneAB34),
  area3 = length(geneAB56),
  n12 = ABn12,
  n23 = ABn23,
  n13 = ABn13,
  n123 = ABn123,
  category = c("Control vs STAT5A-KO", "Control vs UCMSC", "STAT5A-KO vs UCMSC"),
  col = c("#CA0000","#0087FF", "#00BA1D"),
  lwd = 4,
  cat.col = c("#CA0000","#0087FF", "#00BA1D"), 
  cat.cex = 2,
  cex = 2.5,
  fontface = "bold" )

pdf("去掉A4转录分析/韦恩图/分化1.pdf", width=12, height=12)
grid.draw(venn.plot)
dev.off()
grid.draw(venn.plot)

#免疫
geneI12 = read_csv("去掉A4转录分析/新/newimmune A B.csv")
geneI34 = read_csv("去掉A4转录分析/新/newimmune A C.csv")
geneI56 = read_csv("去掉A4转录分析/新/newimmune B C.csv")

geneI12_id = geneI12$`Gene id`
geneI34_id = geneI34$`Gene id`
geneI56_id = geneI56$`Gene id`

In12 <- length(intersect(geneI12, geneI34))
In23 <- length(intersect(geneI34, geneI56))
In13 <- length(intersect(geneI12, geneI56))
In123 <- length(intersect(intersect(geneI12, geneI34), geneI56))

#提取余博要求的基因list 分化
In12_id  = intersect(geneI12_id, geneI34_id)
In123_id = intersect(intersect(geneI12_id, geneI34_id), geneI56_id)
result_id = setdiff(In12_id, In123_id)
result_id_name = read_csv("去掉A4转录分析/gene_expression.addAnno.csv") |> 
  select(`Gene id`, `Gene Name`) |> 
  filter(`Gene id` %in% result_id)
write_csv(result_id_name, "去掉A4转录分析/韦恩图/免疫差异基因韦恩图余博要求基因列表_20240703_1ed.csv")

venn.plot <- draw.triple.venn(
  area1 = length(geneI12),
  area2 = length(geneI34),
  area3 = length(geneI56),
  n12 = In12,
  n23 = In23,
  n13 = In13,
  n123 = In123,
  category = c("Control vs STAT5A-KO", "Control vs UCMSC", "STAT5A-KO vs UCMSC"),
  col = c("#CA0000","#0087FF", "#00BA1D"),
  lwd = 4,
  cat.col = c("#CA0000","#0087FF", "#00BA1D"), 
  cat.cex = 2,
  cex = 2.5,
  fontface = "bold" )

pdf("去掉A4转录分析/韦恩图/免疫1.pdf")
grid.draw(venn.plot)
dev.off()

#年轻
geneO12 = read_csv("去掉A4转录分析/新/newother A B.csv")
geneO34 = read_csv("去掉A4转录分析/新/newother A C.csv")
geneO56 = read_csv("去掉A4转录分析/新/newother B C.csv")


geneO12_id = geneO12$`Gene id`
geneO34_id = geneO34$`Gene id`
geneO56_id = geneO56$`Gene id`

On12 <- length(intersect(geneO12, geneO34))
On23 <- length(intersect(geneO34, geneO56))
On13 <- length(intersect(geneO12, geneO56))
On123 <- length(intersect(intersect(geneO12, geneO34), geneO56))

#提取余博要求的基因list 分化
On12_id  = intersect(geneO12_id, geneO34_id)
On123_id = intersect(intersect(geneO12_id, geneO34_id), geneO56_id)
result_id = setdiff(On12_id, On123_id)
result_id_name = read_csv("去掉A4转录分析/gene_expression.addAnno.csv") |> 
  select(`Gene id`, `Gene Name`) |> 
  filter(`Gene id` %in% result_id)
write_csv(result_id_name, "去掉A4转录分析/韦恩图/年轻差异基因韦恩图余博要求基因列表_20240703_1ed.csv")

venn.plot <- draw.triple.venn(
  area1 = length(geneO12),
  area2 = length(geneO34),
  area3 = length(geneO56),
  n12 = On12,
  n23 = On23,
  n13 = On13,
  n123 = On123,
  category = c("Control vs STAT5A-KO", "Control vs UCMSC", "STAT5A-KO vs UCMSC"),
  col = c("#CA0000","#0087FF", "#00BA1D"),
  lwd = 4,
  cat.col = c("#CA0000","#0087FF", "#00BA1D"), 
  cat.cex = 2,
  cex = 2.5,
  fontface = "bold" ) 

pdf("去掉A4转录分析/韦恩图/年轻1.pdf")
grid.draw(venn.plot)
dev.off()