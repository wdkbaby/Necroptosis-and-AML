rm(list = ls())
library(limma)
#判断数据集GSE7186是否需要Log2+1
data1<-data.table::fread("GSE7186/GSE7186_series_matrix.csv",data.table = F) ## 读入矩阵文件,warning不用管
rownames(data1) <- data1[,1] ## 第一列存入行名
data1 <- data1[,-1] ## 删除第一列
exprSet <- data1
range(exprSet) ## 查看最大最小值
exprSet <- as.data.frame(exprSet)
ex <- exprSet
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
## 如果LogC=T,就进行循环，如果存在小于0的数，就用非数填入
if (LogC) {
  for(i in 1:ncol(ex)){
    ex[which(ex[,i] <= 0),i] <- NaN
  }
  exprSet <- log2(ex + 1) ## 将ex进行log2+1转化
  print("log2 transform finished")
}else{
  print("log2 transform not needed")
}
data1 <- exprSet
range(data1)

### 标准化----

data1_norm <- normalizeBetweenArrays(data1)
write.csv(data1_norm, "GSE7186/GSE7186_Matrix_norm.csv", row.names = T)

#判断数据集GSE23143是否需要Log2+1
data2<-data.table::fread("GSE23143/GSE23143_series_matrix.csv",data.table = F) ## 读入矩阵文件,warning不用管
rownames(data2) <- data2[,1] ## 第一列存入行名
data2 <- data2[,-1] ## 删除第一列
exprSet <- data2
range(exprSet) ## 查看最大最小值
exprSet = as.data.frame(exprSet)
ex <- exprSet
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
## 如果LogC=T,就进行循环，如果存在小于0的数，就用非数填入
if (LogC) {
  for(i in 1:ncol(ex)){
    ex[which(ex[,i] <= 0),i] <- NaN
  }
  exprSet <- log2(ex + 1) ## 将ex进行log2+1转化
  print("log2 transform finished")
}else{
  print("log2 transform not needed")
}
data2 <- exprSet
range(data2)

### 标准化----

data2_norm <- normalizeBetweenArrays(data2)
write.csv(data2_norm,"GSE23143/GSE23143_Matrix_norm.csv", row.names = T)

# #数据集GSE84334 Log2+1
data3<-data.table::fread("GSE84334/GSE84334_series_matrix.csv",data.table = F) ## 读入矩阵文件,warning不用管
rownames(data3) <- data3[,1] ## 第一列存入行名
data3 <- data3[,-1] ## 删除第一列
exprSet <- data3
range(exprSet) ## 查看最大最小值
exprSet <- as.data.frame(exprSet)
ex <- exprSet
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
# 如果LogC=T,就进行循环，如果存在小于0的数，就用非数填入
if (LogC) {
  for(i in 1:ncol(ex)){
    ex[which(ex[,i] <= 0),i] <- NaN
  }
  exprSet <- log2(ex + 1) ## 将ex进行log2+1转化
  print("log2 transform finished")
} else{
  print("log2 transform not needed")
}
data3 <- exprSet
range(data3)

### 标准化----

data3_norm <- normalizeBetweenArrays(data3)
write.csv(data3_norm, "GSE84334/GSE84334_Matrix_norm.csv", row.names = T)

## GEO数据集合并----
rm(list = ls())
library(magrittr)
library(limma)
group1 <- data.table::fread("GSE7186/GSE7186_group.csv", data.table = F)
group2 <- data.table::fread("GSE23143/GSE23143_group.csv", data.table = F)
group3 <- data.table::fread("GSE84334/GSE84334_group.csv", data.table = F)
data1 <- data.table::fread("GSE7186/GSE7186_Matrix_norm.csv", data.table = F) %>% 
  tibble::column_to_rownames("V1")
data2 <- data.table::fread("GSE23143/GSE23143_Matrix_norm.csv", data.table = F) %>% 
  tibble::column_to_rownames("V1")
data3 <- data.table::fread("GSE84334/GSE84334_Matrix_norm.csv", data.table = F) %>%
  tibble::column_to_rownames("V1")
#data2 <- log2(data2 + 1)

#GSET1分组信息
length(group1$group) 
length(which(group1$group %in% "AML")) 
length(which(group1$group %in% "Normal")) 


#GSET2分组信息
length(group2$group) 
length(which(group2$group %in% "AML")) 
#length(which(group2$group %in% "Control"))


#GSET3分组信息
length(group3$group)
length(which(group3$group %in% "AML"))
#length(which(group3$group %in% "Normal"))
same_gene <- intersect(rownames(data1), rownames(data2)) %>% 
   intersect(rownames(data3))

## 保留数据集中共有基因的数据部分

data1 <- data1[same_gene,]
data2 <- data2[same_gene,]
data3 <- data3[same_gene,]

## 数据集合并

count_matrix <- cbind(data1, data2, data3) %>% 
  na.omit()

### 去除批次效应----

## "1"代表 data1,"2"代表 data2,"3"代表 data3
batch <- c(rep("1", length(group1$group)),
           rep("2", length(group2$group)),
           rep("3", length(group3$group)))

## 去除批次效应
library(sva)
adjusted_counts <- ComBat(count_matrix, batch = batch) 

## 标准化
adjusted_counts <- normalizeBetweenArrays(adjusted_counts)

### 保存合并数据集----

## 分组文件

adjusted_group <- rbind(group1, group2, group3)
#adjusted_group <- rbind(group1, group2)
adjusted_group_sort <- adjusted_group[order(adjusted_group$group, decreasing = T),] ## 注意顺序

length(which(adjusted_group_sort$group %in% "Normal")) 
length(which(adjusted_group_sort$group %in% "AML")) 

## 保存合并数据集的分组信息

write.csv(adjusted_group_sort,"Combined_Datasets_Group.csv", row.names = F)

## 根据合并数据集的分组信息对合并数据集样本进行重新排列

adjusted_counts_sort <- adjusted_counts[,adjusted_group_sort$ID]
write.csv(adjusted_counts_sort, "Combined_Datasets_Matrix.csv")

## 保存疾病样本

adjusted_group_Disease <- adjusted_group_sort[which(adjusted_group_sort$group!="Normal"),]
adjusted_counts_Disease <- adjusted_counts_sort[, adjusted_group_Disease$ID]
write.csv(adjusted_counts_Disease,"Combined_Disease_Matrix.csv")





  ## GEO绘图----
#dir.create("GSE7186")
#dir.create("GSE23143")
#dir.create("GSE84334")
group1 <- data.table::fread("GSE7186/GSE7186_group.csv", data.table = F)
group2 <- data.table::fread("GSE23143/GSE23143_group.csv", data.table = F)
group3 <- data.table::fread("GSE84334/GSE84334_group.csv", data.table = F)
A <- "GSE7186"
B <- "GSE23143"
C <- "GSE84334"

gs <- factor( c(rep(A, length(group1$group)),
                rep(B, length(group2$group)),
                rep(C, length(group3$group))), levels = c(A, B, C))
groups <- make.names(c(A,B,C))
group_col <- c("#1F77B4", "#FF7F0E", "#E64B35")
palette(c("#1F77B4", "#FF7F0E", "#E64B35"))

library(ggplot2)

### theme主题设置----
mytheme <- 
  theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"), # 图像
        plot.title = element_text(size = 7)) + 
  theme(panel.background = element_blank(), # 面板
        panel.grid = element_blank(),
        panel.borde = element_rect(fill = NA, size = 0.75 * 0.47)) + # 添加外框
  theme(axis.line = element_line(size = 0.75 * 0.47), # 坐标轴
        axis.text = element_text(size = 6, color = "black"), # 坐标文字为黑色
        axis.title = element_text(size = 6),
        axis.ticks = element_line(size = 0.75 * 0.47)) +
  theme(legend.key = element_rect(fill = "white"),
        legend.key.size = unit(c(0.3, 0.3), "cm"), # 图注
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.margin = margin(),
        legend.box.margin = margin(),
        legend.box.spacing = unit(0, "cm"),
        legend.background = element_blank(), legend.spacing = unit(0, "cm"),
        legend.box.background = element_blank())

### 矫正前的boxplot图----
dat_before_long <- count_matrix %>% 
  t() %>% 
  data.frame() %>% 
  dplyr::mutate(sample = rownames(.), group = gs) %>% 
  dplyr::mutate(sample = factor(.$sample, levels = .$sample)) %>% 
  tidyr::gather(key = geneid, value, - c(sample, group))

p_before <- 
  ggplot(dat_before_long, aes(sample, value)) + # 映射样本和表达量
  geom_boxplot(aes(fill = group), lwd = 0.1 *0.47, outlier.shape = NA) + # 以分组填充颜色，不要离群点
  # scale_fill_manual(values = c("GSE57345" = "#DA3321", "GSE141910" = "#FC8F57")) +
  scale_fill_manual(values = group_col) +
  labs(title = "Before") + 
  mytheme + # 自己的主题
  theme(legend.position = 'top', # 图注位置（上）
        legend.direction = "horizontal", # 图注方向（水平）
        legend.title = element_blank(), # 隐藏图注标签
        axis.title.x = element_blank(), # 隐藏x轴标签
        axis.title.y = element_blank()) + # 隐藏y轴标签
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1) # 这句注释，下面两句运行，则不显示x轴文本（显示与否取决于样本数量，样本太多不建议展示）
    # axis.text.x = element_blank(), # 隐藏x轴内容
    # axis.ticks.x = element_blank() # 隐藏x轴刻度线
  ) + 
  coord_cartesian(clip = "off") # 解决上右框线看起来浅的问题

p_before ## 看下效果

ggsave("boxplot_before.pdf", plot = p_before, units = "cm",width = 8,height = 5) 

### 矫正后的boxplot图----

dat_after_long <- adjusted_counts %>% 
  t() %>% 
  data.frame() %>% 
  dplyr::mutate(sample = rownames(.), group = gs) %>% 
  dplyr::mutate(sample = factor(.$sample, levels = .$sample)) %>% 
  tidyr::gather(key = geneid, value, - c(sample, group))

p_after <- 
  ggplot(dat_after_long, aes(sample, value)) + # 映射样本和表达量
  geom_boxplot(aes(fill = group), lwd = 0.1 *0.47, outlier.shape = NA) + # 以分组填充颜色，不要离群点
  scale_fill_manual(values = group_col) +
  labs(title = "After") + 
  mytheme + # 自己的主题
  theme(legend.position = 'top', # 图注位置（上）
        legend.direction = "horizontal", # 图注方向（水平）
        legend.title = element_blank(), # 隐藏图注标签
        axis.title.x = element_blank(), # 隐藏x轴标签
        axis.title.y = element_blank()) + # 隐藏y轴标签
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1) # 这句注释，下面两句运行，则不显示x轴文本（显示与否取决于样本数量，样本太多不建议展示）
    # axis.text.x = element_blank(), # 隐藏x轴内容
    # axis.ticks.x = element_blank() # 隐藏x轴刻度线
  ) +   
  coord_cartesian(clip = "off") # 解决上右框线看起来浅的问题

p_after ## 看下效果

ggsave("boxplot_after.pdf", plot = p_after, units = "cm",width = 8,height = 5) 

### PCA矫正前----

count_matrix_pca <- prcomp(t(count_matrix), scale. = T)
count_matrix_pcs <- data.frame(count_matrix_pca$x, group = gs)

## pca1,pca2的百分比
percentage <- round(count_matrix_pca$sdev / sum(count_matrix_pca$sdev) * 100, 2)
percentage <- paste(colnames(count_matrix_pcs),"(", paste(as.character(percentage), "%", ")", sep=""))

pca_before <- ggplot(count_matrix_pcs, aes(x = PC1,y = PC2, color = group)) + 
  geom_point(aes(shape = group)) +
  stat_ellipse(level = 0.95, show.legend = F, geom = "polygon", aes(fill = group),  alpha = 0.2) +
  geom_hline(yintercept = 0, colour = "black", linetype = "longdash", size = 0.75 * 0.47) +
  geom_vline(xintercept = 0, colour = "black", linetype = "longdash", size = 0.75 * 0.47) +
  scale_color_manual(values = group_col) +
  scale_fill_manual(values = group_col) +
  xlab(percentage[1]) +ylab(percentage[2]) +
  labs(title = "Before") +
  mytheme +
  theme(panel.grid = element_line(colour = "grey90", size = 0.75 * 0.47, linetype = 1), 
        plot.title = element_text(size = 7, hjust = 0, vjust = 0.5, margin = unit(c(0.1,0.1,0.1,0.1), "cm")),
        legend.position = 'top', # 图注位置（上）
        legend.direction = "horizontal",
        legend.margin = margin(0, 0, 0, 0, "cm"),
        legend.title  = element_blank()) + # 图注方向（水平）
  coord_cartesian(clip = "off")

pca_before

ggsave("pca_before.pdf", plot = pca_before, units = "cm",width = 6, height = 6) 

### PCA矫正后----

adjust_matrix_pca <- prcomp(t(adjusted_counts), scale. = T)
adjust_matrix_pcs <- data.frame(adjust_matrix_pca$x, group = gs)
percentage <- round(adjust_matrix_pca$sdev / sum(adjust_matrix_pca$sdev) * 100, 2)
percentage <- paste(colnames(adjust_matrix_pcs),"(", paste(as.character(percentage), "%", ")", sep=""))

pca_after <- ggplot(adjust_matrix_pcs, aes(x = PC1,y = PC2, color = group)) + 
  geom_point(aes(shape = group)) +
  stat_ellipse(level = 0.95, show.legend = F, geom = "polygon", aes(fill = group),  alpha = 0.2) +
  geom_hline(yintercept = 0, colour = "black", linetype = "longdash", size = 0.75 * 0.47) +
  geom_vline(xintercept = 0, colour = "black", linetype = "longdash", size = 0.75 * 0.47) +
  scale_color_manual(values = group_col) +
  scale_fill_manual(values = group_col) +
  xlab(percentage[1]) +ylab(percentage[2]) +
  labs(title = "After") +
  mytheme +
  theme(panel.grid = element_line(colour = "grey90", size = 0.75 * 0.47, linetype = 1), 
        plot.title = element_text(size = 7, hjust = 0, vjust = 0.5, margin = unit(c(0.1,0.1,0.1,0.1), "cm")),
        legend.position = 'top', # 图注位置（上）
        legend.direction = "horizontal",
        legend.margin = margin(0, 0, 0, 0, "cm"),
        legend.title  = element_blank()) + # 图注方向（水平）
  coord_cartesian(clip = "off")

pca_after

ggsave("pca_after.pdf", plot = pca_after, units = "cm",width = 6, height = 6) 

save.image() ## 保存


##差异表达分析
rm(list=ls()) 
library("limma")
library('dplyr')

# 准备工作 --------------------------------------------------------------------

expname <- "Combined_Datasets_Matrix.csv" ## 输入文件名
groupname <- "Combined_Datasets_Group.csv" ## 分组文件名
outputname <- "Combined_Datasets_DEGs.csv" ## 输出文件名
groupnameA <- "Normal" ## 正常组名
groupnameB <- "AML" ## 疾病组名

# 数据处理 --------------------------------------------------------------------

GSE_mat <- data.table::fread(expname,data.table=F) %>% 
  tibble::column_to_rownames("V1")## 读入

sampleinfo <- data.table::fread(groupname,data.table=F) ## 读入
sampleinfo<-filter(sampleinfo, sampleinfo$group != "")
sampleinfo <- sampleinfo[order(sampleinfo$group),]
GSE_mat <- GSE_mat[,sampleinfo$ID]

raw_count <- GSE_mat
df <- sampleinfo
group <- sampleinfo[,2]
group <- factor(group) ## 因子化，并强制排序
head(group)
# 实验设计矩阵 ------------------------------------------------------------------

design <- model.matrix(~ 0 + group)
rownames(design) <- df$ID ## 行名用df的ID填入
colnames(design) <- levels(group) ## 列名用group的group填入
head(design) ## 看design的前六行

# 差异分析 --------------------------------------------------------------------

## 线性建模
fit <- lmFit(raw_count,design) ##raw_count,design样本顺序要对应
cont.matrix <- makeContrasts(contrasts = c("AML-Normal"), ## 需疾病组在前
                             levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)
## 经验贝叶斯调整
fit2 <- eBayes(fit2) 
plotSA(fit2)
## 筛选差异基因
dif <- topTable(fit2, coef = 1, n = Inf)
dif <- na.omit(dif)
dif$gene_symbol <- row.names(dif)
head(dif)
dif.up <- dif %>%
  filter(logFC > 0 & adj.P.Val < 0.05)
dif.down <- dif %>%
  filter(logFC < (-0) & adj.P.Val < 0.05)
dif <- dif %>% mutate(group = case_when(
  gene_symbol %in% dif.up$gene_symbol ~ "up",
  gene_symbol %in% dif.down$gene_symbol ~ "down",
  TRUE ~ "no"))
DEG <- dif
write.csv(DEG, outputname)

