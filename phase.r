library(dplyr)

# 设置路径和读取数据
hap1 <- read.table("/home/Mzhou/02.F2/xh/merged_Hap1_counts.txt", sep=" ", header=T, check.names = FALSE)
hap2 <- read.table("/home/Mzhou/02.F2/xh/merged_Hap2_counts.txt", sep=" ", header=T, check.names = FALSE)

# 提取基因名
genenames <- as.character(hap1[,1])

# 提取表达量矩阵并添加小数值避免除以零
hap1 <- hap1 %>% select(-1) %>% mutate_all(as.numeric) %>% as.matrix() + 0.001
hap2 <- hap2 %>% select(-1) %>% mutate_all(as.numeric) %>% as.matrix() + 0.001

# 计算总表达量和过滤低表达基因
total_expr <- hap1 + hap2
idx <- total_expr < 5

# 计算父本比例 (假设Hap1为父本)
paternal_ratio <- hap1 / total_expr
paternal_ratio[idx] <- NA

# 过滤缺失值过多的基因
kept <- apply(paternal_ratio, 1, function(x) { sum(is.na(x)) }) < 0.9 * dim(paternal_ratio)[2]
paternal_ratio <- paternal_ratio[kept, ]
genenames <- genenames[kept]

# 对每个基因进行t检验，检测表达偏倚 (检验均值是否偏离0.5)
mytest <- apply(paternal_ratio, 1, function(x) {
  x <- na.omit(x)
  if(length(x) < 3) {
    return(c(NA, NA, NA, NA))  # 样本量太小无法进行检验
  }
  if (sd(x) < 1e-10) {  # 设置一个极小的阈值
    return(c(mean(x), NA, NA, length(x)))  # 数据几乎恒定，不进行t检验
  }
tmp <- t.test(x, mu=0.5)
  return(c(tmp$estimate, tmp$statistic, tmp$p.value, length(x)))
})

# 整理结果矩阵
mytest <- t(mytest)
colnames(mytest) <- c("Mean_Paternal_Ratio", "T_statistic", "P_value", "Sample_Size")
mytest <- cbind(GeneName = genenames, mytest)

# 读取biased.genes文件并处理父本比例
dge <- read.table('/home/Mzhou/02.F2/xh/biased.genes', header=T)
#dge$Paternal_Ratio <- dge$father_counts / (dge$father_counts + dge$mother_counts + 0.001)  # 直接计算父本比例

# 按基因名排序
dge <- dge[order(dge$Genename), ]

# 自定义t检验函数 (检验均值是否偏离0.5)

mu_test <- function(x) {
  x <- na.omit(x)
  res <- matrix(NA, nrow=1, ncol=4)
  if(length(x) < 3) {
    res[1, ] <- c(length(x), NA, NA, NA)
    return(res)
  } 
  if (sd(x) < 1e-10) {  # 设置一个极小的阈值
    res[1, ] <- c(length(x), mean(x), NA, NA)  # 数据几乎恒定，不进行t检验
    return(res)
  } 
    tmp <- t.test(x, mu=0.5)
    res[1, ] <- c(length(x), tmp$estimate, tmp$statistic, tmp$p.value)
    return(res)  
}


# 按基因分组进行t检验
dge_test <- tapply(dge$ratio, dge$Genename, mu_test)
un_dge_test <- unlist(dge_test)

# 整理结果矩阵
dge_testres <- matrix(NA, nrow=length(dge_test), ncol=4)
dge_testres[, 1] <- un_dge_test[seq(1, length(un_dge_test), 4)]
dge_testres[, 2] <- un_dge_test[seq(2, length(un_dge_test), 4)]
dge_testres[, 3] <- un_dge_test[seq(3, length(un_dge_test), 4)]
dge_testres[, 4] <- un_dge_test[seq(4, length(un_dge_test), 4)]

# 设置行列名
rownames(dge_testres) <- names(dge_test)
colnames(dge_testres) <- c("Nind", "Mean_Paternal_Ratio", "Tvalue", "Pvalue")

# 输出结果
write.table(mytest, file = "/home/Mzhou/02.F2/xh/AllelicImprinting_test.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(dge_testres, file = "/home/Mzhou/02.F2/xh/Imprinting_test.txt", quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
