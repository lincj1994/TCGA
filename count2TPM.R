#Count to FPKM
##1.获得每个基因所有外显子长度之和
txdb <- makeTxDbFromGFF("Homo_sapiens.GRCh38.99.gtf",format="gtf")
exons_gene <- exonsBy(txdb, by = "gene")
exons_gene_lens <- lapply(exons_gene,function(x){sum(width(reduce(x)))}) ###为一个list，记录每个基因所有外显子长度之和
exons_gene_lens2 <- data.frame(do.call(rbind,exons_gene_lens2)) ###将list转换为data.frame

##2.将原表达矩阵中不存在于exons_gene_lens的基因剔除
###注意，表达矩阵行为基因，列为样本(56499 x 1203)
row.num <- which(rownames(exons_gene_lens) %in% rownames(exprSet3))
exons_gene_lens2 <- exons_gene_lens[row.num,]
exprSet4 <- exprSet3[names(exons_gene_lens2),]
all(names(exons_gene_lens2) == rownames(exprSet4))

##3.count to FPKM
mapRead <- colSums2(exprSet4) ###某个样本的所有基因reads总和
exprFPKM <- lapply(1:1203, function(i){
  exprFPKM <- exprSet4[,i]/(exons_gene_lens2[,1] * mapRead[i])  ###表达矩阵每一列÷exons矩阵中对应基因外显子长度÷该列之和
  return(exprFPKM)
})
###注意，针对上条循环内的函数，当相同长度的两列数据相除时，默认同一行的两数据相除,如1:3÷1:3 = 1,1,1
exprFPKM2 <- data.frame(do.call(cbind,exprFPKM))
exprFPKM2 <- exprFPKM2 * 10^9 ###为FPKM值

##4.FPKM to TPM
s <- colSums2(as.matrix(exprFPKM2)) ###计算每一个样本所有FPKM之和
exprTPM <- lapply(1:1203, function(i){
  exprTPM <- exprFPKM2[,i] * 10^6/s[i]  ###换算
  return(exprTPM)
})
exprTPM2 <- data.frame(do.call(cbind,exprTPM))
rownames(exprTPM2) <- rownames(exprSet4)
colnames(exprTPM2) <- colnames(exprSet4)
