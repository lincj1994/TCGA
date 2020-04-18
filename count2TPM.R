#Count to FPKM
##1.���ÿ���������������ӳ���֮��
txdb <- makeTxDbFromGFF("Homo_sapiens.GRCh38.99.gtf",format="gtf")
exons_gene <- exonsBy(txdb, by = "gene")
exons_gene_lens <- lapply(exons_gene,function(x){sum(width(reduce(x)))}) ###Ϊһ��list����¼ÿ���������������ӳ���֮��
exons_gene_lens2 <- data.frame(do.call(rbind,exons_gene_lens2)) ###��listת��Ϊdata.frame

##2.��ԭ��������в�������exons_gene_lens�Ļ����޳�
###ע�⣬���������Ϊ������Ϊ����(56499 x 1203)
row.num <- which(rownames(exons_gene_lens) %in% rownames(exprSet3))
exons_gene_lens2 <- exons_gene_lens[row.num,]
exprSet4 <- exprSet3[names(exons_gene_lens2),]
all(names(exons_gene_lens2) == rownames(exprSet4))

##3.count to FPKM
mapRead <- colSums2(exprSet4) ###ĳ�����������л���reads�ܺ�
exprFPKM <- lapply(1:1203, function(i){
  exprFPKM <- exprSet4[,i]/(exons_gene_lens2[,1] * mapRead[i])  ###�������ÿһ�С�exons�����ж�Ӧ���������ӳ��ȡ¸���֮��
  return(exprFPKM)
})
###ע�⣬�������ѭ���ڵĺ���������ͬ���ȵ������������ʱ��Ĭ��ͬһ�е����������,��1:3��1:3 = 1,1,1
exprFPKM2 <- data.frame(do.call(cbind,exprFPKM))
exprFPKM2 <- exprFPKM2 * 10^9 ###ΪFPKMֵ

##4.FPKM to TPM
s <- colSums2(as.matrix(exprFPKM2)) ###����ÿһ����������FPKM֮��
exprTPM <- lapply(1:1203, function(i){
  exprTPM <- exprFPKM2[,i] * 10^6/s[i]  ###����
  return(exprTPM)
})
exprTPM2 <- data.frame(do.call(cbind,exprTPM))
rownames(exprTPM2) <- rownames(exprSet4)
colnames(exprTPM2) <- colnames(exprSet4)