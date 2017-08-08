# module load R/3.2.5

### load libraries
library(DESeq2)
library(dplyr)
library(annotables)
library(calibrate)

### read in data:
counts = read.table("2017.02.21.BiPSC.gene.counts.tsv",h=T,row.names=1)
gene_ann = counts[,1,drop=F]
counts = counts[,-1]

tpms = read.table("2017.02.21.BiPSC.gene.tpm.tsv",h=T,row.names=1)
gene_ann_tpm = tpms[,1,drop=F]
tpms = tpms[,-1]

### filter the counts and TPMs to only include protein-coding genes

prot_coding = grch37$ensgene[grch37$biotype == "protein_coding"]
ensrownames=sapply(strsplit(rownames(counts),split="\\."), function(x) x[1])
coding_counts = counts[which(ensrownames %in% prot_coding),]
coding_gene_ann = gene_ann[rownames(coding_counts),,drop=F]

ensrownames_tpms=sapply(strsplit(rownames(tpms),split="\\."), function(x) x[1])
coding_tpms = tpms[which(ensrownames_tpms %in% prot_coding),]
coding_tpm_ann = gene_ann_tpm[rownames(coding_tpms),,drop=F]


###################### Differential expression analysis ###############################

stage = substr(sapply(strsplit(colnames(counts), split ="\\."),function(x)x[2]),0,2)
ipsc = substr(sapply(strsplit(colnames(counts), split ="\\."),function(x)x[2]),3,7)
subject = substr(sapply(strsplit(colnames(counts), split ="\\."),function(x)x[3]),1,1)
line = sapply(strsplit(colnames(counts), split ="\\."),function(x)x[3])
rep = rep(1,32)
rep[which(duplicated(paste(sapply(strsplit(colnames(counts), split ="\\."),function(x)x[2]), sapply(strsplit(colnames(counts), split ="\\."),function(x)x[3]),sep="_")))] <- 2

design = data.frame(row.names = colnames(counts), stage = stage, ipsc = ipsc, subject = subject, line = line, rep= rep)

#################### FiPSC vs BiPSC at DE ###########################
design_DE = design[which(design$stage == "DE"),]

dds_DE <- DESeqDataSetFromMatrix(
  countData = coding_counts[,which(design$stage == "DE")],
  colData = design_DE,
  design = ~ipsc)

dds_DE$ipsc <- relevel( dds_DE$ipsc, "FIPSC" )

dds_DE <- DESeq(dds_DE)
res_DE = results(dds_DE)
length(which(res_DE$padj < 0.05))	
### 1247 DEGs at 5% FDR

length(which(res_DE$padj < 0.05 & res_DE$log2FoldChange>0))	
### 567 DEGs at 5% FDR, up-regulated in BiPSC DE

length(which(res_DE$padj < 0.05 & res_DE$log2FoldChange<0))
### 680 DEGs at 5% FDR, up-regulated in BiPSC DE

res_DE = cbind(coding_gene_ann, res_DE)
write.table(res_DE, file = "DE.BiPSC_vs_FiPSC.DEGs.DESeq2.tsv",sep="\t",quote=F,col.names=NA)

#################### FiPSC vs BiPSC at PP ###########################

design_PP = design[which(design$stage == "PP"),]

dds_PP <- DESeqDataSetFromMatrix(
  countData = coding_counts[,which(design$stage == "PP")],
  colData = design_PP,
  design = ~ipsc)

dds_PP$ipsc <- relevel( dds_PP$ipsc, "FIPSC" )

dds_PP <- DESeq(dds_PP)
res_PP = results(dds_PP)
length(which(res_PP$padj < 0.05))	### 607
### 607 DEGs at 5% FDR

length(which(res_PP$padj < 0.05 & res_PP$log2FoldChange>0))	
### 181 DEGs at 5% FDR, up-regulated in BiPSC DE

length(which(res_PP$padj < 0.05 & res_PP$log2FoldChange<0))
### 426 DEGs at 5% FDR, up-regulated in BiPSC DE

res_PP = cbind(coding_gene_ann, res_PP)
write.table(res_PP, file = "PP.BiPSC_vs_FiPSC.coding.DEGs.DESeq2.tsv",sep="\t",quote=F,col.names=NA)

################## Volcano plots highlighting top DEGs #####################
############################## Fig. 3E and 3F ##############################
pdf("DE.volcano_plot.pdf")
col=rep(rgb(0,0,0,50,maxColorValue=255),nrow(res_DE))
col[which(res_DE$padj<0.05)] <- rgb(0,0,100,50,maxColorValue=255)
with(res_DE, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot - DE", xlim=c(-3,3), ylim=c(0,25), col=col))

# Label points with the textxy function from the calibrate plot
with(subset(res_DE, padj<1e-10 & abs(log2FoldChange)>0.5), textxy(log2FoldChange, -log10(pvalue), labs=GeneName, cex=.7, offset=0.5 ))
dev.off()

pdf("PP.volcano_plot.pdf")
col=rep(rgb(0,0,0,50,maxColorValue=255),nrow(res_PP))
col[which(res_PP$padj<0.05)] <- rgb(0,0,100,50,maxColorValue=255)
with(res_PP, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot - PP", xlim=c(-3,6), col=col))
# Label points with the textxy function from the calibrate plot
with(subset(res_PP, padj<1e-10 & abs(log2FoldChange)>0.5), textxy(log2FoldChange, -log10(pvalue), labs=GeneName, cex=.7, offset=0.5 ))

### make the same volcano plot again, but cut-off the Y axis at 25, and add ESR1 manually to the plot in the top right corner (padj=1.12e-101)
### next most significant gene has a padj=2.63e-18
with(res_PP, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot - PP", xlim=c(-3,2), ylim=c(0,25), col=col))

# Label points with the textxy function from the calibrate plot
with(subset(res_PP, padj<1e-10 & abs(log2FoldChange)>0.5), textxy(log2FoldChange, -log10(pvalue), labs=GeneName, cex=.7, offset=0.5 ))

points(2,25,pch=20, col="blue")
text(1.75,25, "ESR1*",cex=.7)
dev.off()


########### Hypergeometric enrichments of Bi-DOCS and Fi-DOCS genes ###############
########## and FOXA2 targets at DE (Wang et al. 2015)  ##########
########## (the genes were obtained using Bi-DOCS and Fi-DOCS boundaries ##########
########## plugged into GREAT using the default settings #########################

b_genes = read.table("ATAC_comparison/BiPSC.DOCS.ann.genes.txt",sep="\t")
b_docs = b_genes$V1
f_genes = read.table("ATAC_comparison/FiPSC.DOCS.ann.genes.txt",sep="\t")
f_docs = f_genes$V1

foxa2_DE_targets = read.table("Other_datasets/FOXA2.dev_stages/FOXA2.DE.genes.txt",sep="\t")
foxa2_DE=foxa2_DE_targets$V1

## intersect the FOXA2 targets with the genes expressed at this stage
## define expressed as TPM>=1 in all samples at the given stage

de_exp_tpm=which(rowSums(coding_tpms[,grep("DE",colnames(coding_tpms))]>=1)==16)
de_exp_genes = coding_tpm_ann[de_exp_tpm,]
foxa2_DE_expressed=intersect(foxa2_DE_targets$V1, de_exp_genes)

## define up- and down-regulated DEGs
de_degs_up = unique(res_DE[which(res_DE$padj<0.05 & res_DE$log2FoldChange>0),]$GeneName)
pp_degs_up = unique(res_PP[which(res_PP$padj<0.05 & res_PP$log2FoldChange>0),]$GeneName)
de_degs_down = unique(res_DE[which(res_DE$padj<0.05 & res_DE$log2FoldChange<0),]$GeneName)
pp_degs_down = unique(res_PP[which(res_PP$padj<0.05 & res_PP$log2FoldChange<0),]$GeneName)

# run hypergeometric enrichment

gene_lists = list(DE_degs_up=de_degs_up, DE_degs_down=de_degs_down, PP_degs_up=pp_degs_up, PP_degs_down=pp_degs_down)
targets=list(B_docs=b_docs, F_docs=f_docs, FOXA2_DE=foxa2_DE_expressed)

enrichment = matrix(, nrow = length(gene_lists), ncol = length(targets))
rownames(enrichment) = names(gene_lists)
colnames(enrichment) = names(targets)

for (i in 1:length(gene_lists)){
	for (j in 1:length(targets)){
		A=length(intersect(unlist(gene_lists[i]), unlist(targets[j]) )) ## overlap
		B=length(unlist(gene_lists[i]))#
		C=length(unique(res_DE$GeneName))
		D=length(which(unlist(targets[j]) %in% unique(res_DE$GeneName)))
		enrichment[i,j] = phyper(A-1,B,C-B,D,lower.tail=F)
	}
}

### correct for multiple testing:
enrichment_adj=matrix(p.adjust(enrichment, method="BH"), nrow=4, ncol=3)
rownames(enrichment_adj)=rownames(enrichment)
colnames(enrichment_adj)=colnames(enrichment)

#                    B_docs       F_docs     FOXA2_DE
# DE_degs_up   3.387110e-07 6.256955e-01 1.419588e-53
# DE_degs_down 9.074325e-02 2.536221e-14 9.999881e-01
# PP_degs_up   9.125416e-04 9.642194e-01 1.560222e-03
# PP_degs_down 1.197380e-02 4.980827e-08 9.999881e-01

#### check for how many Bi-DEGs genes we detect a Bi-DOCS site, and Fi-DEGs a Fi-DOCS site, at DE:
## BiPSC:
length(intersect(unlist(gene_lists$DE_degs_up), unlist(targets$B_docs) ))
## 126
## FiPSC:
length(intersect(unlist(gene_lists$DE_degs_down), unlist(targets$F_docs) ))
## 145


#### check logFC and p-value for CXCR4
res_DE[which(res_DE$GeneName=="CXCR4"),]

save.image("bipsc_fipsc.atac_paper.Rdata")
