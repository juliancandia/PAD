library(edgeR)

rm(list=ls())

PROJECT_DIR = "~/git/PAD-main" # replace with your local path

infile = file.path(PROJECT_DIR,"DATA","metadata_samples.csv")
sample_metadata = read.table(infile,header=T,stringsAsFactors=F,sep=",")
sample_label = paste0("S",sample_metadata[,"Code"])
group_label = sample_metadata[,"PAD.non.PAD"]
group_label[group_label=="non-PAD"] = "CTR"
group_sample = paste0(group_label,"_",sample_label)
sample_metadata = data.frame(sample_label,group_label,group_sample,sample_metadata)

infile = file.path(PROJECT_DIR,"DATA","RawCountFile_RSEM_genes_ANNOT.txt")
gene_metadata = as.matrix(read.table(infile,header=T,stringsAsFactors=F,sep="\t",quote=""))

infile = file.path(PROJECT_DIR,"DATA","RawCountFile_RSEM_genes.txt")
data = as.matrix(read.table(infile,header=T,check.names=F,stringsAsFactors=F,sep="\t"))
sample = colnames(data)[-1]
expr = data[,-1]
expr = matrix(as.numeric(expr),ncol=ncol(expr))
rownames(expr) = gene_metadata[,"ENSG.ver"] # unique labels
colnames(expr) = sample
counts = data.frame(expr)

o = match(sample,sample_metadata[,"sample_label"])
sample_metadata = sample_metadata[o,]
cat("Checkpoint passed:",identical(sample_metadata[,"sample_label"],sample),"\n")

group = factor(sample_metadata[,"group_label"])
sex = factor(sample_metadata[,"Sex"])
age = as.numeric(sample_metadata[,"Age"])

d0 <- DGEList(counts)
d0 <- calcNormFactors(d0)

# using function from edgeR
keep = filterByExpr(d0, group=group) # 21922 genes
d = d0[keep,]
pass_filter = rep("N",length(keep))
pass_filter[keep] = "Y"
output = cbind(gene_metadata,pass_filter,data[,-1])
output = rbind(colnames(output),output)
outfile = file.path(PROJECT_DIR,"DATA","RawCountFile_RSEM_genes_wANNOT.txt")
write(t(output),ncol=ncol(output),file=outfile,sep="\t")
gene_metadata = gene_metadata[keep,]

# How do top.table results change based on the model ("group" vs "group+sex+age"?)
#mm <- model.matrix(~0 + group)
mm <- model.matrix(~0 + group + age + sex)
y0 <- voom(d, mm, plot = T, normalize.method="quantile")

fit <- lmFit(y0, mm)
#head(coef(fit))

contr <- makeContrasts(groupPAD - groupCTR, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
# we save the limma object
#outfile = file.path(PROJECT_DIR,"RESULTS","DEG_limma.Robj")
#save(tmp,file=outfile)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
o = match(rownames(top.table),gene_metadata[,"ENSG.ver"])
gene_metadata = gene_metadata[o,]
cat("Checkpoint passed:",identical(gene_metadata[,"ENSG.ver"],rownames(top.table)),"\n")

# we save the top.table results
tmp = cbind(top.table[,c(1,4)],-log10(top.table[,4]),top.table[,5])
tmp = matrix(signif(unlist(tmp),digits=3),ncol=4)
output = rbind(c(colnames(gene_metadata),"logFC","P.Value","-log10(P.Val)","adj.P.Val"),cbind(gene_metadata,tmp))
outfile = file.path(PROJECT_DIR,"DATA","DEG_topTable.txt")
write(t(output),ncol=ncol(output),file=outfile,sep="\t")

# we save the expression matrix (with more annotations)
output = cbind(gene_metadata[match(rownames(y0$E),gene_metadata[,"ENSG.ver"]),],y0$E)
output = rbind(colnames(output),output)
outfile = file.path(PROJECT_DIR,"DATA","Gene_Expression.txt")
write(t(output),ncol=ncol(output),file=outfile,sep="\t")
