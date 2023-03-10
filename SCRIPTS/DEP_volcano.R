library(RColorBrewer)
require(gplots)
library(edgeR)
library(calibrate)
library(readxl)

rm(list=ls())

PROJECT_DIR = "~/git/PAD-main" # replace with your local path

#infile = file.path(PROJECT_DIR,"DATA","PAD_NonPadResults_Sig267Proteins_Annotated_ForLuigi.xlsx")
infile = file.path(PROJECT_DIR,"DATA","Proteins_multivariable_regression.xlsx")
DEP = as.matrix(read_excel(infile,sheet=2))

outfile = file.path(PROJECT_DIR,"RESULTS","FIG2","DEP_Volcano.pdf")
pdf(outfile,width=7,height=5)

gene = DEP[,"Gene"]
beta = as.numeric(DEP[,"PADBeta"])
pval = as.numeric(DEP[,"PADPvalue"])

x = beta
y = -log10(pval)
# we remove 1 outlier
sel = x>1.5
x = x[!sel]
y = y[!sel]
gene = gene[!sel]
pval = pval[!sel]
x_range = c(-1,1)*max(abs(range(x)))
y_range = range(y)
plot(x_range,y_range,type="n",xlab="Effect size",ylab="-log10(p-value)",cex.lab=1.)
### scatterplot symbol parameters
transparency = 60

PAD_up_col = "firebrick3"
PAD_down_col = "navy"

ns_col = "lightgrey"
###
select = pval>=0.05
lines(x[select],y[select],type="p",pch=1,col=ns_col)
select = (x>0)&(pval<0.05)&(pval>=0.01)
col_par = as.numeric(col2rgb(PAD_up_col))
lines(x[select],y[select],type="p",pch=16,col=rgb(col_par[1],col_par[2],col_par[3],
transparency,maxColorValue=255))
select = (x>0)&(pval<0.01)
lines(x[select],y[select],type="p",pch=16,col=PAD_up_col)
select = (x<0)&(pval<0.05)&(pval>=0.01)
col_par = as.numeric(col2rgb(PAD_down_col))
lines(x[select],y[select],type="p",pch=16,col=rgb(col_par[1],col_par[2],col_par[3],
transparency,maxColorValue=255))
select = (x<0)&(pval<0.01)
lines(x[select],y[select],type="p",pch=16,col=PAD_down_col)
#ref_q01 = -log10(max(pval[pval<0.01]))
#abline(h=ref_q01,lty=2)
#ref_q05 = -log10(max(pval[pval<0.05]))
#abline(h=ref_q05,lty=3)
abline(h=-log10(0.01),lty=2)
abline(h=-log10(0.05),lty=3)

# labels to display in volcano plot
gene_name = c("CD276","LDHA",
"RASA1","PTGR2","PPM1B","DNAJB6","HSPB6","HSP90AA1","G3BP2","PARP1","RBBP7","TFAM","TXNRD2","MT-CO1","RPL27","RPF2")

n_gene = length(gene_name)
for (i_gene in 1:n_gene) {
    index = which(gene==gene_name[i_gene])
    textxy(x[index],y[index],gene_name[i_gene],cex = 0.5, offset = 0.65)
}
dev.off()
