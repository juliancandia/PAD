library(RColorBrewer)
require(gplots)
library(edgeR)
library(calibrate)

rm(list=ls())

PROJECT_DIR = "~/git/PAD-main" # replace with your local path

infile = file.path(PROJECT_DIR,"DATA","DEG_topTable.txt")
DEG = read.table(infile,header=T,stringsAsFactors=F,quote="",sep="\t")

outfile = file.path(PROJECT_DIR,"RESULTS","FIG1","DEG_Volcano_NonCoding.pdf")
pdf(outfile,width=5,height=5)
# we select only non-protein coding genes
index = c(which(DEG[,"type"]!="protein_coding"),which(is.na(DEG[,"type"])))
DE = DEG[index,]
x = DE[,"logFC"]
y = -log10(DE[,"P.Value"])
x_range = range(x)
x_range[2] = 4.17
y_range = range(y)
plot(x_range,y_range,type="n",xlab="log2(PAD/CTR)",ylab="-log10(p-val)",cex.lab=1.)
### scatterplot symbol parameters
transparency = 60
PAD_up_col = "firebrick3"
PAD_down_col = "navy"
ns_col = "lightgrey"
select = DE[,"adj.P.Val"]>=0.05
lines(x[select],y[select],type="p",pch=1,col=ns_col)
select = (x>0)&(DE[,"adj.P.Val"]<0.05)&(DE[,"adj.P.Val"]>=0.01)
col_par = as.numeric(col2rgb(PAD_up_col))
lines(x[select],y[select],type="p",pch=16,col=rgb(col_par[1],col_par[2],col_par[3],
transparency,maxColorValue=255))
select = (x>0)&(DE[,"adj.P.Val"]<0.01)
lines(x[select],y[select],type="p",pch=16,col=PAD_up_col)
select = (x<0)&(DE[,"adj.P.Val"]<0.05)&(DE[,"adj.P.Val"]>=0.01)
col_par = as.numeric(col2rgb(PAD_down_col))
lines(x[select],y[select],type="p",pch=16,col=rgb(col_par[1],col_par[2],col_par[3],
transparency,maxColorValue=255))
select = (x<0)&(DE[,"adj.P.Val"]<0.01)
lines(x[select],y[select],type="p",pch=16,col=PAD_down_col)
ref_q01 = -log10(max(DE[DE[,"adj.P.Val"]<0.01,"P.Value"]))
abline(h=ref_q01,lty=2)
ref_q05 = -log10(max(DE[DE[,"adj.P.Val"]<0.05,"P.Value"]))
abline(h=ref_q05,lty=3)

select = DE[,"name"]=="MIR210HG"
textxy(x[select],y[select],DE[select,"name"],cex = 0.55, offset = 0.65)

axis(4,at=c(ref_q01,ref_q05),labels=paste0("q=",c("0.01","0.05")),cex.axis=0.7)
dev.off()
