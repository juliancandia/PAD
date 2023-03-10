library(readxl)
require(RColorBrewer)

rm(list=ls())

PROJECT_DIR = "~/git/PAD-main" # replace with your local path

pathway_ref = c("BIOCARTA","HALLMARK","KEGG","PID","REACTOME","WP")
pathway_ref_ini = c("BC","HM","KG","PID","RE","WP")

infile = file.path(PROJECT_DIR,"RESULTS","FIG2","Supplementary_Data_4.xlsx")
data = as.matrix(read_excel(infile,sheet=1))
data = data[which(data[,"In plot?"]=="Y"),c(1,3,8)]

tmp = table(data[,2])
tmp = tmp[order(-tmp)]
group = names(tmp)
n_group = length(group)

res = NULL
for (i_group in 1:n_group) {
    tmp = data[data[,2]==group[i_group],,drop=F]
    ES = as.numeric(tmp[,3])
    index = which(ES>0)
    index = index[order(-ES[index])]
    index2 = which(ES<0)
    index2 = index2[order(ES[index2])]
    tmp = tmp[c(index,index2),]
    res = rbind(res,tmp)
}
# add ref annotations and improve labels
pwy = res[,"Pathway"]
n_pwy = length(pwy)
ref_annot = rep(NA,n_pwy)
pwy_label = rep(NA,n_pwy)
for (i_pwy in 1:n_pwy) {
    tmp = strsplit(pwy[i_pwy],"_")[[1]]
    ref_annot[i_pwy] = tmp[1]
    pwy_label[i_pwy] = paste0(tolower(paste(tmp[-1],collapse=" "))," (",pathway_ref_ini[which(pathway_ref==tmp[1])],")")
}
# we save and *manually* rework the pathway labels
outfile = file.path(PROJECT_DIR,"RESULTS","FIG2","Pathways_RAW.txt")
output = cbind(1:n_pwy,pwy,ref_annot,pwy_label,res[,2:3])
write(t(output),ncol=ncol(output),file=outfile,sep="\t")
infile = file.path(PROJECT_DIR,"RESULTS","FIG2","Pathways_EDITED.txt")
pwy = as.matrix(read.table(infile,header=F,sep="\t",comment.char="",quote=""))
# we add suffixes
pwy_label = rep(NA,n_pwy)
for (i_pwy in 1:n_pwy) {
    pwy_label[i_pwy] = paste0(pwy[i_pwy,6]," (",pathway_ref_ini[which(pathway_ref==pwy[i_pwy,2])],")")
}
category = pwy[,4]
significance = as.numeric(pwy[,5])
names(significance) = pwy_label

color_scheme_from_Luigi = c("#feb8f9","#9abbf3","#ffffa2","#c2b2d6","#9dd79d","#fdc897","#b6e4eb","#9ca07f")

outfile = file.path(PROJECT_DIR,"RESULTS","FIG2","GSEA_barplot.pdf")
pdf(outfile,width=8,height=4)

cat_label = unique(rev(category))
cat_col = rep(color_scheme_from_Luigi,10)
col = rep(NA,n_pwy)
for (i_pwy in 1:n_pwy) {
   col[i_pwy] = cat_col[which(cat_label==category[i_pwy])]
}

par(mar=c(3.1, 11.1, 3.1, 1.1)) #par(mar = c(bottom, left, top, right))

barplot(rev(significance),main="",xlab="Enrichment score",horiz=T,las=1,cex.names=0.35,col=rev(col),xlim=c(-0.8,0.8),axes=F)
axis(side = 1, at = seq(-0.8,0.8,by=0.2), labels = seq(-0.8,0.8,by=0.2))
dev.off()
