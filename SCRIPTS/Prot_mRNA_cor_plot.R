library(RColorBrewer)
library(readxl)
library(calibrate)

rm(list=ls())

PROJECT_DIR = "~/git/PAD-main" # replace with your local path

CL = 95

# gene expression vs protein abundance correlations
infile = file.path(PROJECT_DIR,"RESULTS","FIG5",paste0("Prot_mRNA_Cor_CL",CL,"_pval_groups.txt"))
data = as.matrix(read.table(infile,stringsAsFactors=F,header=T,sep="\t",quote="",comment.char=""))
r_CTR = as.numeric(data[,"r_CTR"])
r_PAD = as.numeric(data[,"r_PAD"])
diff = as.numeric(data[,"diff"])
r_CTR_pval = as.numeric(data[,"r_CTR_pval"])
r_PAD_pval = as.numeric(data[,"r_PAD_pval"])

file = c("CTR_mitochondrion","CTR_translation","PAD_extracellular_exosome","PAD_tubulin")
title = c("Mitochondrion (GO-CC)","tRNA aminoacylation (Reactome)","Extracellular exosome (GO-CC)","Gamma-tubulin complex (GO-CC)")
enriched = c("CTR","CTR","PAD","PAD")
col = c("red","cyan2","blue","green2")
n_file = length(file)
sel_target = vector("list",n_file)
for (i_file in 1:n_file) {
    infile = file.path(PROJECT_DIR,"RESULTS","FIG5",paste0(file[i_file],".xlsx"))
    target_gene = as.matrix(read_excel(infile,sheet=1,col_names=F))[,2]
    if (enriched[i_file]=="CTR") {
        sel = data[,"group"]%in%c("CTR_ns_PAD_neg","CTR_pos_PAD_neg","CTR_pos_PAD_ns")
    } else if (enriched[i_file]=="PAD") {
        sel = data[,"group"]%in%c("PAD_ns_CTR_neg","PAD_pos_CTR_neg","PAD_pos_CTR_ns")
    }
    sel_target[[i_file]] = sel&(data[,"ENSG"]%in%target_gene)
}

sel_target_all = sel_target[[1]]
for (i_file in 2:n_file) {
    sel_target_all = sel_target_all|sel_target[[i_file]]
}

outfile = file.path(PROJECT_DIR,"RESULTS","FIG5",paste0("Prot_mRNA_Cor_CL",CL,"_plot.pdf"))
pdf(outfile,width=5,height=5)
plot(c(-1,1),c(-1,1),type="l",lty=2,xlab="protein-mRNA correlation (non-PAD)",ylab="protein-mRNA correlation (PAD)")#,main=paste0("Confidence level = ",CL,"%"))

##################
### scatterplot symbol parameters
case_color = rgb(255,165,0,maxColorValue=255)
control_color = rgb(77,38,0,maxColorValue=255)
ns_col = "lightgrey"
###
sel = ((r_PAD_pval>=0.05)&(r_CTR_pval>=0.05))|(diff==0) # 5528 out of 5851
sel = sel&(!sel_target_all)
lines(r_CTR[sel],r_PAD[sel],type="p",pch=1,col=ns_col)

sel = (diff!=0)&(r_PAD>0)&(r_PAD_pval<0.05)&(r_CTR<0)&(r_CTR_pval<0.05) # 1
sel = sel&(!sel_target_all)
lines(r_CTR[sel],r_PAD[sel],type="p",pch=16,col=case_color)

sel = (diff!=0)&(r_CTR>0)&(r_CTR_pval<0.05)&(r_PAD<0)&(r_PAD_pval<0.05) # 2
sel = sel&(!sel_target_all)
lines(r_CTR[sel],r_PAD[sel],type="p",pch=16,col=control_color)

sel = (diff!=0)&(r_PAD>0)&(r_PAD_pval<0.05)&(r_CTR_pval>=0.05) # 56
sel = sel&(!sel_target_all)
lines(r_CTR[sel],r_PAD[sel],type="p",pch=16,col=case_color)

sel = (diff!=0)&(r_CTR>0)&(r_CTR_pval<0.05)&(r_PAD_pval>=0.05) # 101
sel = sel&(!sel_target_all)
lines(r_CTR[sel],r_PAD[sel],type="p",pch=16,col=control_color)

sel = (diff!=0)&(r_PAD_pval>=0.05)&(r_CTR<0)&(r_CTR_pval<0.05) # 71
sel = sel&(!sel_target_all)
lines(r_CTR[sel],r_PAD[sel],type="p",pch=16,col=case_color)

sel = (diff!=0)&(r_CTR_pval>=0.05)&(r_PAD<0)&(r_PAD_pval<0.05)
sel = sel&(!sel_target_all)
lines(r_CTR[sel],r_PAD[sel],type="p",pch=16,col=control_color) # 92

for (i_file in 1:n_file) {
    lines(r_CTR[sel_target[[i_file]]],r_PAD[sel_target[[i_file]]],type="p",pch=16,col=col[i_file])
}

lines(c(-1,1),c(-1,1),type="l",lty=2)

# labels on target genes
label_gene_ENSG = c("ENSG00000100387","ENSG00000204922","ENSG00000123689")
label_gene_name = c("RBX1","UQCC3","G0S2")
index = match(label_gene_ENSG,data[,"ENSG"])
x = r_CTR[index]
y = r_PAD[index]
textxy(x,y,label_gene_name,cex=0.35)
dev.off()
