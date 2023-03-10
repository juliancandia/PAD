library(RColorBrewer)

rm(list=ls())

PROJECT_DIR = "~/git/PAD-main" # replace with your local path

CL = 95

# gene expression vs protein abundance correlations
infile = file.path(PROJECT_DIR,"RESULTS","FIG5",paste0("Prot_mRNA_Cor_CL",CL,"_pval.txt"))
data = as.matrix(read.table(infile,stringsAsFactors=F,header=T,sep="\t",quote="",comment.char=""))
r_CTR = as.numeric(data[,"r_CTR"])
r_PAD = as.numeric(data[,"r_PAD"])
diff = as.numeric(data[,"diff"])
r_CTR_pval = as.numeric(data[,"r_CTR_pval"])
r_PAD_pval = as.numeric(data[,"r_PAD_pval"])

group = rep("",nrow(data))

sel = (diff!=0)&(r_PAD>0)&(r_PAD_pval<0.05)&(r_CTR<0)&(r_CTR_pval<0.05)
group[sel] = "PAD_pos_CTR_neg"

sel = (diff!=0)&(r_CTR>0)&(r_CTR_pval<0.05)&(r_PAD<0)&(r_PAD_pval<0.05)
group[sel] = "CTR_pos_PAD_neg"

sel = (diff!=0)&(r_PAD>0)&(r_PAD_pval<0.05)&(r_CTR_pval>=0.05)
group[sel] = "PAD_pos_CTR_ns"

sel = (diff!=0)&(r_CTR>0)&(r_CTR_pval<0.05)&(r_PAD_pval>=0.05)
group[sel] = "CTR_pos_PAD_ns"

sel = (diff!=0)&(r_PAD_pval>=0.05)&(r_CTR<0)&(r_CTR_pval<0.05)
group[sel] = "PAD_ns_CTR_neg"

sel = (diff!=0)&(r_CTR_pval>=0.05)&(r_PAD<0)&(r_PAD_pval<0.05)
group[sel] = "CTR_ns_PAD_neg"

output = rbind(c(colnames(data),"group"),cbind(data,group))
outfile = file.path(PROJECT_DIR,"RESULTS","FIG5",paste0("Prot_mRNA_Cor_CL",CL,"_pval_groups.txt"))
write(t(output),ncol=ncol(output),sep="\t",file=outfile)
