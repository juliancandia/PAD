library(readxl)

rm(list=ls())

PROJECT_DIR = "~/git/PAD-main" # replace with your local path

CL = 95

# gene expression vs protein abundance correlations
infile = file.path(PROJECT_DIR,"RESULTS","FIG5",
paste0("Prot_mRNA_Cor_CL",CL,"_pval_groups.txt"))
data = as.matrix(read.table(infile,
stringsAsFactors=F,header=T,sep="\t",quote="",comment.char=""))

file = c("CTR_mitochondrion","CTR_translation","PAD_extracellular_exosome")
title = c("Mitochondrion (GO Cellular Component)","tRNA Aminoacylation (Reactome)","Extracellular exosome (GO Cellular Component)")
enriched = c("CTR","CTR","PAD")
col = c("red","cyan2","blue")
label_shift = c(0.25,0.1,0.25)
n_file = length(file)
for (i_file in 1:n_file) {
    infile = file.path(PROJECT_DIR,"RESULTS","FIG5",paste0(file[i_file],".xlsx"))
    target_gene = as.matrix(read_excel(infile,sheet=1,col_names=F))[,2]
    if (enriched[i_file]=="CTR") {
        sel = data[,"group"]%in%c("CTR_ns_PAD_neg","CTR_pos_PAD_neg","CTR_pos_PAD_ns")
    } else if (enriched[i_file]=="PAD") {
        sel = data[,"group"]%in%c("PAD_ns_CTR_neg","PAD_pos_CTR_neg","PAD_pos_CTR_ns")
    }
    target_gene = data[sel&(data[,"ENSG"]%in%target_gene),]
    x = 1:nrow(target_gene)
    r_CTR = as.numeric(target_gene[,"r_CTR"])
    r_CTR_lo = as.numeric(target_gene[,"r_CTR_lower"])
    r_CTR_hi = as.numeric(target_gene[,"r_CTR_upper"])
    r_PAD = as.numeric(target_gene[,"r_PAD"])
    outfile = file.path(PROJECT_DIR,"RESULTS","FIG5",paste0(file[i_file],".pdf"))
    pdf(outfile,width=8,height=3.5)
    plot(range(x),c(-1,1),type="n",xaxt="n",yaxt="n",xlab="",ylab="protein-mRNA correlation",cex.main=0.85)
    #axis(1,at=x,labels=target_gene[,"name"],cex.axis=0.7,las=2)
    axis(2,at=seq(-1,1,by=0.5),labels=T,cex.axis=0.7)
    axis(1,at=x,labels=F,cex.axis=0.7,las=2)
    text(x+label_shift[i_file], par("usr")[3]-0.35, labels=target_gene[,"name"] , srt=45, pos=2, xpd = TRUE)
    
    lines(x,r_CTR,type="p",pch=15,col=col[i_file])
    lines(x,r_PAD,type="p",pch=4,cex=1.5,lwd=2,col=col[i_file])
    arrows(x,r_CTR_lo,x,r_CTR_hi,col=col[i_file],length=0.025,angle=90,code=3)
    abline(h=0,lty=2)
    dev.off()
}
