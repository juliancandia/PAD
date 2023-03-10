rm(list=ls())

PROJECT_DIR = "~/git/PAD-main" # replace with your local path

infile = file.path(PROJECT_DIR,"DATA","metadata_samples.csv")
metadata = read.table(infile,header=T,stringsAsFactors=F,sep=",")
sample_label = paste0("S",metadata[,"Code"])
n_sample = nrow(metadata)
group_label = metadata[,"PAD.non.PAD"]
group_label[group_label=="non-PAD"] = "CTR"
group_sample = paste0(group_label,"_",sample_label)
metadata = data.frame(sample_label,group_label,group_sample,metadata)

infile = file.path(PROJECT_DIR,"DATA","RSEM.genes.TPM.all_samples.txt")
data = as.matrix(read.table(infile,header=T,stringsAsFactors=F,sep="\t",quote=""))
ENSG = data[,"gene_id"]
gene = data[,"GeneName"]
sample = colnames(data)[-(1:2)]
sample = matrix(unlist(strsplit(sample,"_")),ncol=2,byrow=T)[,1]
expr = data[,-(1:2)]
expr = matrix(as.numeric(expr),ncol=ncol(expr))
rownames(expr) = gene
colnames(expr) = sample
o = match(metadata[,"sample_label"],colnames(expr))
expr = expr[,o]
classA_index = which(metadata[,"group_label"]=="PAD")
classB_index = which(metadata[,"group_label"]=="CTR")

case_color = rgb(255,165,0,maxColorValue=255)
control_color = rgb(77,38,0,maxColorValue=255)
transparency = 70
scatter_span = 0.35

# target gene
target = c("MIR210HG") #,"S100A8","S100A9")
target_label = target
n_target = length(target)
for (i_target in 1:n_target) {
    outfile = file.path(PROJECT_DIR,"RESULTS","FIG1",paste0("logTPM_Boxplot_",target_label[i_target],".pdf"))
    pdf(outfile,width=3.2,height=5)
    i_gene = which(gene==target[i_target])
    class = rep(NA,n_sample)
    class[metadata[,"group_label"]=="PAD"] = "1.PAD"
    class[metadata[,"group_label"]=="CTR"] = "2.non-PAD"
    mydata = data.frame(expr[i_gene,],class)
    colnames(mydata) = c("expr","class")
    wilcox_pval = wilcox.test(expr~class,data=mydata)$p.value
    boxplot(expr~class,data=mydata,xlab="", log="y",ylab=paste0(target_label[i_target]," (TPM)"),outline=F,boxwex=0.8,range=0)
    str = paste0("Wilcoxon p-val: ",signif(wilcox_pval,digits=2))
    #title(str,cex.main=0.85,line=0.5)
    x0 = 0.5
    y0 = max(mydata$expr)-(max(mydata$expr)-min(mydata$expr))*0.05
    #text(x=x0,y=y0,labels=paste0("p=",signif(wilcox_pval,digits=2)),adj=0)
    datapoints = expr[i_gene,class=="1.PAD"]
    scatter = runif(length(datapoints),-scatter_span,scatter_span)
    x = 1+scatter-mean(scatter)
    y = datapoints
    col_par = as.numeric(col2rgb(case_color))
    #points(x,y,cex=1.7,col=rgb(col_par[1],col_par[2],col_par[3],
    #transparency,maxColorValue=255),pch=16)
    points(x,y,cex=1.2,col=case_color,pch=16)
    datapoints = expr[i_gene,class=="2.non-PAD"]
    scatter = runif(length(datapoints),-scatter_span,scatter_span)
    x = 2+scatter-mean(scatter)
    y = datapoints
    col_par = as.numeric(col2rgb(control_color))
    #points(x,y,cex=1.7,col=rgb(col_par[1],col_par[2],col_par[3],
    #transparency,maxColorValue=255),pch=16)
    points(x,y,cex=1.2,col=control_color,pch=16)
    dev.off()
}

