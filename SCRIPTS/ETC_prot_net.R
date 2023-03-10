rm(list=ls())

PROJECT_DIR = "~/git/PAD-main" # replace with your local path

infile = file.path(PROJECT_DIR,"DATA","MitoProteins.txt")
data = read.table(infile,header=T,stringsAsFactors=F,sep="\t")
expr = data[,15:73]
sample = colnames(expr)
prot = data[,"Protein"]
gene = data[,"Gene"]
n_prot = length(prot)
expr = matrix(as.numeric(unlist(expr)),nrow=n_prot)
group = c("CTR","PAD")
n_group = length(group)
for (i_group in 1:n_group) {
    if (group[i_group]=="PAD") {
        index = grep('^PAD',sample)
    } else if (group[i_group]=="CTR") {
        index = grep('^non.PAD',sample)
    }
    res = NULL
    for (i_prot in 1:(n_prot-1)) {
        x = expr[i_prot,index]
        x.not.na = !is.na(x)
        for (j_prot in (i_prot+1):n_prot) {
            y = expr[j_prot,index]
            xy.not.na = x.not.na&(!is.na(y))
            if (sum(xy.not.na)>=5) {
                cor.test = cor.test(x[xy.not.na],y[xy.not.na])# ,method="spearman")
                pval = cor.test$p.value
                r = cor.test$estimate
                res_add = c(gene[i_prot],gene[j_prot],r,pval)
                res = rbind(res,res_add)
            }
        }
    }
    padj = p.adjust(as.numeric(res[,4]),method="fdr")
    output = rbind(c("prot1","prot2","r","p","padj"),cbind(res,padj))
    outfile = file.path(PROJECT_DIR,"RESULTS","FIG3",paste0("ETC_",group[i_group],"_net.txt"))
    write(t(output),ncol=ncol(output),file=outfile,sep="\t")
    
    sel = padj<0.05
    output = rbind(c("prot1","prot2","r"),cbind(res[sel,1:2],abs(as.numeric(res[sel,3]))))
    outfile = file.path(PROJECT_DIR,"RESULTS","FIG3",paste0("ETC_",group[i_group],"_net_filt.txt"))
    write(t(output),ncol=ncol(output),file=outfile,sep="\t")
}
outfile = file.path(PROJECT_DIR,"RESULTS","FIG3","ETC_nodes.txt")
gene = data[,"Gene"]
class = data[,"MitoDetails"]
dna = rep("nDNA",length(gene))
dna[grep("^MT",gene)] = "mtDNA"
output = rbind(c("prot","class","dna"),cbind(gene,class,dna))
write(t(output),ncol=ncol(output),file=outfile,sep="\t")
