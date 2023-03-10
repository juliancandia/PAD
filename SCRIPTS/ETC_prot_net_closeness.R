library(igraph)

rm(list=ls())

PROJECT_DIR = "~/git/PAD-main" # replace with your local path

infile = file.path(PROJECT_DIR,"DATA","MitoProteins.txt")
data = read.table(infile,header=T,stringsAsFactors=F,sep="\t")
expr = data[,15:73]
sample = colnames(expr)
n_sample = length(sample)
prot = data[,"Protein"]
gene = data[,"Gene"]
n_prot = length(prot)
expr = matrix(as.numeric(unlist(expr)),nrow=n_prot)
group = c("CTR","PAD")
n_group = length(group)
prot_class_label = c("Complex I","Complex III","Complex IV","Complex V","Assembly Complex","SC","Non-SC")
n_prot_class = length(prot_class_label)
prot_class = vector("list",n_prot_class)
for (i in 1:5) {
    prot_class[[i]] = which(data[,"MitoDetails"]==prot_class_label[i])
}
prot_class[[6]] = which(data[,"MitoDetails"]%in%c("Complex I","Complex III","Complex IV"))
prot_class[[7]] = which(data[,"MitoDetails"]%in%c("Complex II","Complex V"))
CC = matrix(rep(NA,n_prot_class*n_group),ncol=n_group)
colnames(CC) = group
CC_res = vector("list",n_group)
cor_thres = 5
for (i_group in 1:n_group) {
    CC_res[[i_group]] = vector("list",n_prot_class)
    if (group[i_group]=="PAD") {
        index = grep('^PAD',sample)
    } else if (group[i_group]=="CTR") {
        index = grep('^non.PAD',sample)
    }
    edge = NULL
    for (i_prot in 1:(n_prot-1)) {
        x = expr[i_prot,index]
        x.not.na = !is.na(x)
        for (j_prot in (i_prot+1):n_prot) {
            y = expr[j_prot,index]
            xy.not.na = x.not.na&(!is.na(y))
            if (sum(xy.not.na)>=cor_thres) {
                cor.test = cor.test(x[xy.not.na],y[xy.not.na])# ,method="spearman")
                pval = cor.test$p.value
                r = abs(cor.test$estimate)
                edge_add = c(gene[i_prot],gene[j_prot],r,pval)
                edge = rbind(edge,edge_add)
            }
        }
    }
    #padj = p.adjust(as.numeric(edge[,4]),method="fdr")
    #edge = edge[padj<0.05,]
    for (i in 1:n_prot_class) {
        # we select the subgraph
        sel = (edge[,1]%in%gene[prot_class[[i]]])&(edge[,2]%in%gene[prot_class[[i]]])
        g=graph.edgelist(edge[sel,1:2],directed=F)
        E(g)$weight=sqrt(2*(1-as.numeric(edge[sel,3])))
        CC_res[[i_group]][[i]] = closeness(g,normalized=T)
        CC[i,i_group] = mean(CC_res[[i_group]][[i]])
    }
}
CC_D = CC[,1]-CC[,2] # CTR-PAD
#output = rbind(c("Class",paste0("CC_",group),"CC_diff"),cbind(prot_class_label,CC,CC_D))
#outfile = file.path(PROJECT_DIR,"RESULTS","FIG3","CC_avg.txt")
#write(t(output),ncol=ncol(output),file=outfile,sep="\t")

# randomization:
# re-do by randomizing the expression matrix with fixed margins:
# 1. convert expression matrix to pseudo-counts
count = matrix(round(1.e3*(2**expr)),ncol=ncol(expr))
# expr_reconstructed = matrix(log2(count/1.e3),ncol=ncol(count)) # back-transformation
count[which(is.na(count),arr.ind=T)] = 0 # 7.7% are NA's.

library(vegan)

# mymat = matrix(round(100*runif(100)),ncol=10)
# mymat[2,5] = 0
# perm_f = permatfull(mymat, fixedmar = "both", times=2)
# perm_s = permatswap(mymat, fixedmar = "both", times=2)
n_rdm = 1000
set.seed(123) # for reproducibility
perm_s = permatswap(count, fixedmar = "both", times=n_rdm)
CC_D_rdm = matrix(rep(NA,n_prot_class*n_rdm),ncol=n_rdm)
for (i_rdm in 1:n_rdm) {
    perm_s$perm[[i_rdm]][which(perm_s$perm[[i_rdm]]==0,arr.ind=T)] = NA
    expr_rdm = matrix(log2(perm_s$perm[[i_rdm]]/1.e3),ncol=n_sample)
    CC_rdm = matrix(rep(NA,n_prot_class*n_group),ncol=n_group)
    for (i_group in 1:n_group) {
        if (group[i_group]=="PAD") {
            index = grep('^PAD',sample)
        } else if (group[i_group]=="CTR") {
            index = grep('^non.PAD',sample)
        }
        edge = NULL
        for (i_prot in 1:(n_prot-1)) {
            x = expr_rdm[i_prot,index]
            x.not.na = !is.na(x)
            for (j_prot in (i_prot+1):n_prot) {
                y = expr_rdm[j_prot,index]
                xy.not.na = x.not.na&(!is.na(y))
                if (sum(xy.not.na)>=cor_thres) {
                    cor.test = cor.test(x[xy.not.na],y[xy.not.na])# ,method="spearman")
                    pval = cor.test$p.value
                    r = abs(cor.test$estimate)
                    edge_add = c(gene[i_prot],gene[j_prot],r,pval)
                    edge = rbind(edge,edge_add)
                }
            }
        }
        #padj = p.adjust(as.numeric(edge[,4]),method="fdr")
        #edge = edge[padj<0.05,]
        for (i in 1:n_prot_class) {
            # we select the subgraph
            sel = (edge[,1]%in%gene[prot_class[[i]]])&(edge[,2]%in%gene[prot_class[[i]]])
            g=graph.edgelist(edge[sel,1:2],directed=F)
            E(g)$weight=sqrt(2*(1-as.numeric(edge[sel,3])))
            CC_rdm[i,i_group] = mean(closeness(g,normalized=T))
        }
    }
    CC_D_rdm[,i_rdm] = CC_rdm[,2]-CC_rdm[,1]
}
CC_D_pval = rep(NA,n_prot_class)
for (i in 1:n_prot_class) {
    CC_D_pval[i] = (sum(abs(CC_D_rdm[i,])>abs(CC_D[i]))+1)/(n_rdm+1)
}
output = rbind(c("Class",paste0("CC_",group),"CC_diff","CC_diff_pval"),cbind(prot_class_label,CC,CC_D,CC_D_pval))
outfile = file.path(PROJECT_DIR,"RESULTS","FIG3","CC_avg.txt")
write(t(output),ncol=ncol(output),file=outfile,sep="\t")
