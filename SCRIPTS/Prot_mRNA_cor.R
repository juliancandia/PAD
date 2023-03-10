library(readxl)
library(biomaRt)

rm(list=ls())

PROJECT_DIR = "~/git/PAD-main" # replace with your local path

CL = 95

# protein expression data
#infile = file.path(PROJECT_DIR,"DATA","Ceereena_210804","PAD_NonPad_59_NzldProtein_Annotated_ForJulian_V1.xlsx")
infile = file.path(PROJECT_DIR,"DATA","Protein_abundance.xlsx")
data = as.matrix(read_excel(infile,sheet=1))
prot_sample = colnames(data)[9:67]
prot_annot = data[,1:8]
prot_expr = data[,-(1:8)]
prot_expr = matrix(as.numeric(prot_expr),ncol=ncol(prot_expr))

# gene expression data
infile = file.path(PROJECT_DIR,"DATA","Gene_Expression.txt")
data = as.matrix(read.table(infile,stringsAsFactors=F,header=T,sep="\t",quote=""))
gene_sample = colnames(data)[-(1:5)]
gene = data[,3]
gene_expr = data[,-(1:5)]
gene_expr = matrix(as.numeric(gene_expr),ncol=ncol(gene_expr))

# gene annotations
infile = file.path(PROJECT_DIR,"DATA","RawCountFile_RSEM_genes_wANNOT.txt")
gene_annot = as.matrix(read.table(infile,stringsAsFactors=F,header=T,sep="\t",quote=""))[,1:6]
gene_annot = gene_annot[gene_annot[,"pass_filter"]=="Y",]
cat("Checkpoint passed:",identical(gene_annot[,"name"],gene),"\n")
index_remove = grep("_PAR_Y",gene_annot[,"ENSG.ver"]) # "PAR_Y" suffixes (transcript counts are identical)
gene_expr = gene_expr[-index_remove,]
gene_annot = gene_annot[-index_remove,c(1,3,4,5)]

# ensembl/uniprot annotations mapping
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",dataset = "hsapiens_gene_ensembl",host='https://useast.ensembl.org')
map <- biomaRt::getBM(attributes = c("ensembl_gene_id","uniprot_gn_symbol","uniprot_gn_id"), mart=mart,filters="ensembl_gene_id",values=unique(gene_annot[,"ENSG"])) 
map = map[map[,"uniprot_gn_id"]%in%prot_annot[,"Entry"],]

sel = prot_annot[,"Entry"]%in%map[,"uniprot_gn_id"]
prot_annot = prot_annot[sel,]
prot_expr = prot_expr[sel,]

# we align samples in both datasets
infile = file.path(PROJECT_DIR,"DATA","metadata_samples.csv")
sample_metadata = read.table(infile,header=T,stringsAsFactors=F,sep=",")
sample_label = paste0("S",sample_metadata[,"Code"])
group_label = sample_metadata[,"PAD.non.PAD"]
group_label[group_label=="non-PAD"] = "CTR"
group_sample = paste0(group_label,"_",sample_label)
sample_metadata = data.frame(sample_label,group_label,group_sample,sample_metadata)
sample_metadata[,"Original.ID"] = trimws(sample_metadata[,"Original.ID"])

o = match(prot_sample,sample_metadata[,"Original.ID"])
sample_metadata = sample_metadata[o,]
cat("Checkpoint passed:",identical(sample_metadata[,"Original.ID"],prot_sample),"\n")

o = match(sample_metadata[,"sample_label"],gene_sample)
gene_sample = gene_sample[o]
cat("Checkpoint passed:",identical(sample_metadata[,"sample_label"],gene_sample),"\n")
gene_expr = gene_expr[,o]

sel_CTR = sample_metadata[,"group_label"]=="CTR"
sel_PAD = sample_metadata[,"group_label"]=="PAD"

n_prot = nrow(prot_annot)
res = NULL
res_header = c("r_CTR_lower","r_CTR","r_CTR_upper","r_PAD_lower","r_PAD","r_PAD_upper","diff","r_CTR_pval","r_PAD_pval")
res_header = c(colnames(prot_annot),colnames(gene_annot),res_header)
for (i_prot in 1:n_prot) {
    ENSG_all = unique(map[map[,"uniprot_gn_id"]==prot_annot[i_prot,"Entry"],"ensembl_gene_id"])
    for (ENSG in ENSG_all) {
        i_gene = which(gene_annot[,"ENSG"]==ENSG)
        res_annot = c(prot_annot[i_prot,],gene_annot[i_gene,]) 
        notNA_sel = (!is.na(prot_expr[i_prot,]))&(!is.na(gene_expr[i_gene,]))
        sel = notNA_sel&sel_CTR
        cor_CTR = cor.test(prot_expr[i_prot,sel],gene_expr[i_gene,sel],conf.level=CL/100)
        sel = notNA_sel&sel_PAD
        cor_PAD = cor.test(prot_expr[i_prot,sel],gene_expr[i_gene,sel],conf.level=CL/100)
        res_cor = c(cor_CTR$conf.int[1],cor_CTR$estimate,cor_CTR$conf.int[2],cor_PAD$conf.int[1],cor_PAD$estimate,cor_PAD$conf.int[2],"0",cor_CTR$p.value,cor_PAD$p.value)
        if (cor_PAD$estimate>cor_CTR$conf.int[2]) {
            res_cor[7] = "1"
        } else if (cor_PAD$estimate<cor_CTR$conf.int[1]) {
            res_cor[7] = "-1"
        }
        res = rbind(res,c(res_annot,res_cor))
    }
}
outfile = file.path(PROJECT_DIR,"RESULTS","FIG5",paste0("Prot_mRNA_Cor_CL",CL,"_pval.txt"))
output = rbind(res_header,res)
write(t(output),ncol=ncol(output),file=outfile,sep="\t")
