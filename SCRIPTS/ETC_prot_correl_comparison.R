rm(list=ls())

PROJECT_DIR = "~/git/PAD-main" # replace with your local path

infile = file.path(PROJECT_DIR,"DATA","MitoProteins_pairwise_correlation.csv")
data = read.table(infile,header=T,stringsAsFactors=F,sep=",")

complex = c("AC","C_I","C_III","C_IV","C_V")
n_complex = length(complex)

res = matrix(rep(NA,n_complex*3),ncol=3)
for (i_complex in 1:n_complex) {
    r = data[,2*i_complex-1]
    r_nonPAD = r[!is.na(r)]
    quart = quantile(r_nonPAD,probs=c(0.25,0.5,0.75))
    res[i_complex,1] = paste0(signif(quart[2],digits=3)," (IQR: ",signif(quart[1],digits=3),"-",signif(quart[3],digits=3),")")
    r = data[,2*i_complex]
    r_PAD = r[!is.na(r)]
    quart = quantile(r_PAD,probs=c(0.25,0.5,0.75))
    res[i_complex,2] = paste0(signif(quart[2],digits=3)," (IQR: ",signif(quart[1],digits=3),"-",signif(quart[3],digits=3),")")
    res[i_complex,3] = signif(wilcox.test(r_nonPAD,r_PAD)$p.value,digits=3)
}

output = rbind(c("ETC Protein Class","non-PAD","PAD","p-value"),cbind(complex,res))
outfile = file.path(PROJECT_DIR,"RESULTS","FIG3","ETC_correl_comparison.txt")
write(t(output),ncol=ncol(output),file=outfile,sep="\t")
