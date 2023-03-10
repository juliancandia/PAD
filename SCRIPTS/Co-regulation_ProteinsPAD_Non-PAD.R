require(RColorBrewer)
library(pheatmap)
library("Hmisc")

PROJECT_DIR = "~/git/PAD-main" # replace with your local path

infile = file.path(PROJECT_DIR,"DATA","Co-regulation_Proteins.txt")

# read the co-regulation complex proteins
HeatmapProt<-read.delim(infile,row.names = 1,header=TRUE,check.names = FALSE)

##PAD Co-regulation Analysis
my_data<-HeatmapProt[,!grepl("non-PAD", colnames(HeatmapProt))]
my_data<-my_data[,2:length(my_data)]

res2<-NULL
res2 <- rcorr(as.matrix(t(my_data)),type = c("pearson"))


#annotation colors for rows and columns

ann_colors = list(group=c('Assembly Complex'="#FFFF66",'Complex I'="#66FF66",'Complex II'="#0066FF",'Complex III'="#CC0066",
                          'Complex IV'="#996600",'Complex V'="#FF9900"))

colors<-colorRampPalette(rev(brewer.pal(n=7,name="RdYlBu")))(255)

##Anottation table for mitodetails

col_groups <- HeatmapProt$MitoDetails

# Data frame with column annotations.
mat_col <- data.frame(group = col_groups)

## Figure export
outfile = file.path(PROJECT_DIR,"RESULTS","FIG3","Co-regulation_Proteins_PAD.pdf")
pdf(outfile,width=5,height=5)

pheatmap(res2$r,scale='none',annotation_row = mat_col,annotation_legend=F,annotation_colors = ann_colors,annotation_col    = mat_col,border_color= NA, fontsize_row = 7.8, fontsize_col = 7.8,col=colors,show_rownames = F, show_colnames = F, main = "PAD, 31 Donors, 96 ETC proteins")

dev.off()

## Non-PAD Co-regulation Analysis

my_data<-HeatmapProt[,grepl("non-PAD", colnames(HeatmapProt))]
my_data<-my_data[,2:length(my_data)]

res2<-NULL
res2 <- rcorr(as.matrix(t(my_data)),type = c("pearson"))

#annotation colors for rows and columns
ann_colors = list(group=c('Assembly Complex'="#FFFF66",'Complex I'="#66FF66",'Complex II'="#0066FF",'Complex III'="#CC0066",
                          'Complex IV'="#996600",'Complex V'="#FF9900"))

colors<-colorRampPalette(rev(brewer.pal(n=7,name="RdYlBu")))(255)

col_groups <- HeatmapProt$MitoDetails

# Data frame with column annotations.
mat_col <- data.frame(group = col_groups)

## Figure export
outfile = file.path(PROJECT_DIR,"RESULTS","FIG3","Co-regulation_Proteins_non-PAD.pdf")
pdf(outfile,width=5,height=5)
pheatmap(res2$r,scale='none',annotation_row = mat_col,annotation_legend=F,annotation_colors = ann_colors,annotation_col    = mat_col,border_color= NA, fontsize_row = 7.8, fontsize_col = 7.8,col=colors,show_rownames = F, show_colnames = F, main = "Non-PAD, 28 Donors, 96 ETC proteins")
dev.off()
