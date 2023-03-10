require(RColorBrewer)
library(pheatmap)

PROJECT_DIR = "~/git/PAD-main" # replace with your local path

infile = file.path(PROJECT_DIR,"DATA","ETC_ProteinHeatmap.txt")
infileAnno = file.path(PROJECT_DIR,"DATA","AnnoFor96ECData.txt")

HeatmapProt<-read.delim(infile,row.names = 1,header=TRUE,check.names = FALSE)
## an annotation table for columns
AnnoEC<-read.delim(infileAnno,row.names = 1,header=TRUE,check.names = FALSE)

#annotation for the rows
annotation_row = data.frame(HeatmapProt[1])
#annotation colors for rows and columns
annotation_colors = list(MitoDetails=c('Assembly Complex'="#FFFF66",'Complex I'="#66FF66",'Complex II'="#0066FF",'Complex III'="#CC0066",'Complex IV'="#996600",'Complex V'="#FF9900"),Groups = c(nonPAD=rgb(77,38,0,maxColorValue=255), PAD=rgb(255,165,0,maxColorValue=255))) #"gray "#"red"))
HeatmapProt[is.na(HeatmapProt)] = 0

d1<-subset(HeatmapProt,HeatmapProt$Group=="nonPAD")
d1<-d1[,2:length(d1)]
d2<-subset(HeatmapProt,HeatmapProt$Group=="PAD")
d2<-d2[,2:length(d2)]


## an annotation table for rows
Groups <- c(rep("nonPAD",nrow(d1)), rep("PAD",nrow(d2)))
Groups <- as.data.frame(Groups)
rownames(Groups) <- c(rownames(d1), rownames(d2))

## Figure export
outfile = file.path(PROJECT_DIR,"RESULTS","FIG2","ETC_ProteinHeatmap_barplot.pdf")
pdf(outfile,width=6,height=4)

pheatmap(t(HeatmapProt[,2:length(HeatmapProt)]),cluster_cols=FALSE,type = colnames(Groups)[1],annotation_col = Groups,annotation_row = AnnoEC,scale = "none", conditions="Auto",annRow=Groups,border_color=NA,annotation_colors = annotation_colors,cutree_rows =8,show_rownames=FALSE,show_colnames=FALSE,col = colorRampPalette(c("navy", "white", "firebrick3"))(100))

dev.off()
