library(ggplot2)

PROJECT_DIR = "~/git/PAD-main" # replace with your local path

infile = file.path(PROJECT_DIR, "DATA", "Ribosomal_Proteins_Barplot.txt")
#Read file
RiboAnnotation <-
  read.delim(infile, header = TRUE)
RiboAnnotation$Gene <-
  factor(RiboAnnotation$Gene, levels = RiboAnnotation$Gene)

## Figure export
outfile = file.path(PROJECT_DIR, "RESULTS", "FIG2", "Ribosome_barplot.pdf")
pdf(outfile, width = 6, height = 4)

ggplot(RiboAnnotation, aes(x = Gene, y = PADBeta, fill =
    RiboCategory)) +
  geom_bar(
    stat = "identity",
    position = "identity",
    width = .8,
    alpha = .3,
    colour = "black") +
  labs(x = "Ribosome Genes", y = "PAD Beta") + labs(fill = "A") +
  labs(title = "PAD Ribosome, p<0.05") +
  theme(
    axis.line = element_line(color = "gray38", size = .5),
    axis.ticks = element_line(size = .7),
    axis.ticks.length = unit(0.2, "cm"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    text = element_text(size = 14, face = "plain"),
    axis.text.x = element_text(angle = 90, hjust = 1)
  )
dev.off()
