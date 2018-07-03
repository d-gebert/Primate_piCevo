library(gplots)

# Image parameters
imres = 600
imwid = 7.5
imhei = 10
impoi = 14

# Command line arguments
args = commandArgs(trailingOnly=TRUE)

# Read in data
expr <- read.table(args[1], header=TRUE)
# Save as matrix
mat <- as.matrix(expr)
# Clustering, average, pearson distance
hc <- hclust(as.dist(1-cor(t(mat))), method="average")

# Create image
png("heatmap.png", width = imwid*imres, height = imhei*imres, res = imres, pointsize = impoi)

# Draw heatmap
heatmap.2(mat,
          Rowv=as.dendrogram(hc),
          dendrogram="row",
          trace="none",
          Colv=NA,
          col=bluered(50),
          labRow="",
          scale="row",
          density.info="none",
          main = "",
          lhei=c(0.5,3), lwid=c(1,3)
)

# Close image
dev.off()

