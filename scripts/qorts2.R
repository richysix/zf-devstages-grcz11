library(QoRTs)
res <- read.qc.results.data("", decoder.files="decoder.tsv", calc.DESeq2=TRUE, calc.edgeR=TRUE)
makeMultiPlot.all(res, outfile.dir="output/", plot.device.name="pdf")
get.size.factors(res, outfile="output/sizeFactors.GEO.txt")
