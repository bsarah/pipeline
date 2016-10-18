require(qgraph)
library(igraph)


#myFiles <- list.files(pattern="*.gr")
#for(mf in myFiles){
args = commandArgs(trailingOnly=TRUE)
db=read.table(args[1])
db2 = db[,-3]
g<-graph.data.frame(db2)
undirG <- as.undirected(g,"collapse")
x<-c(args[1],"_plot.pdf")
pdf(paste(x,collapse=""),height=9,width=9)
par(mar=c(12,12,12,12),cex=0.8)
plot.igraph(undirG)
dev.off()
                                        #}

