library(data.table)
library(ggplot2)

# transform your *.tab files using cat hg38.tab | awk '{for (f=1;f<=NF;f++) col[f] = col[f]":"$f} END {for (f=1;f<=NF;f++) print col[f]}'| tr ':' ';' | cut -c 2- > hg38.csv



file <- "/Users/Anne/Desktop/tRNAEvo/rheMac3.csv"

pdf("/Users/Anne/Desktop/tRNAEvo/rheMac3.pdf",width=12, height=5)

dt <- as.data.table(read.csv(file, sep = ";", header = F))
dt.names <- as.character(unlist(dt[,1]))
dt <- as.data.table(t(dt)[-1,])
setnames(dt, dt.names)
dt <- melt(dt, measure = dt.names, na.rm = T)
dt[, value := as.numeric(value)]

p1 <- ggplot(dt, aes(x=variable, y = value, group = variable, colour = variable)) +
  geom_point(alpha = .2, size = 3) +
  labs( y = "absolute value", x = "") +
  theme(legend.position="none", plot.title = element_text(hjust = 0.5))

p2 <- ggplot(dt, aes(x=variable, y = value, group = variable, colour = variable)) +
  geom_boxplot()+
  labs( y = "absolute value", x = "") +
  theme(legend.position="none", plot.title = element_text(hjust = 0.5))

grid.arrange(p1,p2, ncol=2,nrow=1)
dev.off()



#log
pdf("/Users/Anne/Desktop/tRNAEvo/rheMac3_log.pdf",width=12, height=5)

dt <- as.data.table(read.csv(file, sep = ";", header = F))
dt.names <- as.character(unlist(dt[,1]))
dt <- as.data.table(t(dt)[-1,])
setnames(dt, dt.names)
dt <- melt(dt, measure = dt.names, na.rm = T)
dt[, value := as.numeric(value)]
dt[, value := -log(value)] #used for log 


p1 <- ggplot(dt, aes(x=variable, y = value, group = variable, colour = variable)) +
  geom_point(alpha = .2, size = 3) +
  labs( y = "-log(e-value)", x = "") +
  theme(legend.position="none", plot.title = element_text(hjust = 0.5))

p2 <- ggplot(dt, aes(x=variable, y = value, group = variable, colour = variable)) +
  geom_boxplot()+
  labs( y = "-log(e-value)", x = "") +
  theme(legend.position="none", plot.title = element_text(hjust = 0.5))

grid.arrange(p1,p2, ncol=2,nrow=1)
dev.off()

