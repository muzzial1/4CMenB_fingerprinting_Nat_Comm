#Created by Marco Bruttini, 2014. Please do NOT remove this line.

require("ggplot2")
require("reshape2")
require("scales")
require("grid")
require("gridExtra")
require("ggdendro")

setwd("~/4CMenB_fingerprinting_Nat_Comm/")

ss <- read.delim("./data/SS 8+ Paper.txt")
x <- read.delim(file="./data/Tutto.txt")
x <- cbind(ss, x[-1])
x <- x[order(x$Filtro, x$Inizio, - x$Fine), ]
rownames(x) <- NULL
rm(ss)

x <- x[x$Filtro %in% c("NadA", "NHBA", "fHbp 1"), c(3:8, 2:1, 1, 9:ncol(x))]
names(x)[1] <- "Frammento"
names(x)[7] <- "Mediana"
names(x)[9] <- "Delta"
x$Delta <- sapply(as.character(x[[8]]),switch, "fHbp 1"=434, "NadA"=350, "NHBA"=644)
x$Frammento <- factor(x$Frammento, levels=rev(x$Frammento), ordered=T)
x[7] <- round(apply(x[-1:-9], 1, median), 0)
rownames(x) <- NULL

sba <- read.delim(file="SBA Tutto.txt")
l <- max(sba$Valore, na.rm=T) # 2^ceiling(log2(max(b$Valore)))

bb <- sba[sba$Visita != "1S", c(3, 5:7)]
bb$Ceppo[bb$Ceppo == "H44/76-SL"] <- "H44/76"
rownames(bb) <- NULL
rm(sba)

colors <- limits <- values <- breaks <- list()
colors$AGE <- c("violetred1","dodgerblue2","chartreuse2")#(c("#000066", "#0099FF", "#CCFFFF")) #c("white","grey50","black")
colors$MFI <- c("grey90", "yellow", "orange", "red", "red", "darkred") #lightgrey
colors$SBA <- c("grey90", "grey75", "black")#c("grey90", "darkolivegreen1", "darkgreen")
limits$MFI <- c(0, 65535)
limits$SBA <- log2(c(1, l))
values$MFI <- rescale(c(0, 8000, 19000, 40000, 55000, 65535), from=limits$MFI)
values$SBA <- c(0, 1/2, 1)
breaks$MFI <- c(0, 5000, 15000, 30000, 60000, 65535)
breaks$SBA <- c(1, l)

capitalize <- function(string) {
  string <- strsplit(string, " ")[[1]]
  paste(toupper(substring(string, 1, 1)), substring(string, 2), sep = "", collapse = " ")
}
  
#for (metodo in c("pearson", "kendall", "spearman", "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski")) {
for (metodo in c("manhattan")) {
  #for (metodocluster in c("ward.D", "single", "complete", "average", "mcquitty", "median", "centroid", "ward.D2")){
  for (metodocluster in c("complete")) {  
    title <- paste("Cluster Totale", metodo, metodocluster)
    title <- capitalize(title)
    
    if (metodo %in% c("pearson", "kendall", "spearman")) {
      hc <- hclust(as.dist(1-cor(as.matrix(x[-1:-9]), method=metodo)), method=metodocluster)  
    } else {
      if (metodo %in% c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski")) {
        hc <- hclust(dist(t(as.matrix(x[-1:-9])),method=metodo), method=metodocluster)
      } else {
        stop("Ma che metodo hai scritto???", call.=F)    
      }
    }
    
    sieri <- substring(hc$labels, 1, 3)
    col.ord <- order.dendrogram(as.dendrogram(hc))
    sieri <- substring(pvc$hclust$labels, 1, 3)
    col.ord <- order.dendrogram(as.dendrogram(pvc$hclust))
    xx <- x[c(1:9, col.ord+9)]
    bb$Soggetto <- factor(bb$Soggetto, unique(bb$Soggetto)[col.ord], ordered=T)
    gruppi<-cutree(hc,2)[col.ord]
    
    A=xx[c(1:9, which(gruppi==2)+9)]
    B=xx[c(1:9, which(gruppi==1)+9)]
    A$Mediana <- round(apply(A[-1:-9], 1, median), 0)
    B$Mediana <- round(apply(B[-1:-9], 1, median), 0)
    #write.csv2(file="Gruppi.csv", data.frame(A[1],A=A[[7]],B=B[[7]]), row.names=F)
    
    prova<-data.frame(row.names=unique(bb$Soggetto),"5/99"=bb$Valore[bb$Ceppo=="5/99"],"H44/76"=bb$Valore[bb$Ceppo=="H44/76"], gruppi)
    apply(prova,2,function(x) median(x,na.rm=T)) #Tutto
    
    by(prova$X5.99,gruppi,function(x) median(x, na.rm=T))
    by(prova$H44.76,gruppi,function(x) median(x, na.rm=T))
    
    apply(prova[1:72,],2,function(x) median(x,na.rm=T)) #A
    apply(prova[73:145,],2,function(x) median(x,na.rm=T)) #B
    apply(prova[1:12,],2,function(x) median(x,na.rm=T)) #A.1
    apply(prova[13:72,],2,function(x) median(x,na.rm=T)) #A.2
    apply(prova[73:78,],2,function(x) median(x,na.rm=T)) #B.1
    apply(prova[79:95,],2,function(x) median(x,na.rm=T)) #B.2
    apply(prova[96:120,],2,function(x) median(x,na.rm=T)) #B.3
    apply(prova[121:145,],2,function(x) median(x,na.rm=T)) #B.4    
    
    mdf <- melt(xx, id=1:9, variable.name="Soggetto", value.name="MFI")
    mdf$Trial <- factor(substr(mdf$Soggetto, 1, 3), c("INF", "ADO", "ADU"), ordered=T)
    n <- nrow(xx)
    m <- max(xx$Fine)
    
    pD <- ggdendrogram(hc) + theme(plot.margin=unit(c(0,0,-10,0),"mm"), axis.text=element_blank())
    
    pC <- ggplot(mapping=aes(x=1:length(sieri), y=1, fill=sieri[col.ord])) + geom_tile(colour="#336699") +
      scale_fill_manual(values=colors$AGE, limits=c("INF","ADO","ADU")) +
      scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) + theme_grey() + labs(x=NULL, y=NULL) +
      theme(legend.position="none", plot.margin=unit(c(0,0,0,0),"mm"), axis.ticks=element_blank(), axis.text=element_blank())
    
    pQ <- ggplot(mdf, aes(x=Soggetto, y=Frammento, fill=MFI)) + geom_tile(colour="white") +
      scale_fill_gradientn(colours=colors$MFI, values=values$MFI, limits=limits$MFI) +
      scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) + theme_grey(base_size=30) + labs(x=NULL, y=NULL) +
      theme(legend.position="none", plot.margin=unit(c(1,0,1,0),"mm"), axis.ticks=element_blank(), axis.text.x=element_blank(), axis.text.y=element_text(color="black", face="bold"))
    
    pB <- ggplot(bb, aes(x=Soggetto, y=Ceppo, fill=log2(Valore))) + geom_tile(colour="white") + #facet_grid(Ceppo ~ .) + 
      scale_fill_gradientn(limits=limits$SBA, values=values$SBA, colours=colors$SBA, na.value="white", guide="colourbar") +
      geom_text(aes(x=Soggetto, y=Ceppo, label=gsub(" NA", "", paste(Segno, Valore)), color=(Valore<l/5)), bb, size=8, fontface="bold", angle=-90)  + theme(
        plot.margin=unit(c(0,0,0,0),"mm"), legend.position="none", axis.ticks=element_blank(), strip.text=element_text(size=10),
        axis.text=element_blank(), panel.background=element_blank()) + labs(x=NULL, y=NULL) + scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
      scale_color_manual(values=c("white", "black"))  
###
    
    pD <- ggdendrogram(pvc$hclust) + theme(plot.margin=unit(c(0,0,-10,0),"mm"), axis.text=element_blank())
    
    pC <- ggplot(mapping=aes(x=1:length(sieri), y=1, fill=sieri[col.ord])) + geom_tile(colour="white", size=1) +
      scale_fill_manual(values=colors$AGE, limits=c("INF","ADO","ADU")) +
      scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) + theme_grey() + labs(x=NULL, y=NULL) +
      theme(legend.position="none", plot.margin=unit(c(0,0,0,0),"mm"), axis.ticks=element_blank(), axis.text=element_blank())
    
    pQ <- ggplot(mdf, aes(x=Soggetto, y=Frammento, fill=MFI)) + geom_tile(colour="white", size=1) +
      scale_fill_gradientn(colours=colors$MFI, values=values$MFI, limits=limits$MFI) +
      scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) + theme_grey(base_size=30) + labs(x=NULL, y=NULL) +
      theme(legend.position="none", plot.margin=unit(c(5,0,5,0),"mm"), axis.ticks=element_blank(), axis.text.x=element_blank(), axis.text.y=element_text(color="black", face="bold"))
    
    pB <- ggplot(bb, aes(x=Soggetto, y=Ceppo, fill=log2(Valore))) + geom_tile(colour="white", size=1) + #facet_grid(Ceppo ~ .) + 
      scale_fill_gradientn(limits=limits$SBA, values=values$SBA, colours=colors$SBA, na.value="white", guide="colourbar") +
      geom_text(aes(x=Soggetto, y=Ceppo, label=Segno), bb, color="black", size=8, fontface="bold", vjust=.5, hjust=.5)  + theme(
        plot.margin=unit(c(0,0,0,0),"mm"), legend.position="none", axis.ticks=element_blank(), strip.text=element_text(size=10),
        axis.text=element_blank(), panel.background=element_blank()) + labs(x=NULL, y=NULL) + scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) 
    
    gpD <- ggplotGrob(pD)
    gpC <- ggplotGrob(pC)
    gpQ <- ggplotGrob(pQ)
    gpB <- ggplotGrob(pB)
    
    maxWidth = grid::unit.pmax(gpD$widths[2:5], gpC$widths[2:5], gpQ$widths[2:5], gpB$widths[2:5])
    gpD$widths[2:5] <- as.list(maxWidth)
    gpC$widths[2:5] <- as.list(maxWidth)
    gpQ$widths[2:5] <- as.list(maxWidth)
    gpB$widths[2:5] <- as.list(maxWidth)
    
    png(filename="./Figure 3/Figur3 3.png", width=4000, height=4000) #5k x 4k 
    grid.arrange(gpD, gpC, gpQ, gpB, nrow=4, ncol=1, heights=c(200, 30, height=15*n + 30, 70))
    dev.off()
    message(title)
  }
}
