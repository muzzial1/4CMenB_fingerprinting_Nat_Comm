#Created by Marco Bruttini, 2014. Please do NOT remove this line.
setwd("~/4CMenB_fingerprinting_Nat_Comm/")

require("ggplot2")
require("reshape2")
require("scales")
require("gridExtra")
require("ggdendro")
require("qvalue")

ss <- read.delim("./data/SS 8+ Paper.txt", row.names=3)
o <- order(ss$Filtro, ss$Inizio, - ss$Fine)
x <- read.delim(file="./data/Tutto.txt", row.names=1)[o,]
ss <- ss[o,]
rm(o)

for (ag in c("NadA", "NHBA", "fHbp 1")){
  title <- paste0("V72P4 + V72P10 + V72P12E1 ", ag)
  ok <- ss$Filtro == ag
  n <- length(which(ok))
  V72P4 <- x[ok, 1:30]
  V72P10 <- x[ok, 31:76]
  V72P12E1 <- x[ok, 77:145]
  
  P4P10 <- numeric(n)
  P10P12E1 <- numeric(n)
  P12E1P4 <- numeric(n)
  
  for (i in 1:n) {
    P4P10[i] <- wilcox.test(unlist(V72P4[i,]), unlist(V72P10[i,]),alternative="t")$p.value
    P10P12E1[i] <- wilcox.test(unlist(V72P10[i,]), unlist(V72P12E1[i,]),alternative="t")$p.value
    P12E1P4[i] <- wilcox.test(unlist(V72P12E1[i,]), unlist(V72P4[i,]),alternative="t")$p.value
  }
  
  P4P10[!is.finite(P4P10)] <- 1
  P10P12E1[!is.finite(P10P12E1)] <- 1
  P12E1P4[!is.finite(P12E1P4)] <- 1
  
  result <- data.frame(Frammenti=rownames(ss)[ok], P4P10, P10P12E1, P12E1P4)
  
  asterisks <- trunc(-log10(result[-1]))
  asterisks[asterisks>5] <- 5  
  result <- cbind(result, asterisks)
  
  qP4P10=qvalue(P4P10)$qvalues
  qP10P12E1=qvalue(P10P12E1)$qvalues
  qP12E1P4=qvalue(P12E1P4)$qvalues
  result <- cbind(result, qP4P10, qP10P12E1, qP12E1P4)
  
  write.table(result, paste0("./Figure 2/",title, " Wilcoxon.txt"), quote=F, sep="\t", row.names=F)
  message(title, " Wilcoxon")
}
