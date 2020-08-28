#Created by Marco Bruttini, 2020. Please do NOT remove this line.
setwd("~/4CMenB_fingerprinting_Nat_Comm/")

require(ggplot2)
require(reshape2)
load("./data/MFI Tutto.RData")#Proteinchip data
load("./data/SBA Tutto.RData")#SBA data

x <- dcast(mx,Soggetto + Visita ~ Frammento, value.var = "MFI")

Adulti <- mx$Trial=="V72P4"
Adolescenti <- mx$Trial=="V72P10"
Bambini <- mx$Trial=="V72P12E1"
Post <- mx$Visita=="POST"

Scelti_F <- mx$Frammento %in% paste0("A-", formatC(c(1:2, 9:28), width=2, flag="0"))
Scelti_F <- mx$Frammento %in% paste0("B-", formatC(1:41, width=2, flag="0"))
Scelti_F <- mx$Frammento %in% paste0("C-", formatC(1:31, width=2, flag="0"))

load("./data/Risultati.RData")#Subjects clustering
clusters <- c(rep(1,10),0,rep(2,14),rep(3,45),rep(4,75))
Sieri <- data.frame(Siero=bsc$labels[bsc$order], Cluster=clusters, stringsAsFactors=F)
Sieri$Siero <- formatC(as.integer(substr(Sieri$Siero,5,100)), width=6, flag="0")

Scelti_S <- mx$Soggetto %in% Sieri$Siero[Sieri$Cluster==1]
Scelti_S <- mx$Soggetto %in% Sieri$Siero[Sieri$Cluster==2]
Scelti_S <- mx$Soggetto %in% Sieri$Siero[Sieri$Cluster==3]
Scelti_S <- mx$Soggetto %in% Sieri$Siero[Sieri$Cluster==4]
mfi <- mx$MFI[Scelti_S & Scelti_F & Post]
plot(sort(mfi))
mean(mfi)
1/mean(1/(mfi))

plot(log2(1+mfi))
boxplot.stats(mfi)
hist(mfi, breaks=13)

Scelti <- paste0("A-", formatC(c(1:2, 9:28), width=2, flag="0"))
Scelti <- paste0("B-", formatC(1:42, width=2, flag="0"))
Scelti <- paste0("C-", formatC(1:31, width=2, flag="0"))

Scelti_F <- mx$Frammento %in% Scelti
popolazioni <- list()

for (i in 1:4) {
  message("Cluster ", i)
  Scelti_S <- mx$Soggetto %in% Sieri$Siero[Sieri$Cluster==i]
  mfi <- mx$MFI[Scelti_S & Scelti_F & Post]
  popolazioni <- c(popolazioni, list(mfi))
  plot(log2(1+sort(mfi)), col=i)
  #Mediane <- numeric(length(Scelti))
  #for(j in 1:length(Scelti)) {
  #  Scelti_F <- mx$Frammento == Scelti[j]
  #  mfi <- mx$MFI[Scelti_S & Scelti_F & Post]
  #  Mediane[j] <- median(mfi)
  #  plot(sort(mfi))
  #}
  
  #message(median(Mediane))
  #plot(sort(Mediane))
  
  #sigmax <- sqrt(mean((mfi - mean(mfi))^2))
  #message(mean((mfi - mean(mfi))^3/sigmax^3))
  #plot(sort(mfi))
  
  message("\tMediana: ", median(mfi))
  message("\tMedia  : ", mean(mfi), "\n")
}

test <- data.frame(A=c(1,1,1,2,2,3), B=c(2,3,4,3,4,4), p.two.tailed=numeric(6), p.greater=numeric(6), p.lesser=numeric(6))
for (i in 1:nrow(test)) {
  test$p.two.tailed[i] <- wilcox.test(popolazioni[[test$A[i]]], popolazioni[[test$B[i]]])$p.value
  test$p.greater[i] <- wilcox.test(popolazioni[[test$A[i]]], popolazioni[[test$B[i]]], alternative="g")$p.value
  test$p.lesser[i] <- wilcox.test(popolazioni[[test$A[i]]], popolazioni[[test$B[i]]], alternative="l")$p.value
}
test$responso <- paste("\" Cluster", test$A, ifelse(test$p.two.tailed < .05, ifelse(test$p.greater < 0.5, ">", "<"), "="), "Cluster", test$B, "\"")
test


sba.boni <- sba[sba$Visita=="POST" & sba$Ceppo=="H44/76",]
sba.boni <- sba[sba$Visita=="POST" & sba$Ceppo=="5/99",]
sba.boni <- sba.boni[bsc$order,]
for (i in 1:4) {
  message("Cluster ", i)
  Scelti_S <- sba.boni$Soggetto %in% Sieri$Siero[Sieri$Cluster==i]
  tit <- sba.boni$Valore[Scelti_S]
  plot(sort(log2(tit+1)), col=i)
  message(mean(tit, na.rm = T))
  popolazioni[[i]] <- tit
}

test <- data.frame(A=c(1,1,1,2,2,3), B=c(2,3,4,3,4,4), p.two.tailed=numeric(6), p.greater=numeric(6), p.lesser=numeric(6))
for (i in 1:nrow(test)) {
  test$p.two.tailed[i] <- wilcox.test(popolazioni[[test$A[i]]], popolazioni[[test$B[i]]])$p.value
  test$p.greater[i] <- wilcox.test(popolazioni[[test$A[i]]], popolazioni[[test$B[i]]], alternative="g")$p.value
  test$p.lesser[i] <- wilcox.test(popolazioni[[test$A[i]]], popolazioni[[test$B[i]]], alternative="l")$p.value
}
test$responso <- paste("\" Cluster", test$A, ifelse(test$p.two.tailed < .05, ifelse(test$p.greater < 0.5, ">", "<"), "="), "Cluster", test$B, "\"")
test
