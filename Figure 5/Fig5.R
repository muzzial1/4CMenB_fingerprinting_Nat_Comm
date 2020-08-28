#Created by Marco Bruttini, 2020. Please do NOT remove this line.
setwd("~/4CMenB_fingerprinting_Nat_Comm/")

#Libs ----
require("ggplot2")
require("reshape2")
require("scales")
require("gridExtra")
require("ggdendro")
require("rgl")
require("ggfortify")
require("PerformanceAnalytics")

#Params ----
X936 <- 3:8
X936fHbp <- c(1, 3:20, 25:28)
XNadA <- 4
XNadA <- NULL

colors.MFI <- c("lightgray", "yellow", "orange", "red", "red", "darkred")
colors.SBA <- c("darkgreen", "gold", "red")
limits.MFI <- c(0, 65535)
values.MFI <- rescale(c(0, 8000, 19000, 40000, 55000, 65535), from=limits.MFI)
breaks.MFI <- c(0, 5000, 15000, 30000, 60000, 65535)

miocolfun <- function(X, minimum = 0, maximum = 100){
  ScaledX<-NULL
  #maximum<-max(X)
  #minimum<-min(X)
  for (i in 1:length(X))
  {
    ScaledX<-c(ScaledX,(X[i]-minimum)/(maximum-minimum))
  }
  return(rgb(0,ScaledX,1))
}

ss <- read.delim("./data/SS 8+ Paper.txt")
x <- read.delim(file="./data/Tutto.txt")
x <- cbind(ss, x[-1])
x <- x[order(x$Filtro, x$Inizio, - x$Fine), ]
rownames(x) <- NULL
rm(ss)
gruppi <- substr(names(x[-1:-8]),1,3)
Age <- c(rep("Adults",30), rep("Adolescents",46), rep("Infants",69))

sba <- read.delim(file="./data/SBA Tutto.txt")

ag<-"fHbp 1"
title<-"fHbp"
ag<-"NadA"
title<-"NadA"
saveCTRL <- F
if(saveCTRL == TRUE){pdf(file=paste0("./Figure 5/",title,".pdf", sep=''), paper = "a4")}
 xx <- x[x$Filtro == ag, c(3:8, 2, 9:ncol(x))]
if (ag=="fHbp 1") xx <- xx[-X936,]#xx[-X936fHbp,]
if (ag=="NadA") xx <- xx[-XNadA,]
 names(xx)[1] <- "Frammento"
 names(xx)[7] <- "Mediana"
 xx[7] <- round(apply(xx[-1:-7], 1, median), 0)
 names(xx)[-1:-7] <- substring(names(xx)[-1:-7], 5)
 rownames(xx) <- xx[[1]]
 f <- sub("\\[","", sub("].*", "", rownames(xx)))
 s <- names(xx[-1:-7])
n <- nrow(xx)
 xx$Frammento <- with(xx, factor(Frammento, levels=rev(Frammento), ordered=TRUE))
 bb <- sba[sba$Antigene == ag, 3:7]
 bb<-bb[bb$Visita!="1S",]
 bb$Soggetto <- factor(bb$Soggetto, unique(bb$Soggetto), ordered=T)
 rownames(bb) <- NULL

ind<-c(rep("ADU",30),rep("ADO",46),rep("INF",69))
by(bb$Valore,ind,function(x) median(x,na.rm=T))

#bbA<-bb[clusters==2,]
#bbB<-bb[clusters==1,]
#AgeA<-Age[clusters==2]
#AgeB<-Age[clusters==1]

#Rescaling ----
mioxx <- xx[c(-1:-7)]
dim(mioxx)
dim(mioxx)
mioyy <- log2(bb$Valore)
length(log2(bb$Valore))
cor.test(t(mioxx["A-02",]),mioyy, use = "n")
miorim <- which(is.na(mioyy))
mioxx <- mioxx[-miorim]
mioyy <- mioyy[-miorim]
length(mioyy)

#PLS pls ----
## VIP returns all VIP values for all variables and all number of components,
## as a ncomp x nvars matrix.
VIP <- function(object) {
  if (object$method != "oscorespls")
    stop("Only implemented for orthogonal scores algorithm.  Refit with 'method = \"oscorespls\"'")
  if (nrow(object$Yloadings) > 1)
    stop("Only implemented for single-response models")
  
  SS <- c(object$Yloadings)^2 * colSums(object$scores^2)
  Wnorm2 <- colSums(object$loading.weights^2)
  SSW <- sweep(object$loading.weights^2, 2, SS / Wnorm2, "*")
  sqrt(nrow(SSW) * apply(SSW, 1, cumsum) / cumsum(SS))
}


## VIPjh returns the VIP of variable j with h components
VIPjh <- function(object, j, h) {
  if (object$method != "oscorespls")
    stop("Only implemented for orthogonal scores algorithm.  Refit with 'method = \"oscorespls\"'")
  if (nrow(object$Yloadings) > 1)
    stop("Only implemented for single-response models")
  
  b <- c(object$Yloadings)[1:h]
  T <- object$scores[,1:h, drop = FALSE]
  SS <- b^2 * colSums(T^2)
  W <- object$loading.weights[,1:h, drop = FALSE]
  Wnorm2 <- colSums(W^2)
  sqrt(nrow(W) * sum(SS * W[j,]^2 / Wnorm2) / sum(SS))
}


#PLS caret ----
library(caret)
cor.test(t(mioxx["A-02",]),mioyy)

#pls training+test ----
miodf<-data.frame(t(mioxx),mioyy)
inTrain <- createDataPartition(
  y = miodf$mioyy,
  p = .8,
  list = FALSE
)
training <- miodf[ inTrain,]
testing  <- miodf[-inTrain,]
plsFit <- train(
  mioyy ~ .,
  data = training,
  method = "pls", 
  metric = "Rsquared",
  tuneLength = 20
)
summary(plsFit)
plot(plsFit, main = paste(title, "(PLS training regression)"))
textplot(capture.output(summary(plsFit)), halign="center", valign="center",cex=.5)
plot(varImp(plsFit), main = paste(title, "(PLS training regression)"))
varImp(plsFit)
textplot(capture.output(varImp(plsFit)), halign="center", valign="center",cex=.5)
mioheight <- as.data.frame(varImp(plsFit)$importance)$Overall
mionames <- rownames(as.data.frame(varImp(plsFit)$importance))
miocolbarplot <- miocolfun(as.data.frame(varImp(plsFit)$importance)$Overall)
miocolbarplot <- color.scale(100-as.data.frame(varImp(plsFit)$importance)$Overall,c(0,1,1),c(1,1,0),0)
barplot(height = mioheight[length(mioheight):1], 
        names = mionames[length(mionames):1],
        horiz=T , las=1, cex.names = .5, col = miocolbarplot[length(miocolbarplot):1],
        main = paste(title, "(PLS training regression)"),
        xlab = "Importance",
        ylab = paste(title, "fragments")
)
miodata <- data.frame(mx = gsub("`", "", mionames), my = as.numeric(mioheight))
miodata <- miodata[sort(miodata$mx, decreasing = T),]
ggplot(miodata, aes(x=mx, y=as.numeric(my))) +
  #geom_segment( aes(xend=x, yend=0)) +
  geom_segment(aes(x=mx, xend=mx, yend=0, y=as.numeric(my))) +
  geom_point( size=4, color="orange") +
  #geom_text(aes(x=mx, y=as.numeric(my), label=as.numeric(my))) +
  scale_x_discrete(limits = rev(levels(miodata$mx))) +
  scale_y_continuous() +
  coord_flip() +
  theme_bw() +
  xlab(paste(title,"fragment", sep=' ')) +
  ylab("feature importance") +
  labs(title = paste(title,"(PLS training regression)", sep=' '))



#pls ----
miomodpls <-train(t(mioxx), mioyy, method = "pls", metric = "Rsquared", tuneLength = 20
                  , trControl=trainControl(method="LOOCV")
                  )#metric = c('RMSE', 'Rsquared'); log2(mioxx + 1)
summary(miomodpls)
plot(miomodpls, main = paste(title,"(PLS regression)", sep=' '))
textplot(capture.output(summary(miomodpls)), halign="center", valign="center",cex=.5)
plot(varImp(miomodpls), main = paste(title,"(PLS regression)", sep=' '))
varImp(miomodpls)
textplot(capture.output(varImp(miomodpls)), halign="center", valign="center",cex=.5)
mioheight <- as.data.frame(varImp(miomodpls)$importance)$Overall
mionames <- rownames(as.data.frame(varImp(miomodpls)$importance))
miocolbarplot <- miocolfun(as.data.frame(varImp(miomodpls)$importance)$Overall)
miocolbarplot <- color.scale(100-as.data.frame(varImp(miomodpls)$importance)$Overall,c(0,1,1),c(1,1,0),0)
barplot(height = mioheight[length(mioheight):1], 
        names = mionames[length(mionames):1],
        horiz=T , las=1, cex.names = .5, col = miocolbarplot[length(miocolbarplot):1],
        main = paste(title, "(PLS regression)"),
        xlab = "Importance",
        ylab = paste(title, "fragments")
)
#mionames <- mionames[length(mionames):1]
#mioheight <- mioheight[length(mioheight):1]
miodata <- data.frame(mx = gsub("`", "", mionames), my = as.numeric(mioheight))
miodata <- miodata[sort(miodata$mx, decreasing = T),]
ggplot(miodata, aes(x=mx, y=as.numeric(my))) +
  #geom_segment( aes(xend=x, yend=0)) +
  geom_segment(aes(x=mx, xend=mx, yend=0, y=as.numeric(my))) +
  geom_point( size=4, color="orange") +
  #geom_text(aes(x=mx, y=as.numeric(my), label=as.numeric(my))) +
  scale_x_discrete(limits = rev(levels(miodata$mx))) +
  scale_y_continuous() +
  coord_flip() +
  theme_bw() +
  xlab(paste(title,"fragment", sep=' ')) +
  ylab("feature importance") +
  labs(title = paste(title,"(PLS regression)", sep=' '))

miocolorelist <- cbind(
  Age = c("Adolescents","Adults","Infants"),
  Col = c("blue","green","black")
)
miocol <- merge(Age,miocolorelist, by.x = "x" ,by.y = "Age")
pcxlab <- paste0("Comp 1 (",round(0, 1),"%)")
pcylab <- paste0("Comp 2 (",round(0, 1),"%)")
pctit <- paste(title,"functional correlation (PLS regression)")
plot(miomodpls$finalModel$scores, main = paste(title,"(PLS regression)")
     #, col = miocol[-miorim]$Col
     )
plot(log2(bb$Valore[-miorim]) ~ miomodpls$finalModel$scores[,1], xlab = pcxlab, ylab = "log2(SBA titer)", 
     main = pctit
     #, col = miocol[-miorim]$Col
       )
plsC <- ggplot() + geom_point(aes(x=miomodpls$finalModel$scores[,1], y=log2(bb$Valore[-miorim]), fill=Age[-miorim], color=Age[-miorim]), size=3, alpha=.9, shape=21) +
    theme_bw(base_size=16) + labs(x=pcxlab, y="Log2(SBA titer)") + theme(
    legend.position=c(1,0), legend.justification=c(1,0), 
    plot.margin=unit(c(0, 0, 0, 0), "mm"), panel.background=element_blank()) +
  scale_x_continuous(expand = c(.01, 0)) + scale_y_continuous(expand = c(.01, 0)) +
  scale_fill_manual(values=c("Adults"="chartreuse2","Adolescents"="dodgerblue2","Infants"="violetred1"),  breaks=c("Infants","Adolescents","Adults")) +
  scale_color_manual(values=c("Adults"="chartreuse4","Adolescents"="dodgerblue4","Infants"="violetred3"),  breaks=c("Infants","Adolescents","Adults")) +
  ggtitle(pctit)
#scale_fill_manual(values=c("Adults"="#CCFFFF","Adolescents"="#0099FF","Infants"="#000066"),  breaks=c("Infants","Adolescents","Adults")) +
#scale_color_manual(values=c("Adults"="#669999","Adolescents"="#0066CC","Infants"="#000033"),  breaks=c("Infants","Adolescents","Adults"))
print(plsC)

#save ----

if(saveCTRL == TRUE){
  dev.off()
  save.image(file=paste0("./Figure 5/",title,".RData"))
}
