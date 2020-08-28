#Created by Marco Bruttini, 2014. Please do NOT remove this line.
setwd("~/4CMenB_fingerprinting_Nat_Comm/")

require("ggplot2")
require("scales")
require("grid")

ags <- c("fHbp", "NadA", "NHBA")
age <- c("Adults", "Adolescents","Adolescents + OMV", "Infants")
age2 <- paste(age[-4], "II")
agex <- c(age2, paste(c(age, age2), "senza 936"))
titles <- c(paste(rep(ags, 4), rep(age, each=3)), paste(rep(ags[3], 2), age2[2:3]), paste(ags[1], agex))
ags <- c("936-fHbp", "NadA", "NHBA-953")
ref <- data.frame(name=c(rep(ags,4),rep(ags[3],2),rep(ags[1],3),rep("fHbp",7)),from=1,to=c(rep(c(434,350,644),4),rep(644,2),rep(434,3),rep(249,7)),line=c(rep(c(181,NA,471),4),rep(471,2),rep(181,3),rep(NA,7)))
rm(ags, age, age2, agex)
colori <- regmatches(titles,regexpr("Adults|Adolescents|Infants",titles))
colori[colori=="Adults"] <- "chartreuse3"
colori[colori=="Adolescents"] <- "dodgerblue2"
colori[colori=="Infants"] <- "violetred1"
colorini <- colori
colorini[colori=="chartreuse3"] <- "darkolivegreen2"
colorini[colori=="dodgerblue2"] <- "deepskyblue1"
colorini[colori=="violetred1"] <- "pink1"

i <- 0
scelti <- c("fHbp Adolescents + OMV II senza 936", "fHbp Adults II senza 936", "fHbp Infants senza 936", "NadA Adolescents + OMV", "NadA Adults", "NadA Infants", "NHBA Adolescents + OMV", "NHBA Adults", "NHBA Infants")
  
for (t in titles) { 
  i <- i + 1
  if (t %in% scelti) {
    x <- read.delim(file.path("./data/Dati", paste0(t, ".txt")))
    names(x) <- c("Start", "Stop")
    v <- integer(ref$to[i])
    vshort <- v
    
    for (j in 1:nrow(x)){
      v[x$Start[j]:x$Stop[j]] <- v[x$Start[j]:x$Stop[j]] + 1
      if (x$Stop[j]-x$Start[j] <= 100)
        vshort[x$Start[j]:x$Stop[j]] <- vshort[x$Start[j]:x$Stop[j]] + 1
    }
    
    v <- v / nrow(x)
    vshort <- vshort / nrow(x)
    
    vshort <- vshort/max(v)
    v <- v/max(v)
    
    line <- ref$line[i]
    bre <- c(1, 1:(ref$to[i] %/% 100)*100, ref$to[i]) #, line
    name <- as.character(ref$name[i])
    
    p <- ggplot() + geom_vline(xintercept=line, linetype="dashed", size=.25) +
      geom_path(aes(x=1:length(vshort), y=vshort), color="grey75", size=.75, show_guide = F) + #colorini[i]
      geom_path(aes(x=1:length(v), y=v), size=.75, show_guide = F, color="black") + #colori[i]
      theme_bw(base_size=6) + labs(x=NULL, y=NULL) + theme(panel.grid=element_blank(), panel.background=element_blank(), axis.ticks=element_blank(), axis.text=element_text(), plot.background=element_rect(fill="transparent",color=NA)) +
      scale_x_continuous(expand=c(0, 0), breaks=bre) + scale_y_continuous(expand=c(0, .02), limits=c(0, 1))
    
    png(paste0("./Figure 1/", t, ".png"), width=6, height=4, res=1200, units="cm", bg="transparent")
    print(p)
    dev.off()
    message(t, " - ", colori[i])
  }
}

cat("\n--------\n\n")
