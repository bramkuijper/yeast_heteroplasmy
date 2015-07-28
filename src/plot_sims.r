#!/usr/bin/Rscript --vanilla

arg.commands <- commandArgs(trailingOnly=T)

file <- arg.commands[[1]]

data <- read.table(file, header=T, sep=";", nrow=20800)

library("lattice")
library("colorRamps")
pdf(file=paste("graph_",basename(file),".pdf",sep=""))
print(levelplot(sqrt(freq) ~ generation * p_i | spore,
                data=data,
                xlab="generations, t",
                ylab=expression(paste("within-ind freq. ",italic(C)[1])),
                zlab=expression(paste("freq. ",italic(C)[1])),
                strip=function(...,strip.levels) { strip.default(...,strip.levels=T) },
                col.regions=matlab.like(200)))
dev.off()

pdf(file=paste("graph_hist_",basename(file),".pdf",sep=""))

datx <- data[data$generation == 25,c("freq","p_i","spore")]
names(datx) <- c("freq_gen_25",names(datx)[2:length(names(datx))])

datx$freq_gen_5 <- data[data$generation == 5,"freq"]

print(xyplot(freq_gen_25 + freq_gen_5 ~ p_i | spore,
                lwd=c(1.0,0.5),
                col=c("black","red"),
                type="l",
                key=list(
                            text=list(c("generation 25","generation 5")),
                            lines=list(col=c("black","red"))
                ),
                data=datx,
                xlab=expression(paste("within-ind freq. ",italic(C)[1])),
                ylab="population freq.",
                strip=function(...,strip.levels) { strip.default(...,strip.levels=T) }
                ))
dev.off()
