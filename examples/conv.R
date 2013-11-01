# with timeadjustment
setwd("D:\\mschauer\\AppData\\julia\\packages\\SDE")
fname = "reslong"
C = read.csv(paste(fname,".csv", sep=""))
#f = 1./1.96
f = 1.
svg(paste(fname,".svg", sep=""))
#C$N = log((1:length(C$N))*100, 2)
C$N = (1:length(C$N))*round(C$N[1])
C$t1 = C$N
C$t2 = C$N
plot(C$t1, C$px, t="l", xlim=c(0, max(tail(C$t1),tail(C$t2))), ylim = c(0.0,0.07),  col="#0000AA", xlab="N", ylab="p")
lines(C$t1, C$px-f*C$sex,lwd = 0.4,  col="#6666AA")
lines(C$t1, C$px+f*C$sex,lwd = 0.4,  col="#6666AA")

lines(C$t2, C$po, col="red")
lines(C$t2, C$po-f*C$seo,lwd = 0.4, col="#FF6666")
lines(C$t2, C$po+f*C$seo,lwd = 0.4, col="#FF6666")
dev.off()
stop("ok")

#without time adjustment
setwd("D:\\mschauer\\AppData\\julia\\packages\\SDE")
C = read.csv("res6.csv")
#C2 = read.csv("res.csv")
#f = 1./1.96
f = 1.
svg("output6.svg")
#C$N = log((1:length(C$N))*100, 2)
C$N = (1:length(C$N))*round(C$N[1])
plot(C$N, C$px, t="l", ylim = c(0.0,0.07),  col="#0000AA", xlab="N", ylab="p")
lines(C$N, C$px-f*C$sex,lwd = 0.4,  col="#6666AA")
lines(C$N, C$px+f*C$sex,lwd = 0.4,  col="#6666AA")

lines(C$N, C$po, col="red")
lines(C$N, C$po-f*C$seo,lwd = 0.4, col="#FF6666")
lines(C$N, C$po+f*C$seo,lwd = 0.4, col="#FF6666")
lines(C$N, rep(1.0995, length(C$N)), col="#555555")
dev.off()
#lines(C$N, (C2$px)[1:length(C$N)],  col="#00AA00")
#lines(C$N, (C2$px-f*C2$sex)[1:length(C$N)], lt="dotted", col="#00AA00")
#lines(C$N, (C2$px+f*C2$sex)[1:length(C$N)], lt="dotted", col="#00AA00")
#lines(C$N, (C2$po)[1:length(C$N)],  col="#00AA00")
#lines(C$N, (C2$po-f*C2$seo)[1:length(C$N)], lt="dotted", col="#00AA00")
#lines(C$N, (C2$po+f*C2$seo)[1:length(C$N)], lt="dotted", col="#00AA00")
#C = read.csv("res2.csv")
stop("ok")
C = read.csv("res2.csv")
C2 = read.csv("res.csv")
C$N = (1:length(C$N))*100
X11()
plot(C$N, NA*C$N, t="l", ylim = c(0.5,1.5))
lines(C$N, (C2$po)[1:length(C$N)],  col="#00AA00")
lines(C$N, (C2$po-f*C2$seo)[1:length(C$N)], lt="dotted", col="#00AA00")
lines(C$N, (C2$po+f*C2$seo)[1:length(C$N)], lt="dotted", col="#00AA00")

lines(C$N, C$po, col="red")
lines(C$N, C$po-f*C$seo, lt="dotted", col="#FF0000")
lines(C$N, C$po+f*C$seo, lt="dotted", col="#FF0000")
