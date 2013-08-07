C = read.csv("C:\\Users\\mschauer\\Dropbox\\res2.csv")
C2 = read.csv("C:\\Users\\mschauer\\Dropbox\\res.csv")
#f = 1./1.96
f = 1.
#svg("C:\\Users\\mschauer\\Dropbox\\output2.svg")
#C$N = log((1:length(C$N))*100, 2)
x11()
C$N = (1:length(C$N))*100
plot(C$N, C$px, t="l", ylim = c(0.5,1.5))
lines(C$N, C$px-f*C$sex, lt="dotted", col="gray")
lines(C$N, C$px+f*C$sex, lt="dotted", col="gray")

lines(C$N, C$po, col="red")
lines(C$N, C$po-f*C$seo, lt="dotted", col="#FF0000")
lines(C$N, C$po+f*C$seo, lt="dotted", col="#FF0000")
#dev.off()
#lines(C$N, (C2$px)[1:length(C$N)],  col="#00AA00")
#lines(C$N, (C2$px-f*C2$sex)[1:length(C$N)], lt="dotted", col="#00AA00")
#lines(C$N, (C2$px+f*C2$sex)[1:length(C$N)], lt="dotted", col="#00AA00")
#lines(C$N, (C2$po)[1:length(C$N)],  col="#00AA00")
#lines(C$N, (C2$po-f*C2$seo)[1:length(C$N)], lt="dotted", col="#00AA00")
#lines(C$N, (C2$po+f*C2$seo)[1:length(C$N)], lt="dotted", col="#00AA00")C = read.csv("C:\\Users\\mschauer\\Dropbox\\res2.csv")
C = read.csv("C:\\Users\\mschauer\\Dropbox\\res2.csv")
C2 = read.csv("C:\\Users\\mschauer\\Dropbox\\res.csv")
C$N = (1:length(C$N))*100
X11()
plot(C$N, NA*C$N, t="l", ylim = c(0.5,1.5))
lines(C$N, (C2$po)[1:length(C$N)],  col="#00AA00")
lines(C$N, (C2$po-f*C2$seo)[1:length(C$N)], lt="dotted", col="#00AA00")
lines(C$N, (C2$po+f*C2$seo)[1:length(C$N)], lt="dotted", col="#00AA00")

lines(C$N, C$po, col="red")
lines(C$N, C$po-f*C$seo, lt="dotted", col="#FF0000")
lines(C$N, C$po+f*C$seo, lt="dotted", col="#FF0000")