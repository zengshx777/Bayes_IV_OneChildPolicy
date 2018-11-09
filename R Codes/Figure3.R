#Codes Below Produce Figure 3 in the Paper
load("Figure.RData")
#load("Figure3.RData")
#pdf("Firstage_Prob_Quantile.pdf",height=4,width= 8)
par(mfrow=c(1, 1),mar = c(2,4.5,0,0))
#Index to Draw Points Along the Line
index=seq(1,1000,by=100)
#Create a Grid of Quantile for the IV
z.quantile=seq(0,1,length=1000)
plot(z.quantile,newz.p[,1],ylim=range(newz.p),type='l',lty=1,
     lwd=2,xlab="Quantile of the IFPPR",
     ylab="Probability of Being the Only Child",
     main="",
     yaxt = "n", xaxt = "n", axes = F)
axis(1, at = seq(0, 1, 0.1), 
     labels = seq(0, 1, 0.1))   
axis(2, at = c(0, 0.1, 0.3, 0.5, 0.7), 
     labels = c("0", ".1", ".3", ".5", ".7")) 

points(z.quantile[index],newz.p[index,1],pch=1,cex=1)
lines(z.quantile,newz.p[,2],type='l',lty=1,lwd=2)
points(z.quantile[index],newz.p[index,2],pch=2,cex=1)
lines(z.quantile,newz.p[,3],type='l',lty=1,lwd=2)
points(z.quantile[index],newz.p[index,3],pch=5,cex=1)
lines(z.quantile,newz.p[,4],type='l',lty=1,lwd=2)
points(z.quantile[index],newz.p[index,4],pch=6,cex=1)
legend("topleft",legend=Subgroup.name,lty=1,pch=c(1,2,5,6),cex=0.7)

#dev.off()
