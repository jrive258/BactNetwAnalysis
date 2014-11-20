#Plotting example
require(plotrix)
par(mar=c(1,1,1,1)) #Better margins
plot(c(0,1),c(0,1),col="white",axes=FALSE,ann=FALSE)
lines(c(0.3,0.7), c(0.4,0.65), type="l", col="green")
draw.circle(0.3, 0.4, 0.05, col="blue", border="red") 
draw.circle(0.7, 0.65, 0.05, col="blue", border="red") 
text(0.3,0.4, labels = "C.Phyl1")
text(0.7,0.65, labels = "F.Colo1")