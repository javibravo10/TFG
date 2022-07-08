

#a y b hacen referencia a rho_0 y a s 
chi<-function(b,a) 2+a^2+b^2

b<-seq(-2,2,length=500)
a<-seq(-2,2,length=500)


z<-outer(b, a, chi)


par(mar=c(5,4.5,4,2.5)+0.1) #margenes para que se visualice bien

#grafica en 3d de la funcion
persp(b,a,z, xlab='', ylab='', zlab='',axes=F)

#representacion de unas curvas de nivel

image(b,a,z,xlab=expression(a[1]),ylab=expression(a[2]))
contour(b,a,z,add=T)
contour(b,a,z,add=T, levels=c(2.5))

points(-1,0,type='b',pch=21, col='blue', bg='cyan', cex=1 )
text(locator(1), labels=expression(a), col='blue')
arrows(-0.95,0, 1.4, 0, col='blue', lwd=2)
points(1.42,0,type='b',pch=21, col='blue', bg='cyan', cex=1 )
text(locator(1), labels=expression(a+alpha*h), col='blue')

library(latex2exp)
