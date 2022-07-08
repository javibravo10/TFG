
par(mar=c(5,4.5,4,2.5)+0.1) #margenes para que se visualice bien
par(mfrow=c(2,2)) # matriz de 2x2 gráficos (se rellenan por filas)

#diferentes perfiles de densidad
curve(1/(1+(x/2)^2),0,6, xlab=expression(r),ylab=expression(rho(r)) ,main='Perfil isotermo', ylim=c(0,1.1))
curve(1/(x*(1+x)^2),0,6 ,xlab=expression(r),ylab=expression(rho(r)) ,main='Perfil Navarro-Frenk-White')
curve(1/((1+x)*(1+x^2)) ,xlab=expression(r),ylab=expression(rho(r)) ,0,6,main='Perfil de Burket', ylim=c(0,1.1))
curve(exp(-(x)^0.5),0,6 ,xlab=expression(r),ylab=expression(rho(r)) ,main='Perfil Einasto', ylim=c(0,1.1))
  
#velocidad teorica y experimental inventada
par(mar=c(5,4.5,4,2.5)+0.1,mgp=c(1.5,0.5,0)) #margenes para que se visualice bien
curve(1*x,0,1, xlim=c(0,6),ylim=c(0,1.5),xaxt='n', yaxt='n', xlab=expression(r), ylab=expression(V(r)),lwd=2)
curve(1*x,0,1, col='blue',add=T,lwd=2)
curve(1/(sqrt(x)),1,6, add=T,col='blue',lwd=2)
curve(0.895+1/((x+0.5)*(x+1.5)^2),1,6,add=T,lwd=2)
segments(x0=1, x1=1, y0=0, y1=1, lwd=2, lty=2)
axis(1, at = c(1),labels = expression(R[bulbo]))
legend("topright",legend=c(expression(V[teo]),expression(V[exp])),col=c("blue","black"),lty=1,lwd=2)


#velocidad teorica y experimental buena
par(mar=c(5,4.5,4,2.5)+0.1) #margenes para que se visualice bien
curve(sqrt((log((x+1)^2*(x^2+1))- 2*atan(x))/x) ,xlab=expression(r),ylab=expression(v[halo](r)) ,0,10, ylim=c(0,1.1), xaxt='n',yaxt='n',lwd=2)
curve(0.8*x,0,1.161, col='blue',lwd=2,add=T)
curve(1/(sqrt(x)),1.155,10, add=T,col='blue',lwd=2)
segments(x0=1.161, x1=1.161, y0=0, y1=0.95, lwd=2, lty=2)
axis(1, at = c(1.161),labels = expression(R[bulbo]))
legend("topright",legend=c(expression(V[teo]),expression(V[exp])),col=c("blue","black"),lty=1,lwd=2)



#velocidades de los distintos halos

par(mfrow=c(2,2)) # matriz de 2x2 gráficos (se rellenan por filas)
curve(sqrt(1-atan(x)/x), 0,10, xlab=expression(r),ylab=expression(v[halo](r)) ,main='Perfil isotermo', ylim=c(0,1.1),xaxt='n', yaxt='n')
curve(sqrt((log(x+1)+1/(x+1)-1)/x),0,10 ,xlab=expression(r),ylab=expression(v[halo](r)) ,main='Perfil Navarro-Frenk-White',xaxt='n', yaxt='n')
curve(sqrt((log((x+1)^2*(x^2+1))- 2*atan(x))/x) ,xlab=expression(r),ylab=expression(v[halo](r)) ,0,10,main='Perfil de Burket', ylim=c(0,1.1),xaxt='n', yaxt='n')



#efecto del parametro s en las densidades
curve(1/(1+(x/2)^2),0,6, xlab=expression(r),ylab=expression(rho(r)) ,main=expression(s==1/2), ylim=c(0,1.1))
curve(1/(1+(x)^2),0,6, xlab=expression(r),ylab=expression(rho(r)) ,main=expression(s==1), ylim=c(0,1.1))
curve(1/(1+(x*2)^2),0,6, xlab=expression(r),ylab=expression(rho(r)) ,main=expression(s==2), ylim=c(0,1.1))
curve(1/(1+(x*4)^2),0,6, xlab=expression(r),ylab=expression(rho(r)) ,main=expression(s==4), ylim=c(0,1.1))

#efecto del parametro s en las velocidades

par(mfrow=c(2,2)) # matriz de 2x2 gráficos (se rellenan por filas)
curve(sqrt(((log(x/2+1)+1/(x/2+1)-1)/x)/(1/2)^3),0,10 ,xlab=expression(r),ylab=expression(v[halo](r)) ,main=expression(s==1/2),xaxt='n', yaxt='n')
curve(sqrt((log(x+1)+1/(x+1)-1)/x),0,10 ,xlab=expression(r),ylab=expression(v[halo](r)) ,main=expression(s==1),xaxt='n', yaxt='n')
curve(sqrt(((log(2*x+1)+1/(2*x+1)-1)/x)/2^3),0,10 ,xlab=expression(r),ylab=expression(v[halo](r)) ,main=expression(s==2),xaxt='n', yaxt='n')
curve(sqrt(((log(5*x+1)+1/(5*x+1)-1)/x)/5^3),-0.175,10 ,xlab=expression(r),ylab=expression(v[halo](r)) ,main=expression(s==5),xaxt='n', yaxt='n')
