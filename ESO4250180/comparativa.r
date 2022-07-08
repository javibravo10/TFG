##################################################################
########### Comparativa de los ajustes obtenidos ESO
##################################################################

datos<-read.table("data_format_ESO4250180.txt", header=F)

colnames(datos)<-c("R","Vgas", "Vstar", "Vbulb", "Vtot", "Verr")

G<-43*10^(-4)


#queremos representar la componente de la velocidad debida al halo
Vmv<-sqrt(datos$Vgas^2+datos$Vstar^2)
velocidad<-sqrt(datos$Vtot^2-Vmv^2)  
velocidad

#calculamos los errores por propagacion cuadrática
delta<-(1/velocidad)*datos$Vtot*datos$Verr


#Iso 

a<-0.0302062496
b<-0.0002273927

ajuste<-function(r){
  sqrt(4*pi*G*(a/b^2)*(1-((1/(b*r*1000))*atan(b*r*1000))))
}

par(mar=c(5,4.5,4,2.5)+0.1) #margenes para que se visualice bien
curve(ajuste,col="red", 0.01,14.5,ylim=c(0,150),xlab=expression(r~~(kpc)),ylab=expression(V[halo](r)~~(km/s^2)),lwd=2, main='ESO4250180')
points(datos$R,velocidad,pch=21,col="blue",bg="blue")
arrows(x0=datos$R,y0=velocidad-delta, x1=datos$R, y1=velocidad+delta,length = 0.15, code = 3, angle = 90,col="blue",lwd=1)


#Burket

a<-0.0307272891

b<-0.0001299592

ajuste<-function(r){
  sqrt(pi*G*(a/b^2)*((1/(b*r*1000)) 
                     *((log(((b*r*1000)^2+1)*(b*r*1000+1)^2))-2*atan(r*b*1000))))
}

curve(ajuste,col="green",lwd=2,add=T)



#NFW
a<-0.000524764333
b<-8.37*10^(-6)

ajuste<-function(r){
  sqrt(4*pi*G*(a/b^2)*(1/(b*r*1000))*((1/(b*r*1000+1))-1+log(b*r*1000+1)))
}

curve(ajuste, add=T, col="purple",lwd=2)




legend("topleft",legend=c(expression(V[exp]),expression(V[teo]~isotermo), expression(V[teo]~Burket), expression(V[teo]~NFW)),col=c("blue","red","green","purple"),lty=1,lwd=3)
