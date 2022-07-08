

################################################################################
###################### LECTURA DE DATOS EN FICHERO #############################
################################################################################

#En este caso tenemos kpc para los radios y las velocidades km/s
#buscamos obtener los resultados trabajando en pc, para obtener 
#a=rho_0 en Msolares/pc^3 y b=s en 1/pc

#cte de gravitacion en las unidades que lo necesitamos km^2*pc / Msol*s^2
#esto es para obtener s en kpc^-1 y pho en Msol/pc^3

#La otra opcion es G=43*10^-4 y los radios pasarlos a pc 
# datos$R<-datos$R*1000, entonces las sale en pc^-1 y sus unidades son 3 0 mas

G<-430


#Leemos los datos tal y como estan en el fichero,  los metemos en un dataset y
#damos nombre a sus columnas 
datos<-read.table("data_format_DDO47.txt", header=F)
colnames(datos)<-c("R","Vgas", "Vstar", "Vbulb", "Vtot", "Verr")


#Definimos la masa debida a las componentes visibles
Vmv<-sqrt(datos$Vgas^2+datos$Vstar^2)



################################################################################



################################################################################
#################### REPRESENTACIÓN GRÁFICA DE CHI^2 ###########################
################################################################################

# A partir de ahora los parámetros a y b hacen referencia a rho_0 y s cada uno

#Definimos la función CHI^2 para el ISO reducida (dividida por el numero de 
#grados de libertad (6) en este caso)
sum<-0 #sumatoria necesaria

chi<-function(b,a){
  for(j in 1:8) 
    
    sum<-sum + ((datos$Vtot[j]-(sqrt(Vmv[j]^2+pi*G*(a/b^2)*((1/(b*datos$R[j]))*log(((b*datos$R[j])^2+1)*((b*datos$R[j])+1)^2)-(2/(b*datos$R[j]))*atan(b*datos$R[j])))))^2)/(datos$Verr[j]^2*6)
  
  return (sum)
}


#Evaluación de la función donde sabemos que alcanza realmente el mínimo
chi(0.2696720327,0.449345609)


#Asignamos valores a los parámetros en torno al mínimo para que sean los que se
#representen
b<-seq(0.10,0.40,length=700) #parámetro s
a<-seq(0.30,0.60,length=700) #parámetro rho_0

z<-outer(b, a, chi)

#Representacion de las curvas de nivel
par(mar=c(5,4.5,4,2.5)+0.1) #margenes para que se visualice bien
image(b,a,z,main="DDO47 con BURKET",xlab=expression(s ~~(kpc^-1)),ylab=expression(rho[0]~~ (x10^-1)~~(M[sol]/pc^3)))
contour(b,a,z,add=T, levels=c(0.15,0.151,0.152,0.155,0.16,0.18))
contour(b,a,z,add=T, levels=c(0.21,0.25,0.3))
contour(b,a,z,add=T, levels=seq(0.35, 0.65, length.out=4))

contour(b,a,z,add=T,levels=c(0.7,0.4,0.2,0.16,0.15))
contour(b,a,z,add=T)

################################################################################



################################################################################
########################### MÉTODOS ITERATIVOS #################################
################################################################################

## Comun a todos los metodos de busqueda de minimos


#Necesitamos la libreria de caluclo numerico de derivadas
library(numDeriv)  

#Parámetros iniciales del método
theta_0 <- c(0.10, 0.55) #valores iniciales de los parámetros s y rho_0

n <- 60 #numero maximo de iteraciones
epsilon<-0.00001 #criterio de parada

#Hay que escribir la función de forma que los parámetros se recojan en un vector

sum<-0

# A lo que le calcularemos el hessiano y el gradiente, función escalar
chi2<-function(w){
  for(j in 1:8) 
    
    sum<-sum + ((datos$Vtot[j]-(sqrt(Vmv[j]^2+pi*G*(w[2]/w[1]^2)*((1/(w[1]*datos$R[j]))*log(((w[1]*datos$R[j])^2+1)*((w[1]*datos$R[j])+1)^2)-(2/(w[1]*datos$R[j]))*atan(w[1]*datos$R[j])))))^2)/(datos$Verr[j]^2*6)
  
  return (sum)
}

#A lo que le calcularemos el jacobiano, función vectorial 

aux<-vector("numeric",8)
f<-function(w){
  for(j in 1:8) 
    
    aux[j]<- (sqrt(Vmv[j]^2+pi*G*(w[2]/w[1]^2)*((1/(w[1]*datos$R[j]))*log(((w[1]*datos$R[j])^2+1)*((w[1]*datos$R[j])+1)^2)-(2/(w[1]*datos$R[j]))*atan(w[1]*datos$R[j]))))
  
  return (aux)
}

jacobian(f,theta_0)
t(jacobian(f,theta_0))

  
#Matriz de errores

W<-diag(8)
for(j in 1:8)
  W[j,j]<-1/(datos$Verr[j]^2*6)
  

#Vector de residuos
res<-vector("numeric",8)
r<-function(w){
  for(j in 1:8) 
    
    res[j]<- (datos$Vtot[j]-(sqrt(Vmv[j]^2+pi*G*(w[2]/w[1]^2)*((1/(w[1]*datos$R[j]))*log(((w[1]*datos$R[j])^2+1)*((w[1]*datos$R[j])+1)^2)-(2/(w[1]*datos$R[j]))*atan(w[1]*datos$R[j])))))
  
  return (res)
}

#as.matrix(r(theta_0)) convierte a matriz (8,1) el vector de residuos

#Se verifica como es de esperar que 
# -2*t(jacobian(f,theta_0))%*%W%*%as.matrix(r(theta_0))=grad(chi2,theta_0)

#Se verifica tambien que 2*t(jacobian(f,theta_0))%*%W%*%jacobian(f,theta_0) 
# es aprox hessian(chi2,theta_0)


################################################################################


################################################################################
############################# AJUSTE FINAL #####################################
################################################################################


#Se guardan los valores de a y b obtenidos en el algoritmo de Levenberg Marquardt

#queremos representar la componente de la velocidad debida al halo
velocidad<-sqrt(datos$Vtot^2-Vmv^2)  
velocidad

#calculamos los errores por propagacion cuadrática
delta<-(1/velocidad)*datos$Vtot*datos$Verr


ajuste<-function(r){
  sqrt(pi*G*(a/b^2)*((1/(b*r)) 
                     *((log(((b*r)^2+1)*(b*r+1)^2))-2*atan(r*b))))
}

par(mar=c(5,4.5,4,2.5)+0.1) #margenes para que se visualice bien
curve(ajuste,0,3.2,ylim=c(0,80),xlab=expression(r~~(kpc)),ylab=expression(V[halo](r)~~(km/s^2)),col="purple",lwd=2)
points(datos$R,velocidad,pch=21,col="blue",bg="blue")
arrows(x0=datos$R,y0=velocidad-delta, x1=datos$R, y1=velocidad+delta,length = 0.15, code = 3, angle = 90,col="blue",lwd=1)
legend("topleft",legend=c(expression(V[exp]),expression(V[teo])),col=c("blue","purple"),lty=1,lwd=3)

