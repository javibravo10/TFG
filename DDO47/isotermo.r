
#PERFIL ISOTERMO
#GALAXIA DDO47


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
    
    sum<-sum + ((datos$Vtot[j]-(sqrt(Vmv[j]^2+4*pi*G*(a/b^2)*(1-((1/(b*datos$R[j]))*atan(b*datos$R[j]))))))^2)/(datos$Verr[j]^2*6)
  
  return (sum)
}


#Evaluación de la función donde sabemos que alcanza realmente el mínimo
chi(0.4483903,0.40673)


#Asignamos valores a los parámetros en torno al mínimo para que sean los que se
#representen
b<-seq(0.38,0.55,length=700) #parámetro s
a<-seq(0.36,0.47,length=700) #parámetro rho_0

z<-outer(b, a, chi)

#Representacion de las curvas de nivel
par(mar=c(5,4.5,4,2.5)+0.1) #margenes para que se visualice bien
image(b,a,z,main="DDO47 con isotermo",xlab=expression(s ~~(kpc^-1)),ylab=expression(rho[0]~~ (x10^-1)~~(M[sol]/pc^3)))
contour(b,a,z,add=T, levels=c(0.1605, 0.161, 0.163, 0.165, 0.172))
contour(b,a,z,add=T, levels=c(0.20, 0.25, 0.325, 0.4, 0.5))
contour(b,a,z,add=T, levels=seq(0.60, 1, length.out=5))

################################################################################




################################################################################
####################### MALLADO PARA EL PUNTO INICIAL ##########################
################################################################################

################################################################################


b<-seq(0.00001,0.05,length=700) #parámetro s
a<-seq(0.0001,0.02,length=700) #parámetro rho_0


z<-outer(b, a, chi)
str(z)

min(z)
posicion<-which.min(z)
posicion

#Calculo de a y b que dan en el minimo en el mallado

cociente<-posicion%/%700 #cociente
resto<-posicion%%700 #resto

if(resto>0)  cociente<-cociente+1

z[resto,cociente]


bmin<-0.00001+ (resto-1)*(0.05-0.00001)/700
amin<-0.0001+ (cociente-1)*(0.02-0.0001)/700

bmin
amin

theta_0<-c(bmin,amin)

#buena semilla inicial bmin,amin



################################################################################
########################### MÉTODOS ITERATIVOS #################################
################################################################################

## Comun a todos los metodos de busqueda de minimos


#Necesitamos la libreria de caluclo numerico de derivadas
library(numDeriv)  

#Parámetros iniciales del método
theta_0 <- c(0.54, 0.37) #valores iniciales de los parámetros s y rho_0

n <- 100 #numero maximo de iteraciones
epsilon<-0.00001 #criterio de parada

#Hay que escribir la función de forma que los parámetros se recojan en un vector

sum<-0

# A lo que le calcularemos el hessiano y el gradiente, función escalar
chi2<-function(w){
  for(j in 1:8) 
    
    sum<-sum + ((datos$Vtot[j]-(sqrt(Vmv[j]^2+4*pi*G*w[2]*(1/w[1])^2*(1-((1/(w[1]*datos$R[j]))*atan(w[1]*datos$R[j]))))))^2)/(datos$Verr[j]^2*6)
  
  return (sum)
}

#A lo que le calcularemos el jacobiano, función vectorial 

aux<-vector("numeric",8)
f<-function(w){
  for(j in 1:8) 
    
    aux[j]<-(sqrt(Vmv[j]^2+4*pi*G*w[2]*(1/w[1])^2*(1-((1/(w[1]*datos$R[j]))*atan(w[1]*datos$R[j])))))
  
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
    
    res[j]<-datos$Vtot[j]-(sqrt(Vmv[j]^2+4*pi*G*w[2]*(1/w[1])^2*(1-((1/(w[1]*datos$R[j]))*atan(w[1]*datos$R[j])))))
  
  return (res)
}


################################################################################

################################################################################
#################### DESCENSO DEL GRADIENTE LINE SEARCH ########################
################################################################################

gamma <- 1 #tasa de aprendizaje (si toco un 0 mas se me dispara)
#(realmente es gamma y no alpha ya que se normaliza la direccion del gradiente)



#Definimos una matriz donde se metan en dos columnas los resultados de los parametros
# tras cada iteracion, 

thetas <- matrix(NA, ncol=2, nrow=n+1)    
colnames(thetas) <- c('s', 'rho_0')
norma<-vector("numeric", n)
gradiente<-matrix(NA, ncol=2, nrow=n+1)#vector donde almacenaremos el gradiente en cada paso
thetas[1, ] <- theta_0  #le metemos a la matriz los parametros inciales
veces<-vector("numeric", n) #para garantizar salir del if

#Para la busqueda de linea

c<-0.5 #parametro de la busqueda de linea
#se toma un alpha grande inicialmente y luego se va disminuyendo hasta que se cumpla la condicion.

alpha<-0.0001 #paso 

aux<- matrix(NA, ncol=2, nrow=n+1)  #vector auxiliar para antes de dar los pasos
gradaux<- matrix(NA, ncol=2, nrow=n+1)  #vector auxiliar para antes de dar los pasos

tau<- 0.1 #parámetro de reduccion del alpha en caso de no encontrarse un descenso


#Algoritmo del descenso
for (i in 2:(n+1)){
  alpha<-0.01 #paso 
  gradiente[i-1,]<-grad(func=chi2, x=thetas[i-1, ])
  norma[i-1]<-sqrt(gradiente[i-1,1]^2+gradiente[i-1,2]^2)
  aux[i,] <- thetas[i-1,] - (alpha) * gradiente[i-1,]
  gradaux[i,]<-grad(chi2, aux[i,])
  
  
  #Para que en ningun caso nos pasemos al descender en esa dirección 
  if(chi2(aux[i,])>chi2(thetas[i-1,])-c*alpha*norma[i-1]^2 && gradiente[i-1,]*gradaux[i,]>c*norma[i-1]^2&&veces[i]<10){
    
    alpha<-alpha*tau
    aux[i,] <- thetas[i-1,] - (alpha) * gradiente[i-1,]
    gradaux[i,]<-grad(chi2, aux[i,])
    veces[i]<-veces[i]+1
    
  }
    
  else 
  {
    thetas[i,] <- aux[i,]
  }
  
  if(veces[i]>=9) {
    break
    cat('No se consigue una dirección de descenso')

  }
  
  
  #criterio de parada
  if((abs(sqrt(thetas[i,1]^2+thetas[i,2]^2)-sqrt(thetas[i-1,1]^2+thetas[i-1,2]^2)))<epsilon){
    break
  }
}

thetas[1:i,]


points(thetas[,1],thetas[,2],type='b',pch=21, col='blue', bg='cyan', cex=1 )



################################################################################


################################################################################
########################## MÉTODO DE GAUSS-NEWTON ##############################
################################################################################

deltas <- matrix(NA, ncol=2, nrow=n+1)    
colnames(deltas) <- c('s', 'rho_0')
deltas[1, ] <- theta_0  #le metemos a la matriz los parametros inciales


for (i in 2:(n+1)){
  deltas[i,] <- deltas[i-1,]-as.vector(solve(2*t(jacobian(f,deltas[i-1,]))%*%W%*%jacobian(f,deltas[i-1,]))%*%grad(func=chi2, x=deltas[i-1, ]))
  
  
  #criterio de parada
  if((abs(sqrt(deltas[i-1,1]^2+deltas[i-1,2]^2)-sqrt(deltas[i,1]^2+deltas[i,2]^2)))<epsilon){
    break
  }
}

deltas[1:i,]


points(deltas[,1],deltas[,2],type='b',pch=21, col='springgreen', bg='springgreen', cex=1 )


################################################################################


################################################################################
########################### LEVENBERG MARQUARDT  ###############################
################################################################################

gammas <- matrix(NA, ncol=2, nrow=n+1)    #matriz donde se almacenan las iteraciones
colnames(gammas) <- c('s', 'rho_0')
gammas[1, ] <- theta_0  #le metemos a la matriz los parametros inciales

tau<-0.001
mu<-tau*max(diag(2*t(jacobian(f,gammas[1,]))%*%W%*%jacobian(f,gammas[1,])))#estimacion incial del paso
rho<-0 #gain ratio
nu<-2 #parametro auiliar para recalucular mu

gradiente<-matrix(NA, ncol=2, nrow=n+1)#vector donde almacenaremos el gradiente en cada paso
previo<-matrix(NA, ncol=2, nrow=n+1) #vector donde almacenaremos la iteracion previa a guardarla


for (i in 2:(n+1)){
  gradiente[i-1,]<-grad(func=chi2, x=gammas[i-1, ])
  previo[i,] <- gammas[i-1,]-as.vector(solve(2*t(jacobian(f,gammas[i-1,]))%*%W%*%jacobian(f,gammas[i-1,])+mu*diag(2))%*%gradiente[i-1,])
  rho<-2*(chi2(gammas[i-1,])-chi2(previo[i,]))/((t(previo[i,]-gammas[i-1,])%*%(mu*((previo[i,]-gammas[i-1,]))-gradiente[i-1,])))
  
  
  #criterio de parada
  
    if((abs(sqrt(gammas[i-1,1]^2+gammas[i-1,2]^2)-sqrt(previo[i,1]^2+previo[i,2]^2)))<epsilon){
    
    b<-previo[i,1]
    a<previo[i,2]
    break
  }
  
  if(rho>0) {
    mu<-mu*max(1/3, 1-(2*rho-1)^3)
    nu<-2
    gammas[i,]<-previo[i,]
  }
  else {
    mu<-mu*nu
    nu<-2*nu
    gammas[i,]<-gammas[i-1,] #en este caso no es una buena aproximacion luego no nos interesan los valores obtenidos
  }
  
  
}

gammas[1:i-1,]


points(gammas[,1],gammas[,2],type='b',pch=21, col='maroon1', bg='maroon1', cex=1 )
################################################################################


legend("topleft",legend=c('Descenso más rápido', 'Gauss-Newton','Levenberg-Marquardt'),col=c("cyan","springgreen","maroon1"),pch=c(19,19,19), lty=1,lwd=2)


################################################################################
############################# AJUSTE FINAL #####################################
################################################################################

#sacas los parametros del ajuste

  
#queremos representar la componente de la velocidad debida al halo
velocidad<-sqrt(datos$Vtot^2-Vmv^2)  
velocidad

#calculamos los errores por propagacion cuadrática
delta<-(1/velocidad)*datos$Vtot*datos$Verr


ajuste<-function(r){
  sqrt(4*pi*G*(a/b^2)*(1-((1/(b*r))*atan(b*r))))
}

par(mar=c(5,4.5,4,2.5)+0.1) #margenes para que se visualice bien
curve(ajuste,0,3.2,ylim=c(0,80),xlab=expression(r~~(kpc)),ylab=expression(V[halo](r)~~(km/s^2)),col="purple",lwd=2, main='DDO47 con Burket')
points(datos$R,velocidad,pch=21,col="blue",bg="blue")
arrows(x0=datos$R,y0=velocidad-delta, x1=datos$R, y1=velocidad+delta,length = 0.15, code = 3, angle = 90,col="blue",lwd=1)
legend("topleft",legend=c(expression(V[exp]),expression(V[teo])),col=c("blue","purple"),lty=1,lwd=3)

################################################################################
