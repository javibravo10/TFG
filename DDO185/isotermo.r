
#PERFIL ISOTERMO
#GALAXIA DDO185


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

G<-43000


#Leemos los datos tal y como estan en el fichero,  los metemos en un dataset y
#damos nombre a sus columnas 
datos<-read.table("data_format_DDO185.txt", header=F)
colnames(datos)<-c("R","Vgas", "Vstar", "Vbulb", "Vtot", "Verr")


#Definimos la masa debida a las componentes visibles
Vmv<-sqrt(datos$Vgas^2+datos$Vstar^2)
datos$R<-datos$R*0.1 #convertimos los datos a pc


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
chi(0.1995774,0.17465)


#Asignamos valores a los parámetros en torno al mínimo para que sean los que se
#representen



b<-seq(0.01,0.38,length=700) #parámetro s
a<-seq(0.160,0.185,length=700) #parámetro rho_0


z<-outer(b, a, chi)

#Representacion de las curvas de nivel
par(mar=c(5,4.5,4,2.5)+0.1) #margenes para que se visualice bien
image(b,a,z,main="DDO185 con isotermo",xlab=expression(s ~~ (x10^-1)~~(kpc^-1)),ylab=expression(rho[0]~~ (x10^-1)~~(M[sol]/pc^3)))
contour(b,a,z,add=T, levels=c(0.988,0.98651))

contour(b,a,z,add=T, levels=c(0.99, 0.993, 0.996, 1, 1.01, 1.02, 1.03, 1.04))

contour(b,a,z,add=T)
################################################################################




################################################################################
####################### MALLADO PARA EL PUNTO INICIAL ##########################
################################################################################



b<-seq(0.1,5,length=700) #parámetro s
a<-seq(0.1,1,length=700) #parámetro rho_0


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


bmin<-0.1+ (resto-1)*(5-0.1)/700
amin<-0.1+ (cociente-1)*(1-0.1)/700

bmin
amin

theta_0<-c(bmin,amin)

#buena semilla inicial bmin,amin


################################################################################




################################################################################
########################### MÉTODOS ITERATIVOS #################################
################################################################################

## Comun a todos los metodos de busqueda de minimos


#Necesitamos la libreria de caluclo numerico de derivadas
library(numDeriv)  

#Parámetros iniciales del método
theta_0 <- c(0.1, 0.17) #valores iniciales de los parámetros s y rho_0

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

#as.matrix(r(theta_0)) convierte a matriz (8,1) el vector de residuos

#Se verifica como es de esperar que 
# -2*t(jacobian(f,theta_0))%*%W%*%as.matrix(r(theta_0))=grad(chi2,theta_0)

#Se verifica tambien que 2*t(jacobian(f,theta_0))%*%W%*%jacobian(f,theta_0) 
# es aprox hessian(chi2,theta_0)


################################################################################
