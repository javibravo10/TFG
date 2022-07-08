

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
datos<-read.table("data_format_DDO47.txt", header=F)
colnames(datos)<-c("R","Vgas", "Vstar", "Vbulb", "Vtot", "Verr")
datos$R<-datos$R*10^(-3) #convertimos siendo conscientes de que queremos visualizar bien los resultados

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
    
    sum<-sum + ((datos$Vtot[j]-(sqrt(Vmv[j]^2+4*pi*G*(a/b^2)*(1/(b*datos$R[j]))*((1/(b*datos$R[j]+1))-1+log(b*datos$R[j]+1)))))^2)/(datos$Verr[j]^2*6)
  
  return (sum)
}


#Evaluación de la función donde sabemos que alcanza realmente el mínimo
chi(0.2559391564,0.968174821)


#Asignamos valores a los parámetros en torno al mínimo para que sean los que se
#representen
b<-seq(0.20,0.30,length=700) #parámetro s
a<-seq(0.9,0.10,length=700) #parámetro rho_0

z<-outer(b, a, chi)

#Representacion de las curvas de nivel
par(mar=c(5,4.5,4,2.5)+0.1) #margenes para que se visualice bien
image(b,a,z,main="DDO47 con NFW",xlab=expression(s ~~ (x10^-6)~~(kpc^-1)),ylab=expression(rho[0]~~ (x10^-5)~~(M[sol]/pc^3)))
contour(b,a,z,add=T, levels=c(0.337,0.35,0.37,0.41))
contour(b,a,z,add=T,levels=c(0.5,0.6,0.8,0.9,1))
################################################################################


################################################################################
########################### MÉTODOS ITERATIVOS #################################
################################################################################

## Comun a todos los metodos de busqueda de minimos


#Necesitamos la libreria de caluclo numerico de derivadas
library(numDeriv)  

#Parámetros iniciales del método
theta_0 <- c(0.2, 0.9) #valores iniciales de los parámetros s y rho_0

n <- 60 #numero maximo de iteraciones
epsilon<-0.00001 #criterio de parada

#Hay que escribir la función de forma que los parámetros se recojan en un vector

sum<-0

# A lo que le calcularemos el hessiano y el gradiente, función escalar
chi2<-function(w){
  for(j in 1:8) 
    
    sum<-sum + ((datos$Vtot[j]-(sqrt(Vmv[j]^2+4*pi*G*(w[2]/w[1]^2)*(1/(w[1]*datos$R[j]))*((1/(w[1]*datos$R[j]+1))-1+log(w[1]*datos$R[j]+1)))))^2)/(datos$Verr[j]^2*6)
  
  return (sum)
}

#A lo que le calcularemos el jacobiano, función vectorial 

aux<-vector("numeric",8)
f<-function(w){
  for(j in 1:8) 
    
    aux[j]<- (sqrt(Vmv[j]^2+4*pi*G*(w[2]/w[1]^2)*(1/(w[1]*datos$R[j]))*((1/(w[1]*datos$R[j]+1))-1+log(w[1]*datos$R[j]+1))))
  
  return (aux)
}
  
#Matriz de errores

W<-diag(8)
for(j in 1:8)
  W[j,j]<-1/(datos$Verr[j]^2*6)
  

#Vector de residuos
res<-vector("numeric",8)
r<-function(w){
  for(j in 1:8) 
    
    res[j]<- (datos$Vtot[j]-(sqrt(Vmv[j]^2+4*pi*G*(w[2]/w[1]^2)*(1/(w[1]*datos$R[j]))*((1/(w[1]*datos$R[j]+1))-1+log(w[1]*datos$R[j]+1)))))
  
  return (res)
}


################################################################################

