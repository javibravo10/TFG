
#PERFIL BURKET
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

G<-4300000


#Leemos los datos tal y como estan en el fichero,  los metemos en un dataset y
#damos nombre a sus columnas 
datos<-read.table("data_format_DDO185.txt", header=F)
colnames(datos)<-c("R","Vgas", "Vstar", "Vbulb", "Vtot", "Verr")


#Definimos la masa debida a las componentes visibles
Vmv<-sqrt(datos$Vgas^2+datos$Vstar^2)
datos$R<-datos$R*0.01 #convertimos siendo conscientes de que queremos visualizar bien los resultados



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
chi(0.2001052333,0.17498574)


#Asignamos valores a los parámetros en torno al mínimo para que sean los que se
#representen
b<-seq(0.15,0.25,length=700) #parámetro s
a<-seq(0.15,0.20,length=700) #parámetro rho_0

z<-outer(b, a, chi)

#Representacion de las curvas de nivel
par(mar=c(5,4.5,4,2.5)+0.1) #margenes para que se visualice bien
image(b,a,z,main="DDO185 con BURKET",xlab=expression(s ~~(kpc^-1)),ylab=expression(rho[0]~~ (x10^-1)~~(M[sol]/pc^3)))

contour(b,a,z,add=T)
################################################################################

