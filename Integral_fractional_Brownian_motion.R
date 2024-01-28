#------Install Packages

install.packages("pracma")
install.packages("MASS")
install.packages("mvtnorm")
install.packages("optimr")
install.packages("bbmle")
install.packages("mle.tools")
install.packages("plot3D")
install.packages("rworldmap")
install.packages("geosphere")

#----Library Packages

library(pracma)
library(MASS)
library(mvtnorm)
library(optimr)
library(bbmle)
library(mle.tools)
library(plot3D)
library(rworldmap)
library(geosphere)

#------fBm covariance function

c_H<-function(u,v,h){
  return(0.5*(u^(2*h)+v^(2*h)-abs(u-v)^(2*h)))
 }


#-----Covariance funtion y

cov_OUH_zeta_2<-function(T,n,b,h) {
    D<-T/n
    index<-0
    G<-function(u,v){return(exp(-b*(u+v))*abs(v-u+D*index)^(2*h))}
    H<-function(u){return(exp(-b*u)*(D*(index+1)-u)^(2*h))}
    int_2<-rep(0,n)
    int_1<-rep(0,n)
    for(i in 1:n){
    int_2[i]<-integral2(G,0,D,0,D,reltol = 1e-6)$Q
    int_1[i]<-integrate(H,0,D)$val
    index<-i 
    }
    K<-matrix(rep(0,n*n),n,n)
    for(i in 1:n) {
        for(j in 1:i){
              K[i,j]<-0.5*((1-exp(-b*D))/b)*(int_1[i]+int_1[j])-0.5*int_2[i-j+1]
              K[j,i]<-K[i,j]
        }
    }
    return(K)
 }


#-----Simulation velocity given position

simulation_vel<-function(T,s,b,h,mu,index,sim){
    n<-index[length(index)]
    D<-T/n
    K_2<-cov_OUH_zeta_2(T,n,b,h)    
    z_2<-rep(0,n)
    S_B<-matrix(rep(0,n*n),n,n)
    K_1<-matrix(rep(0,n*n),n,n)
    for(i in 1:n){
      for(j in 1:i){
          K_1[i,j]<-s*s*(c_H(j*D,i*D,h)-exp(-b*D)*(c_H(j*D,(i-1)*D,h)+c_H(i*D,(j-1)*D,h))+exp(-2*b*D)*c_H((i-1)*D,(j-1)*D,h))
          K_1[j,i]<-K_1[i,j]
      }
    }
   index_aux<-0
   H<-function(u){return(exp(-b*u)*((index_aux+1)*D-u)^(2*h))}
   int_2<-rep(0,2*n)
   int_1<-rep(0,n)
   for(i in 1:n){
    int_1[i]<-integrate(H,0,D)$val
    index_aux<-i 
   }
   G<-function(u){return(exp(b*(u-D))*abs(index_aux*D-u)^(2*h))}
   for(i in (-n+1):n){
    index_aux<-i
    int_2[n+i]<-integrate(G,0,D)$val 
   }
   K_21<-matrix(rep(0,n*n),n,n)
    for(i in 1:n) {
        for(j in 1:n){
              K_21[i,j]<-( -0.5*s*s*b*(  ( (j*D)^(2*h)       )*( (1-exp(-b*D))/b )+int_1[i]-int_2[j-i+1+n]          )
                           +0.5*s*s*b*exp(-1*b*D)*(  ( ( (j-1)*D )^(2*h) )*( (1-exp(-b*D))/b )+int_1[i]-int_2[j-i+n]           )     )
       }
    }
 
  for(i in 1:n) {
        for(j in 1:i){
             S_B[i,j]<-exp((j-i)*b*D)
       }
    }
  M<-(s*S_B%*%K_2%*%t(s*S_B))
  M_2<-M[index,index]
  M_1<-M[-index,-index]
  M_21<-M[index,-index]
  M_12<-M[-index,index]
  K_12<-matrix(rep(0,n*n),n,n)  
  K_12<-t(K_21)   
  vel_aux<-rep(0,n)
  vel<-rep(0,n)
  mu_com<-rep(0,n)
#---------------Simulation----
  mu_vel_sim<-matrix(rep(0,n*2*sim),2*sim,n)
  for(i in 1:sim){
  mu_vel_sim[2*i-1,-index]<-mvrnorm( 1,rep(mu[1],n-length(index))+M_12%*%solve(M_2)%*%as.vector(mu[-1]-mu[1]), M_1-M_12%*%solve(M_2)%*%M_21)
  mu_vel_sim[2*i-1,index]<-mu[-1]         #--Es sin el punto inicial, para imprimir debemos poner el punto inicial 
  z_2<--b*solve(S_B)%*%(as.vector(mu_vel_sim[2*i-1,]-mu[1])) #---con los cambios creo que es sin exp(-b*D)
  z_1<-mvrnorm( 1,K_12%*%solve(s*s*b*b*K_2)%*%as.vector(z_2), K_1-K_12%*%solve(s*s*b*b*K_2)%*%K_21)
  vel_aux<-z_1+z_2
  mu_vel_sim[2*i,1]<-vel_aux[1]
  for(j in 2:n){
    mu_vel_sim[2*i,j]<-vel_aux[j]+exp(-b*D)*mu_vel_sim[2*i,j-1]
  }

 }
  return(mu_vel_sim)
}





#--------Return sig of number

sig<-function(x){
  if(x<0){return(-1)}
  if(x>=0){return(1)}
}


#----Change longitude and latitude for distance on meters

degree_km<-function(lon,lat,vel_lon,vel_lat){
k<-length(lon)
change<-matrix(rep(0,4*k),4,k) 
change[1,1]<-0
change[2,1]<-0
change[3,1]<-0
change[4,1]<-0
for(i in 2:k){
change[1,i]<-sig(lon[i]-lon[1])*distm(c(lon[1], lat[1]), c(lon[i], lat[1]), fun = distHaversine)
change[2,i]<-sig(lat[i]-lat[1])*distm(c(lon[1], lat[1]), c(lon[1], lat[i]), fun = distHaversine)
change[3,i]<-vel_lon[i]*distm(c(lon[i], lat[i]), c(lon[i]+1, lat[i]), fun = distHaversine)
change[4,i]<-vel_lat[i]*distm(c(lon[i], lat[i]), c(lon[i], lat[i]+1), fun = distHaversine)
}
  return(change)
}


#------Simulation one dimensional trayectory given covariance matriz and parameters

simulation<-function(s,b,D,K,mu_ini){
    n<-length(K[1,])
    sim<-mvrnorm( 1,rep(0,n), K )
    mu<-rep(0,n+1)
    mu[1]<-mu_ini
    for(i in 1:n){
      mu[i+1]=(1-exp(-b*D))*mu[1]+s*sim[i]+exp(-b*D)*mu[i]
     }
    return(mu)
}

#------Log-likelihood function-----

log_verosimilitude<-function(T,b,h,mu){
n<-length(mu)-1
D<-T/n
K<-cov_OUH_zeta_2(T,n,b,h)
S_B<-matrix(rep(0,n*n),n,n)
for(i in 1:n){
  for(j in 1:i){
     S_B[i,j]<-exp((j-i)*b*D)
  } 
}

M<-(S_B%*%K%*%t(S_B))
s<-(t(as.vector(mu[2:(n+1)]-mu[1]))%*%solve(M)%*%as.vector(mu[2:(n+1)]-mu[1]) )[1,1]/n
S<-s*M   #----No hay estabilidad numerico con y por lo tanto consideramos S_B
return(dmvnorm(mu[2:(n+1)],rep(mu[1],n),S,log=TRUE))
}



#------Log-likelihood position with index (for practical case)--------#
                    
log_verosimilitude_practical<-function(T,b,h,mu,index){
k<-length(mu)-1
n<-index[k]                
D<-T/n                     
K<-cov_OUH_zeta_2(T,n,b,h)
S_B<-matrix(rep(0,n*n),n,n)
for(i in 1:n){
  for(j in 1:i){
     S_B[i,j]<-exp((j-i)*D*b)
  } 
}
M<-(S_B%*%K%*%t(S_B))[index,index]
g<-(t(as.vector(mu[2:(k+1)]-mu[1]))%*%solve(M)%*%as.vector(mu[2:(k+1)]-mu[1]) )[1,1]/k
S<-g*M   #----No hay estabilidad numerico con y por lo tanto consideramos S_B
return(dmvnorm(mu[2:(k+1)],rep(mu[1],k),S,log=TRUE))
}


#------Prediction function Given data and MLE 

prediction<-function(T,s,b,h,mu,sim,step){
n<-length(mu)-1    #----sim numero de simulaciones
D<-T/n         #----step numero de pasos a predecir
K<-cov_OUH_zeta_2(T,n+step,b,h)
y<-rep(0,n)
for(i in 1:n){
    y[i]<--((1-exp(-b*D))*mu[1]-mu[i+1]+exp(-b*D)*mu[i])/s
  }

K_11<-K[1:n,1:n]
K_12<-K[1:n,(n+1):(n+step)]
K_21<-t(K[1:n,(n+1):(n+step)])
K_22<-K[(n+1):(n+step),(n+1):(n+step)]
K_11_inv<-solve(K_11)

y_predic<-rmvnorm(sim,K_21%*%K_11_inv%*%y,K_22-K_21%*%K_11_inv%*%K_12)

mu_predic<-matrix(rep(0,step*sim),sim,step)


for(i in 1:sim){
 mu_aux<-c(mu,rep(0,step))
 for(j in 1:step){
  mu_aux[n+j+1]<-(1-exp(-b*D))*mu_aux[1]+exp(-b*D)*mu_aux[n+j]+s*y_predic[i,j]
 }
 mu_predic[i,]<-mu_aux[(n+2):(n+step+1)]
}

return(mu_predic)

}


#-----Log-Likelihood function with data fixed

log_vero<-function(x){
-log_verosimilitude(T,b=x[1],h=x[2],mu=mu)
}






#------------Form of covariance function of position

T<-1000
n<-50
t<-seq(0,T,length=n)
D<-T/n

s_1<-3
b_1<-0.1
h_1<-0.75

K<-cov_OUH_zeta_2(T,n,b_1,h_1)
S_B<-matrix(rep(0,n*n),n,n)
for(i in 1:n){
  for(j in 1:i){
     S_B[i,j]<-exp((j-i)*b_1*D)
  } 
}

M<-(s_1*S_B%*%K%*%t(s_1*S_B))

plot(t,M[1,],type="l",lwd=1,col=1,xlim=c(0,T+0.5),ylim=c(0,28000000),ylab="Q(s,t)",main="Covariance function")
for(i in c(2:(n/5-1))){
	lines(t,M[i*5,],col=i)
   
}
legend(0, 28000000, legend=c(round(D+t[c(1,5*c(2:(n/5-1)))],2)),
       col=c(1:(n/5)), lty=1, cex=0.8)




#------------------------------------Empirical asymptotic study

time <- proc.time()

T<-20
n<-250

s<-2
b<-3
h<-0.8
mu_ini<-2
t<-seq(0,T,length=n)
D<-t[2]-t[1]
K<-cov_OUH_zeta_2(T,n,b,h)

m<-7
max<-matrix(rep(0,m*3),m,3)

#-------------Maximos -------------
y<-rep(0,(n-1))
for(i in 1:m){
 mu<-simulation(s,b,D,K,mu_ini)
 maximos_2<-opm(c(4,0.55),fn=log_vero, lower=c(0.05,0.505), upper=c(15,0.99),method="L-BFGS-B")
 max[i,2]<-maximos_2$p1
 max[i,3]<-maximos_2$p2
 K_2<-cov_OUH_zeta_2(T,n-1,max[i,2],max[i,3])
 for(j in 1:(n-1)){
   y[j]<-(exp(max[i,2]*D)-1)*mu[1]-exp(max[i,2]*D)*mu[j+1]+mu[j]
   }
 max[i,1]<-sqrt(-dmvnorm(sqrt(2/(n-1))*y,rep(0,n-1),K_2,log=TRUE)+dmvnorm(rep(0,n-1),rep(0,n-1),K_2,log=TRUE))
}
write.csv(max,"Data2.csv")
proc.time()-time







#-----------------------Simulation with predictions

T<-10
n<-200

s_1<-sqrt(3)
b_1<-12
h_1<-0.56
mu_ini_1<-15
vel_ini_1<-0


s_2<-sqrt(7)
b_2<-6
h_2<-0.75
mu_ini_2<-10
vel_ini_2<-0


t<-seq(0,T,length=(n+1))
D<-t[2]-t[1]


K_1<-cov_OUH_zeta_2(T,n,b_1,h_1)
K_2<-cov_OUH_zeta_2(T,n,b_2,h_2)
mu_1<-simulation(s_1,b_1,D,K_1,mu_ini_1)
mu_2<-simulation(s_2,b_2,D,K_2,mu_ini_2)





#---Trajectory
m<-60
plot(mu_1,mu_2, type="l",xlab="longitude",ylab="latitude",lwd=2,
main="Simulation Path h=(0.03,0.92)")
max_1<-max(mu_1)
max_2<-max(mu_2)
min_1<-min(mu_1)
min_2<-min(mu_2)
reg_1<-seq(min_1,max_1,length=m)
reg_2<-seq(min_2,max_2,length=m)
for(i in 1:m){
   lines(c(reg_1[i],reg_1[i]),c(min_2,max_2),col="green",type="l")
   lines(c(min_1,max_1),c(reg_2[i],reg_2[i]),col="green",type="l")
}
points(mu_1[1],mu_2[1],col="red",pch=19)
points(mu_1[n],mu_2[n],col="blue",pch=19)


D
b_1

T

R<-1/D
T_2<-R*T



#------------Change to lon data---------------
log_vero<-function(x){
-log_verosimilitude(T_2,x[1],x[2],mu_1)
}


#----------------Imprimirmos la perfil de h y b para tener una idea de donde puede el EMV de h


R_2<-1/R
b_a<-b_1/5
b_b<-b_1*5

#----Profile Likelihood (b,H)

time <- proc.time()
block<-20
reg_b<-seq(b_a*R_2,b_b*R_2,length=block)
reg_h<-seq(0.2,0.9,length=block)
profile_bh<-matrix(rep(0,block*block),block,block)
for(i in 1:block){
 for(j in 1:block){
   profile_bh[i,j]<-log_vero(c(reg_b[i],reg_h[j]))
 }
}
res<-persp(reg_b,reg_h,exp(-profile_bh/3),theta = 0, phi = 45,xlab="b",ylab="H",zlab="log(Lp(b,H))",main="Profile Log-Likelihood (b,H)")


proc.time()-time


#-----MLE calculations

max_lon<-rep(0,3)
max_lat<-rep(0,3)
mu<-mu_1
maximos_lon<-opm(c(4,0.25),fn=log_vero, lower=c(0.05,0.05), upper=c(15,0.95),method="L-BFGS-B")
 max_lon[2]<-maximos_lon$p1
 max_lon[3]<-maximos_lon$p2
 

K_1<-cov_OUH_zeta_2(T_2,n,max_lon[2],max_lon[3])
S_B<-matrix(rep(0,n*n),n,n)
for(i in 1:n){
  for(j in 1:i){
     S_B[i,j]<-exp((j-i)*max_lon[2]*D)
  } 
}
M<-(S_B%*%K_2%*%t(S_B))
max_lon[1]<-sqrt((t(as.vector(mu[2:(n+1)]-mu[1]))%*%solve(M)%*%as.vector(mu[2:(n+1)]-mu[1]) )[1,1]/n)


#----Escale D=1

max_lon[3]
max_lon[2]
max_lon[1]
 
 
 #----MLE calculation in Real Scale
 
K_2<-cov_OUH_zeta_2(T,n,max_lon[2]/D,max_lon[3])
 S_B<-matrix(rep(0,n*n),n,n)
for(i in 1:n){
  for(j in 1:i){
     S_B[i,j]<-exp((j-i)*max_lon[2])
  } 
}
M<-(S_B%*%K_2%*%t(S_B))
max_lon[1]<-sqrt((t(as.vector(mu[2:(n+1)]-mu[1]))%*%solve(M)%*%as.vector(mu[2:(n+1)]-mu[1]) )[1,1]/n)


#------Real scale

max_lon[3]
max_lon[2]/D
max_lon[1]



mu<-mu_2
#------------Change to lon data---------------
log_vero<-function(x){
-log_verosimilitude(T_2,x[1],x[2],mu_2)
}


#---MLE calculation

maximos_lat<-opm(c(4,0.25),fn=log_vero, lower=c(0.05,0.05), upper=c(15,0.95),method="L-BFGS-B")
 max_lat[2]<-maximos_lat$p1
 max_lat[3]<-maximos_lat$p2
 

K_1<-cov_OUH_zeta_2(T_2,n,max_lat[2],max_lat[3])
S_B<-matrix(rep(0,n*n),n,n)
for(i in 1:n){
  for(j in 1:i){
     S_B[i,j]<-exp((j-i)*max_lat[2]*D)
  } 
}
M<-(S_B%*%K_2%*%t(S_B))
max_lat[1]<-sqrt((t(as.vector(mu[2:(n+1)]-mu[1]))%*%solve(M)%*%as.vector(mu[2:(n+1)]-mu[1]) )[1,1]/n)

#---Scale D=1



max_lat[3]
max_lat[2]
max_lat[1]


#---MLE Calculation real scale
 
K_2<-cov_OUH_zeta_2(T,n,max_lat[2]/D,max_lat[3])
 S_B<-matrix(rep(0,n*n),n,n)
for(i in 1:n){
  for(j in 1:i){
     S_B[i,j]<-exp((j-i)*max_lat[2])
  } 
}
M<-(S_B%*%K_2%*%t(S_B))
max_lat[1]<-sqrt((t(as.vector(mu[2:(n+1)]-mu[1]))%*%solve(M)%*%as.vector(mu[2:(n+1)]-mu[1]) )[1,1]/n)



#----Real Scale

max_lat[3]
max_lat[2]/D
max_lat[1]


#-----


max_lon[2]<-max_lon[2]/D
max_lat[2]<-max_lat[2]/D



#-----------Predictions--------------


sim<-100
step<-10

predic_lon<-prediction(T,max_lon[1],max_lon[2],max_lon[3],mu_1,sim,step)
predic_lat<-prediction(T,max_lat[1],max_lat[2],max_lat[3],mu_2,sim,step)


medias_pred_lon<-c(0,step)
medias_pred_lat<-c(0,step)


for(i in 1:step) {
   medias_pred_lon[i]<-mean(predic_lon[,i])
   medias_pred_lat[i]<-mean(predic_lat[,i])
}


max_1<-max(mu_1)+abs(max(predic_lon)-mu_1[n])
max_2<-max(mu_2)+abs(max(predic_lat)-mu_2[n])
min_1<-min(mu_1)-abs(min(predic_lon)-mu_1[n])
min_2<-min(mu_2)-abs(min(predic_lat)-mu_2[n])
plot(mu_1,mu_2,type="l",col="black",xlab="longitude",ylab="latitude",main="Path Prediction",xlim=c(min_1,max_1),ylim=c(min_2,max_2),lwd=2)
for(i in 1:sim){
 lines(c(mu_1[n],predic_lon[i,]),c(mu_2[n],predic_lat[i,]),type="l",col="red")
}

lines(c(mu_1[n],medias_pred_lon),c(mu_2[n],medias_pred_lat),type="l",col="yellow",lwd=2)
reg_1<-seq(min_1,max_1,length=m)
reg_2<-seq(min_2,max_2,length=m)
for(i in 1:m){
   lines(c(reg_1[i],reg_1[i]),c(min_2,max_2),col="green",type="l")
   lines(c(min_1,max_1),c(reg_2[i],reg_2[i]),col="green",type="l")
}                        
points(mu_1[1],mu_2[1],col="red",pch=19)
points(mu_1[n],mu_2[n],col="blue",pch=19)


#------------------------End simulation case---------------#
















#---Practical case: Study of whales

datos<-read.csv(file.choose(),header=TRUE) #--Data

number<-3 #---Number of whale in data

cont<-0
for(j in 1:length(datos[,3*(number-1)+2])){
  if(datos[j,3*(number-1)+2]!="NA"){
  cont<-cont+1
 }
  
}

k<-cont
lon<-datos[1:k,3*(number-1)+2]
lat<-datos[1:k,3*number]
index<-datos[2:k,3*(number-1)+1] #--Index of day informations



D<-1
T<-index[k-1]*D #----initial + k-1 data
t<-datos[,1]*D
n<-T/D


#-----------------Empirical velocity

vel_lon<-rep(0,k-1)
vel_lat<-rep(0,k-1)
for(i in 1:(k-1))
{
  vel_lon[i]<-(lon[i+1]-lon[i])/(t[i+1]-t[i])
  vel_lat[i]<-(lat[i+1]-lat[i])/(t[i+1]-t[i])
}

mat <- matrix(c(1,2,
                3,3), # y tercer gráfico
              nrow = 2, ncol = 2,
              byrow = TRUE)

layout(mat = mat)

plot(t[2:k],vel_lon,type="l",col="blue",xlab="t",ylab="velocity longitude",main="velocity longitude")
plot(t[2:k],vel_lat,type="l",col="blue",xlab="t",ylab="velocity latitude",main="velocity latitude")



#----Plot data

map<-getMap(resolution = "high")
plot(map,xlim=c(-113,-110),ylim=c(22,33),main="Path")
lines(lon,lat,type="l",col="red")
points(lon,lat,pch=19,col="blue")

max_1<-max(lon)
max_2<-max(lat)
min_1<-min(lon)
min_2<-min(lat)
reg_1<-seq(min_1,max_1,length=m)
reg_2<-seq(min_2,max_2,length=m)
#for(i in 1:m){
 #  lines(c(reg_1[i],reg_1[i]),c(min_2,max_2),col="green",type="l")
  # lines(c(min_1,max_1),c(reg_2[i],reg_2[i]),col="green",type="l")
#}
lines(lon,lat,type="l",xlab="longitude",lwd=2,ylab="latitude",main="Path")
points(lon[1],lat[1],col="blue",pch=19,lwd=5)
points(lon[k],lat[k],col="Orange",pch=19,lwd=5)
points(lon[2:(k-1)],lat[2:(k-1)],col="red",pch=19)





#-----MLE for lon data



#------------Change to lon data---------------
log_vero_pract<-function(x){
-log_verosimilitude_practical(T,x[1],x[2],lon,index)
}




#------Profile Likelihood (b,H)

time <- proc.time()
block<-75
reg_b<-seq(0.4852763*24/3,0.4852763*24*3,length=block)
reg_h<-seq(0.65,0.8,length=block)
profile_bh<-matrix(rep(0,block*block),block,block)
for(i in 1:block){
 for(j in 1:block){
   profile_bh[i,j]<-log_vero_pract(c(reg_b[i],reg_h[j]))
 }
}
res<-persp(reg_b,reg_h,exp(-profile_bh*10),theta = 0, phi = 0,xlab="b",ylab="H",zlab="Lp(b,H)",main="Profile Likelihood (b,H)")
#mypoints <- trans3d(max_lon_prac[2],max_lon_prac[3],exp(-log_vero_pract(c(max_lon_prac[2],max_lon_prac[3]))*10),pmat=res)
#points(mypoints, pch=19, col="red")

proc.time()-time


#----We have a unique maximum

max_lon_prac<-rep(0,3)
maximos_lon_prac<-opm(c(5,runif(1,0.1,0.9)),fn=log_vero_pract, lower=c(0.5,0.1), upper=c(20,0.9),method="L-BFGS-B")
max_lon_prac[2]<-maximos_lon_prac$p1
max_lon_prac[3]<-maximos_lon_prac$p2
 K_2<-cov_OUH_zeta_2(T,n,max_lon_prac[2],max_lon_prac[3])
 S_B<-matrix(rep(0,n*n),n,n)
for(i in 1:n){
  for(j in 1:i){
     S_B[i,j]<-exp((j-i)*max_lon_prac[2]*D)
  } 
}
M<-(S_B%*%K_2%*%t(S_B))[index,index]
max_lon_prac[1]<-sqrt((t(as.vector(lon[2:k]-lon[1]))%*%solve(M)%*%as.vector(lon[2:k]-lon[1]) )[1,1]/(k-1))



#----We don't have a unique maximum


maximos_lon_prac<-opm(c(5,runif(1,0.1,0.9)),fn=log_vero_pract, lower=c(0.5,0.1), upper=c(20,0.9),method="L-BFGS-B")
max_lon_prac[3]<-maximos_lon_prac$p2



#----------Profile likelihood  b-----------------#

log_vero_pract<-function(x){
-log_verosimilitude_practical(T,x,max_lon_prac[3],lon,index)
}


n_1<-1000

reg_b<-seq(10,400,length=n_1)
reg_b_ver<-rep(0,length(reg_b))
for(i in 1:length(reg_b)){
reg_b_ver[i]<--log_vero_pract(reg_b[i])
}
plot(reg_b,reg_b_ver,xlab="b",ylab="L_p(b)",type="l",col="blue",main="Likelihood beta")


b_val_max<-400

b_max_flat<-function(x){
return(log_vero_pract(x)-log_vero_pract(b_val_max)-0.01) #---0.01 is tolerance value to likelihood functions
}

max_lon_prac[2]<-uniroot(b_max_flat,lower=2,upper=400,tol=1e-9)$root

#-----Plot s MLE with respect b MLE

log_verosimilitude_practical_s<-function(T,b,h,mu,index){
k<-length(mu)-1
n<-index[k]                
D<-T/n                     
K<-cov_OUH_zeta_2(T,n,b,h)
S_B<-matrix(rep(0,n*n),n,n)
for(i in 1:n){
  for(j in 1:i){
     S_B[i,j]<-exp((j-i)*D*b)
  } 
}
M<-(S_B%*%K%*%t(S_B))[index,index]
g<-(t(as.vector(mu[2:(k+1)]-mu[1]))%*%solve(M)%*%as.vector(mu[2:(k+1)]-mu[1]) )[1,1]/k
return(sqrt(g))
}

log_vero_pract_s<-function(x){
-log_verosimilitude_practical_s(T,x,max_lon_prac[3],lon,index)
}


reg<-seq(0.01,6,length=50)
reg_s<-rep(0,length(reg))
for(i in 1:length(reg)){
reg_s[i]<-log_verosimilitude_practical_s(T,reg[i],max_lon_prac[3],lon,index)
}
plot(reg,reg_s,type="l",ylab="s",xlab="b",main="Maximum sigma")


#-------Calculation s MLE


K_2<-cov_OUH_zeta_2(T,n,max_lon_prac[2],max_lon_prac[3])
 S_B<-matrix(rep(0,n*n),n,n)
for(i in 1:n){
  for(j in 1:i){
     S_B[i,j]<-exp((j-i)*max_lon_prac[2]*D)
  } 
}
M<-(S_B%*%K_2%*%t(S_B))[index,index]
max_lon_prac[1]<-sqrt((t(as.vector(lon[2:k]-lon[1]))%*%solve(M)%*%as.vector(lon[2:k]-lon[1]) )[1,1]/(k-1))




#------------Change to lon data---------------
log_vero_pract<-function(x){
-log_verosimilitude_practical(T,x[1],x[2],lon,index)
}




#---Ratio lokelihood
-2*(-log_vero_pract(c(max_lon_prac[2],0.5))+log_vero_pract(c(max_lon_prac[2],max_lon_prac[3])))






#------------Change to lat data---------------
log_vero_pract<-function(x){
-log_verosimilitude_practical(T,x[1],x[2],lat,index)
}

#------Profile Likelihood (b,H)

time <- proc.time()
block<-50
reg_b<-seq(0.05,5,length=block)
reg_h<-seq(0.4,0.92,length=block)
profile_bh_2<-matrix(rep(0,block*block),block,block)
for(i in 1:block){
 for(j in 1:block){
   profile_bh_2[i,j]<-log_vero_pract(c(reg_b[i],reg_h[j]))
 }
}
res<-persp(reg_b,reg_h,-profile_bh_2,theta = 20, phi = 55,xlab="b",ylab="H",zlab="log(Lp(b,H))",main="Profile Log-Likelihood (b,H)")


proc.time()-time

#-----------Calculation (s,b,H) MLE


max_lat_prac<-rep(0,3)
maximos_lat_prac<-opm(c(5,runif(1,0.6,0.8)),fn=log_vero_pract, lower=c(0.5,0.6), upper=c(20,0.8),method="L-BFGS-B")
max_lat_prac[2]<-maximos_lat_prac$p1
max_lat_prac[3]<-maximos_lat_prac$p2
 K_2<-cov_OUH_zeta_2(T,n,max_lat_prac[2],max_lat_prac[3])
 S_B<-matrix(rep(0,n*n),n,n)
for(i in 1:n){
  for(j in 1:i){
     S_B[i,j]<-exp((j-i)*max_lat_prac[2]*D)
  } 
}
M<-(S_B%*%K_2%*%t(S_B))[index,index]
max_lat_prac[1]<-sqrt((t(as.vector(lat[2:k]-lat[1]))%*%solve(M)%*%as.vector(lat[2:k]-lat[1]) )[1,1]/(k-1))















#-------------Distance calculate---------

sim<-40
vel_lon<-simulation_vel(T,max_lon_prac[1],max_lon_prac[2],max_lon_prac[3],lon,index,sim)
vel_lat<-simulation_vel(T,max_lat_prac[1],max_lat_prac[2],max_lat_prac[3],lat,index,sim)



mat <- matrix(c(1,2,
                3,4), # y tercer gráfico
              nrow = 2, ncol = 2,
              byrow = TRUE)

layout(mat = mat)
#---Simulation mean lon trajectory

media_lon_complet<-rep(0,n)
for(i in 1:n){media_lon_complet[i]<-mean(vel_lon[2*(1:sim)-1,i])}
plot(c((0:n)*D),c(lon[1],vel_lon[1,]),type="l",main="Estimation Position",,col="blue",xlab="t",ylab="Position longitude")
for(i in 2:sim){
lines(c((0:n)*D),c(lon[1],vel_lon[2*i-1,]),type="l",col="blue")
}
points(c(0,index*D),lon,col="red",pch=19)
lines(c((0:n)*D),c(lon[1],media_lon_complet),type="l",col="Green",lwd=2)


#---Simulation mean lat trajectory

media_lat_complet<-rep(0,n)
for(i in 1:n){media_lat_complet[i]<-mean(vel_lat[2*(1:sim)-1,i])}
plot(c((0:n)*D),c(lat[1],vel_lat[1,]),type="l",main="Estimation Position",,col="blue",xlab="t",ylab="Position latitude")
for(i in 2:sim){
lines(c((0:n)*D),c(lat[1],vel_lat[2*i-1,]),type="l",col="blue")
}
points(c(0,index*D),lat,col="red",pch=19)
lines(c((0:n)*D),c(lat[1],media_lat_complet),type="l",col="Green",lwd=2)


#---Simulation mean lon velocity


media_lon_vel_complet<-rep(0,n)
for(i in 1:n){media_lon_vel_complet[i]<-mean(vel_lon[2*(1:sim),i])}
plot(c((1:n)*D),vel_lon[2,],type="l",main="Estimation velocity longitude",,col="gray",xlab="t",ylab="velocity (degree longitude)",ylim=c(-0.03,0.03))
for(i in 2:sim){
lines(c((1:n)*D),vel_lon[2*i,],type="l",col="gray")
}
lines(c((1:n)*D),media_lon_vel_complet,type="l",col="Green",lwd=2)


#---Simulation mean lat velocity


media_lat_vel_complet<-rep(0,n)
for(i in 1:n){media_lat_vel_complet[i]<-mean(vel_lat[2*(1:sim),i])}
plot(c((1:n)*D),vel_lat[2,],type="l",main="Estimation velocity latitude",,col="gray",xlab="t",ylab="velocity (degree latitude)",ylim=c(-0.03,0.04))
for(i in 2:sim){
lines(c((1:n)*D),vel_lat[2*i,],type="l",col="gray")
}
lines(c((1:n)*D),media_lat_vel_complet,type="l",col="Green",lwd=2)


#--------Trejectory on map

map<-getMap(resolution = "high")
plot(map,xlim=c(-113,-110),ylim=c(22,33))
med<-length(lon)
points(lon[1:med],lat[1:med],xlab="longitude",pch=19,col="blue")
points(lon[1],lat[1],col="Green",pch=19)
points(lon[length(lon)],lat[length(lat)],col="Red",pch=19)


lines(c(lon[1],vel_lon[1,]),c(lat[1],vel_lat[1,]),type="l")
for(i in 2:sim){
lines(c(lon[1],vel_lon[2*i-1,]),c(lat[1],vel_lat[2*i-1,]),type="l",col="blue")
}
points(lon,lat,col="red",pch=19)
lines(c(lon[1],media_lon_complet),c(lat[1],media_lat_complet),type="l",col="Green",lwd=2)



#--------------Change degree to km--------------


mat <- matrix(c(1,2,
                3,4), # y tercer gráfico
              nrow = 2, ncol = 2,
              byrow = TRUE)

layout(mat = mat)



estimation_km<-degree_km(c(lon[1],vel_lon[1,]),c(lat[1],vel_lat[1,]),c(0,vel_lon[2,]),c(0,vel_lat[2,]))
plot(c(0:n),estimation_km[1,]/10^3,xlab="Day",ylab="Km",type="l",col="Gray",main="Longitude km/d")
for(i in 2:sim){
estimation_km<-degree_km(c(lon[1],vel_lon[2*i-1,]),c(lat[1],vel_lat[2*i-1,]),c(0,vel_lon[2*i,]),c(0,vel_lat[2*i,]))
lines(c(0:n),estimation_km[1,]/10^3,type="l",col="Gray")
}

estimation_km<-degree_km(c(lon[1],media_lon_complet),c(lat[1],media_lat_complet),c(0,media_lon_vel_complet),c(0,media_lat_vel_complet))
lines(c(0:n),estimation_km[1,]/10^3,type="l",col="Green",lwd="2")
points(c(0:n),estimation_km[1,]/10^3,col="Blue",pch=19)
points(index,estimation_km[1,index+1]/10^3,col="Red",pch=19)
lines(c(48,48),c(-300,100),col="red")
lines(c(58,58),c(-300,100),col="red")
lines(c(129,129),c(-300,100),col="red")
lines(c(158,158),c(-300,100),col="red")
points(c(48:58),estimation_km[1,48:58+1]/10^3,col="Orange",pch=19)
points(c(129:158),estimation_km[1,129:158+1]/10^3,col="Yellow",pch=19)



estimation_km<-degree_km(c(lon[1],vel_lon[1,]),c(lat[1],vel_lat[1,]),c(0,vel_lon[2,]),c(0,vel_lat[2,]))
plot(c(0:n),estimation_km[2,]/10^3,,xlab="Day",ylab="Km",type="l",col="Gray",main="Latitude km/d")
for(i in 2:sim){
estimation_km<-degree_km(c(lon[1],vel_lon[2*i-1,]),c(lat[1],vel_lat[2*i-1,]),c(0,vel_lon[2*i,]),c(0,vel_lat[2*i,]))
lines(c(0:n),estimation_km[2,]/10^3,type="l",col="Gray")
}

estimation_km<-degree_km(c(lon[1],media_lon_complet),c(lat[1],media_lat_complet),c(0,media_lon_vel_complet),c(0,media_lat_vel_complet))
lines(c(0:n),estimation_km[2,]/10^3,type="l",col="Green",lwd="2")
points(c(0:n),estimation_km[2,]/10^3,col="Blue",pch=19)
points(index,estimation_km[2,index+1]/10^3,col="Red",pch=19)
lines(c(48,48),c(-100,500),col="red")
lines(c(58,58),c(-100,500),col="red")
lines(c(129,129),c(-100,500),col="red")
lines(c(158,158),c(-100,500),col="red")
points(c(48:58),estimation_km[2,48:58+1]/10^3,col="Orange",pch=19)
points(c(129:158),estimation_km[2,129:158+1]/10^3,col="Yellow",pch=19)




estimation_km<-degree_km(c(lon[1],vel_lon[1,]),c(lat[1],vel_lat[1,]),c(0,vel_lon[2,]),c(0,vel_lat[2,]))
plot(c(0:n),D*estimation_km[3,]/10^3,,xlab="Day",ylab="Km",type="l",col="Gray",main="Velocity longitude km/d")
for(i in 2:sim){
estimation_km<-degree_km(c(lon[1],vel_lon[2*i-1,]),c(lat[1],vel_lat[2*i-1,]),c(0,vel_lon[2*i,]),c(0,vel_lat[2*i,]))
lines(c(0:n),D*estimation_km[3,]/10^3,type="l",col="Gray")
}

estimation_km<-degree_km(c(lon[1],media_lon_complet),c(lat[1],media_lat_complet),c(0,media_lon_vel_complet),c(0,media_lat_vel_complet))
lines(c(0:n),D*estimation_km[3,]/10^3,type="l",col="Green",lwd="2")
lines(c(48,48),c(-60,40),col="red")
lines(c(58,58),c(-60,40),col="red")
lines(c(129,129),c(-60,40),col="red")
lines(c(158,158),c(-60,40),col="red")



estimation_km<-degree_km(c(lon[1],vel_lon[1,]),c(lat[1],vel_lat[1,]),c(0,vel_lon[2,]),c(0,vel_lat[2,]))
plot(c(0:n),D*estimation_km[4,]/10^3,,xlab="Day",ylab="Km",type="l",col="Gray",main="Velocity latitude km/d")
for(i in 2:sim){
estimation_km<-degree_km(c(lon[1],vel_lon[2*i-1,]),c(lat[1],vel_lat[2*i-1,]),c(0,vel_lon[2*i,]),c(0,vel_lat[2*i,]))
lines(c(0:n),D*estimation_km[4,]/10^3,type="l",col="Gray")
}

estimation_km<-degree_km(c(lon[1],media_lon_complet),c(lat[1],media_lat_complet),c(0,media_lon_vel_complet),c(0,media_lat_vel_complet))
lines(c(0:n),D*estimation_km[4,]/10^3,type="l",col="Green",lwd="2")
lines(c(48,48),c(-50,100),col="red")
lines(c(58,58),c(-50,100),col="red")
lines(c(129,129),c(-50,100),col="red")
lines(c(158,158),c(-50,100),col="red")


