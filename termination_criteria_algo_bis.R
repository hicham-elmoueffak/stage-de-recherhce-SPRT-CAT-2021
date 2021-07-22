D=1.702

########################################################### création des banques

# création bank 1
nbank1=500
a1=rlnorm(nbank1, meanlog = -0.1, sdlog = 0.3)
b1=runif(nbank1,-4,4)
c1=0.5*runif(nbank1,0.1,0.225)+0.5*runif(nbank1,0.075,0.35)

bank1=matrix(c(a1,b1,c1),nrow=3,ncol=nbank1,byrow=T)

# création bank 2
nbank2=500
npeaked2=200
a2=rlnorm(nbank2, meanlog = -0.1, sdlog = 0.3)
b2=sample(c(runif(nbank2-npeaked2,-4,4),rnorm(npeaked2,0,1)))
c2=0.5*runif(nbank2,0.1,0.225)+0.5*runif(nbank2,0.075,0.35)

bank2=matrix(c(a2,b2,c2),nrow=3,ncol=nbank2,byrow=T)

# création bank 3
nbank3=100
a3=rlnorm(nbank3, meanlog = -0.1, sdlog = 0.3)
b3=runif(nbank3,-4,4)
c3=0.5*runif(nbank3,0.1,0.225)+0.5*runif(nbank3,0.075,0.35)

bank3=matrix(c(a3,b3,c3),nrow=3,ncol=nbank3,byrow=T)

# création bank 4
nbank4=100
npeaked4=50
a4=rlnorm(nbank4, meanlog = -0.1, sdlog = 0.3)
b4=sample(c(runif(nbank4-npeaked4,-4,4),rnorm(npeaked4,0,1)))
c4=0.5*runif(nbank4,0.1,0.225)+0.5*runif(nbank4,0.075,0.35)

bank4=matrix(c(a4,b4,c4),nrow=3,ncol=nbank4,byrow=T)

############################################################### fonctions utiles

# vecteur de probas pour toute la banque d'items
presponsevect=function(bank,PL,theta){
  if(PL==2){
    exp(bank[1,]*(theta-bank[2,]))/(1+exp(bank[1,]*(theta-bank[2,])))
  }
  else if(PL==3){
    bank[3,]+(1-bank[3,])*exp(D*bank[1,]*(theta-bank[2,]))/(1+exp(D*bank[1,]*(theta-bank[2,])))
  }
}

# vecteur de probas pour les items donnés
presponse=function(bank,PL,i,theta){
  if(PL==2){
    exp(bank[1,i]*(theta-bank[2,i]))/(1+exp(bank[1,i]*(theta-bank[2,i])))
  }
  else if(PL==3){
    bank[3,i]+(1-bank[3,i])*exp(D*bank[1,i]*(theta-bank[2,i]))/(1+exp(D*bank[1,i]*(theta-bank[2,i])))
  }
}

# vecteur des infos de Fisher pour toute la banque d'items
fisherinfovect=function(bank,PL,theta){
  if(PL==2){
    (bank[1,]^2)*presponsevect(bank,PL,theta)*(1-presponsevect(bank,PL,theta))
  }
  else if(PL==3){
    (D^2*bank[1,]^2*exp(D*bank[1,]*(theta-bank[2,]))^2*(1-bank[3,]))/((1+exp(D*bank[1,]*(theta-bank[2,])))^2*(bank[3,]+exp(D*bank[1,]*(theta-bank[2,]))))
  }
}

# fonction qui renvoie la fonction llh correspondante aux observations et aux items donnés
llh=function(bank,PL,x,i){
  if(PL==2){
    llhfunction=function(theta){
      sum(log(1-presponse(bank,PL,i,theta))+x*(bank[1,i]*(theta-bank[2,i])))
    }
  }
  else if(PL==3){
    llhfunction=function(theta){
      sum(log(1-presponse(bank,PL,i,theta))+x*log((bank[3,i]+exp(D*bank[1,i]*(theta-bank[2,i])))/(1-bank[3,i])))
    }
  }
}

# SEM correspondante aux items donnés
SEM=function(bank,PL,x,i,theta){
  if(PL==2){
    1/(sqrt(sum((bank[1,i]^2)*presponse(bank,PL,i,theta)*(1-presponse(bank,PL,i,theta)))))
  }
  else if(PL==3){
    1/(sqrt(sum(abs((D^2)*(bank[1,i]^2)*exp(D*bank[1,i]*(theta-bank[2,i]))*((1/((1+exp(D*bank[1,i]*(theta-bank[2,i])))^2))-((bank[3,i]*x)/((bank[3,i]+exp(D*bank[1,i]*(theta-bank[2,i])))^2)))))))
  }
}

############################################################ estimation de theta

# nombre d'itérations et valeur de theta estimée en fonction du critère d'arrêt
thetaestimate=function(arret,bank,PL,theta){
  theta_0=0
  i_0=c(which.max(fisherinfovect(bank,PL,theta_0)))
  if (runif(1,0,1)<=presponse(bank,PL,i_0,theta)){
    theta_0=0.5
    x_0=c(1)
  }
  else {
    theta_0=-0.5
    x_0=c(0)
  }
  if(arret==1){
    while(SEM(bank,PL,x_0,i_0,theta_0)>=0.385 & length(x_0)<100){
      i=which.max(replace(fisherinfovect(bank,PL,theta_0),i_0,0))
      i_0=append(i_0,i)
      if (runif(1,0,1)<=presponse(bank,PL,i,theta)){
        x_0=append(x_0,1)
        if (length(x_0[x_0==1])==length(x_0)){
          theta_0=theta_0+0.5
        }
        else{
          theta_0=optimize(llh(bank,PL,x_0,i_0),c(-4,4),maximum=TRUE)$maximum
        }
      }
      else{
        x_0=append(x_0,0)
        if(length(x_0[x_0==0])==length(x_0)){
          theta_0=theta_0-0.5
        }
        else{
          theta_0=optimize(llh(bank,PL,x_0,i_0),c(-4,4),maximum=TRUE)$maximum
        }
      }
    }
    c(length(x_0),theta_0)
  }
  else if(arret==2){
    while(SEM(bank,PL,x_0,i_0,theta_0)>=0.315 & length(x_0)<100){
      i=which.max(replace(fisherinfovect(bank,PL,theta_0),i_0,0))
      i_0=append(i_0,i)
      if (runif(1,0,1)<=presponse(bank,PL,i,theta)){
        x_0=append(x_0,1)
        if (length(x_0[x_0==1])==length(x_0)){
          theta_0=theta_0+0.5
        }
        else{
          theta_0=optimize(llh(bank,PL,x_0,i_0),c(-4,4),maximum=TRUE)$maximum
        }
      }
      else{
        x_0=append(x_0,0)
        if(length(x_0[x_0==0])==length(x_0)){
          theta_0=theta_0-0.5
        }
        else{
          theta_0=optimize(llh(bank,PL,x_0,i_0),c(-4,4),maximum=TRUE)$maximum
        }
      }
    }
    c(length(x_0),theta_0)
  }
  else if(arret==3){
    while(SEM(bank,PL,x_0,i_0,theta_0)>=0.220 & length(x_0)<100){
      i=which.max(replace(fisherinfovect(bank,PL,theta_0),i_0,0))
      i_0=append(i_0,i)
      if (runif(1,0,1)<=presponse(bank,PL,i,theta)){
        x_0=append(x_0,1)
        if (length(x_0[x_0==1])==length(x_0)){
          theta_0=theta_0+0.5
        }
        else{
          theta_0=optimize(llh(bank,PL,x_0,i_0),c(-4,4),maximum=TRUE)$maximum
        }
      }
      else{
        x_0=append(x_0,0)
        if(length(x_0[x_0==0])==length(x_0)){
          theta_0=theta_0-0.5
        }
        else{
          theta_0=optimize(llh(bank,PL,x_0,i_0),c(-4,4),maximum=TRUE)$maximum
        }
      }
    }
    c(length(x_0),theta_0)
  }
  else if(arret==4){
    while((TRUE %in% (replace(fisherinfovect(bank,PL,theta_0),i_0,0)>=0.2))==TRUE & length(x_0)<100){
      i=which.max(replace(fisherinfovect(bank,PL,theta_0),i_0,0))
      i_0=append(i_0,i)
      if (runif(1,0,1)<=presponse(bank,PL,i,theta)){
        x_0=append(x_0,1)
        if (length(x_0[x_0==1])==length(x_0)){
          theta_0=theta_0+0.5
        }
        else{
          theta_0=optimize(llh(bank,PL,x_0,i_0),c(-4,4),maximum=TRUE)$maximum
        }
      }
      else{
        x_0=append(x_0,0)
        if(length(x_0[x_0==0])==length(x_0)){
          theta_0=theta_0-0.5
        }
        else{
          theta_0=optimize(llh(bank,PL,x_0,i_0),c(-4,4),maximum=TRUE)$maximum
        }
      }
    }
    c(length(x_0),theta_0)
  }
  else if(arret==5){
    while((TRUE %in% (replace(fisherinfovect(bank,PL,theta_0),i_0,0)>=0.1))==TRUE & length(x_0)<100){
      i=which.max(replace(fisherinfovect(bank,PL,theta_0),i_0,0))
      i_0=append(i_0,i)
      if (runif(1,0,1)<=presponse(bank,PL,i,theta)){
        x_0=append(x_0,1)
        if (length(x_0[x_0==1])==length(x_0)){
          theta_0=theta_0+0.5
        }
        else{
          theta_0=optimize(llh(bank,PL,x_0,i_0),c(-4,4),maximum=TRUE)$maximum
        }
      }
      else{
        x_0=append(x_0,0)
        if(length(x_0[x_0==0])==length(x_0)){
          theta_0=theta_0-0.5
        }
        else{
          theta_0=optimize(llh(bank,PL,x_0,i_0),c(-4,4),maximum=TRUE)$maximum
        }
      }
    }
    c(length(x_0),theta_0)
  }
  else if(arret==6){
    while((TRUE %in% (replace(fisherinfovect(bank,PL,theta_0),i_0,0)>=0.01))==TRUE & length(x_0)<100){
      i=which.max(replace(fisherinfovect(bank,PL,theta_0),i_0,0))
      i_0=append(i_0,i)
      if (runif(1,0,1)<=presponse(bank,PL,i,theta)){
        x_0=append(x_0,1)
        if (length(x_0[x_0==1])==length(x_0)){
          theta_0=theta_0+0.5
        }
        else{
          theta_0=optimize(llh(bank,PL,x_0,i_0),c(-4,4),maximum=TRUE)$maximum
        }
      }
      else{
        x_0=append(x_0,0)
        if(length(x_0[x_0==0])==length(x_0)){
          theta_0=theta_0-0.5
        }
        else{
          theta_0=optimize(llh(bank,PL,x_0,i_0),c(-4,4),maximum=TRUE)$maximum
        }
      }
    }
    c(length(x_0),theta_0)
  }
  else if(arret==7){
    while(SEM(bank,PL,x_0,i_0,theta_0)>=0.315 & (TRUE %in% (replace(fisherinfovect(bank,PL,theta_0),i_0,0)>=0.1))==TRUE & length(x_0)<100){
      i=which.max(replace(fisherinfovect(bank,PL,theta_0),i_0,0))
      i_0=append(i_0,i)
      if (runif(1,0,1)<=presponse(bank,PL,i,theta)){
        x_0=append(x_0,1)
        if (length(x_0[x_0==1])==length(x_0)){
          theta_0=theta_0+0.5
        }
        else{
          theta_0=optimize(llh(bank,PL,x_0,i_0),c(-4,4),maximum=TRUE)$maximum
        }
      }
      else{
        x_0=append(x_0,0)
        if(length(x_0[x_0==0])==length(x_0)){
          theta_0=theta_0-0.5
        }
        else{
          theta_0=optimize(llh(bank,PL,x_0,i_0),c(-4,4),maximum=TRUE)$maximum
        }
      }
    }
    c(length(x_0),theta_0)
  }
  else if(arret==8){
    while(SEM(bank,PL,x_0,i_0,theta_0)>=0.22 & (TRUE %in% (replace(fisherinfovect(bank,PL,theta_0),i_0,0)>=0.01))==TRUE & length(x_0)<100){
      i=which.max(replace(fisherinfovect(bank,PL,theta_0),i_0,0))
      i_0=append(i_0,i)
      if (runif(1,0,1)<=presponse(bank,PL,i,theta)){
        x_0=append(x_0,1)
        if (length(x_0[x_0==1])==length(x_0)){
          theta_0=theta_0+0.5
        }
        else{
          theta_0=optimize(llh(bank,PL,x_0,i_0),c(-4,4),maximum=TRUE)$maximum
        }
      }
      else{
        x_0=append(x_0,0)
        if(length(x_0[x_0==0])==length(x_0)){
          theta_0=theta_0-0.5
        }
        else{
          theta_0=optimize(llh(bank,PL,x_0,i_0),c(-4,4),maximum=TRUE)$maximum
        }
      }
    }
    c(length(x_0),theta_0)
  }
  else if(arret==9){
    change=0.5
    while(change>=0.05 & length(x_0)<100 | length(x_0)<11){
      i=which.max(replace(fisherinfovect(bank,PL,theta_0),i_0,0))
      i_0=append(i_0,i)
      if (runif(1,0,1)<=presponse(bank,PL,i,theta)){
        x_0=append(x_0,1)
        if (length(x_0[x_0==1])==length(x_0)){
          theta_0=theta_0+0.5
        }
        else{
          thetaopt=optimize(llh(bank,PL,x_0,i_0),c(-4,4),maximum=TRUE)$maximum
          change=abs(theta_0-thetaopt)
          theta_0=thetaopt
        }
      }
      else{
        x_0=append(x_0,0)
        if(length(x_0[x_0==0])==length(x_0)){
          theta_0=theta_0-0.5
        }
        else{
          thetaopt=optimize(llh(bank,PL,x_0,i_0),c(-4,4),maximum=TRUE)$maximum
          change=abs(theta_0-thetaopt)
          theta_0=thetaopt
        }
      }
    }
    c(length(x_0),theta_0)
  }
  else if(arret==10){
    change=0.5
    while(change>=0.02 & length(x_0)<100 | length(x_0)<11){
      i=which.max(replace(fisherinfovect(bank,PL,theta_0),i_0,0))
      i_0=append(i_0,i)
      if (runif(1,0,1)<=presponse(bank,PL,i,theta)){
        x_0=append(x_0,1)
        if (length(x_0[x_0==1])==length(x_0)){
          theta_0=theta_0+0.5
        }
        else{
          thetaopt=optimize(llh(bank,PL,x_0,i_0),c(-4,4),maximum=TRUE)$maximum
          change=abs(theta_0-thetaopt)
          theta_0=thetaopt
        }
      }
      else{
        x_0=append(x_0,0)
        if(length(x_0[x_0==0])==length(x_0)){
          theta_0=theta_0-0.5
        }
        else{
          thetaopt=optimize(llh(bank,PL,x_0,i_0),c(-4,4),maximum=TRUE)$maximum
          change=abs(theta_0-thetaopt)
          theta_0=thetaopt
        }
      }
    }
    c(length(x_0),theta_0)
  }
}

####################################################### grandeurs significatives

# mean bias
meanbias=function(truethetas,estimatedthetas){
  sum(estimatedthetas-truethetas)/length(truethetas)
}

# RMSE
RMSE=function(truethetas,estimatedthetas){
  sqrt(sum((estimatedthetas-truethetas)^2)/length(truethetas))
}

##################################################### pertinence des estimateurs

# performance des CATs pour les différents critères d'arrêt
performance=function(bank,PL){
  thetas=seq(-3,3,0.5)
  truethetas=c()
  for(i in thetas){
    truethetas=append(truethetas,rep(i,1000))
  }
  iterations1=c()
  iterations2=c()
  iterations3=c()
  iterations4=c()
  iterations5=c()
  iterations6=c()
  iterations7=c()
  iterations8=c()
  iterations9=c()
  iterations10=c()
  estimatedthetas1=c()
  estimatedthetas2=c()
  estimatedthetas3=c()
  estimatedthetas4=c()
  estimatedthetas5=c()
  estimatedthetas6=c()
  estimatedthetas7=c()
  estimatedthetas8=c()
  estimatedthetas9=c()
  estimatedthetas10=c()
  for(i in thetas){
    thetas_i=rep(i,1000)
    for(j in thetas_i){
      # données pour la condition 1
      te1=thetaestimate(1,bank,PL,j)
      iterations1=append(iterations1,te1[1])
      estimatedthetas1=append(estimatedthetas1,te1[2])
      # données pour la condition 2
      te2=thetaestimate(2,bank,PL,j)
      iterations2=append(iterations2,te2[1])
      estimatedthetas2=append(estimatedthetas2,te2[2])
      # données pour la condition 3
      te3=thetaestimate(3,bank,PL,j)
      iterations3=append(iterations3,te3[1])
      estimatedthetas3=append(estimatedthetas3,te3[2])
      # données pour la condition 4
      te4=thetaestimate(4,bank,PL,j)
      iterations4=append(iterations4,te4[1])
      estimatedthetas4=append(estimatedthetas4,te4[2])
      # données pour la condition 5
      te5=thetaestimate(5,bank,PL,j)
      iterations5=append(iterations5,te5[1])
      estimatedthetas5=append(estimatedthetas5,te5[2])
      # données pour la condition 6
      te6=thetaestimate(6,bank,PL,j)
      iterations6=append(iterations6,te6[1])
      estimatedthetas6=append(estimatedthetas6,te6[2])
      # données pour la condition 7
      te7=thetaestimate(7,bank,PL,j)
      iterations7=append(iterations7,te7[1])
      estimatedthetas7=append(estimatedthetas7,te7[2])
      # données pour la condition 8
      te8=thetaestimate(8,bank,PL,j)
      iterations8=append(iterations8,te8[1])
      estimatedthetas8=append(estimatedthetas8,te8[2])
      # données pour la condition 9
      te9=thetaestimate(9,bank,PL,j)
      iterations9=append(iterations9,te9[1])
      estimatedthetas9=append(estimatedthetas9,te9[2])
      # données pour la condition 10
      te10=thetaestimate(10,bank,PL,j)
      iterations10=append(iterations10,te10[1])
      estimatedthetas10=append(estimatedthetas10,te10[2])
    }
  }
  c1=c(mean(iterations1),meanbias(truethetas,estimatedthetas1),RMSE(truethetas,estimatedthetas1))
  c2=c(mean(iterations2),meanbias(truethetas,estimatedthetas2),RMSE(truethetas,estimatedthetas2))
  c3=c(mean(iterations3),meanbias(truethetas,estimatedthetas3),RMSE(truethetas,estimatedthetas3))
  c4=c(mean(iterations4),meanbias(truethetas,estimatedthetas4),RMSE(truethetas,estimatedthetas4))
  c5=c(mean(iterations5),meanbias(truethetas,estimatedthetas5),RMSE(truethetas,estimatedthetas5))
  c6=c(mean(iterations6),meanbias(truethetas,estimatedthetas6),RMSE(truethetas,estimatedthetas6))
  c7=c(mean(iterations7),meanbias(truethetas,estimatedthetas7),RMSE(truethetas,estimatedthetas7))
  c8=c(mean(iterations8),meanbias(truethetas,estimatedthetas8),RMSE(truethetas,estimatedthetas8))
  c9=c(mean(iterations9),meanbias(truethetas,estimatedthetas9),RMSE(truethetas,estimatedthetas9))
  c10=c(mean(iterations10),meanbias(truethetas,estimatedthetas10),RMSE(truethetas,estimatedthetas10))
  tableau=data.frame(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,row.names=c("Mean Length","Mean Bias","Mean RMSE"))
  colnames(tableau)=c(1:10)
  tableau
}

########################################## plots des RMSE et biais conditionnels

# plot du RMSE en fonction de la valeur de theta
conditionnalRMSE=function(bank,PL,conditions){
  thetas=seq(-3,3,0.5)
  y=c()
  truethetas=c()
  for(i in thetas){
    truethetas=append(truethetas,rep(i,1000))
  }
  estimatedthetas=c()
  for(i in conditions){
    for(j in truethetas){
      estimatedthetas=append(estimatedthetas,thetaestimate(i,bank,PL,j)[2])
    }
  }
  for(i in 0:(length(conditions)-1)){
    for(j in 0:12){
      y=append(y,RMSE(truethetas[((j*1000)+1):((j+1)*1000)],estimatedthetas[((13000*i)+(j*1000)+1):((13000*i)+(j+1)*1000)]))
    }
  }
  plot(thetas,y[1:13],type="b",col="blue",ylim=c(0,max(y)+0.1),xlab="theta",ylab="RMSE",main="RMSE conditionnel")
  colors=c("red","green","purple","darkgoldenrod1","burlywood3","darkturquoise","black","darksalmon","azure4")
  for(i in 1:(length(conditions)-1)){
    lines(thetas,y[((i*13)+1):((i+1)*13)],type="b",col=colors[i])
  }
  clegend=c("1","2","3","4","5","6","7","8","9","10")
  legend("bottomright",legend=clegend[conditions],col=c(c("blue"),colors[1:(length(conditions)-1)]),lty=rep(1,length(conditions)),cex=0.8)
}

# plot du biais en fonction de la valeur de theta
conditionnalbias=function(bank,PL,conditions){
  thetas=seq(-3,3,0.5)
  y=c()
  truethetas=c()
  for(i in thetas){
    truethetas=append(truethetas,rep(i,1000))
  }
  estimatedthetas=c()
  for(i in conditions){
    for(j in truethetas){
      estimatedthetas=append(estimatedthetas,thetaestimate(i,bank,PL,j)[2])
    }
  }
  for(i in 0:(length(conditions)-1)){
    for(j in 0:12){
      y=append(y,meanbias(truethetas[((j*1000)+1):((j+1)*1000)],estimatedthetas[((13000*i)+(j*1000)+1):((13000*i)+(j+1)*1000)]))
    }
  }
  plot(thetas,y[1:13],type="b",col="blue",ylim=c(min(y)-0.25,max(y)+0.1),xlab="theta",ylab="biais",main="biais conditionnel")
  colors=c("red","green","purple","darkgoldenrod1","burlywood3","darkturquoise","black","darksalmon","azure4")
  for(i in 1:(length(conditions)-1)){
    lines(thetas,y[((i*13)+1):((i+1)*13)],type="b",col=colors[i])
  }
  clegend=c("1","2","3","4","5","6","7","8","9","10")
  legend("bottomright",legend=clegend[conditions],col=c(c("blue"),colors[1:(length(conditions)-1)]),lty=rep(1,length(conditions)),cex=0.8)
}