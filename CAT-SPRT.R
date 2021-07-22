############################################## construction de la banque d'items
nbank5=250
a5=rlnorm(nbank5, meanlog = 1.128171, sdlog = 0.3)
b5=runif(nbank5,-4,4)

bank5=matrix(c(a5,b5),nrow=2,ncol=nbank5,byrow=T)

############################################################################ 2PL

# test sans compter les cas où theta tombe dans la région d'indifférence
estimatedareaSPRT_2PL_1=function(bank,selection,cpoint,delta,alpha,beta,theta){
  theta1=cpoint-delta
  theta2=cpoint+delta
  if(theta<=theta1 | theta>=theta2){
    x_0=c()
    i_0=sample(which(bank[2,]>=-1.602 & bank[2,]<=-0.314),3,replace=FALSE)
    for(i in i_0){
      if(runif(1,0,1)<=presponse(bank,2,i,theta)){
        x_0=append(x_0,1)
      }
      else{
        x_0=append(x_0,0)
      }
    }
    C=sum(log((1+exp(bank[1,i_0]*(theta1-bank[2,i_0])))/(1+exp(bank[1,i_0]*(theta2-bank[2,i_0])))))
    C_0=C
    while(((log(beta/(1-alpha))-C)/(2*delta)<sum(bank[1,i_0]*x_0)) & (sum(bank[1,i_0]*x_0)<(log((1-beta)/alpha)-C)/(2*delta)) & length(x_0)<40){
      if(selection==1){
        theta_0=optimize(llh(bank,2,x_0,i_0),c(-4,4),maximum=TRUE)$maximum
        i=which.max(replace(fisherinfovect(bank,2,theta_0),i_0,0))
        i_0=append(i_0,i)
        if(runif(1,0,1)<=presponse(bank,2,i,theta)){
          x_0=append(x_0,1)
        }
        else{
          x_0=append(x_0,0)
        }
      }
      if(selection==2){
        i=which.max(replace(fisherinfovect(bank,2,cpoint),i_0,0))
        i_0=append(i_0,i)
        if(runif(1,0,1)<=presponse(bank,2,i,theta)){
          x_0=append(x_0,1)
        }
        else{
          x_0=append(x_0,0)
        }
      }
      if(selection==3){
        theta_0=optimize(llh(bank,2,x_0,i_0),c(-4,4),maximum=TRUE)$maximum
        j=which.min(abs(theta1-theta_0),abs(theta2-theta_0))
        c=c(theta1,theta2)
        i=which.max(replace(fisherinfovect(bank,2,c[j]),i_0,0))
        i_0=append(i_0,i)
        if(runif(1,0,1)<=presponse(bank,2,i,theta)){
          x_0=append(x_0,1)
        }
        else{
          x_0=append(x_0,0)
        }
      }
      C=sum(log((1+exp(bank[1,i_0]*(theta1-bank[2,i_0])))/(1+exp(bank[1,i_0]*(theta2-bank[2,i_0])))))
      C_0=append(C_0,C)
    }
    if(sum(bank[1,i_0]*x_0)<=(log(beta/(1-alpha))-C)/(2*delta)){
      c(length(x_0),(theta<=theta1))
    }
    else if(((log((1-beta)/alpha)-C)/(2*delta))<=sum(bank[1,i_0]*x_0)){
      c(length(x_0),(theta>=theta2))
    }
  }
}

plotSPRT_2PL_1=function(bank,selection,cpoint,delta,alpha,beta,theta){
  theta1=cpoint-delta
  theta2=cpoint+delta
  x_0=c()
  i_0=sample(which(bank[2,]>=-1.602 & bank[2,]<=-0.314),3,replace=FALSE)
  for(i in i_0){
    if(runif(1,0,1)<=presponse(bank,2,i,theta)){
      x_0=append(x_0,1)
    }
    else{
      x_0=append(x_0,0)
    }
  }
  C=sum(log((1+exp(bank[1,i_0]*(theta1-bank[2,i_0])))/(1+exp(bank[1,i_0]*(theta2-bank[2,i_0])))))
  C_0=C
  while(((log(beta/(1-alpha))-C)/(2*delta)<sum(bank[1,i_0]*x_0)) & (sum(bank[1,i_0]*x_0)<(log((1-beta)/alpha)-C)/(2*delta)) & length(x_0)<40){
    if(selection==1){
      theta_0=optimize(llh(bank,2,x_0,i_0),c(-4,4),maximum=TRUE)$maximum
      i=which.max(replace(fisherinfovect(bank,2,theta_0),i_0,0))
      i_0=append(i_0,i)
      if(runif(1,0,1)<=presponse(bank,2,i,theta)){
        x_0=append(x_0,1)
      }
      else{
        x_0=append(x_0,0)
      }
    }
    if(selection==2){
      i=which.max(replace(fisherinfovect(bank,2,cpoint),i_0,0))
      i_0=append(i_0,i)
      if(runif(1,0,1)<=presponse(bank,2,i,theta)){
        x_0=append(x_0,1)
      }
      else{
        x_0=append(x_0,0)
      }
    }
    if(selection==3){
      theta_0=optimize(llh(bank,2,x_0,i_0),c(-4,4),maximum=TRUE)$maximum
      j=which.min(abs(theta1-theta_0),abs(theta2-theta_0))
      c=c(theta1,theta2)
      i=which.max(replace(fisherinfovect(bank,2,c[j]),i_0,0))
      i_0=append(i_0,i)
      if(runif(1,0,1)<=presponse(bank,2,i,theta)){
        x_0=append(x_0,1)
      }
      else{
        x_0=append(x_0,0)
      }
    }
    C=sum(log((1+exp(bank[1,i_0]*(theta1-bank[2,i_0])))/(1+exp(bank[1,i_0]*(theta2-bank[2,i_0])))))
    C_0=append(C_0,C)
  }
  x=3:length(x_0)
  y1=(log(beta/(1-alpha))-C_0)/(2*delta)
  y2=(log((1-beta)/alpha)-C_0)/(2*delta)
  y=cumsum(bank[1,i_0]*x_0)[-c(1,2)]
  plot(x,y,type="b",col="blue",ylim=c(min(c(y1,y2,y)),max(c(y1,y2,y))))
  lines(x,y1,type="b")
  lines(x,y2,type="b")
}

# test où on compare le rapport de vraisemblance à 1 si theta tombe dans la région d'indifférence
estimatedareaSPRT_2PL_2=function(bank,selection,cpoint,delta,alpha,beta,theta){
  theta1=cpoint-delta
  theta2=cpoint+delta
  x_0=c()
  i_0=sample(which(bank[2,]>=-1.602 & bank[2,]<=-0.314),3,replace=FALSE)
  for(i in i_0){
    if(runif(1,0,1)<=presponse(bank,2,i,theta)){
      x_0=append(x_0,1)
    }
    else{
      x_0=append(x_0,0)
    }
  }
  C=sum(log((1+exp(bank[1,i_0]*(theta1-bank[2,i_0])))/(1+exp(bank[1,i_0]*(theta2-bank[2,i_0])))))
  while(((log(beta/(1-alpha))-C)/(2*delta)<sum(bank[1,i_0]*x_0)) & (sum(bank[1,i_0]*x_0)<(log((1-beta)/alpha)-C)/(2*delta)) & length(x_0)<40){
    if(selection==1){
      theta_0=optimize(llh(bank,2,x_0,i_0),c(-4,4),maximum=TRUE)$maximum
      i=which.max(replace(fisherinfovect(bank,2,theta_0),i_0,0))
      i_0=append(i_0,i)
      if(runif(1,0,1)<=presponse(bank,2,i,theta)){
        x_0=append(x_0,1)
      }
      else{
        x_0=append(x_0,0)
      }
    }
    if(selection==2){
      i=which.max(replace(fisherinfovect(bank,2,cpoint),i_0,0))
      i_0=append(i_0,i)
      if(runif(1,0,1)<=presponse(bank,2,i,theta)){
        x_0=append(x_0,1)
      }
      else{
        x_0=append(x_0,0)
      }
    }
    if(selection==3){
      theta_0=optimize(llh(bank,2,x_0,i_0),c(-4,4),maximum=TRUE)$maximum
      j=which.min(abs(theta1-theta_0),abs(theta2-theta_0))
      c=c(theta1,theta2)
      i=which.max(replace(fisherinfovect(bank,2,c[j]),i_0,0))
      i_0=append(i_0,i)
      if(runif(1,0,1)<=presponse(bank,2,i,theta)){
        x_0=append(x_0,1)
      }
      else{
        x_0=append(x_0,0)
      }
    }
    C=sum(log((1+exp(bank[1,i_0]*(theta1-bank[2,i_0])))/(1+exp(bank[1,i_0]*(theta2-bank[2,i_0])))))
  }
  if(sum(bank[1,i_0]*x_0)<=(log(beta/(1-alpha))-C)/(2*delta)){
    c(length(x_0),(theta<=cpoint))
  }
  else if(((log((1-beta)/alpha)-C)/(2*delta))<=sum(bank[1,i_0]*x_0)){
    c(length(x_0),(theta>=cpoint))
  }
  else{
    if(prod(((presponse(bank,2,i_0,theta2)/presponse(bank,2,i_0,theta1))^x_0)*(((1-presponse(bank,2,i_0,theta2))/(1-presponse(bank,2,i_0,theta1)))^(1-x_0)))>=1){
      c(length(x_0),(theta>=cpoint))
    }
    else{
      c(length(x_0),(theta<=cpoint))
    }
  }
}

c=c()
m=c()
for(theta in rep(0.25,20000)){
  m_c=estimatedareaSPRT_2PL_1(bank5,2,0.1,0.15,0.05,0.05,theta)
  c=append(c,m_c[2])
  m=append(m,m_c[1])
}

(length(c[c==1])/length(c))*100
mean(m)