rm(list = ls())
library(Rsolnp)
library(glmnetUtils)
library(truncnorm)
library(HDInterval)
library(invgamma)

# paramater info
T<-100
T0<-20
T1<-T-T0
p<-8 #number of total covariates
N<-40 #number of units
n.f<-3 #number of factors
Sim.size<-1000 #Simulation size
MC.size<-6000 #MCMC size
burnin.size<-1000 #burnin period length

theta0<- -2

bias.TE<-array(0,c(Sim.size,T1,7))
bias.ATE<-matrix(0,Sim.size,7)
MCMC.CP<-rep(0,Sim.size)
MCMC.POWER<-rep(0,Sim.size)
set.seed(2204)


for(MC.goh in 1:Sim.size){
  Y<-matrix(0,N,T)
  Y0<-matrix(0,N,T)
  Z<-matrix(0,N,p)
  
  F<-matrix(0,T+1,n.f)
  F[1,]<-rnorm(n.f,0,1)
  for(t in 2:(T+1)){
    F[t,]<-F[t-1,]*0.2+rnorm(n.f,0,sqrt(.25))
  }
  
  mu.i<-runif(N,-1,1)
  c.t<-matrix(runif(T*p,-.2,.2),T,p)
  c.t[,-(1:2)]<-0
  Z<-matrix(rnorm(p*N,1,sqrt(2)),N,p)
  
  for(i in 1:N){
    b.i<-rnorm(n.f,0,sqrt(0.5))
    
    for(t in 1:T){
      D.it<-1*((i==1)&(t>T0))
      F.t<-F[t+1,]
      Y0[i,t]<-mu.i[i]+sqrt(5*t)+rnorm(1,0,sqrt(0.1))+sum(b.i*F.t)+sum(c.t[t,]*Z[i,])
      Y[i,t]<-theta0*(0.5+sqrt(t/2))*D.it+Y0[i,t]
    }
  }
  
  true.Treat.Effect<-Y[1,-(1:T0)]-Y0[1,-(1:T0)]
  
  X1<-c(Y[1,1:T0])
  X0<-cbind(rep(1,T0),t(Y[-1,1:T0]))
  OLS<-function(b){
    sum((X1-X0%*%b)^2)
  }
  hat.w.ols<-as.matrix(optim(rep(0,dim(X0)[2]),OLS,method="BFGS")$par)
  

#Convex hull by Abadie&Gardeazabal
X1<-c(Y[1,1:T0],Z[1,])
X0<-rbind(cbind(rep(1,T0),t(Y[-1,1:T0])),cbind(rep(0,p),t(Z[-1,])))
LS<-function(w){
  sum((X1-as.numeric(X0[,-1]%*%w))^2)
}
sum.constrant<-function(w){
  sum(w)
}
hat.w.scm1<-rep(0,N)
invisible(capture.output(hat.w.scm1[-1]<-solnp(rep(1/(N-1),N-1),fun=LS,eqfun=sum.constrant,eqB=1,LB=rep(0,N-1))$pars))

  

#Conical hull by Li
LS1<-function(w){
  sum((X1-as.numeric(X0%*%w))^2)
}
  
invisible(capture.output(hat.w.scm2<-solnp(rep(1/(N),N),fun=LS1,LB=c(-Inf,rep(0,N-1)))$pars))
  
#Enet by Doudchenko and Imbens 
  
X1<-c(Y[1,1:T0])
X0<-cbind(rep(1,T0),t(Y[-1,1:T0]))
  
cv.enet<-cva.glmnet(X0[,-1],X1,alpha=seq(0,1,0.1)) #Performs elastic net CV for alpha and lambda simultaneously
cvm.alpha<-rep(0,11)
for(i in 1:11){
  cvm.alpha[i]<-min((cv.enet$modlist[[i]])$cvm)
}
  
alpha.best<-which.min(cvm.alpha)
lambda.best<-(cv.enet$modlist[[alpha.best]])$lambda.min
hat.w.enet<-as.numeric(coef(glmnet(X0[,-1],X1,alpha=(seq(0,1,0.1)[alpha.best]),lambda=lambda.best)))


hat.Y.1t.ols<-as.numeric(t(Y[-1,-(1:T0)])%*%hat.w.ols[-1])+rep(hat.w.ols[1],T1)
hat.Y.1t.scm1<-as.numeric(t(Y[-1,-(1:T0)])%*%hat.w.scm1[-1])+rep(hat.w.scm1[1],T1)
hat.Y.1t.scm2<-as.numeric(t(Y[-1,-(1:T0)])%*%hat.w.scm2[-1])+rep(hat.w.scm2[1],T1)
hat.Y.1t.enet<-as.numeric(t(Y[-1,-(1:T0)])%*%hat.w.enet[-1])+rep(hat.w.enet[1],T1)

  
  
################################
#Bayes MAP
################################
#EM algorithm
################################
X1<-c(Y[1,1:T0],Z[1,])
X0<-rbind(cbind(rep(1,T0),t(Y[-1,1:T0])),cbind(rep(0,p),t(Z[-1,])))
  
c0<-1/2   
d0<-1/2 
inv.nu<-1
prob0<-0.5
hat.xi<-rep(1,p)
xi.all<-c(rep(1,T0),hat.xi)
  
#Naive Bayes SCM
  
LS2<-function(w){
  inv.nu*sum(xi.all*as.numeric((X1-as.numeric(X0%*%w))^2))
}
  
sum.constrant2<-function(w){
  sum(w[-1])
}
invisible(capture.output(hat.w.MAP.naive<-solnp(rep(1/(N),N),fun=LS2,eqfun=sum.constrant2,eqB=1,LB=c(-Inf,rep(0,N-1)))$pars))
hat.Y.1t.MAP.naive<-as.numeric(t(Y[-1,-(1:T0)])%*%hat.w.MAP.naive[-1])+rep(hat.w.MAP.naive[1],T1)
  

#Bayes SCM
  
hat.w.MAP<-rep(0,N)
rep.r<-1000
for(r in 1:rep.r){
  #print(r)
  hat.w.MAP0<-hat.w.MAP
  invisible(capture.output(hat.w.MAP<-solnp(rep(1/(N),N),fun=LS2,eqfun=sum.constrant2,eqB=1,LB=c(-Inf,rep(0,N-1)))$pars))
  if(max(abs(hat.w.MAP0-hat.w.MAP))<(0.1)^3){break}
    
#MC-E-Step
Gibbs.rep<-17000
MC.nu<-rep(0,Gibbs.rep)
MC.xi<-matrix(0,Gibbs.rep,T0+p)
for(i in 1:Gibbs.rep){
    nu.a<-(T0+sum(hat.xi))/2+c0
    #nu.a<-(T0)/2+c0
    nu.b<-0.5*sum(xi.all*(as.numeric(X1-X0%*%hat.w.MAP)^2))+d0
    hat.nu<-rinvgamma(1,shape=nu.a,rate=nu.b)
    phi.Z<-dnorm(X1[-(1:T0)],as.numeric(X0[-(1:T0),]%*%hat.w.MAP),sqrt(hat.nu))
    xi.prob<-prob0/((1-prob0)/phi.Z+prob0)
    hat.xi<-rbinom(p,1,xi.prob)
    xi.all<-c(rep(1,T0),hat.xi)
    #prob0<-rbeta(1,(MCMC.xi+1),(p-MCMC.xi+1))
    MC.xi[i,]<-xi.all
    MC.nu[i]<-hat.nu
  }
    xi.all<-apply(MC.xi[-(1:2000),],2,mean)
    inv.nu<-mean(1/MC.nu[-(1:2000)])
    if(r==rep.r){print("Algorithm does NOT reach the convergence")}
  }
  
##########################################
hat.Y.1t.MAP<-as.numeric(t(Y[-1,-(1:T0)])%*%hat.w.MAP[-1])+rep(hat.w.MAP[1],T1)
#Gibbs sampler

  MC.size<-15000 #MCMC size
  burnin.size<-5000 #burnin period length
  
  MCMC.omega<-hat.w.MAP
  MCMC.nu<-1/inv.nu
  MCMC.xi<-1*(xi.all[-(1:T0)]>.5)
  lb<-rep(0,N)
  lb[1]<--Inf
  MCMC.hat.Y.Bayes<-matrix(0,MC.size,T1)
  MCMC.hat.Y.Bayes.v2<-matrix(0,MC.size,T1)
  active.donor<-which(abs(hat.w.MAP)>0.1^5)
  MCMC.omega[-active.donor]<-0
  MCMC.omega[(active.donor[-1])]<-hat.w.MAP[(active.donor[-1])]/sum(hat.w.MAP[(active.donor[-1])])
  M<-max(active.donor)
  D1<-dim(X0)[1]
  D2<-dim(X0)[2]
  MCMC.seq.nu<-rep(0,MC.size)
  MCMC.seq.omega<-matrix(0,MC.size,N)
  
  for(MCMC.goh in 1:MC.size){
    sig2.i<-1/sum((inv.nu*xi.all)*(X0[,1]^2))
    mu.i<-sum(X0[,1]*((inv.nu*xi.all)*(X1-as.numeric(X0[,-1]%*%MCMC.omega[-1]))))*sig2.i
    MCMC.omega[1]<-rnorm(1,mean=mu.i,sd=sqrt(sig2.i))
    
    for(i in active.donor[-1]){
      if(i==M){
        MCMC.omega[i]<-1-sum(MCMC.omega[-c(1,i)])
      }else{
        #sig2.i<-1/sum((inv.nu*xi.all)*(X0[,i]^2))
        #mu.i<-sum((X0[,i] *(inv.nu*xi.all)* (X1-as.numeric(X0[,-i]%*%MCMC.omega[-i]))) )*sig2.i
        X.i.s<-X0[,i]-X0[,M]
        sig2.i<-1/sum((inv.nu*xi.all)*(X.i.s^2))
        mu.i<-sum(X.i.s*((inv.nu*xi.all)*(X1-X0[,M]-MCMC.omega[1]*X0[,1]-as.numeric((X0[,-c(1,i,M)]-matrix(X0[,M],D1,(D2-3)) )%*%MCMC.omega[-c(1,i,M)]))))*sig2.i
        UB.i<-1-sum(MCMC.omega[-c(1,i,M)])
        MCMC.omega[i]<-rtruncnorm(1,a=0,b=UB.i,mean =mu.i,sd=sqrt(sig2.i))}
    }
    
    MCMC.seq.omega[MCMC.goh,]<-MCMC.omega
    #3. xi from Bernoulli
    phi.Z<-dnorm(X1[-(1:T0)],as.numeric(X0[-(1:T0),]%*%MCMC.omega),sqrt(MCMC.nu))
    xi.prob<-prob0/((1-prob0)/phi.Z+prob0)
    MCMC.xi<-rbinom(p,1,xi.prob)
    xi.all<-c(rep(1,T0),MCMC.xi)
    ##############################
    
    #2. nu from Inverse-Gamma
    nu.a<-(T0+sum(MCMC.xi))/2+c0
    nu.b<-0.5*sum(xi.all*(as.numeric(X1-X0%*%MCMC.omega)^2))+d0
    
    MCMC.nu<-rinvgamma(1,shape=nu.a,rate=nu.b)
    inv.nu<-1/MCMC.nu #fixed
    MCMC.seq.nu[MCMC.goh]<-MCMC.nu
    MCMC.hat.Y.Bayes[MCMC.goh,]<-as.numeric(t(Y[-1,-(1:T0)])%*%MCMC.omega[-1])+rep(MCMC.omega[1],T1)
    MCMC.hat.Y.Bayes.v2[MCMC.goh,]<-rnorm(T1,MCMC.hat.Y.Bayes[MCMC.goh,],sqrt(MCMC.nu))
  }
  
######################################################
hat.Y.1t.Bayes<-apply(MCMC.hat.Y.Bayes.v2[-(1:burnin.size),],2,mean)
  

#For Treatment effect
#HPD interval #
  
hat.Y.1t.Bayes.HPD95<-matrix(0,T1,2)
for(t in 1:T1){
  hat.Y.1t.Bayes.HPD95[t,]<-(hdi(MCMC.hat.Y.Bayes.v2[-(1:burnin.size),t],credMass=0.95))[1:2]
}
  
hat.Y.1t.Bayes.UB95<-hat.Y.1t.Bayes.HPD95[,2]
hat.Y.1t.Bayes.LB95<-hat.Y.1t.Bayes.HPD95[,1]
  
  
#For Average Treatment Effect
MCMC.theta<-matrix(Y[1,(T0+1):T],dim(MCMC.hat.Y.Bayes)[1],dim(MCMC.hat.Y.Bayes)[2],byrow=TRUE)-MCMC.hat.Y.Bayes.v2
MAP.ATE<-mean(Y[1,(T0+1):T]-hat.Y.1t.MAP)
  
MCMC.ATE<-apply(MCMC.theta[-(1:burnin.size),],1,mean)
  
#Plot 1#
MCMC.ATE.all<-apply(MCMC.theta,1,mean)
  if(MC.goh==1){
    par(mfrow=c(1,2))
    plot(MCMC.ATE.all[1:5000],type='l',xlab="iterations",ylab="ATE",ylim=c(MAP.ATE-0.4,MAP.ATE+0.4))
    abline(h=MAP.ATE,lwd=1,col=2)
    legend("bottomright","Bayes SCM estimate of ATE",col=2,lwd=1)
    hist(MCMC.ATE,breaks="scott",xlab="ATE",main="",prob=TRUE)
  }
  true.ATE<-mean(true.Treat.Effect)
  
  
  Bayes.HPD.95<-(hdi(MCMC.ATE,credMass=0.95))[1:2]
  Bays.CI.LB<-Bayes.HPD.95[1]
  Bays.CI.UB<-Bayes.HPD.95[2]
  if((Bays.CI.LB<=true.ATE)&(Bays.CI.UB>=true.ATE)){
    MCMC.CP[MC.goh]<-1
  }
  if((Bays.CI.LB>true.ATE)|(Bays.CI.UB<true.ATE)){
    MCMC.POWER[MC.goh]<-1
  }
  
  est.Treat.Effect.ols<-Y[1,(T0+1):T]-hat.Y.1t.ols
  est.Treat.Effect.scm1<-Y[1,(T0+1):T]-hat.Y.1t.scm1
  est.Treat.Effect.scm2<-Y[1,(T0+1):T]-hat.Y.1t.scm2
  est.Treat.Effect.enet<-Y[1,(T0+1):T]-hat.Y.1t.enet
  est.Treat.Effect.bayes<-Y[1,(T0+1):T]-hat.Y.1t.MAP
  est.Treat.Effect.bayes.naive<-Y[1,(T0+1):T]-hat.Y.1t.MAP.naive
  est.Treat.Effect.bayes.mean<-Y[1,(T0+1):T]-hat.Y.1t.Bayes
  
  bias.TE[MC.goh,,1]<-est.Treat.Effect.ols-true.Treat.Effect
  bias.TE[MC.goh,,2]<-est.Treat.Effect.scm1-true.Treat.Effect
  bias.TE[MC.goh,,3]<-est.Treat.Effect.scm2-true.Treat.Effect
  bias.TE[MC.goh,,4]<-est.Treat.Effect.enet-true.Treat.Effect
  bias.TE[MC.goh,,5]<-est.Treat.Effect.bayes-true.Treat.Effect
  bias.TE[MC.goh,,6]<-est.Treat.Effect.bayes.naive-true.Treat.Effect
  bias.TE[MC.goh,,7]<-est.Treat.Effect.bayes.mean-true.Treat.Effect
  
  bias.ATE[MC.goh,1]<-mean(est.Treat.Effect.ols)-mean(true.Treat.Effect)
  bias.ATE[MC.goh,2]<-mean(est.Treat.Effect.scm1)-mean(true.Treat.Effect)
  bias.ATE[MC.goh,3]<-mean(est.Treat.Effect.scm2)-mean(true.Treat.Effect)
  bias.ATE[MC.goh,4]<-mean(est.Treat.Effect.enet)-mean(true.Treat.Effect)
  bias.ATE[MC.goh,5]<-mean(est.Treat.Effect.bayes)-mean(true.Treat.Effect)
  bias.ATE[MC.goh,6]<-mean(est.Treat.Effect.bayes.naive)-mean(true.Treat.Effect)
  bias.ATE[MC.goh,7]<-mean(est.Treat.Effect.bayes.mean)-mean(true.Treat.Effect)
  
  print(paste(MC.goh,"-th simulation is running"))
  
}

table1.theta0_1<-cbind(c(
  mean(bias.TE[,,2]^2),
  mean(bias.TE[,,4]^2),
  mean(bias.TE[,,3]^2),
  mean(bias.TE[,,6]^2),
  mean(bias.TE[,,5]^2)),
  apply(bias.ATE[,c(2,4,3,6,5)]^2,2,mean)
)

#Table 1
colnames(table1.theta0_1)<-c("MSE_TE","MSE_ATE")
rownames(table1.theta0_1)<-c("ADH-SCM","DI-SCM","L-SCM","naive Bayes SCM","Bayes SCM")
round(table1.theta0_1,4)
#Table 2
print(mean(MCMC.CP))
