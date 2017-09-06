f=function(o,e){
  sum((o-e)^2/e)
  
}

ob=c(4,5,6,7)
ex=c(6.25,12.5,15.625,15.625)

f(ob,ex)

##least square
x=c(.9,.76,1.67,1.44,.2,.16,1.12,1.04,0.48,1.33,1.1,1.56,1.15)
y=c(1.8,1.36,2.92,2.61,.42,.49,1.0,2.38,1.24,2.8,2.41,2.8,2.16)

line.fun=function(x,y,n){
  b=(n*sum(x*y)-(sum(x)*sum(y)))/((n*(sum(x^2)))-((sum(x))^2))
  y_bar=mean(y)
  x_bar=mean(x)
  a=y_bar-(b*x_bar)
  c(a,b)
}
line.fun(x,y,13)

##
cluster.fun=function(y,m,N){
  n=length(y)
  ybar=sum(y)/sum(m)
  s2=sum((y-ybar*m)^2)/(n-1)
  Mhat=mean(m)
  varhat=(1-(n/N))*s2/(n*Mhat^2)
  ybar-2*sqrt(varhat)
  ybar+2*sqrt(varhat)
  c(ybar-2*sqrt(varhat),ybar+2*sqrt(varhat))  
}

##estimating popultion total for cluster (w/out M)
cluster.totalM.fun=function(y,m,N){
  n=length(y)
  ybar=mean(y)
  Nybar=N*ybar
  sr2=sum((y-ybar)^2)/(n-1)
  varhatMy=N*N*(1-(n/N))*sr2/n
  c(Nybar,varhatMy,Nybar-2*sqrt(varhatMy),Nybar+2*sqrt(varhatMy))
}
cluster.totalM.fun(cost,saws,N)

###Function cluster stratified
cluster.strat.fun=function(y,m,stratas,N){
  ybars=tapply(y,stratas,mean)
  mbars=tapply(m,stratas,mean)
  y.bar.c=sum(N*ybars)/sum(N*mbars)
  M.est=sum(N*mbars)
  var.terms=y-y.bar.c*m
  s2.c=tapply(var.terms,stratas,var)
  n.st=tapply(y,stratas,length)
  V.ybar.c=(1/M.est^2)*sum(N^2*(1-n.st/N)*(s2.c/n.st))
  c(y.bar.c,2*sqrt(V.ybar.c),y.bar.c-2*sqrt(V.ybar.c),y.bar.c+2*sqrt(V.ybar.c))
}

##Function for estimating pop
pop.est.fun=function(t,s,n){
  N_hat=n*t/s
  var.N_hat=(t^2*n*(n-s))/s^3
  c(N_hat-2*sqrt(var.N_hat), N_hat+2*sqrt(var.N_hat))
}
pop.est.fun(320,91,515)

#N-hat inverse
pop.est.inv.fun=function(t,s,n){
  N_hat=n*t/s
  var.N_hat=(t^2*n*(n-s))/(s^2*(s+1))
  c(N_hat-2*sqrt(var.N_hat), N_hat+2*sqrt(var.N_hat))
}

##Estimator of density lamda
dense.est.fun=function(a,m,n){
  lamda_hat=mean(m)/a
  var_lamda=1/a^2*var(m)/n
  c(lamda_hat-2*sqrt(var_lamda), lamda_hat+2*sqrt(var_lamda))
  
}
hills=c(rep(0,13),rep(1,8),rep(2,12),rep(3,10),rep(4,5),rep(5,2))
dense.est.fun(16,hills,50)

##density of lamda
dense.est2.fun=function(a,m_bar,n){
  lamda_hat=m_bar/a
  var_lamda=lamda_hat/(a*n)
  c(lamda_hat-2*sqrt(var_lamda), lamda_hat+2*sqrt(var_lamda))
}
dense.est2.fun(100,210,15)

##estimate pop denisty and size from stocked quadrants
dense.stocked.fun=function(a,n,y){
  lamda_bar=-(1/a)*log(y/n)
  var_lamda=(1/(n*a^2))*(exp(lamda_bar*a)-1)
  c(lamda_bar-2*sqrt(var_lamda),lamda_bar+2*sqrt(var_lamda))
}
dense.stocked.fun(.5,20,4)

##mean density pre cell
y=c(0,12,1,15,7)
m=c(1,4,114,5,2)

dens.perCell.fun=function(n,y,m,N){
  mu_hat=(1/n)*(sum(y/m))
  s2=var(y/m)
  var.perCell=(1-(n/N))*(s2/n)
  c(mu_hat-2*sqrt(var.perCell),mu_hat+2*sqrt(var.perCell))
  
}
dens.perCell.fun(5,y,m,100)
##random response technique
randomRespon.fun=function(n,theta,n1){
  p=(n1/n)
  p_hat=1/(2*theta-1)*(n1/n)-((1-theta)/(2*theta-1))
  var.p_hat=1/(2*theta-1)^2*(n1/n^2)*(1-(n1/n))
  c(p_hat-2*sqrt(var.p_hat),p_hat+2*sqrt(var.p_hat) )
}

randomRespon.fun(300,.7,105)


ratio.fun=function(x,y,N,n,mu.x){
  r=mean(y)/mean(x)
  mu_y= r*mu.x
  s_squar =sum((y-r*x)^2)/(n-1)
  v =(1-(n/N))*(s_squar/n)
  c(mu_y-2*sqrt(v),mu_y+2*sqrt(v),mu_y)
}