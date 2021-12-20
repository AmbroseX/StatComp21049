## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(StatComp21049)

## -----------------------------------------------------------------------------
plot(1:5, 1:5, xlim = c(0,6), ylim = c (0,6), type = "n")
text(x = 3, y = 3, labels = "R for beginner")

## -----------------------------------------------------------------------------
plot(1:5, 1:5, xlim = c(0,6), ylim = c (0,6), type = "n")
text(x = c(3, 3), y = c(5,3), labels = c("R for beginner", "It's easy to learn"))

## -----------------------------------------------------------------------------
x<- c(1:6)
y<- c(1:6)
plot(x,y,main="test",type="p",xlim = c(0,7),ylim = c(0,7))

## -----------------------------------------------------------------------------
x<- c(1:6)
y<- c(1:6)
plot(x,y,main="test",type="l",xlim = c(0,7),ylim = c(0,7))

## -----------------------------------------------------------------------------
x<- c(1:6)
y<- c(1:6)
plot(x,y,main="test",type="b",xlim = c(0,7),ylim = c(0,7))

## -----------------------------------------------------------------------------
x<- c(1:6)
y<- c(1:6)
plot(x,y,main="test",type="c",xlim = c(0,7),ylim = c(0,7))

## -----------------------------------------------------------------------------
x<- c(1:6)
y<- c(1:6)
plot(x,y,main="test",type="o",xlim = c(0,7),ylim = c(0,7))

## -----------------------------------------------------------------------------
x<- c(1:6)
y<- c(1:6)
plot(x,y,main="test",type="h",xlim = c(0,7),ylim = c(0,7))

## -----------------------------------------------------------------------------
x<- c(1:6)
y<- c(1:6)
plot(x,y,main="test",type="s",xlim = c(0,7),ylim = c(0,7))

## -----------------------------------------------------------------------------
x<- c(1:6)
y<- c(1:6)
plot(x,y,main="test",type="S",xlim = c(0,7),ylim = c(0,7))

## -----------------------------------------------------------------------------
x<- c(1:6)
y<- c(1:6)
plot(x,y,main="test",type="n",xlim = c(0,7),ylim = c(0,7))

## ----pressure,echo=FALSE------------------------------------------------------
x<-c(1,2,3)
y<-c(2,3,5)
table(x,y)

## -----------------------------------------------------------------------------
rm(list=ls())
0.1->sigma;1e5->n;k<-1;numeric(n)->x
while (k<n) {
  runif(1)->u
  x[k]<-sqrt(2*log(1/u))*sigma; #广义均匀分布和单位均匀分布
  k+1->k
}
hist(x,prob=TRUE,main=expression(f(y)==x/(sigma^2)*exp(-x^2/(2*sigma^2))))
t<-seq(0,0.5,0.01)
lines(t,t/(sigma^2)*exp(-t^2/(2*sigma^2)))

## -----------------------------------------------------------------------------
rm(list=ls())
1->sigma;1e5->n;k<-1;numeric(n)->x
while (k<n) {
  runif(1)->u
  x[k]<-sqrt(2*log(1/u))*sigma; #广义均匀分布和单位均匀分布
  k <- k+1
}
hist(x,prob=TRUE,main=expression(f(y)==x/(sigma^2)*exp(-x^2/(2*sigma^2))))
t<-seq(0,6,0.01)
lines(t,t/(sigma^2)*exp(-t^2/(2*sigma^2)))

## -----------------------------------------------------------------------------
rm(list=ls())
10->sigma;1e5->n;k<-1;numeric(n)->x
while (k<n) {
  runif(1)->u
  x[k]<-sqrt(2*log(1/u))*sigma; #广义均匀分布和单位均匀分布
  k <- k+1
}
hist(x,prob=TRUE,main=expression(f(y)==x/(sigma^2)*exp(-x^2/(2*sigma^2)),sigma==1))
t<-seq(0,40,0.01)
lines(t,t/(sigma^2)*exp(-t^2/(2*sigma^2)))

## -----------------------------------------------------------------------------
rm(list=ls())
p<-0.75
n<-1e3;
x1<-rnorm(n,0,1)
x2<-rnorm(n,3,1)
x<-numeric(n)
k=1
while(k<n){
  u<-runif(1)
  if (u<=p) x[k]<-x1[k] else x[k]<-x2[k]
  k+1->k
}
hist(x,prob=TRUE,main=expression(f(x)==p*x1+(1-p)*x2))


## -----------------------------------------------------------------------------
rm(list=ls())
p<-0.6
n<-1e3;
x1<-rnorm(n,0,1)
x2<-rnorm(n,3,1)
x<-numeric(n)
k=1
while(k<n){
  u<-runif(1)
  if (u<=p) x[k]<-x1[k] else x[k]<-x2[k]
  k+1->k
}
hist(x,prob=TRUE,main=expression(f(x)==p*x1+(1-p)*x2))


## -----------------------------------------------------------------------------
rm(list=ls())
p<-0.55
n<-1e3
x1<-rnorm(n,0,1)
x2<-rnorm(n,3,1)
x<-numeric(n)
k=1
while(k<n){
  u<-runif(1)
  if (u<=p) x[k]<-x1[k] else x[k]<-x2[k]
  k+1->k
}
hist(x,prob=TRUE,main=expression(f(x)==p*x1+(1-p)*x2))


## -----------------------------------------------------------------------------
rm(list=ls())
p<-0.5
n<-1e4;
x1<-rnorm(n,0,1)
x2<-rnorm(n,3,1)
x<-numeric(n)
k=1
while(k<n){
  u<-runif(1)
  if (u<=p) x[k]<-x1[k] else x[k]<-x2[k]
  k <- k+1
}
hist(x,prob=TRUE,main=expression(f(x)==p*x1+(1-p)*x2))


## -----------------------------------------------------------------------------
rm(list=ls())
p<-0.45
n<-1e3;
x1<-rnorm(n,0,1)
x2<-rnorm(n,3,1)
x<-numeric(n)
k=1
while(k<n){
  u<-runif(1)
  if (u<=p) x[k]<-x1[k] else x[k]<-x2[k]
  k <- k+1
}
hist(x,prob=TRUE,main=expression(f(x)==p*x1+(1-p)*x2))


## -----------------------------------------------------------------------------
rm(list=ls())
p<-0.4
n<-1e3;
x1<-rnorm(n,0,1)
x2<-rnorm(n,3,1)
x<-numeric(n)
k=1
while(k<n){
  u<-runif(1)
  if (u<=p) x[k]<-x1[k] else x[k]<-x2[k]
  k <- k+1
}
hist(x,prob=TRUE,main=expression(f(x)==p*x1+(1-p)*x2))


## -----------------------------------------------------------------------------
rm(list=ls())
p<-0.35
n<-1e3;
x1<-rnorm(n,0,1)
x2<-rnorm(n,3,1)
x<-numeric(n)
k=1
while(k<n){
  u<-runif(1)
  if (u<=p) x[k]<-x1[k] else x[k]<-x2[k]
  k <- k+1
}
hist(x,prob=TRUE,main=expression(f(x)==p*x1+(1-p)*x2))

## -----------------------------------------------------------------------------
rm(list=ls())
lambda.poison<-2
lambda.gamma<-1
a = 1
n<-1e3
Npoison<-numeric(n)
k<-1
while(k<=n){
  Npoison[k]<-rpois(1,lambda.poison*k)
  k <- k+1
}
maxp<- max(Npoison)
Y<-rgamma(maxp,shape = a,scale =1/lambda.gamma)
x<-numeric(n)
k<-1
while(k<=n){
  x[k]<-sum(Y[1:Npoison[k]])
  k+1->k
}
plot(x)

## -----------------------------------------------------------------------------
rm(list=ls())
lambda.poison<-2
lambda.gamma<-1
a = 1
t=10
n<-1e4
k<-1
x<-numeric(n)
while(k<=n){
  Npoison<-rpois(1,lambda.poison*t)
  Y<-rgamma(Npoison,shape = a,scale =1/lambda.gamma)
  x[k]<-sum(Y[1:Npoison])
  k+1->k
}
paste('E[X(10)]=',mean(x))
paste('Var(X(10)=',var(x))

## -----------------------------------------------------------------------------
rm(list=ls())
lambda.poison<-2
lambda.gamma<-1
a = 1
Ex10<-10*lambda.poison*a/lambda.gamma
Var10<-10*lambda.poison*(a+a^2)/lambda.gamma^2
paste("E[X(10)]=",Ex10,",Var[X(10)]=",Var10)

## -----------------------------------------------------------------------------
rm(list=ls())
lambda.poison<-2
lambda.gamma<-2
a = 1
t=10
n<-1e4
k<-1
x<-numeric(n)
while(k<=n){
  Npoison<-rpois(1,lambda.poison*t)
  Y<-rgamma(Npoison,shape = a,scale =1/lambda.gamma)
  x[k]<-sum(Y[1:Npoison])
  k+1->k
}
paste('E[X(10)]=',mean(x))
paste('Var(X(10)=',var(x))

## -----------------------------------------------------------------------------
rm(list=ls())
lambda.poison<-2
lambda.gamma<-2
a = 1
Ex10<-10*lambda.poison*a/lambda.gamma
Var10<-10*lambda.poison*(a+a^2)/lambda.gamma^2
paste("E[X(10)]=",Ex10,",Var[X(10)]=",Var10)

## -----------------------------------------------------------------------------
rm(list=ls())
lambda.poison<-5
lambda.gamma<-4
a = 2
t=10
n<-1e4
k<-1
x<-numeric(n)
while(k<=n){
  Npoison<-rpois(1,lambda.poison*t)
  Y<-rgamma(Npoison,shape = a,scale =1/lambda.gamma)
  x[k]<-sum(Y[1:Npoison])
  k+1->k
}
paste('E[X(10)]=',mean(x),'Var(X(10)=',var(x))

## -----------------------------------------------------------------------------
rm(list=ls())
lambda.poison<-5
lambda.gamma<-4
a = 2
Ex10<-10*lambda.poison*a/lambda.gamma
Var10<-10*lambda.poison*(a+a^2)/lambda.gamma^2
paste("E[X(10)]=",Ex10,",Var[X(10)]=",Var10)

## -----------------------------------------------------------------------------
m<-1e4;x<-runif(m,min=0,max=1)
theta.hat <- mean(x^2*(1-x)^2)
result = matrix(0,1,2)
result[1,1] = theta.hat
result[1,2] = beta(3,3)
rownames(result) = c("Beta(3,3)")
colnames(result) = c("Estimate", "Theoretical")
result

## -----------------------------------------------------------------------------
rm(list=ls())

pbetaCDF <- function(x,a,b) {
  pbetax<-numeric(length(x))
  m<-1e6;
  y<-runif(m,min=0,max=1)
  Beta_a_b <- mean(y^(a-1)*(1-y)^(b-1))
  for(i in 1:length(x)){
    n<-1e6;
    y<-runif(n,min=0,max = x[i])
    Beta_x <- x[i]*mean(y^(a-1)*(1-y)^(b-1))
    pbetax[i] <- Beta_x/Beta_a_b
    }
  return(pbetax)
}
x<-seq(0.1,1,0.1)
y<-pbetaCDF(x,3,3)
print(y)

## -----------------------------------------------------------------------------
x<-seq(0.1,1,0.1)
result = matrix(0,3,length(x))
result[1,1:length(x)] = x
result[2,1:length(x)] = pbetaCDF(x,3,3)
result[3,1:length(x)] = pbeta(x,3,3)
rownames(result) = c("x","pbetaCDF(3,3)","pbeta(x,3,3)")
colnames(result) = c(x)
result

## -----------------------------------------------------------------------------
rm(list=ls())
m<-1e5
y<-runif(m,min=0,max=1)
x1<-sqrt(1-log(y))
x2<-1-log(y)
f1<-x1/(2*exp(1)*sqrt(2*pi))*exp((x1^2)/2)
f2<-x2^2/(sqrt(2*pi))*exp(x2-(x2^2)/2-1)
result = matrix(0,2,2)
result[1,1] = mean(f1)
result[1,2] = mean(f1)
result[2,1] = sd(f1)
result[2,2] = sd(f2)
rownames(result) = c("est","sd")
colnames(result) = c("f1(x)","f2(x)")
result

## -----------------------------------------------------------------------------
rm(list = ls())
m <- 20
n <- 2
alpha <- 0.05
Y <- rchisq(m,n)
a<-t.test(Y,conf.level = 1-alpha)
paste('均值的 95%t-interval:',a[4])
UCL <- (m-1)*var(Y)/qchisq(alpha,df=m-1)
paste('upper confidence limits UCL =',UCL)

## -----------------------------------------------------------------------------
rm(list= ls())
m <- 20
n <- 2
alpha <- .05
mu0 = n
x <- rchisq(m, df = n)
LCL <- mu0+var(x)*qt(alpha/2, df = m-1)/sqrt(m)
UCL <- mu0+var(x)*qt(1-alpha/2, df = m-1)/sqrt(m)
paste('均值的 95%t-interval:[',LCL,',',UCL,']')

## -----------------------------------------------------------------------------
#计算欧氏距离
Eucdidean <- function(x){
    return(sqrt(sum(x*x)))
}

#计算第k项
fk <- function(a,d,k){
  lgfk <- (2*k+2)*log(Eucdidean(a))+lgamma((d+1)/2)+lgamma(k+3/2)-lgamma(k+d/2+1)-log(factorial(k))-k*log(2)-log(2*k+1)-log(2*k+2)
  return((-1)^k *exp(lgfk))
}

#计算求和
compute.Sum <- function(Fun,i,j,a1,d1){
  #计算函数f的第i到j项
  if(i>j){
    paste("范围不对")
  }
  else{
    s <- 0
    for (r in i:j){
      s<- s+Fun(a1,d1,r)
    }
  }
  return(s)
}

a <- c(1,2)
d <- 10

#计算前100项
compute.Sum(fk,1,100,a,d)
#计算前1000项
compute.Sum(fk,1,1000,a,d)

## -----------------------------------------------------------------------------
T <- c(0.54,0.48,0.33,0.43,0.91,0.21,0.85)
E <- c(1,1,1)
t <- 1

EM <- function(T,E,t,max.it = 10000,eps = 1e-6){
  n <- length(T)
  m <- length(E)
  lambda1 <- 3
  lambda2 <-0.5 #初始值
  while(abs(lambda1-lambda2)>=eps){
    lambda1 <- lambda2
    lambda2 <- (sum(T)+m*(lambda1+t))/(n+m)
  }
  return(lambda2)
}
paste("EM = ",EM(T,E,t,max.it = 10000,eps = 1e-6))
paste("MLE = ",(sum(T)+sum(E))/length(T))

## -----------------------------------------------------------------------------
set.seed(12315)
trims <- c(0, 0.1, 0.2, 0.5)
x <- rcauchy(100)

lapply(trims, function(trim) mean(x, trim = trim))
lapply(trims, mean, x = x)

## -----------------------------------------------------------------------------
lapply(trims, function(trim) mean(x,trim))

## -----------------------------------------------------------------------------
rm(list=ls())

formulas <- list(
mpg ~ disp,
mpg ~ I(1/disp),
mpg ~ disp + wt,
mpg ~ I(1/disp) + wt
)

boot_df <- function(x){
  return(x[sample(nrow(x), replace = TRUE), ])
} 

rsq <- function(mod) {
  summary(mod)$r.squared
}

boot_lm <- function(formu){
  dat <- boot_df(mtcars)
  rsq(lm(formu,data = dat))
}
#用lapply()
#lapply(formulas,lm,data = boot_df(mtcars))

lapply(formulas,boot_lm)


#用for

for (i in 1:length(formulas)){
  print(boot_lm(formulas[[i]]))
}



#因为用了boot_strap所以两个结果不一样

## -----------------------------------------------------------------------------
bootstraps <- lapply(1:10, function(i) {
rows <- sample(1:nrow(mtcars), rep = TRUE)
mtcars[rows, ]
})

boot_mpg_disp <- function(dat){
  rsq(lm(mpg~disp,data = boot_df(dat)))
}
#使用lapply
lapply(bootstraps, boot_mpg_disp)

#使用for循环
for (i in 1:length(bootstraps)){
  dat <- boot_df(bootstraps[[i]])
  print(rsq(lm(mpg~disp,data = dat)))
}

## -----------------------------------------------------------------------------
n <- 100
df <- data.frame(x = 1:n, y = rnorm(n,0,1))
vapply(df,sd,c(1))

## -----------------------------------------------------------------------------
n <- 100
df <- data.frame(x = rchisq(n,1), y = letters[1:10],z = rnorm(n,0,1))

vapply(df,function(dat){
  if(is.numeric(dat)==TRUE){
    return(sd(dat))
  }
  else{
    return(NA)
  }
},c(1.2))
#a <- vapply(df, class, character(1)) == "numeric" 
#vapply(df,sd,c(1))

