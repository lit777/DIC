library(lattice)
library("R2OpenBUGS")

M <- mean(y)
N <- length(y)

#########################
#### BUGS code ##############
#########################

data<-list("y","N","M")
inits<-function(){list(mu<-c(M,M),  tau<-c(1,1) )}

model<-function(){
for(i in 1:N){
        S[i]~dcat(pi[])
    y[i]~dnorm(mu[S[i]], tau[S[i]])
		}

for(j in 1:2){
	tau[j]~dgamma(1,1)
    sigma[j]<-1/tau[j]
    pii[j]<-1
	}

    mu[1]~dnorm(M,0.01)
        mu[2]~dnorm(M,0.01)
 
    pi[1:2]~ddirch(pii[])
}
write.model(model, "model.txt")
model=paste(getwd(),"model.txt",sep="/")

fit2<-bugs(data,inits, model.file=model,parameters=c("pi","mu","sigma","S"),n.chain=1, n.iter=25000, n.burnin=15000,debug=TRUE, DIC=TRUE, n.thin=5)

#########################
#### DIC 1 ##############
#########################

log.score <- function(pi1,pi2,mu1,mu2,sigma1,sigma2){
    rowSums(sapply(y, function(x) log(pmax(0.1^100,pi1*dnorm(x,mu1,sqrt(sigma1))+pi2*dnorm(x,mu2, sqrt(sigma2))))))
}

log.score0 <- function(pi1,pi2,mu1,mu2,sigma1,sigma2){
    sum(sapply(y, function(x) log(max(0.1^100,pi1*dnorm(x,mu1,sqrt(sigma1))+pi2*dnorm(x,mu2, sqrt(sigma2))))))
}

d1 <- fit2$sims.array[1:10000,1,1]
d2 <- fit2$sims.array[1:10000,1,2]
d3 <- fit2$sims.array[1:10000,1,3]
d4 <- fit2$sims.array[1:10000,1,4]
d5 <- fit2$sims.array[1:10000,1,5]
d6 <- fit2$sims.array[1:10000,1,6]

bar.d <- -4*mean(log.score(d1,d2,d3,d4,d5,d6))

dbar <- 2*log.score0(fit2$summary[1,1],fit2$summary[2,1],fit2$summary[3,1],fit2$summary[4,1],fit2$summary[5,1],fit2$summary[6,1])

(DIC1 <- bar.d+dbar); (pd1 <- bar.d/2+dbar)


#########################
#### DIC 2 ##############
#########################

log.score <- function(pi1,pi2,mu1,mu2,sigma1,sigma2){
  rowSums(sapply(y, function(x) log(pmax(0.1^100,pi1*dnorm(x,mu1,sqrt(sigma1))+pi2*dnorm(x,mu2,sqrt(sigma2))))))
}

log.score0<-function(pi1,pi2,mu1,mu2,sigma1,sigma2){
  sum(sapply(y, function(x) log(max(0.1^100,pi1*dnorm(x,mu1,sqrt(sigma1))+pi2*dnorm(x,mu2, sqrt(sigma2))))))
}

d1 <- fit2$sims.array[1:10000,1,1]
d2 <- fit2$sims.array[1:10000,1,2]
d3 <- fit2$sims.array[1:10000,1,3]
d4 <- fit2$sims.array[1:10000,1,4]
d5 <- fit2$sims.array[1:10000,1,5]
d6 <- fit2$sims.array[1:10000,1,6]

bar.d <- -4*mean(log.score(d1,d2,d3,d4,d5,d6))

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

dbar <- 2*log.score0(getmode(d1),getmode(d2),getmode(d3),getmode(d4),getmode(d5),getmode(d6))

(DIC2 <- bar.d+dbar); (pd2 <- bar.d/2+dbar)



#########################
#### DIC 3 ##############
#########################

log.score <- function(pi1,pi2,mu1,mu2,sigma1,sigma2){
    rowSums(sapply(y, function(x) log(pmax(0.1^100,pi1*dnorm(x,mu1,sqrt(sigma1))+pi2*dnorm(x,mu2,sqrt(sigma2))))))
}

d1 <- fit2$sims.array[1:10000,1,1]
d2 <- fit2$sims.array[1:10000,1,2]
d3 <- fit2$sims.array[1:10000,1,3]
d4 <- fit2$sims.array[1:10000,1,4]
d5 <- fit2$sims.array[1:10000,1,5]
d6 <- fit2$sims.array[1:10000,1,6]

bar.d <- -4*mean(log.score(d1,d2,d3,d4,d5,d6))

score <- function(y,pi1,pi2,mu1,mu2,sigma1,sigma2){
    pi1*dnorm(y,mu1,sqrt(sigma1))+pi2*dnorm(y,mu2,sqrt(sigma2))
}

phat <- 2*log(prod(rowMeans(apply(fit2$sims.array[1:10000,1,1:6],1, function(x) score(y,x[1],x[2],x[3],x[4],x[5],x[6])))))

(DIC.star<-bar.d+phat);(pstar<-bar.d/2+phat)



#########################
#### DIC 4 ##############
#########################

log.score <- function(pi1,pi2,mu1,mu2,sigma1,sigma2){
    rowSums(sapply(y, function(x) pi1*dnorm(x,mu1,sqrt(sigma1))/(pi1*dnorm(x,mu1,sqrt(sigma1))+pi2*dnorm(x,mu2, sqrt(sigma2)))*log(pi1*dnorm(x,mu1,sqrt(sigma1)))+pi2*dnorm(x,mu2,sqrt(sigma2))/(pi1*dnorm(x,mu1,sqrt(sigma1))+pi2*dnorm(x,mu2, sqrt(sigma2)))*log(pi2*dnorm(x,mu2, sqrt(sigma2)))))
}

d1 <- fit2$sims.array[1:10000,1,1]
d2 <- fit2$sims.array[1:10000,1,2]
d3 <- fit2$sims.array[1:10000,1,3]
d4 <- fit2$sims.array[1:10000,1,4]
d5 <- fit2$sims.array[1:10000,1,5]
d6 <- fit2$sims.array[1:10000,1,6]

bar.d <- -4*mean(log.score(d1,d2,d3,d4,d5,d6))

z1 <- ifelse(fit2$sims.array[1:10000,1,7:170]==1, 1, 0)
z2 <- ifelse(fit2$sims.array[1:10000,1,7:170]==2, 1, 0)

m1 <- rowSums(z1)
m2 <- rowSums(z2)

K<-2

### need to derive the following quantities based on (hyper-)prior settings ###

hat.pi1 <- (1+m1)/(1*K+164)
hat.pi2 <- (1+m2)/(1*K+164)

hat.mu1 <- (M*0.01+(1/d5)*(z1%*%y))/(0.01+m1/d5)
hat.mu2 <- (M*0.01+(1/d6)*(z2%*%y))/(0.01+m2/d6)

MU1 <- (z1%*%y)/m1
MU2 <- (z2%*%y)/m2

s1 <- colSums(sapply(MU1, function(x) (x-y)^2)*t(z1))
s2 <- colSums(sapply(MU2, function(x) (x-y)^2)*t(z2))

hat.tau1 <- (1+m1/2)/(1+s1/2)
hat.tau2 <- (1+m2/2)/(1+s2/2)

###############################################################################

log.score0 <- function(hat.pi1,hat.pi2,hat.mu1,hat.mu2,hat.tau1,hat.tau2){
    rowSums(sapply(1:164, function(x) log(hat.pi1*dnorm(y[x],hat.mu1,sqrt(1/hat.tau1)))*z1[,x]+log(hat.pi2*dnorm(y[x],hat.mu2, sqrt(1/hat.tau2)))*z2[,x]))
}

dbar<- 2*mean(log.score0(hat.pi1,hat.pi2,hat.mu1,hat.mu2,hat.tau1,hat.tau2))

(DIC4<-bar.d+dbar); (pd2<-bar.d/2+dbar)



###########################
#### DIC_loo ##############
###########################

score <- function(y,pi1,pi2,mu1,mu2,sigma1,sigma2){
    pi1*dnorm(y,mu1,sqrt(sigma1))+pi2*dnorm(y,mu2,sqrt(sigma2))
}
wt <- apply(fit2$sims.array[1:10000,1,1:6],1, function(x) score(y,x[1],x[2],x[3],x[4],x[5],x[6]) )

d1 <- fit2$sims.array[1:10000,1,1]
d2 <- fit2$sims.array[1:10000,1,2]
d3 <- fit2$sims.array[1:10000,1,3]
d4 <- fit2$sims.array[1:10000,1,4]
d5 <- fit2$sims.array[1:10000,1,5]
d6 <- fit2$sims.array[1:10000,1,6]

rs <- 1/rowMeans(1/wt)

bartemp <- sum(sapply(1:164, function(t) -4*mean(rs[t]*log(score(y[t],d1,d2,d3,d4,d5,d6))/wt[t,] )) )

temp.bar <- sum(sapply(1:164, function(t) 2*log(mean(rs[t]*score(y[t],d1,d2,d3,d4,d5,d6)/wt[t,] ))))

(DIC.loo<-bartemp+temp.bar); (ploo<-bartemp/2+temp.bar)


############################
#### DIC_cloo ##############
############################

score <- function(y,pi1,pi2,mu1,mu2,sigma1,sigma2){
    pi1*dnorm(y,mu1,sqrt(sigma1))+pi2*dnorm(y,mu2,sqrt(sigma2))
}

wt <- apply(fit2$sims.array[1:10000,1,1:6],1, function(x) score(y,x[1],x[2],x[3],x[4],x[5],x[6]))

d1 <- fit2$sims.array[1:10000,1,1]
d2 <- fit2$sims.array[1:10000,1,2]
d3 <- fit2$sims.array[1:10000,1,3]
d4 <- fit2$sims.array[1:10000,1,4]
d5 <- fit2$sims.array[1:10000,1,5]
d6 <- fit2$sims.array[1:10000,1,6]

rs <- 1/rowMeans(1/wt)

phat <- 2*log(  prod(rowMeans(apply(fit2$sims.array[1:10000,1,1:6],1, function(x) score(y,x[1],x[2],x[3],x[4],x[5],x[6])))))

temp.bar <- sum(sapply(1:164, function(t) 2*log(mean(rs[t]*score(y[t],d1,d2,d3,d4,d5,d6)/wt[t,] ))))

(DIC.cloo<- -phat + phat - temp.bar); (pcloo<- phat - temp.bar)

