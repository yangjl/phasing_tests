#################################################
# Example :
#################################################
seqLen = 100
stayPrb = rep(.95, 2)
tM = MakeTransMat(stayPrb)
seqMeans = c(0, 1)
seqSDs = c(1, 1)
x = GenSeq(transMat = tM, seqLen = seqLen,
           seqMeans = seqMeans, seqSDs = seqSDs)
fwdPrbs = ComputeFwd(x$O, transMat = tM, seqMeans = seqMeans, seqSDs = seqSDs)
hidSamp = SampleSeq(nSamp = 2, transMat = tM, fwdProbs = fwdPrbs$fwdProbs)
plot(x$O, col="gray", xlab="index", ylab="")
cv = rainbow( ncol(hidSamp) )
lines(x$Q - 1, lwd=5)
for(k in 1:ncol(hidSamp) ){
    lines(hidSamp[,k] + runif(1)*.03 - 1,
          col=cv[k], lwd=3, lty="dashed")
}

#######
seqLen = 1000
stayPrb = rep(.98, 4)
tM = MakeTransMat(stayPrb)
seqMeans = 0:( length(stayPrb) - 1)
seqSDs = runif( length(stayPrb) )
x = GenSeq(transMat = tM, seqLen = seqLen,
           seqMeans = seqMeans, seqSDs = seqSDs)
fwdPrbs = ComputeFwd(x$O, transMat = tM, seqMeans = seqMeans, seqSDs = seqSDs)
hidSamp = SampleSeq(nSamp = 2, transMat = tM, fwdProbs = fwdPrbs$fwdProbs)
plot(x$O, col="gray", xlab="index", ylab="")
cv = rainbow( ncol(hidSamp) )
lines(x$Q - 1, lwd=5)
for(k in 1:ncol(hidSamp) ){
    lines(hidSamp[,k] + runif(1)*.03 - 1,
          col=cv[k], lwd=3, lty="dashed")
}


###################
### Gaussian mixture model

### simulate data from Gaussian mixture model

### model parameters
sigma=2; mu1=-2; mu2=4;
N=5000;
coin=rbinom(N,1, .3);

r1=rnorm(N,mu1,sigma);
r2=rnorm(N,mu2,sigma);

rr=r1;
rr[which(coin==0)]=r2[which(coin==0)];

tt=seq(-10,15,by=.1);
d1=dnorm(tt,mu1,sigma);
d2=dnorm(tt, mu2, sigma);
### plot data;
hist(rr, prob=T, col="lightgray", border="gray", main="");
lines(tt, d1*.3, col=3, lwd=1.5);
lines(tt, d2*.7, col=4, lwd=1.5);
lines(tt, d1*.3+d2*.7, col=2, lwd=2, lty=2);


### Gaussian mixture model

### simulate data from Gaussian mixture model

### model parameters
sigma=2; mu1=-2; mu2=4;
N=5000;
coin=rbinom(N,1, .3);

r1=rnorm(N,mu1,sigma);
r2=rnorm(N,mu2,sigma);

rr=r1;
rr[which(coin==0)]=r2[which(coin==0)];

tt=seq(-10,15,by=.1);
d1=dnorm(tt,mu1,sigma);
d2=dnorm(tt, mu2, sigma);
### plot data;
hist(rr, prob=T, col="lightgray", border="gray", main="");
lines(tt, d1*.3, col=3, lwd=1.5);
lines(tt, d2*.7, col=4, lwd=1.5);
lines(tt, d1*.3+d2*.7, col=2, lwd=2, lty=2);


### now run EM to estimate parameters;
### show parameter estimates at each iteration

m1=-6; m2=7; s=4; pp=.5; ### initial guess
w=runif(N,.2, .8); ## conditional expectation of z

iter=20; ## number of EM iterations; can use convergence criterion

m.mat=matrix(0,nrow=iter+1,ncol=4);
m.mat[1,]=c(m1,m2, s, pp);

for (i in 1:iter) {
    hist(rr, prob=T, col="lightgray", border="gray", main="", xlab="y", nclass=50);
    d1=dnorm(tt,m1,s);
    d2=dnorm(tt, m2, s);
    lines(tt, d1*pp, col=3, lwd=1.5);
    lines(tt, d2*(1-pp), col=4, lwd=1.5);
    lines(tt, d1*pp+d2*(1-pp), col=2, lwd=2, lty=2);
    abline(v=c(m1,m2), col=c(3,4), lty=2)
    points(c(mu1,mu2), rep(.001,2), pch=17, col=c(3,4))
    
    #readline(prompt="next step")
    ## estimate new attribution
    h1=dnorm(rr,m1,s);
    h2=dnorm(rr,m2,s);
    w=pp*h1/(pp*h1+(1-pp)*h2);
    pp=mean(w);
    ## Max mu, s
    m1=sum(w*rr)/sum(w);
    m2=sum((1-w)*rr)/sum((1-w));
    s1=sum(w*(rr-m1)^2);
    s2=sum((1-w)*(rr-m2)^2);
    s=sqrt((s1+s2)/N);
    
    m.mat[i+1,]=c(m1,m2,s,pp);
    Sys.sleep(0.5)
    
    
}