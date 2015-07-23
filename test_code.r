##get command line args
options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)

#Do several sims
## Get functions
source("phasekids.R")

## Params and probabilities
#### Real Parameters 
size.array=10 # size of progeny array
het.error=0.7 # het->hom error
hom.error=0.002 # hom->other error
numloci=4000
win_length=11 # size of window to phase
sims=1
errors.correct=FALSE # can assume we know error rates or not
#freqs.correct=FALSE # can assume we know freqs or not
crossovers=as.numeric(args[1]) # mean expected crossovers per chromosome; 0 = no recombination
job=args[2]

#### Set up the neutral SFS \& Probabilities
x=1:99/100 #0.01 bins of freq.
freq=1/x
sfs=as.numeric(); 
for(i in 1:99){sfs=c(sfs,rep(x[i],100*freq[i]))}
p=sample(sfs,numloci) #get freqs for all loci

# row 1 is true_gen 00, row2 is true_gen 01, row 3 is true_gen 11
# cols are obs. genotype (00,01,11)
# JY modification to add in cbind(1) to incorporate missing data coded as genotype 3
gen_error_mat<-matrix(c(1-hom.error,hom.error/2,hom.error/2,het.error/2,1-het.error,het.error/2,hom.error/2,hom.error/2,1-hom.error),byrow=T,nrow=3,ncol=3)
probs<-vector("list",3)
probs[[1]]<-cbind(gen_error_mat*matrix(c(1, 0, 0), nrow = 3,byrow=F,ncol=3),1)
probs[[2]]<-cbind(gen_error_mat*matrix(c(1/4, 1/2, 1/4), nrow = 3,byrow=F,ncol=3),1)
probs[[3]]<-cbind(gen_error_mat*matrix(c(0, 0, 1), nrow = 3,byrow=F,ncol=3),1)
gen_error_mat<-cbind(gen_error_mat,1)

## Make mom, progeny
mom.phase.errors=as.numeric()
mean.kid.geno.errors=as.numeric()
a1=ran.hap(numloci,p) #make haplotypes
a2=ran.hap(numloci,p)
true_mom=list(a1,a2) #phased 
obs_mom=add_error(a1+a2,hom.error,het.error) #convert to diploid genotype
progeny<-vector("list",size.array)
progeny<-lapply(1:size.array, function(a) kid(true_mom,true_mom,het.error,hom.error,crossovers))
	
#CHANGE error matrix
if(errors.correct==FALSE){
	het.error=rnorm(1,0.7,0.05) # random normal,roughly 0.55-0.85
  	hom.error=10^-(runif(1)*4+1) # random bad guess
}
  
#MOM GENO
estimated_mom=true_mom[[1]]+true_mom[[2]] #assumes mom was imputed with perfection or known from e.g. WGS data
#MOM PHASE
newmom=phase_mom(estimated_mom,progeny,win_length,verbose=TRUE)
estimated_hets=which(estimated_mom==1)
# can't make a phasing error at a site which is not really heterozygous, 
# nor is calling a true het site homozygous a phasing error
#  so we only check phasing error at the intersection of both
# note that this could in theory lead to super low phasing error BECAUSE of high
# genotype error, but we're gonna ignore that for now
phase_sites=intersect(which(true_mom[[1]]+true_mom[[2]]==1),estimated_hets) 
mom.phase.errors=sum(sapply(2:length(phase_sites), function(a) check_phase(estimated_hets,newmom,true_mom,phase_sites[a],phase_sites[a-1])))

#KIDS GENOS
inferred_progeny=list()
for(z in 1:length(progeny)){
	inferred_progeny[[z]]=which_phase_kid(newmom,progeny[[z]][[2]][estimated_hets] )
	mean.kid.geno.errors=mean.kid.geno.errors+(sum(abs(progeny[[z]][[1]][estimated_hets]-inferred_progeny[[z]])))/length(progeny)
}
mean.kid.geno.errors=mean.kid.geno.errors/numloci
results=c(mom.phase.errors,mean.kid.geno.errors)
write.table(file=paste("./out/out.",crossovers,".",job,".txt",sep=""),t(results),quote=F,col.names=FALSE,row.names=FALSE)
