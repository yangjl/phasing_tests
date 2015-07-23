

#MOM PHASE
newmom=phase_mom(estimated_mom,progeny,win_length,verbose=TRUE)
estimated_hets=which(estimated_mom==1)
# can't make a phasing error at a site which is not really heterozygous, 
# nor is calling a true het site homozygous a phasing error
#  so we only check phasing error at the intersection of both
# note that this could in theory lead to super low phasing error BECAUSE of high
# genotype error, but we're gonna ignore that for now
phase_sites=intersect(which(true_mom[[1]]+true_mom[[2]]==1),estimated_hets) 
mom.phase.errors[mysim]=sum(sapply(2:length(phase_sites), 
                                   function(a) check_phase(estimated_hets,newmom,true_mom,phase_sites[a],phase_sites[a-1])))

