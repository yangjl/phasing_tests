#Phasing Code

Uses selfed progeny from a single mom to 

1) impute mom's genotype
2) phase mom
3) impute self progeny genotypes

Currently assumes data is coded as 0: AA 1: Aa 2: aa and 3 for missing data.


------------------
# split JRI's codes and test independently

```
### loading all the functions in folder "lib"
sourceall <-function(rm=FALSE){
    if(rm) rm(list=ls())
    f <- sapply(list.files(pattern="[.]R$", path="lib", full.names=TRUE), source)
}
sourceall(rm=T)
```


------------
# File Formats (PLINK format)

It accepts sequences of the form `0/1/2/3` where `0/2` denote the 
genotypes that are homozygous, `1` denotes a heterozygous genotype, and `3` denotes an unknown genotype.

### The PED file format is as follows (`plink --no-sex, --no-pheno`):    
- tab delimited
- 'family ID' 'plant ID' 'Paternal ID' 'Maternal ID' 'genotype sequence'

### The MAP file format
- 'chr' 'SNPID' 'Base-pair position'

### For example:
```

```



