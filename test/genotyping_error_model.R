library(sydneyPaternity)

nall=4
E1 = 0.23
E2 = 0.46
e1 = E1
e2 = E2

#test behavior for monomorphic
sydneyPaternity:::genotyping_error_model(c(0,0),0,0,1,E1,E2)

#test calculation of probabilities
geno <- cbind(combn(1:nall,2)); for(i in 1:nall) geno <- cbind(geno, i)
geno <- geno[,order(geno[1,],geno[2,])]
out <- c()
for(i in 1:ncol(geno))
  for(j in 1:ncol(geno))
{
  val <- sydneyPaternity:::genotyping_error_model(geno[,i],geno[1,j],geno[2,j],4,E1,E2)
  cls <- sydneyPaternity:::genotyping_error_model_class(geno[,i],geno[1,j],geno[2,j])
  out <- rbind(out, data.frame(prob=val, class=cls, phenotype=paste(geno[,i],collapse="/"), genotype=paste(geno[,j],collapse="/")))
}
sum(out[out$genotype=="1/2",1])
sum(out[out$genotype=="1/1",1])

sim_error_model <- function(genotype, nall, err1, err2, err3)
{
  nall <- length(nall)
  phenotype <- genotype
  if (genotype[1] != genotype[2])
  {
    is_err1 <- runif(1) < 2*err1
    if (is_err1)
    {
      err1_what_allele <- sample(1:2, 1)
      phenotype[] <- phenotype[err1_what_allele]
    }
  }
  for (i in 1:2)
  {
    is_err2 <- runif(1) < err2
    if (is_err2)
    {
      phenotype[i] <- sample( (1:nall)[-phenotype[i]], 1)
    }
  }
  return(sort(phenotype))
}
table(apply(replicate(100000,sim_error_model(c(1,2),1:4,E1,E2,0)),2,paste,collapse="/"))/100000
table(apply(replicate(100000,sim_error_model2(c(1,2),1:4,E1,E2,0)),2,paste,collapse="/"))/100000
out[out$genotype=="1/2",]
table(apply(replicate(100000,sim_error_model(c(1,1),1:4,E1,E2,0)),2,paste,collapse="/"))/100000
table(apply(replicate(100000,sim_error_model2(c(1,1),1:4,E1,E2,0)),2,paste,collapse="/"))/100000
out[out$genotype=="1/1",]

#test event simulation probabilities
sim_error_model_events <- function(genotype, nall, err1, err2)
{
  err_count <- c(0,0)
  phenotype <- genotype
  if (genotype[1] != genotype[2])
  {
    is_err1 <- runif(1) < 2*err1
    if (is_err1)
    {
      err1_what_allele <- sample(1:2, 1)
      phenotype[] <- phenotype[err1_what_allele]
    }
    err_count[1] <- err_count[1] + is_err1
  }
  for (i in 1:2)
  {
    is_err2 <- runif(1) < err2
    if (is_err2)
    {
      phenotype[i] <- sample( (1:nall)[-phenotype[i]], 1)
    }
    err_count[2] <- err_count[2] + is_err2
  }
  cls <- sydneyPaternity:::genotyping_error_model_class(phenotype,genotype[1],genotype[2])
  return(c(cls,paste(err_count,collapse=":")))
}
oof11 <- replicate(100000,sim_error_model_events(c(1,1),4,E1,E2))
prop.table(table(oof11[1,],oof11[2,])/100000,1)
table(apply(replicate(10000, sydneyPaternity:::simulate_genotyping_errors(c(1,1),1,1,4,E1,E2)),3,paste,collapse=":"))/10000
table(apply(replicate(10000, sydneyPaternity:::simulate_genotyping_errors(c(1,2),1,1,4,E1,E2)),3,paste,collapse=":"))/10000
table(apply(replicate(10000, sydneyPaternity:::simulate_genotyping_errors(c(2,2),1,1,4,E1,E2)),3,paste,collapse=":"))/10000
oof12 <- replicate(100000,sim_error_model_events(c(1,2),4,E1,E2))
prop.table(table(oof12[1,],oof12[2,])/100000,1)
table(apply(replicate(10000, sydneyPaternity:::simulate_genotyping_errors(c(1,2),1,2,4,E1,E2)),3,paste,collapse=":"))/10000
table(apply(replicate(10000, sydneyPaternity:::simulate_genotyping_errors(c(2,2),1,2,4,E1,E2)),3,paste,collapse=":"))/10000
table(apply(replicate(10000, sydneyPaternity:::simulate_genotyping_errors(c(1,4),2,3,4,E1,E2)),3,paste,collapse=":"))/10000
table(apply(replicate(10000, sydneyPaternity:::simulate_genotyping_errors(c(1,2),2,3,4,E1,E2)),3,paste,collapse=":"))/10000
