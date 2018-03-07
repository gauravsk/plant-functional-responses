# Models to optimize for fitting competition coefficients
# Inspired by JHRL script to NJBK

### Model 1: No effect of competitors on seed production
compmodel1 <- function(par) {
  mu = par[1] 
  sigma = par[2]
  pred = mu # The predicted seed output is just the mean
  llik = dnorm(lseed, log(pred), sd = sigma, log = TRUE)
  return(sum(-1*llik))
}

### Model 2: A global lambda, that decreases due to competition with neighbors
# write a new model that has a global alpha and a global lambda
compmodel2<-function(par){
  lambda<-par[1]
  alpha<-par[2]
  sigma<-par[3]
  pred<-lambda/(1+alpha*(dens)) #this is the predictive model
  llik<-dnorm(lseed,log(pred), sd=sigma, log=TRUE) #these are the log likelihoods
  return(sum(-1*llik)) #sum of negative log likelihoods
}

### Model 3: Global alpha; Global lambda that gets modified by a site-specific coefficient
# write a new model that has a global alpha and a global lambda that is additively modified by a site-term
compmodel3<-function(par){
  lambda <- par[1]
  alpha  <- par[2]
  sigma  <- par[3]
  site1  <- par[4]
  site2  <- par[5] 
  site3  <- par[6] 
  site4  <- par[7] 
  site5  <- par[8] 
  site6  <- par[9] 
  site7  <- par[10] 
  site8  <- par[11] 
  site9  <- par[12]
  site10 <- par[13] 
  site11 <- par[14] 
  site12 <- par[15] 
  site13 <- par[16] 
  site14 <- par[17] 
  site15 <- par[18] 
  site16 <- par[19] 
  site17 <- par[20] 
  site18 <- par[21] 
  site19 <- par[22] 
  site20 <- par[23] 
  site21 <- par[24] 
  site22 <- par[25] 
  site23 <- par[26] 
  site24 <- par[27] 
  pred<-(lambda + site1*site_matrix[, 1]+ site2*site_matrix[, 2] + site3*site_matrix[, 3] + site4*site_matrix[, 4] + site5*site_matrix[, 5] + site6*site_matrix[, 6] + site7*site_matrix[, 7] + site8*site_matrix[, 8] + site9*site_matrix[, 9] + site10*site_matrix[, 10] + site11*site_matrix[, 11] + site12*site_matrix[, 12] + site13*site_matrix[, 13] + site14*site_matrix[, 14] + site15*site_matrix[, 15] + site16*site_matrix[, 16] + site17*site_matrix[, 17] + site18*site_matrix[, 18] + site19*site_matrix[, 19] + site20*site_matrix[, 20] + site21*site_matrix[, 21] + site22*site_matrix[, 22] + site23*site_matrix[, 23] + site24*site_matrix[, 24])/(1+alpha*(dens)) #this is the predictive model
  llik<-dnorm(lseed,log(pred), sd=sigma, log=TRUE) #these are the log likelihoods
  return(sum(-1*llik)) #sum of negative log likelihoods
}

### Model 4: Global alpha; 24 site-specific lambda terms per species
# write a new model that has a global alpha and 24 site-specific lambda terms per species
compmodel4<-function(par){
  alpha  <- par[1]
  sigma  <- par[2]
  site1  <- par[3]
  site2  <- par[4] 
  site3  <- par[5] 
  site4  <- par[6] 
  site5  <- par[7] 
  site6  <- par[8] 
  site7  <- par[9] 
  site8  <- par[10] 
  site9  <- par[11]
  site10 <- par[12] 
  site11 <- par[13] 
  site12 <- par[14] 
  site13 <- par[15] 
  site14 <- par[16] 
  site15 <- par[17] 
  site16 <- par[18] 
  site17 <- par[19] 
  site18 <- par[20] 
  site19 <- par[21] 
  site20 <- par[22] 
  site21 <- par[23] 
  site22 <- par[24] 
  site23 <- par[25] 
  site24 <- par[26] 
  pred<-(site1*site_matrix[, 1]+ site2*site_matrix[, 2] + site3*site_matrix[, 3] + site4*site_matrix[, 4] + site5*site_matrix[, 5] + site6*site_matrix[, 6] + site7*site_matrix[, 7] + site8*site_matrix[, 8] + site9*site_matrix[, 9] + site10*site_matrix[, 10] + site11*site_matrix[, 11] + site12*site_matrix[, 12] + site13*site_matrix[, 13] + site14*site_matrix[, 14] + site15*site_matrix[, 15] + site16*site_matrix[, 16] + site17*site_matrix[, 17] + site18*site_matrix[, 18] + site19*site_matrix[, 19] + site20*site_matrix[, 20] + site21*site_matrix[, 21] + site22*site_matrix[, 22] + site23*site_matrix[, 23] + site24*site_matrix[, 24])/(1+alpha*(dens)) #this is the predictive model
  llik<-dnorm(lseed,log(pred), sd=sigma, log=TRUE) #these are the log likelihoods
  return(sum(-1*llik)) #sum of negative log likelihoods
}

### Model 5: 24 site-specific lambda terms and 24 site-specific alpha terms per species
compmodel5<-function(par){
  # The sigma is still global
  sigma  <- par[1]
  
  # Site specific alphas 
  alpha1  <- par[2]
  alpha2  <- par[3] 
  alpha3  <- par[4] 
  alpha4  <- par[5] 
  alpha5  <- par[6] 
  alpha6  <- par[7] 
  alpha7  <- par[8] 
  alpha8  <- par[9] 
  alpha9  <- par[10]
  alpha10 <- par[11] 
  alpha11 <- par[12] 
  alpha12 <- par[13] 
  alpha13 <- par[14] 
  alpha14 <- par[15] 
  alpha15 <- par[16] 
  alpha16 <- par[17] 
  alpha17 <- par[18] 
  alpha18 <- par[19] 
  alpha19 <- par[20] 
  alpha20 <- par[21] 
  alpha21 <- par[22] 
  alpha22 <- par[23] 
  alpha23 <- par[24] 
  alpha24 <- par[25] 
  
  # Site specific lambas 
  site1 <- par[26]
  site2 <- par[27]
  site3 <- par[28]
  site4 <- par[29]
  site5 <- par[30]
  site6 <- par[31]
  site7 <- par[32]
  site8 <- par[33]
  site9 <- par[34]
  site10 <- par[35]
  site11 <- par[36]
  site12 <- par[37]
  site13 <- par[38]
  site14 <- par[39]
  site15 <- par[40]
  site16 <- par[41]
  site17 <- par[42]
  site18 <- par[43]
  site19 <- par[44]
  site20 <- par[45]
  site21 <- par[46]
  site22 <- par[47]
  site23 <- par[48]
  site24 <- par[49]
  
  
  # Predictive model
  
  pred<-(site1*site_matrix[, 1]+ site2*site_matrix[, 2] + site3*site_matrix[, 3] + site4*site_matrix[, 4] + 
           site5*site_matrix[, 5] + site6*site_matrix[, 6] + site7*site_matrix[, 7] + site8*site_matrix[, 8] + 
           site9*site_matrix[, 9] + site10*site_matrix[, 10] + site11*site_matrix[, 11] + site12*site_matrix[, 12] + 
           site13*site_matrix[, 13] + site14*site_matrix[, 14] + site15*site_matrix[, 15] + site16*site_matrix[, 16] + 
           site17*site_matrix[, 17] + site18*site_matrix[, 18] + site19*site_matrix[, 19] + site20*site_matrix[, 20] + 
           site21*site_matrix[, 21] + site22*site_matrix[, 22] + site23*site_matrix[, 23] + site24*site_matrix[, 24]) /
    (1+ alpha1*compdens[, 1]+ alpha2*compdens[, 2] + alpha3*compdens[, 3] + alpha4*compdens[, 4] + 
       alpha5*compdens[, 5] + alpha6*compdens[, 6] + alpha7*compdens[, 7] + alpha8*compdens[, 8] + 
       alpha9*compdens[, 9] + alpha10*compdens[, 10] + alpha11*compdens[, 11] + alpha12*compdens[, 12] + 
       alpha13*compdens[, 13] + alpha14*compdens[, 14] + alpha15*compdens[, 15] + alpha16*compdens[, 16] + 
       alpha17*compdens[, 17] + alpha18*compdens[, 18] + alpha19*compdens[, 19] + alpha20*compdens[, 20] + 
       alpha21*compdens[, 21] + alpha22*compdens[, 22] + alpha23*compdens[, 23] + alpha24*compdens[, 24]) #this is the predictive model
  llik<-dnorm(lseed,log(pred), sd=sigma, log=TRUE) #these are the log likelihoods
  return(sum(-1*llik)) #sum of negative log likelihoods
}