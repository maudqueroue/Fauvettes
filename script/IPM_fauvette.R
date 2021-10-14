##########################
#     IPM PASSEREAUX     #
#        Fauvettes       #
##########################

rm(list=ls())
#Package
library(nimble)
library(nimbleEcology)
library(tidyverse)

#################
# Fonction 
#################

dCJS_vv_sum <- nimbleFunction(
  # It is assumed that the individual has already been captured.
  # Therefore, the first entry in x represents the first possible recapture event.
  # probSurvive[t] represents survival from t-1 to t.
  # probCapture[t] represents capture probability at time t.
  run = function(x = double(1),    ## standard name for the "data"
                 probSurvive = double(1),
                 probCapture = double(1),
                 mult = double(0), #! NEWLY ADDED: argument stating number of occurences of same capture history in entire dataset 
                 len = double(0, default = 0),
                 log = integer(0, default = 0) ## required log argument
  ) {
    if (len != 0) {
      if (len != length(x)) stop("Argument len must match length of data, or be 0.")
    }
    if (length(probSurvive) < length(x) - 1)
      stop("Length of probSurvive must be at least length of data minus 1.")
    if (length(x) != length(probCapture)) stop("Length of probCapture does not match length of data.")
    if (x[1] != 1) stop("dCJS requires specifying first capture. x[1] must equal 1.")
    
    ## Note the calculations used here are actually in hidden Markov model form.
    probAliveGivenHistory <- 1
    ## logProbData will be the final answer
    logProbData <- 0
    if (len == 0) {  ## l<1 should not occur, but just in case:
      len <- length(x)
    }
    for (t in 2:len) {
      ## probAlive is P(Alive(t) | x(1)...x(t-1))
      ## probAliveGivenHistory is (Alive(t-1) | x(1)...x(t-1))
      probAlive <- probAliveGivenHistory * probSurvive[t - 1]
      if (!is.na(x[t])) {
        if (x[t] == 1) {
          ## ProbThisObs = P(x(t) | x(1)...x(t-1))
          probThisObs <- probAlive * probCapture[t]
          probAliveGivenHistory <- 1
        } else {
          probAliveNotSeen <- probAlive * (1 - probCapture[t])
          probThisObs <- probAliveNotSeen + (1 - probAlive)
          probAliveGivenHistory <- probAliveNotSeen / probThisObs
        }
        logProbData <- logProbData + log(probThisObs) * mult #! NEWLY ADDED: "mult"
      }
    }
    if (log) {
      return(logProbData)
    }
    return(exp(logProbData))
    returnType(double())
  }
)

rCJS_vv_sum <- nimbleFunction(
  run = function(n = integer(),
                 probSurvive = double(1),
                 probCapture = double(1),
                 mult = double(0), #! NEWLY ADDED: argument stating number of occurences of same capture history in entire dataset 
                 len = double(0, default = 0)) {
    if (n != 1) stop("rCJS only works for n = 1")
    if (len < 2)
      stop("len must be greater than 1.")
    if(length(probSurvive) != len - 1)
      stop("Length of probSurvive is not the same as len - 1.")
    if(length(probCapture) != len)
      stop("Length of probCapture is not the same as len.")
    ans <- numeric(length = len, init = FALSE)
    ans[1] <- 1
    alive <- 1
    if (len <= 0) return(ans)
    for (i in 2:len) {
      if (alive)
        alive <- rbinom(1, size = 1, prob = probSurvive[i - 1])
      if (alive) {
        ans[i] <- rbinom(1, size = 1, prob = probCapture[i])
      } else {
        ans[i] <- 0
      }
    }
    return(ans)
    returnType(double(1))
  }
)
# ----------------------------------------------------------------

#################
#     Data      #
#################

#Histoire de vie individus
load('hvie_sylbor.RData')
load('hvie_sylatr.RData')

# Histoire de vie des sites
load('hvie_ID_PROG_syl.RData')


# Load decalage pheno
load("DP.RData")
dp <- scale(DP)[,1]

load("DP_sylatr.RData")
dp_bc <- scale(c(DP_sylatr,0))[,1]

# Nb ind
N_gw <- dim(hvie_sylbor)[1]
N_bc <- dim(hvie_sylatr)[1]

# Nb years
K <- 19

#############################
#  Regroupement d'individus #
#############################

# on renomme les colonnes
hvie_sylbor <- hvie_sylbor %>%
  dplyr::rename(
    a_2001 = '2001',
    a_2002 = '2002',
    a_2003 = '2003',
    a_2004 = '2004',
    a_2005 = '2005',
    a_2006 = '2006',
    a_2007 = '2007',
    a_2008 = '2008',
    a_2009 = '2009',
    a_2010 = '2010',
    a_2011 = '2011',
    a_2012 = '2012',
    a_2013 = '2013',
    a_2014 = '2014',
    a_2015 = '2015',
    a_2016 = '2016',
    a_2017 = '2017',
    a_2018 = '2018',
    a_2019 = '2019')

# Nouvelle colonne : nb_ind_groupe : nb d individu partageant la meme histoire de capture
hvie_sylbor <- hvie_sylbor %>% 
  dplyr::group_by(a_2001,a_2002,a_2003,a_2004,a_2005,a_2006,a_2007,a_2008,a_2009,
                  a_2010,a_2011,a_2012,a_2013,a_2014,a_2015,a_2016,a_2017,a_2018,
                  a_2019,censure,ID_PROG,cov_ind,vector_ID_PROG) %>%
  dplyr::summarise(nb_ind_groupe=dplyr::n())


# on renomme les colonnes
hvie_sylatr <- hvie_sylatr %>%
  dplyr::rename(
    a_2001 = '2001',
    a_2002 = '2002',
    a_2003 = '2003',
    a_2004 = '2004',
    a_2005 = '2005',
    a_2006 = '2006',
    a_2007 = '2007',
    a_2008 = '2008',
    a_2009 = '2009',
    a_2010 = '2010',
    a_2011 = '2011',
    a_2012 = '2012',
    a_2013 = '2013',
    a_2014 = '2014',
    a_2015 = '2015',
    a_2016 = '2016',
    a_2017 = '2017',
    a_2018 = '2018',
    a_2019 = '2019')

# Nouvelle colonne : nb_ind_groupe : nb d individu partageant la meme histoire de capture
hvie_sylatr <- hvie_sylatr %>% 
  dplyr::group_by(a_2001,a_2002,a_2003,a_2004,a_2005,a_2006,a_2007,a_2008,a_2009,
                  a_2010,a_2011,a_2012,a_2013,a_2014,a_2015,a_2016,a_2017,a_2018,
                  a_2019,censure,ID_PROG,cov_ind,vector_ID_PROG) %>%
  dplyr::summarise(nb_ind_groupe=dplyr::n())

########################################################
# On enleve les hvie pour les individus avec f[i] = l[i]
########################################################

# Histoire de vie des mesanges en 0 et 1 (pas vu/vu)
#-----------------------------
# On remplace par un tous les 1 et 2 en 1 (=vu)
# il ne reste que des 1 (vu) ou des 0 (pas vu)
mydata_gw <- as.matrix(hvie_sylbor[,1:K])
mydata_gw <- apply(mydata_gw,2,as.numeric)
mydata_gw[which(mydata_gw>0)] <- 1

mydata_bc <- as.matrix(hvie_sylatr[,1:K])
mydata_bc <- apply(mydata_bc,2,as.numeric)
mydata_bc[which(mydata_bc>0)] <- 1


# Nb individu
#-------------------------------
N_gw <- nrow(mydata_gw)
N_bc <- nrow(mydata_bc)

# Compute the date of first capture for each individual:
##-------------------------------
f_gw <- NULL
for (i in 1:N_gw){
  temp <- 1:K
  f_gw <- c(f_gw,min(temp[mydata_gw[i,]==1]))}

f_bc <- NULL
for (i in 1:N_bc){
  temp <- 1:K
  f_bc <- c(f_bc,min(temp[mydata_bc[i,]==1]))}


# Compute the last date of capture for each individual with censure:
#-------------------------------
l_gw <- NULL
for (i in 1:N_gw){
  temp <- 1:K
  if (hvie_sylbor$censure[i]==1){l_new <- K}
  else {l_new <- max(temp[mydata_gw[i,]==1])}
  l_gw <- c(l_gw,l_new)
  rm(l_new)}

l_bc <- NULL
for (i in 1:N_bc){
  temp <- 1:K
  if (hvie_sylatr$censure[i]==1){l_new <- K}
  else {l_new <- max(temp[mydata_bc[i,]==1])}
  l_bc <- c(l_bc,l_new)
  rm(l_new)}


# Pour quels individus l[i]=f[i]
#----------------------------
id_gw <- NULL
for(i in 1:N_gw){
  if(l_gw[i]==f_gw[i]) {
    id_gw[i] <- "id"
  }
  else{id_gw[i] <- "ok"}
}

ID_suppr_gw <- which(id_gw=="id")

id_bc <- NULL
for(i in 1:N_bc){
  if(l_bc[i]==f_bc[i]) {
    id_bc[i] <- "id"
  }
  else{id_bc[i] <- "ok"}
}

ID_suppr_bc <- which(id_bc=="id")


# On les supprime
#---------------------
hvie_sylbor <- hvie_sylbor[-ID_suppr_gw,]
hvie_sylatr <- hvie_sylatr[-ID_suppr_bc,]


# Choses inutiles
#----------------------
rm(ID_suppr_gw,l_gw,f_gw,mydata_gw,N_gw,id_gw,
   ID_suppr_bc,l_bc,f_bc,mydata_bc,N_bc,id_bc,
   temp,i,DP, DP_sylatr)

###############################
# Formatage des donnees
###############################


# Histoire de vie des mesanges en 0 et 1 (pas vu/vu)
#-----------------------------
# On remplace par un tous les 1 et 2 en 1 (=vu)
# il ne reste que des 1 (vu) ou des 0 (pas vu)
mydata_gw <- as.matrix(hvie_sylbor[,1:K])
mydata_gw <- apply(mydata_gw,2,as.numeric)
mydata_gw[which(mydata_gw>0)] <- 1

mydata_bc <- as.matrix(hvie_sylatr[,1:K])
mydata_bc <- apply(mydata_bc,2,as.numeric)
mydata_bc[which(mydata_bc>0)] <- 1

# Covariable age 
#----------------------
# On cree la matrice des covariables pour l'age 
cov_age_gw <- as.matrix(hvie_sylbor[,1:K])
cov_age_gw <- apply(cov_age_gw,2,as.numeric)

cov_age_bc <- as.matrix(hvie_sylatr[,1:K])
cov_age_bc <- apply(cov_age_bc,2,as.numeric)

# Dans quel site a ete capturee chaque mesange
#---------------------------
# Vecteur donnant le site dans lequel a ete capture chaque mesange 
ind_site_gw <- hvie_sylbor$vector_ID_PROG
ind_site_bc <- hvie_sylatr$vector_ID_PROG

# vecteur covariable individuelle recapture
# --------------------------
cov_ind_gw <- hvie_sylbor$cov_ind
cov_ind_bc <- hvie_sylatr$cov_ind


# Histoire de capture des sites 
#-----------------------------
# On remplace les histoires de vies des sites en 0 et 1
# quand le site est actif = 1, pas actif = 0

hvie_site <- as.matrix(hvie_ID_PROG_syl[1:K])
hvie_site[which(hvie_site>0)] <- 1


# Nb individu partageant la m?me observation
#--------------------------------
nb_ind_gw <- hvie_sylbor$nb_ind_groupe
nb_ind_bc <- hvie_sylatr$nb_ind_groupe

# Index 
#-------------------------
load('index_sylbor_point.RData')
load('index_sylatr_point.RData')

#mean
counts_gw <- c(NA,index_sylbor_point[[1]])
counts_bc <- c(NA,index_sylatr_point[[1]])

#sd
sd.counts_gw <- c(NA,sqrt(index_sylbor_point[[2]]))
sd.counts_bc <- c(NA,sqrt(index_sylatr_point[[2]]))


#  Dimensions   
#------------------
# Nb of individuals
N_gw <- dim(mydata_gw)[1]
N_bc <- dim(mydata_bc)[1]


# Compute the date of first capture for each individual:
#-------------------------------
f_gw <- NULL
for (i in 1:N_gw){
  temp <- 1:K
  f_gw <- c(f_gw,min(temp[mydata_gw[i,]==1]))}

f_bc <- NULL
for (i in 1:N_bc){
  temp <- 1:K
  f_bc <- c(f_bc,min(temp[mydata_bc[i,]==1]))}

# Compute the last date of capture for each individual with censure:
#-------------------------------
l_gw <- NULL
for (i in 1:N_gw){
  temp <- 1:K
  if (hvie_sylbor$censure[i]==1){l_new <- K}
  else {l_new <- max(temp[mydata_gw[i,]==1])}
  l_gw <- c(l_gw,l_new)
  rm(l_new)}

l_bc <- NULL
for (i in 1:N_bc){
  temp <- 1:K
  if (hvie_sylatr$censure[i]==1){l_new <- K}
  else {l_new <- max(temp[mydata_bc[i,]==1])}
  l_bc <- c(l_bc,l_new)
  rm(l_new)}

# Covariable age 
#------------------------
# un individu ne peut etre juvenile que la premiere annee de capture
# Donc apres la premiere occasion de capture f, on a forcement que des 2 (=adulte). 

for (i in 1:N_gw){
  for(j in 1:K){
    if (j > f_gw[i] & j<=l_gw[i]) {cov_age_gw[i,j] <- 2}
    else{next}
  }
}

for (i in 1:N_bc){
  for(j in 1:K){
    if (j > f_bc[i] & j<=l_bc[i]) {cov_age_bc[i,j] <- 2}
    else{next}
  }
}


cov_age_gw[which(cov_age_gw==0)] <- NA
cov_age_bc[which(cov_age_bc==0)] <- NA


rm(index_sylbor_point,hvie_ID_PROG_syl,hvie_sylbor,
   index_sylatr_point,hvie_sylatr,
   i,j,temp)

################################
# MODEL
################################

code <- nimbleCode({
  
  
  ###############
  #      GW     #
  ############### 
  
  # Initial population sizes
  # -------------------------------
  
  # Juveniles
  nN1_gw ~ dnorm(3900, sd = 1300)
  N1_gw[1] <- round(nN1_gw)
  
  # Adults
  nNad_gw ~ dnorm(1700, sd = 600)
  Nad_gw[1] <- round(nNad_gw)
  
  #Total
  Ntot_gw[1] <- N1_gw[1] + Nad_gw[1]
  
  
  # Productivity
  # -------------------------------
  
  # intercept
  mu.fec_gw ~ dnorm(log(4), sd= 0.5)
  sigma.fec_gw ~ dunif(0,1)
  # res
  for(t in 2:K){
    eps.fec_gw[t] ~ dnorm(0, sd = sigma.fec_gw)
    log(fec_gw[t]) <- mu.fec_gw + (DDIA_fec_gw * (Nad_gw[t]/1400)) + (DDIE_fec_gw * (Nad_bc[t]/4500)) + (Ddp_fec_gw * dp[t]) + (Dint_fec_gw * (dp[t] * (Nad_bc[t]/4500))) + eps.fec_gw[t]
  }#t
  
  # Survival / Recapture
  # -------------------------------
  
  # -- Survival
  # selon age
  
  for(u in 1:2) {
    # mean 
    mean.phi_gw[u] <- exp(mu.phi_gw[u]) /(1 + exp(mu.phi_gw[u]))
    mu.phi_gw[u] ~ dnorm(0, sd = 1)
    # sd pour effet aleatoire temps
    sigma.phi_gw[u] ~ dunif(0,5)
  }#u
  
  # --Recapture
  # Effet cov ind recapture
  gamma.i.p_gw   ~  dnorm(0, sd = 1)
  # mean
  mean.p_gw <- exp(mu.p_gw) /(1 + exp(mu.p_gw)) 
  mu.p_gw ~ dnorm(0,sd=1)
  # sd pour effet aleatoire temps
  sigma.p_gw ~ dunif(0,5)
  
  for (t in 1:K){ 
    # Survival
    logit(eta.phi_gw[1,t]) <- mu.phi_gw[1] + (DDIA_phiJ_gw * (Ntot_gw[t]/5000)) + (DDIE_phiJ_gw * (Ntot_bc[t]/15000)) + (Ddp_phiJ_gw * dp[t]) + (Dint_phiJ_gw * (dp[t] * (Ntot_bc[t]/15000))) + eps.phi_gw[1,t]
    logit(eta.phi_gw[2,t]) <- mu.phi_gw[2] + eps.phi_gw[2,t]
    #Random effect time
    eps.phi_gw[1,t] ~ dnorm(0, sd = sigma.phi_gw[1])
    eps.phi_gw[2,t] ~ dnorm(0, sd = sigma.phi_gw[2])
    
    # Recapture
    eta.p_gw[t] <- mu.p_gw + eps.p_gw[t]
    # Random effect time
    eps.p_gw[t] ~ dnorm(0, sd = sigma.p_gw)
    
  }#t
  
  # Loop on individuals
  
  for (i in 1:N_gw){
    for (t in f_gw[i]:l_gw[i]){ 
      # Survie
      phi_gw[i,t] <- eta.phi_gw[cov_age_gw[i,t],t]
      p_gw[i,t]   <-  1 / (1 + exp(-(eta.p_gw[t] + gamma.i.p_gw * cov_ind_gw[i]))) * (1-step(-hvie_site[ind_site_gw[i],t])) 
    } #t
  } #i 
  
  
  # System process
  # ------------------------------- 
  
  for (t in 2:K) {
    Nad_juv_gw[t] ~ dbin(eta.phi_gw[1,t-1],round(N1_gw[t-1]))
    Nad_ad_gw[t] ~ dbin(eta.phi_gw[2,t-1],round(Nad_gw[t-1]))
    
    # Nombre d'adultes qui sont deja sur le site
    Nres_gw[t] <- Nad_juv_gw[t] + Nad_ad_gw[t] 
    
    Nim_gw[t] <- round((Nad_juv_gw[t] + Nad_ad_gw[t])*0.3)
    
    # Ajout du nombre d'immigrant
    Nad_gw[t] <- Nres_gw[t] + Nim_gw[t]
    
    # Jeunes
    meanN1_gw[t] <- round(Nad_gw[t]*fec_gw[t]*0.5)
    N1_gw[t] ~ dpois(meanN1_gw[t])
    
    # Total
    Ntot_gw[t] <- N1_gw[t] + Nad_gw[t]
  } #t
  
  
  # Count likelihood
  # ------------------------------- 
  for (t in 2:K){
    logNad_gw[t] <- log(Nad_gw[t])
    counts_gw[t] ~ dnorm(logNad_gw[t], sd = sd.counts_gw[t])
  } #t
  
  # CR likelihood
  # ------------------------------- 
  for (i in 1:N_gw) {
    mydata_gw[i,f_gw[i]:l_gw[i]] ~ dCJS_vv_sum(probSurvive = phi_gw[i, f_gw[i]:l_gw[i]],
                                               probCapture = p_gw[i, f_gw[i]:l_gw[i]],
                                               mult = nb_ind_gw[i],
                                               len = l_gw[i]-f_gw[i]+1)
  } #i 
  
  
  ###############
  #      BC     #
  ############### 
  
  # Initial population sizes
  # -------------------------------
  
  #Juveniles
  nN1_bc ~ dnorm(9000, sd = 3000)
  N1_bc[1] <- round(nN1_bc)
  
  # Adults
  nNad_bc ~ dnorm(4000, sd = 1300)
  Nad_bc[1] <- round(nNad_bc)
  
  #Total
  Ntot_bc[1] <- N1_bc[1] + Nad_bc[1]
  
  # Productivity
  # -------------------------------
  
  # intercept
  mu.fec_bc ~ dnorm(log(4), sd = 0.5)
  sigma.fec_bc ~ dunif(0,1)
  
  # res
  for(t in 2:K){
    eps.fec_bc[t] ~ dnorm(0, sd = sigma.fec_bc)
    log(fec_bc[t]) <-  mu.fec_bc + (DDIA_fec_bc * (Nad_bc[t]/4500)) + ((Ddp1_fec_bc * dp_bc[t]) + (Ddp2_fec_bc * dp2_bc[t])) + eps.fec_bc[t]
  }#t
  
  # Survival / Recapture
  # -------------------------------
  
  # -- Survival
  # selon age
  
  for(u in 1:2) {
    # mean 
    mean.phi_bc[u] <- exp(mu.phi_bc[u]) /(1 + exp(mu.phi_bc[u]))
    mu.phi_bc[u] ~ dnorm(0,sd=1)
    # sd pour effet aleatoire temps
    sigma.phi_bc[u] ~ dunif(0,5)
  }#u
  
  # --Recapture
  # Effet cov ind recapture
  gamma.i.p_bc   ~  dnorm(0,sd=1)
  # mean
  mean.p_bc <- exp(mu.p_bc) /(1 + exp(mu.p_bc)) 
  mu.p_bc ~ dnorm(0,sd=1)
  # sd pour effet aleatoire temps
  sigma.p_bc ~ dunif(0,5)
  
  for (t in 1:K){ 
    # Survival
    logit(eta.phi_bc[1,t]) <- mu.phi_bc[1] + (DDIA_phiJ_bc * (Ntot_bc[t]/15000)) + ((Ddp1_phiJ_bc * dp_bc[t]) + (Ddp2_phiJ_bc * dp2_bc[t])) + eps.phi_bc[1,t]
    logit(eta.phi_bc[2,t]) <- mu.phi_bc[2] + eps.phi_bc[2,t]
    #Random effect time
    eps.phi_bc[1,t] ~ dnorm(0, sd = sigma.phi_bc[1])
    eps.phi_bc[2,t] ~ dnorm(0, sd = sigma.phi_bc[2])
    
    # Recapture
    eta.p_bc[t] <- mu.p_bc + eps.p_bc[t]
    # Random effect time
    eps.p_bc[t] ~ dnorm(0, sd = sigma.p_bc)
    
  }#t
  
  # Loop on individuals
  
  for (i in 1:N_bc){
    for (t in f_bc[i]:l_bc[i]){ 
      # Survie
      phi_bc[i,t] <- eta.phi_bc[cov_age_bc[i,t],t]
      p_bc[i,t]   <-  1 / (1 + exp(-(eta.p_bc[t] + gamma.i.p_bc * cov_ind_bc[i]))) * (1-step(-hvie_site[ind_site_bc[i],t])) 
    } #t
  } #i 
  
  # System process
  # ------------------------------- 
  
  for (t in 2:K) {
    Nad_juv_bc[t] ~ dbin(eta.phi_bc[1,t-1],round(N1_bc[t-1]))
    Nad_ad_bc[t] ~ dbin(eta.phi_bc[2,t-1],round(Nad_bc[t-1]))
    
    # Nombre d'adultes qui sont deja sur le site
    Nres_bc[t] <- Nad_juv_bc[t] + Nad_ad_bc[t] 
    
    Nim_bc[t] <- round((Nad_juv_bc[t] + Nad_ad_bc[t])*0.5)
    
    # Ajout du nombre d'immigrant
    Nad_bc[t] <- Nres_bc[t] + Nim_bc[t]
    
    # Jeunes
    meanN1_bc[t] <- round(Nad_bc[t]*fec_bc[t]*0.5)
    N1_bc[t] ~ dpois(meanN1_bc[t])
    
    # Total
    Ntot_bc[t] <- N1_bc[t] + Nad_bc[t]
  } #t
  
  
  # Count likelihood
  # ------------------------------- 
  for (t in 2:K){
    logNad_bc[t] <- log(Nad_bc[t])
    counts_bc[t] ~ dnorm(logNad_bc[t], sd = sd.counts_bc[t])
  } #t
  
  # CR likelihood
  # ------------------------------- 
  for (i in 1:N_bc) {
    mydata_bc[i,f_bc[i]:l_bc[i]] ~ dCJS_vv_sum(probSurvive = phi_bc[i, f_bc[i]:l_bc[i]],
                                               probCapture = p_bc[i, f_bc[i]:l_bc[i]],
                                               mult = nb_ind_bc[i],
                                               len = l_bc[i]-f_bc[i]+1)
  } #i 
  
  
  ###############
  #  All syl    #
  ############### 
  
  for(t in 1:K){
    Nad_ts[t] <- Nad_gw[t] + Nad_bc[t]
    Ntot_ts[t] <- Ntot_gw[t] + Ntot_bc[t]
  }#t
  
  
  Dint_fec_gw ~ dnorm(0, sd = 10)
  Dint_phiJ_gw ~ dnorm(0, sd = 10)
  Ddp_fec_gw ~ dnorm(0, sd = 10)
  Ddp_phiJ_gw ~ dnorm(0, sd = 10)
  
  Ddp1_fec_bc ~ dnorm(0, sd = 10)
  Ddp1_phiJ_bc ~ dnorm(0, sd = 10)
  Ddp2_fec_bc ~ dnorm(0, sd = 10)
  Ddp2_phiJ_bc ~ dnorm(0, sd = 10)
  
  DDIA_phiJ_bc ~ dnorm(0, sd = 2)
  DDIA_fec_bc  ~ dnorm(0, sd = 0.5)
  DDIE_phiJ_gw~ dnorm(0, sd = 2)
  DDIE_fec_gw  ~ dnorm(0, sd = 0.5)
  DDIA_phiJ_gw ~ dnorm(0, sd = 2)
  DDIA_fec_gw  ~ dnorm(0, sd = 0.5)
})

########
# DATA #
########

# Form the list of data
const = list(K=K,
             N_gw=N_gw,
             sd.counts_gw = sd.counts_gw,
             f_gw=f_gw,
             l_gw=l_gw,
             hvie_site=hvie_site,
             ind_site_gw = ind_site_gw,
             cov_age_gw=cov_age_gw,
             cov_ind_gw=cov_ind_gw,
             nb_ind_gw = nb_ind_gw,
             N_bc=N_bc,
             sd.counts_bc = sd.counts_bc,
             f_bc=f_bc,
             l_bc=l_bc,
             ind_site_bc = ind_site_bc,
             cov_age_bc=cov_age_bc,
             cov_ind_bc=cov_ind_bc,
             nb_ind_bc = nb_ind_bc,
             dp = dp,
             dp_bc = dp_bc,
             dp2_bc = (dp_bc*dp_bc))

data = list( mydata_gw = mydata_gw,
             counts_gw = counts_gw,
             mydata_bc = mydata_bc,
             counts_bc = counts_bc)

#########
# INITS #
#########

init   <- list(mu.phi_gw =c (0,0),
               mu.p_gw = 0,
               mu.fec_gw = log(4),
               sigma.phi_gw = c(0.2,0.2),
               sigma.p_gw = 0.2,
               sigma.fec_gw = 0.2,
               gamma.i.p_gw = 0,
               nN1_gw = 3900,
               nNad_gw = 1700,
               Nad_juv_gw = rep(200, 19),
               Nad_ad_gw = rep(300, 19),
               N1_gw = rep(3200,19),
               eps.fec_gw = rep(0,19),
               mu.phi_bc=c(-1.5,-1.5),
               mu.p_bc=-1.2,
               mu.fec_bc = log(4),
               sigma.phi_bc=c(0.2,0.2),
               sigma.p_bc=0.2,
               sigma.fec_bc = 0.2,
               gamma.i.p_bc=0,
               nN1_bc = 9000,
               nNad_bc = 4000,
               Nad_juv_bc = rep(500,19),
               Nad_ad_bc = rep(900,19),
               N1_bc = rep(11000,19),
               DDIA_phiJ_gw = 0,
               DDIA_fec_gw  = 0,
               DDIE_phiJ_gw = 0,
               DDIE_fec_gw  = 0,
               DDIA_phiJ_bc = 0,
               DDIA_fec_bc  = 0,
               Ddp_fec_gw = 0,
               Dint_fec_gw = 0,
               Ddp_phiJ_gw = 0,
               Dint_phiJ_gw = 0,
               Ddp1_fec_bc = 0,
               Ddp2_fec_bc = 0,
               Ddp1_phiJ_bc = 0,
               Ddp2_phiJ_bc = 0)

inits <- list(init,init)

# Specify the parameters to be monitored
parameters <- c("eta.phi_gw","eta.p_gw","sigma.phi_gw","eps.phi_gw",
                "mean.phi_gw","mu.phi_gw",
                "mean.p_gw"  ,"mu.p_gw"  ,"sigma.p_gw" ,"eps.p_gw",
                "gamma.i.p_gw","mu.fec_gw","eps.fec_gw","sigma.fec_gw",
                "Ntot_gw",
                "fec_gw","N1_gw", "nN1_gw", "nNad_gw","Nad_ad_gw","Nad_juv_gw","Nad_gw","Nim_gw",'Nres_gw',
                "eta.phi_bc","eta.p_bc","sigma.phi_bc","eps.phi_bc",
                "mean.phi_bc","mu.phi_bc",
                "mean.p_bc"  ,"mu.p_bc"  ,"sigma.p_bc" ,"eps.p_bc",
                "gamma.i.p_bc","mu.fec_bc","eps.fec_bc","sigma.fec_bc",
                "Ntot_bc",
                "fec_bc","N1_bc", "nN1_bc", "nNad_bc","Nad_ad_bc","Nad_juv_bc","Nad_bc","Nim_bc",'Nres_bc',
                'Nad_ts',"Ntot_ts",
                'DDIA_phiJ_gw','DDIA_fec_gw','DDIE_phiJ_gw','DDIE_fec_gw',
                'DDIA_phiJ_bc','DDIA_fec_bc',
                'Ddp_fec_gw','Dint_fec_gw','Ddp_phiJ_gw','Dint_phiJ_gw',
                'Ddp1_fec_bc','Ddp2_fec_bc','Ddp1_phiJ_bc','Ddp2_phiJ_bc')

rm(K, cov_age_bc,cov_age_gw,hvie_site,mydata_bc,mydata_gw, counts_bc, counts_gw, cov_ind_bc, cov_ind_gw, f_bc, f_gw, ind_site_bc, ind_site_gw, l_bc, l_gw, N_bc, N_gw, nb_ind_bc, nb_ind_gw, sd.counts_bc, sd.counts_gw,dp_bc,dp)

#################
# fit the model #
#################
m <- nimbleModel(code, const, data, init, check = FALSE, calculate = FALSE)
m$calculate()

# Inits simulation
simNodes <- c('eps.phi_gw','eps.p_gw',"eps.fec_gw", 
              'eps.phi_bc','eps.p_bc',"eps.fec_bc")

simNodeScalar <- m$expandNodeNames(simNodes)

allNodes <- m$getNodeNames()

nodesSorted <- allNodes[allNodes %in% simNodeScalar]

for(n in nodesSorted) {
  m$simulate(n)
  depNodes <- m$getDependencies(n)
  m$calculate(depNodes)
}

rm(depNodes,n,nodesSorted,allNodes,simNodes,simNodeScalar)

m$calculate()

n.thin = 20
n.burnin = 100000
n.keep.exact = 20000
n.keep = n.keep.exact * n.thin
n.iter = n.burnin + n.keep
n.chains = 2

conf <- configureMCMC(m, 
                      thin=n.thin, 
                      monitors = parameters, 
                      autoBlock = FALSE)

# Build model
Rmcmc<-buildMCMC(conf,enableWAIC = TRUE)

# Compile model
Cmodel<-compileNimble(m)

# Copile MCMC
Cmcmc<-compileNimble(Rmcmc, project=m)

# Run model
out <- runMCMC(Cmcmc, 
               niter = n.iter, 
               nburnin = n.burnin, 
               nchains = n.chains , 
               inits=inits,
               progressBar = TRUE, 
               summary = FALSE, 
               WAIC = TRUE,
               samplesAsCodaMCMC = TRUE) 

save(out,file='out_IPM_syl.RData')

