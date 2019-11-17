
#------------------------------------------------------#
#---# sampling of the populations generated before #---#
#------------------------------------------------------#

# once we have the artifitial populations forming clusters, we do the sampling 
# This function generates a sample of size 10% of the population size from a stratified plan 

echantillonnage <- function(dataset){  
  # simulates one single sample from the specified population with an stratified plan
  
  N <- dim(dataset)[1]            # population size
  p <- dim(dataset)[2]            # dimension od the data base 
  K <- length(table(dataset[,p])) # number of clusters 
  
  # in order to do the stratified sampling, we need the size of each sub-population 
  N_K <- rep(0,K)
  
  for(i in 1:K){
    N_K[i] <- table(dataset[,p])[i]
  }
  P <- c(1,cumsum(N_K))                  # size of P is K+1
  
  # sample size 
  n  <- N*0.1                            # we adopt this criteria, 10% of the total population   
  pi <- prop.table(table(dataset[,p]))   # proportions to do the stratified sampling 

  # size of the strats 
  n_K <- rep(0,K)
  
  for(i in 1:K){
    n_K[i] <- ceiling(n*pi[i]) # proportions of each cluster in the sample 
  }
  sum(n_K)
  
  # particular for K=3, but 
  if(sum(n_K) > n) n_K[K]   <- floor(n*pi[K]) # last cluster 
  if(sum(n_K) > n) n_K[K-1] <- floor(n*pi[K-1]) 
  if(sum(n_K) > n) n_K[K-2] <- floor(n*pi[K-2]) 
  if(K>3 && sum(n_K) > n) n_K[K-3] <- floor(n*pi[K-3]) 
  if(K>4 && sum(n_K) > n) n_K[K-4] <- floor(n*pi[K-4]) 
  if(K>5 && sum(n_K) > n) n_K[K-5] <- floor(n*pi[K-5]) 
  
  
  s_K <- vector("list", K)  
  for(i in 1:K){
    s_K[[i]] <- sample(x=seq(P[i],P[i+1]), size = n_K[i]) # we want a sample of size n_k from 
  }                                                       # from the k-th strate 
  
  s <- unlist(s_K)             # total sample labels 
  echantillon <- dataset[s,]   # data corresponding to the choosen sample 
  return(echantillon)
}

# now we want M = 100 samples from this same population
# the output of this function is a list, each list contains a different sample 

# M       -- number of wanted samples 
# dataset -- dataset from which we apply the sampling 

generer_echantillons <- function(M, dataset){
  # we generate M samples from the specified population 
  samples <- vector("list", M)
  for(i in 1:M){
    samples[[i]] <- echantillonnage(dataset)   
  }
  return(samples)  
}

M_samples_P2K3_HO <- generer_echantillons(M=100, dataset = population_P2K3_HO)
M_samples_P2K3_MO <- generer_echantillons(M=100, dataset = population_P2K3_MO)
M_samples_P2K3_LO <- generer_echantillons(M=100, dataset = population_P2K3_LO)

# save the samplings in a file 
# we are going to create a directory for each k and each p 

# you must have to change your working directory to use this function

write_file_cov <- function(database){ # the database must content the set of covariables 

M <- dim(summary(database))[1]  # number of samples 
p <- dim(database[[1]])[2]-1    # number of covariates
K <- length(table(database[[1]][,p+1])) # number of clusters
n_data <- dim(database[[1]])[1] # number of observations  

dir.samples <- "C:/Users/Documents/project/samples/"
dir.samples <-  paste( dir.samples,"p",p, "K", K, "/", sep="") # concatena 

for(i in 1:M){
  colnames(database[[i]]) <- c()
  data.vals <- (database[[i]][,-(p+1)])
  file.set.cov = paste( dir.samples, "CovSet", i, ".txt", sep="") # concatena 
  write( c(n_data, p), file= file.set.cov, ncolumns=1 )   # The data (usually a matrix) x are written to file file
  write( t( data.vals ), file= file.set.cov, ncolumns=p, append =TRUE ) # el 2 en c(Ndata,2) quiere decir dos covariables? 
  }
}

write_file_cov(M_samples_P2K3_HO)
write_file_cov(M_samples_P2K3_MO)
write_file_cov(M_samples_P2K3_LO)

# example d'échantillon 

e1_P2K3_HO <- M_samples_P2K3_HO[[1]]
e1_P2K3_MO <- M_samples_P2K3_MO[[1]]
e1_P2K3_LO <- M_samples_P2K3_LO[[1]]

colors <- c("tomato", "lightpink", "lightskyblue", "red", "magenta", 
            "blue", "brown", "green", "Skyblue", "orange", "yellow","pink", "Violet")

x11()
par(mfrow=c(3,1), mar= c(1, 2, 2, 2) )
plot(e1_P2K3_HO[,1:2], col = colors[e1_P2K3_HO[,3]],
     pch=19, cex =0.8, xlab="", ylab="",
     main="chevauchement élevé",  xaxt='n',  yaxt='n')

plot(e1_P2K3_MO[,1:2], col = colors[e1_P2K3_MO[,3]],
     pch=19, cex =0.8, xlab="", ylab="",
     main=" chevauchement modéré",  xaxt='n',  yaxt='n')

plot(e1_P2K3_LO[,1:2], col = colors[e1_P2K3_LO[,3]],
     pch=19, cex =0.8, xlab="", ylab="",
     main="faible chevauchement",  xaxt='n',  yaxt='n')
