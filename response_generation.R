
library("graphics")                    
library("xtable")

#-------------------------------------------#
#---# Génération de la variable réponse #---#
#-------------------------------------------#

# For each sample, we have to estimate "y" with the piecewise cox model, 
# in order to do that, we need to establish prior values for parameters 
# alpha, lambda and beta. 

#----------------------------------------------#
#---# 1. fonction de risque cumulé de base #---#
#----------------------------------------------#

# j -- indique l'intervalle I_j tel que 'y' est en I_j
# y -- temps de défaillance observé 
# S -- intervalles pour partitionner l'axe du temps
# alpha et lambda -- vecteurs associés à la fonction de risque de base 
# alpha = (alpha_1, ..., alpha_J)

H <- function(j, y, S, alpha, lambda){  
  if(j ==1){                            
      H_out <- lambda[j]*(y^(alpha[j]) - S[j]^(alpha[j]))
  } else {
      H_sum <- 0
        for(g in seq(1,(j-1))){
          H_sum <- H_sum + lambda[g]*(S[g+1]^(alpha[g]) - S[g]^(alpha[g]))
        }
    H_out <- lambda[j]*(y^(alpha[j]) - S[j]^(alpha[j])) + H_sum
  }
  return(H_out)
}

#------------------------------------#
#---# 2. fonction de répartition #---#
#------------------------------------#

# constant de normalisation 

normalizedF <- function(S, alpha, lambda, u) {
  GI <- length(S)
  w  <- 1 - exp(-H(GI - 1, S[GI], S, alpha, lambda) * u)
  return(w)
}

F_dist <- function(j, y, S, alpha, lambda, x, beta){
  w <- (1-exp(-H(j, y, S, alpha, lambda)*exp(x%*%beta)))/normalizedF(S, alpha, lambda, exp(x%*%beta)) 
  return(w)
}

#--------------------------------------------#
#---# 3. fonction de répartition inverse #---#
#--------------------------------------------#

# fonction de repartition inverse 

H_inv <- function(j, alpha, lambda, S){
  if(j == 1){
    if(S[1] > 0){
      H_i <- lambda[1]*exp(alpha[1]*log(S[1]))
    } else {
      H_i <- 0
    }
  } else {
    H_sum <- 0
    for(g in seq(1, (j-1))){
      H_sum <- H_sum + lambda[g]*(S[g+1]^alpha[g] - S[g]^alpha[g])
    }
    H_i <- lambda[j]*(S[j]^alpha[j]) - H_sum 
  }
  return(H_i)
}

F_inv <- function(j, alpha, lambda, S, beta, x, U){
  
  uN   <- normalizedF( S, alpha, lambda, exp( x%*% beta) ) 
  logU <- log(1 - U *uN )
  Hj   <- H_inv(j, alpha, lambda, S)*exp(x%*%beta)
  #cat("uN=", uN, " log=", log(uN), " j=", j, " Hj=", Hj, "\n")
  F_i  <- ((Hj-logU)/(lambda[j]*exp(x%*%beta)))^(1/alpha[j])
  
  return(F_i)
}

#-----------------------#
#---# 4. Simuler y  #---#
#-----------------------#

# simulate one single "y"

simulation_one <- function(alpha, lambda, beta, S, x){
  U <- runif(1)
  J <- length(alpha[[1]])
  F_d <- rep(0,J)
  for(j in seq(1,J)){
    F_d[j] <- F_dist(j, S[j+1], S, alpha, lambda, x, beta) 
  }
  #normalized, make sure last term is 1
  F_d <- c(0, F_d) /F_d[J]
  #print(F_d)
  y   <- 0
  if(J > 1){
    for(i in seq(1,J-1)){
      #cat("U=", U, " F[", i, "]=", F_d[i], ", ", F_d[i+1], "\n")
      if(F_d[i] <= U & U < F_d[i+1]){
        #cat("Found between Fd[", i, "]=", F_d[i]," and Fd[", i+1, "]=", F_d[i+1], "\n")
        y <- F_inv(i, alpha, lambda, S, beta, x, U)
      } 
    }
  } else {
    y <- F_inv(1, alpha, lambda, S, beta, x, U)
  }
  return(y)
} 

# gives the y generated from the piecewise weibull model for the specified parameters 

# database_cov -- échantillons des covariables génerés dans le script "sample_generation" 
# alpha.os     -- valeurs de départ de alpha
# lambda.os    -- valeurs de départ de lambda 
# beta.os      -- valeurs de départ de beta
# S.os         -- bornes des intervalles 

generer_y <- function(database_cov, alpha.os, lambda.os, beta.os, S.os){  # the input is a list of M=100 samples from the objectif population 
  
  M <- dim(summary(database_cov))[1]  # number of samples from the population 
  p <- dim(database_cov[[1]])[2]-1    # number of covariates 
  K <- length(table(database_cov[[1]][,p+1])) # number of clusters generated 
  n_data <- dim(database_cov[[1]])[1] # number of observations in each sample 
  
  database_y <- vector("list", M )
  
  y_sim <- rep(0, n_data)
  
  for(j in 1:M){ database_y[[j]] <- y_sim } 
  
  for(j in 1:M){
    for(i in 1: n_data){
      
      K <- database_cov[[j]]$label[i]
      database_y[[j]][i] <- simulation_one(alpha = alpha.os[[K]], lambda = lambda.os[[K]], beta = beta.os[[K]], S = S.os, x = as.matrix(database_cov[[j]][i,-(p+1)]))
      
      if(database_y[[j]][i] > S.os[J+1]){ database_y[[j]][i] <- S.os[J+1] }
    }
  }
  return(database_y)
}

#----------------------------------------------#
#---# 5. tester différents hyperparamètres #---#
#----------------------------------------------#

# configuration :  p = 2 & k =3 

J <- 2 # number of interval
K <- 3 # the number of clusters
S_J2 <- c(0, 20, 40) # intervals of the time axis   
nu1 <- 1/20
nu2 <- 1/20

lambda_J2 <- vector("list", K) # values of lambda for the K = 3 simulated clusters 
alpha_J2  <- vector("list", K)
beta_p2   <- vector("list", K)

# case 1 : try similar values of beta and lambda (in each cluster) #
beta_p2[[1]] <- c(1, 2)      #runif(2, -0.5, 0.5)  # p=2 very similar values 
beta_p2[[2]] <- c(-0.5, 0.5) #runif(2, -0.5, 0.5) 
beta_p2[[3]] <- c(0, 1)      #runif(2, -0.5, 0.5) 

alpha_J2[[1]]  <- alpha_J2[[2]] <- alpha_J2[[3]] <- c(1, 1) 

# in order to determine a pausible value for lambda we will do the next 

# use the data set we generated before 

x <- split(x = population_P2K3_HO, f = population_P2K3_HO$label )
x_1 <- as.matrix(x[[1]][,-3])
x_2 <- as.matrix(x[[2]][,-3])
x_3 <- as.matrix(x[[3]][,-3])

beta_1 <- as.vector(beta_p2[[1]])   # sample, it should not change that much since the sampling is stratified
beta_2 <- as.vector(beta_p2[[2]]) 
beta_3 <- as.vector(beta_p2[[3]]) 

lambda_J2[[1]] <- round(c( (nu1  + 0.01) / mean(exp(x_1%*%beta_1)), 
                           (nu2 - 0.005) / mean(exp(x_1%*%beta_1))), 3) 

lambda_J2[[2]] <- round(c((nu1 + 0.005) / mean(exp(x_2%*%beta_2)), 
                          (nu2 - 0.010) / mean(exp(x_2%*%beta_2))), 3) 

lambda_J2[[3]] <- round(c( (nu1 + 0.015) / mean(exp(x_3%*%beta_3)), 
                           (nu2 - 0.015) / mean(exp(x_3%*%beta_3))), 3) 

# generate the response 

M_samples_y_P2K3_HO <- generer_y( database_cov = M_samples_P2K3_HO,
                                  alpha.os = alpha_J2, 
                                  lambda.os = lambda_J2, 
                                  beta.os = beta_p2,
                                  S = S_J2)

# genera la respuesta para las M muestras (es por eso que toma tiempo)

# save the parameters used 
M <- matrix(c(nu1, nu2, beta_1, mean(exp(x_1%*%beta_1)), lambda_J2[[1]], 
              nu1, nu2, beta_2, mean(exp(x_2%*%beta_2)), lambda_J2[[2]],
              nu1, nu2, beta_3, mean(exp(x_3%*%beta_3)), lambda_J2[[3]]), nrow = 3, byrow=T)
colnames(M) <- c("nu1", "nu2", "beta1", "beta2", "E", "lambda1", "lambda2")
rownames(M) <- c("C1", "C2", "C3")
xtable(M, digits = 3)

# you must have to change your working directory to use this function
# database -- object conteining the generated "y"
# p        -- dimension of the covariates 
# K        -- number of clusters 

write_file_y <- function(database, p, K){ # the database must content the response variables correponding to a set of covariates specified before  
  
  M <- dim(summary(database))[1]
  n_data <- length(database[[1]])
  
  dir.samples <- "C:/Users/Documents/project/"
  dir.samples <-  paste( dir.samples,"p",p, "K", K, "/", sep="") # concatena 
  
  for(i in 1:M){
    
    file.set.resp <- paste( dir.samples, "coxResponseJ2Set", i, ".txt", sep="") # concatena 
    dat.resp <- matrix(1, ncol=2, nrow= n_data )         # matriz llena de unos, dos columnas y Ndata renglones 
    dat.resp[,1] = database[[i]]                         # llenamos esta matriz con los datos simulados con ayuda de la rutina "simulation_one"
    
    write(c(n_data,p), file = file.set.resp, ncolumns=1 )   # The data (usually a matrix) x are written to file file
    write(t(dat.resp), file = file.set.resp, ncolumns=2, append =TRUE )
  }
} 

write_file_y(M_samples_y_P2K3_HO,  p = 2, K = 3)
