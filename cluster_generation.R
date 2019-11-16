
library("MixSim")
library("broman")

#----------------------------------------------------------------------------------#
#---# génération des populations artifitielles issues d'une mélange gaussienne #---#
#----------------------------------------------------------------------------------#

# On a généré artificiellement différentes populations de grande taille issues d’un mélange
# de lois gaussiennes, puis différents hyper-paramètres pour cette simulation ont été considérés.

# p_dim   -- dimension du vecteur de paramètres de régression, (beta)
# K_clust -- nombre de clusters à générer 
# lam     -- paramètre de la transformation de Cox (doit être entre -5 et 5) 

conditions <- function(p_dim, K_clust, lam){   
  x <- 0
  if( -5 > lam || lam > 5 ){
    print("lambda must be between -5 and 5")                   # paramètre à utiliser dans la transformation de Cox 
    x <- x + 1
  }
  if( p_dim!=2 && p_dim!=5 && p_dim!=10){        # pour cette experience on va étudier les cas où la 
    print("choose other value of p")                           # dimension de beta est 2, 5, 10 ou 15  
    x <- x+1 
  }
  if( K_clust!=3 && K_clust!=5 && K_clust!=7 && K_clust!=10 ){ # le nombre de clusters génerées est fixé à 3, 5 7 ou bien 10
    print("choose other value of K")  
    x <- x+1
  }
  if(x!=0){
    return(FALSE)
  } else{
    return(TRUE)
  }
} # conditions pour la géneration de mélanges 

# on essaie d'éviter d'avoir des consentration de données isolées,
# surtout quand le nombre de clusteres, K_clust est grand 
# M       - matrice de chevauchement 
# K_clust - nombre de clusters géneré 
not_isolated_clusters <- function(M, K_clust){ 
  M_ind <- matrix(0, nrow=K_clust, ncol=K_clust)
  for(i in 1 : K_clust){
    for( j in 1 : K_clust ){
      if( M[i,j]==0 ){ 
        M_ind[i,j] <- 1 
      }
    }
  }
  return(sum(M_ind))
}

# cette fonction a pour but de construire un ensemble de données simulées 
# à partir d'un mélange généré avec les paramètres spécifiés

# le paramètre lam est celui utilisé dans la transformation Box-Cox 
# (si on veut des données autres que la normale), il doit être comprise entre 
# -5 et 5 )
# p_dim        - dimension de beta
# K_clust      - nombre de clusters à generer 
# mean_overlap - parametre de chevaucheent entre clusters
# N_size       - taille de la base de données à generer à partir du modèle de mélange 
# lam          - parametre utilisé dans la transformation de Cox 

# cette fonction est propre aux besoins de cette experience, mais on presente par la suite 
# un fonction plus generale 
MixData <- function(p_dim, K_clust, mean_overlap, N_size, lam){
  if(conditions(p_dim, K_clust, lam)==FALSE){
    stop( " conditions in the parameters not satisfied " ) 
  } else {
    pi_min  <- 1/K_clust- 0.05 # proportion des donnés au sein de chaque cluster 
                               # nous ne voulons pas qu'elles soient exactement égaux,
                               # mais assez semblables pour que nos donnés restent équilibrés.
    
    if(K_clust == 3){ 
      max_ovelp <- mean_overlap*(K_clust*(K_clust-1)/2) - (mean_overlap/2)
      
      # chevauchement élevé (10 fois plus grand par rapport au moderé)
      Mix_HO <- MixSim(BarOmega = mean_overlap*10, MaxOmega = max_ovelp*10, K = K_clust, p = p_dim, PiLow = pi_min, resN=1000) 
      
      # chevauchement modèrée 
      Mix_MO <- MixSim(BarOmega = mean_overlap,    MaxOmega = max_ovelp,    K = K_clust, p = p_dim, PiLow = pi_min) 
    }      # critéres de chevaucheant pour 3 clusters 
    
    if(K_clust == 5){
      repeat{
        
        max_ovelp <- mean_overlap*(K_clust*(K_clust-1)/2) - (12/2)*mean_overlap
        
        Mix_HO    <- MixSim(BarOmega = mean_overlap*(10/2),  MaxOmega = max_ovelp*(15/2), K = K_clust, p = p_dim, PiLow = pi_min, resN=1000) 
        Mix_MO    <- MixSim(BarOmega = mean_overlap,         MaxOmega = max_ovelp,        K = K_clust, p = p_dim, PiLow = pi_min, resN=1000) 
        
        # a partir de k = 5 on s'asure qu'il n'ai pas de clusters isolées 
        # le parametre de chevauchement est sensible au nombre de clusters 
        # verification  
        
        if( not_isolated_clusters(M = Mix_MO$OmegaMap, K_clust = K_clust) > (K_clust^2-K_clust)/2){
          mean_overlap <- mean_overlap + 0.001  # si la condition n'est pas satisfaite, on augmente le 
        } else {                                # chevaichement moyen 
          break
        }
      }
    }        # critéres de chevaucheant pour 5 clusters
    
    if(K_clust == 7){
      repeat{
        max_ovelp <- mean_overlap*(K_clust*(K_clust-1)/2) - (20/2)*mean_overlap
        
        Mix_HO <- MixSim(BarOmega = mean_overlap*(6/2),  MaxOmega = max_ovelp*(4/2), K = K_clust, p = p_dim, PiLow = pi_min, resN=1000) 
        
        Mix_MO <- MixSim(BarOmega = mean_overlap,        MaxOmega = max_ovelp,       K = K_clust, p = p_dim, PiLow = pi_min, resN=1000) 
        
        # verification  
        if(not_isolated_clusters(M = Mix_MO$OmegaMap, K_clust) > (K_clust^2-K_clust)/2){
          mean_overlap <- mean_overlap + 0.001
        } else {
          break
        }
      }
    }        # critéres de chevaucheant pour 7 clusters
    
    if(K_clust == 10){
      repeat{
        max_ovelp <- mean_overlap*(K_clust*(K_clust-1)/2) - (50/2)*mean_overlap
        
        Mix_HO <- MixSim(BarOmega = mean_overlap*(5/2),  MaxOmega = max_ovelp*(5/4), K = K_clust, p = p_dim, PiLow = pi_min, resN=1000) 
        Mix_MO <- MixSim(BarOmega = mean_overlap*(3/2),  MaxOmega = max_ovelp,       K = K_clust, p = p_dim, PiLow = pi_min, resN=1000) 
        
        # verification 
        if(not_isolated_clusters(M = Mix_MO$OmegaMap, K_clust) > (K_clust^2-K_clust)/2){
          mean_overlap <- mean_overlap + 0.001
        } else {
          break
        }
      }
    }       # critéres de chevaucheant pour 10 clusters
    
    #formation des données a partir du modèle de mélange géneré avant  
    dataset_HO <- simdataset(n = N_size, Pi = Mix_HO$Pi, Mu = Mix_HO$Mu, S = Mix_HO$S, lambda = rep(lam, p_dim)) 
    dataset_MO <- simdataset(n = N_size, Pi = Mix_MO$Pi, Mu = Mix_MO$Mu, S = Mix_MO$S, lambda = rep(lam, p_dim))
    
    # si la dimention du vecteur beta est deux, on peut effectuer une vérification graphique
    if(p_dim==2){
      colors <- c("red", "magenta", "blue", "brown", "green",
              "Skyblue", "orange", "yellow","pink", "Violet")
      
      x11()
      par(mfrow=c(2,1), mar= c(1, 2, 2, 2) )
      plot(dataset_HO$X, col = colors[dataset_HO$id], pch=19, cex =0.8,
           xlab="", ylab="", main="High overlap",  xaxt='n',  yaxt='n')
      
      plot(dataset_MO$X, col = colors[dataset_MO$id], pch=19, cex =0.8,
           xlab="", ylab="", main="Moderated overlap",  xaxt='n', yaxt='n')
    }
    
    # creation du data frame avec la donnée géneré 
    dataset_H <- data.frame(dataset_HO$X, dataset_HO$id)
    dataset_M <- data.frame(dataset_MO$X, dataset_MO$id)
    
    # creation des labels 
    x_p <- rep(0, p_dim)
    
    for(i in 1:p_dim){
      x_p[i] <- paste("x", i, sep = "_")
    }
    
    colnames(dataset_H) <- colnames(dataset_M) <- c(x_p, "label")
    
    mix_data <- vector("list", 2)
    mix_data[[1]] <- dataset_H
    mix_data[[2]] <- dataset_M
    
    return(mix_data)  # la fonction returne la base de données
  }
}


# on n'impose pas des conditions de verification, mais on s'attend à que 
# les hyperparametres rentrés par l'utilisateur sont appropriés 

MixData_generation <- function(p_dim, K_clust, mean_overlap, N_size, lam){
  
  # proportion des donnés au sein de chaque cluster 
  pi_min <- 1/K_clust- 0.05 
  
  # paramètre de chevauchement maximale 
  max_ovelp <- mean_overlap*(K_clust*(K_clust-1)/2) - (mean_overlap/2)
  
  # generation des mélanges 
  # chevauchement élevé 
  Mix_HO <- MixSim(BarOmega = mean_overlap*10, MaxOmega = max_ovelp*10,
                   K = K_clust, p = p_dim, PiLow = pi_min, resN=1000) 
  
  # chevauchement modèrée 
  Mix_MO <- MixSim(BarOmega = mean_overlap,    MaxOmega = max_ovelp, 
                   K = K_clust, p = p_dim, PiLow = pi_min) 
  
  #formation des données a partir du modèle de mélange géneré avant  
  # chevauchement élevé 
  dataset_HO <- simdataset(n = N_size, Pi = Mix_HO$Pi, Mu = Mix_HO$Mu,
                           S = Mix_HO$S, lambda = rep(lam, p_dim)) 
  # chevauchement modèrée 
  dataset_MO <- simdataset(n = N_size, Pi = Mix_MO$Pi, Mu = Mix_MO$Mu,
                           S = Mix_MO$S, lambda = rep(lam, p_dim))
  
  # si la dimention du vecteur beta est deux, on peut effectuer une validation graphique
  if(p_dim==2){
    colors <- c("tomato", "lightpink", "lightskyblue", "red", "magenta", 
                "blue", "brown", "green", "Skyblue", "orange", "yellow","pink", "Violet")
    
    x11()
    par(mfrow=c(2,1), mar= c(1, 2, 2, 2) )
    plot(dataset_HO$X, col = colors[dataset_HO$id], pch=19, cex =0.8,
         xlab="", ylab="", main="High overlap",  xaxt='n',  yaxt='n')
    
    plot(dataset_MO$X, col = colors[dataset_MO$id], pch=19, cex =0.8,
         xlab="", ylab="", main="Moderated overlap",  xaxt='n', yaxt='n')
  }
  
  # creation du data frame avec la donnée géneré 
  dataset_H <- data.frame(dataset_HO$X, dataset_HO$id)
  dataset_M <- data.frame(dataset_MO$X, dataset_MO$id)
  
  # creation des labels 
  x_p <- rep(0, p_dim)
  
  for(i in 1:p_dim){
    x_p[i] <- paste("x", i, sep = "_")
  }
  
  colnames(dataset_H) <- colnames(dataset_M) <- c(x_p, "label")
  
  mix_data <- vector("list", 2)
  mix_data[[1]] <- dataset_H
  mix_data[[2]] <- dataset_M
  
  return(mix_data)  # la fonction returne une liste avec les deux bases de données
}


# Test de la fonction 
P2_K3_datasets <- MixData_generation(p_dim = 2, K_clust = 3, mean_overlap = 0.005,
                                     N_size = 10000, lam = 0.5)

population_P2_K3_HO <- P2_K3_datasets[[1]] # base de donnés avec chevauchement elevé 
population_P2_K3_MO <- P2_K3_datasets[[2]] # base de donnés avec chevauchement moderé

head(population_P2_K3_HO)
head(population_P2_K3_MO)

