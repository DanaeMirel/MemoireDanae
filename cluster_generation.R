
library("MixSim")
library("broman")

#----------------------------------------------------------------------------------#
#---# génération des populations artifitielles issues d'une mélange gaussienne #---#
#----------------------------------------------------------------------------------#

# On a généré artificiellement différentes populations de grande taille issues d’un mélange
# de lois gaussiennes, puis différents hyper-paramètres pour cette simulation ont été considérés.

#-----------------------------------------------------#
#---# 1. conditions préalables pour la génération #---#
#-----------------------------------------------------#

# paramètres de la fonction "conditions"
# p_dim   -- dimension du vecteur de paramètres de régression, (beta)
# K_clust -- nombre de clusters à générer 
# lam     -- paramètre de la transformation de Cox (doit être entre -5 et 5) 

conditions <- function(p_dim, K_clust, lam){   
  x <- 0
  if( -5 > lam || lam > 5 ){
      print("lambda must be between -5 and 5")               # paramètre à utiliser dans la transformation de Cox 
      x <- x + 1
  }
  if( p_dim!=2 & p_dim!=5 & p_dim!=10){                      # pour cette experience on va étudier les cas où la 
      print("choose other value of p")                       # dimension de beta est 2, 5 ou 10  
      x <- x+1 
  }
  if( K_clust!=3 & K_clust!=6 & K_clust!=10 ){               # le nombre de clusters génerées est fixé à 3, 6 ou 10
      print("choose other value of K")  
      x <- x+1
  }
  if(x!=0){
    return(FALSE)
  } else{
    return(TRUE)
  }
} # conditions pour la géneration de mélanges 

# on essaie d'éviter des concentrations de données isolées, 
# surtout quand le nombre de clusters, "K_clust" est grand 

# paramètres de la fonction "not_isolated_clusters"

# M       -- matrice de chevauchement 
# K_clust -- nombre de clusters génerés 

not_isolated_clusters <- function(M, K_clust){ 
  M_ind <- matrix(0, nrow = K_clust, ncol = K_clust)
  for( i in 1 : K_clust ){
    for( j in 1 : K_clust ){
      if( M[i,j]==0 ){ 
        M_ind[i,j] <- 1 
      }
    }
  }
  return(sum(M_ind))
}

#-----------------------#
#---# 2. génération #---# 
#-----------------------#

# la fonction suivante a pour but de construire un ensemble de données simulées 
# à partir d'un mélange généré avec les paramètres spécifiés par l'utilisateur

# le paramètre lam est celui utilisé dans la transformation Box-Cox 
# (si on veut générer des données autres que de la normale)
# -5 < lam < 5

# p_dim        -- dimension de beta
# K_clust      -- nombre de clusters à générer 
# mean_overlap -- paramètre de chevauchement entre clusters
# N_size       -- taille de la base de données à générer à partir du modèle de mélange 
# lam          -- paramètre utilisé dans la transformation de Cox 

# cette fonction est propre aux besoins de cette experience, 
# mais les critères peuvent être changés 

MixData <- function(p_dim, K_clust, mean_overlap, N_size, lam){
  
  if(conditions(p_dim, K_clust, lam) == FALSE){
    stop( " conditions in the parameters not satisfied " ) 
  } else {
    
    pi_min  <- 1/K_clust - 0.05 # proportion des donnés au sein de chaque cluster 
                                # nous ne voulons pas qu'elles soient exactement égaux,
                                # mais assez semblables pour que nos jeux de données restent équilibrés.
    
    if(K_clust == 3){ 
      
      max_ovelp <- mean_overlap*(K_clust*(K_clust-1)/2) - (mean_overlap/2)
      
      # chevauchement élevée 
      Mix_HO <- MixSim(BarOmega = mean_overlap*10,
                       MaxOmega = max_ovelp*10, 
                       K = K_clust, 
                       p = p_dim, 
                       PiLow = pi_min, 
                       resN=1000) 
      
      # chevauchement modérée 
      Mix_MO <- MixSim(BarOmega = mean_overlap,    
                       MaxOmega = max_ovelp,    
                       K = K_clust, 
                       p = p_dim, 
                       PiLow = pi_min) 
    
      # faible chevauchement  
      Mix_LO <- MixSim(BarOmega = mean_overlap/10,    
                       MaxOmega = max_ovelp/10,    
                       K = K_clust, 
                       p = p_dim, 
                       PiLow = pi_min) 
      }      # critéres de chevaucheant pour 3 clusters 
    
    if(K_clust == 6){
      repeat{
        
        max_ovelp <- mean_overlap*(K_clust*(K_clust-1)/2) - (12/2)*mean_overlap
        
        Mix_HO <- MixSim(BarOmega = mean_overlap*(10/2),
                         MaxOmega = max_ovelp*(15/2),
                         K = K_clust,
                         p = p_dim,
                         PiLow = pi_min,
                         resN=1000) 
        
        Mix_MO <- MixSim(BarOmega = mean_overlap,
                         MaxOmega = max_ovelp, 
                         K = K_clust,
                         p = p_dim,
                         PiLow = pi_min, 
                         resN=1000) 
        
        Mix_LO <- MixSim(BarOmega = mean_overlap/10,
                         MaxOmega = max_ovelp/10,
                         K = K_clust,
                         p = p_dim,
                         PiLow = pi_min,
                         resN=1000) 
        
        # à partir de k = 5 on s'assure qu'il n'a pas de clusters isolés.
        # le paramètre de chevauchement est sensible au nombre de clusters vérification 
        
        if( not_isolated_clusters(M = Mix_MO$OmegaMap, K_clust = K_clust) > (K_clust^2-K_clust)/2){
          mean_overlap <- mean_overlap + 0.001  # si les conditions ne sont pas satisfaites, on augmente le chevauchement moyen
        } else {                               
          break
        }
      }
    } # critéres de chevaucheant pour 6 clusters
    
    if(K_clust == 10){
      repeat{
        
        max_ovelp <- mean_overlap*(K_clust*(K_clust-1)/2) - (50/2)*mean_overlap
        
        Mix_HO <- MixSim(BarOmega = mean_overlap*(6/2),
                         MaxOmega = max_ovelp*(4/2),
                         K = K_clust, 
                         p = p_dim,
                         PiLow = pi_min,
                         resN=1000) 
        
        Mix_MO <- MixSim(BarOmega = mean_overlap,
                         MaxOmega = max_ovelp, 
                         K = K_clust,
                         p = p_dim, 
                         PiLow = pi_min, 
                         resN=1000) 
        
        Mix_LO <- MixSim(BarOmega = mean_overlap/10,
                         MaxOmega = max_ovelp/10, 
                         K = K_clust,
                         p = p_dim, 
                         PiLow = pi_min, 
                         resN=1000) 
        
        # verification  
        if(not_isolated_clusters(M = Mix_MO$OmegaMap, K_clust) > (K_clust^2-K_clust)/2){
          mean_overlap <- mean_overlap + 0.001
        } else {
          break
        }
      }
    }        # critéres de chevaucheant pour 10 clusters
    
    # formation des données à partir du modèle de mélange généré dans l'étape précédente 
    
    dataset_HO <- simdataset(n = N_size, Pi = Mix_HO$Pi, Mu = Mix_HO$Mu, S = Mix_HO$S, lambda = rep(lam, p_dim)) 
    dataset_MO <- simdataset(n = N_size, Pi = Mix_MO$Pi, Mu = Mix_MO$Mu, S = Mix_MO$S, lambda = rep(lam, p_dim))
    dataset_LO <- simdataset(n = N_size, Pi = Mix_LO$Pi, Mu = Mix_LO$Mu, S = Mix_LO$S, lambda = rep(lam, p_dim))
    
    # si la dimention du vecteur beta est deux, on peut effectuer une vérification graphique
    if(p_dim == 2){
      colors <- c("red", "magenta", "blue", "brown", "green",
                  "Skyblue", "orange", "yellow","pink", "Violet")
      
      x11()
      par(mfrow=c(3,1), mar= c(1, 2, 2, 2) )
      plot(dataset_HO$X, col = colors[dataset_HO$id], pch=19, cex =0.8,
           xlab="", ylab="", main="High overlap", xaxt='n',  yaxt='n')
      
      plot(dataset_MO$X, col = colors[dataset_MO$id], pch=19, cex =0.8,
           xlab="", ylab="", main="Moderated overlap", xaxt='n', yaxt='n')
      
      plot(dataset_LO$X, col = colors[dataset_LO$id], pch=19, cex =0.8,
           xlab="", ylab="", main="Low overlap", xaxt='n', yaxt='n')
      
    }
    
    # Création du data frame avec la donnée générée 
    dataset_H <- data.frame(dataset_HO$X, dataset_HO$id)
    dataset_M <- data.frame(dataset_MO$X, dataset_MO$id)
    dataset_L <- data.frame(dataset_LO$X, dataset_LO$id)
    
    # creation des labels 
    x_p <- rep(0, p_dim)
    
    for(i in 1:p_dim){
      x_p[i] <- paste("x", i, sep = "_")
    }
    
    colnames(dataset_H) <- colnames(dataset_M) <- colnames(dataset_L) <- c(x_p, "label")
    
    mix_data <- vector("list", 3)
    mix_data[[1]] <- dataset_H
    mix_data[[2]] <- dataset_M
    mix_data[[3]] <- dataset_L
    
    return(mix_data)  #La fonction retourne une liste avec les 3 jeux de données générées
  }
}

# fonction alternative qui n'impose pas des conditions particulaires

MixData_generation <- function(p_dim, K_clust, mean_overlap, N_size, lam){
  
  # proportion des donnés au sein de chaque cluster 
  pi_min <- 1/K_clust- 0.05 
  
  # paramètre de chevauchement maximale 
  max_ovelp <- mean_overlap*(K_clust*(K_clust-1)/2) - (mean_overlap/2)
  
  # generation des mélanges 
  # chevauchement élevé 
  Mix_HO <- MixSim(BarOmega = mean_overlap*10,
                   MaxOmega = max_ovelp*10,
                   K = K_clust, 
                   p = p_dim,
                   PiLow = pi_min, 
                   resN=1000) 
  
  # chevauchement modèré 
  Mix_MO <- MixSim(BarOmega = mean_overlap,
                   MaxOmega = max_ovelp, 
                   K = K_clust,
                   p = p_dim, 
                   PiLow = pi_min) 
  
  # faible chevauchement  
  Mix_LO <- MixSim(BarOmega = mean_overlap/10,
                   MaxOmega = max_ovelp/10, 
                   K = K_clust,
                   p = p_dim, 
                   PiLow = pi_min) 
  
  # formation des jeux de données à partir du modèle de mélange générées avant 

  dataset_HO <- simdataset(n = N_size, Pi = Mix_HO$Pi, Mu = Mix_HO$Mu, S = Mix_HO$S, lambda = rep(lam, p_dim)) 
  dataset_MO <- simdataset(n = N_size, Pi = Mix_MO$Pi, Mu = Mix_MO$Mu, S = Mix_MO$S, lambda = rep(lam, p_dim))
  dataset_LO <- simdataset(n = N_size, Pi = Mix_LO$Pi, Mu = Mix_LO$Mu, S = Mix_LO$S, lambda = rep(lam, p_dim))
  
  # si la dimention du vecteur beta est deux, on peut effectuer une validation graphique
  if(p_dim==2){
    colors <- c("tomato", "lightpink", "lightskyblue", "red", "magenta", 
                "blue", "brown", "green", "Skyblue", "orange", "yellow","pink", "Violet")
    
    x11()
    par(mfrow=c(3,1), mar= c(1, 2, 2, 2) )
    plot(dataset_HO$X, col = colors[dataset_HO$id], pch=19, cex =0.8,
         xlab="", ylab="", main="High overlap",  xaxt='n',  yaxt='n')
    
    plot(dataset_MO$X, col = colors[dataset_MO$id], pch=19, cex =0.8,
         xlab="", ylab="", main="Moderated overlap",  xaxt='n', yaxt='n')
    
    plot(dataset_LO$X, col = colors[dataset_LO$id], pch=19, cex =0.8,
         xlab="", ylab="", main="Low overlap",  xaxt='n', yaxt='n')
  }
  
  # creation du data frame avec la donnée géneré 
  dataset_H <- data.frame(dataset_HO$X, dataset_HO$id)
  dataset_M <- data.frame(dataset_MO$X, dataset_MO$id)
  dataset_L <- data.frame(dataset_LO$X, dataset_LO$id)
  
  # creation des labels 
  x_p <- rep(0, p_dim)
  
  for(i in 1:p_dim){
    x_p[i] <- paste("x", i, sep = "_")
  }
  
  colnames(dataset_H) <- colnames(dataset_M) <- colnames(dataset_L) <- c(x_p, "label")
  
  mix_data <- vector("list", 3)
  mix_data[[1]] <- dataset_H
  mix_data[[2]] <- dataset_M
  mix_data[[3]] <- dataset_L
  
  return(mix_data)  # la fonction returne une liste avec les trois jeux de données
}

# Test de la fonction 
set.seed(123)
# no conditions particulaires 
P2K3_datasets <- MixData_generation(p_dim = 2, 
                                    K_clust = 3,
                                    mean_overlap = 0.005,
                                    N_size = 10000, 
                                    lam = 0.5)
# ajustements manuels
P2K3_datasets_cond <- MixData(p_dim = 2, 
                              K_clust = 3,
                              mean_overlap = 0.005,
                              N_size = 10000, 
                              lam = 0.5)

population_P2K3_HO <- P2K3_datasets[[1]] # base de donnés avec chevauchement elevé 
population_P2K3_MO <- P2K3_datasets[[2]] # base de donnés avec chevauchement moderé
population_P2K3_LO <- P2K3_datasets[[3]] # base de donnés avec faible chevauchement 

# suvgarder la donnée generée au besoin 

setwd("C:/Users/Documents/project/")

write.table(x=population_P2K3_HO, 
            file="P2K3_HO.txt", 
            sep = " ",
            dec = ".",
            row.names = FALSE, 
            col.names = TRUE)
