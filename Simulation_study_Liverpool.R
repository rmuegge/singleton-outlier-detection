library(sf)
library(spdep)
library(MASS) #mvrnorm
library(cluster)
library(dplyr)
library(tidyr)
library(INLA)

sf_use_s2(FALSE)

#set working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#the LSOA shapefiles are available at
#https://data.cambridgeshireinsight.org.uk/dataset/output-areas/resource/3bc4faef-38c7-417d-88a9-f302ad845ebe
Bdry <- st_read("Lower_Layer_Super_Output_Areas_December_2011_Generalised_Clipped__Boundaries_in_England_and_Wales.shp")
Bdry2 <- st_transform(x=Bdry, crs='+proj=longlat +datum=WGS84 +no_defs')

#the mapping from LSOA to LAD is available at
#https://geoportal.statistics.gov.uk/documents/d75d81ae855d470f8c2bb2425003457f/about
lads <- read_excel("LSOA11_WD21_LAD21_EW_LU.xlsx")
lads <- unique(lads[,c(1,2,5,6)])
colnames(lads) <- c("area","area_name","LAD_code","LAD_name")
lads <- lads[which(lads$LAD_name=="Liverpool"),]

#load RDS file created in "Application_asthma_England.R"
Asthma.Bdry <- readRDS("Asthma_Bdry_RDOS.rds")
prob_est <- median(Asthma.Bdry$proportion)
logit_prob_est <- log(prob_est/(1-prob_est))
Liv.asthma <- Asthma.Bdry[which(Asthma.Bdry$lsoa11cd %in% lads$area),c("lsoa11cd","population")]
Liv.asthma <- as.data.frame(Liv.asthma[,c("lsoa11cd","population")])
Liv.asthma <- as.data.frame(Liv.asthma[,c("lsoa11cd","population")])
rownames(Liv.asthma) <- c()
lads <- merge(x=lads,y=Liv.asthma,by.x="area",by.y="lsoa11cd")

#this is the geography used in the simulation study
Liv.Bdry <- merge(x=Bdry2,y=lads,by.x="lsoa11cd",by.y="area")

#create KNN neighbourhood matrix to generate data (W must be symmetric)
coords <- st_centroid(st_geometry(Liv.Bdry))
knn <- 10
nb.test <- knn2nb(knearneigh(coords, k=knn), row.names=Liv.Bdry$lsoa11cd)
W.nb <- make.sym.nb(nb.test)
summary(W.nb)
W <- nb2mat(W.nb, style = "B")
K <- nrow(W)
#compute W1 (number of neighbours for each area) 
W1 <- rowSums(W)

W.list <- nb2listw(W.nb) #for local Moran's I computation

#use adjacency matrix to generate outliers (ensures no two outliers share a border)
W.nbadj <- poly2nb(Liv.Bdry,row.names=Liv.Bdry$InterZone)
W_adj.list <- nb2listw(W.nbadj)
W_adj <- nb2mat(W.nbadj,style = "B")

#compute the Q matrix from the Leroux CAR model
rho <- 0.999 #use high correlation to ensure smoothness
Q <- rho*(diag(W1)-W)+(1-rho)*diag(rep(1,K))
Qinv <- solve(Q) #compute the inverse of Q

tau2 <- 0.1 #variance parameter

S <- list() #structure matrix used to compute RDOS
for(i in 1:K){
  S[[i]] <- which(W[i,]==1) 
}

#define the local neighbourhoods for outliers
Sadj <- list()
for(i in 1:K){
  Sadj[[i]] <- which(W_adj[i,]==1) 
}

#Gaussian kernel
Gk <- function(x,S,h){
  1/sqrt(2*pi)*(1/h)*exp(-1/2*((x-S)^2)/(h^2))
}

#modelling, kNN neighbours for k=6
coords <- st_centroid(st_geometry(Liv.Bdry))
knn_INLA <- 6 #number of neighbours to use in the INLA model
nb.test_INLA <- knn2nb(knearneigh(coords, k=knn_INLA), row.names=Liv.Bdry$lsoa11cd)
W.INLA_symm <- make.sym.nb(nb.test_INLA)
nb2INLA("map.6nn_all",W.INLA_symm)
W_all <- inla.read.graph(filename="map.6nn_all")

#create a function to run the simulation study
simulation_study <- function(m,v){
  Liv.Bdry$phi <- mvrnorm(n=1,mu=rep(0,K),Sigma=tau2*Qinv) #random effects
  Liv.Bdry$phi <- Liv.Bdry$phi-mean(Liv.Bdry$phi) #sum to zero constraint
  outliers <- c() #this will become a vector containing the outlier indices
  areas <- 1:K
  areasleft <- areas #areasleft contains all areas that do not border an outlier
  remove <- c() #areas that are outliers or sharing a border with an outlier
  for(i in 1:m){
    j <- sample(areasleft,size=1) #sample one area from the available ones
    outliers <- c(outliers,j) #make that area an outlier
    remove <- unique(c(remove,j,Sadj[[j]])) #the outlier and its bordering neighbours can't become outliers
    areasleft <- areas[-remove] #update remaining areas that could be outliers
  }
  
  u <- rep(0,K)
  for(i in 1:m){
    u[outliers[i]] <- sign(Liv.Bdry$phi[outliers[i]]-mean(Liv.Bdry$phi[Sadj[[outliers[i]]]]))*v[i] #shifted in direction opposite to mean of border sharing areas
  }
  Liv.Bdry$u <- u
  
  Liv.Bdry$prob <- exp(logit_prob_est+Liv.Bdry$phi+Liv.Bdry$u)/(1+exp(logit_prob_est+Liv.Bdry$phi+Liv.Bdry$u))
  
  Liv.Bdry$Y <- rbinom(n=K,size=Liv.Bdry$population,prob=Liv.Bdry$prob)
  
  Liv.Bdry$outl <- ifelse(1:K %in% outliers,1,0) #outlier ID
  
  Liv.Bdry$proportion <- Liv.Bdry$Y/Liv.Bdry$population
  
  absdiff <- c()
  for(i in 1:K){
    absdiff <- c(absdiff,median(abs(Liv.Bdry$proportion[S[[i]]]-Liv.Bdry$proportion[i])))
  }
  
  h1 <- median(absdiff)
  
  outlier_LSOAs <- Liv.Bdry[outliers,]
  inlier_LSOAs <- Liv.Bdry[-outliers,]
  
  h <- h1
  est_dens_context <- rep(0,K) #kernel density
  for(i in 1:K){
    est_dens_context[i] <- 1/(length(S[[i]])+1)*(sum(Gk(x=Liv.Bdry$proportion[i],S=Liv.Bdry$proportion[S[[i]]],h=h))+Gk(x=Liv.Bdry$proportion[i],S=Liv.Bdry$proportion[i],h=h))
  }
  
  RDOS_context <- rep(0,K) #RDOS values
  for(i in 1:K){
    RDOS_context[i] <- mean(est_dens_context[S[[i]]])/(est_dens_context[i])
  }
  Liv.Bdry$RDOS_context <- RDOS_context
  
  #increase the bandwidth incrementally until Kendall's rank coefficient is greater than 0.99
  tau <- 0
  c <- 1
  while(tau<0.99){
    c <- c+0.1
    h <- c*h1
    est_dens_context <- rep(0,K) #kernel density
    for(i in 1:K){
      est_dens_context[i] <- 1/(length(S[[i]])+1)*(sum(Gk(x=Liv.Bdry$proportion[i],S=Liv.Bdry$proportion[S[[i]]],h=h))+Gk(x=Liv.Bdry$proportion[i],S=Liv.Bdry$proportion[i],h=h))
    }
    
    RDOS_context <- rep(0,K) #RDOS values
    for(i in 1:K){
      RDOS_context[i] <- mean(est_dens_context[S[[i]]])/(est_dens_context[i])
    }
    Liv.Bdry$RDOS_context_new <- RDOS_context
    
    #data in decreasing order of the previous RDOS value
    rdos_data_context <- as.data.frame(Liv.Bdry[order(Liv.Bdry$RDOS_context,decreasing=TRUE),])
    rdos_data_context$ranking <- 1:K 
    
    #data in decreasing order of the current RDOS value 
    rdos_data_context_new <- as.data.frame(rdos_data_context[order(rdos_data_context$RDOS_context_new,decreasing=TRUE),])
    rdos_data_context_new$ranking_new <- 1:K
    
    #compute Kendall's rank coefficient tau
    kendall_test <- cor.test(rdos_data_context_new$ranking,rdos_data_context_new$ranking_new,method="kendall")
    tau <- kendall_test$estimate
    Liv.Bdry$RDOS_context <- Liv.Bdry$RDOS_context_new
  }
  
  rdos_data_context <- as.data.frame(Liv.Bdry[order(Liv.Bdry$RDOS_context,decreasing=TRUE),])
  rdos_data_context$ranking <- 1:K 
  rdos_data_context$inl <- 1-rdos_data_context$outl
  
  #compute the TPR and FPR after adding each observation to the outlier set by decreasing RDOS value
  TPR_rdos <- c(0,rep(NA,K))
  FPR_rdos <- c(0,rep(NA,K))
  for(i in 1:K){
    TPR_rdos[i+1] <- sum(rdos_data_context$outl[1:i])/m
    FPR_rdos[i+1] <- sum(rdos_data_context$inl[1:i])/(K-m)
  }
  
  #use local Moran's I for comparison
  lisaRslt <- localmoran(Liv.Bdry$Y,W.list)
  lRdat <- as.data.frame(lisaRslt)
  lRdat$lsoa11cd <- Liv.Bdry$lsoa11cd
  lRdat$outl <- Liv.Bdry$outl
  lRdat$inl <- 1-lRdat$outl
  lRdat <- lRdat[order(lRdat$Ii,decreasing=FALSE),] #by increasing local Moran's I 
  lRdat_abs <- lRdat[order(abs(lRdat$Ii),decreasing=TRUE),] #decreasing absolute value
  
  #TPR and FPR for local Moran's I in increasing value
  TPR_locM <- c(0,rep(NA,K))
  FPR_locM <- c(0,rep(NA,K))
  for(i in 1:K){
    TPR_locM[i+1] <- sum(lRdat$outl[1:i])/m
    FPR_locM[i+1] <- sum(lRdat$inl[1:i])/(K-m)
  }
  
  #TPR and FPR for local Moran's I in decreasing absolute value
  TPR_locM_abs <- c(0,rep(NA,K))
  FPR_locM_abs <- c(0,rep(NA,K))
  for(i in 1:K){
    TPR_locM_abs[i+1] <- sum(lRdat_abs$outl[1:i])/m
    FPR_locM_abs[i+1] <- sum(lRdat_abs$inl[1:i])/(K-m)
  }
  
  Liv.Bdry$idarea <- 1:K
  
  #fit the conventional BYM2 model
  res_no_outl <- inla(Y ~ f(idarea,model="bym2",graph=W_all,scale.model=TRUE),
                      family="binomial",
                      data=Liv.Bdry,
                      Ntrials=population,
                      control.compute=list(dic=TRUE,waic=TRUE,cpo=TRUE),
                      control.fixed = list(prec.intercept=0.001))
  
  phi_mean_no_outl <- res_no_outl$summary.hyperpar$mean[2]
  phi_CrI_no_outl <- res_no_outl$summary.hyperpar$`0.975quant`[2]-res_no_outl$summary.hyperpar$`0.025quant`[2]
  
  RMSE_no_outl <- sqrt(mean((Liv.Bdry$prob-res_no_outl$summary.fitted.values$mean)^2))
  RMSE_no_outl_inliers <- sqrt(mean((Liv.Bdry$prob[which(Liv.Bdry$outl==0)]-res_no_outl$summary.fitted.values$mean[which(Liv.Bdry$outl==0)])^2))
  RMSE_no_outl_outliers <- sqrt(mean((Liv.Bdry$prob[which(Liv.Bdry$outl==1)]-res_no_outl$summary.fitted.values$mean[which(Liv.Bdry$outl==1)])^2))
  Coverage_no_outl <- sum(Liv.Bdry$prob>res_no_outl$summary.fitted.values$`0.025quant` & Liv.Bdry$prob<res_no_outl$summary.fitted.values$`0.975quant`)/K
  Coverage_no_outl_inliers <- sum(Liv.Bdry$prob[which(Liv.Bdry$outl==0)]>res_no_outl$summary.fitted.values$`0.025quant`[which(Liv.Bdry$outl==0)] & Liv.Bdry$prob[which(Liv.Bdry$outl==0)]<res_no_outl$summary.fitted.values$`0.975quant`[which(Liv.Bdry$outl==0)])/(K-m)
  Coverage_no_outl_outliers <- sum(Liv.Bdry$prob[which(Liv.Bdry$outl==1)]>res_no_outl$summary.fitted.values$`0.025quant`[which(Liv.Bdry$outl==1)] & Liv.Bdry$prob[which(Liv.Bdry$outl==1)]<res_no_outl$summary.fitted.values$`0.975quant`[which(Liv.Bdry$outl==1)])/m
  MeanCrIwidth_no_outl <- mean(res_no_outl$summary.fitted.values$`0.975quant`-res_no_outl$summary.fitted.values$`0.025quant`) 
  MeanCrIwidth_no_outl_inliers <- mean(res_no_outl$summary.fitted.values$`0.975quant`[which(Liv.Bdry$outl==0)]-res_no_outl$summary.fitted.values$`0.025quant`[which(Liv.Bdry$outl==0)]) 
  MeanCrIwidth_no_outl_outliers <- mean(res_no_outl$summary.fitted.values$`0.975quant`[which(Liv.Bdry$outl==1)]-res_no_outl$summary.fitted.values$`0.025quant`[which(Liv.Bdry$outl==1)]) 
  MAE_no_outl <- median(abs(Liv.Bdry$prob-res_no_outl$summary.fitted.values$mean))
  MAE_no_outl_inliers <- median(abs(Liv.Bdry$prob[which(Liv.Bdry$outl==0)]-res_no_outl$summary.fitted.values$mean[which(Liv.Bdry$outl==0)]))
  MAE_no_outl_outliers <- median(abs(Liv.Bdry$prob[which(Liv.Bdry$outl==1)]-res_no_outl$summary.fitted.values$mean[which(Liv.Bdry$outl==1)]))
  WAIC_no_outl <- res_no_outl$waic$waic
  DIC_no_outl <- res_no_outl$dic$dic
  
  detect_pam <- function(oinit){
    #partitioning around medoids
    #two clusters
    clust <- pam(rdos_data_context$RDOS_context[1:oinit],2,metric="manhattan",nstart=15,cluster.only=TRUE)
    #if cluster 2 contains the outliers (large RDOS values)
    if(mean(rdos_data_context$RDOS_context[which(clust==1)])<mean(rdos_data_context$RDOS_context[which(clust==2)])){
      rdos_data_context$cluster <- c(clust,rep(1,K-oinit))
      rdos_data_context$cluster <- ifelse(rdos_data_context$cluster==1,"not a contextual outlier","contextual outlier")
    }else{#if cluster 1 contains the outliers (large RDOS values)
      rdos_data_context$cluster <- c(clust,rep(2,K-oinit))
      rdos_data_context$cluster <- ifelse(rdos_data_context$cluster==2,"not a contextual outlier","contextual outlier")
    }
    r_inl_max <- max(rdos_data_context$RDOS_context[which(rdos_data_context$cluster=="not a contextual outlier")])
    r_outl_min <- min(rdos_data_context$RDOS_context[which(rdos_data_context$cluster=="contextual outlier")])
    
    if(sum(Liv.Bdry$RDOS_context>r_inl_max)>0){
      top_n <- Liv.Bdry[which(Liv.Bdry$RDOS_context>r_inl_max),]
      n_context <- nrow(top_n) #the number of possible outliers identified using the RDOS values
      top_n <- top_n[order(top_n$RDOS_context,decreasing=TRUE),]
      not_context <- Liv.Bdry$lsoa11cd[-which(Liv.Bdry$lsoa11cd %in% top_n$lsoa11cd)]
      pprec <- sum(outlier_LSOAs$lsoa11cd %in% top_n$lsoa11cd)/(n_context)
      psens <- sum(outlier_LSOAs$lsoa11cd %in% top_n$lsoa11cd)/m
      pspec <- sum(inlier_LSOAs$lsoa11cd %in% not_context)/(K-m)
    }else{
      n_context <- 0
      pprec <- 0
      psens <- 0
      pspec <- 1
    }
    n_outl <- nrow(top_n)
    
    #compute the neighbourhood matrix for the inliers
    Liv.Bdry_inliers <- Liv.Bdry[which(Liv.Bdry$lsoa11cd %in% rdos_data_context$lsoa11cd[(n_outl+1):K]),]
    inlier_ids <- which(Liv.Bdry$lsoa11cd %in% rdos_data_context$lsoa11cd[(n_outl+1):K])
    outlier_ids <- which(Liv.Bdry$lsoa11cd %in% rdos_data_context$lsoa11cd[1:n_outl])
    coords <- st_centroid(st_geometry(Liv.Bdry_inliers))
    knn_INLA <- 6 #number of neighbours to use in the INLA model
    nb.test_INLA_inliers <- knn2nb(knearneigh(coords, k=knn_INLA), row.names=Liv.Bdry_inliers$lsoa11cd)
    W.INLA_symm_inliers <- make.sym.nb(nb.test_INLA_inliers)
    nb2INLA("map.6nn_inliers",W.INLA_symm_inliers)
    W_inliers <- inla.read.graph(filename="map.6nn_inliers")
    
    Liv.Bdry_inliers$idinliers <- 1:nrow(Liv.Bdry_inliers)
    
    Liv.Bdry$idinliers <- NA
    for(i in 1:K){
      Liv.Bdry$idinliers[i] <- ifelse(Liv.Bdry$lsoa11cd[i] %in% Liv.Bdry_inliers$lsoa11cd,Liv.Bdry_inliers$idinliers[which(Liv.Bdry_inliers$lsoa11cd==Liv.Bdry$lsoa11cd[i])],NA)
    }
    
    j <- 1
    Liv.Bdry$idoutliers <- NA
    for(i in 1:K){
      if(is.na(Liv.Bdry$idinliers[i])){
        Liv.Bdry$idoutliers[i] <- j
        j <- j+1
      }
    }
    
    #modified smoothing model BYM2-O
    res_outl_bym2 <- inla(Y ~ f(idinliers,model="bym2",graph=W_inliers,scale.model=TRUE) +
                            f(idoutliers,model="iid"),
                          family="binomial",
                          data=Liv.Bdry,
                          Ntrials=population,
                          control.compute=list(dic=TRUE,waic=TRUE,cpo=TRUE),
                          control.fixed = list(prec.intercept=0.001))
    RMSE_outl_bym2 <- sqrt(mean((Liv.Bdry$prob-res_outl_bym2$summary.fitted.values$mean)^2))
    RMSE_outl_bym2_inliers <- sqrt(mean((Liv.Bdry$prob[which(Liv.Bdry$outl==0)]-res_outl_bym2$summary.fitted.values$mean[which(Liv.Bdry$outl==0)])^2))
    RMSE_outl_bym2_outliers <- sqrt(mean((Liv.Bdry$prob[which(Liv.Bdry$outl==1)]-res_outl_bym2$summary.fitted.values$mean[which(Liv.Bdry$outl==1)])^2))
    RMSE_outl_bym2_identified_inliers <- sqrt(mean((Liv.Bdry$prob[inlier_ids]-res_outl_bym2$summary.fitted.values$mean[inlier_ids])^2))
    RMSE_outl_bym2_identified_outliers <- sqrt(mean((Liv.Bdry$prob[outlier_ids]-res_outl_bym2$summary.fitted.values$mean[outlier_ids])^2))
    Coverage_outl_bym2 <- sum(Liv.Bdry$prob>res_outl_bym2$summary.fitted.values$`0.025quant` & Liv.Bdry$prob<res_outl_bym2$summary.fitted.values$`0.975quant`)/K
    Coverage_outl_bym2_inliers <- sum(Liv.Bdry$prob[which(Liv.Bdry$outl==0)]>res_outl_bym2$summary.fitted.values$`0.025quant`[which(Liv.Bdry$outl==0)] & Liv.Bdry$prob[which(Liv.Bdry$outl==0)]<res_outl_bym2$summary.fitted.values$`0.975quant`[which(Liv.Bdry$outl==0)])/(K-m)
    Coverage_outl_bym2_outliers <- sum(Liv.Bdry$prob[which(Liv.Bdry$outl==1)]>res_outl_bym2$summary.fitted.values$`0.025quant`[which(Liv.Bdry$outl==1)] & Liv.Bdry$prob[which(Liv.Bdry$outl==1)]<res_outl_bym2$summary.fitted.values$`0.975quant`[which(Liv.Bdry$outl==1)])/m
    Coverage_outl_bym2_identified_inliers <- sum(Liv.Bdry$prob[inlier_ids]>res_outl_bym2$summary.fitted.values$`0.025quant`[inlier_ids] & Liv.Bdry$prob[inlier_ids]<res_outl_bym2$summary.fitted.values$`0.975quant`[inlier_ids])/(K-n_outl)
    Coverage_outl_bym2_identified_outliers <- sum(Liv.Bdry$prob[outlier_ids]>res_outl_bym2$summary.fitted.values$`0.025quant`[outlier_ids] & Liv.Bdry$prob[outlier_ids]<res_outl_bym2$summary.fitted.values$`0.975quant`[outlier_ids])/n_outl
    MeanCrIwidth_outl_bym2 <- mean(res_outl_bym2$summary.fitted.values$`0.975quant`-res_outl_bym2$summary.fitted.values$`0.025quant`) 
    MeanCrIwidth_outl_bym2_inliers <- mean(res_outl_bym2$summary.fitted.values$`0.975quant`[which(Liv.Bdry$outl==0)]-res_outl_bym2$summary.fitted.values$`0.025quant`[which(Liv.Bdry$outl==0)]) 
    MeanCrIwidth_outl_bym2_outliers <- mean(res_outl_bym2$summary.fitted.values$`0.975quant`[which(Liv.Bdry$outl==1)]-res_outl_bym2$summary.fitted.values$`0.025quant`[which(Liv.Bdry$outl==1)]) 
    MeanCrIwidth_outl_bym2_identified_inliers <- mean(res_outl_bym2$summary.fitted.values$`0.975quant`[inlier_ids]-res_outl_bym2$summary.fitted.values$`0.025quant`[inlier_ids]) 
    MeanCrIwidth_outl_bym2_identified_outliers <- mean(res_outl_bym2$summary.fitted.values$`0.975quant`[outlier_ids]-res_outl_bym2$summary.fitted.values$`0.025quant`[outlier_ids]) 
    MAE_outl_bym2 <- median(abs(Liv.Bdry$prob-res_outl_bym2$summary.fitted.values$mean))
    MAE_outl_bym2_inliers <- median(abs(Liv.Bdry$prob[which(Liv.Bdry$outl==0)]-res_outl_bym2$summary.fitted.values$mean[which(Liv.Bdry$outl==0)]))
    MAE_outl_bym2_outliers <- median(abs(Liv.Bdry$prob[which(Liv.Bdry$outl==1)]-res_outl_bym2$summary.fitted.values$mean[which(Liv.Bdry$outl==1)]))
    MAE_outl_bym2_identified_inliers <- median(abs(Liv.Bdry$prob[inlier_ids]-res_outl_bym2$summary.fitted.values$mean[inlier_ids]))
    MAE_outl_bym2_identified_outliers <- median(abs(Liv.Bdry$prob[outlier_ids]-res_outl_bym2$summary.fitted.values$mean[outlier_ids]))
    WAIC_outl_bym2 <- res_outl_bym2$waic$waic
    DIC_outl_bym2 <- res_outl_bym2$dic$dic
    phi_mean <- res_outl_bym2$summary.hyperpar$mean[2]
    phi_CrI <- res_outl_bym2$summary.hyperpar$`0.975quant`[2]-res_outl_bym2$summary.hyperpar$`0.025quant`[2]
    
    return(list(oinit,n_context,pprec,psens,pspec,phi_mean,phi_CrI,
                RMSE_outl_bym2,RMSE_outl_bym2_inliers,RMSE_outl_bym2_outliers,
                RMSE_outl_bym2_identified_inliers,RMSE_outl_bym2_identified_outliers,
                Coverage_outl_bym2,Coverage_outl_bym2_inliers,Coverage_outl_bym2_outliers,
                Coverage_outl_bym2_identified_inliers,Coverage_outl_bym2_identified_outliers,
                MeanCrIwidth_outl_bym2,MeanCrIwidth_outl_bym2_inliers,MeanCrIwidth_outl_bym2_outliers,
                MeanCrIwidth_outl_bym2_identified_inliers,MeanCrIwidth_outl_bym2_identified_outliers,
                MAE_outl_bym2,MAE_outl_bym2_inliers,MAE_outl_bym2_outliers,
                MAE_outl_bym2_identified_inliers,MAE_outl_bym2_identified_outliers,
                WAIC_outl_bym2,DIC_outl_bym2))
  }
  
  #apply PAM to 10% largest (30), 20% largest (60), and all RDOS values
  res_pam_30 <- detect_pam(oinit=30)
  res_pam_60 <- detect_pam(oinit=60)
  res_pam_all <- detect_pam(oinit=K)
  
  #best case scenario
  Liv.Bdry_inliers <- Liv.Bdry[which(Liv.Bdry$outl==0),]
  coords <- st_centroid(st_geometry(Liv.Bdry_inliers))
  knn_INLA <- 6 #number of neighbours to use in the INLA model
  nb.test_INLA_inliers <- knn2nb(knearneigh(coords, k=knn_INLA), row.names=Liv.Bdry_inliers$lsoa11cd)
  W.INLA_symm_inliers <- make.sym.nb(nb.test_INLA_inliers)
  nb2INLA("map.6nn_inliers",W.INLA_symm_inliers)
  W_inliers <- inla.read.graph(filename="map.6nn_inliers")
  
  Liv.Bdry_inliers$idinliers <- 1:nrow(Liv.Bdry_inliers)
  
  Liv.Bdry$idinliers <- NA
  for(i in 1:K){
    Liv.Bdry$idinliers[i] <- ifelse(Liv.Bdry$lsoa11cd[i] %in% Liv.Bdry_inliers$lsoa11cd,Liv.Bdry_inliers$idinliers[which(Liv.Bdry_inliers$lsoa11cd==Liv.Bdry$lsoa11cd[i])],NA)
  }
  
  j <- 1
  Liv.Bdry$idoutliers <- NA
  for(i in 1:K){
    if(is.na(Liv.Bdry$idinliers[i])){
      Liv.Bdry$idoutliers[i] <- j
      j <- j+1
    }
  }
  
  #fit the BYM2-O model under the best case scenario
  res_outl_bym2 <- inla(Y ~ f(idinliers,model="bym2",graph=W_inliers,scale.model=TRUE) +
                          f(idoutliers,model="iid"),
                        family="binomial",
                        data=Liv.Bdry,
                        Ntrials=population,
                        control.compute=list(dic=TRUE,waic=TRUE,cpo=TRUE),
                        control.fixed = list(prec.intercept=0.001))
  phi_mean_outl_bcs <- res_outl_bym2$summary.hyperpar$mean[2]
  phi_CrI_outl_bcs <- res_outl_bym2$summary.hyperpar$`0.975quant`[2]-res_no_outl$summary.hyperpar$`0.025quant`[2]
  RMSE_outl_bcs <- sqrt(mean((Liv.Bdry$prob-res_outl_bym2$summary.fitted.values$mean)^2))
  RMSE_outl_bcs_inliers <- sqrt(mean((Liv.Bdry$prob[which(Liv.Bdry$outl==0)]-res_outl_bym2$summary.fitted.values$mean[which(Liv.Bdry$outl==0)])^2))
  RMSE_outl_bcs_outliers <- sqrt(mean((Liv.Bdry$prob[which(Liv.Bdry$outl==1)]-res_outl_bym2$summary.fitted.values$mean[which(Liv.Bdry$outl==1)])^2))
  Coverage_outl_bcs <- sum(Liv.Bdry$prob>res_outl_bym2$summary.fitted.values$`0.025quant` & Liv.Bdry$prob<res_outl_bym2$summary.fitted.values$`0.975quant`)/K
  Coverage_outl_bcs_inliers <- sum(Liv.Bdry$prob[which(Liv.Bdry$outl==0)]>res_outl_bym2$summary.fitted.values$`0.025quant`[which(Liv.Bdry$outl==0)] & Liv.Bdry$prob[which(Liv.Bdry$outl==0)]<res_outl_bym2$summary.fitted.values$`0.975quant`[which(Liv.Bdry$outl==0)])/(K-m)
  Coverage_outl_bcs_outliers <- sum(Liv.Bdry$prob[which(Liv.Bdry$outl==1)]>res_outl_bym2$summary.fitted.values$`0.025quant`[which(Liv.Bdry$outl==1)] & Liv.Bdry$prob[which(Liv.Bdry$outl==1)]<res_outl_bym2$summary.fitted.values$`0.975quant`[which(Liv.Bdry$outl==1)])/m
  MeanCrIwidth_outl_bcs <- mean(res_outl_bym2$summary.fitted.values$`0.975quant`-res_outl_bym2$summary.fitted.values$`0.025quant`) 
  MeanCrIwidth_outl_bcs_inliers <- mean(res_outl_bym2$summary.fitted.values$`0.975quant`[which(Liv.Bdry$outl==0)]-res_outl_bym2$summary.fitted.values$`0.025quant`[which(Liv.Bdry$outl==0)]) 
  MeanCrIwidth_outl_bcs_outliers <- mean(res_outl_bym2$summary.fitted.values$`0.975quant`[which(Liv.Bdry$outl==1)]-res_outl_bym2$summary.fitted.values$`0.025quant`[which(Liv.Bdry$outl==1)]) 
  MAE_outl_bcs <- median(abs(Liv.Bdry$prob-res_outl_bym2$summary.fitted.values$mean))
  MAE_outl_bcs_inliers <- median(abs(Liv.Bdry$prob[which(Liv.Bdry$outl==0)]-res_outl_bym2$summary.fitted.values$mean[which(Liv.Bdry$outl==0)]))
  MAE_outl_bcs_outliers <- median(abs(Liv.Bdry$prob[which(Liv.Bdry$outl==1)]-res_outl_bym2$summary.fitted.values$mean[which(Liv.Bdry$outl==1)]))
  WAIC_outl_bcs <- res_outl_bym2$waic$waic
  DIC_outl_bcs <- res_outl_bym2$dic$dic
  
  return(list(c(c,h,TPR_rdos,FPR_rdos,TPR_locM,FPR_locM,TPR_locM_abs,FPR_locM_abs,
                phi_mean_no_outl,phi_CrI_no_outl,
                RMSE_no_outl,RMSE_no_outl_inliers,RMSE_no_outl_outliers,
                Coverage_no_outl,Coverage_no_outl_inliers,Coverage_no_outl_outliers,
                MeanCrIwidth_no_outl,MeanCrIwidth_no_outl_inliers,MeanCrIwidth_no_outl_outliers,
                MAE_no_outl,MAE_no_outl_inliers,MAE_no_outl_outliers,WAIC_no_outl,DIC_no_outl,res_pam_30,res_pam_60,res_pam_all,
                phi_mean_outl_bcs,phi_CrI_outl_bcs,RMSE_outl_bcs,RMSE_outl_bcs_inliers,RMSE_outl_bcs_outliers,
                Coverage_outl_bcs,Coverage_outl_bcs_inliers,Coverage_outl_bcs_outliers,
                MeanCrIwidth_outl_bcs,MeanCrIwidth_outl_bcs_inliers,MeanCrIwidth_outl_bcs_outliers,
                MAE_outl_bcs,MAE_outl_bcs_inliers,MAE_outl_bcs_outliers,WAIC_outl_bcs,DIC_outl_bcs)))
}

#Here, one example for the simulations with 15 large outliers
nsim <- 100
trials_rdos <- as.data.frame(matrix(data=NA,nrow=nsim,ncol=599))
colnames(trials_rdos) <- c("sim",paste0("TPR_n",0:K),paste0("FPR_n",0:K))
trials_lM <- as.data.frame(matrix(data=NA,nrow=nsim,ncol=599))
colnames(trials_lM) <- c("sim",paste0("TPR_locM_n",0:K),paste0("FPR_locM_n",0:K))
trials_lM_abs <- as.data.frame(matrix(data=NA,nrow=nsim,ncol=599))
colnames(trials_lM_abs) <- c("sim",paste0("TPR_locM_abs_n",0:K),paste0("FPR_locM_abs_n",0:K))
model_no_outl <- as.data.frame(matrix(data=NA,nrow=nsim,ncol=17))
colnames(model_no_outl) <- c("sim","phi_mean_no_outl","phi_CrI_no_outl","RMSE_no_outl","RMSE_no_outl_inliers","RMSE_no_outl_outliers",
                             "Coverage_no_outl","Coverage_no_outl_inliers","Coverage_no_outl_outliers",
                             "MeanCrIwidth_no_outl","MeanCrIwidth_no_outl_inliers","MeanCrIwidth_no_outl_outliers",
                             "MAE_no_outl","MAE_no_outl_inliers","MAE_no_outl_outliers","WAIC_no_outl","DIC_no_outl")
model_rdos_30 <- as.data.frame(matrix(data=NA,nrow=nsim,ncol=30))
colnames(model_rdos_30) <- c("sim","oinit","n_outl","pprec","psens","pspec","phi_mean","phi_CrI",
                             "RMSE_outl_bym2","RMSE_outl_bym2_inliers","RMSE_outl_bym2_outliers",
                             "RMSE_outl_bym2_identified_inliers","RMSE_outl_bym2_identified_outliers",
                             "Coverage_outl_bym2","Coverage_outl_bym2_inliers","Coverage_outl_bym2_outliers",
                             "Coverage_outl_bym2_identified_inliers","Coverage_outl_bym2_identified_outliers",
                             "MeanCrIwidth_outl_bym2","MeanCrIwidth_outl_bym2_inliers","MeanCrIwidth_outl_bym2_outliers",
                             "MeanCrIwidth_outl_bym2_identified_inliers","MeanCrIwidth_outl_bym2_identified_outliers",
                             "MAE_outl_bym2","MAE_outl_bym2_inliers","MAE_outl_bym2_outliers",
                             "MAE_outl_bym2_identified_inliers","MAE_outl_bym2_identified_outliers",
                             "WAIC_outl_bym2","DIC_outl_bym2")
model_rdos_60 <- as.data.frame(matrix(data=NA,nrow=nsim,ncol=30))
colnames(model_rdos_60) <- c("sim","oinit","n_outl","pprec","psens","pspec","phi_mean","phi_CrI",
                             "RMSE_outl_bym2","RMSE_outl_bym2_inliers","RMSE_outl_bym2_outliers",
                             "RMSE_outl_bym2_identified_inliers","RMSE_outl_bym2_identified_outliers",
                             "Coverage_outl_bym2","Coverage_outl_bym2_inliers","Coverage_outl_bym2_outliers",
                             "Coverage_outl_bym2_identified_inliers","Coverage_outl_bym2_identified_outliers",
                             "MeanCrIwidth_outl_bym2","MeanCrIwidth_outl_bym2_inliers","MeanCrIwidth_outl_bym2_outliers",
                             "MeanCrIwidth_outl_bym2_identified_inliers","MeanCrIwidth_outl_bym2_identified_outliers",
                             "MAE_outl_bym2","MAE_outl_bym2_inliers","MAE_outl_bym2_outliers",
                             "MAE_outl_bym2_identified_inliers","MAE_outl_bym2_identified_outliers",
                             "WAIC_outl_bym2","DIC_outl_bym2")
model_rdos_all <- as.data.frame(matrix(data=NA,nrow=nsim,ncol=30))
colnames(model_rdos_all) <- c("sim","oinit","n_outl","pprec","psens","pspec","phi_mean","phi_CrI",
                              "RMSE_outl_bym2","RMSE_outl_bym2_inliers","RMSE_outl_bym2_outliers",
                              "RMSE_outl_bym2_identified_inliers","RMSE_outl_bym2_identified_outliers",
                              "Coverage_outl_bym2","Coverage_outl_bym2_inliers","Coverage_outl_bym2_outliers",
                              "Coverage_outl_bym2_identified_inliers","Coverage_outl_bym2_identified_outliers",
                              "MeanCrIwidth_outl_bym2","MeanCrIwidth_outl_bym2_inliers","MeanCrIwidth_outl_bym2_outliers",
                              "MeanCrIwidth_outl_bym2_identified_inliers","MeanCrIwidth_outl_bym2_identified_outliers",
                              "MAE_outl_bym2","MAE_outl_bym2_inliers","MAE_outl_bym2_outliers",
                              "MAE_outl_bym2_identified_inliers","MAE_outl_bym2_identified_outliers",
                              "WAIC_outl_bym2","DIC_outl_bym2")
model_outl_bcs <- as.data.frame(matrix(data=NA,nrow=nsim,ncol=17))
colnames(model_outl_bcs) <- c("sim","phi_mean_outl_bcs","phi_CrI_outl_bcs",
                              "RMSE_outl_bcs","RMSE_outl_bcs_inliers","RMSE_outl_bcs_outliers",
                              "Coverage_outl_bcs","Coverage_outl_bcs_inliers","Coverage_outl_bcs_outliers",
                              "MeanCrIwidth_outl_bcs","MeanCrIwidth_outl_bcs_inliers","MeanCrIwidth_outl_bcs_outliers",
                              "MAE_outl_bcs","MAE_outl_bcs_inliers","MAE_outl_bcs_outliers","WAIC_outl_bcs","DIC_outl_bcs")
h_vec <- rep(NA,nsim)
c_vec <- rep(NA,nsim)

for(i in 1:nsim){
  res <- unlist(simulation_study(m=15,v=c(rep(0.7,15))))
  c_vec[i] <- res[1]
  h_vec[i] <- res[2]
  trials_rdos[i,] <- c(i,res[3:600])
  trials_lM[i,] <- c(i,res[601:1198])
  trials_lM_abs[i,] <- c(i,res[1199:1796])
  model_no_outl[i,] <- c(i,res[1797:1812])
  model_rdos_30[i,] <- c(i,res[1813:1841])
  model_rdos_60[i,] <- c(i,res[1842:1870])
  model_rdos_all[i,] <- c(i,res[1871:1899])
  model_outl_bcs[i,] <- c(i,res[1900:1915]) 
}

saveRDS(trials_rdos,"ROC_rdos_m15_v7_k10symm.rds")
saveRDS(trials_lM,"ROC_lM_m15_v7_k10symm.rds")
saveRDS(trials_lM_abs,"ROC_lM_abs_m15_v7_k10symm.rds")
saveRDS(model_no_outl,"model_no_outl_m15_v7.rds")
saveRDS(model_rdos_30,"model_rdos_m15_v7_30.rds")
saveRDS(model_rdos_60,"model_rdos_m15_v7_60.rds")
saveRDS(model_rdos_all,"model_rdos_m15_v7_all.rds")
saveRDS(model_outl_bcs,"model_outl_bcs_m15_v7.rds")
saveRDS(h_vec,"h_vec_m15_v7.rds")
saveRDS(c_vec,"c_vec_m15_v7.rds")

#local Moran's I, increasing
trials_lM <- readRDS("ROC_lM_m15_v7_k10symm.rds")
lM_ROC_TPR <- trials_lM[,1:300]
lM_ROC_TPR_long <- pivot_longer(lM_ROC_TPR,cols=2:300)
lM_ROC_TPR_long$name <- gsub(".*_","",lM_ROC_TPR_long$name)
colnames(lM_ROC_TPR_long) <- c("sim","rank","TPR")

lM_ROC_FPR <- trials_lM[,c(1,301:599)]
lM_ROC_FPR_long <- pivot_longer(lM_ROC_FPR,cols=2:300)
lM_ROC_FPR_long$name <- gsub(".*_","",lM_ROC_FPR_long$name)
colnames(lM_ROC_FPR_long) <- c("sim","rank","FPR")

lM_ROC <- cbind(lM_ROC_TPR_long,lM_ROC_FPR_long$FPR)
colnames(lM_ROC) <- c(colnames(lM_ROC_TPR_long),"FPR")

compute_AUC <- function(ROC_dat){
  max_FPRs <- ROC_dat %>% group_by(sim,FPR) %>% summarise(TPR=max(TPR))
  pivots <- max_FPRs %>% group_by(sim,TPR) %>% summarise(FPR=max(FPR))
  pivots$prev <- 0
  for(i in 2:nrow(pivots)){
    if(pivots$sim[i]==pivots$sim[i-1]){
      pivots$prev[i] <- pivots$FPR[i-1]
    }else{
      pivots$prev[i] <- 0
    }
  }
  pivots$diff <- pivots$FPR-pivots$prev
  AUCs <- pivots %>% group_by(sim) %>% summarise(AUC=sum(TPR*diff))
  return(AUCs)
}
AUCs <- compute_AUC(lM_ROC)
AUCs <- AUCs[order(AUCs$AUC),]
sim_median <- AUCs$sim[50] 
sim_lower <- AUCs$sim[3]
sim_upper <- AUCs$sim[98]
AUCs$AUC[50]
AUCs$AUC[3]
AUCs$AUC[98]

#Figure 2 (a), middle row 
ggplot() + geom_line(data=lM_ROC[which(lM_ROC$sim==sim_median),],aes(FPR,TPR)) +
  geom_line(data=lM_ROC[which(lM_ROC$sim==sim_lower),],aes(FPR,TPR),linetype="dashed") +
  geom_line(data=lM_ROC[which(lM_ROC$sim==sim_upper),],aes(FPR,TPR),linetype="dashed") +
  ggtitle("Local Moran's I") + xlab("FPR") + ylab("TPR") +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=13),
        title=element_text(size=13)) +
  annotate("label",x=0.55,y=0.1,label="AUC=0.55 (0.35,0.73)",size=4.5)

#local Moran's I, decreasing absolute value
trials_lM_abs <- readRDS("ROC_lM_abs_m15_v7_k10symm.rds")
lM_ROC_TPR_abs <- trials_lM_abs[,1:300]
lM_ROC_TPR_abs_long <- pivot_longer(lM_ROC_TPR_abs,cols=2:300)
lM_ROC_TPR_abs_long$name <- gsub(".*_","",lM_ROC_TPR_abs_long$name)
colnames(lM_ROC_TPR_abs_long) <- c("sim","rank","TPR")

lM_ROC_FPR_abs <- trials_lM_abs[,c(1,301:599)]
lM_ROC_FPR_abs_long <- pivot_longer(lM_ROC_FPR_abs,cols=2:300)
lM_ROC_FPR_abs_long$name <- gsub(".*_","",lM_ROC_FPR_abs_long$name)
colnames(lM_ROC_FPR_abs_long) <- c("sim","rank","FPR")

lM_ROC_abs <- cbind(lM_ROC_TPR_abs_long,lM_ROC_FPR_abs_long$FPR)
colnames(lM_ROC_abs) <- c(colnames(lM_ROC_TPR_abs_long),"FPR")

AUCs_abs <- compute_AUC(lM_ROC_abs)
AUCs_abs <- AUCs_abs[order(AUCs_abs$AUC),]
sim_median_abs <- AUCs_abs$sim[50] 
sim_lower_abs <- AUCs_abs$sim[3]
sim_upper_abs <- AUCs_abs$sim[98]
AUCs_abs$AUC[50]
AUCs_abs$AUC[3]
AUCs_abs$AUC[98]

#Figure 2 (a), bottom row 
ggplot() + geom_line(data=lM_ROC_abs[which(lM_ROC_abs$sim==sim_median_abs),],aes(FPR,TPR)) +
  geom_line(data=lM_ROC_abs[which(lM_ROC_abs$sim==sim_lower_abs),],aes(FPR,TPR),linetype="dashed") +
  geom_line(data=lM_ROC_abs[which(lM_ROC_abs$sim==sim_upper_abs),],aes(FPR,TPR),linetype="dashed") +
  ggtitle("Abs. local Moran's I") + xlab("FPR") + ylab("TPR") + 
  theme(axis.text=element_text(size=12),axis.title=element_text(size=13),
        title=element_text(size=13)) +
  annotate("label",x=0.55,y=0.1,label="AUC=0.80 (0.70,0.88)",size=4.5)

#RDOS
trials_rdos <- readRDS("ROC_rdos_m15_v7_k10symm.rds")
RDOS_ROC <- trials_rdos
RDOS_ROC_TPR <- RDOS_ROC[,c(1,2:300)]
RDOS_ROC_TPR_long <- pivot_longer(RDOS_ROC_TPR,cols=2:300)
RDOS_ROC_TPR_long$name <- gsub(".*_","",RDOS_ROC_TPR_long$name)
colnames(RDOS_ROC_TPR_long) <- c("sim","rank","TPR")

RDOS_ROC_FPR <- RDOS_ROC[,c(1,301:599)]
RDOS_ROC_FPR_long <- pivot_longer(RDOS_ROC_FPR,cols=2:300)
RDOS_ROC_FPR_long$name <- gsub(".*_","",RDOS_ROC_FPR_long$name)
colnames(RDOS_ROC_FPR_long) <- c("sim","rank","FPR")

RDOS_ROC_long <- cbind(RDOS_ROC_TPR_long,RDOS_ROC_FPR_long$FPR)
colnames(RDOS_ROC_long) <- c(colnames(RDOS_ROC_TPR_long),"FPR")
AUCs <- compute_AUC(RDOS_ROC_long)
AUCs <- AUCs[order(AUCs$AUC),]
sim_median <- AUCs$sim[50] 
sim_lower <- AUCs$sim[3]
sim_upper <- AUCs$sim[98]
AUCs$AUC[50]
AUCs$AUC[3]
AUCs$AUC[98]

#Figure 2 (a), top row 
ggplot() + geom_line(data=RDOS_ROC_long[which(RDOS_ROC_long$sim==sim_median),],aes(FPR,TPR)) +
  geom_line(data=RDOS_ROC_long[which(RDOS_ROC_long$sim==sim_lower),],aes(FPR,TPR),linetype="dashed") +
  geom_line(data=RDOS_ROC_long[which(RDOS_ROC_long$sim==sim_upper),],aes(FPR,TPR),linetype="dashed") +
  ggtitle("RDOS") +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=13),
        title=element_text(size=13)) + xlab("FPR") + ylab("TPR") +
  annotate("label",x=0.55,y=0.1,label="AUC=1.00 (1.00,1.00)",size=4.5)

#bandwidth summaries presented in table S1
c_vec_m15_v7 <- readRDS("c_vec_m15_v7.rds")
paste(format(round(mean(c_vec_m15_v7),2),nsmall=2),paste0("(",format(round(quantile(c_vec_m15_v7,probs=c(0.025)),2),nsmall=2),",",format(round(quantile(c_vec_m15_v7,probs=0.975),2),nsmall=2),")"),sep=" ")
c_vec_m15_v3 <- readRDS("c_vec_m15_v3.rds")
paste(format(round(mean(c_vec_m15_v3),2),nsmall=2),paste0("(",format(round(quantile(c_vec_m15_v3,probs=c(0.025)),2),nsmall=2),",",format(round(quantile(c_vec_m15_v3,probs=0.975),2),nsmall=2),")"),sep=" ")
c_vec_m3_v7 <- readRDS("c_vec_m3_v7.rds")
paste(format(round(mean(c_vec_m3_v7),2),nsmall=2),paste0("(",format(round(quantile(c_vec_m3_v7,probs=c(0.025)),2),nsmall=2),",",format(round(quantile(c_vec_m3_v7,probs=0.975),2),nsmall=2),")"),sep=" ")
c_vec_m3_v3 <- readRDS("c_vec_m3_v3.rds")
paste(format(round(mean(c_vec_m3_v3),2),nsmall=2),paste0("(",format(round(quantile(c_vec_m3_v3,probs=c(0.025)),2),nsmall=2),",",format(round(quantile(c_vec_m3_v3,probs=0.975),2),nsmall=2),")"),sep=" ")

#PAM for RDOS and model comparisons
#As an example, for 15 large outliers
sims <- readRDS("model_rdos_m15_v7_30.rds") #10% largest RDOS values
#sims <- readRDS("model_rdos_m15_v7_60.rds") #20% largest RDOS values
#sims <- readRDS("model_rdos_m15_v7_all.rds") #all RDOS values

summ <- as.data.frame(rbind(paste(format(round(mean(sims$n_outl),2),nsmall=2),paste0("(",paste(round(quantile(sims$n_outl,probs=0.025)),round(quantile(sims$n_outl,probs=0.975)),sep=","),")"),sep=" "),
                            paste(format(round(mean(sims$pprec),2),nsmall=2),paste0("(",paste(format(round(quantile(sims$pprec,probs=0.025),2),nsmall=2),format(round(quantile(sims$pprec,probs=0.975),2),nsmall=2),sep=","),")"),sep=" "),
                            paste(format(round(mean(sims$psens),2),nsmall=2),paste0("(",paste(format(round(quantile(sims$psens,probs=0.025),2),nsmall=2),format(round(quantile(sims$psens,probs=0.975),2),nsmall=2),sep=","),")"),sep=" "),
                            paste(format(round(mean(sims$pspec),2),nsmall=2),paste0("(",paste(format(round(quantile(sims$pspec,probs=0.025),2),nsmall=2),format(round(quantile(sims$pspec,probs=0.975),2),nsmall=2),sep=","),")"),sep=" ")))

#used for building up Table S2
print(xtable(summ), include.rownames=FALSE)

#model performance, RMSE
model_rdos_m15_v7_30 <- readRDS("model_rdos_m15_v7_30.rds")
res_rdos_m15_v7_30 <- model_rdos_m15_v7_30[,c("oinit","phi_mean","phi_CrI",
                                              "RMSE_outl_bym2","RMSE_outl_bym2_inliers","RMSE_outl_bym2_outliers",
                                              "RMSE_outl_bym2_identified_inliers","RMSE_outl_bym2_identified_outliers",
                                              "Coverage_outl_bym2","Coverage_outl_bym2_inliers","Coverage_outl_bym2_outliers",
                                              "Coverage_outl_bym2_identified_inliers","Coverage_outl_bym2_identified_outliers",
                                              "MeanCrIwidth_outl_bym2","MeanCrIwidth_outl_bym2_inliers","MeanCrIwidth_outl_bym2_outliers",
                                              "MeanCrIwidth_outl_bym2_identified_inliers","MeanCrIwidth_outl_bym2_identified_outliers",
                                              "MAE_outl_bym2","MAE_outl_bym2_inliers","MAE_outl_bym2_outliers",
                                              "MAE_outl_bym2_identified_inliers","MAE_outl_bym2_identified_outliers",
                                              "WAIC_outl_bym2","DIC_outl_bym2")]
colnames(res_rdos_m15_v7_30) <- c("method","phi_mean","phi_CrI","RMSE","RMSE_inliers","RMSE_outliers",
                                  "RMSE_identified_inliers","RMSE_identified_outliers",
                                  "Coverage","Coverage_inliers","Coverage_outliers",
                                  "Coverage_identified_inliers","Coverage_identified_outliers",
                                  "MeanCrIwidth","MeanCrIwidth_inliers","MeanCrIwidth_outliers",
                                  "MeanCrIwidth_identified_inliers","MeanCrIwidth_identified_outliers",
                                  "MAE","MAE_inliers","MAE_outliers","MAE_identified_inliers","MAE_identified_outliers",
                                  "WAIC","DIC")
model_rdos_m15_v7_60 <- readRDS("model_rdos_m15_v7_60.rds")
res_rdos_m15_v7_60 <- model_rdos_m15_v7_60[,c("oinit","phi_mean","phi_CrI",
                                              "RMSE_outl_bym2","RMSE_outl_bym2_inliers","RMSE_outl_bym2_outliers",
                                              "RMSE_outl_bym2_identified_inliers","RMSE_outl_bym2_identified_outliers",
                                              "Coverage_outl_bym2","Coverage_outl_bym2_inliers","Coverage_outl_bym2_outliers",
                                              "Coverage_outl_bym2_identified_inliers","Coverage_outl_bym2_identified_outliers",
                                              "MeanCrIwidth_outl_bym2","MeanCrIwidth_outl_bym2_inliers","MeanCrIwidth_outl_bym2_outliers",
                                              "MeanCrIwidth_outl_bym2_identified_inliers","MeanCrIwidth_outl_bym2_identified_outliers",
                                              "MAE_outl_bym2","MAE_outl_bym2_inliers","MAE_outl_bym2_outliers",
                                              "MAE_outl_bym2_identified_inliers","MAE_outl_bym2_identified_outliers",
                                              "WAIC_outl_bym2","DIC_outl_bym2")]
colnames(res_rdos_m15_v7_60) <- c("method","phi_mean","phi_CrI","RMSE","RMSE_inliers","RMSE_outliers",
                                  "RMSE_identified_inliers","RMSE_identified_outliers",
                                  "Coverage","Coverage_inliers","Coverage_outliers",
                                  "Coverage_identified_inliers","Coverage_identified_outliers",
                                  "MeanCrIwidth","MeanCrIwidth_inliers","MeanCrIwidth_outliers",
                                  "MeanCrIwidth_identified_inliers","MeanCrIwidth_identified_outliers",
                                  "MAE","MAE_inliers","MAE_outliers","MAE_identified_inliers","MAE_identified_outliers",
                                  "WAIC","DIC")
model_rdos_m15_v7_all <- readRDS("model_rdos_m15_v7_all.rds")
res_rdos_m15_v7_all <- model_rdos_m15_v7_all[,c("oinit","phi_mean","phi_CrI",
                                                "RMSE_outl_bym2","RMSE_outl_bym2_inliers","RMSE_outl_bym2_outliers",
                                                "RMSE_outl_bym2_identified_inliers","RMSE_outl_bym2_identified_outliers",
                                                "Coverage_outl_bym2","Coverage_outl_bym2_inliers","Coverage_outl_bym2_outliers",
                                                "Coverage_outl_bym2_identified_inliers","Coverage_outl_bym2_identified_outliers",
                                                "MeanCrIwidth_outl_bym2","MeanCrIwidth_outl_bym2_inliers","MeanCrIwidth_outl_bym2_outliers",
                                                "MeanCrIwidth_outl_bym2_identified_inliers","MeanCrIwidth_outl_bym2_identified_outliers",
                                                "MAE_outl_bym2","MAE_outl_bym2_inliers","MAE_outl_bym2_outliers",
                                                "MAE_outl_bym2_identified_inliers","MAE_outl_bym2_identified_outliers",
                                                "WAIC_outl_bym2","DIC_outl_bym2")]
colnames(res_rdos_m15_v7_all) <- c("method","phi_mean","phi_CrI","RMSE","RMSE_inliers","RMSE_outliers",
                                   "RMSE_identified_inliers","RMSE_identified_outliers",
                                   "Coverage","Coverage_inliers","Coverage_outliers",
                                   "Coverage_identified_inliers","Coverage_identified_outliers",
                                   "MeanCrIwidth","MeanCrIwidth_inliers","MeanCrIwidth_outliers",
                                   "MeanCrIwidth_identified_inliers","MeanCrIwidth_identified_outliers",
                                   "MAE","MAE_inliers","MAE_outliers","MAE_identified_inliers","MAE_identified_outliers",
                                   "WAIC","DIC")

res_rdos_m15_v7 <- as.data.frame(rbind(res_rdos_m15_v7_30,res_rdos_m15_v7_60,res_rdos_m15_v7_all))
res_rdos_m15_v7$method <- factor(res_rdos_m15_v7$method,levels=c(30,60,298),labels=c("BYM2-O, PAM(10%)","BYM2-O, PAM(20%)","BYM2-O, PAM(all)"))
RMSE_rdos_m15_v7 <- res_rdos_m15_v7[,c("method","RMSE")]

model_no_outl_m15_v7 <- readRDS("model_no_outl_m15_v7.rds")
res_no_outl_m15_v7 <- as.data.frame(model_no_outl_m15_v7[,c("phi_mean_no_outl","phi_CrI_no_outl",
                                                            "RMSE_no_outl","RMSE_no_outl_inliers","RMSE_no_outl_outliers",
                                                            "Coverage_no_outl","Coverage_no_outl_inliers","Coverage_no_outl_outliers",
                                                            "MeanCrIwidth_no_outl","MeanCrIwidth_no_outl_inliers","MeanCrIwidth_no_outl_outliers",
                                                            "MAE_no_outl","MAE_no_outl_inliers","MAE_no_outl_outliers","WAIC_no_outl","DIC_no_outl")])
res_no_outl_m15_v7$method <- "BYM2"
colnames(res_no_outl_m15_v7) <- c("phi_mean","phi_CrI","RMSE","RMSE_inliers","RMSE_outliers",
                                  "Coverage","Coverage_inliers","Coverage_outliers",
                                  "MeanCrIwidth","MeanCrIwidth_inliers","MeanCrIwidth_outliers",
                                  "MAE","MAE_inliers","MAE_outliers","WAIC","DIC","method")

RMSE_no_outl_m15_v7 <- res_no_outl_m15_v7[,c("method","RMSE")]

model_outl_bcs_m15_v7 <- readRDS("model_outl_bcs_m15_v7.rds")
summary(model_outl_bcs_m15_v7$phi_mean_outl_bcs)
res_outl_bcs_m15_v7 <- as.data.frame(model_outl_bcs_m15_v7[,c("phi_mean_outl_bcs","phi_CrI_outl_bcs",
                                                              "RMSE_outl_bcs","RMSE_outl_bcs_inliers","RMSE_outl_bcs_outliers",
                                                              "Coverage_outl_bcs","Coverage_outl_bcs_inliers","Coverage_outl_bcs_outliers",
                                                              "MeanCrIwidth_outl_bcs","MeanCrIwidth_outl_bcs_inliers","MeanCrIwidth_outl_bcs_outliers",
                                                              "MAE_outl_bcs","MAE_outl_bcs_inliers","MAE_outl_bcs_outliers","WAIC_outl_bcs","DIC_outl_bcs")])
res_outl_bcs_m15_v7$method <- "BYM2-O, best case"
colnames(res_outl_bcs_m15_v7) <- c("phi_mean","phi_CrI","RMSE","RMSE_inliers","RMSE_outliers",
                                   "Coverage","Coverage_inliers","Coverage_outliers",
                                   "MeanCrIwidth","MeanCrIwidth_inliers","MeanCrIwidth_outliers",
                                   "MAE","MAE_inliers","MAE_outliers","WAIC","DIC","method")
RMSE_res_outl_bcs_m15_v7 <- res_outl_bcs_m15_v7[,c("method","RMSE")]

RMSE_comp_m15_v7 <- as.data.frame(rbind(RMSE_rdos_m15_v7,RMSE_no_outl_m15_v7,RMSE_res_outl_bcs_m15_v7))
RMSE_comp_m15_v7$method <- factor(RMSE_comp_m15_v7$method,levels=c("BYM2","BYM2-O, PAM(10%)","BYM2-O, PAM(20%)","BYM2-O, PAM(all)","BYM2-O, best case"),
                                  labels=c("BYM2","BYM2-O, PAM(10%)","BYM2-O, PAM(20%)","BYM2-O, PAM(all)","BYM2-O, best case"))

Coverage_rdos_m15_v7 <- res_rdos_m15_v7[,c("method","Coverage")]
Coverage_no_outl_m15_v7 <- res_no_outl_m15_v7[,c("method","Coverage")]
Coverage_res_outl_bcs_m15_v7 <- res_outl_bcs_m15_v7[,c("method","Coverage")]
Coverage_comp_m15_v7 <- as.data.frame(rbind(Coverage_rdos_m15_v7,Coverage_no_outl_m15_v7,Coverage_res_outl_bcs_m15_v7))
Coverage_comp_m15_v7$method <- factor(Coverage_comp_m15_v7$method,levels=c("BYM2","BYM2-O, PAM(10%)","BYM2-O, PAM(20%)","BYM2-O, PAM(all)","BYM2-O, best case"),
                                      labels=c("BYM2","BYM2-O, PAM(10%)","BYM2-O, PAM(20%)","BYM2-O, PAM(all)","BYM2-O, best case"))

MeanCrIwidth_rdos_m15_v7 <- res_rdos_m15_v7[,c("method","MeanCrIwidth")]
MeanCrIwidth_no_outl_m15_v7 <- res_no_outl_m15_v7[,c("method","MeanCrIwidth")]
MeanCrIwidth_res_outl_bcs_m15_v7 <- res_outl_bcs_m15_v7[,c("method","MeanCrIwidth")]

MeanCrIwidth_comp_m15_v7 <- as.data.frame(rbind(MeanCrIwidth_rdos_m15_v7,MeanCrIwidth_no_outl_m15_v7,MeanCrIwidth_res_outl_bcs_m15_v7))
MeanCrIwidth_comp_m15_v7$method <- factor(MeanCrIwidth_comp_m15_v7$method,levels=c("BYM2","BYM2-O, PAM(10%)","BYM2-O, PAM(20%)","BYM2-O, PAM(all)","BYM2-O, best case"),
                                          labels=c("BYM2","BYM2-O, PAM(10%)","BYM2-O, PAM(20%)","BYM2-O, PAM(all)","BYM2-O, best case"))

MAE_rdos_m15_v7 <- res_rdos_m15_v7[,c("method","MAE")]
MAE_no_outl_m15_v7 <- res_no_outl_m15_v7[,c("method","MAE")]
MAE_res_outl_bcs_m15_v7 <- res_outl_bcs_m15_v7[,c("method","MAE")]
MAE_comp_m15_v7 <- as.data.frame(rbind(MAE_rdos_m15_v7,MAE_no_outl_m15_v7,MAE_res_outl_bcs_m15_v7))
MAE_comp_m15_v7$method <- factor(MAE_comp_m15_v7$method,levels=c("BYM2","BYM2-O, PAM(10%)","BYM2-O, PAM(20%)","BYM2-O, PAM(all)","BYM2-O, best case"),
                                 labels=c("BYM2","BYM2-O, PAM(10%)","BYM2-O, PAM(20%)","BYM2-O, PAM(all)","BYM2-O, best case"))

#repeat the above for m15_v3, m3_v7, and m3_v3
#combine the results, as shown below
RMSE_comp_m15_v7$group <- "m15v7"
RMSE_comp_m15_v3$group <- "m15v3"
RMSE_comp_m3_v7$group <- "m3v7"
RMSE_comp_m3_v3$group <- "m3v3"
RMSE_comp <- rbind(RMSE_comp_m15_v7,RMSE_comp_m15_v3,RMSE_comp_m3_v7,RMSE_comp_m3_v3)
RMSE_comp$group <- factor(RMSE_comp$group,levels=c("m15v7","m15v3","m3v7","m3v3"),
                          labels=c("15 large outliers","15 small outliers","3 large outliers","3 small outliers"))
#Figure 3 (a)
ggplot(data=RMSE_comp, aes(y=RMSE, fill=factor(method))) +
  geom_boxplot() + scale_x_continuous(breaks=NULL) +  
  facet_grid(~group) +  # Shrink the x-axis space
  scale_fill_manual(values=c("#CC79A7","#238b45","#74c476","#bae4b3","#1E88E5")) +
  theme(axis.title=element_text(size=14), axis.text=element_text(size=13),
        legend.title=element_text(size=14), legend.text=element_text(size=13),
        strip.text=element_text(size=13), legend.position="bottom",
        legend.margin=margin(t=-10),  # Reduce space between plots and legend
        legend.box.margin=margin(t=-10), # Further reduce the margin
        panel.spacing = unit(0.3, "lines"),  # Reduce facet spacing
        panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank()) +
  labs(fill="Model") +
  guides(fill=guide_legend(nrow=1)) +
  ylab("RMSE") + xlab("")

RMSE_comp %>% group_by(method,group) %>% summarise(mean=mean(RMSE))

Coverage_comp_m15_v7$group <- "m15v7"
Coverage_comp_m15_v3$group <- "m15v3"
Coverage_comp_m3_v7$group <- "m3v7"
Coverage_comp_m3_v3$group <- "m3v3"
Coverage_comp <- rbind(Coverage_comp_m15_v7,Coverage_comp_m15_v3,Coverage_comp_m3_v7,Coverage_comp_m3_v3)
Coverage_comp$group <- factor(Coverage_comp$group,levels=c("m15v7","m15v3","m3v7","m3v3"),
                              labels=c("15 large outliers","15 small outliers","3 large outliers","3 small outliers"))
#Figure 3 (c)
ggplot(data=Coverage_comp, aes(y=Coverage, fill=factor(method))) +
  geom_hline(yintercept=0.95,colour="black",linewidth=1) +
  geom_boxplot() + scale_x_continuous(breaks=NULL) +  
  facet_grid(~group) +  # Shrink the x-axis space
  scale_fill_manual(values=c("#CC79A7","#238b45","#74c476","#bae4b3","#1E88E5")) +
  theme(axis.title=element_text(size=14), axis.text=element_text(size=13),
        legend.title=element_text(size=14), legend.text=element_text(size=13),
        strip.text=element_text(size=13), legend.position="bottom",
        legend.margin=margin(t=-10),  # Reduce space between plots and legend
        legend.box.margin=margin(t=-10), # Further reduce the margin
        panel.spacing = unit(0.3, "lines"),  # Reduce facet spacing
        panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank())  +
  labs(fill="Model") +
  guides(fill=guide_legend(nrow=1)) +
  ylab("Coverage") + xlab("")

Coverage_comp %>% group_by(method,group) %>% summarise(median=median(Coverage))

MAE_comp_m15_v7$group <- "m15v7"
MAE_comp_m15_v3$group <- "m15v3"
MAE_comp_m3_v7$group <- "m3v7"
MAE_comp_m3_v3$group <- "m3v3"
MAE_comp <- rbind(MAE_comp_m15_v7,MAE_comp_m15_v3,MAE_comp_m3_v7,MAE_comp_m3_v3)
MAE_comp$group <- factor(MAE_comp$group,levels=c("m15v7","m15v3","m3v7","m3v3"),
                         labels=c("15 large outliers","15 small outliers","3 large outliers","3 small outliers"))
#Figure 3 (b)
ggplot(data=MAE_comp, aes(y=MAE, fill=factor(method))) +
  geom_boxplot() + scale_x_continuous(breaks=NULL) +  
  facet_grid(~group) +  # Shrink the x-axis space
  scale_fill_manual(values=c("#CC79A7","#238b45","#74c476","#bae4b3","#1E88E5")) +
  theme(axis.title=element_text(size=14), axis.text=element_text(size=13),
        legend.title=element_text(size=14), legend.text=element_text(size=13),
        strip.text=element_text(size=13), legend.position="bottom",
        legend.margin=margin(t=-10),  # Reduce space between plots and legend
        legend.box.margin=margin(t=-10), # Further reduce the margin
        panel.spacing = unit(0.3, "lines"),  # Reduce facet spacing
        panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank()) +
  labs(fill="Model") +
  guides(fill=guide_legend(nrow=1)) +
  ylab("MAE") + xlab("")

MeanCrIwidth_comp_m15_v7$group <- "m15v7"
MeanCrIwidth_comp_m15_v3$group <- "m15v3"
MeanCrIwidth_comp_m3_v7$group <- "m3v7"
MeanCrIwidth_comp_m3_v3$group <- "m3v3"
MeanCrIwidth_comp <- rbind(MeanCrIwidth_comp_m15_v7,MeanCrIwidth_comp_m15_v3,MeanCrIwidth_comp_m3_v7,MeanCrIwidth_comp_m3_v3)
MeanCrIwidth_comp$group <- factor(MeanCrIwidth_comp$group,levels=c("m15v7","m15v3","m3v7","m3v3"),
                                  labels=c("15 large outliers","15 small outliers","3 large outliers","3 small outliers"))
#Figure 3 (d)
ggplot(data=MeanCrIwidth_comp, aes(y=MeanCrIwidth, fill=factor(method))) +
  geom_boxplot() + scale_x_continuous(breaks=NULL) +  
  facet_grid(~group) +  # Shrink the x-axis space
  scale_fill_manual(values=c("#CC79A7","#238b45","#74c476","#bae4b3","#1E88E5")) +
  theme(axis.title=element_text(size=14), axis.text=element_text(size=13),
        legend.title=element_text(size=14), legend.text=element_text(size=13),
        strip.text=element_text(size=13), legend.position="bottom",
        legend.margin=margin(t=-10),  # Reduce space between plots and legend
        legend.box.margin=margin(t=-10), # Further reduce the margin
        panel.spacing = unit(0.3, "lines"),  # Reduce facet spacing
        panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank()) +
  labs(fill="Model") +
  guides(fill=guide_legend(nrow=1)) +
  ylab("Mean credible interval width") + xlab("")

