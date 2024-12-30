#load the relevant packages
library(ggplot2)
library(GGally)
library(sf)
library(readxl)
library(dplyr)
library(leaflet)
library(RColorBrewer)
library(tidyr)
library(spdep)
library(gridExtra)
library(cluster)
library(purrr)
library(broom)
library(xtable)
library(INLA)

sf_use_s2(FALSE)

#set working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#the data are available at 
#https://pldr.org/dataset/2noyv/small-area-mental-health-index-samhi
load("Data - England LSOA 4 diseases.RData") 

asthma <- dat[which(dat$Year=="2017"),]
asthma <- asthma[,c(2:4)]
colnames(asthma) <- c("area","population","incidence")
asthma$proportion <- asthma$incidence/asthma$population

#Figure 1
ggplot(data=asthma,aes(x=proportion,group=1)) + geom_histogram(fill="#33a02c",colour="black",alpha=0.7) +
  ylab("Count") + xlab("Observed prevalence") +
  theme(axis.title=element_text(size=12),axis.text=element_text(size=11))

#numerical summaries
mean(asthma$proportion)
median(asthma$proportion)
sd(asthma$proportion)

#the LSOA shapefiles are available at
#https://data.cambridgeshireinsight.org.uk/dataset/output-areas/resource/3bc4faef-38c7-417d-88a9-f302ad845ebe
Bdry <- st_read("Lower_Layer_Super_Output_Areas_December_2011_Generalised_Clipped__Boundaries_in_England_and_Wales.shp")
Bdry2 <- st_transform(x=Bdry, crs='+proj=longlat +datum=WGS84 +no_defs')
Asthma.Bdry <- merge(x=Bdry2,y=asthma,by.x="lsoa11cd",by.y="area",all.x=FALSE)

#derive neighbourhood matrix for bordering/adjacent areas
W.nb_adj <- poly2nb(Asthma.Bdry,row.names=Asthma.Bdry$lsoa11cd)
W.list_adj <- nb2listw(W.nb_adj,style="B")
summary(W.nb_adj) #average number of neighbours is 5.87

K <- nrow(Asthma.Bdry) #number of areas in the study region

#Moran's I test
coords <- st_centroid(st_geometry(Asthma.Bdry))
knn_moran <- 6 #number of neighbours to compute the Moran's I statistic
nb.test_moran <- knn2nb(knearneigh(coords, k=knn_moran), row.names=Asthma.Bdry$lsoa11cd)
W.moran_symm <- make.sym.nb(nb.test_moran) #add additional neighbours
W6nn_list_symm <- nb2listw(W.moran_symm,style="B")

mori <- moran.mc(x=as.numeric(Asthma.Bdry$proportion),listw=W6nn_list_symm,nsim=10000)
mori$statistic #0.89
mori$p.value #less than 0.0001

#create neighbourhoods for RDOS
knn_RDOS <- 10 #number of neighbours to compute the RDOS values
nb.test_RDOS <- knn2nb(knearneigh(coords, k=knn_RDOS), row.names=Asthma.Bdry$lsoa11cd)
W.RDOS_symm <- make.sym.nb(nb.test_RDOS) #add additional neighbours
W10nn_symm <- nb2mat(W.RDOS_symm, style = "B")
W10nn_list_symm <- nb2listw(W.RDOS_symm,style="B")

S_10_symm <- list()
for(i in 1:K){
  S_10_symm[[i]] <- which(W10nn_symm[i,]==1) #10-nearest neighbours, symmetric
}

#Gaussian kernel
Gk <- function(x,S,h){
  1/sqrt(2*pi)*(1/h)*exp(-1/2*((x-S)^2)/(h^2))
}

#matrix used to compute RDOS values
S <- S_10_symm
absdiff <- c() #local median absolute difference
for(i in 1:K){
  absdiff <- c(absdiff,median(abs(Asthma.Bdry$proportion[S[[i]]]-Asthma.Bdry$proportion[i])))
}

#initial bandwidth value, to be multiplied by constant c
h1 <- median(absdiff)

h <- h1
#initial density estimates for h1
est_dens_context <- rep(0,K) #estimated density
for(i in 1:K){
  est_dens_context[i] <- 1/(length(S[[i]])+1)*(sum(Gk(x=Asthma.Bdry$proportion[i],S=Asthma.Bdry$proportion[S[[i]]],h=h))+Gk(x=Asthma.Bdry$proportion[i],S=Asthma.Bdry$proportion[i],h=h))
}

#initial RDOS values for h1
RDOS_context <- rep(0,K) #RDOS values
for(i in 1:K){
  RDOS_context[i] <- mean(est_dens_context[S[[i]]])/(est_dens_context[i])
}
Asthma.Bdry$RDOS_context <- RDOS_context

#increase c until Kendall's rank coefficient exceeds 0.99
tau <- 0
c <- 1
while(tau<0.99){
  c <- c+0.1
  h <- c*h1
  est_dens_context <- rep(0,K) #estimated density
  for(i in 1:K){
    est_dens_context[i] <- 1/(length(S[[i]])+1)*(sum(Gk(x=Asthma.Bdry$proportion[i],S=Asthma.Bdry$proportion[S[[i]]],h=h))+Gk(x=Asthma.Bdry$proportion[i],S=Asthma.Bdry$proportion[i],h=h))
  }
  
  RDOS_context <- rep(0,K) #RDOS values
  for(i in 1:K){
    RDOS_context[i] <- mean(est_dens_context[S[[i]]])/(est_dens_context[i])
  }
  Asthma.Bdry$RDOS_context_new <- RDOS_context
  
  #all observations in decreasing order of RDOS value
  rdos_data_context <- as.data.frame(Asthma.Bdry[order(Asthma.Bdry$RDOS_context,decreasing=TRUE),])
  rdos_data_context$ranking <- 1:K 
  
  #ranks under the new RDOS values
  rdos_data_context_new <- as.data.frame(rdos_data_context[order(rdos_data_context$RDOS_context_new,decreasing=TRUE),])
  rdos_data_context_new$ranking_new <- 1:K
  
  #compare the ranks for the two iterations using Kendall's rank coefficient
  kendall_test <- cor.test(rdos_data_context_new$ranking,rdos_data_context_new$ranking_new,method="kendall")
  tau <- kendall_test$estimate
  Asthma.Bdry$RDOS_context <- Asthma.Bdry$RDOS_context_new
}

saveRDS(Asthma.Bdry,"Asthma_Bdry_RDOS.rds")
c #the scalar c takes on a value of 3.2

#rank the data in decreasing order of RDOS
rdos_data <- as.data.frame(Asthma.Bdry[order(Asthma.Bdry$RDOS_context,decreasing=TRUE),])
rdos_data$ranking <- 1:K

max(rdos_data_context$RDOS_context) #max value of 10.75
#Figure 5 (a)
ggplot(rdos_data, aes(x = as.numeric(ranking), y = RDOS_context)) + 
  geom_point() + xlab("Rank") + ylab("RDOS") +
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16)) +
  scale_x_continuous(breaks = c(1, seq(4000, 32000, 4000)), 
                     minor_breaks = c(1, seq(2000, 32000, 2000)), 
                     expand = c(0.02, 0.02)) +
  scale_y_continuous(breaks = c(1:11))

#Figure 5 (b)
ggplot(rdos_data[1:2000,],aes(x=factor(ranking),y=RDOS_context)) + geom_point() +
  xlab("Rank") + ylab("RDOS") +
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16)) +
  scale_x_discrete(breaks=c(1,seq(200,2000,200)),expand = c(0.02, 0.02)) +
  scale_y_continuous(breaks=c(1:11))

#Figure 5 (c)
ggplot(rdos_data[1:400,],aes(x=factor(ranking),y=RDOS_context)) + 
  geom_vline(xintercept=37.5,linetype="dashed",linewidth=0.8,alpha=0.5) +
  geom_vline(xintercept=116.5,linetype="dashed",linewidth=0.8,alpha=0.5) +
  geom_point() +
  xlab("Rank") + ylab("RDOS") + 
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16)) +
  scale_x_discrete(breaks=c(1,seq(20,400,20)),
                   expand = c(0.02, 0.02)) +
  scale_y_continuous(breaks=c(1:11))

#Figure 5 (d)
ggplot(rdos_data[20:200,], aes(x=factor(ranking), y=RDOS_context)) + 
  geom_vline(xintercept=which(levels(factor(rdos_data$ranking[20:200])) == 37)+0.3, 
             linetype="dashed", linewidth=0.8, alpha=0.5) +
  geom_vline(xintercept=which(levels(factor(rdos_data$ranking[20:200])) == 116)+0.3, 
             linetype="dashed", linewidth=0.8, alpha=0.5) +
  geom_point() + xlab("Rank") + ylab("RDOS") + 
  theme(axis.text=element_text(size=14), axis.title=element_text(size=16)) +
  scale_x_discrete(breaks=seq(20, 200, 20), expand=c(0.02, 0.02)) +
  scale_y_continuous(breaks=1:11)

#Rug density plot (supplementary material)
ggplot(data=Asthma.Bdry, aes(x=RDOS_context)) +
  geom_density() +  
  geom_rug(sides="b", alpha=0.5) +
  xlab("RDOS") + ylab("Density") +
  theme(axis.title=element_text(size=13),axis.text=element_text(size=12)) +
  geom_segment(aes(x=3.828,xend=3.828,y=0,yend=10),linetype = "dashed") +
  geom_segment(aes(x=2.695,xend=2.695,y=0,yend=10),linetype = "dashed") +
  scale_x_continuous(breaks=seq(1,10,1))

#check RDOS vs observed proportions for 37 identified outliers
n_outl <- 37
Asthma.Bdry$nb_comp <- NA
for(i in 1:K){
  if(Asthma.Bdry$proportion[i] < mean(Asthma.Bdry$proportion[S_10_symm[[i]]])){
    Asthma.Bdry$nb_comp[i] <- "#009E73"
  }else{
    Asthma.Bdry$nb_comp[i] <- "#E69F00"
  }
}

rdos_data <- as.data.frame(Asthma.Bdry[order(Asthma.Bdry$RDOS_context,decreasing=TRUE),])
rdos_data$nb_comp[38:K] <- "black"
rdos_data$nb_comp <- factor(rdos_data$nb_comp,levels=c("black","#009E73","#E69F00"))
#Figure 6
ggplot() + geom_point(data=rdos_data,aes(x=proportion,y=RDOS_context,colour=nb_comp)) +
  scale_colour_manual(values=c("black","#009E73","#E69F00"),labels=c("Inlier","Low local outlier","High local outlier")) +
  xlab("Observed prevalence") + ylab("RDOS") +
  theme(legend.title=element_blank(),axis.text=element_text(size=12),axis.title=element_text(size=13),
        legend.text=element_text(size=12),title=element_text(size=13)) + 
  scale_y_continuous(breaks=seq(from=1,to=11,by=2)) 

#Modelling
coords <- st_centroid(st_geometry(Asthma.Bdry))
knn_INLA <- 6 #number of neighbours to use in the INLA model
nb.test_INLA <- knn2nb(knearneigh(coords, k=knn_INLA), row.names=Asthma.Bdry$lsoa11cd)
W.INLA_symm <- make.sym.nb(nb.test_INLA)
nb2INLA("map.6nn_all",W.INLA_symm)
W_all <- inla.read.graph(filename="map.6nn_all")

Asthma.Bdry$idarea <- 1:K
#conventional BYM2 model
res_no_outl <- inla(incidence ~ f(idarea,model="bym2",graph=W_all,scale.model=TRUE),
                    family="binomial",
                    data=Asthma.Bdry,
                    Ntrials=population,
                    control.compute=list(dic=TRUE,waic=TRUE,cpo=TRUE),
                    control.fixed = list(prec.intercept=0.001))

#modified BYM2-O model with 37 identified outliers
n_outl <- 37
W.INLA_symm_inl <- subset(W.INLA_symm,!(1:length(W.INLA_symm) %in% which(Asthma.Bdry$lsoa11cd %in% rdos_data_context$lsoa11cd[1:n_outl])))
nb2INLA("map.6nn_inl_37out",W.INLA_symm_inl)
W_inl_37 <- inla.read.graph(filename="map.6nn_inl_37out")

j <- 1
id_inliers <- c(rep(NA,K))
for(i in 1:K){
  id_inliers[i] <- ifelse(i %in% which(Asthma.Bdry$lsoa11cd %in% rdos_data_context$lsoa11cd[1:n_outl]),NA,j)
  if(i %in% which(Asthma.Bdry$lsoa11cd %in% rdos_data_context$lsoa11cd[(n_outl+1):K])){
    j <- j+1
  }
}

j <- 1 
id_outliers <- c(rep(NA,K))
for(i in 1:K){
  id_outliers[i] <- ifelse(i %in% which(Asthma.Bdry$lsoa11cd %in% rdos_data_context$lsoa11cd[1:n_outl]),j,NA)
  if(i %in% which(Asthma.Bdry$lsoa11cd %in% rdos_data_context$lsoa11cd[1:n_outl])){
    j <- j+1
  }
}
Asthma.Bdry$id_inliers <- id_inliers 
Asthma.Bdry$id_outliers <- id_outliers

res_outl_bym2_37outl <- inla(incidence ~ f(id_inliers,model="bym2",graph=W_inl_37,scale.model=TRUE) +
                               f(id_outliers,model="iid"),
                             family="binomial",
                             data=Asthma.Bdry,
                             Ntrials=population,
                             control.compute=list(dic=TRUE,waic=TRUE,cpo=TRUE),
                             control.fixed = list(prec.intercept=0.001))

