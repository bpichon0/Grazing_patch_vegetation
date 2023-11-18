rm(list=ls())
source("./Structure_grazing_function.R")

# ---------------------- Step 0: Sensitivity of change in cover ---- 
## >> 1) How do spatial statistics change with uncertainty of cover ----

d=tibble();thresh=0.01
n_site=100

for (id in 1:n_site){ # Id of the site
    
  for (n_clust in c(2,3,4)){ # number of kmeans clusters
    
    for (binary_thresh in 1:(n_clust-1)){ # binary thresholds
      
      #getting the images and the kept matrix = "true" one 
      cover_site=d_biocom_old$Cover[id]
      img1=Get_png_empirical_site_biocom(id)
      kept_land=Get_empirical_site_biocom(id)
      
      #performing the kmean
      kmean_land=k_means_RGB(img1,n_clust)
      
      #defining the thresholds to binarize the matrix
      if (n_clust>2){
        val_scale=get_cut_grayscale_values(n_clust)[[binary_thresh]]
        site_land=binarize(kmean_land,val_scale[[1]],val_scale[[2]])
      }else {
        if (kept_land[1,1]!=kmean_land[1,1]){ #we order K-means classification
          kmean_land[kmean_land==0]=2
          kmean_land[kmean_land==1]=0
          kmean_land[kmean_land==2]=1
          site_land=kmean_land
        }
      }
      cover_img_binary=sum(site_land==1)/(nrow(site_land)**2) #cover of the binarized image
      
      if ((cover_img_binary-d_biocom_old$Cover[id])**2<thresh){ #we only keep the transformations close in term of cover
        d=rbind(d,Get_sumstat(site_land==1)%>%
                  add_column(., Nclust=n_clust,binary_thresh=binary_thresh,ID_site=id,Type="Changed"))
      } 
    }
  }
}
write.table(d,"../Data/Sensitivity_cover/Varying_threshold.csv",sep=";")

d_tot=rbind(d,d_biocom_old[1:100,14:24]%>%add_column(Nclust=0,binary_thresh=0,
                                                    ID_site=1:100,Type="Kept"))%>%
  arrange(., ID_site,Nclust)%>%
  add_column(., Error_cover=sapply(1:nrow(.),function(x){
    return((.$rho_p[x]-d_biocom_old$rho_p[.$ID_site[x]])**2)
}))


#Doing a PCA

sumstat_name=colnames(d_tot)[1:11]
res.comp=imputePCA(d_tot[,which(colnames(d_tot) %in% sumstat_name)],ncp=3,scale = T) 

if ("completeObs" %in% names(res.comp)){
  res.pca=PCA(res.comp$completeObs, ncp = 3,  graph=F)
}else {
  res.pca=PCA(res.comp, ncp = 3,  graph=F)
}

axes_for_plot=tibble(x=c(1,1,2),y=c(2,3,3))

for (i in 1:3){
  assign(paste0("p",i),
         d_tot%>%
           mutate(., Error_cover=unlist(sapply(1:nrow(.),function(x){
             if (.$Error_cover[x]==0){
               return(NA)
             }else { 
               return(.$Error_cover[x])}
           })))%>%
           add_column(., PC1=res.pca$ind$coord[,axes_for_plot$x[i]],PC2=res.pca$ind$coord[,axes_for_plot$y[i]])%>%
           ggplot(.) +
           geom_hline(yintercept = 0, lty = 2) +
           geom_vline(xintercept = 0, lty = 2) +
           geom_point(aes(x = PC1, y = PC2, color = Error_cover,fill=Error_cover,shape=Type))+
           geom_line(aes(x = PC1, y = PC2, group=ID_site),color="gray",alpha=.5)+
           scale_color_viridis_c(na.value = "red")+
           scale_fill_viridis_c(na.value = "red")+
           scale_shape_manual(values = c("Kept"=25,"Changed"=21))+
           labs(x=paste0("PC ",axes_for_plot$x[i]," (",round(res.pca$eig[axes_for_plot$x[i],2], 1)," %)"),
                y=paste0("PC ",axes_for_plot$y[i]," (",round(res.pca$eig[axes_for_plot$y[i],2], 1)," %)"),
                color="Quadratic error on cover",fill="")+
           ggtitle("")+guides()+
           theme_classic()+theme(legend.position = "bottom")+
           guides(fill="none",shape="none")+
           theme(legend.box = "vertical")
  )
}

p=ggarrange(ggarrange(p1+theme(legend.position = "none"),
                      p2+theme(legend.position = "none"),
                      p3+theme(legend.position = "none"),
                      ncol=3,align = "hv"),ggarrange(ggplot()+theme_void(),get_legend(p2),ggplot()+theme_void(),ncol=3,widths = c(.3,1,.3)),
            nrow=2,heights = c(1,.2))

ggsave("../Figures/PCA_error.pdf",p,width = 11,height = 5)





## >> 2) Change in parameter inferred ----

NA_kept=100
id=1
d_biocom_old=read.table("../Data/Sensitivity_cover/Varying_threshold.csv",sep=";")
d_sim=read.table("../Data_new/All_new_sim2.csv",sep=";")%>%
  dplyr::relocate(., Pooling,.after =q )%>%
  #filter(., PL_expo>0,!is.na(PLR),ID %in% c(1:15))%>%
  filter(., PL_expo>0,!is.na(PLR))%>%
  dplyr::select(., -ID)

rownames(d_sim)=1:nrow(d_sim)


for (id_plot in id){
  `%!in%` = Negate(`%in%`)
  n_param=3
  
  
  d_param_infer_NN=array(0,c(NA_kept,nrow(d_biocom_old),3))
  d_param_infer_rej=array(0,c(NA_kept,nrow(d_biocom_old),3))
  
  d_NRMSE_sumstat=x_y_stat=tibble()
  
  name_plot=c("all","no_cv","no_cv_fmax","no_PL_PLR","not_the_4","only_fmax","no_PL","no_PLR")
  
  sumstat_kept=4:14
  
  
  
  for (empirical_id in 1:nrow(d_biocom_old)){
    print(empirical_id)
    target=d_biocom_old[empirical_id,sumstat_kept-2]
    matrix_param=d_sim[,1:3]
    
    mat_sumstat=rbind(d_sim[,sumstat_kept],target)
    
    #1) Boxcox
    
    for (x in 1:ncol(mat_sumstat)){
      if (any(is.na(target)) & x %!in% which(is.na(target))){
        if (colnames(mat_sumstat)[x] %in% c("skewness","moran_I","fmax_psd")){
          
          
          
          b=boxcox(lm(mat_sumstat[,x]+abs(min(mat_sumstat[,x]))+.5 ~ 1),plotit = F,eps = .05)     #Working with positive values
          lambda_x=b$x[which.max(b$y)]
          if (lambda_x !=0){ #to avoid errors
            mat_sumstat[,x] = (exp(mat_sumstat[,x]*(lambda_x)) -1)/(lambda_x)
          }
          
        }else {
          b=boxcox(lm(mat_sumstat[,x]+.5 ~ 1),plotit = F,eps = .05)
          lambda_x=b$x[which.max(b$y)]
          if (lambda_x !=0){ #to avoid errors
            mat_sumstat[,x] = (mat_sumstat[,x]^(lambda_x) -1)/(lambda_x)
          }
        }
      }
    }
    
    
    #2) Scaling
    
    for (x in 1:ncol(mat_sumstat)) mat_sumstat[,x] = (mat_sumstat[,x]-mean(mat_sumstat[,x],na.rm = T))/sd(mat_sumstat[,x],na.rm = T)
    
    if (any(is.na(mat_sumstat[nrow(mat_sumstat),]))){
      
      which_na=which(is.na(mat_sumstat[nrow(mat_sumstat),]))
      
      cross_valid=abc(target = mat_sumstat[nrow(mat_sumstat),-which_na],
                      param = matrix_param[-nrow(mat_sumstat),],sumstat = mat_sumstat[-nrow(mat_sumstat),-which_na], #removing the target data
                      tol = 1000/nrow(matrix_param),method = "rejection") #we keep the 1000 closest simulations for the first step
      
    }else {
      cross_valid=abc(target = mat_sumstat[nrow(mat_sumstat),],
                      param = matrix_param[-nrow(mat_sumstat),],sumstat = mat_sumstat[-nrow(mat_sumstat),], #removing the target data
                      tol = 1000/nrow(matrix_param),method = "rejection") #we keep the 1000 closest simulations for the first step
    }
    
    
    #Keeping 1000 simulations and doing the same steps again: normality, scaling and PLS
    
    mat_sumstat_step1=d_sim[as.numeric(rownames(cross_valid$ss)),sumstat_kept] #we keep information with the true values
    mat_sumstat_step1=rbind(mat_sumstat_step1,target)
    
    #again, first box cox
    which_na=which(is.na(mat_sumstat[nrow(mat_sumstat),]))
    
    for (x in 1:ncol(mat_sumstat_step1)){
      
      if (any(which_na) & x %!in% which_na){
        
        if (colnames(mat_sumstat_step1)[x] %in% c("skewness","moran_I","fmax_psd")){
          
          b=boxcox(lm(mat_sumstat_step1[,x]+abs(min(mat_sumstat_step1[,x]))+.5 ~ 1),plotit = F,eps = .05)     #Working with positive values
          lambda_x=b$x[which.max(b$y)]
          
          if (lambda_x !=0){ #to avoid errors
            mat_sumstat_step1[,x] = (exp(mat_sumstat_step1[,x]*(lambda_x)) -1)/(lambda_x)
          }
          
        }else {
          b=boxcox(lm(mat_sumstat_step1[,x] ~ 1),plotit = F,eps = .05)
          lambda_x=b$x[which.max(b$y)]
          if (lambda_x !=0){ #to avoid errors
            mat_sumstat_step1[,x] = (mat_sumstat_step1[,x]^(lambda_x) -1)/(lambda_x)
          }
        }
      }
    }
    
    #and normalization
    for (x in 1:ncol(mat_sumstat_step1)) {
      if (length(unique(mat_sumstat_step1[,x])) != 1){
        mat_sumstat_step1[,x] = (mat_sumstat_step1[,x]-mean(mat_sumstat_step1[,x],na.rm = T))/sd(mat_sumstat_step1[,x],na.rm = T)
      }
    }
    if (any(is.na(mat_sumstat_step1[nrow(mat_sumstat_step1),]))){
      
      which_na=which(is.na(mat_sumstat_step1[nrow(mat_sumstat_step1),]))
      
      cross_valid=abc(target = mat_sumstat_step1[nrow(mat_sumstat_step1),-which_na],
                      param = cross_valid$unadj.values,
                      sumstat = mat_sumstat_step1[-nrow(mat_sumstat_step1),-which_na], #removing the target data
                      tol = NA_kept/nrow(mat_sumstat_step1),method = "rejection",transf = c(rep("logit",2),"none"), #as parameters are proba, we perform logit regression
                      logit.bounds = matrix(c(0,1),3,2,byrow = T),
                      numnet = 10,sizenet = 15)
      
      cross_valid$ss=d_sim[as.numeric(rownames(cross_valid$ss)),sumstat_kept] #we keep information with the true values
      
      if(name_plot[id_plot]=="no_PLR"){
        
        x_y_stat=rbind(x_y_stat,target[-which_na]%>%
                         add_column(., PL_expo=NA, Site_ID=empirical_id,Method="rejection",Type="Obs")%>%
                         relocate(., PL_expo,.after =Spectral_ratio ))
        
        x_y_stat=rbind(x_y_stat,as_tibble(t(colMeans(cross_valid$ss)[-which_na]))%>%
                         add_column(.,PL_expo=NA, Site_ID=empirical_id,Method="rejection",Type="Sim")%>%
                         relocate(., PL_expo,.after =Spectral_ratio ))
      }else{
        
        x_y_stat=rbind(x_y_stat,target[-which_na]%>%
                         add_column(., PL_expo=NA, Site_ID=empirical_id,Method="rejection",Type="Obs")%>%
                         relocate(., PL_expo,.after =PLR ))
        
        x_y_stat=rbind(x_y_stat,as_tibble(t(colMeans(cross_valid$ss)[-which_na]))%>%
                         add_column(.,PL_expo=NA, Site_ID=empirical_id,Method="rejection",Type="Sim")%>%
                         relocate(., PL_expo,.after =PLR ))
        
        
      }
      
      
    }else {
      cross_valid=abc(target = mat_sumstat_step1[nrow(mat_sumstat_step1),],
                      param = cross_valid$unadj.values,
                      sumstat = mat_sumstat_step1[-nrow(mat_sumstat_step1),], #removing the target data
                      tol = NA_kept/nrow(mat_sumstat_step1),method = "rejection",transf = c(rep("logit",2),"none"), #as parameters are proba, we perform logit regression
                      logit.bounds = matrix(c(0,1),3,2,byrow = T),
                      numnet = 10,sizenet = 15)
      
      cross_valid$ss=d_sim[as.numeric(rownames(cross_valid$ss)),sumstat_kept] #we keep information with the true values
      
      x_y_stat=rbind(x_y_stat,target%>%
                       add_column(., Site_ID=empirical_id,Method="rejection",Type="Obs"))
      
      x_y_stat=rbind(x_y_stat,as_tibble(t(colMeans(cross_valid$ss)))%>%
                       add_column(.,Site_ID=empirical_id,Method="rejection",Type="Sim"))
      
    }
    
    
    cross_valid$ss=d_sim[as.numeric(rownames(cross_valid$ss)),sumstat_kept] #we keep information with the true values
    
    mat_sumstat=d_sim[,sumstat_kept]
    
    cross_valid$adj.values=cross_valid$adj.values
    
    #NRMSE for the summary statistics observed
    RMSE = sapply(1:ncol(cross_valid$ss),function(x){
      sqrt(sum((cross_valid$ss[,x]-target[,x])**2,na.rm = T)/nrow(cross_valid$ss) )
    }
    )
    
    RMSE_prior=sapply(1:ncol(mat_sumstat),function(x){
      sqrt(sum((mat_sumstat[,x]-target[,x])**2,na.rm = T)/nrow(mat_sumstat) )
    }
    )
    NRMSE = RMSE/RMSE_prior
    
    d_NRMSE_sumstat=rbind(d_NRMSE_sumstat,as_tibble(t(NRMSE)))
    
    d_param_infer_rej[,empirical_id,1]=cross_valid$unadj.values[,1] # we keep the whole distribution for p
    d_param_infer_rej[,empirical_id,2]=cross_valid$unadj.values[,2] # for q
    d_param_infer_rej[,empirical_id,3]=cross_valid$unadj.values[,3] # for the scale of observation
    
    
  }
  
  write.table(d_NRMSE_sumstat,paste0("../Data/Sensitivity_cover/NRMSE_sumstat.csv"),sep=";")
  
  write.table(x_y_stat,paste0("../Data/Sensitivity_cover/x_y_stat.csv"),sep=";")
  
  write.table(d_param_infer_rej,paste0("../Data/Sensitivity_cover/param_inferred.csv"),sep=";")
}


nsite=100
param_infered=read.table("../Data/Sensitivity_cover/param_inferred.csv",sep=";")
param_infered_biocom=read.table("../Data/Sensitivity_cover/posterior_param.csv",sep=";")
d_threshold=read.table("../Data/Sensitivity_cover/Varying_threshold.csv",sep=";")

d_inferrence=tibble(p=c(colMeans(param_infered_biocom[1:nsite]),colMeans(param_infered)[1:nrow(d_threshold)]),
                    q=c(colMeans(param_infered_biocom[c(1:nsite)+345]),colMeans(param_infered)[(nrow(d_threshold)+1):(2*nrow(d_threshold))]),
                    Site=c(1:nsite,d_threshold$ID_site),
                    Type=c(rep(1,nsite),d_threshold$Nclust))


p=ggplot(d_inferrence%>%melt(.,measure.vars=c("p","q")))+
  geom_line(aes(Type,value,group=Site),color="gray",lwd=.5,alpha=.5)+
  facet_wrap(.~variable)+
  scale_x_continuous(labels = c("True","Modif_1","Modif_2","Modif_3"))+
  the_theme


ggsave("../Figures/Change_parameters_no_info_cover.pdf",p,width = 8,height = 4)



# ---------------------- Step 1: Computing the spatial metrics ----
## >> 1) Transforming images into binary matrices ----


dir.create("../Data/Landscapes/Binary_landscapes/",showWarnings = F)

Extract_binary_matrix=function(id){  
  
  info_kmean=read.table("../Data/Landscapes/Kmean_clust_info.csv",sep=",",header = T)%>%
    filter(., Size!=200,Dataset=="biodesert") #keeping the kept sites
  
  # we load the landscape
  img=readJPEG(paste0("../Data/Landscapes/",info_kmean$Dataset[id],
                      "/",info_kmean$Size[id],"/",
                      info_kmean$Site[id],"_",info_kmean$Image[id],".jpeg")) 
  
  #and binarize the kmean output it with Josquin thresholds
  kmean_img=k_means_RGB(img,info_kmean$nclust[id])
  
  cats=get_cut_grayscale_values(info_kmean$nclust[id])[[info_kmean$cut[id]]]
  mat=kmean_img %>% binarize(cats[[1]], cats[[2]])

  #saving the binary matrix
  write.table(mat,paste0("../Data/Landscapes/Binary_landscapes/",info_kmean$Dataset[id],
                         "_",info_kmean$Size[id],"_",
                         info_kmean$Site[id],"_",info_kmean$Image[id],".csv"),
              row.names = F,col.names = F,sep=",")
}

library(parallel)
mclapply(Extract_binary_matrix,1:978,mc.cores = 25)

## >> 2) Computing the metrics on the binary landscapes (of Biodesert only) ----

## Computing the metrics

dir.create("../Data/Metrics",showWarnings = F)

Compute_all_metrics_biodesert=function(k){
  
  landscape=as.matrix(read.table(paste0("../Data/Landscapes/Binary_landscapes/",k),sep=","))
  name=gsub(".csv","",k)
  name_dataset=strsplit(k,"_")[[1]][1]
  
  if (any(grep("cropped",k))){
    id_site=as.numeric(strsplit(k,"_")[[1]][4])
    sub_id=as.numeric(gsub(".csv","",(strsplit(k,"_")[[1]][5])))
  }else{
    id_site=as.numeric(strsplit(k,"_")[[1]][3])
    sub_id=as.numeric(gsub(".csv","",(strsplit(k,"_")[[1]][4])))
  }
  
  slope=d_biodesert$`SLOPE-ALOS30`[which(d_biodesert$ID==id_site)]
  
  d=Get_sumstat(landscape,slope = slope)%>% #Computing all summary statistics
    add_column(.,Full_name=name,Resolution=Get_spatial_resolution(landscape), #full name & spatial resolution
               Site_ID=id_site, # ID of the site
               Sub_id=sub_id) # id of the sublandscape taken at a given site
  
  write.table(d,paste0("../Data/Metrics/Metric_",name,".csv"),sep=";")  
}

list_f=list.files("../Data/Landscapes/Binary_landscapes/",".csv")


library(parallel)
mclapply(list_f[grep(list_f,pattern = "biodesert")],Compute_all_metrics_biodesert,mc.cores = 5)


## aggregating the metrics into a df

info_kmean=read.table("../Data/Landscapes/Kmean_clust_info.csv",sep=",",header = T)%>%
  filter(., Size!=200,Dataset=="biodesert",status=="kept") #keeping the kept sites

d=tibble()
for (k in paste0("Metric_",info_kmean$Site_ID,"_",info_kmean$Image,".csv")){
  d=rbind(d,read.table(paste0("../Data/Metrics/",k),sep=";",header = T))}

d=d%>%
  tibble::add_column(.,
             Grazing=sapply(1:nrow(.),function(x){return(d_biodesert$GRAZ[which(d_biodesert$ID==.$Site_ID[x])])}),
             Aridity=sapply(1:nrow(.),function(x){return(d_biodesert$ARIDITY[which(d_biodesert$ID==.$Site_ID[x])])}),
             Nurse=sapply(1:nrow(.),function(x){return(d_biodesert$Nurse[which(d_biodesert$ID==.$Site_ID[x])])}),
             Slope=sapply(1:nrow(.),function(x){return(d_biodesert$`SLOPE-ALOS30`[which(d_biodesert$ID==.$Site_ID[x])])}),
             Elevation=sapply(1:nrow(.),function(x){return(d_biodesert$`ELE-ALOS30`[which(d_biodesert$ID==.$Site_ID[x])])}),
             Sp_richness=sapply(1:nrow(.),function(x){return(d_biodesert$SR[which(d_biodesert$ID==.$Site_ID[x])])}),
             Lattitude=sapply(1:nrow(.),function(x){return(d_biodesert$Lat_decimal[which(d_biodesert$ID==.$Site_ID[x])])}),
             Longitude=sapply(1:nrow(.),function(x){return(d_biodesert$Long_decimal[which(d_biodesert$ID==.$Site_ID[x])])}),
             Type_veg=sapply(1:nrow(.),function(x){return(d_biodesert$VEG[which(d_biodesert$ID==.$Site_ID[x])])}),
             Herbivores=sapply(1:nrow(.),function(x){return(d_biodesert$Dominant_Livestock[which(d_biodesert$ID==.$Site_ID[x])])}),
             Org_C=sapply(1:nrow(.),function(x){return(d_biodesert$`ORC veg`[which(d_biodesert$ID==.$Site_ID[x])]- #difference between plant and bare soil for organic carbon
                                                         d_biodesert$`ORC b`[which(d_biodesert$ID==.$Site_ID[x])])}),
             Sand=sapply(1:nrow(.),function(x){return(d_biodesert$`SAC veg`[which(d_biodesert$ID==.$Site_ID[x])])}),
             Shannon_div=sapply(1:nrow(.),function(x){return(d_biodesert$DIV_SHQ[which(d_biodesert$ID==.$Site_ID[x])])}),
             Simpson_div=sapply(1:nrow(.),function(x){return(d_biodesert$DIV_SIQ[which(d_biodesert$ID==.$Site_ID[x])])}),
             Woody=sapply(1:nrow(.),function(x){return(d_biodesert$RWCQ[which(d_biodesert$ID==.$Site_ID[x])])})
  )%>%
  dplyr::mutate(., Type_veg=recode_factor(Type_veg,"1"="Grassland","2"="Shrubland","3"="Forest","1_2"="Grass_Shrub"))


#Then, we extract summarized climatic variables using a PCA on all climatic variables

clim_variables=colnames(d_biodesert)[38:58]
res.comp=imputePCA(d_biodesert[,which(colnames(d_biodesert) %in% clim_variables)],ncp=4,scale = T) 

if ("completeObs" %in% names(res.comp)){
  res.pca=PCA(res.comp$completeObs, ncp = 4,  graph=F)
}else {
  res.pca=PCA(res.comp, ncp = 4,  graph=F)
}

fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 50))
p1=fviz_pca_var(res.pca,col.var=c(rep("Temperature",11),
                                  rep("Precipitation",8),
                                  rep("Temperature",2)),
                axes = c(1,2))+ #the ones that are the most imporant
  the_theme+
  scale_color_manual(values=c("#4761D0","#EF5454"))+
  theme(legend.position = "none")+
  labs(x=paste0("PC 1 (",round(res.pca$eig[1,2],1),")"),
       y=paste0("PC 2 (",round(res.pca$eig[2,2],1),")"))+
  ggtitle("")

p2=fviz_pca_var(res.pca,col.var=c(rep("Temperature",11),
                                  rep("Precipitation",8),
                                  rep("Temperature",2)),
                axes = c(3,4))+ #the ones that are the most imporant
  the_theme+
  scale_color_manual(values=c("#4761D0","#EF5454"))+
  theme(legend.position = "none")+
  labs(x=paste0("PC 3 (",round(res.pca$eig[3,2],1),")"),
       y=paste0("PC 4 (",round(res.pca$eig[4,2],1),")"))+
  ggtitle("")

ggsave("../Figures/Final_figs/SI/PCA_climatic_variables.pdf",ggarrange(p1,p2),
       width = 14,height = 5)

#We extract the first 4 ones

d=d%>%
  add_column(.,
             Clim1=sapply(1:nrow(.),function(x){return(res.pca$ind$coord[which(d_biodesert$ID==.$Site_ID[x]),1])}),
             Clim2=sapply(1:nrow(.),function(x){return(res.pca$ind$coord[which(d_biodesert$ID==.$Site_ID[x]),2])}),
             Clim3=sapply(1:nrow(.),function(x){return(res.pca$ind$coord[which(d_biodesert$ID==.$Site_ID[x]),3])}),
             Clim4=sapply(1:nrow(.),function(x){return(res.pca$ind$coord[which(d_biodesert$ID==.$Site_ID[x]),4])})
             )

#we add longitude (cos & sin) as well as functioning of the soil
d=d%>%
  add_column(., Long_sin=sin(.$Longitude),Long_cos=cos(.$Longitude))%>%
  add_column(., 
             Aromatic=as.numeric(sapply(1:nrow(.),function(x){return(d_biodesert$`ARO veg`[which(d_biodesert$ID==.$Site_ID[x])]
                                                                     # -
                                                                     #   d_biodesert$`ARO b`[which(d_biodesert$ID==.$Site_ID[x])]
             )})),
             Org_C_v=as.numeric(sapply(1:nrow(.),function(x){return(d_biodesert$`ORC veg`[which(d_biodesert$ID==.$Site_ID[x])]
                                                                  # -
                                                                  #   d_biodesert$`ORC b`[which(d_biodesert$ID==.$Site_ID[x])]
             )})),
             Beta_gluco=as.numeric(sapply(1:nrow(.),function(x){return(d_biodesert$`BGL veg`[which(d_biodesert$ID==.$Site_ID[x])]
                                                                       # -
                                                                       #   d_biodesert$`BGL b`[which(d_biodesert$ID==.$Site_ID[x])]
             )})),
             Phosphatase=as.numeric(sapply(1:nrow(.),function(x){return(d_biodesert$`FOS veg`[which(d_biodesert$ID==.$Site_ID[x])]
                                                                        # -
                                                                        #   d_biodesert$`FOS b`[which(d_biodesert$ID==.$Site_ID[x])]
             )})),
             Total_N=as.numeric(sapply(1:nrow(.),function(x){return(d_biodesert$`TON veg`[which(d_biodesert$ID==.$Site_ID[x])]
                                                                    # -
                                                                    #   d_biodesert$`TON b`[which(d_biodesert$ID==.$Site_ID[x])]
             )})),
             Amonium=as.numeric(sapply(1:nrow(.),function(x){return(d_biodesert$`AMO veg`[which(d_biodesert$ID==.$Site_ID[x])]
                                                                    -
                                                                      d_biodesert$`AMO b`[which(d_biodesert$ID==.$Site_ID[x])]
             )})),
             Nitrate=as.numeric(sapply(1:nrow(.),function(x){return(d_biodesert$`NIT veg`[which(d_biodesert$ID==.$Site_ID[x])]
                                                                    # -
                                                                    #   d_biodesert$`NIT b`[which(d_biodesert$ID==.$Site_ID[x])]
             )})),
             Hexose=as.numeric(sapply(1:nrow(.),function(x){return(d_biodesert$`HEX veg`[which(d_biodesert$ID==.$Site_ID[x])]
                                                                   # -
                                                                   #   d_biodesert$`HEX b`[which(d_biodesert$ID==.$Site_ID[x])]
             )})),
             N_transfo=as.numeric(sapply(1:nrow(.),function(x){return(d_biodesert$`NTR veg`[which(d_biodesert$ID==.$Site_ID[x])]
                                                                      # -
                                                                      #   d_biodesert$`NTR b`[which(d_biodesert$ID==.$Site_ID[x])]
             )})),
             N_minera=as.numeric(sapply(1:nrow(.),function(x){return(d_biodesert$`MIN veg`[which(d_biodesert$ID==.$Site_ID[x])]
                                                                     # -
                                                                     #   d_biodesert$`MIN b`[which(d_biodesert$ID==.$Site_ID[x])]
             )})),
             Total_P=as.numeric(sapply(1:nrow(.),function(x){return(as.numeric(c(d_biodesert$Total_P_vegetation[which(d_biodesert$ID==.$Site_ID[x])]))
                                                                    # -
                                                                    #   as.numeric(d_biodesert$Total_P_b[which(d_biodesert$ID==.$Site_ID[x])])
                                                                      )})),
             lnTotal_P=as.numeric(sapply(1:nrow(.),function(x){return(as.numeric(c(d_biodesert$Total_P_vegetation[which(d_biodesert$ID==.$Site_ID[x])]))
                                                                     -
                                                                       as.numeric(d_biodesert$Total_P_b[which(d_biodesert$ID==.$Site_ID[x])])
             )})),
             lnTotal_N=as.numeric(sapply(1:nrow(.),function(x){return(as.numeric(c(d_biodesert$`TON veg`[which(d_biodesert$ID==.$Site_ID[x])]))
                                                                      -
                                                                        as.numeric(d_biodesert$`TON b`[which(d_biodesert$ID==.$Site_ID[x])])
             )}))
  )

#transforming into true Site ID
d$Site_ID=sapply(1:nrow(d),function(x){return(d_biodesert$SITE_ID[which(d_biodesert$ID==d$Site_ID[x])])})

write.table(d,"../Data/Spatial_structure_grazing.csv",sep=";")

## >> 3) Spatial metrics on simulations ----

#simulating and computing spatial statistics
param_space=rbind(tibble(b=seq(0,.8,.03),g0=0,Driver="Aridity"),
                  tibble(b=.8,g0=seq(0,0.5,.02),Driver="Grazing"))

d_sim=tibble()
param=Get_classical_param(g0=0)

for (k in 1:nrow(param_space)){
  print(k)
  
  param$g0=param_space$g0[k]
  param$b=param_space$b[k]
  
  ini_land=Get_initial_lattice()
  out_model=Run_spatial_model(params = param,ini=ini_land)
  
  d_sim=rbind(d_sim,Get_sumstat(out_model$State==1)%>%
                add_column(., g0=param_space$g0[k],b=param_space$b[k],Driver=param_space$Driver[k]))
}
write.table(d_sim,"../Data/Spatial_structure_simulations.csv",sep=";")


# ---------------------- Step 2: Direct & interactive effects of grazing  ----


## >> 1) Running mixed-effect models  ----

d_data=read.table("../Data/Spatial_structure_grazing.csv",sep=";")%>%
  Closer_to_normality(.)

d_data[,c(1,9:29,36,38:42,45:ncol(d_data))] = apply(d_data[,c(1,9:29,36,38:42,45:ncol(d_data))],2,z_tranform)


#correlations between predictors ?
p=ggplot(cor(d_data[,c("Sand","Clim1","Clim2","Clim3","Clim4","Org_C","Woody",
                     "Grazing","rho_p","Long_cos","Long_sin","Lattitude","Slope","Elevation")],
           use =  "na.or.complete")%>%
         melt(.))+geom_tile(aes(Var1,Var2,fill=value))+
  scale_fill_gradient2(low="red","white","blue")+
  labs(x="",y="")+
  theme(axis.text.x=element_text(angle=60,hjust=1))
  
ggsave("../Figures/Step1_Understanding_grazing/Correlation_btw_predictors.pdf",p,width = 5,height = 4)




### Running the models with all climatic variables ----

dir.create("../Data/Step1_Understanding_grazing/VIF",showWarnings = F)
dir.create("../Data/Step1_Understanding_grazing/Auto_corr",showWarnings = F)
dir.create("../Data/Step1_Understanding_grazing/Importance",showWarnings = F)
dir.create("../Data/Step1_Understanding_grazing/Estimator",showWarnings = F)
dir.create("../Data/Step1_Understanding_grazing/Keep_models",showWarnings = F)
dir.create("../Data/Step1_Understanding_grazing/Keep_data",showWarnings = F)


Run_model_importance=function(id){
  
  list_mod=expand.grid(with_cover=c(T,F),
                       Stats=c("perim_area_scaling","fmax_psd","Cond_H","PL_expo","Spectral_ratio",
                               "core_area_land","division","fractal_dim","contig","core_area",
                               "flow_length","PLR","KS_dist",
                               "Struct1","Struct2")) #We also add the case where the stats correspond to the first two components of the PCA
  
  d_all=d_all2=tibble()#for keeping all information
  
  d_data=read.table("../Data/Spatial_structure_grazing.csv",sep=";")%>%
    Closer_to_normality(.)
  
  d_data[,c(1,9:29,36,38:42,45:ncol(d_data))] = 
    apply(d_data[,c(1,9:29,36,38:42,45:ncol(d_data))],2,z_tranform)
  
  d_data=Perform_PCA_spatial_struc(d_data)
  
  stat=list_mod$Stats[id]
  
  d_data_mod=d_data%>%melt(., measure.vars=stat)%>%
    filter(., !is.na(value))
  
  with_cover=list_mod$with_cover[id]
  
  if (with_cover){
    
    formula_mod=formula(formula_(paste("value ~ Long_cos + Long_sin + Lattitude + Slope + Elevation 
      + Clim1 + Clim2 + Clim3 + Clim4 + Grazing 
      + Woody + Woody*Grazing
      + Clim1*Grazing + Clim2*Grazing + Clim3*Grazing + Clim4*Grazing  
      + rho_p + rho_p*Grazing + Type_veg + Type_veg*Grazing
      + Sand + Sand*Grazing + Org_C + Org_C*Grazing
      + (1|Site_ID)")))
    
  } else{
    
    formula_mod=formula(formula_(paste("value ~ Long_cos + Long_sin + Lattitude + Slope + Elevation 
      + Clim1 + Clim2 + Clim3 + Clim4 + Grazing 
      + Woody + Woody*Grazing
      + Clim1*Grazing + Clim2*Grazing + Clim3*Grazing + Clim4*Grazing  
      + Type_veg + Type_veg*Grazing
      + Sand + Sand*Grazing + Org_C + Org_C*Grazing
      + (1|Site_ID)")))
  }
  
  model_spa_stat  = lmer(formula_mod, d_data_mod,
                         na.action = na.fail ,REML ="FALSE")

  #we remove potential outliers
  rm.outliers = romr.fnc(model_spa_stat, d_data_mod, trim=2.5)
  d_data_out = rm.outliers$data
  
  #saving data
  write.table(d_data_out,paste0("../Data/Step1_Understanding_grazing/Keep_data/Data_",stat,"_",with_cover,".csv"))
  
  model_spa_stat = lmer(formula_mod, data = d_data_out, 
                        na.action = na.fail,REML ="FALSE")
  
  # #do some model selection
  select_model=dredge(model_spa_stat, subset = ~ Slope & Elevation &
                        Long_cos & Long_sin & Lattitude &
                        dc(Woody & Grazing, Woody : Grazing) &
                        dc(Clim1 & Grazing, Clim1 : Grazing) &
                        dc(Clim2 & Grazing, Clim2 : Grazing) &
                        dc(Clim3 & Grazing, Clim3 : Grazing) &
                        dc(Clim4 & Grazing, Clim4 : Grazing) &
                        dc(rho_p & Grazing, rho_p : Grazing) &
                        dc(Org_C & Grazing, Org_C : Grazing) &
                        dc(Sand  & Grazing, Sand  : Grazing) &
                        dc(Type_veg & Grazing, Type_veg : Grazing),
                      extra=c(R2m=function(x) r.squaredGLMM(x)[1],
                              R2c=function(x) r.squaredGLMM(x)[2]),
                      options(na.action = "na.fail") )
  
  #extract the result of model selection
  result_select=model.avg(select_model, subset = delta < 2)
  
  #Get the importance of each metric
  importance_mod=sw(result_select)
  
  saveRDS(model_spa_stat,paste0("../Data/Step1_Understanding_grazing/Keep_models/Mod_",stat,"_",with_cover,".rds"))
  
  R2=select_model%>%filter(., AICc<min(AICc)+2)
  
  d_all=rbind(d_all,cbind(tibble(Sp_stat=stat,
                                 N_outliers=rm.outliers$n.removed,
                                 R2m=mean(R2$R2m),R2C=mean(R2$R2c)),
                          Aggregate_importance(importance_mod))%>%
                add_column(., With_cover=with_cover))
  
  
  
  summary_coef=confint(result_select)
  
  #Merge in a df
  d_all2=tibble(Median=confint(result_select,level = 1e-9)[-1,1],
                q1=summary_coef[-1,1], 
                q3=summary_coef[-1,2],
                pvalue=summary(result_select)$coefmat.full[-1,5],
                term=rownames(summary_coef)[-1],
                Stat=stat,
                R2m=mean(R2$R2m),
                R2c=mean(R2$R2c))%>%
    add_column(., With_cover=with_cover)
  
  

  #Checking multicolinearity using VIF 
  
  vif_model=vif(model_spa_stat)
  write.table(vif_model,paste0("../Data/Step1_Understanding_grazing/VIF/VIF_",stat,"_cover_",with_cover,".csv"),sep=";")
  
  
  #Checking spatial auto-correlation at different spatial scales (10, 30, 50)
  
  d_auto_corr=tibble()
  for (spatial_scales in c(10,20,50)){
    
    #Aggregating at different scales
    spatial_coords = cbind(X=as.numeric(d_data_out$Longitude), 
                           Y=as.numeric(d_data_out$Lattitude))
    
    KNN = knearneigh(spatial_coords, k=spatial_scales) 
    nb_KNN=knn2nb(KNN)
    
    #getting residuals
    resids_mod=residuals(model_spa_stat)
    
    #Testing for auto-correlation
    Moran_test=moran.test(resids_mod, nb2listw(nb_KNN))
    
    d_auto_corr=rbind(d_auto_corr,tibble(K=spatial_scales,
                                         Moran_stat=Moran_test$statistic,
                                         p_val=Moran_test$p.value))
    
  }
  
  # Saving autocorrelation checks
  write.table(d_auto_corr,paste0("../Data/Step1_Understanding_grazing/Auto_corr/Auto_corr_",stat,"_cover_",with_cover,".csv"),sep=";")
  
  #Saving each DF
  write.table(d_all,paste0("../Data/Step1_Understanding_grazing/Importance/Importance_",stat,"_",with_cover,".csv"),sep=";")
  write.table(d_all2,paste0("../Data/Step1_Understanding_grazing/Estimator/Estimators_model_",stat,"_",with_cover,".csv"),sep=";")
  
}

library(parallel)

mclapply(1:30,Run_model_importance,mc.cores = 30)


d_all=d_all2=tibble()
for (k in list.files("../Data/Step1_Understanding_grazing/Importance","E.csv")){
  d_all=rbind(d_all,read.table(paste0("../Data/Step1_Understanding_grazing/Importance/",k),sep=";"))
}
for (k in list.files("../Data/Step1_Understanding_grazing/Estimator","E.csv")){
  d_all2=rbind(d_all2,read.table(paste0("../Data/Step1_Understanding_grazing/Estimator/",k),sep=";"))
}
write.table(d_all,"../Data/Step1_Understanding_grazing/Importance.csv",sep=";")
write.table(d_all2,"../Data/Step1_Understanding_grazing/Estimators_model.csv",sep=";")



Run_model_importance_no_inter=function(id){
  
  list_mod=expand.grid(with_cover=c(F,T),
                       Stats=c("perim_area_scaling","fmax_psd","Cond_H","PL_expo","Spectral_ratio",
                               "core_area_land","division","fractal_dim","contig","core_area",
                               "flow_length","PLR",
                               "Struct1","Struct2")) #We also add the case where the stats correspond to the first two components of the PCA
  
  d_all=d_all2=tibble()#for keeping all information
  
  d_data=read.table("../Data/Spatial_structure_grazing.csv",sep=";")%>%
    Closer_to_normality(.)
  
  d_data[,c(1,9:29,36,38:42,45:ncol(d_data))] = 
    apply(d_data[,c(1,9:29,36,38:42,45:ncol(d_data))],2,z_tranform)
  
  d_data=Perform_PCA_spatial_struc(d_data)
  
  stat=list_mod$Stats[id]
  
  d_data_mod=d_data%>%melt(., measure.vars=stat)%>%
    filter(., !is.na(value))
  
  with_cover=list_mod$with_cover[id]
  
  if (with_cover){
    
    formula_mod=formula(formula_(paste("value ~ Long_cos + Long_sin + Lattitude + Slope + Elevation 
      + Clim1 + Clim2 + Clim3 + Clim4 + Grazing 
      + Woody + Woody*Grazing
      + Type_veg + rho_p       + Sand + Org_C       + (1|Site_ID)")))
    
  } else{
    
    formula_mod=formula(formula_(paste("value ~ Long_cos + Long_sin + Lattitude + Slope + Elevation 
      + Clim1 + Clim2 + Clim3 + Clim4 + Grazing 
      + Woody + Woody*Grazing
      + Type_veg  + Sand + Org_C
      + (1|Site_ID)")))
  }
  
  model_spa_stat  = lmer(formula_mod, d_data_mod,
                         na.action = na.fail ,REML ="FALSE")
  
  #we remove potential outliers
  rm.outliers = romr.fnc(model_spa_stat, d_data_mod, trim=2.5)
  d_data_out = rm.outliers$data
  
  
  
  model_spa_stat = lmer(formula_mod, data = d_data_out, 
                        na.action = na.fail,REML ="FALSE")
  
  # #do some model selection
  select_model=dredge(model_spa_stat, subset = ~ Slope & Elevation &
                        Long_cos & Long_sin & Lattitude &
                        dc(Woody & Grazing, Woody : Grazing) &
                        dc(Clim1 & Grazing, Clim1 : Grazing) &
                        dc(Clim2 & Grazing, Clim2 : Grazing) &
                        dc(Clim3 & Grazing, Clim3 : Grazing) &
                        dc(Clim4 & Grazing, Clim4 : Grazing) &
                        dc(rho_p & Grazing, rho_p : Grazing) &
                        dc(Org_C & Grazing, Org_C : Grazing) &
                        dc(Sand  & Grazing, Sand  : Grazing) &
                        dc(Type_veg & Grazing, Type_veg : Grazing),
                      extra=c(R2m=function(x) r.squaredGLMM(x)[1],
                              R2c=function(x) r.squaredGLMM(x)[2]),
                      options(na.action = "na.fail") )
  
  
  
  #extract the result of model selection
  result_select=summary(model.avg(select_model, subset = delta < 2))
  
  
  #Get the importance of each metric
  importance_mod=sw(result_select)
  R2=select_model%>%filter(., AICc<min(AICc)+2)
  
  d_all=rbind(d_all,cbind(tibble(Sp_stat=stat,
                                 N_outliers=rm.outliers$n.removed,
                                 R2m=mean(R2$R2m),R2C=mean(R2$R2c)),
                          Aggregate_importance(importance_mod))%>%
                add_column(., With_cover=with_cover))
  
  write.table(d_all,paste0("../Data/Step1_Understanding_grazing/Importance/Importance_",stat,"_",with_cover,"nointer.csv"),sep=";")
  
}

library(parallel)

mclapply(1:30,Run_model_importance_no_inter,mc.cores = 30)


d_all=d_all2=tibble()
for (k in list.files("../Data/Step1_Understanding_grazing/Importance","nointer")){
  d_all=rbind(d_all,read.table(paste0("../Data/Step1_Understanding_grazing/Importance/",k),sep=";"))
}
for (k in list.files("../Data/Step1_Understanding_grazing/Estimator","nointer")){
  d_all2=rbind(d_all2,read.table(paste0("../Data/Step1_Understanding_grazing/Estimator/",k),sep=";"))
}
write.table(d_all,"../Data/Step1_Understanding_grazing/Importance_no_inter.csv",sep=";")
write.table(d_all2,"../Data/Step1_Understanding_grazing/Estimators_model_no_inter.csv",sep=";")


### Running the models with only aridity and grazing ----

Run_model_importance_aridity=function(id){
  
  list_mod=expand.grid(with_cover=c(T,F),
                       Stats=c("perim_area_scaling","fmax_psd","Cond_H","PL_expo","Spectral_ratio",
                               "core_area_land","division","fractal_dim","contig","core_area",
                               "flow_length","PLR","KS_dist",
                               "Struct1","Struct2")) #We also add the case where the stats correspond to the first two components of the PCA
  
  d_all=d_all2=tibble()#for keeping all information
  
  d_data=read.table("../Data/Spatial_structure_grazing.csv",sep=";")%>%
    Closer_to_normality(.)
  
  d_data[,c(1,9:29,36,38:42,45:ncol(d_data))] = 
    apply(d_data[,c(1,9:29,36,38:42,45:ncol(d_data))],2,z_tranform)
  
  d_data=Perform_PCA_spatial_struc(d_data)
  
  stat=list_mod$Stats[id]
  
  d_data_mod=d_data%>%melt(., measure.vars=stat)%>%
    filter(., !is.na(value))
  
  with_cover=list_mod$with_cover[id]
  
  if (with_cover){
    
    formula_mod=formula(formula_(paste("value ~ Long_cos + Long_sin + Lattitude + Slope + Elevation 
      + Aridity + Grazing 
      + Woody + Woody*Grazing
      + Aridity*Grazing 
      + rho_p + rho_p*Grazing + Type_veg + Type_veg*Grazing
      + Sand + Sand*Grazing + Org_C + Org_C*Grazing
      + (1|Site_ID)")))
    
  } else{
    
    formula_mod=formula(formula_(paste("value ~ Long_cos + Long_sin + Lattitude + Slope + Elevation 
      + Aridity + Grazing 
      + Woody + Woody*Grazing
      + Aridity*Grazing 
      + Type_veg + Type_veg*Grazing
      + Sand + Sand*Grazing + Org_C + Org_C*Grazing
      + (1|Site_ID)")))
  }
  
  model_spa_stat  = lmer(formula_mod, d_data_mod,
                         na.action = na.fail ,REML ="FALSE")
  
  #we remove potential outliers
  rm.outliers = romr.fnc(model_spa_stat, d_data_mod, trim=2.5)
  d_data_out = rm.outliers$data
  
  #saving data
  write.table(d_data_out,paste0("../Data/Step1_Understanding_grazing/Keep_data/Data_",stat,"_",with_cover,"_aridity.csv"))
  
  model_spa_stat = lmer(formula_mod, data = d_data_out, 
                        na.action = na.fail,REML ="FALSE")
  
  # #do some model selection
  select_model=dredge(model_spa_stat, subset = ~ Slope & Elevation &
                        Long_cos & Long_sin & Lattitude &
                        dc(Woody & Grazing, Woody : Grazing) &
                        dc(Aridity & Grazing, Aridity : Grazing) &
                        dc(rho_p & Grazing, rho_p : Grazing) &
                        dc(Org_C & Grazing, Org_C : Grazing) &
                        dc(Sand  & Grazing, Sand  : Grazing) &
                        dc(Type_veg & Grazing, Type_veg : Grazing),
                      extra=c(R2m=function(x) r.squaredGLMM(x)[1],
                              R2c=function(x) r.squaredGLMM(x)[2]),
                      options(na.action = "na.fail") )
  
  #extract the result of model selection
  result_select=model.avg(select_model, subset = delta < 2)
  
  #Get the importance of each metric
  importance_mod=sw(result_select)
  
  saveRDS(model_spa_stat,paste0("../Data/Step1_Understanding_grazing/Keep_models/Mod_",stat,"_",with_cover,"_aridity.rds"))
  
  R2=select_model%>%filter(., AICc<min(AICc)+2)
  
  d_all=rbind(d_all,cbind(tibble(Sp_stat=stat,
                                 N_outliers=rm.outliers$n.removed,
                                 R2m=mean(R2$R2m),R2C=mean(R2$R2c)),
                          Aggregate_importance(importance_mod,T))%>%
                add_column(., With_cover=with_cover))
  
  
  
  summary_coef=confint(result_select)
  
  #Merge in a df
  d_all2=tibble(Median=confint(result_select,level = 1e-9)[-1,1],
                q1=summary_coef[-1,1], 
                q3=summary_coef[-1,2],
                pvalue=summary(result_select)$coefmat.full[-1,5],
                term=rownames(summary_coef)[-1],
                Stat=stat,
                R2m=mean(R2$R2m),
                R2c=mean(R2$R2c))%>%
    add_column(., With_cover=with_cover)
  
  
  
  #Checking multicolinearity using VIF 
  
  vif_model=vif(model_spa_stat)
  write.table(vif_model,paste0("../Data/Step1_Understanding_grazing/VIF/VIF_",stat,"_cover_",with_cover,"_aridity.csv"),sep=";")
  
  
  #Checking spatial auto-correlation at different spatial scales (10, 30, 50)
  
  d_auto_corr=tibble()
  for (spatial_scales in c(10,20,50)){
    
    #Aggregating at different scales
    spatial_coords = cbind(X=as.numeric(d_data_out$Longitude), 
                           Y=as.numeric(d_data_out$Lattitude))
    
    KNN = knearneigh(spatial_coords, k=spatial_scales) 
    nb_KNN=knn2nb(KNN)
    
    #getting residuals
    resids_mod=residuals(model_spa_stat)
    
    #Testing for auto-correlation
    Moran_test=moran.test(resids_mod, nb2listw(nb_KNN))
    
    d_auto_corr=rbind(d_auto_corr,tibble(K=spatial_scales,
                                         Moran_stat=Moran_test$statistic,
                                         p_val=Moran_test$p.value))
    
  }
  
  # Saving autocorrelation checks
  write.table(d_auto_corr,paste0("../Data/Step1_Understanding_grazing/Auto_corr/Auto_corr_",stat,"_cover_",with_cover,"_aridity.csv"),sep=";")
  
  #Saving each DF
  write.table(d_all,paste0("../Data/Step1_Understanding_grazing/Importance/Importance_",stat,"_",with_cover,"_aridity.csv"),sep=";")
  write.table(d_all2,paste0("../Data/Step1_Understanding_grazing/Estimator/Estimators_model_",stat,"_",with_cover,"_aridity.csv"),sep=";")
  
}

library(parallel)

mclapply(1:30,Run_model_importance_aridity,mc.cores = 30)


d_all=d_all2=tibble()
for (k in list.files("../Data/Step1_Understanding_grazing/Importance","aridity")){
  d_all=rbind(d_all,read.table(paste0("../Data/Step1_Understanding_grazing/Importance/",k),sep=";"))
}
for (k in list.files("../Data/Step1_Understanding_grazing/Estimator","aridity")){
  d_all2=rbind(d_all2,read.table(paste0("../Data/Step1_Understanding_grazing/Estimator/",k),sep=";"))
}
write.table(d_all,"../Data/Step1_Understanding_grazing/Importance_aridity.csv",sep=";")
write.table(d_all2,"../Data/Step1_Understanding_grazing/Estimators_model_aridity.csv",sep=";")


Run_model_importance_aridity_no_inter=function(id){
  
  list_mod=expand.grid(with_cover=c(T,F),
                       Stats=c("perim_area_scaling","fmax_psd","Cond_H","PL_expo","Spectral_ratio",
                               "core_area_land","division","fractal_dim","contig","core_area",
                               "flow_length","PLR","KS_dist",
                               "Struct1","Struct2")) #We also add the case where the stats correspond to the first two components of the PCA
  
  d_all=d_all2=tibble()#for keeping all information
  
  d_data=read.table("../Data/Spatial_structure_grazing.csv",sep=";")%>%
    Closer_to_normality(.)
  
  d_data[,c(1,9:29,36,38:42,45:ncol(d_data))] = 
    apply(d_data[,c(1,9:29,36,38:42,45:ncol(d_data))],2,z_tranform)
  
  d_data=Perform_PCA_spatial_struc(d_data)
  
  stat=list_mod$Stats[id]
  
  d_data_mod=d_data%>%melt(., measure.vars=stat)%>%
    filter(., !is.na(value))
  
  with_cover=list_mod$with_cover[id]
  
  if (with_cover){
    
    formula_mod=formula(formula_(paste("value ~ Long_cos + Long_sin + Lattitude + Slope + Elevation 
      + Aridity + Grazing 
      + Woody + rho_p + Type_veg 
      + Sand + Org_C 
      + (1|Site_ID)")))
    
  } else{
    
    formula_mod=formula(formula_(paste("value ~ Long_cos + Long_sin + Lattitude + Slope + Elevation 
      + Aridity + Grazing 
      + Woody
      + Type_veg
      + Sand + Org_C 
      + (1|Site_ID)")))
  }
  
  model_spa_stat  = lmer(formula_mod, d_data_mod,
                         na.action = na.fail ,REML ="FALSE")
  
  #we remove potential outliers
  rm.outliers = romr.fnc(model_spa_stat, d_data_mod, trim=2.5)
  d_data_out = rm.outliers$data
  
  #saving data
  write.table(d_data_out,paste0("../Data/Step1_Understanding_grazing/Keep_data/Data_",stat,"_",with_cover,"_aridity_no_inter.csv"))
  
  model_spa_stat = lmer(formula_mod, data = d_data_out, 
                        na.action = na.fail,REML ="FALSE")
  
  # #do some model selection
  select_model=dredge(model_spa_stat, subset = ~ Slope & Elevation &
                        Long_cos & Long_sin & Lattitude &
                        dc(Woody & Grazing, Woody : Grazing) &
                        dc(Aridity & Grazing, Aridity : Grazing) &
                        dc(rho_p & Grazing, rho_p : Grazing) &
                        dc(Org_C & Grazing, Org_C : Grazing) &
                        dc(Sand  & Grazing, Sand  : Grazing) &
                        dc(Type_veg & Grazing, Type_veg : Grazing),
                      extra=c(R2m=function(x) r.squaredGLMM(x)[1],
                              R2c=function(x) r.squaredGLMM(x)[2]),
                      options(na.action = "na.fail") )
  
  #extract the result of model selection
  result_select=model.avg(select_model, subset = delta < 2)
  
  #Get the importance of each metric
  importance_mod=sw(result_select)
  
  saveRDS(model_spa_stat,paste0("../Data/Step1_Understanding_grazing/Keep_models/Mod_",stat,"_",with_cover,"_aridity_no_inter.rds"))
  
  R2=select_model%>%filter(., AICc<min(AICc)+2)
  
  d_all=rbind(d_all,cbind(tibble(Sp_stat=stat,
                                 N_outliers=rm.outliers$n.removed,
                                 R2m=mean(R2$R2m),R2C=mean(R2$R2c)),
                          Aggregate_importance(importance_mod,T))%>%
                add_column(., With_cover=with_cover))
  
  
  
  summary_coef=confint(result_select)
  
  #Merge in a df
  d_all2=tibble(Median=confint(result_select,level = 1e-9)[-1,1],
                q1=summary_coef[-1,1], 
                q3=summary_coef[-1,2],
                pvalue=summary(result_select)$coefmat.full[-1,5],
                term=rownames(summary_coef)[-1],
                Stat=stat,
                R2m=mean(R2$R2m),
                R2c=mean(R2$R2c))%>%
    add_column(., With_cover=with_cover)
  
  
  
  #Checking multicolinearity using VIF 
  
  vif_model=vif(model_spa_stat)
  write.table(vif_model,paste0("../Data/Step1_Understanding_grazing/VIF/VIF_",stat,"_cover_",with_cover,"_aridity_no_inter.csv"),sep=";")
  
  
  #Checking spatial auto-correlation at different spatial scales (10, 30, 50)
  
  d_auto_corr=tibble()
  for (spatial_scales in c(10,20,50)){
    
    #Aggregating at different scales
    spatial_coords = cbind(X=as.numeric(d_data_out$Longitude), 
                           Y=as.numeric(d_data_out$Lattitude))
    
    KNN = knearneigh(spatial_coords, k=spatial_scales) 
    nb_KNN=knn2nb(KNN)
    
    #getting residuals
    resids_mod=residuals(model_spa_stat)
    
    #Testing for auto-correlation
    Moran_test=moran.test(resids_mod, nb2listw(nb_KNN))
    
    d_auto_corr=rbind(d_auto_corr,tibble(K=spatial_scales,
                                         Moran_stat=Moran_test$statistic,
                                         p_val=Moran_test$p.value))
    
  }
  
  # Saving autocorrelation checks
  write.table(d_auto_corr,paste0("../Data/Step1_Understanding_grazing/Auto_corr/Auto_corr_",stat,"_cover_",with_cover,"_aridity_no_inter.csv"),sep=";")
  
  #Saving each DF
  write.table(d_all,paste0("../Data/Step1_Understanding_grazing/Importance/Importance_",stat,"_",with_cover,"_aridity_no_inter.csv"),sep=";")
  write.table(d_all2,paste0("../Data/Step1_Understanding_grazing/Estimator/Estimators_model_",stat,"_",with_cover,"_aridity_no_inter.csv"),sep=";")
  
}

library(parallel)

mclapply(1:30,Run_model_importance_aridity_no_inter,mc.cores = 30)


d_all=d_all2=tibble()
for (k in list.files("../Data/Step1_Understanding_grazing/Importance","aridity_no_inter")){
  d_all=rbind(d_all,read.table(paste0("../Data/Step1_Understanding_grazing/Importance/",k),sep=";"))
}
for (k in list.files("../Data/Step1_Understanding_grazing/Estimator","aridity_no_inter")){
  d_all2=rbind(d_all2,read.table(paste0("../Data/Step1_Understanding_grazing/Estimator/",k),sep=";"))
}
write.table(d_all,"../Data/Step1_Understanding_grazing/Importance_aridity_no_inter.csv",sep=";")
write.table(d_all2,"../Data/Step1_Understanding_grazing/Estimators_model_aridity_no_inter.csv",sep=";")

## >> 2) Analyzing the residuals ----

### 1) Model with all climatic variables ----

d_data=read.table("../Data/Spatial_structure_grazing.csv",sep=";")%>%
  Closer_to_normality(.)
d_data[,c(1,9:29,36,38:42,45:ncol(d_data))] = 
  apply(d_data[,c(1,9:29,36,38:42,45:ncol(d_data))],2,z_tranform)
d_data=Perform_PCA_spatial_struc(d_data)
save=d_data

d_slope=tibble()
for (k in c("perim_area_scaling","fmax_psd","Cond_H","PL_expo","Spectral_ratio",
            "core_area_land","division","contig","core_area",
            "flow_length","PLR","KS_dist",
            "Struct1","Struct2")){
  
  #for each we plot the slope against the partial residuals with and without cover
  
  
  d_data_out=read.table(paste0("../Data/Step1_Understanding_grazing/Keep_data/Data_",k,"_TRUE.csv"),sep=" ")
  model_spa_stat=readRDS(paste0("../Data/Step1_Understanding_grazing/Keep_models/Mod_",k,"_TRUE.rds"))
  
  resid_mod=visreg::visreg(fit = model_spa_stat,xvar="Grazing",plot=F) 
  
  mod_cov=lm(visregRes~Grazing,resid_mod$res)
  
  d_slope=rbind(d_slope,
                tibble(pval=summary(mod_cov)$coefficient[2,4],
                       slope=summary(mod_cov)$coefficient[2,1],
                       Low_int=confint(mod_cov)[2,1],
                       High_int=confint(mod_cov)[2,2],
                       With_cover=T,
                       Stat=k
                ))
  
  d_data_out=read.table(paste0("../Data/Step1_Understanding_grazing/Keep_data/Data_",k,"_FALSE.csv"),sep=" ")
  model_spa_stat=readRDS(paste0("../Data/Step1_Understanding_grazing/Keep_models/Mod_",k,"_FALSE.rds"))
  
  resid_mod=visreg::visreg(fit = model_spa_stat,xvar="Grazing",plot=F) 
  
  mod_no_cov=lm(visregRes~Grazing,resid_mod$res)
  
  d_slope=rbind(d_slope,
                tibble(pval=summary(mod_no_cov)$coefficient[2,4],
                       slope=summary(mod_no_cov)$coefficient[2,1],
                       Low_int=confint(mod_no_cov)[2,1],
                       High_int=confint(mod_no_cov)[2,2],
                       With_cover=F,
                       Stat=k
                ))
}

write.table(d_slope,"../Data/Step1_Understanding_grazing/Slope_partial_residuals.csv",sep=";")

d_data=read.table("../Data/Spatial_structure_grazing.csv",sep=";")%>%
  Closer_to_normality(.)
d_data[,c(1,9:29,36,38:42,45:ncol(d_data))] = 
  apply(d_data[,c(1,9:29,36,38:42,45:ncol(d_data))],2,z_tranform)
d_data=Perform_PCA_spatial_struc(d_data)
save=d_data

d_slope=tibble()
for (k in c("perim_area_scaling","fmax_psd","Cond_H","PL_expo","Spectral_ratio",
            "core_area_land","division","contig","core_area",
            "flow_length","PLR","KS_dist",
            "Struct1","Struct2")){
  
  #for each we plot the slope against the partial residuals with and without cover
  
  
  d_data_out=read.table(paste0("../Data/Step1_Understanding_grazing/Keep_data/Data_",k,"_TRUE.csv"),sep=" ")
  model_spa_stat=readRDS(paste0("../Data/Step1_Understanding_grazing/Keep_models/Mod_",k,"_TRUE.rds"))
  
  resid_mod=visreg::visreg(fit = model_spa_stat,xvar="Grazing",plot=F) 
  
  mod_cov=lm(visregRes~Grazing,resid_mod$res)
  
  d_slope=rbind(d_slope,
                tibble(pval=summary(mod_cov)$coefficient[2,4],
                       slope=summary(mod_cov)$coefficient[2,1],
                       Low_int=confint(mod_cov)[2,1],
                       High_int=confint(mod_cov)[2,2],
                       With_cover=T,
                       Stat=k
                ))
  
  d_data_out=read.table(paste0("../Data/Step1_Understanding_grazing/Keep_data/Data_",k,"_FALSE.csv"),sep=" ")
  model_spa_stat=readRDS(paste0("../Data/Step1_Understanding_grazing/Keep_models/Mod_",k,"_FALSE.rds"))
  
  resid_mod=visreg::visreg(fit = model_spa_stat,xvar="Grazing",plot=F) 
  
  mod_no_cov=lm(visregRes~Grazing,resid_mod$res)
  
  d_slope=rbind(d_slope,
                tibble(pval=summary(mod_no_cov)$coefficient[2,4],
                       slope=summary(mod_no_cov)$coefficient[2,1],
                       Low_int=confint(mod_no_cov)[2,1],
                       High_int=confint(mod_no_cov)[2,2],
                       With_cover=F,
                       Stat=k
                ))
}

write.table(d_slope,"../Data/Step1_Understanding_grazing/Slope_partial_residuals.csv",sep=";")



#Residuals without cover in the model


list_title=c("PLR","Power-law exp. of the PSD",
             "log (largest patch)","Bare soil connectivity",
             "Fractal scaling area, perim.","Mean % of core in patches",
             "Contiguity",
             "% landscape covered by core")

list_name_x=c("Grazing intensity","Facilitation",
              'PC1 climatic variables',"PC2 climatic variables",
              "PC3 climatic variables","PC4 climatic variables","Sand","% woody cover")

abrev_metric=c("","F","C1","C2","C3","C4","S","W")


list_labels_plot=c(
  "'Y = ' * alpha[1] * 'G*' * + alpha[2] * 'F'* + alpha[3] *  'F x G'",
  "'Y = ' * alpha[1] * 'G*' * + alpha[2] * 'C2***'* + alpha[3] *  'C2 x G'",
  "'Y = ' * alpha[1] * 'G*' * + alpha[2] * 'C4'* + alpha[3] *  'C4 x G*'",
  "'Y = ' * alpha[1] * 'G**' * + alpha[2] * 'C2***'* + alpha[3] *  'C2 x G**'",
  "'Y = ' * alpha[1] * 'G**' * + alpha[2] * 'C3***'* + alpha[3] *  'C3 x G'",
  "'Y = ' * alpha[1] * 'G*' * + alpha[2] * 'C4***'* + alpha[3] *  'C4 x G***'",
  "'Y = ' * alpha[1] * 'G**' * + alpha[2] * 'W?'* + alpha[3] *  'W x G'",
  "'Y = ' * alpha[1] * 'G' * + alpha[2] * 'C2***'* + alpha[3] *  'C2 x G***'",
  "'Y = ' * alpha[1] * 'G' * + alpha[2] * 'C3***'* + alpha[3] *  'C3 x G?'",
  "'Y = ' * alpha[1] * 'G' * + alpha[2] * 'C4*'* + alpha[3] *  'C4 x G'",
  "'Y = ' * alpha[1] * 'G' * + alpha[2] * 'W**'* + alpha[3] *  'W x G'",
  "'Y = ' * alpha[1] * 'G*' * + alpha[2] * 'F***'* + alpha[3] *  'F x G***'",
  "'Y = ' * alpha[1] * 'G***' * + alpha[2] * 'C1***'* + alpha[3] *  'C1 x G***'",
  "'Y = ' * alpha[1] * 'G***' * + alpha[2] * 'C3***'* + alpha[3] *  'C3 x G'",
  "'Y = ' * alpha[1] * 'G***' * + alpha[2] * 'C4***'* + alpha[3] *  'C3 x G'",
  "'Y = ' * alpha[1] * 'G***' * + alpha[2] * 'C4***'* + alpha[3] *  'C3 x G'",
  "'Y = ' * alpha[1] * 'G***' * + alpha[2] * 'S'* + alpha[3] *  'S x G'",
  "'Y = ' * alpha[1] * 'G' * + alpha[2] * 'F***'* + alpha[3] *  'F x G***'",
  "'Y = ' * alpha[1] * 'G' * + alpha[2] * 'C1'* + alpha[3] *  'C1 x G***'",
  "'Y = ' * alpha[1] * 'G' * + alpha[2] * 'C2***'* + alpha[3] *  'C1 x G**'",
  "'Y = ' * alpha[1] * 'G' * + alpha[2] * 'F***'* + alpha[3] *  'F x G***'",
  "'Y = ' * alpha[1] * 'G' * + alpha[2] * 'C2***'* + alpha[3] *  'C2 x G'",
  "'Y = ' * alpha[1] * 'G' * + alpha[2] * 'F***'* + alpha[3] *  'F x G***'",
  "'Y = ' * alpha[1] * 'G***' * + alpha[2] * 'C1***'* + alpha[3] *  'C1 x G***'",
  "'Y = ' * alpha[1] * 'G***' * + alpha[2] * 'C2***'* + alpha[3] *  'C2 x G'",
  "'Y = ' * alpha[1] * 'G***' * + alpha[2] * 'C4***'* + alpha[3] *  'C4 x G'" )




ID_label=ID_title=1

for (stat in c("PLR","PL_expo","fmax_psd","flow_length",
               "perim_area_scaling","core_area","contig",
               "core_area_land",
               "KS_dist")){
  
  title=list_title[ID_title]  
  
  ID_plot=1
  list_plot=list()
  ID_x=1
  
  for (metric in c("Type_veg","Org_C","Clim1","Clim2","Clim3","Clim4","Sand","Woody")){
    
    name_x=list_name_x[ID_x]
    
    d_data_out=  read.table(paste0("../Data/Step1_Understanding_grazing/Keep_data/Data_",stat,"_FALSE.csv"),sep=";")
    if (ncol(d_data_out)==1){
      d_data_out=  read.table(paste0("../Data/Step1_Understanding_grazing/Keep_data/Data_",stat,"_FALSE.csv"),sep=" ")
    }
    model_spa_stat=readRDS(paste0("../Data/Step1_Understanding_grazing/Keep_models/Mod_",stat,"_FALSE.rds"))
    
    
    if (metric=="Type_veg" & "Type_vegGrassland" %in% rownames(summary(model_spa_stat)$coefficients)){
      
      resid_mod=visreg::visreg(fit = model_spa_stat,xvar = "Grazing",by="Type_veg",plot=F) #extracting residuals
      
      p_val=tibble(pval=c(summary(lm(visregRes~Grazing,resid_mod$res%>%filter(.,Type_veg=="Grass_Shrub")))$coefficients[2,4],
                          summary(lm(visregRes~Grazing,resid_mod$res%>%filter(.,Type_veg=="Forest")))$coefficients[2,4],
                          summary(lm(visregRes~Grazing,resid_mod$res%>%filter(.,Type_veg=="Grassland")))$coefficients[2,4],
                          summary(lm(visregRes~Grazing,resid_mod$res%>%filter(.,Type_veg=="Shrubland")))$coefficients[2,4]),
                   Habitat=c("Grass_Shrub","Forest","Grassland","Shrubland")
      )
      
      if (all(p_val$pval<.05)){
        label_pval="P slopes < 0.05"
      }else if (all(p_val$pval>.05)){
        label_pval="NS"
      }else {
        label_pval=paste0(    sapply(which(p_val$pval>.05),function(x){return(
          paste0(p_val$Habitat[x]," = NS")
        )}),collapse = ", "
        )
      }
      
      annot_mod=annotate("text", x =1.5 , y = min(resid_mod$res$visregRes),
                         label = label_pval,size=5)
      
      list_plot[[ID_plot]]=ggplot(resid_mod$res%>%
                                    mutate(., Type_veg=as.character(Type_veg)))+
        geom_jitter(aes(Grazing,visregRes),color="black",alpha=.2,width = .1)+
        geom_smooth(aes(Grazing,visregRes,fill=Type_veg,color=Type_veg),method = "lm",level=.95)+
        scale_x_continuous(labels = c("Ungrazed","Low","Medium","High"))+
        annot_mod+
        the_theme+
        scale_color_manual(values=c("#05401E","#D7414F","#41D77C","#6367C3"),
                           labels=c("Forest","Grass & Shrubs","Grassland","Shrubland"))+
        scale_fill_manual(values=c("#05401E","#D7414F","#41D77C","#6367C3"),
                          labels=c("Forest","Grass & Shrubs","Grassland","Shrubland"))+
        labs(x=name_x,y=title,fill="",color="")+
        theme(title=element_text(size=13),axis.title = element_text(size=15))
      
      
      ID_plot=ID_plot+1
      
      
    }else if (metric  %in% rownames(summary(model_spa_stat)$coefficients)){
      
      resid_mod=visreg::visreg(fit = model_spa_stat,xvar = metric,by="Grazing",plot=F) #extracting residuals
      resid_mod$res$Grazing=rep(0:3,as.vector(table(d_data_out$Grazing))) #correcting vigreg output
      
      p_val=sapply(summary(lm(formula(paste("visregRes~Grazing+",metric,"+",metric,"*Grazing")),resid_mod$res))$coefficients[-1,4],
                   function(x){return(is_signif(x))})
      
      
      annot_mod=annotate("text", x = (max(resid_mod$res[,metric])+min(resid_mod$res[,metric]))/2,
                         y = min(resid_mod$res$visregRes), parse = TRUE,
                         label = list_labels_plot[ID_label],size=5)
      
      ID_label=ID_label+1
      
      
      
      list_plot[[ID_plot]]=ggplot(resid_mod$res%>%
                                    mutate(., Grazing=as.character(Grazing))%>%
                                    dplyr::select(., -value)%>%
                                    melt(., measure.vars=metric)%>%arrange(., Grazing,decreasing=T))+
        geom_point(aes(value,visregRes),color="gray30",alpha=.2)+
        geom_smooth(aes(value,visregRes,fill=Grazing,color=Grazing),method = "lm",level=.95)+
        annot_mod+
        the_theme+
        scale_color_manual(values=c("0"="#FFC170",'1'="#43FFF4","2"="#5196F9","3"="#012E71"),
                           labels=c("Ungrazed","Low","Medium","High"))+
        scale_fill_manual(values=c("0"="#FFC170",'1'="#43FFF4","2"="#5196F9","3"="#012E71"),
                          labels=c("Ungrazed","Low","Medium","High"))+
        labs(x=name_x,y=title,fill="",color="")+
        theme(title=element_text(size=13),axis.title = element_text(size=15))
      
      ID_plot=ID_plot+1
      
    }
    
    
    ID_x=ID_x+1
  }
  
  
  if (stat=="perim_area_scaling"){
    p=ggarrange(plotlist=list_plot,ncol = length(list_plot),nrow = 1,legend = "bottom",labels = c("i","ii","iii","iv","v")[1:length(list_plot)])
  }else{
    p=ggarrange(plotlist=list_plot,ncol = length(list_plot),nrow = 1,legend = "none",labels = c("i","ii","iii","iv","v")[1:length(list_plot)])
  }
  
  
  ggsave(paste0("../Figures/Step1_Understanding_grazing/",
                "Interaction_no_cover/All_interactions_",stat,".pdf"),p,width = 4*length(list_plot),height = 4)
  
  ID_title=ID_title+1
  
}








### 2) Model with aridity and grazing ----

d_data=read.table("../Data/Spatial_structure_grazing.csv",sep=";")%>%
  Closer_to_normality(.)
d_data[,c(1,9:29,36,38:42,45:ncol(d_data))] = 
  apply(d_data[,c(1,9:29,36,38:42,45:ncol(d_data))],2,z_tranform)
d_data=Perform_PCA_spatial_struc(d_data)
save=d_data

d_slope=tibble()
for (k in c("perim_area_scaling","fmax_psd","Cond_H","PL_expo","Spectral_ratio",
            "core_area_land","division","contig","core_area",
            "flow_length","PLR","KS_dist",
            "Struct1","Struct2")){
  
  #for each we plot the slope against the partial residuals with and without cover
  
  
  d_data_out=read.table(paste0("../Data/Step1_Understanding_grazing/Keep_data/Data_",k,"_TRUE_aridity.csv"),sep=" ")
  model_spa_stat=readRDS(paste0("../Data/Step1_Understanding_grazing/Keep_models/Mod_",k,"_TRUE_aridity.rds"))
  
  resid_mod=visreg::visreg(fit = model_spa_stat,xvar="Grazing",plot=F) 
  
  mod_cov=lm(visregRes~Grazing,resid_mod$res)
  
  d_slope=rbind(d_slope,
                tibble(pval=summary(mod_cov)$coefficient[2,4],
                       slope=summary(mod_cov)$coefficient[2,1],
                       Low_int=confint(mod_cov)[2,1],
                       High_int=confint(mod_cov)[2,2],
                       With_cover=T,
                       Stat=k,Driver="Grazing"
                ))
  
  resid_mod=visreg::visreg(fit = model_spa_stat,xvar="Aridity",plot=F) 
  
  mod_cov=lm(visregRes~Aridity,resid_mod$res)
  
  d_slope=rbind(d_slope,
                tibble(pval=summary(mod_cov)$coefficient[2,4],
                       slope=summary(mod_cov)$coefficient[2,1],
                       Low_int=confint(mod_cov)[2,1],
                       High_int=confint(mod_cov)[2,2],
                       With_cover=T,
                       Stat=k,Driver="Aridity"
                ))
  
  
  d_data_out=read.table(paste0("../Data/Step1_Understanding_grazing/Keep_data/Data_",k,"_FALSE_aridity.csv"),sep=" ")
  model_spa_stat=readRDS(paste0("../Data/Step1_Understanding_grazing/Keep_models/Mod_",k,"_FALSE_aridity.rds"))
  
  resid_mod=visreg::visreg(fit = model_spa_stat,xvar="Grazing",plot=F) 
  
  mod_cov=lm(visregRes~Grazing,resid_mod$res)
  
  d_slope=rbind(d_slope,
                tibble(pval=summary(mod_cov)$coefficient[2,4],
                       slope=summary(mod_cov)$coefficient[2,1],
                       Low_int=confint(mod_cov)[2,1],
                       High_int=confint(mod_cov)[2,2],
                       With_cover=F,
                       Stat=k,Driver="Grazing"
                ))
  
  resid_mod=visreg::visreg(fit = model_spa_stat,xvar="Aridity",plot=F) 
  
  mod_cov=lm(visregRes~Aridity,resid_mod$res)
  
  d_slope=rbind(d_slope,
                tibble(pval=summary(mod_cov)$coefficient[2,4],
                       slope=summary(mod_cov)$coefficient[2,1],
                       Low_int=confint(mod_cov)[2,1],
                       High_int=confint(mod_cov)[2,2],
                       With_cover=F,
                       Stat=k,Driver="Aridity"
                ))
}

write.table(d_slope,"../Data/Step1_Understanding_grazing/Slope_partial_residuals_aridity.csv",sep=";")

# ---------------------- Step 3: Indirect effect of grazing  ----

d_data=read.table("../Data/Spatial_structure_grazing.csv",sep=";")%>%
  Closer_to_normality(.)

d_data[,c(1,9:29,36,38:42,45:ncol(d_data))] = apply(d_data[,c(1,9:29,36,38:42,45:ncol(d_data))],2,z_tranform)

#Getting the multifunctionality index
d_data$MF=z_tranform(rowSums(d_data[,56:(ncol(d_data))],na.rm = T))

#Adding the first two components from the PCA
d_data=Perform_PCA_spatial_struc(d_data)

dir.create("../Figures/Step1_Understanding_grazing/SEM/",showWarnings = F)

param_list=expand.grid(Stat=c("Struct1","Struct2",
                              "perim_area_scaling","PL_expo","fmax_psd",
                              "PLR","flow_length","core_area_land","core_area",
                              "contig","KS_dist"),
                       MF=c(T,F)) #MF instead of organic carbon

d_indirect=tibble()

for (ID in 1:nrow(param_list)){
  
  dir.create(paste0("../Figures/Step1_Understanding_grazing/SEM/",param_list$Stat[ID]),showWarnings = F)
  
  k=as.character(param_list$Stat[ID])
  
  model_lmer=lm("value ~ Long_cos + Long_sin + Lattitude + Slope + Elevation +
                        Clim1 + Clim2 + Clim3 + Clim4 + Sand + Type_veg",
                data = d_data%>%melt(., measure.vars=k)%>%filter(., !is.na(value)),
                na.action = na.omit)
  
  #we remove potential outliers
  rm.outliers = romr.fnc(model_lmer, d_data%>%melt(., measure.vars=k)%>%filter(., !is.na(value)),
                         trim=2.5)
  d_data_out = rm.outliers$data
  
  #We first control for all covariates and extract the residuals
  model_lmer=lm("value ~ Long_cos + Long_sin + Lattitude + Slope + Elevation +
                        Clim1 + Clim2 + Clim3 + Clim4 + Sand + Type_veg",
                data = d_data_out,
                na.action = na.fail)
  
  resid_model=residuals(model_lmer) #extract residuals
  
  save=d_data[as.numeric(names(resid_model)),]%>% #add it to the dataframe 
    add_column(., Resid_mod=resid_model)
  
  if (param_list$MF[ID]){ save$Org_C=save$MF}
  
  #DOING the SEMs
  
  #SEM with all grazing intensity
  
  d_sem=save
  all_d=summary(psem(
    lmer(Resid_mod ~ (1|Site_ID) + Grazing + Woody + rho_p + Org_C, d_sem),
    lmer(Woody ~ (1|Site_ID) + Grazing + Org_C , d_sem),
    lmer(rho_p ~ (1|Site_ID) + Grazing + Org_C   , d_sem),
    lmer(Org_C ~ (1|Site_ID) + Grazing, d_sem)
  ))
  
  #SEM with low grazing intensity
  
  d_sem=filter(save,Grazing %in% 0:1)
  low_graz=summary(psem(
    lmer(Resid_mod ~ (1|Site_ID) + Grazing + Woody + rho_p + Org_C, d_sem),
    lmer(Woody ~ (1|Site_ID) + Grazing + Org_C , d_sem),
    lmer(rho_p ~ (1|Site_ID) + Grazing + Org_C   , d_sem),
    lmer(Org_C ~ (1|Site_ID) + Grazing, d_sem)
  ))
  
  #SEM with high grazing intensity
  
  d_sem=filter(save,Grazing %in% 2:3)
  high_graz=summary(psem(
    lmer(Resid_mod ~ (1|Site_ID) + Grazing + Woody + rho_p + Org_C, d_sem),
    lmer(Woody ~ (1|Site_ID) + Grazing + Org_C , d_sem),
    lmer(rho_p ~ (1|Site_ID) + Grazing + Org_C  , d_sem),
    lmer(Org_C ~ (1|Site_ID) + Grazing, d_sem)
  ))
  
  Plot_SEM_with_cover(summary_sem = all_d,pdf_ = T,
                      name = paste0("../Figures/Step1_Understanding_grazing/SEM/",param_list$Stat[ID],"/",
                                    "SEM_type_veg_",param_list$Stat[ID],"_",param_list$Cov[ID],"_MF_",param_list$MF[ID],"_all"),
                      title_ = paste0("SEM_",k,"_all"),
                      name_var = k,MF = param_list$MF[ID])
  
  Plot_SEM_with_cover(summary_sem = low_graz,pdf_ = T,
                      name = paste0("../Figures/Step1_Understanding_grazing/SEM/",param_list$Stat[ID],"/",
                                    "SEM_type_veg_",param_list$Stat[ID],"_",param_list$Cov[ID],"_MF_",param_list$MF[ID],"_low"),
                      title_ = paste0("SEM_",k,"_low"),
                      name_var = k,MF = param_list$MF[ID])
  
  Plot_SEM_with_cover(summary_sem = high_graz,pdf_ = T,
                      name = paste0("../Figures/Step1_Understanding_grazing/SEM/",param_list$Stat[ID],"/",
                                    "SEM_type_veg_",param_list$Stat[ID],"_",param_list$Cov[ID],"_MF_",param_list$MF[ID],"_high"),
                      title_ = paste0("SEM_",k,"_high"),
                      name_var = k,MF = param_list$MF[ID])
  
  d_indirect=rbind(d_indirect,Get_indirect_effects_grazing(all_d)%>%add_column(., Stat=k,Grazing="All"))
  d_indirect=rbind(d_indirect,Get_indirect_effects_grazing(low_graz)%>%add_column(., Stat=k,Grazing="Low"))
  d_indirect=rbind(d_indirect,Get_indirect_effects_grazing(high_graz)%>%add_column(., Stat=k,Grazing="High"))
  
}







# ---------------------- Step 4: Effect of spatial resolution  ----



#To evaluate the importance of spatial resolution, we compare the best models selected before
#with and without spatial resolution as a predictor


d_AIC=tibble()
list_mod=expand.grid(with_cover=c(T),
                     Stats=c("perim_area_scaling","fmax_psd","Cond_H","PL_expo","Spectral_ratio",
                             "core_area_land","division","fractal_dim","contig","core_area",
                             "flow_length","PLR",
                             "Struct1","Struct2")) #We also add the case where the stats correspond to the first two components of the PCA

for (id in 1:nrow(list_mod)){  

  stat=list_mod$Stats[id]

  d_data_out=  read.table(paste0("../Data/Step1_Understanding_grazing/Keep_data/Data_",stat,"_TRUE.csv"),sep=";")
  if (ncol(d_data_out)==1){
    d_data_out=  read.table(paste0("../Data/Step1_Understanding_grazing/Keep_data/Data_",stat,"_TRUE.csv"),sep=" ")
  }
  model_spa_stat=readRDS(paste0("../Data/Step1_Understanding_grazing/Keep_models/Mod_",stat,"_TRUE.rds"))
  
  mod_with_res=lmer(formula = formula(
    paste(paste(formula(model_spa_stat))[2],paste(formula(model_spa_stat))[1],paste(paste(formula(model_spa_stat))[3],"+ Resolution")
  )),data = d_data_out,REML = F)
  
  d_AIC=rbind(d_AIC,tibble(AIC_no_res=AICc(mod_with_res,model_spa_stat)[2,2],
                           AIC_res=AICc(mod_with_res,model_spa_stat)[1,2],
                           Stat=stat))
}
d_AIC%>%Filter_relevant_stats(.)%>%Rename_spatial_statistics(.)%>%
  add_column(., Delta_AIC=.$AIC_no_res-.$AIC_res)

# ---------------------- Step 5: Effect on the stability and resilience ----
## >> 1) Inference parameters ----

NA_kept=100;id_plot=1    
d_biodesert=read.table("../Data/Spatial_structure_grazing.csv",sep=";")[,1:11]%>%
  add_column(., ID_sim=1:nrow(.))

d_sim=read.table("../Data/All_sim.csv",sep=";")%>%
  dplyr::relocate(., Pooling,.after =q )%>%
  filter(., PL_expo>0,!is.na(PLR))%>%
  dplyr::select(., -ID)

rownames(d_sim)=1:nrow(d_sim)
`%!in%` = Negate(`%in%`)
n_param=3

d_param_infer_NN=array(0,c(NA_kept,nrow(d_biodesert),3))
d_param_infer_rej=array(0,c(NA_kept,nrow(d_biodesert),3))

d_NRMSE_sumstat=x_y_stat=tibble()


sumstat_kept=1:11

for (empirical_id in which(d_biodesert$rho_p>.05)){
  
  print(empirical_id)
  target=d_biodesert[empirical_id,sumstat_kept]
  matrix_param=d_sim[,1:3]
  
  mat_sumstat=rbind(d_sim[,3+sumstat_kept],target)
  
  #1) Boxcox
  
  for (x in 1:ncol(mat_sumstat)){
    if (any(is.na(target)) & x %!in% which(is.na(target))){
      if (colnames(mat_sumstat)[x] %in% c("skewness","moran_I","fmax_psd")){
        
        
        
        b=boxcox(lm(mat_sumstat[,x]+abs(min(mat_sumstat[,x]))+.5 ~ 1),plotit = F,eps = .05)     #Working with positive values
        lambda_x=b$x[which.max(b$y)]
        if (lambda_x !=0){ #to avoid errors
          mat_sumstat[,x] = (exp(mat_sumstat[,x]*(lambda_x)) -1)/(lambda_x)
        }
        
      }else {
        b=boxcox(lm(mat_sumstat[,x]+.5 ~ 1),plotit = F,eps = .05)
        lambda_x=b$x[which.max(b$y)]
        if (lambda_x !=0){ #to avoid errors
          mat_sumstat[,x] = (mat_sumstat[,x]^(lambda_x) -1)/(lambda_x)
        }
      }
    }
  }
  
  
  #2) Scaling
  
  for (x in 1:ncol(mat_sumstat)) mat_sumstat[,x] = (mat_sumstat[,x]-mean(mat_sumstat[,x],na.rm = T))/sd(mat_sumstat[,x],na.rm = T)
  
  if (any(is.na(mat_sumstat[nrow(mat_sumstat),]))){
    
    which_na=which(is.na(mat_sumstat[nrow(mat_sumstat),]))
    
    cross_valid=abc(target = mat_sumstat[nrow(mat_sumstat),-which_na],
                    param = matrix_param[-nrow(mat_sumstat),],sumstat = mat_sumstat[-nrow(mat_sumstat),-which_na], #removing the target data
                    tol = 1000/nrow(matrix_param),method = "rejection") #we keep the 1000 closest simulations for the first step
    
  }else {
    cross_valid=abc(target = mat_sumstat[nrow(mat_sumstat),],
                    param = matrix_param[-nrow(mat_sumstat),],sumstat = mat_sumstat[-nrow(mat_sumstat),], #removing the target data
                    tol = 1000/nrow(matrix_param),method = "rejection") #we keep the 1000 closest simulations for the first step
  }
  
  
  #Keeping 1000 simulations and doing the same steps again: normality, scaling and PLS
  
  mat_sumstat_step1=d_sim[as.numeric(rownames(cross_valid$ss)),3+sumstat_kept] #we keep information with the true values
  mat_sumstat_step1=rbind(mat_sumstat_step1,target)
  
  #again, first box cox
  which_na=which(is.na(mat_sumstat[nrow(mat_sumstat),]))
  
  for (x in 1:ncol(mat_sumstat_step1)){
    
    if (any(which_na) & x %!in% which_na){
      
      if (colnames(mat_sumstat_step1)[x] %in% c("skewness","moran_I","fmax_psd")){
        
        b=boxcox(lm(mat_sumstat_step1[,x]+abs(min(mat_sumstat_step1[,x]))+.5 ~ 1),plotit = F,eps = .05)     #Working with positive values
        lambda_x=b$x[which.max(b$y)]
        
        if (lambda_x !=0){ #to avoid errors
          mat_sumstat_step1[,x] = (exp(mat_sumstat_step1[,x]*(lambda_x)) -1)/(lambda_x)
        }
        
      }else {
        b=boxcox(lm(mat_sumstat_step1[,x] ~ 1),plotit = F,eps = .05)
        lambda_x=b$x[which.max(b$y)]
        if (lambda_x !=0){ #to avoid errors
          mat_sumstat_step1[,x] = (mat_sumstat_step1[,x]^(lambda_x) -1)/(lambda_x)
        }
      }
    }
  }
  
  #and normalization
  for (x in 1:ncol(mat_sumstat_step1)) {
    if (length(unique(mat_sumstat_step1[,x])) != 1){
      mat_sumstat_step1[,x] = (mat_sumstat_step1[,x]-mean(mat_sumstat_step1[,x],na.rm = T))/sd(mat_sumstat_step1[,x],na.rm = T)
    }
  }
  if (any(is.na(mat_sumstat_step1[nrow(mat_sumstat_step1),]))){
    
    which_na=which(is.na(mat_sumstat_step1[nrow(mat_sumstat_step1),]))
    
    cross_valid=abc(target = mat_sumstat_step1[nrow(mat_sumstat_step1),-which_na],
                    param = cross_valid$unadj.values,
                    sumstat = mat_sumstat_step1[-nrow(mat_sumstat_step1),-which_na], #removing the target data
                    tol = NA_kept/nrow(mat_sumstat_step1),method = "neuralnet",transf = c(rep("logit",2),"none"), #as parameters are proba, we perform logit regression
                    logit.bounds = matrix(c(0,1),3,2,byrow = T),
                    numnet = 10,sizenet = 15)
    
    cross_valid$ss=d_sim[as.numeric(rownames(cross_valid$ss)),3+sumstat_kept] #we keep information with the true values
    
    if(name_plot[id_plot]=="no_PLR"){
      
      x_y_stat=rbind(x_y_stat,target[-which_na]%>%
                       add_column(., PL_expo=NA, Site_ID=empirical_id,Method="rejection",Type="Obs")%>%
                       relocate(., PL_expo,.after =Spectral_ratio ))
      
      x_y_stat=rbind(x_y_stat,as_tibble(t(colMeans(cross_valid$ss)[-which_na]))%>%
                       add_column(.,PL_expo=NA, Site_ID=empirical_id,Method="rejection",Type="Sim")%>%
                       relocate(., PL_expo,.after =Spectral_ratio ))
    }else{
      
      x_y_stat=rbind(x_y_stat,target[-which_na]%>%
                       add_column(., PL_expo=NA, Site_ID=empirical_id,Method="rejection",Type="Obs")%>%
                       relocate(., PL_expo,.after =PLR ))
      
      x_y_stat=rbind(x_y_stat,as_tibble(t(colMeans(cross_valid$ss)[-which_na]))%>%
                       add_column(.,PL_expo=NA, Site_ID=empirical_id,Method="rejection",Type="Sim")%>%
                       relocate(., PL_expo,.after =PLR ))
      
      
    }
    
    
  }else {
    cross_valid=abc(target = mat_sumstat_step1[nrow(mat_sumstat_step1),],
                    param = cross_valid$unadj.values,
                    sumstat = mat_sumstat_step1[-nrow(mat_sumstat_step1),], #removing the target data
                    tol = NA_kept/nrow(mat_sumstat_step1),method = "neuralnet",transf = c(rep("logit",2),"none"), #as parameters are proba, we perform logit regression
                    logit.bounds = matrix(c(0,1),3,2,byrow = T),
                    numnet = 10,sizenet = 15)
    
    cross_valid$ss=d_sim[as.numeric(rownames(cross_valid$ss)),3+sumstat_kept] #we keep information with the true values
    
    x_y_stat=rbind(x_y_stat,target%>%
                     add_column(., Site_ID=empirical_id,Method="rejection",Type="Obs"))
    
    x_y_stat=rbind(x_y_stat,as_tibble(t(colMeans(cross_valid$ss)))%>%
                     add_column(.,Site_ID=empirical_id,Method="rejection",Type="Sim"))
    
  }
  
  
  cross_valid$ss=d_sim[as.numeric(rownames(cross_valid$ss)),3+sumstat_kept] #we keep information with the true values
  
  mat_sumstat=d_sim[,3+sumstat_kept]
  
  if (names(cross_valid)[1]=="unadj.values")names(cross_valid)[1] = "adj.values"
  
  cross_valid$adj.values=cross_valid$adj.values
  
  #NRMSE for the summary statistics observed
  RMSE = sapply(1:ncol(cross_valid$ss),function(x){
    sqrt(sum((cross_valid$ss[,x]-target[,x])**2,na.rm = T)/nrow(cross_valid$ss) )
  }
  )
  
  RMSE_prior=sapply(1:ncol(mat_sumstat),function(x){
    sqrt(sum((mat_sumstat[,x]-target[,x])**2,na.rm = T)/nrow(mat_sumstat) )
  }
  )
  NRMSE = RMSE/RMSE_prior
  
  d_NRMSE_sumstat=rbind(d_NRMSE_sumstat,as_tibble(t(NRMSE)))
  
  d_param_infer_NN[,empirical_id,1]=cross_valid$adj.values[,1] # we keep the whole distribution for p
  d_param_infer_NN[,empirical_id,2]=cross_valid$adj.values[,2] # for q
  d_param_infer_NN[,empirical_id,3]=cross_valid$adj.values[,3] # for the scale of observation
  
  d_param_infer_rej[,empirical_id,1]=cross_valid$unadj.values[,1] # we keep the whole distribution for p
  d_param_infer_rej[,empirical_id,2]=cross_valid$unadj.values[,2] # for q
  d_param_infer_rej[,empirical_id,3]=cross_valid$unadj.values[,3] # for the scale of observation
  
  
}

write.table(d_NRMSE_sumstat,paste0("./NRMSE_sumstat_all_",NA_kept,".csv"),sep=";")
write.table(x_y_stat,paste0("./x_y_stat_all",NA_kept,".csv"),sep=";")
write.table(d_param_infer_rej,paste0("./param_inferred.csv"),sep=";")

## >> 2) Selecting sites ----

library(diptest)
post_param=read.table("../Data/Inferrence/param_inferred.csv",sep=";")

d_biodesert$bimod=sapply(1:nrow(d_biodesert),function(x){
  if (dip.test(post_param[,x])$p.value<.05 | dip.test(post_param[,x+504])$p.value<.05){
    return("bimod")
  }else {
    return("unimod")
  }
})
write.table(which(d_biodesert$bimod!="bimod"),"../Data/Inferrence/Keeping_sites.csv",sep=";")




## >> 3) Computing stability metrics ----

# See ABC_stability.jl file

## >> 4) Post-processing stability metrics ----

d=tibble();step_size=0.005
for (site in list.files("../Data/Inferrence/Prediction/","Dist")){
  
  site_id=as.numeric(gsub(".csv","",strsplit(site,"_")[[1]][3]))
  
  pred=read.table(paste0("../Data/Inferrence/Prediction/",site),sep=",")%>%
    filter(., V1 != 0)
  colnames(pred)=c("p","q","cover")
  
  if (any(pred$cover>.05)){
    index=0;pred$ID_sim=NA
    for (x in 1:nrow(pred)){
      if (pred$p[x]==0.005){
        index=index+1
      }
      pred$ID_sim[x]=index
    }
    p_desert=sapply(unique(pred$ID_sim),function(x){
      d_fil=filter(pred,ID_sim==x)
      if (any(d_fil$cover>0)){
        return(d_fil$p[min(which(d_fil$cover !=0))]-step_size)
      }else {
        return(NA)
      }
    }) 
    
    p_infer=sapply(unique(pred$ID_sim),function(x){
      d_fil=filter(pred,ID_sim==x)
      return(d_fil$p[nrow(d_fil)])
    }) 
    
    q_infer=sapply(unique(pred$ID_sim),function(x){
      d_fil=filter(pred,ID_sim==x)
      return(d_fil$q[nrow(d_fil)])
    }) 
    
    size_tipping=sapply(unique(pred$ID_sim),function(x){
      d_fil=filter(pred,ID_sim==x)
      if (any(d_fil$cover>0)){
        return(d_fil$cover[which(d_fil$p==d_fil$p[min(which(d_fil$cover !=0))])])
      }else {
        return(NA)
      }
      
    }) 
    
    d=rbind(d,tibble(ID_sim=1:length(p_desert),pcrit=p_desert,pinfer=p_infer,qinfer=q_infer,Size_tipping=size_tipping,
                     Site=site_id))
  }
  
  print(site)
}

write.table(d,"../Data/Inferrence/Prediction/Raw_stability_metrics.csv",sep=";")


## >> 5) Running mixed-effect models ----

d=read.table("../Data/Inferrence/param_inferred.csv",sep=";",header=T)
d_data=read.table("../Data/Spatial_structure_grazing.csv",sep=";")%>%
  Closer_to_normality(.)
n_sites=504

keep_sites=read.table("../Data/Inferrence/Keeping_sites.csv",sep=";")$x

d=tibble(Site=1:n_sites,mean_p=apply(d[,(1:n_sites)],2,mean),sd_p=apply(d[,(1:n_sites)],2,sd),
         mean_q=apply(d[,(n_sites+1):(2*n_sites)],2,mean),sd_q=apply(d[,(n_sites+1):(2*n_sites)],2,sd),
         Plot_n=d_data$Site_ID,
         Aridity=d_data$Aridity,
         Clim1=d_data$Clim1,
         Clim2=d_data$Clim2,
         Clim3=d_data$Clim3,
         Clim4=d_data$Clim4,
         Sp_richness=d_data$Sp_richness,
         Type_veg=d_data$Type_veg,
         Org_C=d_data$Org_C,
         Sand=d_data$Sand,
         Grazing=d_data$Grazing,
         Lattitude=d_data$Lattitude,
         Longitude=d_data$Longitude,
         Elevation=d_data$Elevation,
         Long_sin=d_data$Long_sin,
         Long_cos=d_data$Long_cos,
         Woody=d_data$Woody,
         Slope=d_data$Slope,
         Cover=d_data$rho_p)%>%
  filter(., Site %in% keep_sites)

d2=read.table("../Data/Inferrence/Prediction/Raw_stability_metrics.csv",sep=";")%>%
  group_by(., Site)%>%
  dplyr::summarise(., .groups = "keep",
                   abs_mean=mean(pinfer-pcrit,na.rm = T),
                   abs_sd=sd(pinfer-pcrit,na.rm = T),
                   relativ_mean=mean((pinfer-pcrit)/pinfer,na.rm = T),
                   relativ_sd=sd((pinfer-pcrit)/pinfer,na.rm = T),
                   Size_mean=mean(Size_tipping,na.rm = T),
                   Size_sd=sd(Size_tipping,na.rm = T))%>%
  filter(., Site %in% keep_sites)

d=cbind(d%>%filter(., Site %in% d2$Site),d2)

#as the parameters have uncertainty -> Monte Carlo approach

nsim=1000
d_mod=d_partial=tibble()

for (k in 1:nsim){
  
  d2=tibble(p=scale(sapply(1:nrow(d),function(x){return(rnorm(1,d$mean_p[x],d$sd_p[x]))}))[,1],
            q=scale(sapply(1:nrow(d),function(x){return(rnorm(1,d$mean_q[x],d$sd_q[x]))}))[,1],
            Size_tipping=scale(sapply(1:nrow(d),function(x){return(rnorm(1,d$Size_mean[x],d$Size_sd[x]))}))[,1],
            abs_dist=scale(sapply(1:nrow(d),function(x){return(rnorm(1,d$abs_mean[x],d$abs_sd[x]))}))[,1],
            rela_dist=scale(sapply(1:nrow(d),function(x){return(rnorm(1,d$relativ_mean[x],d$relativ_sd[x]))}))[,1],
            
            #site related variables
            Sand=(d$Sand-mean(d$Sand,na.rm=T))/sd(d$Sand,na.rm = T),
            Sp_richness=(d$Sp_richness-mean(d$Sp_richness,na.rm=T))/sd(d$Sp_richness,na.rm = T),
            Grazing=(d$Grazing-mean(d$Grazing,na.rm=T))/sd(d$Grazing,na.rm = T),
            Site=d$Site,
            Plot_n=d$Plot_n,
            Org_C=(d$Org_C-mean(d$Org_C,na.rm=T))/sd(d$Org_C,na.rm = T),
            Cover=(d$Cover-mean(d$Cover,na.rm=T))/sd(d$Cover,na.rm = T),
            Woody=(d$Woody-mean(d$Woody,na.rm=T))/sd(d$Woody,na.rm = T),
            
            #climatic variables
            Clim1=(d$Clim1-mean(d$Clim1,na.rm=T))/sd(d$Clim1,na.rm = T),
            Clim2=(d$Clim2-mean(d$Clim2,na.rm=T))/sd(d$Clim2,na.rm = T),
            Clim3=(d$Clim3-mean(d$Clim3,na.rm=T))/sd(d$Clim3,na.rm = T),
            Clim4=(d$Clim4-mean(d$Clim4,na.rm=T))/sd(d$Clim4,na.rm = T),
            Aridity=(d$Aridity-mean(d$Aridity,na.rm=T))/sd(d$Aridity,na.rm = T),
            
            #covariates
            Lat=(d$Lattitude -mean(d$Lattitude ,na.rm=T))/sd(d$Lattitude ,na.rm = T),
            Long_cos=(d$Long_cos-mean(d$Long_cos,na.rm=T))/sd(d$Long_cos,na.rm = T),
            Long_sin=(d$Long_sin-mean(d$Long_sin,na.rm=T))/sd(d$Long_sin,na.rm = T),
            Elevation=(d$Elevation-mean(d$Elevation,na.rm=T))/sd(d$Elevation,na.rm = T),
            Slope=(d$Slope-mean(d$Slope,na.rm=T))/sd(d$Slope,na.rm = T))
  
  mod_predictors=gsub("\n     ","","Clim1 + Clim2 + Grazing + Sand + Cover + Sp_richness + Org_C +
      Lat + Long_cos + Long_sin + Slope + Elevation + Grazing * Clim1 + Grazing * Clim2 + Grazing * Org_C + ( 1 | Plot_n)")
  
  model_p=lmer(formula = paste("p ~ ",mod_predictors), data = d2,na.action = na.fail ,REML ="FALSE")
  model_q=lmer(formula = paste("q ~ ",mod_predictors), data = d2,na.action = na.fail ,REML ="FALSE")
  model_abs=lmer(formula = paste("abs_dist ~ ",mod_predictors), data = d2,na.action = na.fail ,REML ="FALSE")
  model_rela=lmer(formula = paste("rela_dist ~ ",mod_predictors), data = d2,na.action = na.fail ,REML ="FALSE")
  model_size=lmer(formula = paste("Size_tipping ~ ",mod_predictors), data = d2,na.action = na.fail ,REML ="FALSE")
  
  #Getting partial prediction
  
  #Full grazing pressure
  resid_mod_size=visreg::visreg(fit = model_size,xvar="Grazing",plot=F) 
  mod_size=lm(visregRes~Grazing,resid_mod_size$res)
  
  resid_mod_rela=visreg::visreg(fit = model_rela,xvar="Grazing",plot=F) 
  mod_rela=lm(visregRes~Grazing,resid_mod_rela$res)
  
  resid_mod_q=visreg::visreg(fit = model_q,xvar="Grazing",plot=F) 
  mod_q=lm(visregRes~Grazing,resid_mod_q$res)
  
  d_partial=rbind(d_partial,
                  tibble(slope=c(summary(mod_size)$coefficient[2,1],
                                 summary(mod_rela)$coefficient[2,1],
                                 summary(mod_q)$coefficient[2,1]),
                         Stat=c("Distance to tipping","Size tipping","q (Spatial structure)"),
                         Type="Full"))

  #Low grazing pressure
  Sub=d2%>%filter(., Grazing %in% c(0,1))
  
  model_p=lmer(formula = paste("p ~ ",mod_predictors), data = Sub,na.action = na.fail ,REML ="FALSE")
  model_q=lmer(formula = paste("q ~ ",mod_predictors), data = Sub,na.action = na.fail ,REML ="FALSE")
  model_abs=lmer(formula = paste("abs_dist ~ ",mod_predictors), data = Sub,na.action = na.fail ,REML ="FALSE")
  model_rela=lmer(formula = paste("rela_dist ~ ",mod_predictors), data = Sub,na.action = na.fail ,REML ="FALSE")
  model_size=lmer(formula = paste("Size_tipping ~ ",mod_predictors), data = Sub,na.action = na.fail ,REML ="FALSE")
  
  resid_mod_size=visreg::visreg(fit = model_size,xvar="Grazing",plot=F) 
  mod_size=lm(visregRes~Grazing,resid_mod_size$res)
  
  resid_mod_rela=visreg::visreg(fit = model_rela,xvar="Grazing",plot=F) 
  mod_rela=lm(visregRes~Grazing,resid_mod_rela$res)
  
  resid_mod_q=visreg::visreg(fit = model_q,xvar="Grazing",plot=F) 
  mod_q=lm(visregRes~Grazing,resid_mod_q$res)
  
  d_partial=rbind(d_partial,
                  tibble(slope=c(summary(mod_size)$coefficient[2,1],
                                 summary(mod_rela)$coefficient[2,1],
                                 summary(mod_q)$coefficient[2,1]),
                         Stat=c("Distance to tipping","Size tipping","q (Spatial structure)"),
                         Type="Low"))

  #High grazing pressure
  Sub = d2%>%filter(., Grazing %in% c(2,3))
  
  model_p=lmer(formula = paste("p ~ ",mod_predictors), data = Sub,na.action = na.fail ,REML ="FALSE")
  model_q=lmer(formula = paste("q ~ ",mod_predictors), data = Sub,na.action = na.fail ,REML ="FALSE")
  model_abs=lmer(formula = paste("abs_dist ~ ",mod_predictors), data = Sub,na.action = na.fail ,REML ="FALSE")
  model_rela=lmer(formula = paste("rela_dist ~ ",mod_predictors), data = Sub,na.action = na.fail ,REML ="FALSE")
  model_size=lmer(formula = paste("Size_tipping ~ ",mod_predictors), data = Sub,na.action = na.fail ,REML ="FALSE")
  
  resid_mod_size=visreg::visreg(fit = model_size,xvar="Grazing",plot=F) 
  mod_size=lm(visregRes~Grazing,resid_mod_size$res)
  
  resid_mod_rela=visreg::visreg(fit = model_rela,xvar="Grazing",plot=F) 
  mod_rela=lm(visregRes~Grazing,resid_mod_rela$res)
  
  resid_mod_q=visreg::visreg(fit = model_q,xvar="Grazing",plot=F) 
  mod_q=lm(visregRes~Grazing,resid_mod_q$res)
  
  d_partial=rbind(d_partial,
                  tibble(slope=c(summary(mod_size)$coefficient[2,1],
                                 summary(mod_rela)$coefficient[2,1],
                                 summary(mod_q)$coefficient[2,1]),
                         Stat=c("Distance to tipping","Size tipping","q (Spatial structure)"),
                         Type="High"))
  
    
  model_p=summary(model_p)
  model_q=summary(model_q)
  model_abs=summary(model_abs)
  model_rela=summary(model_rela)
  model_size=summary(model_size)
  
  d_mod=rbind(d_mod,tibble(Clim1=c(model_p$coefficients[2,1],model_q$coefficients[2,1],model_abs$coefficients[2,1],
                                     model_rela$coefficients[2,1],model_size$coefficients[2,1]),
                           Clim2=c(model_p$coefficients[3,1],model_q$coefficients[3,1],model_abs$coefficients[3,1],
                                                model_rela$coefficients[3,1],model_size$coefficients[3,1]),
                           Grazing=c(model_p$coefficients[4,1],model_q$coefficients[4,1],model_abs$coefficients[4,1],
                                  model_rela$coefficients[4,1],model_size$coefficients[4,1]),
                           Sand=c(model_p$coefficients[5,1],model_q$coefficients[5,1],model_abs$coefficients[5,1],
                                   model_rela$coefficients[5,1],model_size$coefficients[5,1]),
                           Cover=c(model_p$coefficients[6,1],model_q$coefficients[6,1],model_abs$coefficients[6,1],
                                model_rela$coefficients[6,1],model_size$coefficients[6,1]),
                           Sp_richness=c(model_p$coefficients[7,1],model_q$coefficients[7,1],model_abs$coefficients[7,1],
                                   model_rela$coefficients[7,1],model_size$coefficients[7,1]),
                           Org_C=c(model_p$coefficients[8,1],model_q$coefficients[8,1],model_abs$coefficients[8,1],
                                 model_rela$coefficients[8,1],model_size$coefficients[8,1]),
                           Lat=c(model_p$coefficients[9,1],model_q$coefficients[9,1],model_abs$coefficients[9,1],
                                      model_rela$coefficients[9,1],model_size$coefficients[9,1]),
                           Long_cos=c(model_p$coefficients[10,1],model_q$coefficients[10,1],model_abs$coefficients[10,1],
                                      model_rela$coefficients[10,1],model_size$coefficients[10,1]),
                           Long_sin=c(model_p$coefficients[11,1],model_q$coefficients[11,1],model_abs$coefficients[11,1],
                                   model_rela$coefficients[11,1],model_size$coefficients[11,1]),
                           Slope=c(model_p$coefficients[12,1],model_q$coefficients[12,1],model_abs$coefficients[12,1],
                                   model_rela$coefficients[12,1],model_size$coefficients[12,1]),
                           Elevation=c(model_p$coefficients[13,1],model_q$coefficients[13,1],model_abs$coefficients[13,1],
                                   model_rela$coefficients[13,1],model_size$coefficients[13,1]),
                           Clim1_Grazing=c(model_p$coefficients[14,1],model_q$coefficients[14,1],model_abs$coefficients[14,1],
                                   model_rela$coefficients[14,1],model_size$coefficients[14,1]),
                           Clim2_Grazing=c(model_p$coefficients[15,1],model_q$coefficients[15,1],model_abs$coefficients[15,1],
                                           model_rela$coefficients[15,1],model_size$coefficients[15,1]),
                           Grazing_Org_C=c(model_p$coefficients[16,1],model_q$coefficients[16,1],model_abs$coefficients[16,1],
                                           model_rela$coefficients[16,1],model_size$coefficients[16,1]),
                           Param=c("p","q","Absolute distance","Relative distance","Size tipping")),
              )
  
  
  print(k)
  
}

write.table(d_mod,"../Data/Inferrence/Drivers_stability_metrics.csv",sep=";")
write.table(d_partial,"../Data/Inferrence/Partial_residuals_grazing.csv",sep=";")





## >> 6) Running mixed-effect models without cover ----

d=read.table("../Data/Inferrence/param_inferred.csv",sep=";",header=T)
d_data=read.table("../Data/Spatial_structure_grazing.csv",sep=";")%>%
  Closer_to_normality(.)
n_sites=504


keep_sites=read.table("../Data/Inferrence/Keeping_sites.csv",sep=";")$x

d=tibble(Site=1:n_sites,mean_p=apply(d[,(1:n_sites)],2,mean),sd_p=apply(d[,(1:n_sites)],2,sd),
         mean_q=apply(d[,(n_sites+1):(2*n_sites)],2,mean),sd_q=apply(d[,(n_sites+1):(2*n_sites)],2,sd),
         Plot_n=d_data$Site_ID,
         Aridity=d_data$Aridity,
         Clim1=d_data$Clim1,
         Clim2=d_data$Clim2,
         Clim3=d_data$Clim3,
         Clim4=d_data$Clim4,
         Sp_richness=d_data$Sp_richness,
         Type_veg=d_data$Type_veg,
         Org_C=d_data$Org_C,
         Sand=d_data$Sand,
         Grazing=d_data$Grazing,
         Lattitude=d_data$Lattitude,
         Longitude=d_data$Longitude,
         Elevation=d_data$Elevation,
         Long_sin=d_data$Long_sin,
         Long_cos=d_data$Long_cos,
         Woody=d_data$Woody,
         Slope=d_data$Slope,
         Cover=d_data$rho_p)%>%
  filter(., Site %in% keep_sites)

d2=read.table("../Data/Inferrence/Prediction/Raw_stability_metrics.csv",sep=";")%>%
  group_by(., Site)%>%
  dplyr::summarise(., .groups = "keep",
                   abs_mean=mean(pinfer-pcrit,na.rm = T),
                   abs_sd=sd(pinfer-pcrit,na.rm = T),
                   relativ_mean=mean((pinfer-pcrit)/pinfer,na.rm = T),
                   relativ_sd=sd((pinfer-pcrit)/pinfer,na.rm = T),
                   Size_mean=mean(Size_tipping,na.rm = T),
                   Size_sd=sd(Size_tipping,na.rm = T))%>%
  filter(., Site %in% keep_sites)

d=cbind(d%>%filter(., Site %in% d2$Site),d2)

#as the parameters have uncertainty -> Monte Carlo approach

nsim=1000
d_mod=d_partial=tibble()

for (k in 1:nsim){
  
  d2=tibble(p=scale(sapply(1:nrow(d),function(x){return(rnorm(1,d$mean_p[x],d$sd_p[x]))}))[,1],
            q=scale(sapply(1:nrow(d),function(x){return(rnorm(1,d$mean_q[x],d$sd_q[x]))}))[,1],
            Size_tipping=scale(sapply(1:nrow(d),function(x){return(rnorm(1,d$Size_mean[x],d$Size_sd[x]))}))[,1],
            abs_dist=scale(sapply(1:nrow(d),function(x){return(rnorm(1,d$abs_mean[x],d$abs_sd[x]))}))[,1],
            rela_dist=scale(sapply(1:nrow(d),function(x){return(rnorm(1,d$relativ_mean[x],d$relativ_sd[x]))}))[,1],
            
            #site related variables
            Sand=(d$Sand-mean(d$Sand,na.rm=T))/sd(d$Sand,na.rm = T),
            Sp_richness=(d$Sp_richness-mean(d$Sp_richness,na.rm=T))/sd(d$Sp_richness,na.rm = T),
            Grazing=d$Grazing,
            Site=d$Site,
            Plot_n=d$Plot_n,
            Org_C=(d$Org_C-mean(d$Org_C,na.rm=T))/sd(d$Org_C,na.rm = T),
            Cover=(d$Cover-mean(d$Cover,na.rm=T))/sd(d$Cover,na.rm = T),
            Woody=(d$Woody-mean(d$Woody,na.rm=T))/sd(d$Woody,na.rm = T),
            
            #climatic variables
            Clim1=(d$Clim1-mean(d$Clim1,na.rm=T))/sd(d$Clim1,na.rm = T),
            Clim2=(d$Clim2-mean(d$Clim2,na.rm=T))/sd(d$Clim2,na.rm = T),
            Clim3=(d$Clim3-mean(d$Clim3,na.rm=T))/sd(d$Clim3,na.rm = T),
            Clim4=(d$Clim4-mean(d$Clim4,na.rm=T))/sd(d$Clim4,na.rm = T),
            Aridity=(d$Aridity-mean(d$Aridity,na.rm=T))/sd(d$Aridity,na.rm = T),
            
            #covariates
            Lat=(d$Lattitude -mean(d$Lattitude ,na.rm=T))/sd(d$Lattitude ,na.rm = T),
            Long_cos=(d$Long_cos-mean(d$Long_cos,na.rm=T))/sd(d$Long_cos,na.rm = T),
            Long_sin=(d$Long_sin-mean(d$Long_sin,na.rm=T))/sd(d$Long_sin,na.rm = T),
            Elevation=(d$Elevation-mean(d$Elevation,na.rm=T))/sd(d$Elevation,na.rm = T),
            Slope=(d$Slope-mean(d$Slope,na.rm=T))/sd(d$Slope,na.rm = T))
  
  mod_predictors=gsub("\n     ","","Clim1 + Clim2 + Grazing + Sand + Sp_richness + Org_C +
      Lat + Long_cos + Long_sin + Slope + Elevation + Grazing * Clim1 + Grazing * Clim2 + Grazing * Org_C + ( 1 | Plot_n)")
  
  
  model_p=lmer(formula = paste("p ~ ",mod_predictors), data = d2,na.action = na.fail ,REML ="FALSE")
  model_q=lmer(formula = paste("q ~ ",mod_predictors), data = d2,na.action = na.fail ,REML ="FALSE")
  model_abs=lmer(formula = paste("abs_dist ~ ",mod_predictors), data = d2,na.action = na.fail ,REML ="FALSE")
  model_rela=lmer(formula = paste("rela_dist ~ ",mod_predictors), data = d2,na.action = na.fail ,REML ="FALSE")
  model_size=lmer(formula = paste("Size_tipping ~ ",mod_predictors), data = d2,na.action = na.fail ,REML ="FALSE")
  
  #Getting partial prediction
  
  #Full grazing pressure
  resid_mod_size=visreg::visreg(fit = model_size,xvar="Grazing",plot=F) 
  mod_size=lm(visregRes~Grazing,resid_mod_size$res)
  
  resid_mod_rela=visreg::visreg(fit = model_rela,xvar="Grazing",plot=F) 
  mod_rela=lm(visregRes~Grazing,resid_mod_rela$res)
  
  resid_mod_q=visreg::visreg(fit = model_q,xvar="Grazing",plot=F) 
  mod_q=lm(visregRes~Grazing,resid_mod_q$res)
  
  d_partial=rbind(d_partial,
                  tibble(slope=c(summary(mod_size)$coefficient[2,1],
                                 summary(mod_rela)$coefficient[2,1],
                                 summary(mod_q)$coefficient[2,1]),
                         Stat=c("Distance to tipping","Size tipping","q (Spatial structure)"),
                         Type="Full"))
  
  model_p=summary(model_p)
  model_q=summary(model_q)
  model_abs=summary(model_abs)
  model_rela=summary(model_rela)
  model_size=summary(model_size)
  
  
  d_mod=rbind(d_mod,tibble(Clim1=c(model_p$coefficients[2,1],model_q$coefficients[2,1],model_abs$coefficients[2,1],
                                   model_rela$coefficients[2,1],model_size$coefficients[2,1]),
                           Clim2=c(model_p$coefficients[3,1],model_q$coefficients[3,1],model_abs$coefficients[3,1],
                                   model_rela$coefficients[3,1],model_size$coefficients[3,1]),
                           Grazing=c(model_p$coefficients[4,1],model_q$coefficients[4,1],model_abs$coefficients[4,1],
                                     model_rela$coefficients[4,1],model_size$coefficients[4,1]),
                           Sand=c(model_p$coefficients[5,1],model_q$coefficients[5,1],model_abs$coefficients[5,1],
                                  model_rela$coefficients[5,1],model_size$coefficients[5,1]),
                           # Cover=c(model_p$coefficients[6,1],model_q$coefficients[6,1],model_abs$coefficients[6,1],
                           #         model_rela$coefficients[6,1],model_size$coefficients[6,1]),
                           Sp_richness=c(model_p$coefficients[6,1],model_q$coefficients[6,1],model_abs$coefficients[6,1],
                                         model_rela$coefficients[6,1],model_size$coefficients[6,1]),
                           Org_C=c(model_p$coefficients[7,1],model_q$coefficients[7,1],model_abs$coefficients[7,1],
                                   model_rela$coefficients[7,1],model_size$coefficients[7,1]),
                           Lat=c(model_p$coefficients[8,1],model_q$coefficients[8,1],model_abs$coefficients[8,1],
                                 model_rela$coefficients[8,1],model_size$coefficients[8,1]),
                           Long_cos=c(model_p$coefficients[9,1],model_q$coefficients[9,1],model_abs$coefficients[9,1],
                                      model_rela$coefficients[9,1],model_size$coefficients[9,1]),
                           Long_sin=c(model_p$coefficients[10,1],model_q$coefficients[10,1],model_abs$coefficients[10,1],
                                      model_rela$coefficients[10,1],model_size$coefficients[10,1]),
                           Slope=c(model_p$coefficients[11,1],model_q$coefficients[11,1],model_abs$coefficients[11,1],
                                   model_rela$coefficients[11,1],model_size$coefficients[11,1]),
                           Elevation=c(model_p$coefficients[12,1],model_q$coefficients[12,1],model_abs$coefficients[12,1],
                                       model_rela$coefficients[12,1],model_size$coefficients[12,1]),
                           Clim1_Grazing=c(model_p$coefficients[13,1],model_q$coefficients[13,1],model_abs$coefficients[13,1],
                                           model_rela$coefficients[13,1],model_size$coefficients[13,1]),
                           Clim2_Grazing=c(model_p$coefficients[14,1],model_q$coefficients[14,1],model_abs$coefficients[14,1],
                                           model_rela$coefficients[14,1],model_size$coefficients[14,1]),
                           Grazing_Org_C=c(model_p$coefficients[15,1],model_q$coefficients[15,1],model_abs$coefficients[15,1],
                                           model_rela$coefficients[15,1],model_size$coefficients[15,1]),
                           Param=c("p","q","Absolute distance","Relative distance","Size tipping"),
                           Type_grazing="All"))
  
  
  #High grazing pressure
  
  Sub=d2%>%filter(., Grazing %in% c(2,3))
  
  model_p=lmer(formula = paste("p ~ ",mod_predictors), data = Sub,na.action = na.fail ,REML ="FALSE")
  model_q=lmer(formula = paste("q ~ ",mod_predictors), data = Sub,na.action = na.fail ,REML ="FALSE")
  model_abs=lmer(formula = paste("abs_dist ~ ",mod_predictors), data = Sub,na.action = na.fail ,REML ="FALSE")
  model_rela=lmer(formula = paste("rela_dist ~ ",mod_predictors), data = Sub,na.action = na.fail ,REML ="FALSE")
  model_size=lmer(formula = paste("Size_tipping ~ ",mod_predictors), data = Sub,na.action = na.fail ,REML ="FALSE")
  
  
  
  resid_mod_size=visreg::visreg(fit = model_size,xvar="Grazing",plot=F) 
  mod_size=lm(visregRes~Grazing,resid_mod_size$res)
  
  resid_mod_rela=visreg::visreg(fit = model_rela,xvar="Grazing",plot=F) 
  mod_rela=lm(visregRes~Grazing,resid_mod_rela$res)
  
  resid_mod_q=visreg::visreg(fit = model_q,xvar="Grazing",plot=F) 
  mod_q=lm(visregRes~Grazing,resid_mod_q$res)
  
  d_partial=rbind(d_partial,
                  tibble(slope=c(summary(mod_size)$coefficient[2,1],
                                 summary(mod_rela)$coefficient[2,1],
                                 summary(mod_q)$coefficient[2,1]),
                         Stat=c("Distance to tipping","Size tipping","q (Spatial structure)"),
                         Type="High"))
  
  
  
  
  model_p=summary(model_p)
  model_q=summary(model_q)
  model_abs=summary(model_abs)
  model_rela=summary(model_rela)
  model_size=summary(model_size)
  
  d_mod=rbind(d_mod,tibble(Clim1=c(model_p$coefficients[2,1],model_q$coefficients[2,1],model_abs$coefficients[2,1],
                                   model_rela$coefficients[2,1],model_size$coefficients[2,1]),
                           Clim2=c(model_p$coefficients[3,1],model_q$coefficients[3,1],model_abs$coefficients[3,1],
                                   model_rela$coefficients[3,1],model_size$coefficients[3,1]),
                           Grazing=c(model_p$coefficients[4,1],model_q$coefficients[4,1],model_abs$coefficients[4,1],
                                     model_rela$coefficients[4,1],model_size$coefficients[4,1]),
                           Sand=c(model_p$coefficients[5,1],model_q$coefficients[5,1],model_abs$coefficients[5,1],
                                  model_rela$coefficients[5,1],model_size$coefficients[5,1]),
                           # Cover=c(model_p$coefficients[6,1],model_q$coefficients[6,1],model_abs$coefficients[6,1],
                           #         model_rela$coefficients[6,1],model_size$coefficients[6,1]),
                           Sp_richness=c(model_p$coefficients[6,1],model_q$coefficients[6,1],model_abs$coefficients[6,1],
                                         model_rela$coefficients[6,1],model_size$coefficients[6,1]),
                           Org_C=c(model_p$coefficients[7,1],model_q$coefficients[7,1],model_abs$coefficients[7,1],
                                   model_rela$coefficients[7,1],model_size$coefficients[7,1]),
                           Lat=c(model_p$coefficients[8,1],model_q$coefficients[8,1],model_abs$coefficients[8,1],
                                 model_rela$coefficients[8,1],model_size$coefficients[8,1]),
                           Long_cos=c(model_p$coefficients[9,1],model_q$coefficients[9,1],model_abs$coefficients[9,1],
                                      model_rela$coefficients[9,1],model_size$coefficients[9,1]),
                           Long_sin=c(model_p$coefficients[10,1],model_q$coefficients[10,1],model_abs$coefficients[10,1],
                                      model_rela$coefficients[10,1],model_size$coefficients[10,1]),
                           Slope=c(model_p$coefficients[11,1],model_q$coefficients[11,1],model_abs$coefficients[11,1],
                                   model_rela$coefficients[11,1],model_size$coefficients[11,1]),
                           Elevation=c(model_p$coefficients[12,1],model_q$coefficients[12,1],model_abs$coefficients[12,1],
                                       model_rela$coefficients[12,1],model_size$coefficients[12,1]),
                           Clim1_Grazing=c(model_p$coefficients[13,1],model_q$coefficients[13,1],model_abs$coefficients[13,1],
                                           model_rela$coefficients[13,1],model_size$coefficients[13,1]),
                           Clim2_Grazing=c(model_p$coefficients[14,1],model_q$coefficients[14,1],model_abs$coefficients[14,1],
                                           model_rela$coefficients[14,1],model_size$coefficients[14,1]),
                           Grazing_Org_C=c(model_p$coefficients[15,1],model_q$coefficients[15,1],model_abs$coefficients[15,1],
                                           model_rela$coefficients[15,1],model_size$coefficients[15,1]),
                           Param=c("p","q","Absolute distance","Relative distance","Size tipping"),
                           Type_grazing="High"))
  
  
  
  #Low grazing pressure
  
  Sub=d2%>%filter(., Grazing %in% c(0,1))
  
  model_p=lmer(formula = paste("p ~ ",mod_predictors), data = Sub,na.action = na.fail ,REML ="FALSE")
  model_q=lmer(formula = paste("q ~ ",mod_predictors), data = Sub,na.action = na.fail ,REML ="FALSE")
  model_abs=lmer(formula = paste("abs_dist ~ ",mod_predictors), data = Sub,na.action = na.fail ,REML ="FALSE")
  model_rela=lmer(formula = paste("rela_dist ~ ",mod_predictors), data = Sub,na.action = na.fail ,REML ="FALSE")
  model_size=lmer(formula = paste("Size_tipping ~ ",mod_predictors), data = Sub,na.action = na.fail ,REML ="FALSE")
  
  
  
  resid_mod_size=visreg::visreg(fit = model_size,xvar="Grazing",plot=F) 
  mod_size=lm(visregRes~Grazing,resid_mod_size$res)
  
  resid_mod_rela=visreg::visreg(fit = model_rela,xvar="Grazing",plot=F) 
  mod_rela=lm(visregRes~Grazing,resid_mod_rela$res)
  
  resid_mod_q=visreg::visreg(fit = model_q,xvar="Grazing",plot=F) 
  mod_q=lm(visregRes~Grazing,resid_mod_q$res)
  
  d_partial=rbind(d_partial,
                  tibble(slope=c(summary(mod_size)$coefficient[2,1],
                                 summary(mod_rela)$coefficient[2,1],
                                 summary(mod_q)$coefficient[2,1]),
                         Stat=c("Distance to tipping","Size tipping","q (Spatial structure)"),
                         Type="Low"))
  
  
  
  
  model_p=summary(model_p)
  model_q=summary(model_q)
  model_abs=summary(model_abs)
  model_rela=summary(model_rela)
  model_size=summary(model_size)
  
  d_mod=rbind(d_mod,tibble(Clim1=c(model_p$coefficients[2,1],model_q$coefficients[2,1],model_abs$coefficients[2,1],
                                   model_rela$coefficients[2,1],model_size$coefficients[2,1]),
                           Clim2=c(model_p$coefficients[3,1],model_q$coefficients[3,1],model_abs$coefficients[3,1],
                                   model_rela$coefficients[3,1],model_size$coefficients[3,1]),
                           Grazing=c(model_p$coefficients[4,1],model_q$coefficients[4,1],model_abs$coefficients[4,1],
                                     model_rela$coefficients[4,1],model_size$coefficients[4,1]),
                           Sand=c(model_p$coefficients[5,1],model_q$coefficients[5,1],model_abs$coefficients[5,1],
                                  model_rela$coefficients[5,1],model_size$coefficients[5,1]),
                           # Cover=c(model_p$coefficients[6,1],model_q$coefficients[6,1],model_abs$coefficients[6,1],
                           #         model_rela$coefficients[6,1],model_size$coefficients[6,1]),
                           Sp_richness=c(model_p$coefficients[6,1],model_q$coefficients[6,1],model_abs$coefficients[6,1],
                                         model_rela$coefficients[6,1],model_size$coefficients[6,1]),
                           Org_C=c(model_p$coefficients[7,1],model_q$coefficients[7,1],model_abs$coefficients[7,1],
                                   model_rela$coefficients[7,1],model_size$coefficients[7,1]),
                           Lat=c(model_p$coefficients[8,1],model_q$coefficients[8,1],model_abs$coefficients[8,1],
                                 model_rela$coefficients[8,1],model_size$coefficients[8,1]),
                           Long_cos=c(model_p$coefficients[9,1],model_q$coefficients[9,1],model_abs$coefficients[9,1],
                                      model_rela$coefficients[9,1],model_size$coefficients[9,1]),
                           Long_sin=c(model_p$coefficients[10,1],model_q$coefficients[10,1],model_abs$coefficients[10,1],
                                      model_rela$coefficients[10,1],model_size$coefficients[10,1]),
                           Slope=c(model_p$coefficients[11,1],model_q$coefficients[11,1],model_abs$coefficients[11,1],
                                   model_rela$coefficients[11,1],model_size$coefficients[11,1]),
                           Elevation=c(model_p$coefficients[12,1],model_q$coefficients[12,1],model_abs$coefficients[12,1],
                                       model_rela$coefficients[12,1],model_size$coefficients[12,1]),
                           Clim1_Grazing=c(model_p$coefficients[13,1],model_q$coefficients[13,1],model_abs$coefficients[13,1],
                                           model_rela$coefficients[13,1],model_size$coefficients[13,1]),
                           Clim2_Grazing=c(model_p$coefficients[14,1],model_q$coefficients[14,1],model_abs$coefficients[14,1],
                                           model_rela$coefficients[14,1],model_size$coefficients[14,1]),
                           Grazing_Org_C=c(model_p$coefficients[15,1],model_q$coefficients[15,1],model_abs$coefficients[15,1],
                                           model_rela$coefficients[15,1],model_size$coefficients[15,1]),
                           Param=c("p","q","Absolute distance","Relative distance","Size tipping"),
                           Type_grazing="Low"))
  
  
  print(k)
  
}

write.table(d_mod,"../Data/Inferrence/Drivers_stability_metrics_no_cover.csv",sep=";")
write.table(d_partial,"../Data/Inferrence/Partial_residuals_grazing_no_cover.csv",sep=";")




# ---------------------- Step 6: Grazing as a category ----



## >> 1) Run importance ----



Run_model_importance=function(id){
  
  list_mod=expand.grid(with_cover=c(T,F),
                       Stats=c("perim_area_scaling","fmax_psd","Cond_H","PL_expo","Spectral_ratio",
                               "core_area_land","division","fractal_dim","contig","core_area",
                               "flow_length","PLR",
                               "Struct1","Struct2")) #We also add the case where the stats correspond to the first two components of the PCA
  
  d_all=d_all2=tibble()#for keeping all information
  
  d_data=read.table("../Data/Spatial_structure_grazing.csv",sep=";")%>%
    Closer_to_normality(.)
  
  d_data[,c(1,9:29,35,37:41,44:ncol(d_data))] = 
    apply(d_data[,c(1,9:29,35,37:41,44:ncol(d_data))],2,z_tranform)
  
  d_data=Perform_PCA_spatial_struc(d_data)
  d_data$Grazing=as.factor(d_data$Grazing)
  stat=list_mod$Stats[id]
  
  d_data_mod=d_data%>%melt(., measure.vars=stat)%>%
    filter(., !is.na(value))
  
  with_cover=list_mod$with_cover[id]
  
  if (with_cover){
    
    formula_mod=formula(formula_(paste("value ~ Long_cos + Long_sin + Lattitude + Slope + Elevation 
      + Clim1 + Clim2 + Clim3 + Clim4 + Grazing 
      + Woody + Woody*Grazing
      + Clim1*Grazing + Clim2*Grazing + Clim3*Grazing + Clim4*Grazing  
      + rho_p + rho_p*Grazing + Type_veg + Type_veg*Grazing
      + Sand + Sand*Grazing + Org_C + Org_C*Grazing
      + (1|Site_ID)")))
    
  } else{
    
    formula_mod=formula(formula_(paste("value ~ Long_cos + Long_sin + Lattitude + Slope + Elevation 
      + Clim1 + Clim2 + Clim3 + Clim4 + Grazing 
      + Woody + Woody*Grazing
      + Clim1*Grazing + Clim2*Grazing + Clim3*Grazing + Clim4*Grazing  
      + Type_veg + Type_veg*Grazing
      + Sand + Sand*Grazing + Org_C + Org_C*Grazing
      + (1|Site_ID)")))
  }
  
  model_spa_stat  = lmer(formula_mod, d_data_mod,
                         na.action = na.fail ,REML ="FALSE")
  
  #we remove potential outliers
  rm.outliers = romr.fnc(model_spa_stat, d_data_mod, trim=2.5)
  d_data_out = rm.outliers$data
  
  #saving data
  write.table(d_data_out,paste0("../Data/Step1_Understanding_grazing/Factor/Keep_data/Data_",stat,"_",with_cover,".csv"))
  
  model_spa_stat = lmer(formula_mod, data = d_data_out, 
                        na.action = na.fail,REML ="FALSE")
  
  # #do some model selection
  select_model=dredge(model_spa_stat, subset = ~ Slope & Elevation &
                        Long_cos & Long_sin & Lattitude &
                        dc(Woody & Grazing, Woody : Grazing) &
                        dc(Clim1 & Grazing, Clim1 : Grazing) &
                        dc(Clim2 & Grazing, Clim2 : Grazing) &
                        dc(Clim3 & Grazing, Clim3 : Grazing) &
                        dc(Clim4 & Grazing, Clim4 : Grazing) &
                        dc(rho_p & Grazing, rho_p : Grazing) &
                        dc(Org_C & Grazing, Org_C : Grazing) &
                        dc(Sand  & Grazing, Sand  : Grazing) &
                        dc(Type_veg & Grazing, Type_veg : Grazing),
                      options(na.action = "na.fail") )
  
  #extract the result of model selection
  result_select=model.avg(select_model, subset = delta < 2)
  
  #Get the importance of each metric
  importance_mod=sw(result_select)
  
  saveRDS(model_spa_stat,paste0("../Data/Step1_Understanding_grazing/Factor/Keep_models/Mod_",stat,"_",with_cover,".rds"))
  
  R2=select_model%>%filter(., AICc<min(AICc)+2)
  
  d_all=rbind(d_all,cbind(tibble(Sp_stat=stat,
                                 N_outliers=rm.outliers$n.removed,
                                 R2m=mean(R2$R2m),R2C=mean(R2$R2c)),
                          Aggregate_importance(importance_mod))%>%
                add_column(., With_cover=with_cover))
  
  
  
  summary_coef=confint(result_select)
  
  #Merge in a df
  d_all2=tibble(Median=confint(result_select,level = 1e-9)[-1,1],
                q1=summary_coef[-1,1], 
                q3=summary_coef[-1,2],
                pvalue=summary(result_select)$coefmat.full[-1,5],
                term=rownames(summary_coef)[-1],
                Stat=stat,
                R2m=mean(R2$R2m),
                R2c=mean(R2$R2c))%>%
    add_column(., With_cover=with_cover)
  
  write.table(d_all,paste0("../Data/Step1_Understanding_grazing/Factor/Importance/Importance_",stat,"_",with_cover,".csv"),sep=";")
  write.table(d_all2,paste0("../Data/Step1_Understanding_grazing/Factor/Estimator/Estimators_model_",stat,"_",with_cover,".csv"),sep=";")
  
}

library(parallel)

mclapply(1:28,Run_model_importance,mc.cores = 28)


d_all=d_all2=tibble()
for (k in list.files("../Data/Step1_Understanding_grazing/Grazing_factor/Importance","TRUE")){
  d_all=rbind(d_all,read.table(paste0("../Data/Step1_Understanding_grazing/Grazing_factor/Importance/",k),sep=";"))
}
for (k in list.files("../Data/Step1_Understanding_grazing/Grazing_factor/Estimator","TRUE")){
  d_all2=rbind(d_all2,read.table(paste0("../Data/Step1_Understanding_grazing/Grazing_factor/Estimator/",k),sep=";"))
}
write.table(d_all,"../Data/Step1_Understanding_grazing/Grazing_factor/Importance.csv",sep=";")
write.table(d_all2,"../Data/Step1_Understanding_grazing/Grazing_factor/Estimators_model.csv",sep=";")



d_all=d_all2=tibble()
for (k in list.files("../Data/Step1_Understanding_grazing/Grazing_factor/Importance","FALSE")){
  d_all=rbind(d_all,read.table(paste0("../Data/Step1_Understanding_grazing/Grazing_factor/Importance/",k),sep=";"))
}
for (k in list.files("../Data/Step1_Understanding_grazing/Grazing_factor/Estimator","FALSE")){
  d_all2=rbind(d_all2,read.table(paste0("../Data/Step1_Understanding_grazing/Grazing_factor/Estimator/",k),sep=";"))
}
write.table(d_all,"../Data/Step1_Understanding_grazing/Grazing_factor/Importance_no_cov.csv",sep=";")
write.table(d_all2,"../Data/Step1_Understanding_grazing/Grazing_factor/Estimators_model_no_cov.csv",sep=";")


## >> 2) Comparing AIC ----


d_data=read.table("../Data/Spatial_structure_grazing.csv",sep=";")%>%
  Closer_to_normality(.)
d_data[,c(1,9:29,35,37:41,44:ncol(d_data))] = 
  apply(d_data[,c(1,9:29,35,37:41,44:ncol(d_data))],2,z_tranform)
d_data=Perform_PCA_spatial_struc(d_data)

d_AIC=tibble()

for (k in c("perim_area_scaling","fmax_psd","PL_expo","Spectral_ratio",
            "core_area_land","division","fractal_dim","contig","core_area",
            "flow_length","PLR",
            "Struct1","Struct2")){
  
  #for each we plot the slope against the partial residuals with and without cover
  
  d_data$Grazing=as.factor(d_data$Grazing)
  
  d_data_out=read.table(paste0("../Data/Step1_Understanding_grazing/Grazing_factor/Keep_data/Data_",k,"_TRUE.csv"),sep=" ")
  model_spa_stat=readRDS(paste0("../Data/Step1_Understanding_grazing/Grazing_factor/Keep_models/Mod_",k,"_TRUE.rds"))
  d_data$Grazing=as.numeric(d_data$Grazing)
  d_data_out=read.table(paste0("../Data/Step1_Understanding_grazing/Keep_data/Data_",k,"_TRUE.csv"),sep=" ")
  model_spa_stat2=readRDS(paste0("../Data/Step1_Understanding_grazing/Keep_models/Mod_",k,"_TRUE.rds"))
  d_AIC=rbind(d_AIC,tibble(AIC_factor=AIC(model_spa_stat,model_spa_stat2)[1,2],AIC_num=AIC(model_spa_stat,model_spa_stat2)[2,2])%>%
                add_column(., Stat=k))
  
}
d_AIC

## >> 3) Analyse residuals ----

pdf("../Figures/Step1_Understanding_grazing/Grazing_factor/Test_factor_with_cov.pdf",width = 6,height = 4)


d_data=read.table("../Data/Spatial_structure_grazing.csv",sep=";")%>%
  Closer_to_normality(.)
d_data[,c(1,9:29,35,37:41,44:ncol(d_data))] = 
  apply(d_data[,c(1,9:29,35,37:41,44:ncol(d_data))],2,z_tranform)
d_data=Perform_PCA_spatial_struc(d_data)

d_resid=tibble()
for (k in c("perim_area_scaling","fmax_psd","PL_expo","Spectral_ratio",
            "core_area_land","division","fractal_dim","contig","core_area",
            "flow_length","PLR",
            "Struct1","Struct2")){
  
  #for each we plot the slope against the partial residuals with and without cover
  
  d_data$Grazing=as.factor(d_data$Grazing)
  
  d_data_out=read.table(paste0("../Data/Step1_Understanding_grazing/Grazing_factor/Keep_data/Data_",k,"_TRUE.csv"),sep=" ")
  d_data_out$Grazing=as.factor(d_data_out$Grazing)
  model_spa_stat=readRDS(paste0("../Data/Step1_Understanding_grazing/Grazing_factor/Keep_models/Mod_",k,"_TRUE.rds"))
  
  
  resid_mod=visreg::visreg(fit = model_spa_stat,xvar="Grazing",plot=F,ylab=as.character(k)) 
  
  mod_res=lm("visregRes ~ Grazing",resid_mod$res)
  
  d_resid=rbind(d_resid,tibble(Low=confint(mod_res)[-1,1],High=confint(mod_res)[-1,2],Stat=k,Name=c("Low","Medium","High")))
  
  if ("Type_vegGrassland" %in% rownames(summary(model_spa_stat)$coefficients)){
    resid_mod=visreg::visreg(fit = model_spa_stat,xvar="Grazing",by="Type_veg",plot=T,ylab=as.character(k))
  }
  if ("Clim1" %in% rownames(summary(model_spa_stat)$coefficients)){
    metric="Clim1"
    resid_mod=visreg::visreg(fit = model_spa_stat,xvar="Clim1",by="Grazing",plot=F,ylab=as.character(k))
    print(ggplot(resid_mod$res%>%
                   mutate(., Grazing=as.character(Grazing))%>%
                   dplyr::select(., -value)%>%
                   melt(., measure.vars=metric)%>%arrange(., Grazing,decreasing=T))+
            geom_point(aes(value,visregRes),color="gray30",alpha=.2)+
            geom_smooth(aes(value,visregRes,fill=Grazing,color=Grazing),method = "lm",level=.95)+
            the_theme+
            scale_color_manual(values=c("0"="#FFC170",'1'="#43FFF4","2"="#5196F9","3"="#012E71"),
                               labels=c("Ungrazed","Low","Medium","High"))+
            scale_fill_manual(values=c("0"="#FFC170",'1'="#43FFF4","2"="#5196F9","3"="#012E71"),
                              labels=c("Ungrazed","Low","Medium","High"))+
            labs(x=metric))
    
  }
  if ("Clim2" %in% rownames(summary(model_spa_stat)$coefficients)){
    metric="Clim2"
    resid_mod=visreg::visreg(fit = model_spa_stat,xvar="Clim2",by="Grazing",plot=F,ylab=as.character(k))
    print(ggplot(resid_mod$res%>%
                   mutate(., Grazing=as.character(Grazing))%>%
                   dplyr::select(., -value)%>%
                   melt(., measure.vars=metric)%>%arrange(., Grazing,decreasing=T))+
            geom_point(aes(value,visregRes),color="gray30",alpha=.2)+
            geom_smooth(aes(value,visregRes,fill=Grazing,color=Grazing),method = "lm",level=.95)+
            the_theme+
            scale_color_manual(values=c("0"="#FFC170",'1'="#43FFF4","2"="#5196F9","3"="#012E71"),
                               labels=c("Ungrazed","Low","Medium","High"))+
            scale_fill_manual(values=c("0"="#FFC170",'1'="#43FFF4","2"="#5196F9","3"="#012E71"),
                              labels=c("Ungrazed","Low","Medium","High"))+
            labs(x=metric))
    
  }
  if ("Clim3" %in% rownames(summary(model_spa_stat)$coefficients)){
    metric="Clim3"
    resid_mod=visreg::visreg(fit = model_spa_stat,xvar="Clim3",by="Grazing",plot=F,ylab=as.character(k))
    print(ggplot(resid_mod$res%>%
                   mutate(., Grazing=as.character(Grazing))%>%
                   dplyr::select(., -value)%>%
                   melt(., measure.vars=metric)%>%arrange(., Grazing,decreasing=T))+
            geom_point(aes(value,visregRes),color="gray30",alpha=.2)+
            geom_smooth(aes(value,visregRes,fill=Grazing,color=Grazing),method = "lm",level=.95)+
            the_theme+
            scale_color_manual(values=c("0"="#FFC170",'1'="#43FFF4","2"="#5196F9","3"="#012E71"),
                               labels=c("Ungrazed","Low","Medium","High"))+
            scale_fill_manual(values=c("0"="#FFC170",'1'="#43FFF4","2"="#5196F9","3"="#012E71"),
                              labels=c("Ungrazed","Low","Medium","High"))+
            labs(x=metric))
    
  }
  if ("Clim4" %in% rownames(summary(model_spa_stat)$coefficients)){
    metric="Clim4"
    resid_mod=visreg::visreg(fit = model_spa_stat,xvar="Clim4",by="Grazing",plot=F,ylab=as.character(k))
    print(ggplot(resid_mod$res%>%
                   mutate(., Grazing=as.character(Grazing))%>%
                   dplyr::select(., -value)%>%
                   melt(., measure.vars=metric)%>%arrange(., Grazing,decreasing=T))+
            geom_point(aes(value,visregRes),color="gray30",alpha=.2)+
            geom_smooth(aes(value,visregRes,fill=Grazing,color=Grazing),method = "lm",level=.95)+
            the_theme+
            scale_color_manual(values=c("0"="#FFC170",'1'="#43FFF4","2"="#5196F9","3"="#012E71"),
                               labels=c("Ungrazed","Low","Medium","High"))+
            scale_fill_manual(values=c("0"="#FFC170",'1'="#43FFF4","2"="#5196F9","3"="#012E71"),
                              labels=c("Ungrazed","Low","Medium","High"))+
            labs(x=metric))
    
  }
  if ("Org_C" %in% rownames(summary(model_spa_stat)$coefficients)){
    metric="Org_C"
    resid_mod=visreg::visreg(fit = model_spa_stat,xvar="Org_C",by="Grazing",plot=F,ylab=as.character(k))
    print(ggplot(resid_mod$res%>%
                   mutate(., Grazing=as.character(Grazing))%>%
                   dplyr::select(., -value)%>%
                   melt(., measure.vars=metric)%>%arrange(., Grazing,decreasing=T))+
            geom_point(aes(value,visregRes),color="gray30",alpha=.2)+
            geom_smooth(aes(value,visregRes,fill=Grazing,color=Grazing),method = "lm",level=.95)+
            the_theme+
            scale_color_manual(values=c("0"="#FFC170",'1'="#43FFF4","2"="#5196F9","3"="#012E71"),
                               labels=c("Ungrazed","Low","Medium","High"))+
            scale_fill_manual(values=c("0"="#FFC170",'1'="#43FFF4","2"="#5196F9","3"="#012E71"),
                              labels=c("Ungrazed","Low","Medium","High"))+
            labs(x=metric))
    
  }
  
  
  
}
dev.off()


pdf("../Figures/Step1_Understanding_grazing/Grazing_factor/Test_factor_no_cov.pdf",width = 6,height = 4)


d_data=read.table("../Data/Spatial_structure_grazing.csv",sep=";")%>%
  Closer_to_normality(.)
d_data[,c(1,9:29,35,37:41,44:ncol(d_data))] = 
  apply(d_data[,c(1,9:29,35,37:41,44:ncol(d_data))],2,z_tranform)
d_data=Perform_PCA_spatial_struc(d_data)

save=d_data

for (k in c("perim_area_scaling","fmax_psd","PL_expo","Spectral_ratio",
            "core_area_land","division","fractal_dim","contig","core_area",
            "flow_length","PLR",
            "Struct1","Struct2")){
  
  #for each we plot the slope against the partial residuals with and without cover
  
  
  d_data_out=read.table(paste0("../Data/Step1_Understanding_grazing/Grazing_factor/Keep_data/Data_",k,"_FALSE.csv"),sep=" ")
  d_data_out$Grazing=as.factor(d_data_out$Grazing)
  model_spa_stat=readRDS(paste0("../Data/Step1_Understanding_grazing/Grazing_factor/Keep_models/Mod_",k,"_FALSE.rds"))
  
  
  resid_mod=visreg::visreg(fit = model_spa_stat,xvar="Grazing",plot=T,ylab=as.character(k))
  
  
  if ("Type_vegGrassland" %in% rownames(summary(model_spa_stat)$coefficients)){
    resid_mod=visreg::visreg(fit = model_spa_stat,xvar="Grazing",by="Type_veg",plot=T,ylab=as.character(k))
  }
  if ("Clim1" %in% rownames(summary(model_spa_stat)$coefficients)){
    metric="Clim1"
    resid_mod=visreg::visreg(fit = model_spa_stat,xvar="Clim1",by="Grazing",plot=F,ylab=as.character(k))
    print(ggplot(resid_mod$res%>%
                   mutate(., Grazing=as.character(Grazing))%>%
                   dplyr::select(., -value)%>%
                   melt(., measure.vars=metric)%>%arrange(., Grazing,decreasing=T))+
            geom_point(aes(value,visregRes),color="gray30",alpha=.2)+
            geom_smooth(aes(value,visregRes,fill=Grazing,color=Grazing),method = "lm",level=.95)+
            the_theme+
            scale_color_manual(values=c("0"="#FFC170",'1'="#43FFF4","2"="#5196F9","3"="#012E71"),
                               labels=c("Ungrazed","Low","Medium","High"))+
            scale_fill_manual(values=c("0"="#FFC170",'1'="#43FFF4","2"="#5196F9","3"="#012E71"),
                              labels=c("Ungrazed","Low","Medium","High"))+
            labs(x=metric))
    
  }
  if ("Clim2" %in% rownames(summary(model_spa_stat)$coefficients)){
    metric="Clim2"
    resid_mod=visreg::visreg(fit = model_spa_stat,xvar="Clim2",by="Grazing",plot=F,ylab=as.character(k))
    print(ggplot(resid_mod$res%>%
                   mutate(., Grazing=as.character(Grazing))%>%
                   dplyr::select(., -value)%>%
                   melt(., measure.vars=metric)%>%arrange(., Grazing,decreasing=T))+
            geom_point(aes(value,visregRes),color="gray30",alpha=.2)+
            geom_smooth(aes(value,visregRes,fill=Grazing,color=Grazing),method = "lm",level=.95)+
            the_theme+
            scale_color_manual(values=c("0"="#FFC170",'1'="#43FFF4","2"="#5196F9","3"="#012E71"),
                               labels=c("Ungrazed","Low","Medium","High"))+
            scale_fill_manual(values=c("0"="#FFC170",'1'="#43FFF4","2"="#5196F9","3"="#012E71"),
                              labels=c("Ungrazed","Low","Medium","High"))+
            labs(x=metric))
    
  }
  if ("Clim3" %in% rownames(summary(model_spa_stat)$coefficients)){
    metric="Clim3"
    resid_mod=visreg::visreg(fit = model_spa_stat,xvar="Clim3",by="Grazing",plot=F,ylab=as.character(k))
    print(ggplot(resid_mod$res%>%
                   mutate(., Grazing=as.character(Grazing))%>%
                   dplyr::select(., -value)%>%
                   melt(., measure.vars=metric)%>%arrange(., Grazing,decreasing=T))+
            geom_point(aes(value,visregRes),color="gray30",alpha=.2)+
            geom_smooth(aes(value,visregRes,fill=Grazing,color=Grazing),method = "lm",level=.95)+
            the_theme+
            scale_color_manual(values=c("0"="#FFC170",'1'="#43FFF4","2"="#5196F9","3"="#012E71"),
                               labels=c("Ungrazed","Low","Medium","High"))+
            scale_fill_manual(values=c("0"="#FFC170",'1'="#43FFF4","2"="#5196F9","3"="#012E71"),
                              labels=c("Ungrazed","Low","Medium","High"))+
            labs(x=metric))
    
  }
  if ("Clim4" %in% rownames(summary(model_spa_stat)$coefficients)){
    metric="Clim4"
    resid_mod=visreg::visreg(fit = model_spa_stat,xvar="Clim4",by="Grazing",plot=F,ylab=as.character(k))
    print(ggplot(resid_mod$res%>%
                   mutate(., Grazing=as.character(Grazing))%>%
                   dplyr::select(., -value)%>%
                   melt(., measure.vars=metric)%>%arrange(., Grazing,decreasing=T))+
            geom_point(aes(value,visregRes),color="gray30",alpha=.2)+
            geom_smooth(aes(value,visregRes,fill=Grazing,color=Grazing),method = "lm",level=.95)+
            the_theme+
            scale_color_manual(values=c("0"="#FFC170",'1'="#43FFF4","2"="#5196F9","3"="#012E71"),
                               labels=c("Ungrazed","Low","Medium","High"))+
            scale_fill_manual(values=c("0"="#FFC170",'1'="#43FFF4","2"="#5196F9","3"="#012E71"),
                              labels=c("Ungrazed","Low","Medium","High"))+
            labs(x=metric))
    
  }
  if ("Org_C" %in% rownames(summary(model_spa_stat)$coefficients)){
    metric="Org_C"
    resid_mod=visreg::visreg(fit = model_spa_stat,xvar="Org_C",by="Grazing",plot=F,ylab=as.character(k))
    print(ggplot(resid_mod$res%>%
                   mutate(., Grazing=as.character(Grazing))%>%
                   dplyr::select(., -value)%>%
                   melt(., measure.vars=metric)%>%arrange(., Grazing,decreasing=T))+
            geom_point(aes(value,visregRes),color="gray30",alpha=.2)+
            geom_smooth(aes(value,visregRes,fill=Grazing,color=Grazing),method = "lm",level=.95)+
            the_theme+
            scale_color_manual(values=c("0"="#FFC170",'1'="#43FFF4","2"="#5196F9","3"="#012E71"),
                               labels=c("Ungrazed","Low","Medium","High"))+
            scale_fill_manual(values=c("0"="#FFC170",'1'="#43FFF4","2"="#5196F9","3"="#012E71"),
                              labels=c("Ungrazed","Low","Medium","High"))+
            labs(x=metric))
    
  }
  
  
  
}
dev.off()



## >> 4) Ploting residuals ----


d_data=read.table("../Data/Spatial_structure_grazing.csv",sep=";")%>%
  Closer_to_normality(.)
d_data[,c(1,9:29,35,37:41,44:ncol(d_data))] = 
  apply(d_data[,c(1,9:29,35,37:41,44:ncol(d_data))],2,z_tranform)
d_data=Perform_PCA_spatial_struc(d_data)

d_resid=tibble()
for (k in c("perim_area_scaling","fmax_psd","PL_expo","Spectral_ratio",
            "core_area_land","division","fractal_dim","contig","core_area",
            "flow_length","PLR",
            "Struct1","Struct2")){
  
  #for each we plot the slope against the partial residuals with and without cover
  
  d_data$Grazing=as.factor(d_data$Grazing)
  
  d_data_out=read.table(paste0("../Data/Step1_Understanding_grazing/Grazing_factor/Keep_data/Data_",k,"_TRUE.csv"),sep=" ")
  d_data_out$Grazing=as.factor(d_data_out$Grazing)
  model_spa_stat=readRDS(paste0("../Data/Step1_Understanding_grazing/Grazing_factor/Keep_models/Mod_",k,"_TRUE.rds"))
  
  
  resid_mod=visreg::visreg(fit = model_spa_stat,xvar="Grazing",plot=F,ylab=as.character(k)) 
  
  mod_res=lm("visregRes ~ Grazing",resid_mod$res)
  
  d_resid=rbind(d_resid,tibble(Low=confint(mod_res)[-1,1],
                               High=confint(mod_res)[-1,2],
                               Stat=k,Name=c("Low","Medium","High"),Cover=T))
  
  
  d_data$Grazing=as.factor(d_data$Grazing)
  
  d_data_out=read.table(paste0("../Data/Step1_Understanding_grazing/Grazing_factor/Keep_data/Data_",k,"_FALSE.csv"),sep=" ")
  d_data_out$Grazing=as.factor(d_data_out$Grazing)
  model_spa_stat=readRDS(paste0("../Data/Step1_Understanding_grazing/Grazing_factor/Keep_models/Mod_",k,"_FALSE.rds"))
  
  
  resid_mod=visreg::visreg(fit = model_spa_stat,xvar="Grazing",plot=F,ylab=as.character(k)) 
  
  mod_res=lm("visregRes ~ Grazing",resid_mod$res)
  
  d_resid=rbind(d_resid,tibble(Low=confint(mod_res)[-1,1],
                               High=confint(mod_res)[-1,2],
                               Stat=k,Name=c("Low","Medium","High"),Cover=F))
  
}

p=ggplot(d_resid%>%add_column(., median=(.$Low+.$High)/2)%>%Filter_relevant_stats(.)%>%Rename_spatial_statistics(.))+
  geom_pointrange(aes(y=Stat,xmin=Low,xmax=High,x=median,color=Name))+
  geom_vline(xintercept = 0,linetype=9)+
  scale_color_manual(values=c("red","orange","green"))+
  the_theme+
  facet_wrap(.~Cover,labeller = label_bquote(cols = "Cover"==.(as.numeric(Cover))))+
  labs(color="Grazing intensity",y="")

ggsave("../Figures/Step1_Understanding_grazing/Grazing_factor/Residuals_with_factor_grazing.pdf",p,width = 7,height = 5)


