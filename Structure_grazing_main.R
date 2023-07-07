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

NA_kept=250
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
  
  write.table(d_param_infer_rej,paste0("../Data/Sensitivity_cover/param_rej.csv"),sep=";")
}


nsite=100
param_infered=read.table("../Data/Sensitivity_cover/param_rej.csv",sep=";")
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

list_f=list.files("../Data/Metrics/",".csv")
d=tibble()
for (k in list_f){d=rbind(d,read.table(paste0("../Data/Metrics/",k),sep=";",header = T))}

d=d%>%
  add_column(.,
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
  mutate(., Type_veg=recode_factor(Type_veg,"1"="Grassland","2"="Shrubland","3"="Forest","1_2"="Grass_Shrub"))


#Then, we extract summarized climatic variables using a PCA on all climatic variables

clim_variables=colnames(d_biodesert)[38:58]
res.comp=imputePCA(d_biodesert[,which(colnames(d_biodesert) %in% clim_variables)],ncp=4,scale = T) 

if ("completeObs" %in% names(res.comp)){
  res.pca=PCA(res.comp$completeObs, ncp = 4,  graph=F)
}else {
  res.pca=PCA(res.comp, ncp = 4,  graph=F)
}

fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 50))
fviz_pca_var(res.pca,col.var="#2D569A",axes = c(1,2))+ #the ones that are the most imporant
  the_theme

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
             Org_C=as.numeric(sapply(1:nrow(.),function(x){return(d_biodesert$`ORC veg`[which(d_biodesert$ID==.$Site_ID[x])]
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
                                                                      )}))
  )



#transforming into true Site ID
d$Site_ID=sapply(1:nrow(d),function(x){return(d_biodesert$SITE_ID[which(d_biodesert$ID==d$Site_ID[x])])})

write.table(d%>%dplyr::rename(., fractal_dim=factal_dim),
            "../Data/Spatial_structure_grazing.csv",sep=";")

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

d_data=read.table("../Data/Spatial_structure_grazing.csv",sep=";")%>%
  Closer_to_normality(.)

d_data[,c(1,9:29,35,37:41,44:ncol(d_data))] = apply(d_data[,c(1,9:29,35,37:41,44:ncol(d_data))],2,z_tranform)


#correlations between predictors ?
p=ggplot(cor(d_data[,c("Sand","Clim1","Clim2","Clim3","Clim4","Org_C","Woody",
                     "Grazing","rho_p","Long_cos","Long_sin","Lattitude","Slope","Elevation")],
           use =  "na.or.complete")%>%
         melt(.))+geom_tile(aes(Var1,Var2,fill=value))+
  scale_fill_gradient2(low="red","white","blue")+
  labs(x="",y="")+
  theme(axis.text.x=element_text(angle=60,hjust=1))
  
ggsave("../Figures/Step1_Understanding_grazing/Correlation_btw_predictors.pdf",p,width = 5,height = 4)




## Runnning the models

dir.create("../Data/Step1_Understanding_grazing/VIF",showWarnings = F)
dir.create("../Data/Step1_Understanding_grazing/Auto_corr",showWarnings = F)
dir.create("../Data/Step1_Understanding_grazing/Importance",showWarnings = F)
dir.create("../Data/Step1_Understanding_grazing/Estimator",showWarnings = F)
dir.create("../Data/Step1_Understanding_grazing/Keep_models",showWarnings = F)


Run_model_importance=function(id){
  
  list_mod=expand.grid(with_cover=c(T),
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
  result_select=summary(model.avg(select_model, subset = delta < 2))
  
  #Get the importance of each metric
  importance_mod=sw(result_select)
  
  #Get R² of the best model
  model_spa_stat=MuMIn::get.models(select_model, subset = delta == 0)[[1]]
  
  saveRDS(model_spa_stat,paste0("../Data/Step1_Understanding_grazing/Keep_models/Mod_",stat,"_",with_cover,".rds"))
  
  R2_mod=r.squaredGLMM(model_spa_stat)
  
  d_all=rbind(d_all,cbind(tibble(Sp_stat=stat,
                                 N_outliers=rm.outliers$n.removed,
                                 R2m=R2_mod[1],R2C=R2_mod[2]),
                          Aggregate_importance(importance_mod))%>%
                add_column(., With_cover=with_cover))
  
  #bootstraping slopes
  boot_mod=bootstrap(model_spa_stat,.f=fixef,type = "parametric",500)
  
  #Merge in a df
  d_all2=rbind(d_all2,
               tibble(Median=apply(boot_mod$replicates[,-1],2,median),
                      q1=apply(boot_mod$replicates[,-1],2,quantile,.025), 
                      q3=apply(boot_mod$replicates[,-1],2,quantile,.975),
                      pvalue=apply(boot_mod$replicates[,-1],2,twoside_pvalue),
                      term=colnames(boot_mod$replicates)[-1])%>%
                 add_column(Stat=stat)%>%
                 add_column(., With_cover=with_cover))
  
  
  
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

mclapply(1:14,Run_model_importance,mc.cores = 14)


d_all=d_all2=tibble()
for (k in list.files("../Data/Step1_Understanding_grazing/Importance")){
  d_all=rbind(d_all,read.table(paste0("../Data/Step1_Understanding_grazing/Importance/",k),sep=";"))
}
for (k in list.files("../Data/Step1_Understanding_grazing/Estimator")){
  d_all2=rbind(d_all2,read.table(paste0("../Data/Step1_Understanding_grazing/Estimator/",k),sep=";"))
}
write.table(d_all,"../Data/Step1_Understanding_grazing/Importance.csv",sep=";")
write.table(d_all2,"../Data/Step1_Understanding_grazing/Estimators_model.csv",sep=";")


Run_model_importance_Grazing_square=function(id){
  
  list_mod=expand.grid(with_cover=c(T),
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
  
  stat=list_mod$Stats[id]
  
  d_data_mod=d_data%>%melt(., measure.vars=stat)%>%
    filter(., !is.na(value))
  
  with_cover=list_mod$with_cover[id]
  
  if (with_cover){
    
    formula_mod=formula(formula_(paste("value ~ Long_cos + Long_sin + Lattitude + Slope + Elevation 
      + Clim1 + Clim2 + Clim3 + Clim4 + Grazing 
      + Woody + Woody*Grazing + I(Grazing^2)+
      + Clim1*Grazing + Clim2*Grazing + Clim3*Grazing + Clim4*Grazing  
      + rho_p + rho_p*Grazing + Type_veg + Type_veg*Grazing
      + Sand + Sand*Grazing + Org_C + Org_C*Grazing
      + (1|Site_ID)")))
    
  } else{
    
    formula_mod=formula(formula_(paste("value ~ Long_cos + Long_sin + Lattitude + Slope + Elevation 
      + Clim1 + Clim2 + Clim3 + Clim4 + Grazing 
      + Woody + Woody*Grazing + I(Grazing^2)+
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
                        dc(Grazing, I(Grazing^2)) &
                        dc(Type_veg & Grazing, Type_veg : Grazing),
                      options(na.action = "na.fail") )
  
  #extract the result of model selection
  result_select=summary(model.avg(select_model, subset = delta < 2))
  
  #Get the importance of each metric
  importance_mod=sw(result_select)
  
  #Get R² of the best model
  model_spa_stat=MuMIn::get.models(select_model, subset = delta == 0)[[1]]
  
  saveRDS(model_spa_stat,paste0("../Data/Step1_Understanding_grazing/Keep_models/Mod_",stat,"_",with_cover,"Grazing_square.rds"))
  
  R2_mod=r.squaredGLMM(model_spa_stat)
  
  d_all=rbind(d_all,cbind(tibble(Sp_stat=stat,
                                 N_outliers=rm.outliers$n.removed,
                                 R2m=R2_mod[1],R2C=R2_mod[2]),
                          Aggregate_importance(importance_mod))%>%
                add_column(., With_cover=with_cover))
  
  #bootstraping slopes
  boot_mod=bootstrap(model_spa_stat,.f=fixef,type = "parametric",500)
  
  #Merge in a df
  d_all2=rbind(d_all2,
               tibble(Median=apply(boot_mod$replicates[,-1],2,median),
                      q1=apply(boot_mod$replicates[,-1],2,quantile,.025), 
                      q3=apply(boot_mod$replicates[,-1],2,quantile,.975),
                      pvalue=apply(boot_mod$replicates[,-1],2,twoside_pvalue),
                      term=colnames(boot_mod$replicates)[-1])%>%
                 add_column(Stat=stat)%>%
                 add_column(., With_cover=with_cover))
  
  
  #Saving each DF
  write.table(d_all,paste0("../Data/Step1_Understanding_grazing/Importance/Importance_",stat,"_",with_cover,"grazing_square.csv"),sep=";")
  write.table(d_all2,paste0("../Data/Step1_Understanding_grazing/Estimator/Estimators_model_",stat,"_",with_cover,"grazing_square.csv"),sep=";")
  
}

library(parallel)

mclapply(1:14,Run_model_importance_Grazing_square,mc.cores = 14)


d_all=d_all2=tibble()
for (k in list.files("../Data/Step1_Understanding_grazing/Importance","grazing")){
  d_all=rbind(d_all,read.table(paste0("../Data/Step1_Understanding_grazing/Importance/",k),sep=";"))
}
for (k in list.files("../Data/Step1_Understanding_grazing/Estimator","grazing")){
  d_all2=rbind(d_all2,read.table(paste0("../Data/Step1_Understanding_grazing/Estimator/",k),sep=";"))
}
write.table(d_all,"../Data/Step1_Understanding_grazing/Importance_grazing_square.csv",sep=";")
write.table(d_all2,"../Data/Step1_Understanding_grazing/Estimators_model_grazing_square.csv",sep=";")





