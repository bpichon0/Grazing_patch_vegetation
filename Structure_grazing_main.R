rm(list=ls())
source("./Structure_grazing_function.R")

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
  
  #and binarize the kmean output 
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

## >> 2) Computing the metrics on the binary landscapes ----

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
mclapply(list_f[grep(list_f,pattern = "biodesert")],Compute_all_metrics_biodesert,mc.cores = 50)


# 
# dir.create("../Data/Metric_null",showWarnings = F)
# 
# Compute_all_metrics_biodesert_null=function(k){
# 
#   landscape=as.matrix(read.table(paste0("../Data/Landscapes/Binary_landscapes/",k),sep=","))
#   name=gsub(".csv","",k)
#   name_dataset=strsplit(k,"_")[[1]][1]
# 
#   if (any(grep("cropped",k))){
#     id_site=as.numeric(strsplit(k,"_")[[1]][4])
#     sub_id=as.numeric(gsub(".csv","",(strsplit(k,"_")[[1]][5])))
#   }else{
#     id_site=as.numeric(strsplit(k,"_")[[1]][3])
#     sub_id=as.numeric(gsub(".csv","",(strsplit(k,"_")[[1]][4])))
#   }
# 
#   slope=d_biodesert$`SLOPE-ALOS30`[which(d_biodesert$ID==id_site)]
#   
#   d_null=tibble()
#   
#   for (k in 1:199){ #each random landscape
#     null_mat = matrix(sample(landscape), nrow = nrow(landscape), ncol = ncol(landscape))
#     spatial_metrics=Get_sumstat(null_mat,slope=slope,log_ = T)
#     
#     d_null=rbind(d_null,spatial_metrics%>%add_column(., N_random=k))
#   }
# 
#   write.table(d_null,paste0("../Data/Metrics_null/Metric_",name,".csv"),sep=";")
# }
# 
# info_kmean=read.table("../Data/Landscapes/Kmean_clust_info.csv",sep=",",header = T)%>%
#   filter(., Size!=200,Dataset=="biodesert",status=="kept") #keeping the kept sites
# 
# 
# library(parallel)
# mclapply(paste0(info_kmean$Site_ID,"_",info_kmean$Image,".csv"),Compute_all_metrics_biodesert_null,mc.cores = 70)
# 

## aggregating the metrics into a df

info_kmean=read.table("../Data/Landscapes/Kmean_clust_info.csv",sep=",",header = T)%>%
  filter(., Size!=200,Dataset=="biodesert",status=="kept") #keeping the kept sites

d=tibble()
for (k in paste0("Metric_",info_kmean$Site_ID,"_",info_kmean$Image,".csv")){
  d=rbind(d,read.table(paste0("../Data/Metrics/",k),sep=";",header = T))}


## aggregating the metrics into a df

info_kmean=read.table("../Data/Landscapes/Kmean_clust_info.csv",sep=",",header = T)%>%
  filter(., Size!=200,Dataset=="biodesert",status=="kept") #keeping the kept sites

d=tibble()
for (k in paste0("Metric_",info_kmean$Site_ID,"_",info_kmean$Image,".csv")){
  d=rbind(d,read.table(paste0("../Data/Metrics/",k),sep=";",header = T))
}

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

ggsave("../Figures/SI/PCA_climatic_variables.pdf",ggarrange(p1,p2),
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
             Org_C_tot=as.numeric(sapply(1:nrow(.),function(x){return(d_biodesert$`TOC veg`[which(d_biodesert$ID==.$Site_ID[x])]
                                                                      -
                                                                        d_biodesert$`TOC b`[which(d_biodesert$ID==.$Site_ID[x])]
             )})),
             Org_C_tot_v=as.numeric(sapply(1:nrow(.),function(x){return(d_biodesert$`TOC veg`[which(d_biodesert$ID==.$Site_ID[x])]
                                                                      #-
                                                                      #  d_biodesert$`TOC b`[which(d_biodesert$ID==.$Site_ID[x])]
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
                                                                    #-
                                                                    #  d_biodesert$`AMO b`[which(d_biodesert$ID==.$Site_ID[x])]
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
             )})),
             lnNitrate=as.numeric(sapply(1:nrow(.),function(x){return(d_biodesert$`NIT veg`[which(d_biodesert$ID==.$Site_ID[x])]
                                                                     -
                                                                       d_biodesert$`NIT b`[which(d_biodesert$ID==.$Site_ID[x])]
             )})),
             lnAmonium=as.numeric(sapply(1:nrow(.),function(x){return(d_biodesert$`AMO veg`[which(d_biodesert$ID==.$Site_ID[x])]
                                                                    -
                                                                      d_biodesert$`AMO b`[which(d_biodesert$ID==.$Site_ID[x])]
             )}))
             
  )

#transforming into true Site ID
d$Site_ID=sapply(1:nrow(d),function(x){return(d_biodesert$SITE_ID[which(d_biodesert$ID==d$Site_ID[x])])})

write.table(d,"../Data/Spatial_structure_grazing_8_neigh.csv",sep=";")

# ---------------------- Step 2: Effects of grazing on the spatial structure  ----
## >> 1) Running mixed-effect models ---- 

Run_model_importance_aridity=function(id){
  
  list_mod=rbind(expand.grid(with_cover=c(T,F),
                             Stats=c("perim_area_scaling","fmax_psd","PL_expo",
                                     "core_area_land","core_area",
                                     "flow_length","PLR","KS_dist","Shape_metric",
                                     "mean_psd")),
                 tibble(with_cover=F,
                        Stats=c("rho_p")))
  
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
      + Aridity*Grazing 
      + rho_p + Type_veg 
      + Sand + Org_C_v + Grazing * Org_C_v
      + (1|Site_ID)")))
    
  } else{
    
    formula_mod=formula(formula_(paste("value ~ Long_cos + Long_sin + Lattitude + Slope + Elevation 
      + Aridity + Grazing 
      + Aridity*Grazing 
      + Type_veg 
      + Sand + Org_C_v + Grazing * Org_C_v
      + (1|Site_ID)")))
  }

  save=d_data_mod
  
  d_data_mod=save
  grazing_intensity="all"
  
  d_data_mod$Grazing=as.factor(d_data_mod$Grazing)
  d_data_mod$Grazing = relevel(d_data_mod$Grazing,ref = "0")
  
  model_spa_stat  = lmer(formula_mod, d_data_mod,
                         na.action = na.fail ,REML ="FALSE")
  
  #we remove potential outliers
  rm.outliers = romr.fnc(model_spa_stat, d_data_mod, trim=2.5)
  d_data_out = rm.outliers$data
  
  #saving data
  write.table(d_data_out,paste0("../Data/Linear_models_factor/Keep_data/Data_",stat,"_",with_cover,"_aridity_",grazing_intensity,".csv"))
  
  model_spa_stat = lmer(formula_mod, data = d_data_out, 
                        na.action = na.fail,REML ="FALSE")
  
  # #do some model selection
  select_model=dredge(model_spa_stat, subset = ~ Slope & Elevation &
                        Long_cos & Long_sin & Lattitude &
                        dc(Aridity & Grazing, Aridity : Grazing) &
                        dc(Org_C_v & Grazing, Org_C_v : Grazing),
                      extra=c(R2m=function(x) r.squaredGLMM(x)[1],
                              R2c=function(x) r.squaredGLMM(x)[2]),
                      options(na.action = "na.fail"))
  
  #extract the result of model selection
  if (dim(select_model%>%filter(., delta<2))[1]==1){ #one model
    
    result_select=model_spa_stat
    R2=select_model%>%filter(., AICc<min(AICc)+2)
    
    importance_mod=tibble(Interactions=1,
                          Sand =1,
                          Aridity=1,
                          Org_C=1,
                          Type_veg=1,
                          Grazing=1,
                          Cover=1,
                          Woody=1)
    
    d_all=rbind(d_all,cbind(tibble(Sp_stat=stat,
                                   N_outliers=rm.outliers$n.removed,
                                   R2m=mean(R2$R2m),R2C=mean(R2$R2c),
                                   With_cover=with_cover,
                                   Grazing_intensity=grazing_intensity),
                            importance_mod))
    
    
  }else{
    
    result_select=model.avg(select_model, subset = delta < 2)
    
    #Get the importance of each metric
    importance_mod=sw(result_select)
    R2=select_model%>%filter(., AICc<min(AICc)+2)
    
    d_all=rbind(d_all,cbind(tibble(Sp_stat=stat,
                                   N_outliers=rm.outliers$n.removed,
                                   R2m=mean(R2$R2m),R2C=mean(R2$R2c),
                                   With_cover=with_cover,
                                   Grazing_intensity=grazing_intensity),
                            Aggregate_importance(importance_mod,T)))
    
  }
  
  saveRDS(MuMIn::get.models(select_model,subset = delta<2)[[1]],
          paste0("../Data/Linear_models_factor/Keep_models/Mod_",stat,"_",
                 with_cover,"_aridity_",grazing_intensity,".rds"))
  
  summary_coef=confint(result_select)
  
  #Merge in a df
  d_all2=rbind(d_all2,tibble(Median=confint(result_select,level = 1e-9)[-1,1],
                             q1=summary_coef[-1,1], 
                             q3=summary_coef[-1,2],
                             # pvalue=summary(result_select)$coefmat.full[-1,5],
                             term=rownames(summary_coef)[-1],
                             Stat=stat,
                             R2m=mean(R2$R2m),
                             R2c=mean(R2$R2c),
                             Grazing_intensity=grazing_intensity)%>%
                 add_column(., With_cover=with_cover))
  
  
  
  #Saving each DF
  write.table(d_all,paste0("../Data/Linear_models_factor/Importance/Importance_",stat,"_",with_cover,"_aridity.csv"),sep=";")
  write.table(d_all2,paste0("../Data/Linear_models_factor/Estimator/Estimators_model_",stat,"_",with_cover,"_aridity.csv"),sep=";")
  
}

library(parallel)

mclapply(1:23,Run_model_importance_aridity,mc.cores = 23)

d_all=d_all2=tibble()
for (k in list.files("../Data/Linear_models_factor/Importance",paste0("aridity.csv"))){
  d=read.table(paste0("../Data/Linear_models_factor/Importance/",k),sep=";")
  if (nrow(d)!=0){
    d_all=rbind(d_all,d)
  }
}
for (k in list.files("../Data/Linear_models_factor/Estimator",paste0("aridity.csv"))){
  d2=read.table(paste0("../Data/Linear_models_factor/Estimator/",k),sep=";")
  if (nrow(d2)!=0){
    d_all2=rbind(d_all2,d2)
  }
}
write.table(d_all,paste0("../Data/Linear_models_factor/Importance_aridity_factor.csv"),sep=";")
write.table(d_all2,paste0("../Data/Linear_models_factor/Estimators_model_aridity_factor.csv"),sep=";")


Run_model_importance_aridity_no_inter=function(id){
  
  list_mod=rbind(expand.grid(with_cover=c(T,F),
                             Stats=c("perim_area_scaling","fmax_psd","PL_expo",
                                     "core_area_land","core_area",
                                     "flow_length","PLR","KS_dist","Shape_metric",
                                     "mean_psd")),
                 tibble(with_cover=F,
                        Stats=c("rho_p")))
  
  
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
      + rho_p + Type_veg 
      + Sand + Org_C_v 
      + (1|Site_ID)")))
    
  } else{
    
    formula_mod=formula(formula_(paste("value ~ Long_cos + Long_sin + Lattitude + Slope + Elevation 
      + Aridity + Grazing 
      + Type_veg
      + Sand + Org_C_v 
      + (1|Site_ID)")))
  }
  
  save=d_data_mod
  grazing_intensity='all'
  d_data_mod=save
  
  d_data_mod$Grazing=as.factor(d_data_mod$Grazing)
  d_data_mod$Grazing = relevel(d_data_mod$Grazing,ref = "0")
  
  model_spa_stat  = lmer(formula_mod, d_data_mod,
                         na.action = na.fail ,REML ="FALSE")
  
  #we remove potential outliers
  rm.outliers = romr.fnc(model_spa_stat, d_data_mod, trim=2.5)
  d_data_out = rm.outliers$data
  
  
  #saving data
  write.table(d_data_out,paste0("../Data/Linear_models_factor/Keep_data/Data_",stat,"_",with_cover,"_aridity_no_inter_",grazing_intensity,".csv"))
  
  model_spa_stat = lmer(formula_mod, data = d_data_out, 
                        na.action = na.fail,REML ="FALSE")
  
  # #do some model selection
  select_model=dredge(model_spa_stat, subset = ~ Slope & Elevation &
                        Long_cos & Long_sin & Lattitude &
                        dc(Woody & Grazing, Woody : Grazing) &
                        dc(Aridity & Grazing, Aridity : Grazing) &
                        dc(rho_p & Grazing, rho_p : Grazing) &
                        dc(Org_C_v & Grazing, Org_C_v : Grazing) &
                        dc(Sand  & Grazing, Sand  : Grazing) &
                        dc(Type_veg & Grazing, Type_veg : Grazing),
                      extra=c(R2m=function(x) r.squaredGLMM(x)[1],
                              R2c=function(x) r.squaredGLMM(x)[2]),
                      options(na.action = "na.fail") )
  
  if (dim(select_model%>%filter(., delta<2))[1]==1){ #one model
    
    result_select=model_spa_stat
    R2=select_model%>%filter(., AICc<min(AICc)+2)
    importance_mod=tibble(Interactions=1,
                          Sand =1,
                          Aridity=1,
                          Org_C=1,
                          Type_veg=1,
                          Grazing=1,
                          Cover=1,
                          Woody=1)
    
    d_all=rbind(d_all,cbind(tibble(Sp_stat=stat,
                                   N_outliers=rm.outliers$n.removed,
                                   R2m=mean(R2$R2m),R2C=mean(R2$R2c),
                                   With_cover=with_cover,
                                   Grazing_intensity=grazing_intensity),
                            importance_mod))
    
  }else{
    
    result_select=model.avg(select_model, subset = delta < 2)
    
    #Get the importance of each metric
    importance_mod=sw(result_select)
    R2=select_model%>%filter(., AICc<min(AICc)+2)
    
    d_all=rbind(d_all,cbind(tibble(Sp_stat=stat,
                                   N_outliers=rm.outliers$n.removed,
                                   R2m=mean(R2$R2m),R2C=mean(R2$R2c)),
                            Aggregate_importance(importance_mod,T))%>%
                  add_column(., With_cover=with_cover,
                             Grazing_intensity=grazing_intensity))
    
  }
  
  saveRDS(MuMIn::get.models(select_model,subset = delta<2)[[1]],
          paste0("../Data/Linear_models_factor/Keep_models/Mod_",stat,
                 "_",with_cover,"_aridity_no_inter_",grazing_intensity,".rds"))
  
  summary_coef=confint(result_select)
  
  
  #Merge in a df
  d_all2=rbind(d_all2,tibble(Median=confint(result_select,level = 1e-9)[-1,1],
                             q1=summary_coef[-1,1], 
                             q3=summary_coef[-1,2],
                             # pvalue=summary(result_select)$coefmat.full[-1,5],
                             term=rownames(summary_coef)[-1],
                             Stat=stat,
                             R2m=mean(R2$R2m),
                             R2c=mean(R2$R2c),
                             Grazing_intensity=grazing_intensity)%>%
                 add_column(., With_cover=with_cover))

  #Saving each DF
  write.table(d_all,paste0("../Data/Linear_models_factor/Importance/Importance_",stat,"_",with_cover,"_aridity_no_inter.csv"),sep=";")
  write.table(d_all2,paste0("../Data/Linear_models_factor/Estimator/Estimators_model_",stat,"_",with_cover,"_aridity_no_inter.csv"),sep=";")
  
}

library(parallel)

mclapply(1:23,Run_model_importance_aridity_no_inter,mc.cores = 23)

d_all=d_all2=tibble()
for (k in list.files("../Data/Linear_models_factor/Importance",paste0("aridity_no_inter"))){
  d=read.table(paste0("../Data/Linear_models_factor/Importance/",k),sep=";")
  if (nrow(d)!=0){
    d_all=rbind(d_all,d)
  }
}
for (k in list.files("../Data/Linear_models_factor/Estimator",paste0("aridity_no_inter"))){
  d2=read.table(paste0("../Data/Linear_models_factor/Estimator/",k),sep=";")
  if (nrow(d)!=0){
    d_all2=rbind(d_all2,d2)
  }
}
write.table(d_all,paste0("../Data/Linear_models_factor/Importance_aridity_no_inter_factor.csv"),sep=";")
write.table(d_all2,paste0("../Data/Linear_models_factor/Estimators_model_aridity_no_inter_factor.csv"),sep=";")


## >> 2) Analyzing the residuals ---- 



d_data=read.table("../Data/Spatial_structure_grazing.csv",sep=";")%>%
  Closer_to_normality(.)
d_data[,c(1,9:29,36,38:42,45:ncol(d_data))] = 
  apply(d_data[,c(1,9:29,36,38:42,45:ncol(d_data))],2,z_tranform)
d_data=Perform_PCA_spatial_struc(d_data)
d_data$Grazing=as.factor(d_data$Grazing)
save=d_data



boot_function_lm_graz = function(formula, data, indices) {
  d = data[indices,] 
  fit = lm(formula, data=d) 
  return(summary(fit)$coefficient[2:4,1])
}
boot_function_lm_arid = function(formula, data, indices) {
  d = data[indices,] 
  fit = lm(formula, data=d) 
  return(summary(fit)$coefficient[2,1])
}

grazing_intensity="all"
for (with_interaction in c(T,F)){
  
  d_slope=tibble()
  
  for (k in c("perim_area_scaling","fmax_psd","PL_expo",
              "core_area_land",#"core_area", since the mean % of core pixels in patches has no been selected for grazing, we remove it
              "flow_length","PLR","KS_dist","Shape_metric","mean_psd","rho_p")){
    
    #for each we plot the slope against the partial residuals with and without cover
    
    if (k !="rho_p"){
      
      d_data_out=read.table(paste0("../Data/Linear_models_factor/Keep_data/Data_",k,"_TRUE_aridity_",
                                   ifelse(with_interaction,"","no_inter_"),grazing_intensity,".csv"),sep=" ")
      model_spa_stat=readRDS(paste0("../Data/Linear_models_factor/Keep_models/Mod_",k,"_TRUE_aridity_",
                                    ifelse(with_interaction,"","no_inter_"),grazing_intensity,".rds"))
      
      d_data_out$Grazing = as.factor(d_data_out$Grazing)
      d_data_out$Grazing = relevel(d_data_out$Grazing,ref = "0")
      
      resid_mod=visreg::visreg(fit = model_spa_stat,xvar="Grazing",plot=F) 
      
      mod_cov = boot(data=resid_mod$res, statistic=boot_function_lm_graz, R=1000, formula=visregRes~Grazing)
      
      d_slope=rbind(d_slope,
                    tibble(pval=apply(mod_cov$t,2,twoside_pvalue),
                           q2=apply((mod_cov$t),2,median),
                           q1=apply(mod_cov$t,2,quantile,.025),
                           q3=apply(mod_cov$t,2,quantile,.975),
                           q1_90=apply(mod_cov$t,2,quantile,.05),
                           q3_90=apply(mod_cov$t,2,quantile,.95),
                           With_cover=T,Grazing_intensity=grazing_intensity,
                           ID_grazing=c("1","2","3"),
                           Stat=k,Driver="Grazing"
                    ))
      
    }
    
    d_data_out=read.table(paste0("../Data/Linear_models_factor/Keep_data/Data_",k,"_FALSE_aridity_",
                                 ifelse(with_interaction,"","no_inter_"),grazing_intensity,".csv"),sep=" ")
    model_spa_stat=readRDS(paste0("../Data/Linear_models_factor/Keep_models/Mod_",k,"_FALSE_aridity_",
                                  ifelse(with_interaction,"","no_inter_"),grazing_intensity,".rds"))
    
    d_data_out$Grazing = as.factor(d_data_out$Grazing)
    d_data_out$Grazing = relevel(d_data_out$Grazing,ref = "0")
    
    resid_mod=visreg::visreg(fit = model_spa_stat,xvar="Grazing",plot=F) 
    
    mod_cov = boot(data=resid_mod$res, statistic=boot_function_lm_graz, R=1000, formula=visregRes~Grazing)
    
    d_slope=rbind(d_slope,
                  tibble(pval=apply(mod_cov$t,2,twoside_pvalue),
                         q2=apply((mod_cov$t),2,median),
                         q1=apply(mod_cov$t,2,quantile,.025),
                         q3=apply(mod_cov$t,2,quantile,.975),
                         q1_90=apply(mod_cov$t,2,quantile,.05),
                         q3_90=apply(mod_cov$t,2,quantile,.95),
                         With_cover=F,Grazing_intensity=grazing_intensity,
                         ID_grazing=c("1","2","3"),
                         Stat=k,Driver="Grazing"
                  ))
    
    
  }
  write.table(d_slope,paste0("../Data/Linear_models_factor/Slope_partial_residuals_aridity",
                             ifelse(with_interaction,"","_no_inter"),
                             ".csv"),sep=";")
  
}





# ---------------------- Step 3: Island of fertility ----

## >> 1) Running mixed-effect models ---- 

Run_model_fertility_nutrients=function(id){
  
  list_mod=expand.grid(with_cover=c(T,F),
                       Stats=c("Org_C","Org_C_v","Total_N","lnTotal_N","Nitrate","lnNitrate"))
  
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
      + rho_p + Type_veg 
      + Sand 
      + (1|Site_ID)")))
    
  } else{
    
    formula_mod=formula(formula_(paste("value ~ Long_cos + Long_sin + Lattitude + Slope + Elevation 
      + Aridity + Grazing 
      + Type_veg 
      + Sand 
      + (1|Site_ID)")))
  }
  
  
  save=d_data_mod
  grazing_intensity="all"
  
  d_data_mod$Grazing=as.factor(d_data_mod$Grazing)
  d_data_mod$Grazing = relevel(d_data_mod$Grazing,ref = "0")
  
  model_spa_stat  = lmer(formula_mod, d_data_mod,
                         na.action = na.fail ,REML ="FALSE")
  
  #we remove potential outliers
  rm.outliers = romr.fnc(model_spa_stat, d_data_mod, trim=2.5)
  d_data_out = rm.outliers$data
  d_data_out$Grazing=as.factor(d_data_out$Grazing)
  d_data_out$Grazing = relevel(d_data_out$Grazing,ref = "0")
  
  #saving data
  write.table(d_data_out,paste0("../Data/Linear_models_fertility_factor/Keep_data/Data_",stat,"_",with_cover,"_aridity_",grazing_intensity,".csv"))
  
  model_spa_stat = lmer(formula_mod, data = d_data_out, 
                        na.action = na.fail,REML ="FALSE")
  
  # #do some model selection
  select_model=dredge(model_spa_stat, subset = ~ Slope & Elevation &
                        Long_cos & Long_sin & Lattitude &
                        dc(Aridity & Grazing, Aridity : Grazing),
                      extra=c(R2m=function(x) r.squaredGLMM(x)[1],
                              R2c=function(x) r.squaredGLMM(x)[2]),
                      options(na.action = "na.fail") )
  
  #extract the result of model selection
  if (dim(select_model%>%filter(., delta<2))[1]==1){ #one model
    
    result_select=model_spa_stat
    R2=select_model%>%filter(., AICc<min(AICc)+2)
    
    importance_mod=tibble(Interactions=1,
                          Sand =1,
                          Aridity=1,
                          Org_C=1,
                          Type_veg=1,
                          Grazing=1,
                          Cover=1,
                          Woody=1)
    
    d_all=rbind(d_all,cbind(tibble(Sp_stat=stat,
                                   N_outliers=rm.outliers$n.removed,
                                   R2m=mean(R2$R2m),R2C=mean(R2$R2c),
                                   With_cover=with_cover,
                                   Grazing_intensity=grazing_intensity),
                            importance_mod))
    
  }else{
    
    result_select=model.avg(select_model, subset = delta < 2)
    
    #Get the importance of each metric
    importance_mod=sw(result_select)
    R2=select_model%>%filter(., AICc<min(AICc)+2)
    
    d_all=rbind(d_all,cbind(tibble(Sp_stat=stat,
                                   N_outliers=rm.outliers$n.removed,
                                   R2m=mean(R2$R2m),R2C=mean(R2$R2c),
                                   With_cover=with_cover,
                                   Grazing_intensity=grazing_intensity),
                            Aggregate_importance(importance_mod,T)))
    
  }
  
  saveRDS(model_spa_stat,paste0("../Data/Linear_models_fertility_factor/Keep_models/Mod_",stat,"_",with_cover,"_aridity_",grazing_intensity,".rds"))
  
  summary_coef=confint(result_select)
  
  #Merge in a df
  d_all2=rbind(d_all2,tibble(Median=confint(result_select,level = 1e-9)[-1,1],
                             q1=summary_coef[-1,1], 
                             q3=summary_coef[-1,2],
                             # pvalue=summary(result_select)$coefmat.full[-1,5],
                             term=rownames(summary_coef)[-1],
                             Stat=stat,
                             R2m=mean(R2$R2m),
                             R2c=mean(R2$R2c),
                             Grazing_intensity=grazing_intensity)%>%
                 add_column(., With_cover=with_cover))
  
  
  
  #Saving each DF
  write.table(d_all,paste0("../Data/Linear_models_fertility_factor/Importance/Importance_",stat,"_",with_cover,"_aridity.csv"),sep=";")
  write.table(d_all2,paste0("../Data/Linear_models_fertility_factor/Estimator/Estimators_model_",stat,"_",with_cover,"_aridity.csv"),sep=";")
  
}

library(parallel)

mclapply(1:12,Run_model_fertility_nutrients,mc.cores = 1)

d_all=d_all2=tibble()
for (k in list.files("../Data/Linear_models_fertility_factor/Importance",paste0("aridity.csv"))){
  d=read.table(paste0("../Data/Linear_models_fertility_factor/Importance/",k),sep=";")
  if (nrow(d)!=0){
    d_all=rbind(d_all,d)
  }
}
for (k in list.files("../Data/Linear_models_fertility_factor/Estimator",paste0("aridity.csv"))){
  d2=read.table(paste0("../Data/Linear_models_fertility_factor/Estimator/",k),sep=";")
  if (nrow(d2)!=0){
    d_all2=rbind(d_all2,d2)
  }
}
write.table(d_all,paste0("../Data/Linear_models_fertility_factor/Importance_nutrients.csv"),sep=";")
write.table(d_all2,paste0("../Data/Linear_models_fertility_factor/Estimators_model_nutrients.csv"),sep=";")


## >> 2) Analyzing the residuals ---- 


boot_function_lm_graz = function(formula, data, indices) {
  d = data[indices,] 
  d$Grazing=as.factor(d$Grazing)
  fit = lm(formula, data=d) 
  return(summary(fit)$coefficient[2:4,1])
}
boot_function_lm_arid = function(formula, data, indices) {
  d = data[indices,] 
  fit = lm(formula, data=d) 
  return(summary(fit)$coefficient[2,1])
}


d_data=read.table("../Data/Spatial_structure_grazing.csv",sep=";")%>%
  Closer_to_normality(.)
d_data[,c(1,9:29,36,38:42,45:ncol(d_data))] = 
  apply(d_data[,c(1,9:29,36,38:42,45:ncol(d_data))],2,z_tranform)
d_data=Perform_PCA_spatial_struc(d_data)
save=d_data
d_slope=tibble()

for (k in c("Org_C","Org_C_v","Total_N","lnTotal_N",
            "Total_P","lnTotal_P","Nitrate","lnNitrate")){
  
  grazing_intensity="all"
  #for each we plot the slope against the partial residuals with and without cover
  
  d_data_out=read.table(paste0("../Data/Linear_models_fertility_factor/Keep_data/Data_",k,
                               "_TRUE_aridity_",grazing_intensity,".csv"),sep=" ")
  model_spa_stat=readRDS(paste0("../Data/Linear_models_fertility_factor/Keep_models/Mod_",k,
                                "_TRUE_aridity_",grazing_intensity,".rds"))
  
  d_data_out$Grazing = as.factor(d_data_out$Grazing)
  d_data_out$Grazing = relevel(d_data_out$Grazing,ref = "0")
  
  resid_mod=visreg::visreg(fit = model_spa_stat,xvar="Grazing",plot=F) 
  
  mod_cov = boot(data=resid_mod$res, statistic=boot_function_lm_graz, R=1000, formula=visregRes~Grazing)
  
  
  d_slope=rbind(d_slope,
                tibble(pval=apply(mod_cov$t,2,twoside_pvalue),
                       q2=apply((mod_cov$t),2,median),
                       q1=apply(mod_cov$t,2,quantile,.025),
                       q3=apply(mod_cov$t,2,quantile,.975),
                       q1_90=apply(mod_cov$t,2,quantile,.05),
                       q3_90=apply(mod_cov$t,2,quantile,.95),
                       With_cover=T,Grazing_intensity=grazing_intensity,
                       ID_grazing=c("1","2","3"),
                       Stat=k,Driver="Grazing"
                ))
  
  
  resid_mod=visreg::visreg(fit = model_spa_stat,xvar="Aridity",plot=F) 
  
  mod_cov = boot(data=resid_mod$res, statistic=boot_function_lm_arid, R=1000, formula=visregRes~Aridity)
  
  
  d_slope=rbind(d_slope,
                tibble(pval=apply(mod_cov$t,2,twoside_pvalue),
                       q2=apply((mod_cov$t),2,median),
                       q1=apply(mod_cov$t,2,quantile,.025),
                       q3=apply(mod_cov$t,2,quantile,.975),
                       q1_90=apply(mod_cov$t,2,quantile,.05),
                       q3_90=apply(mod_cov$t,2,quantile,.95),
                       With_cover=T,Grazing_intensity=grazing_intensity,
                       ID_grazing="none",
                       Stat=k,Driver="Aridity"
                ))
  
  
  
  
  grazing_intensity="all"
  #for each we plot the slope against the partial residuals with and without cover
  
  d_data_out=read.table(paste0("../Data/Linear_models_fertility_factor/Keep_data/Data_",k,
                               "_FALSE_aridity_",grazing_intensity,".csv"),sep=" ")
  model_spa_stat=readRDS(paste0("../Data/Linear_models_fertility_factor/Keep_models/Mod_",k,
                                "_FALSE_aridity_",grazing_intensity,".rds"))
  
  d_data_out$Grazing = as.factor(d_data_out$Grazing)
  d_data_out$Grazing = relevel(d_data_out$Grazing,ref = "0")
  
  resid_mod=visreg::visreg(fit = model_spa_stat,xvar="Grazing",plot=F) 
  
  mod_cov = boot(data=resid_mod$res, statistic=boot_function_lm_graz, R=1000, formula=visregRes~Grazing)
  
  
  d_slope=rbind(d_slope,
                tibble(pval=apply(mod_cov$t,2,twoside_pvalue),
                       q2=apply((mod_cov$t),2,median),
                       q1=apply(mod_cov$t,2,quantile,.025),
                       q3=apply(mod_cov$t,2,quantile,.975),
                       q1_90=apply(mod_cov$t,2,quantile,.05),
                       q3_90=apply(mod_cov$t,2,quantile,.95),
                       With_cover=F,Grazing_intensity=grazing_intensity,
                       ID_grazing=c("1","2","3"),
                       Stat=k,Driver="Grazing"
                ))
  
  
  resid_mod=visreg::visreg(fit = model_spa_stat,xvar="Aridity",plot=F) 
  
  mod_cov = boot(data=resid_mod$res, statistic=boot_function_lm_arid, R=1000, formula=visregRes~Aridity)
  
  
  d_slope=rbind(d_slope,
                tibble(pval=apply(mod_cov$t,2,twoside_pvalue),
                       q2=apply((mod_cov$t),2,median),
                       q1=apply(mod_cov$t,2,quantile,.025),
                       q3=apply(mod_cov$t,2,quantile,.975),
                       q1_90=apply(mod_cov$t,2,quantile,.05),
                       q3_90=apply(mod_cov$t,2,quantile,.95),
                       With_cover=F,Grazing_intensity=grazing_intensity,
                       ID_grazing="none",
                       Stat=k,Driver="Aridity"
                ))
}

write.table(d_slope,"../Data/Linear_models_fertility_factor/Slope_partial_residuals_aridity_grazing.csv",sep=";")




# ---------------------- Step 4: SEM ----

for (k in c("fmax_psd","PL_expo","perim_area_scaling","core_area_land","core_area")){
  
  save=Get_data_resid_SEM(k,"grazed")
  
  d_sem=save#%>%filter(., Grazing %in% c(0,1))
  SEM_distance=summary((psem(
    lmer(Resid_mod ~  (1|Site_ID) + Org_C_v+Grazing+Sand, d_sem),
    lmer(Nitrate ~  (1|Site_ID) + Grazing +Aridity+Sand, d_sem),
    lmer(Org_C_v ~  (1|Site_ID) + Nitrate  + Aridity+Sand+Grazing, d_sem),
    Sand%~~%Aridity
  )))
  
  Plot_SEM(SEM_distance,pdf_ = F,title_ = paste0("",
                                                 "\n Goodness of fit: Fisher C stat = ",
                                                 round(SEM_distance$Cstat,2)[1],
                                                 ",df = ",SEM_distance$Cstat[2],
                                                 ", Pval = ",SEM_distance$Cstat[3]),
           name_var = "% landscape   \n covered by core    ",type_N = "Nitrate",
           label_cex = 1.1,edge_cex = 1.2)
  
  
}




# ---------------------- Step 5: Moving average aridity ----


Run_model_moving_average_aridity=function(id){
  
  list_mod=expand.grid(with_cover=c(T,F),
                       Stats=c("perim_area_scaling","fmax_psd","PL_expo",
                               "core_area_land","core_area",
                               "flow_length","PLR","KS_dist"))
  
  d_slope=tibble()#for keeping all information
  
  d_data=read.table("../Data/Spatial_structure_grazing.csv",sep=";")%>%
    Closer_to_normality(.)
  
  d_data[,c(1,9:29,38:42,45:ncol(d_data))] = 
    apply(d_data[,c(1,9:29,38:42,45:ncol(d_data))],2,z_tranform)
  d_data$Grazing=as.factor(d_data$Grazing)
  d_data$Grazing = relevel(d_data$Grazing,ref = "0")
  
  d_data=Perform_PCA_spatial_struc(d_data)
  
  stat=list_mod$Stats[id]
  
  d_data_mod=d_data%>%melt(., measure.vars=stat)%>%
    filter(., !is.na(value))
  
  with_cover=list_mod$with_cover[id]
  
  if (with_cover){
    
    formula_mod=formula(formula_(paste("value ~ Long_cos + Long_sin + Lattitude + Slope + Elevation 
      + Aridity + Grazing 
      + rho_p + Type_veg 
      + Sand + Org_C_v "))) #removing site effect to avoid singularity issues due to non-complete dataset
    
  } else{
    
    formula_mod=formula(formula_(paste("value ~ Long_cos + Long_sin + Lattitude + Slope + Elevation 
      + Aridity + Grazing 
      + Type_veg
      + Sand + Org_C_v "))) #removing site effect to avoid singularity issues due to non-complete dataset
  }
  
  save=d_data_mod
  grazing_intensity="all"
  
  for (arid_thresh in seq(.73,max(save$Aridity),length.out=100)){  #keeping at least 90 sites
    
    d_data_mod=save%>%filter(., Aridity<arid_thresh)
    d_data_mod$Aridity=scale(d_data_mod$Aridity)[,1]
    
    model_spa_stat  = lm(formula_mod, d_data_mod,
                         na.action = na.fail )
    
    resid_mod=visreg::visreg(fit = model_spa_stat,xvar="Aridity",plot=F) 
    
    mod_cov = boot(data=resid_mod$res, statistic=boot_function_lm, R=1000, formula=visregRes~Aridity)
    
    #Merge in a df
    d_slope=rbind(d_slope,
                  tibble(pval=twoside_pvalue(mod_cov$t),
                         q2=median(mod_cov$t),
                         q1=quantile(mod_cov$t,.025),
                         q3=quantile(mod_cov$t,.975),
                         With_cover=with_cover,Aridity_thresh=arid_thresh,
                         Stat=stat,Driver="Aridity"
                  ))
  }
  
  write.table(d_slope,paste0("../Data/Moving_average_aridity/Moving_average_aridity_",stat,"_",with_cover,".csv"),sep=";")
  print("--")
}

library(parallel)

mclapply(1:18,Run_model_moving_average_aridity,mc.cores = 18)


d_slope=tibble()
for (k in list.files("../Data/Moving_average_aridity",paste0("Moving"))){
  d=read.table(paste0("../Data/Moving_average_aridity/",k),sep=";")
  d_slope=rbind(d_slope,d)
}
write.table(d_slope,paste0("../Data/Moving_average_aridity/Moving_average_aridity.csv"),sep=";")







# ---------------------- Step 6: Effect on the resilience ----
## >> 1) Inference parameters ----

NA_kept=100;id_plot=1    
d_biodesert=read.table("../Data/Spatial_structure_grazing.csv",sep=";")[,1:11]%>%
  add_column(., ID_sim=1:nrow(.))

d_sim=read.table("../Data/Simulations.csv",sep=";")%>%
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
  
  d_param_infer_rej[,empirical_id,1]=cross_valid$unadj.values[,1] # we keep the whole distribution for p
  d_param_infer_rej[,empirical_id,2]=cross_valid$unadj.values[,2] # for q
  d_param_infer_rej[,empirical_id,3]=cross_valid$unadj.values[,3] # for the scale of observation
}

write.table(d_NRMSE_sumstat,paste0("../Data/Inferrence/NRMSE_sumstat_all_",NA_kept,".csv"),sep=";")
write.table(x_y_stat,paste0("../Data/Inferrence/x_y_stat_all",NA_kept,".csv"),sep=";")
write.table(d_param_infer_rej,paste0("../Data/Inferrence/param_inferred.csv"),sep=";")

## >> 2) Selecting sites ----

library(diptest)
post_param=read.table("../Data/Inferrence/param_inferred.csv",sep=";")

d_data=read.table("../Data/Spatial_structure_grazing.csv",sep=";")
d_data$bimod=sapply(1:nrow(d_data),function(x){
  if (dip.test(post_param[,x])$p.value<.05 | dip.test(post_param[,x+504])$p.value<.05){
    return("bimod")
  }else {
    return("unimod")
  }
})
write.table(which(d_data$bimod!="bimod"),"../Data/Inferrence/Keeping_sites.csv",sep=";")




## >> 3) Computing resilience metrics ----

# See ABC_stability.jl file

## >> 4) Post-processing resilience metrics ----

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

write.table(d,"../Data/Inferrence/Prediction/Resilience_metrics.csv",sep=";")


## >> 5) Running mixed-effect models ----

d=read.table(paste0("../Data/Inferrence/param_inferred.csv"),sep=";",header=T)
d_data=read.table(paste0("../Data/Spatial_structure_grazing.csv"),sep=";")%>%
  Closer_to_normality(.)
n_sites=504

keep_sites=read.table(paste0("../Data/Inferrence/Keeping_sites.csv"),sep=";")$x

d=tibble(Site=1:n_sites,mean_p=apply(d[,(1:n_sites)],2,mean),sd_p=apply(d[,(1:n_sites)],2,sd),
         mean_q=apply(d[,(n_sites+1):(2*n_sites)],2,mean),sd_q=apply(d[,(n_sites+1):(2*n_sites)],2,sd),
         median_q=apply(d[,(n_sites+1):(2*n_sites)],2,median),
         Plot_n=d_data$Site_ID,
         Aridity=d_data$Aridity,
         Sp_richness=d_data$Sp_richness,
         Type_veg=d_data$Type_veg,
         Org_C=d_data$Org_C,
         Org_C_v=d_data$Org_C_v,
         Type_veg=d$Type_veg,
         Nitrate=(d_data$Nitrate),
         Amonium=(d_data$Amonium),
         Total_N=(d_data$Total_N),
         Total_P=(d_data$Total_P),
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

d2=read.table(paste0("../Data/Inferrence/Resilience_metrics.csv"),sep=";")%>%
  group_by(., Site)%>%
  dplyr::summarise(., .groups = "keep",
                   abs_mean=mean(pinfer-pcrit,na.rm = T),
                   abs_median=median(pinfer-pcrit,na.rm = T),
                   abs_sd=sd(pinfer-pcrit,na.rm = T),
                   relativ_mean=mean((pinfer-pcrit)/pinfer,na.rm = T),
                   relativ_sd=sd((pinfer-pcrit)/pinfer,na.rm = T),
                   relativ_median=median((pinfer-pcrit)/pinfer,na.rm = T),
                   Size_mean=mean(Size_tipping,na.rm = T),
                   Size_sd=sd(Size_tipping,na.rm = T))%>%
  filter(., Site %in% keep_sites)

d=cbind(d%>%filter(., Site %in% d2$Site),d2)

d2=tibble(q=scale(logit(d$median_q))[,1],
          abs_dist=scale(log(d$abs_median))[,1], #for normality purposes
          rela_dist=scale(d$relativ_median)[,1],
          
          #site related variables
          Sand=(d$Sand-mean(d$Sand,na.rm=T))/sd(d$Sand,na.rm = T),
          Sp_richness=(d$Sp_richness-mean(d$Sp_richness,na.rm=T))/sd(d$Sp_richness,na.rm = T),
          Grazing=as.factor(d$Grazing),
          Site=d$Site,
          Plot_n=d$Plot_n,
          Org_C=(d$Org_C-mean(d$Org_C,na.rm=T))/sd(d$Org_C,na.rm = T),
          Org_C_v=(d$Org_C_v-mean(d$Org_C_v,na.rm=T))/sd(d$Org_C_v,na.rm = T),
          Nitrate=(d$Nitrate-mean(d$Nitrate,na.rm=T))/sd(d$Nitrate,na.rm=T),
          Amonium=(d$Amonium-mean(d$Amonium,na.rm=T))/sd(d$Amonium,na.rm=T),
          Total_N=(d$Total_N-mean(d$Total_N,na.rm=T))/sd(d$Total_N,na.rm=T),
          Total_P=(d$Total_P-mean(d$Total_P,na.rm=T))/sd(d$Total_P,na.rm=T),
          Cover=(d$Cover-mean(d$Cover,na.rm=T))/sd(d$Cover,na.rm = T),
          Woody=(d$Woody-mean(d$Woody,na.rm=T))/sd(d$Woody,na.rm = T),
          
          Aridity=(d$Aridity-mean(d$Aridity,na.rm=T))/sd(d$Aridity,na.rm = T),
          
          #covariates
          Type_veg=d$Type_veg,
          Lat=(d$Lattitude -mean(d$Lattitude ,na.rm=T))/sd(d$Lattitude ,na.rm = T),
          Long_cos=(d$Long_cos-mean(d$Long_cos,na.rm=T))/sd(d$Long_cos,na.rm = T),
          Long_sin=(d$Long_sin-mean(d$Long_sin,na.rm=T))/sd(d$Long_sin,na.rm = T),
          Elevation=(d$Elevation-mean(d$Elevation,na.rm=T))/sd(d$Elevation,na.rm = T),
          Slope=(d$Slope-mean(d$Slope,na.rm=T))/sd(d$Slope,na.rm = T))

d2$Grazing=relevel(d2$Grazing,ref = "0")

d_partial=d_mod=tibble()
Sub=d2
Grazing_intensity="Full"

for (with_interaction in c(F)){
  
  if (with_interaction){
    mod_predictors=gsub("\n     ","","Aridity+Aridity*Grazing + Grazing + Sand + Sp_richness+ Org_C_v*Grazing + Type_veg + 
          Lat + Long_cos + Long_sin + Slope + Elevation + ( 1 | Plot_n)")
  }else{
    mod_predictors=gsub("\n     ","","Aridity + Grazing + Sand + Sp_richness+ Org_C_v + Type_veg +
          Lat + Long_cos + Long_sin + Slope + Elevation + ( 1 | Plot_n)")
  }
  
  model_q=lmer(formula = paste("q ~ ",mod_predictors), data = Sub,na.action = na.fail ,REML ="FALSE")
  model_abs=lmer(formula = paste("abs_dist ~ ",mod_predictors), data = Sub,na.action = na.fail ,REML ="FALSE")
  model_rela=lmer(formula = paste("rela_dist ~ ",mod_predictors), data = Sub,na.action = na.fail ,REML ="FALSE")
  model_cover=lmer(formula = paste("Cover ~ ",mod_predictors), data = Sub,na.action = na.fail ,REML ="FALSE")
  
  #Getting partial prediction
  
  boot_function_lm_graz = function(formula, data, indices) {
    d = data[indices,] 
    fit = lm(formula, data=d) 
    return(summary(fit)$coefficient[2:4,1])
  }
  boot_function_lm_arid = function(formula, data, indices) {
    d = data[indices,] 
    fit = lm(formula, data=d) 
    return(summary(fit)$coefficient[2,1])
  }
  
  #GRAZING
  
  resid_mod_abs=visreg::visreg(fit = model_abs,xvar="Grazing",plot=F) 
  model_abs_boot=boot(data=resid_mod_abs$res,statistic=boot_function_lm_graz,R=500, formula=visregRes~Grazing)$t[,1:3]
  
  resid_mod_rela=visreg::visreg(fit = model_rela,xvar="Grazing",plot=F) 
  model_rela_boot=boot(data=resid_mod_rela$res,statistic=boot_function_lm_graz,R=500, formula=visregRes~Grazing)$t[,1:3]
  
  resid_mod_q=visreg::visreg(fit = model_q,xvar="Grazing",plot=F) 
  model_q_boot=boot(data=resid_mod_q$res,statistic=boot_function_lm_graz,R=500, formula=visregRes~Grazing)$t[,1:3]
  
  resid_mod_cover=visreg::visreg(fit = model_cover,xvar="Grazing",plot=F) 
  model_cover_boot=boot(data=resid_mod_cover$res,statistic=boot_function_lm_graz,R=500, formula=visregRes~Grazing)$t[,1:3]
  
  name_grazing_id=c("1","2","3")
  for (k in 1:3){
    d_mod=rbind(d_mod,data.frame(Stat=model_abs_boot[,k])%>%
                  add_column(., Param="Absolute distance",Type=Grazing_intensity,Metric="Grazing",Grazing_id=name_grazing_id[k],
                             With_inter=with_interaction))
    d_mod=rbind(d_mod,data.frame(Stat=model_rela_boot[,k])%>%
                  add_column(., Param="Relative distance",Type=Grazing_intensity,Metric="Grazing",Grazing_id=name_grazing_id[k],
                             With_inter=with_interaction))
    d_mod=rbind(d_mod,data.frame(Stat=model_q_boot[,k])%>%
                  add_column(., Param="q",Type=Grazing_intensity,Metric="Grazing",Grazing_id=name_grazing_id[k],
                             With_inter=with_interaction))
    d_mod=rbind(d_mod,data.frame(Stat=model_cover_boot[,k])%>%
                  add_column(., Param="Cover",Type=Grazing_intensity,Metric="Grazing",Grazing_id=name_grazing_id[k],
                             With_inter=with_interaction))
  }
  
  
  #ARIDITY
  resid_mod_abs=visreg::visreg(fit = model_abs,xvar="Aridity",plot=F) 
  model_abs_boot=boot(data=resid_mod_abs$res,statistic=boot_function_lm_arid,R=500, formula=visregRes~Aridity)$t[,1]
  
  resid_mod_rela=visreg::visreg(fit = model_rela,xvar="Aridity",plot=F) 
  model_rela_boot=boot(data=resid_mod_rela$res,statistic=boot_function_lm_arid,R=500, formula=visregRes~Aridity)$t[,1]
  
  resid_mod_q=visreg::visreg(fit = model_q,xvar="Aridity",plot=F) 
  model_q_boot=boot(data=resid_mod_q$res,statistic=boot_function_lm_arid,R=500, formula=visregRes~Aridity)$t[,1]
  
  resid_mod_cover=visreg::visreg(fit = model_cover,xvar="Aridity",plot=F) 
  model_cover_boot=boot(data=resid_mod_cover$res,statistic=boot_function_lm_arid,R=500, formula=visregRes~Aridity)$t[,1]
  
  d_mod=rbind(d_mod,data.frame(Stat=model_abs_boot)%>%
                add_column(., Param="Absolute distance",Type=Grazing_intensity,Metric="Aridity",Grazing_id="none",
                           With_inter=with_interaction))
  d_mod=rbind(d_mod,data.frame(Stat=model_rela_boot)%>%
                add_column(., Param="Relative distance",Type=Grazing_intensity,Metric="Aridity",Grazing_id="none",
                           With_inter=with_interaction))
  d_mod=rbind(d_mod,data.frame(Stat=model_q_boot)%>%
                add_column(., Param="q",Type=Grazing_intensity,Metric="Aridity",Grazing_id="none",
                           With_inter=with_interaction))
  d_mod=rbind(d_mod,data.frame(Stat=model_cover_boot)%>%
                add_column(., Param="Cover",Type=Grazing_intensity,Metric="Aridity",Grazing_id="none",
                           With_inter=with_interaction))

}
write.table(d_mod,paste0("../Data/Inferrence/Grazing_on_resilience_factor.csv"),sep=";")

## >> 6) SEM on the resilience of drylands ----

d=read.table(paste0("../Data/Inferrence/param_inferred.csv"),sep=";",header=T)
d_data=read.table(paste0("../Data/Spatial_structure_grazing.csv"),sep=";")%>%
  Closer_to_normality(.)
n_sites=504
keep_sites=read.table(paste0("../Data/Inferrence/Keeping_sites.csv"),sep=";")$x

d=tibble(Site=1:n_sites,mean_p=apply(d[,(1:n_sites)],2,mean),sd_p=apply(d[,(1:n_sites)],2,sd),
         mean_q=apply(d[,(n_sites+1):(2*n_sites)],2,mean),sd_q=apply(d[,(n_sites+1):(2*n_sites)],2,sd),
         median_q=apply(d[,(n_sites+1):(2*n_sites)],2,median),
         Site_ID=d_data$Site_ID,
         Aridity=d_data$Aridity,
         Sp_richness=d_data$Sp_richness,
         Type_veg=d_data$Type_veg,
         Org_C=d_data$Org_C,
         Org_C_v=d_data$Org_C_v,
         Org_C_tot=d_data$Org_C_tot,
         Nitrate=d_data$Nitrate,
         Amonium=d_data$Amonium,
         Total_N=d_data$Total_N,
         Total_P=d_data$Total_P,
         lnNitrate=d_data$lnNitrate,
         lnAmonium=d_data$lnAmonium,
         lnTotal_N=d_data$lnTotal_N,
         lnTotal_P=d_data$lnTotal_P,
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

d2=read.table(paste0("../Data/Inferrence/Resilience_metrics.csv"),sep=";")%>%
  group_by(., Site)%>%
  dplyr::summarise(., .groups = "keep",
                   abs_mean=mean(pinfer-pcrit,na.rm = T),
                   abs_median=median(pinfer-pcrit,na.rm = T),
                   abs_sd=sd(pinfer-pcrit,na.rm = T),
                   relativ_mean=mean((pinfer-pcrit)/pinfer,na.rm = T),
                   relativ_sd=sd((pinfer-pcrit)/pinfer,na.rm = T),
                   relativ_median=median((pinfer-pcrit)/pinfer,na.rm = T),
                   Size_mean=mean(Size_tipping,na.rm = T),
                   Size_sd=sd(Size_tipping,na.rm = T))%>%
  filter(., Site %in% keep_sites)

d=cbind(d%>%filter(., Site %in% d2$Site),d2)

d2=tibble(q=scale(logit(d$median_q))[,1],
          abs_dist=scale(log(d$abs_median))[,1],
          rela_dist=scale(log(d$relativ_median))[,1],
          
          #site related variables
          Sand=(d$Sand-mean(d$Sand,na.rm=T))/sd(d$Sand,na.rm = T),
          Sp_richness=(d$Sp_richness-mean(d$Sp_richness,na.rm=T))/sd(d$Sp_richness,na.rm = T),
          Grazing=d$Grazing,
          Site=d$Site,
          Site_ID=d$Site_ID,
          Org_C=(d$Org_C-mean(d$Org_C,na.rm=T))/sd(d$Org_C,na.rm=T),
          Org_C_v=(d$Org_C_v-mean(d$Org_C_v,na.rm=T))/sd(d$Org_C_v,na.rm=T),
          Org_C_tot=(d$Org_C_tot-mean(d$Org_C_tot,na.rm=T))/sd(d$Org_C_tot,na.rm=T),
          Nitrate=(d$Nitrate-mean(d$Nitrate,na.rm=T))/sd(d$Nitrate,na.rm=T),
          Amonium=(d$Amonium-mean(d$Amonium,na.rm=T))/sd(d$Amonium,na.rm=T),
          Total_N=(d$Total_N-mean(d$Total_N,na.rm=T))/sd(d$Total_N,na.rm=T),
          Total_P=(d$Total_P-mean(d$Total_P,na.rm=T))/sd(d$Total_P,na.rm=T),
          lnTotal_N=(d$lnTotal_N-mean(d$lnTotal_N,na.rm=T))/sd(d$lnTotal_N,na.rm=T),
          lnTotal_P=(d$lnTotal_P-mean(d$lnTotal_P,na.rm=T))/sd(d$lnTotal_P,na.rm=T),
          lnNitrate=(d$lnNitrate-mean(d$lnNitrate,na.rm=T))/sd(d$lnNitrate,na.rm=T),
          lnAmonium=(d$lnAmonium-mean(d$lnAmonium,na.rm=T))/sd(d$lnAmonium,na.rm=T),
          Cover=(d$Cover-mean(d$Cover,na.rm=T))/sd(d$Cover,na.rm = T),
          Woody=(d$Woody-mean(d$Woody,na.rm=T))/sd(d$Woody,na.rm = T),
          Aridity=(d$Aridity-mean(d$Aridity,na.rm=T))/sd(d$Aridity,na.rm = T),
          
          #covariates
          Type_veg=d$Type_veg,
          Lat=(d$Lattitude -mean(d$Lattitude ,na.rm=T))/sd(d$Lattitude ,na.rm = T),
          Long_cos=(d$Long_cos-mean(d$Long_cos,na.rm=T))/sd(d$Long_cos,na.rm = T),
          Long_sin=(d$Long_sin-mean(d$Long_sin,na.rm=T))/sd(d$Long_sin,na.rm = T),
          Elevation=(d$Elevation-mean(d$Elevation,na.rm=T))/sd(d$Elevation,na.rm = T),
          Slope=(d$Slope-mean(d$Slope,na.rm=T))/sd(d$Slope,na.rm = T))



k="abs_dist"

#We first control for all covariates and extract the residuals
model_lmer=lm("value ~ Long_cos + Long_sin + Lat + Elevation",
              data = d2%>%melt(., measure.vars=k)%>%filter(., !is.na(value)),
              na.action = na.fail)

resid_model=residuals(model_lmer) #extract residuals of the distance to the tipping point after controlling for all covariates

save=d2[as.numeric(names(resid_model)),]%>% #add it to the dataframe 
  add_column(., Resid_mod=resid_model)

boxplot(save$abs_dist~save$Grazing) #linear changes of resilience along grazing gradient
boxplot(save$q~save$Grazing) #linear changes of resilience along grazing gradient
boxplot(save$Cover~save$Grazing) #linear changes of resilience along grazing gradient
d_sem=save

SEM_distance=psem(
  lm(Resid_mod ~ Cover + q, d_sem),
  lmer(Cover ~ (1|Site_ID)   + Org_C_v + Grazing+Nitrate, d_sem),
  lmer(q ~ (1|Site_ID) +   Nitrate + Org_C_v + Grazing, d_sem),
  lmer(Nitrate ~ (1|Site_ID)   + Grazing, d_sem),
  lmer(Org_C_v ~ (1|Site_ID)   + Grazing+Nitrate, d_sem),
  q%~~%Cover
)
saveRDS(SEM_distance,"../Data/SEM/SEM_resilience.rds")

SEM_distance_boot = bootEff(SEM_distance, R = 1000, seed = 13, parallel = "snow",ran.eff = "Site_ID",ncpus = 50)
saveRDS(SEM_distance_boot,"../Data/SEM/SEM_boot_resilience.rds")



#bootstrapped direct and total effects
SEM_resilience=semEff(readRDS("../Data/SEM/SEM_boot_resilience.rds"))

Direct_effects=tibble()
for (predictors in c("Resid.mod","q","Org.C.v","Cover","Nitrate")){
  effect_pred=SEM_resilience$`Bootstrapped Effects`[[predictors]]
  Direct_effects=rbind(Direct_effects,
                       tibble(q1=apply(effect_pred$Direct,2,quantile,.025)[-1],
                              q3=apply(effect_pred$Direct,2,quantile,.975)[-1],
                              q2=apply(effect_pred$Direct,2,median)[-1],
                              pval=apply(effect_pred$Direct,2,get_bootstrapped_pval)[-1],
                              Term=names(apply(effect_pred$Direct,2,quantile,.025)[-1]),
                              Response=predictors))
  
}


effect_on_q=SEM_resilience$`Bootstrapped Effects`$q
effect_on_cov=SEM_resilience$`Bootstrapped Effects`$Cover
effect_on_dist=SEM_resilience$`Bootstrapped Effects`$Resid.mod



Total_effects=rbind(tibble(q1=apply(effect_on_q$Total[,-1],2,quantile,.025),
                           q3=apply(effect_on_q$Total[,-1],2,quantile,.975),
                           q2=apply(effect_on_q$Total[,-1],2,median),
                           pval=apply(effect_on_q$Total[,-1],2,get_bootstrapped_pval),
                           Term=names(apply(effect_on_q$Total[,-1],2,quantile,.025)),
                           Response="q"),
                    tibble(q1=apply(effect_on_cov$Total[,-1],2,quantile,.025),
                           q3=apply(effect_on_cov$Total[,-1],2,quantile,.975),
                           q2=apply(effect_on_cov$Total[,-1],2,median),
                           pval=apply(effect_on_cov$Total[,-1],2,get_bootstrapped_pval),
                           Term=names(apply(effect_on_cov$Total[,-1],2,quantile,.025)),
                           Response="Cover"),
                    tibble(q1=apply(effect_on_dist$Total[,-1],2,quantile,.025),
                           q3=apply(effect_on_dist$Total[,-1],2,quantile,.975),
                           q2=apply(effect_on_dist$Total[,-1],2,median),
                           pval=apply(effect_on_dist$Total[,-1],2,get_bootstrapped_pval),
                           Term=names(apply(effect_on_dist$Total[,-1],2,quantile,.025)),
                           Response="Distance"))


write.table(Total_effects,"../Data/SEM/Total_effects_SEM_resilience.csv",sep=";")
write.table(Direct_effects,"../Data/SEM/Direct_effects_SEM_resilience.csv",sep=";")

# ---------------------- Step 7: Changing grazing reference level ----
## >> 1) Changing to grazing intensity 1 ----

dir.create("./Linear_models_factor_ref1",showWarnings = F)
dir.create("./Linear_models_factor_ref1/Importance",showWarnings = F)
dir.create("./Linear_models_factor_ref1/Estimator",showWarnings = F)
dir.create("./Linear_models_factor_ref1/Keep_data",showWarnings = F)
dir.create("./Linear_models_factor_ref1/Keep_models",showWarnings = F)

Run_model_importance_aridity_no_inter=function(id){
  
  list_mod=rbind(expand.grid(with_cover=c(T,F),
                             Stats=c("perim_area_scaling","fmax_psd","PL_expo",
                                     "core_area_land","core_area","flow_length","PLR","KS_dist","Shape_metric")))
  
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
      + rho_p + Type_veg 
      + Sand + Org_C_v 
      + (1|Site_ID)")))
    
  } else{
    
    formula_mod=formula(formula_(paste("value ~ Long_cos + Long_sin + Lattitude + Slope + Elevation 
      + Aridity + Grazing 
      + Type_veg
      + Sand + Org_C_v 
      + (1|Site_ID)")))
  }
  
  save=d_data_mod
  grazing_intensity='all'
  d_data_mod=save
  
  d_data_mod$Grazing=as.factor(d_data_mod$Grazing)
  d_data_mod$Grazing = relevel(d_data_mod$Grazing,ref = "1")
  
  model_spa_stat  = lmer(formula_mod, d_data_mod,
                         na.action = na.fail ,REML ="FALSE")
  
  #we remove potential outliers
  rm.outliers = romr.fnc(model_spa_stat, d_data_mod, trim=2.5)
  d_data_out = rm.outliers$data
  
  #saving data
  write.table(d_data_out,paste0("../Data/Linear_models_factor_ref1/Keep_data/Data_",stat,"_",with_cover,"_aridity_no_inter_",grazing_intensity,".csv"))
  
  model_spa_stat = lmer(formula_mod, data = d_data_out, 
                        na.action = na.fail,REML ="FALSE")
  
  # #do some model selection
  select_model=dredge(model_spa_stat, subset = ~ Slope & Elevation &
                        Long_cos & Long_sin & Lattitude &
                        dc(Woody & Grazing, Woody : Grazing) &
                        dc(Aridity & Grazing, Aridity : Grazing) &
                        dc(rho_p & Grazing, rho_p : Grazing) &
                        dc(Org_C_v & Grazing, Org_C_v : Grazing) &
                        dc(Sand  & Grazing, Sand  : Grazing) &
                        dc(Type_veg & Grazing, Type_veg : Grazing),
                      extra=c(R2m=function(x) r.squaredGLMM(x)[1],
                              R2c=function(x) r.squaredGLMM(x)[2]),
                      options(na.action = "na.fail") )
  
  if (dim(select_model%>%filter(., delta<2))[1]==1){ #one model
    
    result_select=model_spa_stat
    R2=select_model%>%filter(., AICc<min(AICc)+2)
    importance_mod=tibble(Interactions=1,
                          Sand =1,
                          Aridity=1,
                          Org_C=1,
                          Type_veg=1,
                          Grazing=1,
                          Cover=1,
                          Woody=1)
    
    d_all=rbind(d_all,cbind(tibble(Sp_stat=stat,
                                   N_outliers=rm.outliers$n.removed,
                                   R2m=mean(R2$R2m),R2C=mean(R2$R2c),
                                   With_cover=with_cover,
                                   Grazing_intensity=grazing_intensity),
                            importance_mod))
    
  }else{
    
    result_select=model.avg(select_model, subset = delta < 2)
    
    #Get the importance of each metric
    importance_mod=sw(result_select)
    R2=select_model%>%filter(., AICc<min(AICc)+2)
    
    d_all=rbind(d_all,cbind(tibble(Sp_stat=stat,
                                   N_outliers=rm.outliers$n.removed,
                                   R2m=mean(R2$R2m),R2C=mean(R2$R2c)),
                            Aggregate_importance(importance_mod,T))%>%
                  add_column(., With_cover=with_cover,
                             Grazing_intensity=grazing_intensity))
    
  }
  
  saveRDS(model_spa_stat,paste0("../Data/Linear_models_factor_ref1/Keep_models/Mod_",stat,"_",with_cover,"_aridity_no_inter_",grazing_intensity,".rds"))
  
  summary_coef=confint(result_select)
  
  
  #Merge in a df
  d_all2=rbind(d_all2,tibble(Median=confint(result_select,level = 1e-9)[-1,1],
                             q1=summary_coef[-1,1], 
                             q3=summary_coef[-1,2],
                             # pvalue=summary(result_select)$coefmat.full[-1,5],
                             term=rownames(summary_coef)[-1],
                             Stat=stat,
                             R2m=mean(R2$R2m),
                             R2c=mean(R2$R2c),
                             Grazing_intensity=grazing_intensity)%>%
                 add_column(., With_cover=with_cover))
  
  
  
  
  
  #Saving each DF
  write.table(d_all,paste0("../Data/Linear_models_factor_ref1/Importance/Importance_",stat,"_",with_cover,"_aridity_no_inter.csv"),sep=";")
  write.table(d_all2,paste0("../Data/Linear_models_factor_ref1/Estimator/Estimators_model_",stat,"_",with_cover,"_aridity_no_inter.csv"),sep=";")
  
}

library(parallel)

mclapply(1:18,Run_model_importance_aridity_no_inter,mc.cores = 18)

d_all=d_all2=tibble()
for (k in list.files("../Data/Linear_models_factor_ref1/Importance",paste0("aridity_no_inter"))){
  d=read.table(paste0("../Data/Linear_models_factor_ref1/Importance/",k),sep=";")
  if (nrow(d)!=0){
    d_all=rbind(d_all,d)
  }
}
for (k in list.files("../Data/Linear_models_factor_ref1/Estimator",paste0("aridity_no_inter"))){
  d2=read.table(paste0("../Data/Linear_models_factor_ref1/Estimator/",k),sep=";")
  if (nrow(d)!=0){
    d_all2=rbind(d_all2,d2)
  }
}
write.table(d_all,paste0("../Data/Linear_models_factor_ref1/Importance_aridity_no_inter_factor_ref1.csv"),sep=";")
write.table(d_all2,paste0("../Data/Linear_models_factor_ref1/Estimators_model_aridity_no_inter_factor_ref1.csv"),sep=";")


d_data=read.table("../Data/Spatial_structure_grazing.csv",sep=";")%>%
  Closer_to_normality(.)
d_data[,c(1,9:29,36,38:42,45:ncol(d_data))] = 
  apply(d_data[,c(1,9:29,36,38:42,45:ncol(d_data))],2,z_tranform)
d_data=Perform_PCA_spatial_struc(d_data)
d_data$Grazing=as.factor(d_data$Grazing)
save=d_data



boot_function_lm_graz = function(formula, data, indices) {
  d = data[indices,] 
  fit = lm(formula, data=d) 
  return(summary(fit)$coefficient[2:4,1])
}
boot_function_lm_arid = function(formula, data, indices) {
  d = data[indices,] 
  fit = lm(formula, data=d) 
  return(summary(fit)$coefficient[2,1])
}

grazing_intensity="all"
for (with_interaction in c(T,F)){
  
  d_slope=tibble()
  
  for (k in c("perim_area_scaling","fmax_psd","PL_expo",
              "core_area_land","core_area",
              "flow_length","PLR","KS_dist","Shape_metric")){
    
    #for each we plot the slope against the partial residuals with and without cover
    
    if (k !="rho_p"){
      
      d_data_out=read.table(paste0("../Data/Linear_models_factor_ref1/Keep_data/Data_",k,"_TRUE_aridity_",
                                   ifelse(with_interaction,"","no_inter_"),grazing_intensity,".csv"),sep=" ")
      model_spa_stat=readRDS(paste0("../Data/Linear_models_factor_ref1/Keep_models/Mod_",k,"_TRUE_aridity_",
                                    ifelse(with_interaction,"","no_inter_"),grazing_intensity,".rds"))
      
      d_data_out$Grazing = as.factor(d_data_out$Grazing)
      d_data_out$Grazing = relevel(d_data_out$Grazing,ref = "1")
      
      resid_mod=visreg::visreg(fit = model_spa_stat,xvar="Grazing",plot=F) 
      
      mod_cov = boot(data=resid_mod$res, statistic=boot_function_lm_graz, R=1000, formula=visregRes~Grazing)
      
      d_slope=rbind(d_slope,
                    tibble(pval=apply(mod_cov$t,2,twoside_pvalue),
                           q2=apply((mod_cov$t),2,median),
                           q1=apply(mod_cov$t,2,quantile,.025),
                           q3=apply(mod_cov$t,2,quantile,.975),
                           With_cover=T,Grazing_intensity=grazing_intensity,
                           ID_grazing=c("0","2","3"),
                           Stat=k,Driver="Grazing"
                    ))
      
      resid_mod=visreg::visreg(fit = model_spa_stat,xvar="Aridity",plot=F) 
      
      mod_cov = boot(data=resid_mod$res, statistic=boot_function_lm_arid, R=1000, formula=visregRes~Aridity)
      
      d_slope=rbind(d_slope,
                    tibble(pval=apply(mod_cov$t,2,twoside_pvalue),
                           q2=apply((mod_cov$t),2,median),
                           q1=apply(mod_cov$t,2,quantile,.025),
                           q3=apply(mod_cov$t,2,quantile,.975),
                           With_cover=T,Grazing_intensity=grazing_intensity,
                           ID_grazing=c("0","2","3"),
                           Stat=k,Driver="Aridity"
                    ))
    }
    
    d_data_out=read.table(paste0("../Data/Linear_models_factor_ref1/Keep_data/Data_",k,"_FALSE_aridity_",
                                 ifelse(with_interaction,"","no_inter_"),grazing_intensity,".csv"),sep=" ")
    model_spa_stat=readRDS(paste0("../Data/Linear_models_factor_ref1/Keep_models/Mod_",k,"_FALSE_aridity_",
                                  ifelse(with_interaction,"","no_inter_"),grazing_intensity,".rds"))
    
    d_data_out$Grazing = as.factor(d_data_out$Grazing)
    d_data_out$Grazing = relevel(d_data_out$Grazing,ref = "1")
    
    resid_mod=visreg::visreg(fit = model_spa_stat,xvar="Grazing",plot=F) 
    
    mod_cov = boot(data=resid_mod$res, statistic=boot_function_lm_graz, R=1000, formula=visregRes~Grazing)
    
    d_slope=rbind(d_slope,
                  tibble(pval=apply(mod_cov$t,2,twoside_pvalue),
                         q2=apply((mod_cov$t),2,median),
                         q1=apply(mod_cov$t,2,quantile,.025),
                         q3=apply(mod_cov$t,2,quantile,.975),
                         With_cover=F,Grazing_intensity=grazing_intensity,
                         ID_grazing=c("0","2","3"),
                         Stat=k,Driver="Grazing"
                  ))
    
    
    resid_mod=visreg::visreg(fit = model_spa_stat,xvar="Aridity",plot=F) 
    
    mod_cov = boot(data=resid_mod$res, statistic=boot_function_lm_arid, R=1000, formula=visregRes~Aridity)
    
    d_slope=rbind(d_slope,
                  tibble(pval=apply(mod_cov$t,2,twoside_pvalue),
                         q2=apply((mod_cov$t),2,median),
                         q1=apply(mod_cov$t,2,quantile,.025),
                         q3=apply(mod_cov$t,2,quantile,.975),
                         q1=apply(mod_cov$t,2,quantile,.025),
                         q3=apply(mod_cov$t,2,quantile,.975),
                         With_cover=F,Grazing_intensity=grazing_intensity,
                         ID_grazing=c("0","2","3"),
                         Stat=k,Driver="Aridity"
                  ))
  }
  write.table(d_slope,paste0("../Data/Linear_models_factor_ref1/Slope_partial_residuals_aridity",
                             ifelse(with_interaction,"","_no_inter"),
                             "_ref1.csv"),sep=";")
  
}

## >> 2) Changing to grazing intensity 2 ----
dir.create("./Linear_models_factor_ref2",showWarnings = F)
dir.create("./Linear_models_factor_ref2/Importance",showWarnings = F)
dir.create("./Linear_models_factor_ref2/Estimator",showWarnings = F)
dir.create("./Linear_models_factor_ref2/Keep_data",showWarnings = F)
dir.create("./Linear_models_factor_ref2/Keep_models",showWarnings = F)

Run_model_importance_aridity_no_inter=function(id){
  
  list_mod=rbind(expand.grid(with_cover=c(T,F),
                             Stats=c("perim_area_scaling","fmax_psd","PL_expo",
                                     "core_area_land","core_area","flow_length","PLR","KS_dist","Shape_metric")))
  
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
      + rho_p + Type_veg 
      + Sand + Org_C_v 
      + (1|Site_ID)")))
    
  } else{
    
    formula_mod=formula(formula_(paste("value ~ Long_cos + Long_sin + Lattitude + Slope + Elevation 
      + Aridity + Grazing 
      + Type_veg
      + Sand + Org_C_v 
      + (1|Site_ID)")))
  }
  
  save=d_data_mod
  grazing_intensity='all'
  d_data_mod=save
  
  d_data_mod$Grazing=as.factor(d_data_mod$Grazing)
  d_data_mod$Grazing = relevel(d_data_mod$Grazing,ref = "2")
  
  model_spa_stat  = lmer(formula_mod, d_data_mod,
                         na.action = na.fail ,REML ="FALSE")
  
  #we remove potential outliers
  rm.outliers = romr.fnc(model_spa_stat, d_data_mod, trim=2.5)
  d_data_out = rm.outliers$data
  
  #saving data
  write.table(d_data_out,paste0("../Data/Linear_models_factor_ref2/Keep_data/Data_",stat,"_",with_cover,"_aridity_no_inter_",grazing_intensity,".csv"))
  
  model_spa_stat = lmer(formula_mod, data = d_data_out, 
                        na.action = na.fail,REML ="FALSE")
  
  # #do some model selection
  select_model=dredge(model_spa_stat, subset = ~ Slope & Elevation &
                        Long_cos & Long_sin & Lattitude &
                        dc(Woody & Grazing, Woody : Grazing) &
                        dc(Aridity & Grazing, Aridity : Grazing) &
                        dc(rho_p & Grazing, rho_p : Grazing) &
                        dc(Org_C_v & Grazing, Org_C_v : Grazing) &
                        dc(Sand  & Grazing, Sand  : Grazing) &
                        dc(Type_veg & Grazing, Type_veg : Grazing),
                      extra=c(R2m=function(x) r.squaredGLMM(x)[1],
                              R2c=function(x) r.squaredGLMM(x)[2]),
                      options(na.action = "na.fail") )
  
  if (dim(select_model%>%filter(., delta<2))[1]==1){ #one model
    
    result_select=model_spa_stat
    R2=select_model%>%filter(., AICc<min(AICc)+2)
    importance_mod=tibble(Interactions=1,
                          Sand =1,
                          Aridity=1,
                          Org_C=1,
                          Type_veg=1,
                          Grazing=1,
                          Cover=1,
                          Woody=1)
    
    d_all=rbind(d_all,cbind(tibble(Sp_stat=stat,
                                   N_outliers=rm.outliers$n.removed,
                                   R2m=mean(R2$R2m),R2C=mean(R2$R2c),
                                   With_cover=with_cover,
                                   Grazing_intensity=grazing_intensity),
                            importance_mod))
    
  }else{
    
    result_select=model.avg(select_model, subset = delta < 2)
    
    #Get the importance of each metric
    importance_mod=sw(result_select)
    R2=select_model%>%filter(., AICc<min(AICc)+2)
    
    d_all=rbind(d_all,cbind(tibble(Sp_stat=stat,
                                   N_outliers=rm.outliers$n.removed,
                                   R2m=mean(R2$R2m),R2C=mean(R2$R2c)),
                            Aggregate_importance(importance_mod,T))%>%
                  add_column(., With_cover=with_cover,
                             Grazing_intensity=grazing_intensity))
    
  }
  
  saveRDS(model_spa_stat,paste0("../Data/Linear_models_factor_ref2/Keep_models/Mod_",stat,"_",with_cover,"_aridity_no_inter_",grazing_intensity,".rds"))
  
  summary_coef=confint(result_select)
  
  
  #Merge in a df
  d_all2=rbind(d_all2,tibble(Median=confint(result_select,level = 1e-9)[-1,1],
                             q1=summary_coef[-1,1], 
                             q3=summary_coef[-1,2],
                             # pvalue=summary(result_select)$coefmat.full[-1,5],
                             term=rownames(summary_coef)[-1],
                             Stat=stat,
                             R2m=mean(R2$R2m),
                             R2c=mean(R2$R2c),
                             Grazing_intensity=grazing_intensity)%>%
                 add_column(., With_cover=with_cover))
  
  
  
  
  
  #Saving each DF
  write.table(d_all,paste0("../Data/Linear_models_factor_ref2/Importance/Importance_",stat,"_",with_cover,"_aridity_no_inter.csv"),sep=";")
  write.table(d_all2,paste0("../Data/Linear_models_factor_ref2/Estimator/Estimators_model_",stat,"_",with_cover,"_aridity_no_inter.csv"),sep=";")
  
}

library(parallel)

mclapply(1:18,Run_model_importance_aridity_no_inter,mc.cores = 18)

d_all=d_all2=tibble()
for (k in list.files("../Data/Linear_models_factor_ref2/Importance",paste0("aridity_no_inter"))){
  d=read.table(paste0("../Data/Linear_models_factor_ref2/Importance/",k),sep=";")
  if (nrow(d)!=0){
    d_all=rbind(d_all,d)
  }
}
for (k in list.files("../Data/Linear_models_factor_ref2/Estimator",paste0("aridity_no_inter"))){
  d2=read.table(paste0("../Data/Linear_models_factor_ref2/Estimator/",k),sep=";")
  if (nrow(d)!=0){
    d_all2=rbind(d_all2,d2)
  }
}
write.table(d_all,paste0("../Data/Linear_models_factor_ref2/Importance_aridity_no_inter_factor_ref2.csv"),sep=";")
write.table(d_all2,paste0("../Data/Linear_models_factor_ref2/Estimators_model_aridity_no_inter_factor_ref2.csv"),sep=";")


d_data=read.table("../Data/Spatial_structure_grazing.csv",sep=";")%>%
  Closer_to_normality(.)
d_data[,c(1,9:29,36,38:42,45:ncol(d_data))] = 
  apply(d_data[,c(1,9:29,36,38:42,45:ncol(d_data))],2,z_tranform)
d_data=Perform_PCA_spatial_struc(d_data)
d_data$Grazing=as.factor(d_data$Grazing)
save=d_data



boot_function_lm_graz = function(formula, data, indices) {
  d = data[indices,] 
  fit = lm(formula, data=d) 
  return(summary(fit)$coefficient[2:4,1])
}
boot_function_lm_arid = function(formula, data, indices) {
  d = data[indices,] 
  fit = lm(formula, data=d) 
  return(summary(fit)$coefficient[2,1])
}

grazing_intensity="all"
for (with_interaction in c(T,F)){
  
  d_slope=tibble()
  
  for (k in c("perim_area_scaling","fmax_psd","PL_expo",
              "core_area_land","core_area",
              "flow_length","PLR","KS_dist","Shape_metric")){
    
    #for each we plot the slope against the partial residuals with and without cover
    
    if (k !="rho_p"){
      
      d_data_out=read.table(paste0("../Data/Linear_models_factor_ref2/Keep_data/Data_",k,"_TRUE_aridity_",
                                   ifelse(with_interaction,"","no_inter_"),grazing_intensity,".csv"),sep=" ")
      model_spa_stat=readRDS(paste0("../Data/Linear_models_factor_ref2/Keep_models/Mod_",k,"_TRUE_aridity_",
                                    ifelse(with_interaction,"","no_inter_"),grazing_intensity,".rds"))
      
      d_data_out$Grazing = as.factor(d_data_out$Grazing)
      d_data_out$Grazing = relevel(d_data_out$Grazing,ref = "2")
      
      resid_mod=visreg::visreg(fit = model_spa_stat,xvar="Grazing",plot=F) 
      
      mod_cov = boot(data=resid_mod$res, statistic=boot_function_lm_graz, R=1000, formula=visregRes~Grazing)
      
      d_slope=rbind(d_slope,
                    tibble(pval=apply(mod_cov$t,2,twoside_pvalue),
                           q2=apply((mod_cov$t),2,median),
                           q1=apply(mod_cov$t,2,quantile,.025),
                           q3=apply(mod_cov$t,2,quantile,.975),
                           With_cover=T,Grazing_intensity=grazing_intensity,
                           ID_grazing=c("0","1","3"),
                           Stat=k,Driver="Grazing"
                    ))
      
      resid_mod=visreg::visreg(fit = model_spa_stat,xvar="Aridity",plot=F) 
      
      mod_cov = boot(data=resid_mod$res, statistic=boot_function_lm_arid, R=1000, formula=visregRes~Aridity)
      
      d_slope=rbind(d_slope,
                    tibble(pval=apply(mod_cov$t,2,twoside_pvalue),
                           q2=apply((mod_cov$t),2,median),
                           q1=apply(mod_cov$t,2,quantile,.025),
                           q3=apply(mod_cov$t,2,quantile,.975),
                           With_cover=T,Grazing_intensity=grazing_intensity,
                           ID_grazing=c("0","1","3"),
                           Stat=k,Driver="Aridity"
                    ))
    }
    
    d_data_out=read.table(paste0("../Data/Linear_models_factor_ref2/Keep_data/Data_",k,"_FALSE_aridity_",
                                 ifelse(with_interaction,"","no_inter_"),grazing_intensity,".csv"),sep=" ")
    model_spa_stat=readRDS(paste0("../Data/Linear_models_factor_ref2/Keep_models/Mod_",k,"_FALSE_aridity_",
                                  ifelse(with_interaction,"","no_inter_"),grazing_intensity,".rds"))
    
    d_data_out$Grazing = as.factor(d_data_out$Grazing)
    d_data_out$Grazing = relevel(d_data_out$Grazing,ref = "2")
    
    resid_mod=visreg::visreg(fit = model_spa_stat,xvar="Grazing",plot=F) 
    
    mod_cov = boot(data=resid_mod$res, statistic=boot_function_lm_graz, R=1000, formula=visregRes~Grazing)
    
    d_slope=rbind(d_slope,
                  tibble(pval=apply(mod_cov$t,2,twoside_pvalue),
                         q2=apply((mod_cov$t),2,median),
                         q1=apply(mod_cov$t,2,quantile,.025),
                         q3=apply(mod_cov$t,2,quantile,.975),
                         With_cover=F,Grazing_intensity=grazing_intensity,
                         ID_grazing=c("0","1","3"),
                         Stat=k,Driver="Grazing"
                  ))
    
    
    resid_mod=visreg::visreg(fit = model_spa_stat,xvar="Aridity",plot=F) 
    
    mod_cov = boot(data=resid_mod$res, statistic=boot_function_lm_arid, R=1000, formula=visregRes~Aridity)
    
    d_slope=rbind(d_slope,
                  tibble(pval=apply(mod_cov$t,2,twoside_pvalue),
                         q2=apply((mod_cov$t),2,median),
                         q1=apply(mod_cov$t,2,quantile,.025),
                         q3=apply(mod_cov$t,2,quantile,.975),
                         With_cover=F,Grazing_intensity=grazing_intensity,
                         ID_grazing=c("0","1","3"),
                         Stat=k,Driver="Aridity"
                  ))
  }
  write.table(d_slope,paste0("../Data/Linear_models_factor_ref2/Slope_partial_residuals_aridity",
                             ifelse(with_interaction,"","_no_inter"),
                             "_ref2.csv"),sep=";")
  
}


## >> 3) Changing to grazing intensity 3 ----
dir.create("./Linear_models_factor_ref3",showWarnings = F)
dir.create("./Linear_models_factor_ref3/Importance",showWarnings = F)
dir.create("./Linear_models_factor_ref3/Estimator",showWarnings = F)
dir.create("./Linear_models_factor_ref3/Keep_data",showWarnings = F)
dir.create("./Linear_models_factor_ref3/Keep_models",showWarnings = F)

Run_model_importance_aridity_no_inter=function(id){
  
  list_mod=rbind(expand.grid(with_cover=c(T,F),
                             Stats=c("perim_area_scaling","fmax_psd","PL_expo",
                                     "core_area_land","core_area","flow_length","PLR","KS_dist","Shape_metric")))
  
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
      + rho_p + Type_veg 
      + Sand + Org_C_v 
      + (1|Site_ID)")))
    
  } else{
    
    formula_mod=formula(formula_(paste("value ~ Long_cos + Long_sin + Lattitude + Slope + Elevation 
      + Aridity + Grazing 
      + Type_veg
      + Sand + Org_C_v 
      + (1|Site_ID)")))
  }
  
  save=d_data_mod
  grazing_intensity='all'
  d_data_mod=save
  
  d_data_mod$Grazing=as.factor(d_data_mod$Grazing)
  d_data_mod$Grazing = relevel(d_data_mod$Grazing,ref = "3")
  
  model_spa_stat  = lmer(formula_mod, d_data_mod,
                         na.action = na.fail ,REML ="FALSE")
  
  #we remove potential outliers
  rm.outliers = romr.fnc(model_spa_stat, d_data_mod, trim=2.5)
  d_data_out = rm.outliers$data
  
  #saving data
  write.table(d_data_out,paste0("../Data/Linear_models_factor_ref2/Keep_data/Data_",stat,"_",with_cover,"_aridity_no_inter_",grazing_intensity,".csv"))
  
  model_spa_stat = lmer(formula_mod, data = d_data_out, 
                        na.action = na.fail,REML ="FALSE")
  
  # #do some model selection
  select_model=dredge(model_spa_stat, subset = ~ Slope & Elevation &
                        Long_cos & Long_sin & Lattitude &
                        dc(Woody & Grazing, Woody : Grazing) &
                        dc(Aridity & Grazing, Aridity : Grazing) &
                        dc(rho_p & Grazing, rho_p : Grazing) &
                        dc(Org_C_v & Grazing, Org_C_v : Grazing) &
                        dc(Sand  & Grazing, Sand  : Grazing) &
                        dc(Type_veg & Grazing, Type_veg : Grazing),
                      extra=c(R2m=function(x) r.squaredGLMM(x)[1],
                              R2c=function(x) r.squaredGLMM(x)[2]),
                      options(na.action = "na.fail") )
  
  if (dim(select_model%>%filter(., delta<2))[1]==1){ #one model
    
    result_select=model_spa_stat
    R2=select_model%>%filter(., AICc<min(AICc)+2)
    importance_mod=tibble(Interactions=1,
                          Sand =1,
                          Aridity=1,
                          Org_C=1,
                          Type_veg=1,
                          Grazing=1,
                          Cover=1,
                          Woody=1)
    
    d_all=rbind(d_all,cbind(tibble(Sp_stat=stat,
                                   N_outliers=rm.outliers$n.removed,
                                   R2m=mean(R2$R2m),R2C=mean(R2$R2c),
                                   With_cover=with_cover,
                                   Grazing_intensity=grazing_intensity),
                            importance_mod))
    
  }else{
    
    result_select=model.avg(select_model, subset = delta < 2)
    
    #Get the importance of each metric
    importance_mod=sw(result_select)
    R2=select_model%>%filter(., AICc<min(AICc)+2)
    
    d_all=rbind(d_all,cbind(tibble(Sp_stat=stat,
                                   N_outliers=rm.outliers$n.removed,
                                   R2m=mean(R2$R2m),R2C=mean(R2$R2c)),
                            Aggregate_importance(importance_mod,T))%>%
                  add_column(., With_cover=with_cover,
                             Grazing_intensity=grazing_intensity))
    
  }
  
  saveRDS(model_spa_stat,paste0("../Data/Linear_models_factor_ref3/Keep_models/Mod_",stat,"_",with_cover,"_aridity_no_inter_",grazing_intensity,".rds"))
  
  summary_coef=confint(result_select)
  
  
  #Merge in a df
  d_all2=rbind(d_all2,tibble(Median=confint(result_select,level = 1e-9)[-1,1],
                             q1=summary_coef[-1,1], 
                             q3=summary_coef[-1,2],
                             # pvalue=summary(result_select)$coefmat.full[-1,5],
                             term=rownames(summary_coef)[-1],
                             Stat=stat,
                             R2m=mean(R2$R2m),
                             R2c=mean(R2$R2c),
                             Grazing_intensity=grazing_intensity)%>%
                 add_column(., With_cover=with_cover))
  
  
  
  
  
  #Saving each DF
  write.table(d_all,paste0("../Data/Linear_models_factor_ref3/Importance/Importance_",stat,"_",with_cover,"_aridity_no_inter.csv"),sep=";")
  write.table(d_all2,paste0("../Data/Linear_models_factor_ref3/Estimator/Estimators_model_",stat,"_",with_cover,"_aridity_no_inter.csv"),sep=";")
  
}

library(parallel)

mclapply(1:18,Run_model_importance_aridity_no_inter,mc.cores = 18)

d_all=d_all2=tibble()
for (k in list.files("../Data/Linear_models_factor_ref3/Importance",paste0("aridity_no_inter"))){
  d=read.table(paste0("../Data/Linear_models_factor_ref3/Importance/",k),sep=";")
  if (nrow(d)!=0){
    d_all=rbind(d_all,d)
  }
}
for (k in list.files("../Data/Linear_models_factor_ref3/Estimator",paste0("aridity_no_inter"))){
  d2=read.table(paste0("../Data/Linear_models_factor_ref3/Estimator/",k),sep=";")
  if (nrow(d)!=0){
    d_all2=rbind(d_all2,d2)
  }
}
write.table(d_all,paste0("../Data/Linear_models_factor_ref3/Importance_aridity_no_inter_factor_ref2.csv"),sep=";")
write.table(d_all2,paste0("../Data/Linear_models_factor_ref3/Estimators_model_aridity_no_inter_factor_ref2.csv"),sep=";")


d_data=read.table("../Data/Spatial_structure_grazing.csv",sep=";")%>%
  Closer_to_normality(.)
d_data[,c(1,9:29,36,38:42,45:ncol(d_data))] = 
  apply(d_data[,c(1,9:29,36,38:42,45:ncol(d_data))],2,z_tranform)
d_data=Perform_PCA_spatial_struc(d_data)
d_data$Grazing=as.factor(d_data$Grazing)
save=d_data



boot_function_lm_graz = function(formula, data, indices) {
  d = data[indices,] 
  fit = lm(formula, data=d) 
  return(summary(fit)$coefficient[2:4,1])
}
boot_function_lm_arid = function(formula, data, indices) {
  d = data[indices,] 
  fit = lm(formula, data=d) 
  return(summary(fit)$coefficient[2,1])
}

grazing_intensity="all"
for (with_interaction in c(T,F)){
  
  d_slope=tibble()
  
  for (k in c("perim_area_scaling","fmax_psd","PL_expo",
              "core_area_land","core_area",
              "flow_length","PLR","KS_dist","Shape_metric")){
    
    #for each we plot the slope against the partial residuals with and without cover
    
    if (k !="rho_p"){
      
      d_data_out=read.table(paste0("../Data/Linear_models_factor_ref3/Keep_data/Data_",k,"_TRUE_aridity_",
                                   ifelse(with_interaction,"","no_inter_"),grazing_intensity,".csv"),sep=" ")
      model_spa_stat=readRDS(paste0("../Data/Linear_models_factor_ref3/Keep_models/Mod_",k,"_TRUE_aridity_",
                                    ifelse(with_interaction,"","no_inter_"),grazing_intensity,".rds"))
      
      d_data_out$Grazing = as.factor(d_data_out$Grazing)
      d_data_out$Grazing = relevel(d_data_out$Grazing,ref = "3")
      
      resid_mod=visreg::visreg(fit = model_spa_stat,xvar="Grazing",plot=F) 
      
      mod_cov = boot(data=resid_mod$res, statistic=boot_function_lm_graz, R=1000, formula=visregRes~Grazing)
      
      d_slope=rbind(d_slope,
                    tibble(pval=apply(mod_cov$t,2,twoside_pvalue),
                           q2=apply((mod_cov$t),2,median),
                           q1=apply(mod_cov$t,2,quantile,.025),
                           q3=apply(mod_cov$t,2,quantile,.975),
                           With_cover=T,Grazing_intensity=grazing_intensity,
                           ID_grazing=c("0","1","2"),
                           Stat=k,Driver="Grazing"
                    ))
      
      resid_mod=visreg::visreg(fit = model_spa_stat,xvar="Aridity",plot=F) 
      
      mod_cov = boot(data=resid_mod$res, statistic=boot_function_lm_arid, R=1000, formula=visregRes~Aridity)
      
      d_slope=rbind(d_slope,
                    tibble(pval=apply(mod_cov$t,2,twoside_pvalue),
                           q2=apply((mod_cov$t),2,median),
                           q1=apply(mod_cov$t,2,quantile,.025),
                           q3=apply(mod_cov$t,2,quantile,.975),
                           With_cover=T,Grazing_intensity=grazing_intensity,
                           ID_grazing=c("0","1","2"),
                           Stat=k,Driver="Aridity"
                    ))
    }
    
    d_data_out=read.table(paste0("../Data/Linear_models_factor_ref3/Keep_data/Data_",k,"_FALSE_aridity_",
                                 ifelse(with_interaction,"","no_inter_"),grazing_intensity,".csv"),sep=" ")
    model_spa_stat=readRDS(paste0("../Data/Linear_models_factor_ref3/Keep_models/Mod_",k,"_FALSE_aridity_",
                                  ifelse(with_interaction,"","no_inter_"),grazing_intensity,".rds"))
    
    d_data_out$Grazing = as.factor(d_data_out$Grazing)
    d_data_out$Grazing = relevel(d_data_out$Grazing,ref = "3")
    
    resid_mod=visreg::visreg(fit = model_spa_stat,xvar="Grazing",plot=F) 
    
    mod_cov = boot(data=resid_mod$res, statistic=boot_function_lm_graz, R=1000, formula=visregRes~Grazing)
    
    d_slope=rbind(d_slope,
                  tibble(pval=apply(mod_cov$t,2,twoside_pvalue),
                         q2=apply((mod_cov$t),2,median),
                         q1=apply(mod_cov$t,2,quantile,.025),
                         q3=apply(mod_cov$t,2,quantile,.975),
                         With_cover=F,Grazing_intensity=grazing_intensity,
                         ID_grazing=c("0","1","2"),
                         Stat=k,Driver="Grazing"
                  ))
    
    
    resid_mod=visreg::visreg(fit = model_spa_stat,xvar="Aridity",plot=F) 
    
    mod_cov = boot(data=resid_mod$res, statistic=boot_function_lm_arid, R=1000, formula=visregRes~Aridity)
    
    d_slope=rbind(d_slope,
                  tibble(pval=apply(mod_cov$t,2,twoside_pvalue),
                         q2=apply((mod_cov$t),2,median),
                         q1=apply(mod_cov$t,2,quantile,.025),
                         q3=apply(mod_cov$t,2,quantile,.975),
                         With_cover=F,Grazing_intensity=grazing_intensity,
                         ID_grazing=c("0","1","2"),
                         Stat=k,Driver="Aridity"
                  ))
  }
  write.table(d_slope,paste0("../Data/Linear_models_factor_ref3/Slope_partial_residuals_aridity",
                             ifelse(with_interaction,"","_no_inter"),
                             "_ref2.csv"),sep=";")
  
}


