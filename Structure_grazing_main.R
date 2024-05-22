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
    add_column(.,
               Full_name=name,Resolution=Get_spatial_resolution(landscape), #full name & spatial resolution
               Site_ID=id_site, # ID of the site
               Sub_id=sub_id) # id of the sublandscape taken at a given site
  
  write.table(d,paste0("../Data/Metrics/Metric_",name,".csv"),sep=";")  
}

list_f=list.files("../Data/Landscapes/Binary_landscapes/",".csv")


library(parallel)
mclapply(list_f[grep(list_f,pattern = "biodesert")],Compute_all_metrics_biodesert,mc.cores = 10)


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
             Herbivores=sapply(1:nrow(.),function(x){return(d_biodesert$`Herbivore 1`[which(d_biodesert$ID==.$Site_ID[x])])}),
             Org_C=sapply(1:nrow(.),function(x){return(d_biodesert$`ORC veg`[which(d_biodesert$ID==.$Site_ID[x])]- #difference between plant and bare soil for organic carbon
                                                         d_biodesert$`ORC b`[which(d_biodesert$ID==.$Site_ID[x])])}),
             Sand=sapply(1:nrow(.),function(x){return(d_biodesert$`SAC veg`[which(d_biodesert$ID==.$Site_ID[x])])}),
             Shannon_div=sapply(1:nrow(.),function(x){return(d_biodesert$DIV_SHQ[which(d_biodesert$ID==.$Site_ID[x])])}),
             Simpson_div=sapply(1:nrow(.),function(x){return(d_biodesert$DIV_SIQ[which(d_biodesert$ID==.$Site_ID[x])])}),
             Fertility=(sapply(1:nrow(.),function(x){return(d_biodesert$Fertility[which(d_biodesert$ID==.$Site_ID[x])])})),
             Forage_Quality=(sapply(1:nrow(.),function(x){return(d_biodesert$Forage_Quality[which(d_biodesert$ID==.$Site_ID[x])])})),
             C_stock=(sapply(1:nrow(.),function(x){return(d_biodesert$C_stock[which(d_biodesert$ID==.$Site_ID[x])])})),
             Productivity=(sapply(1:nrow(.),function(x){return(d_biodesert$Productivity[which(d_biodesert$ID==.$Site_ID[x])])})),
             Woody=sapply(1:nrow(.),function(x){return(d_biodesert$RWCQ[which(d_biodesert$ID==.$Site_ID[x])])}),
             Herb=sapply(1:nrow(.),function(x){return(d_biodesert$RHCQ[which(d_biodesert$ID==.$Site_ID[x])])}),
             Grass=sapply(1:nrow(.),function(x){return(d_biodesert$RGCQ[which(d_biodesert$ID==.$Site_ID[x])])}),
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

d_biodesert$Type_herb=sapply(1:nrow(d_biodesert),function(x){
  ifelse((as.numeric(d_biodesert$Dung_natives_kg_ha)/
            as.numeric(d_biodesert$Dung_all_herbivores_kg_ha))[x]==1,"Native",
         ifelse((as.numeric(d_biodesert$Dung_natives_kg_ha)/
                   as.numeric(d_biodesert$Dung_all_herbivores_kg_ha))[x]==0,"Livestock","Mixed"))
})

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

write.table(d,"../Data/Spatial_structure_grazing.csv",sep=";")


# Changing the spatial statistics to controll for vegetation cover

d=read.table("../Data/Spatial_structure_grazing.csv",sep=";")

for (x in c(2:31)){
  d_lm=melt(d,measure.vars=colnames(d)[x])
  d_lm[,c("rho_p","value")]=apply(d_lm[,c("rho_p","value")],2,scale)
  reg_cover=lm(value~rho_p,d_lm)
  d[as.numeric(names(residuals(reg_cover))),colnames(d)[x]]=residuals(reg_cover)
}

write.table(d,"../Data/Spatial_structure_grazing_control_cover.csv",sep=";")




# ---------------------- Step 2: Effects of grazing on the spatial structure  ----
## >> 1) Running mixed-effect models ---- 

Run_model_importance_aridity=function(id){
  
  list_mod=rbind(expand.grid(with_cover=c(T,F),
                             Stats=c("perim_area_scaling","fmax_psd","PL_expo",
                                     "core_area_land","core_area","Small_patches",
                                     "flow_length","mean_psd","moran_I")),
                 tibble(with_cover=F,
                        Stats=c("rho_p")))
  
  d_all=d_all2=tibble()#for keeping all information
  
  d_data=read.table("../Data/Spatial_structure_grazing.csv",sep=";")%>%
    Closer_to_normality(.)
  
  d_data[,c(1,9:29,38,40:44,47:ncol(d_data))] = 
    apply(d_data[,c(1,9:29,38,40:44,47:ncol(d_data))],2,z_tranform)
  
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
      + Sand + Org_C_v 
      + (1|Site_ID)")))
    
  } else{
    
    formula_mod=formula(formula_(paste("value ~ Long_cos + Long_sin + Lattitude + Slope + Elevation 
      + Aridity + Grazing 
      + Aridity*Grazing 
      + Type_veg 
      + Sand + Org_C_v 
      + (1|Site_ID)")))
  }
  
  save=d_data_mod
  
  d_data_mod=save
  
  
  for (grazing_intensity in c("binary","all")){
    
    if (grazing_intensity =="binary"){
      d_data_mod=save
      d_data_mod$Grazing= as.numeric(d_data_mod$Grazing!=0)
    } else{
      d_data_mod=save
    }
    
    
    d_data_mod$Grazing=as.factor(d_data_mod$Grazing)
    d_data_mod$Grazing = relevel(d_data_mod$Grazing,ref = "0")
    
    model_spa_stat  = lmer(formula_mod, d_data_mod,
                           na.action = na.fail ,REML ="FALSE")
    
    #we remove potential outliers
    rm.outliers = romr.fnc(model_spa_stat, d_data_mod, trim=2.5)
    d_data_out = rm.outliers$data
    
    #saving data
    write.table(d_data_out,paste0("../Data/Linear_models_factor/Keep_data/Data_",
                                  stat,"_",with_cover,"_aridity_",grazing_intensity,".csv"))
    
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
    
    saveRDS(model_spa_stat,#MuMIn::get.models(select_model,subset = delta<2)[[1]],
            paste0("../Data/Linear_models_factor/Keep_models/Mod_",stat,"_",
                   with_cover,"_aridity_",grazing_intensity,".rds"))
    
    summary_coef=confint(result_select)
    
    #Merge in a df
    d_all2=rbind(d_all2,tibble(Median=confint(result_select,level = 1e-9)[-1,1],
                               q1=summary_coef[-1,1], 
                               q3=summary_coef[-1,2],
                               term=rownames(summary_coef)[-1],
                               Stat=stat,
                               R2m=mean(R2$R2m),
                               R2c=mean(R2$R2c),
                               Grazing_intensity=grazing_intensity)%>%
                   add_column(., With_cover=with_cover))
    
  }  
  
  #Saving each DF
  write.table(d_all,paste0("../Data/Linear_models_factor/Importance/Importance_",
                           stat,"_",with_cover,"_aridity_",grazing_intensity,".csv"),sep=";")
  
  write.table(d_all2,paste0("../Data/Linear_models_factor/Estimator/Estimators_model_",stat
                            ,"_",with_cover,"_aridity_",grazing_intensity,".csv"),sep=";")
  
}

library(parallel)

mclapply(1:19,Run_model_importance_aridity,mc.cores = 19)

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
                                     "core_area_land","core_area","Small_patches",
                                     "flow_length","mean_psd","moran_I")),
                 tibble(with_cover=F,
                        Stats=c("rho_p")))
  
  
  d_all=d_all2=tibble()#for keeping all information
  
  d_data=read.table("../Data/Spatial_structure_grazing.csv",sep=";")%>%
    Closer_to_normality(.)
  
  d_data[,c(1,9:29,38,40:44,47:ncol(d_data))] = 
    apply(d_data[,c(1,9:29,38,40:44,47:ncol(d_data))],2,z_tranform)
  
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
  
  for (grazing_intensity in c("binary","all")){
    
    if (grazing_intensity =="binary"){
      d_data_mod=save
      d_data_mod$Grazing= as.numeric(d_data_mod$Grazing!=0)
    } else{
      d_data_mod=save
    }
    
    d_data_mod=save
    
    d_data_mod$Grazing=as.factor(d_data_mod$Grazing)
    d_data_mod$Grazing = relevel(d_data_mod$Grazing,ref = "0")
    
    model_spa_stat  = lmer(formula_mod, d_data_mod,
                           na.action = na.fail ,REML ="FALSE")
    
    #we remove potential outliers
    rm.outliers = romr.fnc(model_spa_stat, d_data_mod, trim=2.5)
    d_data_out = rm.outliers$data
    
    
    #saving data
    write.table(d_data_out,paste0("../Data/Linear_models_factor/Keep_data/Data_",
                                  stat,"_",with_cover,"_aridity_no_inter_",grazing_intensity,".csv"))
    
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
    
    saveRDS(model_spa_stat,#MuMIn::get.models(select_model,subset = delta<2)[[1]],
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
    
  }
  #Saving each DF
  write.table(d_all,paste0("../Data/Linear_models_factor/Importance/Importance_",
                           stat,"_",with_cover,"_aridity_no_inter_",grazing_intensity,".csv"),sep=";")
  write.table(d_all2,paste0("../Data/Linear_models_factor/Estimator/Estimators_model_",
                            stat,"_",with_cover,"_aridity_no_inter_",grazing_intensity,".csv"),sep=";")
  
}

library(parallel)

mclapply(1:19,Run_model_importance_aridity_no_inter,mc.cores = 19)

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
d_data[,c(1,9:29,38,40:44,47:ncol(d_data))] = 
  apply(d_data[,c(1,9:29,38,40:44,47:ncol(d_data))],2,z_tranform)
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
boot_function_lm_graz_binary = function(formula, data, indices) {
  d = data[indices,] 
  fit = lm(formula, data=d) 
  return(summary(fit)$coefficient[2,1])
}

for (with_interaction in c(T,F)){
  d_slope=tibble()
  
  for (grazing_intensity in c("binary","all")){
    
    for (k in c("perim_area_scaling","fmax_psd","PL_expo",
                "core_area_land","core_area","Small_patches",
                "flow_length","mean_psd","moran_I","rho_p")){
      
      #for each we plot the slope against the partial residuals with and without cover
      
      if (k !="rho_p"){
        
        d_data_out=read.table(paste0("../Data/Linear_models_factor/Keep_data/Data_",k,"_TRUE_aridity_",
                                     ifelse(with_interaction,"","no_inter_"),grazing_intensity,".csv"),sep=" ")
        model_spa_stat=readRDS(paste0("../Data/Linear_models_factor/Keep_models/Mod_",k,"_TRUE_aridity_",
                                      ifelse(with_interaction,"","no_inter_"),grazing_intensity,".rds"))
        
        d_data_out$Grazing = as.factor(d_data_out$Grazing)
        d_data_out$Grazing = relevel(d_data_out$Grazing,ref = "0")
        
        resid_mod=visreg::visreg(fit = model_spa_stat,xvar="Grazing",plot=F) 
        
        
        if (grazing_intensity=="binary"){
          mod_cov = boot(data=resid_mod$res, statistic=boot_function_lm_graz_binary, R=1000, formula=visregRes~Grazing)
          
          d_slope=rbind(d_slope,
                        tibble(pval=twoside_pvalue(mod_cov$t),
                               q2=median(mod_cov$t),
                               q1=quantile(mod_cov$t,.025),
                               q3=quantile(mod_cov$t,.975),
                               q1_90=quantile(mod_cov$t,.05),
                               q3_90=quantile(mod_cov$t,.95),
                               With_cover=T,Grazing_intensity=grazing_intensity,
                               ID_grazing=c("grazed"),
                               Stat=k,Driver="Grazing"
                        ))
        } else{
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
        
        resid_mod=visreg::visreg(fit = model_spa_stat,xvar="Aridity",plot=F) 
        
        
        mod_cov = boot(data=resid_mod$res, statistic=boot_function_lm_arid, R=1000, formula=visregRes~Aridity)
        
        d_slope=rbind(d_slope,
                      tibble(pval=twoside_pvalue(mod_cov$t),
                             q2=median(mod_cov$t),
                             q1=quantile(mod_cov$t,.025),
                             q3=quantile(mod_cov$t,.975),
                             q1_90=quantile(mod_cov$t,.05),
                             q3_90=quantile(mod_cov$t,.95),
                             With_cover=T,Grazing_intensity=grazing_intensity,
                             ID_grazing=c("all"),
                             Stat=k,Driver="Aridity"
                      ))
        
        
      }
      
      d_data_out=read.table(paste0("../Data/Linear_models_factor/Keep_data/Data_",k,"_FALSE_aridity_",
                                   ifelse(with_interaction,"","no_inter_"),grazing_intensity,".csv"),sep=" ")
      model_spa_stat=readRDS(paste0("../Data/Linear_models_factor/Keep_models/Mod_",k,"_FALSE_aridity_",
                                    ifelse(with_interaction,"","no_inter_"),grazing_intensity,".rds"))
      
      d_data_out$Grazing = as.factor(d_data_out$Grazing)
      d_data_out$Grazing = relevel(d_data_out$Grazing,ref = "0")
      
      resid_mod=visreg::visreg(fit = model_spa_stat,xvar="Grazing",plot=F) 
      

      if (grazing_intensity=="binary"){
        
        mod_cov = boot(data=resid_mod$res, statistic=boot_function_lm_graz_binary, R=1000, formula=visregRes~Grazing)
        
        d_slope=rbind(d_slope,
                      tibble(pval=twoside_pvalue(mod_cov$t),
                             q2=median(mod_cov$t),
                             q1=quantile(mod_cov$t,.025),
                             q3=quantile(mod_cov$t,.975),
                             q1_90=quantile(mod_cov$t,.05),
                             q3_90=quantile(mod_cov$t,.95),
                             With_cover=F,Grazing_intensity=grazing_intensity,
                             ID_grazing=c("grazed"),
                             Stat=k,Driver="Grazing"
                      ))
      } else{
        
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
      
      resid_mod=visreg::visreg(fit = model_spa_stat,xvar="Aridity",plot=F) 
      
      
      mod_cov = boot(data=resid_mod$res, statistic=boot_function_lm_arid, R=1000, formula=visregRes~Aridity)
      
      d_slope=rbind(d_slope,
                    tibble(pval=twoside_pvalue(mod_cov$t),
                           q2=median(mod_cov$t),
                           q1=quantile(mod_cov$t,.025),
                           q3=quantile(mod_cov$t,.975),
                           q1_90=quantile(mod_cov$t,.05),
                           q3_90=quantile(mod_cov$t,.95),
                           With_cover=F,Grazing_intensity=grazing_intensity,
                           ID_grazing=c("all"),
                           Stat=k,Driver="Aridity"
                    ))
      
      
      
    }
  }
  
  write.table(d_slope,paste0("../Data/Linear_models_factor/Slope_partial_residuals_aridity",
                             ifelse(with_interaction,"","_no_inter"),
                             ".csv"),sep=";")
}

# ---------------------- Step 2bis: Effects of grazing on the spatial structure cover control residuals ----
## >> 1) Running mixed-effect model ----

Run_model_importance_aridity=function(id){
  
  list_mod=rbind(expand.grid(with_cover=c(F),
                             Stats=c("perim_area_scaling","fmax_psd","PL_expo",
                                     "core_area_land","core_area","Small_patches",
                                     "flow_length","mean_psd","moran_I")),
                 tibble(with_cover=F,
                        Stats=c("rho_p")))
  
  d_all=d_all2=tibble()#for keeping all information
  
  d_data=read.table("../Data/Spatial_structure_grazing.csv",sep=";")%>%
    Closer_to_normality(.)
  
  d_data[,c(1,9:29,38,40:44,47:ncol(d_data))] = 
    apply(d_data[,c(1,9:29,38,40:44,47:ncol(d_data))],2,z_tranform)
  
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
      + Sand + Org_C_v
      + (1|Site_ID)")))
    
  } else{
    
    formula_mod=formula(formula_(paste("value ~ Long_cos + Long_sin + Lattitude + Slope + Elevation 
      + Aridity + Grazing 
      + Aridity*Grazing 
      + Type_veg 
      + Sand + Org_C_v 
      + (1|Site_ID)")))
  }
  
  save=d_data_mod
  
  d_data_mod=save
  
  
  for (grazing_intensity in c("binary","all")){
    
    if (grazing_intensity =="binary"){
      d_data_mod=save
      d_data_mod$Grazing= as.numeric(d_data_mod$Grazing!=0)
    } else{
      d_data_mod=save
    }
    
    
    #controlling for cover 
    model_spa_stat  = lm(value~rho_p, d_data_mod)
    
    d_data_mod$value=residuals(model_spa_stat)
    
    d_data_mod$Grazing=as.factor(d_data_mod$Grazing)
    d_data_mod$Grazing = relevel(d_data_mod$Grazing,ref = "0")
    
    model_spa_stat  = lmer(formula_mod, d_data_mod,
                           na.action = na.fail ,REML ="FALSE")
    
    #we remove potential outliers
    rm.outliers = romr.fnc(model_spa_stat, d_data_mod, trim=2.5)
    d_data_out = rm.outliers$data
    
    #saving data
    write.table(d_data_out,paste0("../Data/Linear_models_factor_cover_control/Keep_data/Data_",
                                  stat,"_",with_cover,"_aridity_",grazing_intensity,".csv"))
    
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
    
    saveRDS(model_spa_stat,#MuMIn::get.models(select_model,subset = delta<2)[[1]],
            paste0("../Data/Linear_models_factor_cover_control/Keep_models/Mod_",stat,"_",
                   with_cover,"_aridity_",grazing_intensity,".rds"))
    
    summary_coef=confint(result_select)
    
    #Merge in a df
    d_all2=rbind(d_all2,tibble(Median=confint(result_select,level = 1e-9)[-1,1],
                               q1=summary_coef[-1,1], 
                               q3=summary_coef[-1,2],
                               term=rownames(summary_coef)[-1],
                               Stat=stat,
                               R2m=mean(R2$R2m),
                               R2c=mean(R2$R2c),
                               Grazing_intensity=grazing_intensity)%>%
                   add_column(., With_cover=with_cover))
    
  }  
  
  #Saving each DF
  write.table(d_all,paste0("../Data/Linear_models_factor_cover_control/Importance/Importance_",
                           stat,"_",with_cover,"_aridity_",grazing_intensity,".csv"),sep=";")
  
  write.table(d_all2,paste0("../Data/Linear_models_factor_cover_control/Estimator/Estimators_model_",stat
                            ,"_",with_cover,"_aridity_",grazing_intensity,".csv"),sep=";")
  
}

library(parallel)

mclapply(1:10,Run_model_importance_aridity,mc.cores = 10)

d_all=d_all2=tibble()
for (k in list.files("../Data/Linear_models_factor_cover_control/Importance",paste0("aridity.csv"))){
  d=read.table(paste0("../Data/Linear_models_factor_cover_control/Importance/",k),sep=";")
  if (nrow(d)!=0){
    d_all=rbind(d_all,d)
  }
}
for (k in list.files("../Data/Linear_models_factor_cover_control/Estimator",paste0("aridity.csv"))){
  d2=read.table(paste0("../Data/Linear_models_factor_cover_control/Estimator/",k),sep=";")
  if (nrow(d2)!=0){
    d_all2=rbind(d_all2,d2)
  }
}
write.table(d_all,paste0("../Data/Linear_models_factor_cover_control/Importance_aridity_factor.csv"),sep=";")
write.table(d_all2,paste0("../Data/Linear_models_factor_cover_control/Estimators_model_aridity_factor.csv"),sep=";")


Run_model_importance_aridity_no_inter=function(id){
  
  list_mod=expand.grid(with_cover=c(F),
                       Stats=c("perim_area_scaling","fmax_psd","PL_expo",
                               "core_area_land","core_area","Small_patches",
                               "flow_length","mean_psd","moran_I"))
  
  d_all=d_all2=tibble()#for keeping all information
  
  d_data=read.table("../Data/Spatial_structure_grazing.csv",sep=";")%>%
    Closer_to_normality(.)
  
  d_data[,c(1,9:29,38,40:44,47:ncol(d_data))] = 
    apply(d_data[,c(1,9:29,38,40:44,47:ncol(d_data))],2,z_tranform)
  
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
  
  for (grazing_intensity in c("binary","all")){
    
    if (grazing_intensity =="binary"){
      d_data_mod=save
      d_data_mod$Grazing= as.numeric(d_data_mod$Grazing!=0)
    } else{
      d_data_mod=save
    }
    
    d_data_mod=save
    
    #controlling for cover 
    model_spa_stat  = lm(value~rho_p, d_data_mod)
    
    d_data_mod$value=residuals(model_spa_stat)
    
    
    d_data_mod$Grazing=as.factor(d_data_mod$Grazing)
    d_data_mod$Grazing = relevel(d_data_mod$Grazing,ref = "0")
    
    model_spa_stat  = lmer(formula_mod, d_data_mod,
                           na.action = na.fail ,REML ="FALSE")
    
    #we remove potential outliers
    rm.outliers = romr.fnc(model_spa_stat, d_data_mod, trim=2.5)
    d_data_out = rm.outliers$data
    
    
    #saving data
    write.table(d_data_out,paste0("../Data/Linear_models_factor_cover_control/Keep_data/Data_",
                                  stat,"_",with_cover,"_aridity_no_inter_",grazing_intensity,".csv"))
    
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
    
    saveRDS(model_spa_stat,#MuMIn::get.models(select_model,subset = delta<2)[[1]],
            paste0("../Data/Linear_models_factor_cover_control/Keep_models/Mod_",stat,
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
    
  }
  #Saving each DF
  write.table(d_all,paste0("../Data/Linear_models_factor_cover_control/Importance/Importance_",
                           stat,"_",with_cover,"_aridity_no_inter_",grazing_intensity,".csv"),sep=";")
  write.table(d_all2,paste0("../Data/Linear_models_factor_cover_control/Estimator/Estimators_model_",
                            stat,"_",with_cover,"_aridity_no_inter_",grazing_intensity,".csv"),sep=";")
  
}

library(parallel)

mclapply(1:9,Run_model_importance_aridity_no_inter,mc.cores = 9)

d_all=d_all2=tibble()
for (k in list.files("../Data/Linear_models_factor_cover_control/Importance",paste0("aridity_no_inter"))){
  d=read.table(paste0("../Data/Linear_models_factor_cover_control/Importance/",k),sep=";")
  if (nrow(d)!=0){
    d_all=rbind(d_all,d)
  }
}
for (k in list.files("../Data/Linear_models_factor_cover_control/Estimator",paste0("aridity_no_inter"))){
  d2=read.table(paste0("../Data/Linear_models_factor_cover_control/Estimator/",k),sep=";")
  if (nrow(d)!=0){
    d_all2=rbind(d_all2,d2)
  }
}
write.table(d_all,paste0("../Data/Linear_models_factor_cover_control/Importance_aridity_no_inter_factor.csv"),sep=";")
write.table(d_all2,paste0("../Data/Linear_models_factor_cover_control/Estimators_model_aridity_no_inter_factor.csv"),sep=";")


## >> 2) Analyzing the residuals ---- 

d_data=read.table("../Data/Spatial_structure_grazing.csv",sep=";")%>%
  Closer_to_normality(.)
d_data[,c(1,9:29,38,40:44,47:ncol(d_data))] = 
  apply(d_data[,c(1,9:29,38,40:44,47:ncol(d_data))],2,z_tranform)
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
boot_function_lm_graz_binary = function(formula, data, indices) {
  d = data[indices,] 
  fit = lm(formula, data=d) 
  return(summary(fit)$coefficient[2,1])
}

for (with_interaction in c(F,T)){
  d_slope=tibble()
  
  for (grazing_intensity in c("binary","all")){
    
    for (k in c("perim_area_scaling","fmax_psd","PL_expo",
                "core_area_land","core_area","Small_patches",
                "flow_length","mean_psd","moran_I")){
      
      #for each we plot the slope against the partial residuals with and without cover
      
      d_data_out=read.table(paste0("../Data/Linear_models_factor_cover_control/Keep_data/Data_",k,"_FALSE_aridity_",
                                   ifelse(with_interaction,"","no_inter_"),grazing_intensity,".csv"),sep=" ")
      model_spa_stat=readRDS(paste0("../Data/Linear_models_factor_cover_control/Keep_models/Mod_",k,"_FALSE_aridity_",
                                    ifelse(with_interaction,"","no_inter_"),grazing_intensity,".rds"))
      
      d_data_out$Grazing = as.factor(d_data_out$Grazing)
      d_data_out$Grazing = relevel(d_data_out$Grazing,ref = "0")
      
      resid_mod=visreg::visreg(fit = model_spa_stat,xvar="Grazing",plot=F) 

      if (grazing_intensity=="binary"){
        
        mod_cov = boot(data=resid_mod$res, statistic=boot_function_lm_graz_binary, R=1000, formula=visregRes~Grazing)
        
        d_slope=rbind(d_slope,
                      tibble(pval=twoside_pvalue(mod_cov$t),
                             q2=median(mod_cov$t),
                             q1=quantile(mod_cov$t,.025),
                             q3=quantile(mod_cov$t,.975),
                             q1_90=quantile(mod_cov$t,.05),
                             q3_90=quantile(mod_cov$t,.95),
                             With_cover=F,Grazing_intensity=grazing_intensity,
                             ID_grazing=c("grazed"),
                             Stat=k,Driver="Grazing"
                      ))
      } else{
        
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
      
      
      
      resid_mod=visreg::visreg(fit = model_spa_stat,xvar="Aridity",plot=F) 
      
      
      mod_cov = boot(data=resid_mod$res, statistic=boot_function_lm_arid, R=1000, formula=visregRes~Aridity)
      
      d_slope=rbind(d_slope,
                    tibble(pval=twoside_pvalue(mod_cov$t),
                           q2=median(mod_cov$t),
                           q1=quantile(mod_cov$t,.025),
                           q3=quantile(mod_cov$t,.975),
                           q1_90=quantile(mod_cov$t,.05),
                           q3_90=quantile(mod_cov$t,.95),
                           With_cover=F,Grazing_intensity=grazing_intensity,
                           ID_grazing=c("all"),
                           Stat=k,Driver="Aridity"
                    ))
      
      
    }
  }
  
  write.table(d_slope,paste0("../Data/Linear_models_factor_cover_control/Slope_partial_residuals_aridity",
                             ifelse(with_interaction,"","_no_inter"),
                             ".csv"),sep=";")
}




# ---------------------- Step 3: Processing traits----
## >> 1) Processing cover quadrats: writing a file per site ----

dir.create("../Data/Traits",showWarnings = F)
dir.create("../Data/Traits/Cover_quadrats",showWarnings = F)

list_sheets=readxl::excel_sheets("../Data/Traits/Cover_quadrats_Biodesert.xlsx")
list_cover_quadrat=list()

index=1
for (k in 1:length(list_sheets)){
  
  name_sheet_k=list_sheets[k] #name country
  
  quadrats_k=readxl::read_xlsx("../Data/Traits/Cover_quadrats_Biodesert.xlsx",sheet = name_sheet_k)
  
  rows_quadrat_sites=as_tibble(matrix(sort(c(1,which(is.na(quadrats_k[,1])),
                                             which(is.na(quadrats_k[,1]))+1,nrow(quadrats_k))),
                                      ncol=2,byrow = T))
  
  for (site_k in 1:nrow(rows_quadrat_sites)){ #for each site/grazing pressure, extract the cover
    
    quadrats_k_site = quadrats_k[rows_quadrat_sites$V1[site_k]:rows_quadrat_sites$V2[site_k],] #taking a site
    quadrats_k_site=quadrats_k_site[which(!is.na(quadrats_k_site[,1])),1:101] #removing the NA (space in the sheet between sites) and the sums of speces cover
    
    if (site_k>1){
      colnames(quadrats_k_site)=paste0(quadrats_k_site[1,]) # changing columns names
      quadrats_k_site=quadrats_k_site[-1,] #removing the column names from the data frame
    }
    
    Plot_number=as.numeric(gsub("\\D", "", colnames(quadrats_k_site)[1]))
    Site=gsub(paste0(" ",Plot_number),"",gsub(paste0(" ",Plot_number),"",colnames(quadrats_k_site)[1]))
    Country = name_sheet_k
    
    if (!is.na(Plot_number)){
      write.table(quadrats_k_site,paste0("../Data/Traits/Cover_quadrats/",Country,"_",Site,"_",Plot_number,".csv"),sep=";")
    }
    
    list_cover_quadrat[[index]]=quadrats_k_site
    names(list_cover_quadrat)[index]= paste0(Country,"_",Site,"_",Plot_number)
    index=index+1
  }
}
saveRDS(list_cover_quadrat[-grep("NA", names(list_cover_quadrat))],"../Data/Traits/All_quadrats.rds")


## >> 2) Computing CWM, CW variance and FD (functional divergence) ----

d_CWT_FD=d_traits%>%
  filter(., !is.na(Cover))%>%#remove absent species
  mutate(., Cover=Cover/sum(Cover))%>%#relative species abundance
  melt(., measure.vars=colnames(d_traits)[7:16])%>%#for each trait
  dplyr::group_by(.,Country,Site,Plot,variable)%>% 
  dplyr::summarise(., .groups = "keep",
                   CWM = sum(Cover*value,na.rm = T), # community weighted mean
                   CW_Var = sum(Cover*((value-mean(value,na.rm=T))**2),na.rm = T), # community weighted variance
                   FD = sum(Cover * (abs(value-sum(Cover*value,na.rm = T))/
                                       (sum(abs(value-sum(Cover*value,na.rm = T)),na.rm = T))),na.rm = T) # functional divergence
                   
  )

d_biodesert$Complete_name=paste0(d_biodesert$COU,"_",d_biodesert$SITE,"_",d_biodesert$PLOT)
d_CWT_FD$Complete_name=paste0(d_CWT_FD$Country,"_",d_CWT_FD$Site,"_",d_CWT_FD$Plot)

d_CWT_FD=d_CWT_FD%>%
  add_column(.,
             Graz=unlist(sapply(1:nrow(.),function(x){
               if (length(d_biodesert$GRAZ[which(d_biodesert$Complete_name==.$Complete_name[x])])==0){
                 return(NA)
                 print(x)
               }else{
                 return(d_biodesert$GRAZ[which(d_biodesert$Complete_name==.$Complete_name[x])])
               }
             })),
             ID=unlist(sapply(1:nrow(.),function(x){
               if (length(d_biodesert$ID[which(d_biodesert$Complete_name==.$Complete_name[x])])==0){
                 return(NA)
               }else{
                 return(d_biodesert$ID[which(d_biodesert$Complete_name==.$Complete_name[x])])
               }
             })))

write.table(d_CWT_FD,"../Data/Traits/CWM_FD_sites.csv",sep=";")



## >> 3) Within quadrat scale: trait dispersion or aggregation (for each trait) ----

dir.create("../Data/Traits/PwD",showWarnings = F)
dir.create("../Data/Traits/QWT_random",showWarnings = F)

Run_pair_weighted_trait_distance=function(plot_id){
  
  #Prior to analyses: standardize each trait 
  d_traits_scaled=d_traits%>%
    filter(., !is.na(Cover))%>%#remove absent species
    dplyr::group_by(.,Country,Site,Plot)%>% 
    dplyr::mutate(., LL=(LL-mean(LL,na.rm=T))/(max(LL,na.rm=T)-min(LL,na.rm = T)),
                  SLA=(SLA-mean(SLA,na.rm=T))/(max(SLA,na.rm=T)-min(SLA,na.rm = T)),
                  LDMC=(LDMC-mean(LDMC,na.rm=T))/(max(LDMC,na.rm=T)-min(LDMC,na.rm = T)),
                  LA=(LA-mean(LA,na.rm=T))/(max(LA,na.rm=T)-min(LA,na.rm = T)),
                  MaxH=(MaxH-mean(MaxH,na.rm=T))/(max(MaxH,na.rm=T)-min(MaxH,na.rm = T)),
                  MaxLS=(MaxLS-mean(MaxLS,na.rm=T))/(max(MaxLS,na.rm=T)-min(MaxLS,na.rm = T)),
                  Maxvolume=(Maxvolume-mean(Maxvolume,na.rm=T))/(max(Maxvolume,na.rm=T)-min(Maxvolume,na.rm = T)),
                  Phenolics=(Phenolics-mean(Phenolics,na.rm=T))/(max(Phenolics,na.rm=T)-min(Phenolics,na.rm = T)),
                  LNC=(LNC-mean(LNC,na.rm=T))/(max(LNC,na.rm=T)-min(LNC,na.rm = T)),
                  LCC=(LCC-mean(LCC,na.rm=T))/(max(LCC,na.rm=T)-min(LCC,na.rm = T)),
                  Cover=Cover,
                  Genus=Genus,
                  Species=Species
    )
  
  
  list_quadrat_sites=list.files("../Data/Traits/Cover_quadrats/")
  n_perm=199
  
  quadrat_k=read.table(paste0("../Data/Traits/Cover_quadrats/",list_quadrat_sites[plot_id]),sep=";")
  
  country=strsplit(list_quadrat_sites[plot_id],"_")[[1]][1]
  site=strsplit(list_quadrat_sites[plot_id],"_")[[1]][2]
  plot=as.numeric(gsub(".csv","",strsplit(list_quadrat_sites[plot_id],"_")[[1]][3]))
  
  mat_quadrat=as.matrix(quadrat_k[,-1])
  colnames(mat_quadrat)=paste0("Quadrat_",1:100)
  rownames(mat_quadrat)=quadrat_k[,1]
  
  traits_site=d_traits_scaled%>%
    dplyr::filter(.,Country==country,Site==site,Plot==plot)
  
  #Randomization of the quadrats
  Randomization_quadrat=vegan::permatfull(mat_quadrat,times = n_perm)
  
  d_PwD=tibble()
  for (random_id in 1:length(Randomization_quadrat$perm)){ # for each of the randomization
    
    quadrat_k_id=Randomization_quadrat$perm[[random_id]]
    
    d_QWT=tibble()
    for (column_id in 1:ncol(quadrat_k_id)){# for each quadrat, compute the Quadrat weighted trait
      d_QWT=rbind(d_QWT,Get_CWT_traits(abundance = quadrat_k_id[,column_id],traits_site = traits_site)%>%
                    add_column(., Quadrat=column_id,Randomization_id=random_id))
    }
    
    d_PwD=rbind(d_PwD,d_QWT%>%  #for the randomized QWT compute pairwise trait distance
      dplyr::group_by(., variable,Country,Site,Plot,Randomization_id)%>%
      dplyr::do(.,Compute_pairwise_trait_distance(.$CWM)))
    
  write.table(d_QWT,paste0("../Data/Traits/QWT_random/QWT_random_Randomization_",random_id,"_",list_quadrat_sites[plot_id]),sep=";")
  }
  
  #Observed traits metrics
  
  d_QWT_obs=tibble()
  for (column_id in 1:ncol(mat_quadrat)){# for each quadrat, compute the Quadrat weighted trait (here on observed species abundance)
    d_QWT_obs=rbind(d_QWT_obs,Get_CWT_traits(abundance = mat_quadrat[,column_id],traits_site = traits_site)%>%
                      add_column(., Quadrat=column_id,Randomization_id="Observed"))  
  }
  
  PwD_obs=d_QWT_obs%>% #for the observed QWT compute pairwise trait distance
    dplyr::group_by(., variable,Country,Site,Plot,Randomization_id)%>%
    dplyr::do(.,Compute_pairwise_trait_distance(.$CWM))
  
  PwD_random=d_PwD%>%
    dplyr::ungroup(.)%>%  #then extract the mean and confidence interval across randomization 
    dplyr::group_by(., variable,Country,Site,Plot)%>%
    dplyr::summarise(., .groups = "keep",
                     mean_PwD=mean(PwD),
                     q025_PwD=quantile(PwD,.025),
                     q975_PwD=quantile(PwD,.975),
                     q05_PwD=quantile(PwD,.05),
                     q95_PwD=quantile(PwD,.95)
    )
  
  PwD_tot=PwD_obs%>%
    add_column(.,
               mean_PwD_random=PwD_random$mean_PwD,
               q025_PwD_random=PwD_random$q025_PwD,
               q975_PwD_random=PwD_random$q975_PwD,
               q05_PwD_random=PwD_random$q05_PwD,
               q95_PwD_random=PwD_random$q95_PwD)%>%
    dplyr::rename(., PwD_obs=PwD)
  
  
  write.table(PwD_tot,paste0("../Data/Traits/PwD_",list_quadrat_sites[plot_id]),sep=";")
}

list_quadrat_sites=list.files("../Data/Traits/Cover_quadrats/")

library(parallel)
mclapply(1:length(list_quadrat_sites),Run_pair_weighted_trait_distance,mc.cores = 20)




# ---------------------- Step 4: Adding community composition, facilitation and fertility effects ----

#Adding these variables to the data

d_sem=Perform_PCA_spatial_struc(d_data,plot=T)
k="Struct1"
#We first control for all covariates and extract the residuals
model_lmer=lm("value ~ Long_cos + Long_sin + Lattitude + Elevation",
              data = d_sem%>%melt(., measure.vars=k)%>%filter(., !is.na(value)),
              na.action = na.fail)

resid_model=residuals(model_lmer) #extract residuals of the distance to the tipping point after controlling for all covariates

d_sem=d_sem[as.numeric(names(resid_model)),]%>% #add it to the dataframe 
  add_column(., Resid_mod=resid_model)

SEM_distance=psem(
  lmer(Resid_mod ~ Grazing + Woody + Forage_Quality + Fertility +Aridity+ (1|Site_ID), d_sem),
  lm(Woody ~ Grazing+Aridity, d_sem),
  lm(Forage_Quality ~ Aridity + Grazing, d_sem),
  lm(Fertility ~ Grazing+Aridity, d_sem))

summary(SEM_distance)


# ---------------------- Step 5: Changing grazing reference level ----
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
  
  d_data[,c(1,9:29,38,40:44,47:ncol(d_data))] = 
    apply(d_data[,c(1,9:29,38,40:44,47:ncol(d_data))],2,z_tranform)
  
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
d_data[,c(1,9:29,38,40:44,47:ncol(d_data))] = 
  apply(d_data[,c(1,9:29,38,40:44,47:ncol(d_data))],2,z_tranform)
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
  
  d_data[,c(1,9:29,38,40:44,47:ncol(d_data))] = 
    apply(d_data[,c(1,9:29,38,40:44,47:ncol(d_data))],2,z_tranform)
  
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
d_data[,c(1,9:29,38,40:44,47:ncol(d_data))] = 
  apply(d_data[,c(1,9:29,38,40:44,47:ncol(d_data))],2,z_tranform)
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
  
  d_data[,c(1,9:29,38,40:44,47:ncol(d_data))] = 
    apply(d_data[,c(1,9:29,38,40:44,47:ncol(d_data))],2,z_tranform)
  
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
d_data[,c(1,9:29,38,40:44,47:ncol(d_data))] = 
  apply(d_data[,c(1,9:29,38,40:44,47:ncol(d_data))],2,z_tranform)
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



# ---------------------- Step 6: Using the Schneider model with different spatial structures ----
library(chouca)

#high cover, high aggregation
d=tibble()
CA=ca_library("aridvege",
              parms = list(r=0,
                           f=.9,
                           delta=.1,
                           c=.2,
                           d=.1,
                           m0=0.001,
                           g0=.3,
                           b=.49,
                           pr=1
              ))
im = generate_initmat(CA, c(0.2, 0.3,.5), nrow = 100, ncol = 100)
out=run_camodel(CA, im, times = seq(0,2000),control = list(save_snapshots_every=1))
landscape_sim11=matrix(out$output$snapshots[[which(out$output$covers[,4]==.55)[1]]]=="VEGE",100,100)
image(landscape_sim11)
d=rbind(d,Get_sumstat(landscape_sim11,compute_KS = F))

#high cover, medium aggregation
CA=ca_library("aridvege",
              parms = list(r=0,
                           f=.9,
                           delta=.1,
                           c=.2,
                           d=.1,
                           m0=0.05,
                           g0=.14,
                           b=.46,
                           pr=1
              ))
im = generate_initmat(CA, c(0.2, 0.3,.5), nrow = 100, ncol = 100)
out=run_camodel(CA, im, times = seq(0,4000),control = list(save_snapshots_every=1))
landscape_sim21=matrix(out$output$snapshots[[which(out$output$covers[,4]==.55)[length(which(out$output$covers[,4]==.55))]]]=="VEGE",100,100)
image(landscape_sim21)
d=rbind(d,Get_sumstat(landscape_sim21,compute_KS = F))


#high cover, low aggregation

CA=ca_library("aridvege",
              parms = list(r=0,
                           f=.9,
                           delta=.1,
                           c=.2,
                           d=.1,
                           m0=0.1,
                           g0=0,
                           b=.42,
                           pr=1
              ))
im = generate_initmat(CA, c(0.2, 0.3,.5), nrow = 100, ncol = 100)
out=run_camodel(CA, im, times = seq(0,2000),control = list(save_snapshots_every=1))
landscape_sim31=matrix(out$output$snapshots[[which(out$output$covers[,4]==.55)[length(which(out$output$covers[,4]==.55))]]]=="VEGE",100,100)
image(landscape_sim31)
d=rbind(d,Get_sumstat(landscape_sim31,compute_KS = F))


#medium cover, high aggregation

CA=ca_library("aridvege",
              parms = list(r=0,
                           f=.9,
                           delta=.1,
                           c=.2,
                           d=.1,
                           m0=0.001,
                           g0=.319,
                           b=.5,
                           pr=1
              ))
im = generate_initmat(CA, c(0.2, 0.3,.5), nrow = 100, ncol = 100)
out=run_camodel(CA, im, times = seq(0,4000),control = list(save_snapshots_every=1))
landscape_sim12=matrix(out$output$snapshots[[which(out$output$covers[,4]==.35)[length(which(out$output$covers[,4]==.35))]]]=="VEGE",100,100)
image(landscape_sim12)
d=rbind(d,Get_sumstat(landscape_sim12,compute_KS = F))

#medium cover, medium aggregation
CA=ca_library("aridvege",
              parms = list(r=0,
                           f=.9,
                           delta=.1,
                           c=.2,
                           d=.1,
                           m0=0.05,
                           g0=.202,
                           b=.5,
                           pr=1
              ))
im = generate_initmat(CA, c(0.2, 0.3,.5), nrow = 100, ncol = 100)
out=run_camodel(CA, im, times = seq(0,4000),control = list(save_snapshots_every=1))
landscape_sim22=matrix(out$output$snapshots[[which(out$output$covers[,4]==.35)[length(which(out$output$covers[,4]==.35))]]]=="VEGE",100,100)
image(landscape_sim22)
d=rbind(d,Get_sumstat(landscape_sim22,compute_KS = F))

#medium cover, low aggregation
CA=ca_library("aridvege",
              parms = list(r=0,
                           f=.9,
                           delta=.1,
                           c=.2,
                           d=.1,
                           m0=0.17,
                           g0=0,
                           b=.5,
                           pr=1
              ))
im = generate_initmat(CA, c(0.2, 0.3,.5), nrow = 100, ncol = 100)
out=run_camodel(CA, im, times = seq(0,2000),control = list(save_snapshots_every=1))
landscape_sim32=matrix(out$output$snapshots[[which(out$output$covers[,4]==.35)[length(which(out$output$covers[,4]==.35))]]]=="VEGE",100,100)
image(landscape_sim32)
d=rbind(d,Get_sumstat(landscape_sim32,compute_KS = F))


#low cover, high aggregation

CA=ca_library("aridvege",
              parms = list(r=0,
                           f=.9,
                           delta=.1,
                           c=.2,
                           d=.1,
                           m0=0.001,
                           g0=.325,
                           b=.5,
                           pr=1
              ))
im = generate_initmat(CA, c(0.2, 0.3,.5), nrow = 100, ncol = 100)
out=run_camodel(CA, im, times = seq(0,4000),control = list(save_snapshots_every=1))
landscape_sim13=matrix(out$output$snapshots[[which(out$output$covers[,4]==.17)[length(which(out$output$covers[,4]==.17))]]]=="VEGE",100,100)
image(landscape_sim13)
d=rbind(d,Get_sumstat(landscape_sim13,compute_KS = F))

# #low cover, medium aggregation
CA=ca_library("aridvege",
              parms = list(r=0,
                           f=.9,
                           delta=.1,
                           c=.2,
                           d=.1,
                           m0=0.05,
                           g0=.20802,
                           b=.5,
                           pr=1
              ))
im = generate_initmat(CA, c(0.2, 0.3,.5), nrow = 100, ncol = 100)
out=run_camodel(CA, im, times = seq(0,4000),control = list(save_snapshots_every=1))
landscape_sim23=matrix(out$output$snapshots[[which(out$output$covers[,4] > .1695 & out$output$covers[,4]<.1705)[length(which(out$output$covers[,4] > .1695 & out$output$covers[,4]<.1705))]]]=="VEGE",100,100)
image(landscape_sim23)
d=rbind(d,Get_sumstat(landscape_sim23,compute_KS = F))

#low cover, low aggregation
CA=ca_library("aridvege",
              parms = list(r=0,
                           f=.9,
                           delta=.1,
                           c=.2,
                           d=.1,
                           m0=0.18,
                           g0=0,
                           b=.5,
                           pr=1
              ))
im = generate_initmat(CA, c(0.2, 0.3,.5), nrow = 100, ncol = 100)
out=run_camodel(CA, im, times = seq(0,2000),control = list(save_snapshots_every=1))
landscape_sim33=matrix(out$output$snapshots[[which(out$output$covers[,4]==.17)[length(which(out$output$covers[,4]==.17))]]]=="VEGE",100,100)
image(landscape_sim33)
d=rbind(d,Get_sumstat(landscape_sim33,compute_KS = F))



list_landscape=list(landscape_sim11,landscape_sim12,landscape_sim13,
                    landscape_sim21,landscape_sim22,landscape_sim23,
                    landscape_sim31,landscape_sim32,landscape_sim33)

saveRDS(object = list_landscape,"../Data/Simulations/Minimal_examples_stats_landscapes.rds")

write.table(d,"../Data/Simulations/Minimal_examples_stats.csv",sep=";")


# ---------------------- Step 7: Sensitivity spatial resolution ----
## >> 1) Running mixed-effect model ----

Run_model_importance_aridity=function(id){
  
  list_mod=rbind(expand.grid(with_cover=c(T,F),
                             Stats=c("perim_area_scaling","fmax_psd","PL_expo",
                                     "core_area_land","core_area","Small_patches",
                                     "flow_length","mean_psd","moran_I")),
                 tibble(with_cover=F,
                        Stats=c("rho_p")))
  
  d_all=d_all2=tibble()#for keeping all information
  
  d_data=read.table("../Data/Spatial_structure_grazing.csv",sep=";")%>%
    Closer_to_normality(.)
  
  d_data[,c(1,9:29,38,40:44,47:ncol(d_data))] = 
    apply(d_data[,c(1,9:29,38,40:44,47:ncol(d_data))],2,z_tranform)
  
  d_data=Perform_PCA_spatial_struc(d_data)
  
  stat=list_mod$Stats[id]
  
  d_data_mod=d_data%>%melt(., measure.vars=stat)%>%
    filter(., !is.na(value))
  
  with_cover=list_mod$with_cover[id]
  
  if (with_cover){
    
    formula_mod=formula(formula_(paste("value ~ Long_cos + Long_sin + Lattitude + Slope + Elevation 
      + Aridity + Grazing + Resolution
      + Aridity*Grazing 
      + rho_p + Type_veg 
      + Sand + Org_C_v 
      + (1|Site_ID)")))
    
  } else{
    
    formula_mod=formula(formula_(paste("value ~ Long_cos + Long_sin + Lattitude + Slope + Elevation 
      + Aridity + Grazing + Resolution
      + Aridity*Grazing 
      + Type_veg 
      + Sand + Org_C_v 
      + (1|Site_ID)")))
  }
  
  save=d_data_mod
  
  d_data_mod=save
  
  
  for (grazing_intensity in c("binary","all")){
    
    if (grazing_intensity =="binary"){
      d_data_mod=save
      d_data_mod$Grazing= as.numeric(d_data_mod$Grazing!=0)
    } else{
      d_data_mod=save
    }
    
    
    #controlling for cover 
    model_spa_stat  = lm(value~rho_p, d_data_mod)
    
    d_data_mod$value=residuals(model_spa_stat)
    
    d_data_mod$Grazing=as.factor(d_data_mod$Grazing)
    d_data_mod$Grazing = relevel(d_data_mod$Grazing,ref = "0")
    
    model_spa_stat  = lmer(formula_mod, d_data_mod,
                           na.action = na.fail ,REML ="FALSE")
    
    #we remove potential outliers
    rm.outliers = romr.fnc(model_spa_stat, d_data_mod, trim=2.5)
    d_data_out = rm.outliers$data
    
    #saving data
    write.table(d_data_out,paste0("../Data/Linear_models_resolution/Keep_data/Data_",
                                  stat,"_",with_cover,"_aridity_",grazing_intensity,".csv"))
    
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
    
    saveRDS(model_spa_stat,#MuMIn::get.models(select_model,subset = delta<2)[[1]],
            paste0("../Data/Linear_models_resolution/Keep_models/Mod_",stat,"_",
                   with_cover,"_aridity_",grazing_intensity,".rds"))
    
    summary_coef=confint(result_select)
    
    #Merge in a df
    d_all2=rbind(d_all2,tibble(Median=confint(result_select,level = 1e-9)[-1,1],
                               q1=summary_coef[-1,1], 
                               q3=summary_coef[-1,2],
                               term=rownames(summary_coef)[-1],
                               Stat=stat,
                               R2m=mean(R2$R2m),
                               R2c=mean(R2$R2c),
                               Grazing_intensity=grazing_intensity)%>%
                   add_column(., With_cover=with_cover))
    
  }  
  
  #Saving each DF
  write.table(d_all,paste0("../Data/Linear_models_resolution/Importance/Importance_",
                           stat,"_",with_cover,"_aridity_",grazing_intensity,".csv"),sep=";")
  
  write.table(d_all2,paste0("../Data/Linear_models_resolution/Estimator/Estimators_model_",stat
                            ,"_",with_cover,"_aridity_",grazing_intensity,".csv"),sep=";")
  
}

library(parallel)

mclapply(1:19,Run_model_importance_aridity,mc.cores = 19)

d_all=d_all2=tibble()
for (k in list.files("../Data/Linear_models_resolution/Importance",paste0("aridity.csv"))){
  d=read.table(paste0("../Data/Linear_models_resolution/Importance/",k),sep=";")
  if (nrow(d)!=0){
    d_all=rbind(d_all,d)
  }
}
for (k in list.files("../Data/Linear_models_resolution/Estimator",paste0("aridity.csv"))){
  d2=read.table(paste0("../Data/Linear_models_resolution/Estimator/",k),sep=";")
  if (nrow(d2)!=0){
    d_all2=rbind(d_all2,d2)
  }
}
write.table(d_all,paste0("../Data/Linear_models_resolution/Importance_aridity_factor.csv"),sep=";")
write.table(d_all2,paste0("../Data/Linear_models_resolution/Estimators_model_aridity_factor.csv"),sep=";")


Run_model_importance_aridity_no_inter=function(id){
  
  list_mod=expand.grid(with_cover=c(T,F),
                       Stats=c("perim_area_scaling","fmax_psd","PL_expo",
                               "core_area_land","core_area","Small_patches",
                               "flow_length","mean_psd","moran_I"))
  
  
  d_all=d_all2=tibble()#for keeping all information
  
  d_data=read.table("../Data/Spatial_structure_grazing.csv",sep=";")%>%
    Closer_to_normality(.)
  
  d_data[,c(1,9:29,38,40:44,47:ncol(d_data))] = 
    apply(d_data[,c(1,9:29,38,40:44,47:ncol(d_data))],2,z_tranform)
  
  d_data=Perform_PCA_spatial_struc(d_data)
  
  stat=list_mod$Stats[id]
  
  d_data_mod=d_data%>%melt(., measure.vars=stat)%>%
    filter(., !is.na(value))
  
  with_cover=list_mod$with_cover[id]
  
  if (with_cover){
    
    formula_mod=formula(formula_(paste("value ~ Long_cos + Long_sin + Lattitude + Slope + Elevation 
      + Aridity + Grazing + Resolution
      + rho_p + Type_veg 
      + Sand + Org_C_v 
      + (1|Site_ID)")))
    
  } else{
    
    formula_mod=formula(formula_(paste("value ~ Long_cos + Long_sin + Lattitude + Slope + Elevation 
      + Aridity + Grazing + Resolution
      + Type_veg
      + Sand + Org_C_v 
      + (1|Site_ID)")))
  }
  
  save=d_data_mod
  
  for (grazing_intensity in c("binary","all")){
    
    if (grazing_intensity =="binary"){
      d_data_mod=save
      d_data_mod$Grazing= as.numeric(d_data_mod$Grazing!=0)
    } else{
      d_data_mod=save
    }
    
    d_data_mod=save
    
    #controlling for cover 
    model_spa_stat  = lm(value~rho_p, d_data_mod)
    
    d_data_mod$value=residuals(model_spa_stat)
    
    
    d_data_mod$Grazing=as.factor(d_data_mod$Grazing)
    d_data_mod$Grazing = relevel(d_data_mod$Grazing,ref = "0")
    
    model_spa_stat  = lmer(formula_mod, d_data_mod,
                           na.action = na.fail ,REML ="FALSE")
    
    #we remove potential outliers
    rm.outliers = romr.fnc(model_spa_stat, d_data_mod, trim=2.5)
    d_data_out = rm.outliers$data
    
    
    #saving data
    write.table(d_data_out,paste0("../Data/Linear_models_resolution/Keep_data/Data_",
                                  stat,"_",with_cover,"_aridity_no_inter_",grazing_intensity,".csv"))
    
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
    
    saveRDS(model_spa_stat,#MuMIn::get.models(select_model,subset = delta<2)[[1]],
            paste0("../Data/Linear_models_resolution/Keep_models/Mod_",stat,
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
    
  }
  #Saving each DF
  write.table(d_all,paste0("../Data/Linear_models_resolution/Importance/Importance_",
                           stat,"_",with_cover,"_aridity_no_inter_",grazing_intensity,".csv"),sep=";")
  write.table(d_all2,paste0("../Data/Linear_models_resolution/Estimator/Estimators_model_",
                            stat,"_",with_cover,"_aridity_no_inter_",grazing_intensity,".csv"),sep=";")
  
}

library(parallel)

mclapply(1:19,Run_model_importance_aridity_no_inter,mc.cores = 19)

d_all=d_all2=tibble()
for (k in list.files("../Data/Linear_models_resolution/Importance",paste0("aridity_no_inter"))){
  d=read.table(paste0("../Data/Linear_models_resolution/Importance/",k),sep=";")
  if (nrow(d)!=0){
    d_all=rbind(d_all,d)
  }
}
for (k in list.files("../Data/Linear_models_resolution/Estimator",paste0("aridity_no_inter"))){
  d2=read.table(paste0("../Data/Linear_models_resolution/Estimator/",k),sep=";")
  if (nrow(d)!=0){
    d_all2=rbind(d_all2,d2)
  }
}
write.table(d_all,paste0("../Data/Linear_models_resolution/Importance_aridity_no_inter_factor.csv"),sep=";")
write.table(d_all2,paste0("../Data/Linear_models_resolution/Estimators_model_aridity_no_inter_factor.csv"),sep=";")


## >> 2) Analyzing the residuals ---- 

d_data=read.table("../Data/Spatial_structure_grazing.csv",sep=";")%>%
  Closer_to_normality(.)
d_data[,c(1,9:29,38,40:44,47:ncol(d_data))] = 
  apply(d_data[,c(1,9:29,38,40:44,47:ncol(d_data))],2,z_tranform)
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
boot_function_lm_graz_binary = function(formula, data, indices) {
  d = data[indices,] 
  fit = lm(formula, data=d) 
  return(summary(fit)$coefficient[2,1])
}

for (with_interaction in c(F)){
  d_slope=tibble()
  
  for (grazing_intensity in c("binary","all")){
    
    for (k in c("perim_area_scaling","fmax_psd","PL_expo",
                "core_area_land","core_area","Small_patches",
                "flow_length","mean_psd","moran_I")){
      
      #for each we plot the slope against the partial residuals with and without cover
      
      d_data_out=read.table(paste0("../Data/Linear_models_resolution/Keep_data/Data_",k,"_FALSE_aridity_",
                                   ifelse(with_interaction,"","no_inter_"),grazing_intensity,".csv"),sep=" ")
      model_spa_stat=readRDS(paste0("../Data/Linear_models_resolution/Keep_models/Mod_",k,"_FALSE_aridity_",
                                    ifelse(with_interaction,"","no_inter_"),grazing_intensity,".rds"))
      
      d_data_out$Grazing = as.factor(d_data_out$Grazing)
      d_data_out$Grazing = relevel(d_data_out$Grazing,ref = "0")
      
      resid_mod=visreg::visreg(fit = model_spa_stat,xvar="Grazing",plot=F) 
      
      
      if (grazing_intensity=="binary"){
        
        mod_cov = boot(data=resid_mod$res, statistic=boot_function_lm_graz_binary, R=1000, formula=visregRes~Grazing)
        
        d_slope=rbind(d_slope,
                      tibble(pval=twoside_pvalue(mod_cov$t),
                             q2=median(mod_cov$t),
                             q1=quantile(mod_cov$t,.025),
                             q3=quantile(mod_cov$t,.975),
                             q1_90=quantile(mod_cov$t,.05),
                             q3_90=quantile(mod_cov$t,.95),
                             With_cover=F,Grazing_intensity=grazing_intensity,
                             ID_grazing=c("grazed"),
                             Stat=k,Driver="Grazing"
                      ))
      } else{
        
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
      
      
      
      resid_mod=visreg::visreg(fit = model_spa_stat,xvar="Aridity",plot=F) 
      
      
      mod_cov = boot(data=resid_mod$res, statistic=boot_function_lm_arid, R=1000, formula=visregRes~Aridity)
      
      d_slope=rbind(d_slope,
                    tibble(pval=twoside_pvalue(mod_cov$t),
                           q2=median(mod_cov$t),
                           q1=quantile(mod_cov$t,.025),
                           q3=quantile(mod_cov$t,.975),
                           q1_90=quantile(mod_cov$t,.05),
                           q3_90=quantile(mod_cov$t,.95),
                           With_cover=F,Grazing_intensity=grazing_intensity,
                           ID_grazing=c("all"),
                           Stat=k,Driver="Aridity"
                    ))
      
      
    }
  }
  
  write.table(d_slope,paste0("../Data/Linear_models_resolution/Slope_partial_residuals_aridity",
                             ifelse(with_interaction,"","_no_inter"),
                             ".csv"),sep=";")
}





# ---------------------- Step 8: Including the information on the herbivore species ----

Run_model_importance_aridity=function(id){
  
  list_mod=rbind(expand.grid(with_cover=c(T,F),
                             Stats=c("perim_area_scaling","fmax_psd","PL_expo",
                                     "core_area_land","core_area","Small_patches",
                                     "flow_length","mean_psd","moran_I")),
                 tibble(with_cover=F,
                        Stats=c("rho_p")))
  
  d_all=d_all2=tibble()#for keeping all information
  
  d_data=read.table("../Data/Spatial_structure_grazing.csv",sep=";")%>%
    Closer_to_normality(.)
  
  d_data[,c(1,9:29,38,40:44,47:ncol(d_data))] = 
    apply(d_data[,c(1,9:29,38,40:44,47:ncol(d_data))],2,z_tranform)
  
  d_data=Perform_PCA_spatial_struc(d_data)
  
  stat=list_mod$Stats[id]
  
  d_data_mod=d_data%>%melt(., measure.vars=stat)%>%
    filter(., !is.na(value))
  
  with_cover=list_mod$with_cover[id]
  
  if (with_cover){
    
    formula_mod=formula(formula_(paste("value ~ Long_cos + Long_sin + Lattitude + Slope + Elevation 
      + Aridity + Grazing + Herbivores + Herbivores*Grazing
      + Aridity*Grazing 
      + rho_p + Type_veg 
      + Sand + Org_C_v 
      + (1|Site_ID)")))
    
  } else{
    
    formula_mod=formula(formula_(paste("value ~ Long_cos + Long_sin + Lattitude + Slope + Elevation 
      + Aridity + Grazing + Herbivores + Herbivores*Grazing
      + Aridity*Grazing 
      + Type_veg 
      + Sand + Org_C_v 
      + (1|Site_ID)")))
  }
  
  save=d_data_mod
  save$Herbivores[save$Herbivores=="Donkey"]="Horse"
  
  d_data_mod=save
  
  #controlling for cover 
  model_spa_stat  = lm(value~rho_p, d_data_mod)
  
  d_data_mod$value=residuals(model_spa_stat)
  
  d_data_mod$Grazing=as.factor(d_data_mod$Grazing)
  d_data_mod$Grazing = relevel(d_data_mod$Grazing,ref = "0")
  
  model_spa_stat  = lmer(formula_mod, d_data_mod,
                         na.action = na.fail ,REML ="FALSE")
  
  #we remove potential outliers
  rm.outliers = romr.fnc(model_spa_stat, d_data_mod, trim=2.5)
  d_data_out = rm.outliers$data
  
  #saving data
  summary_coef=confint(model_spa_stat)
  
  #Merge in a df
  d_all2=tibble(Median=confint(model_spa_stat,level = 1e-9)[-1,1],
                q1=summary_coef[-1,1], 
                q3=summary_coef[-1,2],
                term=rownames(summary_coef)[-1],
                Stat=stat)%>%
    add_column(., with_cover=with_cover)
  
  write.table(d_all2,paste0("../Data/Linear_models_herbivores/Estimator/Estimators_model_",stat
                            ,"_",with_cover,"_aridity.csv"),sep=";")
  
  saveRDS(model_spa_stat,#MuMIn::get.models(select_model,subset = delta<2)[[1]],
          paste0("../Data/Linear_models_herbivores/Keep_models/Mod_",stat,"_",
                 with_cover,"_aridity.rds"))
  write.table(d_data_out,paste0("../Data/Linear_models_herbivores/Keep_data/Data_",
                                stat,"_",with_cover,"_aridity.csv"))
  
}

library(parallel)

mclapply(1:19,Run_model_importance_aridity,mc.cores = 19)


d_all2=tibble()
for (k in list.files("../Data/Linear_models_herbivores/Estimator",paste0("aridity.csv"))){
  d2=read.table(paste0("../Data/Linear_models_herbivores/Estimator/",k),sep=";")
  if (nrow(d2)!=0){
    d_all2=rbind(d_all2,d2)
  }
}
write.table(d_all2,paste0("../Data/Linear_models_herbivores/Estimators_model_aridity_factor.csv"),sep=";")






Run_model_importance_aridity_no_inter=function(id){
  
  list_mod=expand.grid(with_cover=c(T,F),
                       Stats=c("perim_area_scaling","fmax_psd","PL_expo",
                               "core_area_land","core_area","Small_patches",
                               "flow_length","mean_psd","moran_I"))
  
  
  d_all=d_all2=tibble()#for keeping all information
  
  d_data=read.table("../Data/Spatial_structure_grazing.csv",sep=";")%>%
    Closer_to_normality(.)
  
  d_data[,c(1,9:29,38,40:44,47:ncol(d_data))] = 
    apply(d_data[,c(1,9:29,38,40:44,47:ncol(d_data))],2,z_tranform)
  
  d_data=Perform_PCA_spatial_struc(d_data)
  
  stat=list_mod$Stats[id]
  
  d_data_mod=d_data%>%melt(., measure.vars=stat)%>%
    filter(., !is.na(value))
  
  with_cover=list_mod$with_cover[id]
  
  if (with_cover){
    
    formula_mod=formula(formula_(paste("value ~ Long_cos + Long_sin + Lattitude + Slope + Elevation 
      + Aridity + Grazing + Herbivores + Grazing * Herbivores
      + rho_p + Type_veg 
      + Sand + Org_C_v 
      + (1|Site_ID)")))
    
  } else{
    
    formula_mod=formula(formula_(paste("value ~ Long_cos + Long_sin + Lattitude + Slope + Elevation 
      + Aridity + Grazing + Herbivores + Grazing * Herbivores
      + Type_veg
      + Sand + Org_C_v 
      + (1|Site_ID)")))
  }
  
  save=d_data_mod%>%filter(., !is.na(Herbivores))
  save$Herbivores[save$Herbivores=="Donkey"]="Horse"
  
  
  d_data_mod=save
  
  #controlling for cover 
  model_spa_stat  = lm(value~rho_p, d_data_mod)
  
  d_data_mod$value=residuals(model_spa_stat)
  
  d_data_mod$Grazing=as.factor(d_data_mod$Grazing)
  d_data_mod$Grazing = relevel(d_data_mod$Grazing,ref = "0")
  
  model_spa_stat  = lmer(formula_mod, d_data_mod,
                         na.action = na.fail ,REML ="FALSE")
  
  #we remove potential outliers
  rm.outliers = romr.fnc(model_spa_stat, d_data_mod, trim=2.5)
  d_data_out = rm.outliers$data
  
  
  #saving data
  summary_coef=confint(model_spa_stat)
  
  #Merge in a df
  d_all2=tibble(Median=confint(model_spa_stat,level = 1e-9)[-1,1],
                q1=summary_coef[-1,1], 
                q3=summary_coef[-1,2],
                term=rownames(summary_coef)[-1],
                Stat=stat)%>%
    add_column(., with_cover=with_cover)
  
  write.table(d_all2,paste0("../Data/Linear_models_herbivores/Estimator/Estimators_model_",stat
                            ,"_",with_cover,"_aridity_no_inter.csv"),sep=";")  
  
  saveRDS(model_spa_stat,#MuMIn::get.models(select_model,subset = delta<2)[[1]],
          paste0("../Data/Linear_models_herbivores/Keep_models/Mod_",stat,"_",
                 with_cover,"_aridity_no_inter.rds"))
  write.table(d_data_out,paste0("../Data/Linear_models_herbivores/Keep_data/Data_",
                                stat,"_",with_cover,"_aridity_no_inter.csv"))
  
}

library(parallel)

mclapply(1:19,Run_model_importance_aridity_no_inter,mc.cores = 19)


d_all2=tibble()
for (k in list.files("../Data/Linear_models_herbivores/Estimator",paste0("aridity_no_inter"))){
  d2=read.table(paste0("../Data/Linear_models_herbivores/Estimator/",k),sep=";")
  if (nrow(d2)!=0){
    d_all2=rbind(d_all2,d2)
  }
}
write.table(d_all2,paste0("../Data/Linear_models_herbivores/Estimators_model_aridity_no_inter_factor.csv"),sep=";")






k="fmax_psd"
with_interaction=F
d_data_out=read.table(paste0("../Data/Linear_models_herbivores/Keep_data/Data_",k,"_FALSE_aridity",
                             ifelse(with_interaction,"","_no_inter"),".csv"),sep=" ")
model_spa_stat=readRDS(paste0("../Data/Linear_models_herbivores/Keep_models/Mod_",k,"_FALSE_aridity",
                              ifelse(with_interaction,"","_no_inter"),".rds"))

d_data_out$Grazing = as.factor(d_data_out$Grazing)
d_data_out$Grazing = relevel(d_data_out$Grazing,ref = "0")

resid_mod=visreg::visreg(fit = model_spa_stat,xvar="Grazing",plot=T) 
resid_mod=visreg::visreg(fit = model_spa_stat,xvar="Grazing",by="Herbivores",plot=T) 
