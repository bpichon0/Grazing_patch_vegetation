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
             Dung_herbivores=sapply(1:nrow(.),function(x){return(d_biodesert$`Dung_all_herbivores_g/m2`[which(d_biodesert$ID==.$Site_ID[x])])}),
  )%>%
  dplyr::mutate(., Type_veg=recode_factor(Type_veg,"1.0"="Grassland","2.0"="Shrubland","3.0"="Forest","1_2"="Grass_Shrub"))


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

#Multifunctionality
MF = read.table("../Data/Multifunctionality.csv",sep=";")

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
d$ID=sapply(1:nrow(d),function(x){
  return(d_biodesert$ID[which(d_biodesert$SITE_ID==d$Site_ID[x] & 
                                d_biodesert$GRAZ==d$Grazing[x])])
})

d$MF_CNP=sapply(1:nrow(d),function(x){
  return(MF$MF_CNP[which(MF$Full_name==paste0(d$ID[x],"_",d$Sub_id[x]))])
})
d$MF_all=sapply(1:nrow(d),function(x){
  return(MF$MF_all[which(MF$Full_name==paste0(d$ID[x],"_",d$Sub_id[x]))])
})


write.table(d,"../Data/Spatial_structure_grazing.csv",sep=";")



# ---------------------- Step 2: Effects of grazing on the spatial structure  ----
## >> 1) Running mixed-effect models ---- 

Run_model_importance_aridity=function(id){
  
  list_mod=rbind(expand.grid(with_cover=c(F),
                             Stats=c("perim_area_scaling","fmax_psd","PL_expo",
                                     "core_area_land","core_area","Small_patches",
                                     "flow_length","mean_psd","moran_I")),
                 tibble(with_cover=F,
                        Stats=c("rho_p")))
  
  
  d_data=read.table("../Data/Spatial_structure_grazing.csv",sep=";")%>%
    Closer_to_normality(.)
  
  d_data[,c(1,9:29,38,40:44,47:80,82:ncol(d_data))] = 
    apply(d_data[,c(1,9:29,38,40:44,47:80,82:ncol(d_data))],2,z_tranform)
  
  d_data=Perform_PCA_spatial_struc(d_data)
  
  stat=list_mod$Stats[id]
  
  d_data_mod=d_data%>%melt(., measure.vars=stat)%>%
    filter(., !is.na(value),!is.na(Herbivores))
  
  with_cover=list_mod$with_cover[id]
  
  if (with_cover){
    
    formula_mod=formula(formula_(paste("value ~ Long_cos + Long_sin + Lattitude + Slope + Elevation 
      + Aridity + Grazing  + Herbivores
      + Aridity*Grazing 
      + rho_p + Type_veg 
      + Sand + Org_C_v 
      + (1|Site_ID)")))
    
  } else{
    
    formula_mod=formula(formula_(paste("value ~ Long_cos + Long_sin + Lattitude + Slope + Elevation 
      + Aridity + Grazing  + Herbivores
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

    
    saveRDS(model_spa_stat,#MuMIn::get.models(select_model,subset = delta<2)[[1]],
            paste0("../Data/Linear_models_factor/Keep_models/Mod_",stat,"_",
                   with_cover,"_aridity_",grazing_intensity,".rds"))
    
    
  }  
  
  
}

library(parallel)

mclapply(1:19,Run_model_importance_aridity,mc.cores = 19)


Run_model_importance_herbivores=function(id){
  
  list_mod=rbind(expand.grid(with_cover=c(T,F),
                             Stats=c("perim_area_scaling","fmax_psd","PL_expo",
                                     "core_area_land","core_area","Small_patches",
                                     "flow_length","mean_psd","moran_I")),
                 tibble(with_cover=F,
                        Stats=c("rho_p")))
  
  
  d_all=d_all2=tibble()#for keeping all information
  
  d_data=read.table("../Data/Spatial_structure_grazing.csv",sep=";")%>%
    Closer_to_normality(.)
  
  d_data[,c(1,9:29,38,40:44,47:80,82:ncol(d_data))] = 
    apply(d_data[,c(1,9:29,38,40:44,47:80,82:ncol(d_data))],2,z_tranform)
  
  d_data=Perform_PCA_spatial_struc(d_data)
  
  stat=list_mod$Stats[id]
  
  d_data_mod=d_data%>%melt(., measure.vars=stat)%>%
    filter(., !is.na(value),!is.na(Herbivores))
  
  with_cover=list_mod$with_cover[id]
  
  if (with_cover){
    
    formula_mod=formula(formula_(paste("value ~ Long_cos + Long_sin + Lattitude + Slope + Elevation 
      + Aridity + Grazing  + Herbivores
      + rho_p + Type_veg 
      + Sand + Org_C_v 
      + (1|Site_ID)")))
    
  } else{
    
    formula_mod=formula(formula_(paste("value ~ Long_cos + Long_sin + Lattitude + Slope + Elevation 
      + Aridity + Grazing  + Herbivores
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
    
    saveRDS(model_spa_stat,#MuMIn::get.models(select_model,subset = delta<2)[[1]],
            paste0("../Data/Linear_models_factor/Keep_models/Mod_",stat,
                   "_",with_cover,"_aridity_no_inter_",grazing_intensity,".rds"))
    
    
  }
  
}

library(parallel)

mclapply(1:19,Run_model_importance_aridity_no_inter,mc.cores = 19)

## >> 2) Analyzing the residuals ---- 

d_data=read.table("../Data/Spatial_structure_grazing.csv",sep=";")%>%
  Closer_to_normality(.)
d_data[,c(1,9:29,38,40:44,47:80,82:ncol(d_data))] = 
  apply(d_data[,c(1,9:29,38,40:44,47:80,82:ncol(d_data))],2,z_tranform)
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
      
      #for each we plot the slope against the partial residuals
      
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
      
      
      
      resid_mod=visreg::visreg(fit = model_spa_stat,xvar="Herbivores",plot=F) 
      
      
      mod_cov = boot(data=resid_mod$res, statistic=boot_function_lm_graz, R=1000, formula=visregRes~Herbivores)
      
      d_slope=rbind(d_slope,
                    tibble(pval=apply(mod_cov$t,2,twoside_pvalue),
                           q2=apply((mod_cov$t),2,median),
                           q1=apply(mod_cov$t,2,quantile,.025),
                           q3=apply(mod_cov$t,2,quantile,.975),
                           q1_90=apply(mod_cov$t,2,quantile,.05),
                           q3_90=apply(mod_cov$t,2,quantile,.95),
                           With_cover=F,Grazing_intensity=grazing_intensity,
                           ID_grazing=c("Goat","Horse","Sheep"),
                           Stat=k,Driver="Herbivores"
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
                                     "flow_length","mean_psd","moran_I","Struct1","Struct2")),
                 tibble(with_cover=F,
                        Stats=c("rho_p")))
  
  d_all=d_all2=tibble()#for keeping all information
  
  d_data=read.table("../Data/Spatial_structure_grazing.csv",sep=";")%>%
    Closer_to_normality(.)%>%Perform_PCA_spatial_struc(.)
  
  d_data[,c(1,9:29,38,40:44,47:80,82:ncol(d_data))] = 
    apply(d_data[,c(1,9:29,38,40:44,47:80,82:ncol(d_data))],2,z_tranform)
  
  d_data=Perform_PCA_spatial_struc(d_data)
  
  stat=list_mod$Stats[id]
  
  d_data_mod=d_data%>%melt(., measure.vars=stat)%>%
    dplyr::filter(., !is.na(value),!is.na(Herbivores))
  
  with_cover=list_mod$with_cover[id]
  
  if (with_cover){
    
    formula_mod=formula(formula_(paste("value ~ Long_cos + Long_sin + Lattitude + Slope + Elevation 
      + Aridity + Grazing + Herbivores
      + Aridity*Grazing 
      + rho_p + Type_veg 
      + Sand + Org_C_v
      + (1|Site_ID)")))
    
  } else{
    
    formula_mod=formula(formula_(paste("value ~ Long_cos + Long_sin + Lattitude + Slope + Elevation 
      + Aridity + Grazing + Herbivores
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
    
    
    saveRDS(model_spa_stat,#MuMIn::get.models(select_model,subset = delta<2)[[1]],
            paste0("../Data/Linear_models_factor_cover_control/Keep_models/Mod_",stat,"_",
                   with_cover,"_aridity_",grazing_intensity,".rds"))
    
    
  }  
  
}

library(parallel)

mclapply(1:12,Run_model_importance_aridity,mc.cores = 12)

Run_model_importance_aridity_no_inter=function(id){
  
  list_mod=expand.grid(with_cover=c(F),
                       Stats=c("perim_area_scaling","fmax_psd","PL_expo",
                               "core_area_land","core_area","Small_patches",
                               "flow_length","mean_psd","moran_I","Struct1","Struct2"))
  
  d_all=d_all2=tibble()#for keeping all information
  
  d_data=read.table("../Data/Spatial_structure_grazing.csv",sep=";")%>%
    Closer_to_normality(.)%>%Perform_PCA_spatial_struc(.)
  
  d_data[,c(1,9:29,38,40:44,47:80,82:ncol(d_data))] = 
    apply(d_data[,c(1,9:29,38,40:44,47:80,82:ncol(d_data))],2,z_tranform)
  
  d_data=Perform_PCA_spatial_struc(d_data)
  
  stat=list_mod$Stats[id]
  
  d_data_mod=d_data%>%melt(., measure.vars=stat)%>%
    filter(., !is.na(value),!is.na(Herbivores))
  
  with_cover=list_mod$with_cover[id]
  
  if (with_cover){
    
    formula_mod=formula(formula_(paste("value ~ Long_cos + Long_sin + Lattitude + Slope + Elevation 
      + Aridity + Grazing  + Herbivores
      + rho_p + Type_veg 
      + Sand + Org_C_v 
      + (1|Site_ID)")))
    
  } else{
    
    formula_mod=formula(formula_(paste("value ~ Long_cos + Long_sin + Lattitude + Slope + Elevation 
      + Aridity + Grazing  + Herbivores
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
    
    
    saveRDS(model_spa_stat,#MuMIn::get.models(select_model,subset = delta<2)[[1]],
            paste0("../Data/Linear_models_factor_cover_control/Keep_models/Mod_",stat,
                   "_",with_cover,"_aridity_no_inter_",grazing_intensity,".rds"))
    
    
  }
  
}

library(parallel)

mclapply(1:11,Run_model_importance_aridity_no_inter,mc.cores = 11)

## >> 2) Analyzing the residuals ---- 

d_data=read.table("../Data/Spatial_structure_grazing.csv",sep=";")%>%
  Closer_to_normality(.)
d_data[,c(1,9:29,38,40:44,47:80,82:ncol(d_data))] = 
  apply(d_data[,c(1,9:29,38,40:44,47:80,82:ncol(d_data))],2,z_tranform)
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
                "flow_length","mean_psd","moran_I","Struct1","Struct2")){
      
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
      
      
      resid_mod=visreg::visreg(fit = model_spa_stat,xvar="Herbivores",plot=F) 
      
      
      mod_cov = boot(data=resid_mod$res, statistic=boot_function_lm_graz, R=1000, formula=visregRes~Herbivores)
      
      d_slope=rbind(d_slope,
                    tibble(pval=apply(mod_cov$t,2,twoside_pvalue),
                           q2=apply((mod_cov$t),2,median),
                           q1=apply(mod_cov$t,2,quantile,.025),
                           q3=apply(mod_cov$t,2,quantile,.975),
                           q1_90=apply(mod_cov$t,2,quantile,.05),
                           q3_90=apply(mod_cov$t,2,quantile,.95),
                           With_cover=F,Grazing_intensity=grazing_intensity,
                           ID_grazing=c("Goat","Horse","Sheep"),
                           Stat=k,Driver="Herbivores"
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
d_traits$Cover_quadrats=d_traits$Cover
d_traits$Pos=1:nrow(d_traits)
d_traits$Complete_name=paste0(d_traits$Genus," ", d_traits$Species)
d_traits=dplyr::arrange(d_traits,Country,Site,Plot,Complete_name)

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
    if (nrow(quadrats_k_site)>0){ #other it means that there was no recorded species
      Plot_number=as.numeric(gsub("\\D", "", colnames(quadrats_k_site)[1]))
      Site=gsub(" ","",gsub(paste0(" ",Plot_number),"",
                gsub(paste0("",Plot_number),"",
                     (colnames(quadrats_k_site)[1])))) #making sure that no plot number is present and homogenize sites names
      Site=gsub(" ","",Homogeneize_sites_name(Site))
      
      Country = Homogeneize_country_name(name_sheet_k)
      quadrats_k_site=as.data.frame(quadrats_k_site)
      quadrats_k_site[,-1]=apply(quadrats_k_site[,-1],2,as.numeric) #Making numerical
          
      #Getting species cover using the quadrats
      cover_species=tibble(Sp_name=quadrats_k_site[,1],Cover=rowSums(quadrats_k_site[,-1],na.rm = T))
      colnames(cover_species)[1]="Sp_name"
      
      
      #Getting species cover in dryfun
      species_traits_DB=d_traits[which(d_traits$Country==Country & gsub(" ","",d_traits$Site)==Site & d_traits$Plot==Plot_number),]%>%
        dplyr::arrange(.,Complete_name)
      
      test_cover=species_traits_DB$Complete_name %!in% cover_species$Sp_name  #testing which species in drypop is not in the quadrats
      test_cover2=cover_species$Sp_name %!in% species_traits_DB$Complete_name  #testing which species in the quadrats is not in drypop
      
      if (any(test_cover) & any(test_cover2)==F){ #if a species in drypop has no measured cover in the data, we keep its cover in drypop
        species_traits_DB$Cover_quadrats[-which(test_cover)] = cover_species$Cover
        if (length(species_traits_DB$Cover_quadrats[-which(test_cover)])!=length( cover_species$Cover)){
          print(paste0(Country,"_",Site,"_",Plot_number))
        }
      }else if (any(test_cover2) & any(test_cover)==F){ #if a species in the quadrats has no measured cover in drypop, we keep its cover in the quadrats
        cover_species=cover_species[-which(test_cover2),]
        species_traits_DB$Cover_quadrats=cover_species$Cover
      }else if (any(test_cover) & any(test_cover2)){
        cover_species=cover_species[-which(test_cover2),]
        species_traits_DB=species_traits_DB[-which(test_cover),]
        species_traits_DB$Cover_quadrats=cover_species$Cover
      }else{
        species_traits_DB$Cover_quadrats=cover_species$Cover
      }
      
      d_traits[which(d_traits$Pos %in% species_traits_DB$Pos),]=species_traits_DB
      
      if (!is.na(Plot_number)){
        write.table(quadrats_k_site,paste0("../Data/Traits/Cover_quadrats/",Country,"_",Site,"_",Plot_number,".csv"),sep=";")
      }
      
      list_cover_quadrat[[index]]=quadrats_k_site
      names(list_cover_quadrat)[index]= paste0(Country,"_",Site,"_",Plot_number)
      index=index+1
  }
  }
}
saveRDS(list_cover_quadrat[-grep("NA", names(list_cover_quadrat))],"../Data/Traits/All_quadrats.rds")
write.table(d_traits,"../Data/Traits/Drypop_biodesert_cover.csv",sep=";")

## >> 2) Computing CWM  ----

d_traits=read.table("../Data/Traits/Drypop_biodesert_cover.csv",sep=";")
d_CWT_FD=Add_PCA_traits(d_traits)%>%
  #mutate(., Cover=Cover_quadrats)%>%
  dplyr::filter(., !is.na(Cover))%>%#remove absent species
  melt(., measure.vars=colnames(.)[c(8:18,grep("PC",colnames(.)))])%>%#for each trait
  dplyr::group_by(.,Country,Site,Plot,variable)%>%
  dplyr::summarise(., .groups = "keep",
                   CWM = sum((Cover/sum(Cover))*value,na.rm = T), # community weighted mean
                   CW_Var = sum((Cover/sum(Cover))*((value-sum((Cover/sum(Cover))*value,na.rm = T))**2),na.rm = T), # community weighted variance
                   CW_skew = sum((Cover/sum(Cover))*((value-sum((Cover/sum(Cover))*value,na.rm = T))**3)/(sum((Cover/sum(Cover))*((value-sum((Cover/sum(Cover))*value,na.rm = T))**2),na.rm = T)**(3/2)),na.rm = T), # community weighted skewness
                   CW_kurt = sum((Cover/sum(Cover))*((value-sum((Cover/sum(Cover))*value,na.rm = T))**4)/(sum((Cover/sum(Cover))*((value-sum((Cover/sum(Cover))*value,na.rm = T))**2),na.rm = T)**2),na.rm = T), # community weighted kurtosis
                   FD = sum((Cover/sum(Cover)) * (abs(value-sum((Cover/sum(Cover))*value,na.rm = T))/
                                                    (sum(abs(value-sum((Cover/sum(Cover))*value,na.rm = T)),na.rm = T))),na.rm = T) # functional divergence
  )

d_CWM_life_form=d_traits%>%
  #mutate(., Cover=Cover_quadrats)%>%
  dplyr::filter(., !is.na(Cover))%>%#remove absent species
  dplyr::group_by(.,Country,Site,Plot)%>%
  dplyr::do(., Life_form_quadrat(.$Life_form,.$Cover))

#Diversity
d_diversity_index=d_traits%>%
  dplyr::filter(., !is.na(Cover))%>%#remove absent species
  dplyr::group_by(.,Country,Site,Plot)%>%
  dplyr::do(., Diversity_community_level(.))

d_diversity_index[,4:ncol(d_diversity_index)]=
  apply(d_diversity_index[,4:ncol(d_diversity_index)],
        2,
        function(x){
          if (any(x==1e9)){
            x[which(x==1e9)]=NA
          }
          return(x)
        })
  
d_biodesert$Complete_name=gsub(" ","",paste0(d_biodesert$COU,"_",d_biodesert$SITE,"_",d_biodesert$PLOT))
d_CWT_FD$Complete_name=sub(" ","",paste0(d_CWT_FD$Country,"_",gsub(" ","",d_CWT_FD$Site),"_",d_CWT_FD$Plot))
d_CWM_life_form$Complete_name=sub(" ","",paste0(d_CWM_life_form$Country,"_",gsub(" ","",d_CWM_life_form$Site),"_",d_CWM_life_form$Plot))
d_diversity_index$Complete_name=sub(" ","",paste0(d_diversity_index$Country,"_",gsub(" ","",d_diversity_index$Site),"_",d_diversity_index$Plot))

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

d_CWM_life_form=d_CWM_life_form%>%
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

d_diversity_index=d_diversity_index%>%
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
write.table(d_CWM_life_form,"../Data/Traits/CWM_life_forms_sites.csv",sep=";")
write.table(d_diversity_index,"../Data/Traits/Diversity_indices_whole_sites.csv",sep=";")

## >> 3) Within quadrat scale: trait dispersion, aggregation,  (for each trait) ----

dir.create("../Data/Traits/PwD",showWarnings = F)
dir.create("../Data/Traits/Diversity",showWarnings = F)
dir.create("../Data/Traits/QWT",showWarnings = F)

Run_diversity_spatialstruct_traits_quadrats=function(plot_id){
  
  print(plot_id)
  
  d_traits=read.table("../Data/Traits/Drypop_biodesert_cover.csv",sep=";")
  
  #Prior to analyses: standardize each trait 
  d_traits_scaled=Add_PCA_traits(d_traits)%>%
    #mutate(., Cover=Cover_quadrats)%>%
    dplyr::filter(., !is.na(Cover))%>%#remove absent species
    dplyr::group_by(.,Country,Site,Plot)%>% 
    dplyr::mutate(., LL=(LL-mean(LL,na.rm=T))/(max(LL,na.rm=T)-min(LL,na.rm = T)),
                  SLA=(SLA-mean(SLA,na.rm=T))/(max(SLA,na.rm=T)-min(SLA,na.rm = T)),
                  LDMC=(LDMC-mean(LDMC,na.rm=T))/(max(LDMC,na.rm=T)-min(LDMC,na.rm = T)),
                  LA=(LA-mean(LA,na.rm=T))/(max(LA,na.rm=T)-min(LA,na.rm = T)),
                  MaxH=(MaxH-mean(MaxH,na.rm=T))/(max(MaxH,na.rm=T)-min(MaxH,na.rm = T)),
                  MaxLS=(MaxLS-mean(MaxLS,na.rm=T))/(max(MaxLS,na.rm=T)-min(MaxLS,na.rm = T)),
                  d15N=(d15N-mean(d15N,na.rm=T))/(max(d15N,na.rm=T)-min(d15N,na.rm = T)),
                  d13C=(d13C-mean(d13C,na.rm=T))/(max(d13C,na.rm=T)-min(d13C,na.rm = T)),
                  # Maxvolume=(Maxvolume-mean(Maxvolume,na.rm=T))/(max(Maxvolume,na.rm=T)-min(Maxvolume,na.rm = T)),
                  Phenolics=(Phenolics-mean(Phenolics,na.rm=T))/(max(Phenolics,na.rm=T)-min(Phenolics,na.rm = T)),
                  LNC=(LNC-mean(LNC,na.rm=T))/(max(LNC,na.rm=T)-min(LNC,na.rm = T)),
                  LCC=(LCC-mean(LCC,na.rm=T))/(max(LCC,na.rm=T)-min(LCC,na.rm = T)),
                  PC1_traits=(PC1_traits-mean(PC1_traits,na.rm=T))/(max(PC1_traits,na.rm=T)-min(PC1_traits,na.rm = T)),
                  PC2_traits=(PC2_traits-mean(PC2_traits,na.rm=T))/(max(PC2_traits,na.rm=T)-min(PC2_traits,na.rm = T)),
                  PC3_traits=(PC3_traits-mean(PC3_traits,na.rm=T))/(max(PC3_traits,na.rm=T)-min(PC3_traits,na.rm = T)),
                  PC4_traits=(PC4_traits-mean(PC4_traits,na.rm=T))/(max(PC4_traits,na.rm=T)-min(PC4_traits,na.rm = T)),
                  Cover=Cover,
                  Genus=Genus,
                  Species=Species)
  
  d_traits_not_scaled=Add_PCA_traits(d_traits)%>%
    #mutate(., Cover=Cover_quadrats)%>%
    dplyr::filter(., !is.na(Cover))
  
  
  list_quadrat_sites=list.files("../Data/Traits/Cover_quadrats/")
  n_perm=100
  quadrat_k=read.table(paste0("../Data/Traits/Cover_quadrats/",list_quadrat_sites[plot_id]),sep=";")
  
  country=gsub(" ","",strsplit(list_quadrat_sites[plot_id],"_")[[1]][1])
  site=gsub(" ","",gsub('[0-9]+', '', strsplit(list_quadrat_sites[plot_id],"_")[[1]][2])) #in case 
  plot=as.numeric(gsub(".csv","",strsplit(list_quadrat_sites[plot_id],"_")[[1]][3]))
  
  mat_quadrat=as.matrix(quadrat_k[,-1]) #build the quadrat-abundance matrix
  colnames(mat_quadrat)=paste0("Quadrat_",1:100)
  rownames(mat_quadrat)=quadrat_k[,1]
  mat_quadrat=ceiling(mat_quadrat) #because some observations are integer, we round all values
  
  if (any(is.na(as.numeric(mat_quadrat))==T)){ # in case there is a NA, put a 0 instead
    mat_quadrat[is.na(mat_quadrat)]=0
  }
  if (any(colSums(mat_quadrat)==0)){ # in case there is a NA, put a 0 instead
    mat_quadrat=mat_quadrat[,-which(colSums(mat_quadrat)==0)]
  }
  
  traits_site=d_traits_scaled%>% #Get the traits
    mutate(., Site = gsub(" ","",Site))%>%
    dplyr::filter(.,Country==country,Site==site,Plot==plot)
  
  traits_site_ns=d_traits_not_scaled%>% #Get the traits
    mutate(., Site = gsub(" ","",Site))%>%
    dplyr::filter(.,Country==country,Site==site,Plot==plot)
  
  if (!is.matrix(mat_quadrat)){
    mat_quadrat=matrix(mat_quadrat,ncol = 100)
    colnames(mat_quadrat)=paste0("Quadrat_",1:100)
  }
  
  if (nrow(mat_quadrat)>1 & nrow(traits_site_ns)>1){ #At least two species and traits that have been measured
    
    #Randomization of the quadrats
    Randomization_quadrat=vegan::permatfull(mat_quadrat,times = n_perm)
    
    d_PwD=stats_QWT=d_diversity=tibble()
    for (random_id in 1:length(Randomization_quadrat$perm)){ # for each of the randomization
      
      quadrat_k_id=Randomization_quadrat$perm[[random_id]]
      
      # 1.1-Random) Spatial structure of traits within quadrats
      
      d_QWT=d_QWT_not_scaled=tibble()
      for (column_id in 1:ncol(quadrat_k_id)){# for each quadrat, compute the Quadrat weighted trait
        d_QWT=rbind(d_QWT,Get_CWT_traits(abundance = quadrat_k_id[,column_id],traits_site = traits_site)%>%
                      ungroup(.)%>%
                      add_column(., Quadrat=column_id,Randomization_id=random_id))
        d_QWT_not_scaled=rbind(d_QWT_not_scaled,Get_CWT_traits(abundance = quadrat_k_id[,column_id],traits_site = traits_site_ns)%>%
                                 ungroup(.)%>%
                                 add_column(., Quadrat=column_id,Randomization_id=random_id))
      }
      
      d_PwD=rbind(d_PwD,d_QWT%>%  #for the randomized QWT compute pairwise trait distance
                    mutate(., Randomization_id=as.character(Randomization_id))%>%
                    dplyr::group_by(., variable,Country,Site,Plot,Randomization_id)%>%
                    dplyr::do(.,Compute_trait_variance(.$CWM))%>%
                    add_column(.,PwD=d_QWT%>%  #for the randomized QWT compute pairwise trait distance
                                 mutate(., Randomization_id=as.character(Randomization_id))%>%
                                 dplyr::group_by(., variable,Country,Site,Plot,Randomization_id)%>%
                                 dplyr::do(.,Compute_pairwise_trait_distance(.$CWM))%>%
                                 dplyr::ungroup(.)%>%dplyr::select(., PwD)%>%dplyr::pull(.)))
      
      
      # 1.2-Random) mean and variance of the dominant traits in quadrats
      stats_QWT=rbind(stats_QWT,d_QWT_not_scaled%>%  #for the randomized QWT compute pairwise trait distance
                        mutate(., Randomization_id=as.character(Randomization_id))%>%
                        dplyr::group_by(., variable,Country,Site,Plot,Randomization_id)%>%
                        dplyr::summarise(., .groups = "keep",
                                         mean_QWT=mean(CWM),var_QWT=var(CWM),
                                         mean_QWT_var=mean(CWM_var)))
      
      
      # 1.3-Random) Diversity of quadrats (alpha, beta)
      d_diversity=rbind(d_diversity,Get_diversity_indices(abundance_mat = quadrat_k_id,
                                                          traits_site = traits_site_ns)%>%
                          add_column(., Randomization_id=random_id))
    }
    
    
    d_QWT_obs=d_QWT_not_scaled_obs=tibble()
    for (column_id in 1:ncol(mat_quadrat)){# for each quadrat, compute the Quadrat weighted trait mean (here on observed species abundance)
      d_QWT_obs=rbind(d_QWT_obs,Get_CWT_traits(abundance = mat_quadrat[,column_id],traits_site = traits_site)%>%
                        add_column(., Quadrat=column_id,Randomization_id="Observed"))  
      d_QWT_not_scaled_obs=rbind(d_QWT_not_scaled_obs,Get_CWT_traits(abundance = mat_quadrat[,column_id],
                                                                     traits_site = traits_site_ns)%>%
                                   add_column(., Quadrat=column_id,Randomization_id="Observed"))  
    }
    
    
    # 2.1-Observed) Spatial structure of traits within quadrats
    
    PwD_obs=d_QWT_obs%>% #for the observed QWT compute pairwise trait distance
      mutate(., Randomization_id=as.character(Randomization_id))%>%
      dplyr::group_by(., variable,Country,Site,Plot,Randomization_id)%>%
      dplyr::do(.,Compute_trait_variance(.$CWM))%>%
      add_column(.,PwD=d_QWT_obs%>% #for the observed QWT compute pairwise trait distance
                   mutate(., Randomization_id=as.character(Randomization_id))%>%
                   dplyr::group_by(., variable,Country,Site,Plot,Randomization_id)%>%
                   dplyr::do(.,Compute_pairwise_trait_distance(.$CWM))%>%
                   dplyr::ungroup(.)%>%dplyr::select(., PwD)%>%dplyr::pull(.))
    
    PwD_random=d_PwD%>%
      dplyr::ungroup(.)%>%  #then extract the mean and confidence interval across randomization 
      dplyr::group_by(., variable,Country,Site,Plot)%>%
      dplyr::summarise(., .groups = "keep",
                       mean_PwD=mean(PwD),
                       q025_PwD=quantile(PwD,.025),
                       q975_PwD=quantile(PwD,.975),
                       q05_PwD=quantile(PwD,.05),
                       q95_PwD=quantile(PwD,.95),
                       mean_Var=mean(Variance),
                       mean_sd=mean(sqrt(Variance)),
                       q025_Var=quantile(Variance,.025),
                       q975_Var=quantile(Variance,.975),
                       q05_Var=quantile(Variance,.05),
                       q95_Var=quantile(Variance,.95))
    
    PwD_tot=PwD_obs%>%
      add_column(.,
                 mean_PwD_random=PwD_random$mean_PwD,
                 q025_PwD_random=PwD_random$q025_PwD,
                 q975_PwD_random=PwD_random$q975_PwD,
                 q05_PwD_random=PwD_random$q05_PwD,
                 q95_PwD_random=PwD_random$q95_PwD,
                 mean_Var_random=PwD_random$mean_Var,
                 mean_sd_random=PwD_random$mean_sd,
                 q025_Var_random=PwD_random$q025_Var,
                 q975_Var_random=PwD_random$q975_Var,
                 q05_Var_random=PwD_random$q05_Var,
                 q95_Var_random=PwD_random$q95_Var)%>%
      dplyr::rename(., PwD_obs=PwD)%>%
      add_column(., sd_obs=sqrt(.$Variance))
    
    # 2.2-Observed) mean and variance of the dominant traits in quadrats
    
    stats_QWT_random=stats_QWT%>%
      dplyr::group_by(., variable,Country,Site,Plot)%>%
      dplyr::do(., Compute_mean_var_qwt(.$mean_QWT,.$var_QWT,.$mean_QWT_var))
    
    stats_QWT_obs=d_QWT_not_scaled_obs%>% #for the observed QWT compute pairwise trait distance
      mutate(., Randomization_id=as.character(Randomization_id))%>%
      dplyr::group_by(., variable,Country,Site,Plot,Randomization_id)%>%
      dplyr::summarise(., .groups = "keep",
                       mean_QWT=mean(CWM),
                       var_QWT=var(CWM),
                       mean_QWTvar=mean(CWM_var))
    
    
    stats_QWT_tot=stats_QWT_random%>%
      dplyr::ungroup(.)%>%  #then extract the mean and confidence interval across randomization 
      add_column(.,
                 mean_QWT_obs=stats_QWT_obs$mean_QWT,
                 var_QWT_obs=stats_QWT_obs$var_QWT,
                 mean_QWTvar_obs=stats_QWT_obs$mean_QWTvar)
    
    
    # 2.3-Observed) Diversity of quadrats (alpha, beta)
    
    d_diversity_obs=Get_diversity_indices(abundance_mat=mat_quadrat,
                                          traits_site=traits_site_ns)%>%
      add_column(., Randomization_id="Observed")
    
    Diversity_tot=d_diversity%>%
      melt(.,id.vars = "Randomization_id")%>%
      add_column(., 
                 Country=country,
                 Site=site,
                 Plot=plot)%>%
      dplyr::group_by(., variable,Country,Site,Plot)%>%
      dplyr::summarise(., .groups = "keep",
                       mean_div_random=mean(value,na.rm = T),
                       q025_div_random=quantile(value,.025,na.rm = T),
                       q975_div_random=quantile(value,.975,na.rm = T),
                       q05_div_random=quantile(value,.05,na.rm = T),
                       q95_div_random=quantile(value,.95,na.rm = T))
    
    Diversity_tot=Diversity_tot%>%
      dplyr::ungroup(.)%>%
      add_column(.,Div_observed=d_diversity_obs%>%reshape2::melt(.,id.vars="Randomization_id")%>%pull(., value))
    
    Diversity_tot[Diversity_tot==1e9]=0 #Since we use 1e9 as a numeric to avoid problems with dplyr
    
    write.table(PwD_tot,paste0("../Data/Traits/PwD/PwD_",list_quadrat_sites[plot_id]),sep=";")
    write.table(Diversity_tot,paste0("../Data/Traits/Diversity/Diversity_",list_quadrat_sites[plot_id]),sep=";")
    write.table(stats_QWT_tot,paste0("../Data/Traits/QWT/QWT_",list_quadrat_sites[plot_id]),sep=";")
  }
}
list_quadrat_sites=list.files("../Data/Traits/Cover_quadrats/")
library(parallel)
mclapply(1:length(list_quadrat_sites),Run_diversity_spatialstruct_traits_quadrats,mc.cores = 70)


for (type_d in c("PwD","Diversity","QWT")){
  
  #Aggregating data for all sites
  list_data_sites=list.files(paste0("../Data/Traits/",type_d))
  d_pwd=tibble()
  for (k in list_data_sites){
    trait_table=read.table(paste0("../Data/Traits/",type_d,"/",k),sep=";")
    if (trait_table[1,1]!="variable"){
      d_pwd=rbind(d_pwd,trait_table)
    }
  }
  d_biodesert$Complete_name=gsub(" ","",paste0(d_biodesert$COU,"_",d_biodesert$SITE,"_",d_biodesert$PLOT))
  d_pwd$Complete_name=gsub(" ","",paste0(d_pwd$Country,"_",gsub(" ","",d_pwd$Site),"_",d_pwd$Plot))
  
  d_pwd=d_pwd%>%
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
               })),
               Aridity=unlist(sapply(1:nrow(.),function(x){
                 if (length(d_biodesert$ARIDITY[which(d_biodesert$Complete_name==.$Complete_name[x])])==0){
                   return(NA)
                 }else{
                   return(d_biodesert$ARIDITY[which(d_biodesert$Complete_name==.$Complete_name[x])])
                 }
               })))
  
  write.table(d_pwd,paste0("../Data/Traits/",type_d,"_all_sites.csv"),sep=";")
}

## >> 4) Creating final datasets ----

#With the images
d_pwd=read.table("../Data/Traits/Pwd_all_sites.csv",sep=";")%>%
  add_column(., 
             Site_ID=sapply(1:nrow(.),function(x){
               return(d_biodesert$SITE_ID[which(d_biodesert$ID==.$ID[x])])
             }),
             Deviation_sd=.$PwD_obs-.$mean_PwD_random,
             Deviationgross=.$sd_obs-.$mean_sd_random,
             Deviation_var=.$Variance-.$mean_Var_random)

d_diversity=read.table("../Data/Traits/Diversity_all_sites.csv",sep=";")%>%
  add_column(., 
             Site_ID=sapply(1:nrow(.),function(x){
               return(d_biodesert$SITE_ID[which(d_biodesert$ID==.$ID[x])])
             })
  )%>%
  dplyr::filter(.,variable %!in% c("alpha_mean_star","alpha_var_star"))

d_QWT=read.table("../Data/Traits/QWT_all_sites.csv",sep=";")%>%
  add_column(., 
             Site_ID=sapply(1:nrow(.),function(x){
               return(d_biodesert$SITE_ID[which(d_biodesert$ID==.$ID[x])])
             })
  )



d_CWT_FD=read.table("../Data/Traits/CWM_FD_sites.csv",sep=";")

for (file_name in list.files("../Data/","Spatial_structure_grazing")){
  
  d_data=read.table(paste0("../Data/",file_name),sep=";")
  
  if (!any(grep("control_cover",file_name))){
    d_data=d_data%>%
      Closer_to_normality(.)
  }
  
  
  d_data=d_data%>%
    add_column(.,
               Dev_MaxLS=NA,
               CWM_MaxLS=NA
    )
  
  traits=grep("Dev",colnames(d_data))
  all_traits_name=unique(d_pwd$variable)
  
  #Adding the deviation of each trait to the global dataset 
  d_data[,traits]=as.data.frame(
    matrix(
      unlist(
        lapply(traits,function(x){ #Getting the deviation from randomness for each trait x and for each vegetation image y
          
          return(sapply(1:nrow(d_data),function(y){
            if (d_data$ID[y] %in% unique(d_pwd$ID)){
              return(d_pwd$Deviation_sd[which(d_pwd$variable==all_traits_name[grep(gsub("Dev_","",colnames(d_data)[x]),all_traits_name)] &
                                             d_pwd$ID==d_data$ID[y] & d_pwd$Graz==d_data$Grazing[y])])
            }else{
              return(NA)
            }
          })
          )
        })),
      nrow = nrow(d_data),length(traits)))
  
  traits=grep("CWM_",colnames(d_data)) 
  all_traits_name=unique(d_CWT_FD$variable)
  
  #Adding community weighted means to  to the global dataset 
  d_data[,traits]=as.data.frame(
    matrix(
      unlist(
        lapply(traits,function(x){ #Getting the deviation from randomness for each trait x and for each vegetation image y
          
          return(sapply(1:nrow(d_data),function(y){
            if (d_data$ID[y] %in% unique(d_CWT_FD$ID)){
              return(d_CWT_FD$CWM[which(d_CWT_FD$variable==all_traits_name[grep(gsub("CWM_","",colnames(d_data)[x]),all_traits_name)] &
                                          d_CWT_FD$ID==d_data$ID[y] & d_CWT_FD$Graz==d_data$Grazing[y])])
            }else{
              return(NA)
            }
          })
          )
        })),
      nrow = nrow(d_data),length(traits)))

  write.table(d_data,paste0("../Data/",gsub(".csv","",gsub("_grazing","",file_name)),"_with_traits.csv"),sep=";")
  
}


d=readxl::read_xlsx("../Data/Final_biodesert.xlsx")[,-c(112:319)]%>%
  add_column(.,Dev_MaxLS=NA,
             CWM_MaxLS=NA)%>%
  dplyr::rename(., 
                "Grazing"="GRAZ",
                "Slope"=`SLOPE-ALOS30`,
                "Elevation"=`ELE-ALOS30`,
                "Lattitude"="Lat_decimal",
                "Longitude"="Long_decimal",
                "Aridity"="ARIDITY",
                "Site_ID"="SITE_ID",
                "MAT"="AMT",
                "MAP"="RAI",
                "Herb"="RHCQ",
                "Woody"="RWCQ",
                "Grass"="RGCQ",
                "Species_richness"="SR",
                "Evenness_quadrat"="SEQ")%>%
  add_column(., 
             Long_cos=cos(.$Longitude),
             Long_sin=sin(.$Longitude)
  )%>%as.data.frame(.)

d_pwd=read.table("../Data/Traits/Pwd_all_sites.csv",sep=";")%>%
  add_column(., 
             Site_ID=sapply(1:nrow(.),function(x){
               return(d_biodesert$SITE_ID[which(d_biodesert$ID==.$ID[x])])
             }),
             Deviation_sd=.$sd_obs-.$mean_sd_random
  )


d_CWT_FD=read.table("../Data/Traits/CWM_FD_sites.csv",sep=";")

traits=grep("Dev",colnames(d))
all_traits_name=unique(d_pwd$variable)

#Adding the deviation of each trait to the global dataset 
d[,traits]=as.data.frame(
  matrix(
    unlist(
      lapply(traits,function(x){ #Getting the deviation from randomness for each trait x and for each vegetation image y
        
        return(sapply(1:nrow(d),function(y){
          if (d$ID[y] %in% unique(d_pwd$ID)){
            return(d_pwd$Deviation_sd[which(d_pwd$variable==all_traits_name[grep(gsub("Dev_","",colnames(d)[x]),all_traits_name)] &
                                              d_pwd$ID==d$ID[y] & d_pwd$Graz==d$Grazing[y])])
          }else{
            return(NA)
          }
        })
        )
      })),
    nrow = nrow(d),length(traits)))


#Community traits weighted by abundance
traits=grep("CWM_",colnames(d)) 
all_traits_name=unique(d_CWT_FD$variable)

#Adding community weighted means to  to the global dataset 
d[,traits]=as.data.frame(
  matrix(
    unlist(
      lapply(traits,function(x){ #Getting the deviation from randomness for each trait x and for each vegetation image y
        
        return(sapply(1:nrow(d),function(y){
          if (d$ID[y] %in% unique(d_CWT_FD$ID)){
            return(d_CWT_FD$CWM[which(d_CWT_FD$variable==all_traits_name[grep(gsub("CWM_","",colnames(d)[x]),all_traits_name)] &
                                        d_CWT_FD$ID==d$ID[y] & d_CWT_FD$Graz==d$Grazing[y])])
          }else{
            return(NA)
          }
        })
        )
      })),
    nrow = nrow(d),length(traits)))

#normality assumptions

d=d%>%mutate(., 
             CWM_MaxLS=log(CWM_MaxLS+1),
             Slope=log(Slope+1))

d$ID=as.character(d$ID)
d$Site_ID=as.character(d$Site_ID)
write.table(as.data.frame(d),"../Data/Biodesert_all_sites_traits.csv",sep=";")

# ---------------------- Step 4: Path analysis with traits ----

#Data of spatial structure

d_data=Perform_PCA_spatial_struc(
  Closer_to_normality_CWM(
    read.table("../Data/Spatial_structure_control_cover_with_traits.csv",sep=";")))%>%
  mutate(., Grazing=as.factor(Grazing),ID=as.character(ID))%>%
  mutate(., Dev_MaxLS=-Dev_MaxLS,#positive values = facilitation
         Grazing=scale(as.numeric(Grazing))[,1])

d_data[,colnames(dplyr::select_if(d_data, is.numeric))]=apply(d_data[,colnames(dplyr::select_if(d_data, is.numeric))],
                                                    2,
                                                    function(x){return((x-mean(x,na.rm=T))/sd(x,na.rm = T))})

d=read.table("../Data/Biodesert_all_sites_traits.csv",sep=";",quote = "")
colnames(d)=gsub("X.","",colnames(d))

current_names=c("CWM_MaxLS.","Dev_MaxLS.","Woody.","Lattitude.",
                "Long_cos.","Long_sin.","Elevation.",
                "Grazing.","Aridity.","Site_ID.")
changed_names=c("CWM_MaxLS","Dev_MaxLS","Woody","Lattitude",
                "Long_cos","Long_sin","Elevation",
                "Grazing","Aridity","Site_ID")
for (k in 1:length(current_names)){
  colnames(d)[which(colnames(d)==current_names[k])]=changed_names[k]
}
d[,colnames(dplyr::select_if(d, is.numeric))]=apply(d[,colnames(dplyr::select_if(d, is.numeric))],
                                                    2,
                                                    function(x){return((x-mean(x,na.rm=T))/sd(x,na.rm = T))})

d_mod1=dplyr::filter(d,!is.na(CWM_MaxLS))%>%
  mutate(., Grazing=as.numeric(Grazing))

mod1=lmer(CWM_MaxLS~Aridity+I(Aridity^2)+(1|Site_ID)+Grazing+Woody+Lattitude+Long_cos+Long_sin+Elevation,data=d_mod1,na.action = na.fail)

d_mod3=dplyr::filter(d,!is.na(Woody))%>%
  mutate(., Grazing=as.numeric(Grazing))

mod3=lmer(Woody~Aridity+I(Aridity^2)+Grazing+(1|Site_ID)+Lattitude+Long_cos+Long_sin+Elevation,data=d_mod3,na.action = na.fail)


#Since woody and CWM correlated, we perform PLS
#First for DEV_MAXLS

d_mod2=d%>%
  mutate(., Grazing=as.numeric(Grazing))

d_mod2=d_mod2[,c("Dev_MaxLS","CWM_MaxLS","Aridity","Woody","Grazing","Lattitude","Long_cos","Long_sin","Elevation")]%>%
  dplyr::filter(., !is.na(Dev_MaxLS),!is.na(CWM_MaxLS))%>%
  dplyr::mutate(., Dev_MaxLS= -Dev_MaxLS) #positive = facilitation

cv.modpls = cv.plsR(Dev_MaxLS~.,data=d_mod2,K=10,nt=10,modele = "pls",
                    grouplist = createFolds(d_mod2[,1], k = 10, list = F, returnTrain = FALSE),NK=200,verbose = F)

#We check number of components
res_cv=cvtable(summary(cv.modpls,MClassed=T))
#5 or 3 components according to CV

weights_mod=res_cv$CVPress[c(3,5)] #98 and 46

mod_pls = plsR(Dev_MaxLS~.,data=d_mod2,nt=10,pvals.expli=TRUE,typeVC = "adaptative")

#2 or 3 components according to AIC 


#We perform the two PLS
mod_pls = plsR(Dev_MaxLS~.,data=d_mod2,nt=3,pvals.expli=TRUE,typeVC = "adaptative")
boot_pls_mod_3= bootpls(mod_pls,typeboot="fmodel_np",R=2000,verbose = T)

mod_pls = plsR(Dev_MaxLS~.,data=d_mod2,nt=5,pvals.expli=TRUE,typeVC = "adaptative")
boot_pls_mod_5= bootpls(mod_pls,typeboot="fmodel_np",R=2000,verbose = T)
#R2 = 0.39

#And coefficient by the relative support of the two PLS
d_result_mod2 = data.frame(predictor=c("CWM_MaxLS","Aridity","Woody","Grazing","Lattitude","Long_cos","Long_sin","Elevation"),
                          q2=cbind(boot_pls_mod_3$t0[-1,],boot_pls_mod_5$t0[-1,]) %*% prop.table(weights_mod), #weighted mean
                          mean=cbind(apply(boot_pls_mod_3$t[,-1],2,mean),
                                     apply(boot_pls_mod_5$t[,-1],2,mean))%*% prop.table(weights_mod),
                          q1=apply(cbind(apply(boot_pls_mod_3$t[,-1],2,quantile,.025),
                                         apply(boot_pls_mod_5$t[,-1],2,quantile,.025)),1,
                                   function(x){return(min((x)))}),#return largest incertitude
                          q3=apply(cbind(apply(boot_pls_mod_3$t[,-1],2,quantile,.975),
                                         apply(boot_pls_mod_5$t[,-1],2,quantile,.975)),1,
                                   function(x){return(max((x)))}),#return largest incertitude
                          sd=apply(cbind(apply(boot_pls_mod_3$t[,-1],2,sd),
                                         apply(boot_pls_mod_5$t[,-1],2,sd)),1,
                                   function(x){return(max((x)))})#return largest incertitude
)


d_mod4=d_data[,c("Struct1","Dev_MaxLS","CWM_MaxLS","Woody","Lattitude","Long_cos","Long_sin","Elevation")]

cv.modpls = cv.plsR(Struct1~.,data=d_mod4,K=10,nt=10,modele = "pls",
                    grouplist = createFolds(d_mod4[,1], k = 10, list = F, returnTrain = FALSE),NK=200,verbose = F)

#We check number of components
res_cv=cvtable(summary(cv.modpls,MClassed=T))
weights_mod=res_cv$CVPress[c(4,5,6)]
#4 5 6 for CV. 43,65,82

mod_pls = plsR(Struct1~.,data=d_mod4,nt=10,pvals.expli=TRUE,typeVC = "adaptative")
mod_pls$AIC
#R2 =.38
#3 or 4 for AIC components according to AIC and cross validation


#We perform the three PLS
mod_pls = plsR(Struct1~.,data=d_mod4,nt=4,pvals.expli=TRUE,typeVC = "adaptative")
boot_pls_mod_4= bootpls(mod_pls,typeboot="fmodel_np",R=2000,verbose = T)

mod_pls = plsR(Struct1~.,data=d_mod4,nt=5,pvals.expli=TRUE,typeVC = "adaptative")
boot_pls_mod_5= bootpls(mod_pls,typeboot="fmodel_np",R=2000,verbose = T)

mod_pls = plsR(Struct1~.,data=d_mod4,nt=6,pvals.expli=TRUE,typeVC = "adaptative")
boot_pls_mod_6= bootpls(mod_pls,typeboot="fmodel_np",R=2000,verbose = T)


#And coefficient by the relative support of the two PLS
d_result_mod4 = data.frame(predictor=c("Dev_MaxLS","CWM_MaxLS","Woody","Lattitude","Long_cos","Long_sin","Elevation"),
                          q2=cbind(boot_pls_mod_4$t0[-1,],boot_pls_mod_5$t0[-1,],boot_pls_mod_6$t0[-1,]) %*% prop.table(weights_mod), #weighted mean
                          mean=cbind(apply(boot_pls_mod_4$t[,-1],2,mean),
                                     apply(boot_pls_mod_5$t[,-1],2,mean),
                                     apply(boot_pls_mod_6$t[,-1],2,mean)) %*% prop.table(weights_mod),
                          q1=apply(cbind(apply(boot_pls_mod_4$t[,-1],2,quantile,.025),
                                         apply(boot_pls_mod_5$t[,-1],2,quantile,.025),
                                         apply(boot_pls_mod_6$t[,-1],2,quantile,.025)),1,
                                   function(x){return(min((x)))}),#return largest incertitude
                          q3=apply(cbind(apply(boot_pls_mod_4$t[,-1],2,quantile,.975),
                                         apply(boot_pls_mod_5$t[,-1],2,quantile,.975),
                                         apply(boot_pls_mod_6$t[,-1],2,quantile,.975)),1,
                                   function(x){return(max((x)))}),#return largest incertitude
                          sd=apply(cbind(apply(boot_pls_mod_4$t[,-1],2,sd),
                                         apply(boot_pls_mod_5$t[,-1],2,sd),
                                         apply(boot_pls_mod_6$t[,-1],2,sd)),1,
                                   function(x){return(max((x)))})#return largest incertitude
)


r.squaredGLMM(mod1)
r.squaredGLMM(mod3)


boot_mod1 = bootstrap(mod1, .f = fixef, type = "parametric", B = 1000)
boot_mod3 = bootstrap(mod3, .f = fixef, type = "parametric", B = 1000)

d_all=tibble(
  Response=c(rep("CWM_MaxLS",dim(boot_mod1$replicates)[2]),
             rep("Dev_MaxLS",nrow(d_result_mod2)),
             rep("Woody",dim(boot_mod3$replicates)[2]),
             rep("Struct_PC1",nrow(d_result_mod4))),
  Pval=c(apply(boot_mod1$replicates,2,twoside_pvalue),
         apply(boot_pls_mod_3$t[,-1],2,twoside_pvalue),#apply(boot_pls_mod_8$t[,-1],2,twoside_pvalue) give the same
         apply(boot_mod3$replicates,2,twoside_pvalue),
         apply(boot_pls_mod_4$t[,-1],2,twoside_pvalue)),
  Variable=c(colnames(boot_mod1$replicates),
             d_result_mod2$predictor,
             colnames(boot_mod3$replicates),
             d_result_mod4$predictor),
  q2=c(apply(boot_mod1$replicates,2,median),
       d_result_mod2$q2,
       apply(boot_mod3$replicates,2,median),
       d_result_mod4$q2),
  mean=c(apply(boot_mod1$replicates,2,mean),
         d_result_mod2$mean,
         apply(boot_mod3$replicates,2,mean),
         d_result_mod4$mean),
  sd=c(apply(boot_mod1$replicates,2,sd),
       d_result_mod2$sd,
       apply(boot_mod3$replicates,2,sd),
       d_result_mod4$sd),
  q1=c(apply(boot_mod1$replicates,2,quantile,.025),
       d_result_mod2$q1,
       apply(boot_mod3$replicates,2,quantile,.025),
       d_result_mod4$q1),
  q3=c(apply(boot_mod1$replicates,2,quantile,.975),
       d_result_mod2$q3,
       apply(boot_mod3$replicates,2,quantile,.975),
       d_result_mod4$q3)
)

dir.create("../Data/Dsep",showWarnings = F)
write.table(d_all,"../Data/Dsep/CI_links_dsep.csv",sep=";")

# ---------------------- Step 5: Using the Schneider model with different spatial structures ----
## >> 1) Making simulations ----
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




## >> 2) Effect size parameter ----

parameters=matrix(c(0.001,.3,0.05,.14,0.1,0,.001,.319,0.05,.202,0.17,0,0.001,.325,0.05,.20802,0.18,0),9,2,byrow = T)
d=read.table("../Data/Simulations/Minimal_examples_stats.csv",sep=";")%>%
  add_column(., Sim_ID=rep(1:3,each=3))%>%
  add_column(., ID=rep(c("High","Medium","Low"),3),
             Cover=rep(c("High","Medium","Low"),each=3),
             m=parameters[,1],
             g0=parameters[,2])

d_effect=tibble()
# par(mfrow=c(3,3))
for (stat in c("PL_expo","fmax_psd","flow_length",
               "moran_I","core_area","Small_patches",
               "mean_psd")){
  for (cover_id in c("Low","Medium","High")){
    d2=d%>%melt(., measure.vars=stat)%>%filter(., Cover==cover_id)
    # plot(d2$g0,d2$value,ylab=stat)
    model_stat=lm(scale(value)~scale(g0),data = d2)
    
    d_effect=rbind(d_effect,tibble(q2=confint(model_stat,level=0)[2,1],
                                   Stat=stat,Cover=cover_id))
  }
}

write.table(d_effect,"../Data/Simulations/Model_coefficients.csv",sep=";")




d=read.table("./All_sim.csv",sep=";")
d_effect=tibble()
par(mfrow=c(3,3))
for (stat in c("PL_expo","fmax_psd","flow_length",
               "moran_I","core_area","Small_patches",
               "mean_psd")){
    d2=d%>%melt(., measure.vars=stat)%>%filter(.,!is.na(value))
    if (stat!="rho_p"){
      d2=d2[which(d2$rho_p>.05 & d2$rho_p<.7),]
    }
    d2$value=residuals(lm(value~rho_p,data = d2))
    model_stat=lm(scale(value)~scale(g0),data = d2)
    
    d_effect=rbind(d_effect,tibble(q2=confint(model_stat,level=0)[2,1],
                                   q1=confint(model_stat,level=.95)[2,1],
                                   q3=confint(model_stat,level=.95)[2,2],
                                   Stat=stat))
    
    plot(d2$rho_p,d2$value,ylab=stat)
    boxplot(d2$value~d2$g0,ylab=stat)
    # print(ggplot(d2)+geom_point(aes(x=g0,value))+geom_line(aes(x=g0,value,group=ID))+
    #   geom_smooth(aes(x=g0,value)))
    
    # plot(d2$value,d2$rho_p)

}
ggplot(d2)+geom_boxplot(aes(x=g0,value,group=g0))
write.table(d_effect,"../Data/Model_coefficients.csv",sep=";")




for (stat in c("PL_expo","fmax_psd","flow_length",
               "moran_I","core_area","Small_patches",
               "mean_psd")){
  
  for (unique_id in unique(d$ID)){
    d2=d%>%melt(., measure.vars=stat)%>%filter(.,!is.na(value),ID==unique_id)
    if (stat!="rho_p"){
      d2=d2[which(d2$rho_p>.05 & d2$rho_p<.7),]
    }
    # d2$value=residuals(lm(value~rho_p,data = d2))
    model_stat=lm(scale(value)~scale(g0)+scale(rho_p),data = d2)
    
    d_effect=rbind(d_effect,tibble(q2=confint(model_stat,level=0)[2,1],
                                   q1=confint(model_stat,level=.95)[2,1],
                                   q3=confint(model_stat,level=.95)[2,2],
                                   Stat=stat,
                                   ID=unique_id))
    
    # plot(d2$rho_p,d2$value,ylab=stat)
    # boxplot(d2$value~d2$g0,ylab=stat)
    # print(ggplot(d2)+geom_point(aes(x=g0,value))+geom_line(aes(x=g0,value,group=ID))+
    #   geom_smooth(aes(x=g0,value)))
    
    # plot(d2$value,d2$rho_p)
  }
  
}

## >> 3) 2D plot ----

# Run_PA_schneider=function(id){
#   d=tibble()
#   k=1
#   index=1
#   g0_=seq(0,.5,length.out=150)[id]
#   
#   for (b_ in rev(seq(0,1,length.out=150))){
#     if (index>1){
#       if (d$Cover[nrow(d)]==0 & d$Cover[nrow(d)-1]==0){#both desert
#         d=rbind(d,tibble(g0=g0_,
#                          b=b_,
#                          Branch="Restoration",
#                          Cover=0))
#         d=rbind(d,tibble(g0=param_list$g0[k],
#                          b=param_list$b[k],
#                          Branch="Degradation",
#                          Cover=0))
#         
#       }else{
#         test=Compute_ode(param = Get_params_model(r=.01,d=.1,b=b_,g0=g0_),method_ode = "ode45",ini_cover = c(.02,.49,.49),time_max = 700)
#         test2=Compute_ode(param = Get_params_model(r=.01,d=.1,b=b_,g0=g0_),method_ode = "ode45",ini_cover = c(.8,.1,.1),time_max = 700)
#         d=rbind(d,tibble(g0=g0_,
#                          b=b_,
#                          Branch="Restoration",
#                          Cover=round(test$Rho_p[nrow(test)-10],3)))
#         d=rbind(d,tibble(g0=g0_,
#                          b=b_,
#                          Branch="Degradation",
#                          Cover=round(test2$Rho_p[nrow(test)-10],3)))
#         
#       }
#     }else{
#       test=Compute_ode(param = Get_params_model(r=.01,d=.1,b=b_,g0=g0_),method_ode = "ode45",ini_cover = c(.02,.49,.49),time_max = 700)
#       test2=Compute_ode(param = Get_params_model(r=.01,d=.1,b=b_,g0=g0_),method_ode = "ode45",ini_cover = c(.8,.1,.1),time_max = 700)
#       d=rbind(d,tibble(g0=g0_,
#                        b=b_,
#                        Branch="Restoration",
#                        Cover=round(test$Rho_p[nrow(test)-10],3)))
#       d=rbind(d,tibble(g0=g0_,
#                        b=b_,
#                        Branch="Degradation",
#                        Cover=round(test2$Rho_p[nrow(test)-10],3)))
#       
#     }
#     
#     index=index+1
#     print(k)
#     k=k+1
#     
#   }
#   
#   write.table(d,paste0("./TEST/PA_",id,".csv"),sep=";")
#   
# }
# library(parallel)
# mclapply(1:150,Run_PA_schneider,mc.cores=75)
# d=tibble()
# for (k in 1:150){
#   d=rbind(d,read.table(paste0("./TEST/PA_",k,".csv"),sep=";"))
# }
# write.table(d,"./All_d.csv",sep=";")


Run_PA_schneider=function(id){
  d=tibble()
  k=1
  index=1
  g0_=seq(0,.3,length.out=50)[id]
  m_=rev(seq(0.001,.18,length.out=50))[id]
  
  for (b_ in rev(seq(0,1,length.out=50))){
    if (index>1){
      if (d$Cover[nrow(d)]==0 & d$Cover[nrow(d)-1]==0){#both desert
        d=rbind(d,tibble(g0=g0_,
                         b=b_,
                         Branch="Restoration",
                         Cover=0))
        d=rbind(d,tibble(g0=param_list$g0[k],
                         b=param_list$b[k],
                         Branch="Degradation",
                         Cover=0))
        
      }else{
        test=Compute_ode(param = Get_params_model(r=.01,d=.1,b=b_,g0=g0_,m=m_),method_ode = "ode45",ini_cover = c(.02,.49,.49),time_max = 700)
        test2=Compute_ode(param = Get_params_model(r=.01,d=.1,b=b_,g0=g0_,m=m_),method_ode = "ode45",ini_cover = c(.8,.1,.1),time_max = 700)
        d=rbind(d,tibble(g0=g0_,
                         b=b_,
                         Branch="Restoration",
                         Cover=round(test$Rho_p[nrow(test)-10],3)))
        d=rbind(d,tibble(g0=g0_,
                         b=b_,
                         Branch="Degradation",
                         Cover=round(test2$Rho_p[nrow(test)-10],3)))
        
      }
    }else{
      test=Compute_ode(param = Get_params_model(r=.01,d=.1,b=b_,g0=g0_,m=m_),method_ode = "ode45",ini_cover = c(.02,.49,.49),time_max = 700)
      test2=Compute_ode(param = Get_params_model(r=.01,d=.1,b=b_,g0=g0_,m=m_),method_ode = "ode45",ini_cover = c(.8,.1,.1),time_max = 700)
      d=rbind(d,tibble(g0=g0_,
                       b=b_,
                       Branch="Restoration",
                       Cover=round(test$Rho_p[nrow(test)-10],3)))
      d=rbind(d,tibble(g0=g0_,
                       b=b_,
                       Branch="Degradation",
                       Cover=round(test2$Rho_p[nrow(test)-10],3)))
      
    }
    
    index=index+1
    print(k)
    k=k+1
    
  }
  
  write.table(d,paste0("../Data/Simulations/PA_",id,".csv"),sep=";")
  
}
library(parallel)
mclapply(1:50,Run_PA_schneider,mc.cores=50)


d=tibble()
for (k in 1:50){
  d=rbind(d,read.table(paste0("../Data/Simulations/PA_",k,".csv"),sep=";"))
}
write.table(d,"../Data/Simulations/All_simulations_PA.csv",sep=";")

# ---------------------- Step 6: Sensitivity spatial resolution ----
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
  
  d_data[,c(1,9:29,38,40:44,47:80,82:ncol(d_data))] = 
    apply(d_data[,c(1,9:29,38,40:44,47:80,82:ncol(d_data))],2,z_tranform)
  
  d_data=Perform_PCA_spatial_struc(d_data)
  
  stat=list_mod$Stats[id]
  
  d_data_mod=d_data%>%melt(., measure.vars=stat)%>%
    filter(., !is.na(value),!is.na(Herbivores))
  
  with_cover=list_mod$with_cover[id]
  
  if (with_cover){
    
    formula_mod=formula(formula_(paste("value ~ Long_cos + Long_sin + Lattitude + Slope + Elevation 
      + Aridity + Grazing + Resolution + Herbivores
      + Aridity*Grazing 
      + rho_p + Type_veg 
      + Sand + Org_C_v 
      + (1|Site_ID)")))
    
  } else{
    
    formula_mod=formula(formula_(paste("value ~ Long_cos + Long_sin + Lattitude + Slope + Elevation 
      + Aridity + Grazing + Resolution + Herbivores
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

mclapply(1:10,Run_model_importance_aridity,mc.cores = 10)


Run_model_importance_aridity_no_inter=function(id){
  
  list_mod=rbind(expand.grid(with_cover=c(F),
                             Stats=c("perim_area_scaling","fmax_psd","PL_expo",
                                     "core_area_land","core_area","Small_patches",
                                     "flow_length","mean_psd","moran_I")),
                 tibble(with_cover=F,
                        Stats=c("rho_p")))
  
  
  d_all=d_all2=tibble()#for keeping all information
  
  d_data=read.table("../Data/Spatial_structure_grazing.csv",sep=";")%>%
    Closer_to_normality(.)
  
  d_data[,c(1,9:29,38,40:44,47:80,82:ncol(d_data))] = 
    apply(d_data[,c(1,9:29,38,40:44,47:80,82:ncol(d_data))],2,z_tranform)
  
  d_data=Perform_PCA_spatial_struc(d_data)
  
  stat=list_mod$Stats[id]
  
  d_data_mod=d_data%>%melt(., measure.vars=stat)%>%
    filter(., !is.na(value),!is.na(Herbivores))
  
  with_cover=list_mod$with_cover[id]
  
  if (with_cover){
    
    formula_mod=formula(formula_(paste("value ~ Long_cos + Long_sin + Lattitude + Slope + Elevation 
      + Aridity + Grazing + Resolution + Herbivores 
      + rho_p + Type_veg 
      + Sand + Org_C_v 
      + (1|Site_ID)")))
    
  } else{
    
    formula_mod=formula(formula_(paste("value ~ Long_cos + Long_sin + Lattitude + Slope + Elevation 
      + Aridity + Grazing + Resolution + Herbivores
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

mclapply(1:10,Run_model_importance_aridity_no_inter,mc.cores = 10)

## >> 2) Analyzing the residuals ---- 

d_data=read.table("../Data/Spatial_structure_grazing.csv",sep=";")%>%
  Closer_to_normality(.)
d_data[,c(1,9:29,38,40:44,47:80,82:ncol(d_data))] = 
  apply(d_data[,c(1,9:29,38,40:44,47:80,82:ncol(d_data))],2,z_tranform)
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







