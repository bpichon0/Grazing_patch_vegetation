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
mclapply(list_f[grep(list_f,pattern = "biodesert")],Compute_all_metrics_biodesert,mc.cores = 5)



dir.create("../Data/Metric_null",showWarnings = F)

Compute_all_metrics_biodesert_null=function(k){
  
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
  
  
  all_dist = (lapply(seq.int(199), function(i) {
    null_mat = matrix(sample(landscape), nrow = nrow(landscape), ncol = ncol(landscape))
    psd=spatialwarnings::patchdistr_sews(null_mat>0)
    max_patchsize=max(psd$psd_obs)
    cv_patch=sd(psd$psd_obs)/mean(psd$psd_obs)
    mean_nb_neigh = mean(simecol::neighbors(x =null_mat,state = 1, wdist =  matrix( c(0, 1, 0,1, 0, 1, 0, 1, 0), nrow = 3),bounds = 1)[which(null_mat == 1)])
    var_spa=spatialwarnings::raw_cg_variance(null_mat>0)
    
    return(data.frame(max_patchsize=log(max_patchsize),
                  cv_patch=cv_patch,
                  mean_nb_neigh=mean_nb_neigh,
                  var_spa=var_spa))
    
  }))
  
  d2=matrix(unlist(all_dist),199,4,byrow = T)
  
  #observed
  psd=spatialwarnings::patchdistr_sews(landscape>0)
  max_patchsize=max(psd$psd_obs)
  cv_patch=sd(psd$psd_obs)/mean(psd$psd_obs)
  mean_nb_neigh = mean(simecol::neighbors(x =landscape,state = 1, wdist =  matrix( c(0, 1, 0,1, 0, 1, 0, 1, 0), nrow = 3),bounds = 1)[which(landscape == 1)])
  var_spa=spatialwarnings::raw_cg_variance(landscape>0)
  
  d_obs=data.frame(max_patchsize=log(max_patchsize),
         cv_patch=cv_patch,
         mean_nb_neigh=mean_nb_neigh,
         var_spa=var_spa)
  
  for (k in 1:4){d2[,k]=abs(d2[,k]-d_obs[,k])}
  
  d=as_tibble(t(colMeans(abs(d2))))
  colnames(d)=colnames(all_dist[[1]])
  
  d=d%>% #Computing all summary statistics
    add_column(.,Full_name=name,Resolution=Get_spatial_resolution(landscape), #full name & spatial resolution
               Site_ID=id_site, # ID of the site
               Sub_id=sub_id) # id of the sublandscape taken at a given site
  
  write.table(d,paste0("../Data/Metrics_null/Metric_",name,".csv"),sep=";")  
}

info_kmean=read.table("../Data/Landscapes/Kmean_clust_info.csv",sep=",",header = T)%>%
  filter(., Size!=200,Dataset=="biodesert",status=="kept") #keeping the kept sites


library(parallel)
mclapply(paste0(info_kmean$Site_ID,"_",info_kmean$Image,".csv"),Compute_all_metrics_biodesert_null,mc.cores = 70)


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

write.table(d,"../Data/Spatial_structure_grazing.csv",sep=";")

# ---------------------- Step 2: Running mixed-effect models  ----



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
  
ggsave("../Figures/Linear_models/Correlation_btw_predictors.pdf",p,width = 5,height = 4)





Run_model_importance_aridity=function(id){
  
  list_mod=rbind(expand.grid(with_cover=c(T,F),
                             Stats=c("perim_area_scaling","fmax_psd","PL_expo","Spectral_ratio",
                                     "core_area_land","division","fractal_dim","contig","core_area",
                                     "flow_length","PLR","KS_dist",
                                     "Struct1","Struct2")), #We also add the case where the stats correspond to the first two components of the PCA
                 tibble(with_cover=F,Stats="rho_p"))
  
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
  
  for(grazing_intensity in c("low","high","all")){
    
    if (grazing_intensity =="low"){
      d_data_mod=save%>%filter(., Grazing %in% c(0,1))
    } else if (grazing_intensity =="high"){
      d_data_mod=save%>%filter(., Grazing %in% c(2,3))
    } else{
      d_data_mod=save
    }
    
    d_data_mod$Grazing=scale(d_data_mod$Grazing)[,1]
    
    model_spa_stat  = lmer(formula_mod, d_data_mod,
                           na.action = na.fail ,REML ="FALSE")
    
    #we remove potential outliers
    rm.outliers = romr.fnc(model_spa_stat, d_data_mod, trim=2.5)
    d_data_out = rm.outliers$data
    
    #saving data
    write.table(d_data_out,paste0("../Data/Linear_models/Keep_data/Data_",stat,"_",with_cover,"_aridity_",grazing_intensity,".csv"))
    
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
    
    saveRDS(model_spa_stat,paste0("../Data/Linear_models/Keep_models/Mod_",stat,"_",with_cover,"_aridity_",grazing_intensity,".rds"))
    
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
  write.table(d_all,paste0("../Data/Linear_models/Importance/Importance_",stat,"_",with_cover,"_aridity.csv"),sep=";")
  write.table(d_all2,paste0("../Data/Linear_models/Estimator/Estimators_model_",stat,"_",with_cover,"_aridity.csv"),sep=";")
  
}

library(parallel)

mclapply(1:30,Run_model_importance_aridity,mc.cores = 1)

d_all=d_all2=tibble()
for (k in list.files("../Data/Linear_models/Importance",paste0("aridity.csv"))){
  d=read.table(paste0("../Data/Linear_models/Importance/",k),sep=";")
  if (nrow(d)!=0){
    d_all=rbind(d_all,d)
  }
}
for (k in list.files("../Data/Linear_models/Estimator",paste0("aridity.csv"))){
  d2=read.table(paste0("../Data/Linear_models/Estimator/",k),sep=";")
  if (nrow(d2)!=0){
    d_all2=rbind(d_all2,d2)
  }
}
write.table(d_all,paste0("../Data/Linear_models/Importance_","aridity.csv"),sep=";")
write.table(d_all2,paste0("../Data/Linear_models/Estimators_model_","aridity.csv"),sep=";")


Run_model_importance_aridity_no_inter=function(id){
  
  list_mod=rbind(expand.grid(with_cover=c(T,F),
                       Stats=c("perim_area_scaling","fmax_psd","PL_expo","Spectral_ratio",
                               "core_area_land","division","fractal_dim","contig","core_area",
                               "flow_length","PLR","KS_dist",
                               "Struct1","Struct2")), #We also add the case where the stats correspond to the first two components of the PCA
                tibble(with_cover=F,Stats="rho_p"))
  
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
  
  for(grazing_intensity in c("low","high","all")){
    
    if (grazing_intensity =="low"){
      d_data_mod=save%>%filter(., Grazing %in% c(0,1))
    } else if (grazing_intensity =="high"){
      d_data_mod=save%>%filter(., Grazing %in% c(2,3))
    } else{
      d_data_mod=save
    }
    
    d_data_mod$Grazing=scale(d_data_mod$Grazing)[,1]
    
    model_spa_stat  = lmer(formula_mod, d_data_mod,
                           na.action = na.fail ,REML ="FALSE")
    
    #we remove potential outliers
    rm.outliers = romr.fnc(model_spa_stat, d_data_mod, trim=2.5)
    d_data_out = rm.outliers$data
    
    #saving data
    write.table(d_data_out,paste0("../Data/Linear_models/Keep_data/Data_",stat,"_",with_cover,"_aridity_no_inter_",grazing_intensity,".csv"))
    
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
    
    if (dim(select_model%>%filter(., delta<2))[1]==1){ #one model
      
      result_select=model_spa_stat
      R2=select_model%>%filter(., AICc<min(AICc)+2)
      
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
    
    saveRDS(model_spa_stat,paste0("../Data/Linear_models/Keep_models/Mod_",stat,"_",with_cover,"_aridity_no_inter_",grazing_intensity,".rds"))
    
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
  write.table(d_all,paste0("../Data/Linear_models/Importance/Importance_",stat,"_",with_cover,"_aridity_no_inter.csv"),sep=";")
  write.table(d_all2,paste0("../Data/Linear_models/Estimator/Estimators_model_",stat,"_",with_cover,"_aridity_no_inter.csv"),sep=";")
  
}

library(parallel)

mclapply(1:30,Run_model_importance_aridity_no_inter,mc.cores = 1)

d_all=d_all2=tibble()
for (k in list.files("../Data/Linear_models/Importance",paste0("aridity_no_inter"))){
  d=read.table(paste0("../Data/Linear_models/Importance/",k),sep=";")
  if (nrow(d)!=0){
    d_all=rbind(d_all,d)
  }
}
for (k in list.files("../Data/Linear_models/Estimator",paste0("aridity_no_inter"))){
  d2=read.table(paste0("../Data/Linear_models/Estimator/",k),sep=";")
  if (nrow(d)!=0){
    d_all2=rbind(d_all2,d2)
  }
}
write.table(d_all,paste0("../Data/Linear_models/Importance_","aridity_no_inter.csv"),sep=";")
write.table(d_all2,paste0("../Data/Linear_models/Estimators_model_","aridity_no_inter.csv"),sep=";")



# ---------------------- Step 3: Analyzing the residuals  ----
## >> 1) Model with all climatic variables ----

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
  
  
  d_data_out=read.table(paste0("../Data/Linear_models/Keep_data/Data_",k,"_TRUE.csv"),sep=" ")
  model_spa_stat=readRDS(paste0("../Data/Linear_models/Keep_models/Mod_",k,"_TRUE.rds"))
  
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
  
  d_data_out=read.table(paste0("../Data/Linear_models/Keep_data/Data_",k,"_FALSE.csv"),sep=" ")
  model_spa_stat=readRDS(paste0("../Data/Linear_models/Keep_models/Mod_",k,"_FALSE.rds"))
  
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

write.table(d_slope,"../Data/Linear_models/Slope_partial_residuals.csv",sep=";")

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
  
  
  d_data_out=read.table(paste0("../Data/Linear_models/Keep_data/Data_",k,"_TRUE.csv"),sep=" ")
  model_spa_stat=readRDS(paste0("../Data/Linear_models/Keep_models/Mod_",k,"_TRUE.rds"))
  
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
  
  d_data_out=read.table(paste0("../Data/Linear_models/Keep_data/Data_",k,"_FALSE.csv"),sep=" ")
  model_spa_stat=readRDS(paste0("../Data/Linear_models/Keep_models/Mod_",k,"_FALSE.rds"))
  
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

write.table(d_slope,"../Data/Linear_models/Slope_partial_residuals.csv",sep=";")



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
    
    d_data_out=  read.table(paste0("../Data/Linear_models/Keep_data/Data_",stat,"_FALSE.csv"),sep=";")
    if (ncol(d_data_out)==1){
      d_data_out=  read.table(paste0("../Data/Linear_models/Keep_data/Data_",stat,"_FALSE.csv"),sep=" ")
    }
    model_spa_stat=readRDS(paste0("../Data/Linear_models/Keep_models/Mod_",stat,"_FALSE.rds"))
    
    
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
  
  
  ggsave(paste0("../Figures/Linear_models/",
                "Interaction_no_cover/All_interactions_",stat,".pdf"),p,width = 4*length(list_plot),height = 4)
  
  ID_title=ID_title+1
  
}








## >> 2) Model with aridity and grazing ----

d_data=read.table("../Data/Spatial_structure_grazing.csv",sep=";")%>%
  Closer_to_normality(.)
d_data[,c(1,9:29,36,38:42,45:ncol(d_data))] = 
  apply(d_data[,c(1,9:29,36,38:42,45:ncol(d_data))],2,z_tranform)
d_data=Perform_PCA_spatial_struc(d_data)
save=d_data
d_slope=tibble()

for (k in c("perim_area_scaling","fmax_psd","PL_expo","Spectral_ratio",
            "core_area_land","division","contig","core_area",
            "flow_length","PLR","KS_dist",
            "Struct1","Struct2","rho_p")){
  
  for (grazing_intensity in c("low","high","all")){
    
    #for each we plot the slope against the partial residuals with and without cover
    
    if (k !="rho_p"){
      
      d_data_out=read.table(paste0("../Data/Linear_models/Keep_data/Data_",k,"_TRUE_aridity_",grazing_intensity,".csv"),sep=" ")
      model_spa_stat=readRDS(paste0("../Data/Linear_models/Keep_models/Mod_",k,"_TRUE_aridity_",grazing_intensity,".rds"))
      
      resid_mod=visreg::visreg(fit = model_spa_stat,xvar="Grazing",plot=F) 
      
      mod_cov = boot(data=resid_mod$res, statistic=boot_function_lm, R=1000, formula=visregRes~Grazing)
      
      d_slope=rbind(d_slope,
                    tibble(pval=twoside_pvalue(mod_cov$t),
                           q2=median(mod_cov$t),
                           q1=quantile(mod_cov$t,.025),
                           q3=quantile(mod_cov$t,.975),
                           With_cover=T,Grazing_intensity=grazing_intensity,
                           Stat=k,Driver="Grazing"
                    ))
      
      resid_mod=visreg::visreg(fit = model_spa_stat,xvar="Aridity",plot=F) 
      
      mod_cov = boot(data=resid_mod$res, statistic=boot_function_lm, R=1000, formula=visregRes~Aridity)
      
      d_slope=rbind(d_slope,
                    tibble(pval=twoside_pvalue(mod_cov$t),
                           q2=median(mod_cov$t),
                           q1=quantile(mod_cov$t,.025),
                           q3=quantile(mod_cov$t,.975),
                           With_cover=T,Grazing_intensity=grazing_intensity,
                           Stat=k,Driver="Aridity"
                    ))
    }
    
    d_data_out=read.table(paste0("../Data/Linear_models/Keep_data/Data_",k,"_FALSE_aridity_",grazing_intensity,".csv"),sep=" ")
    model_spa_stat=readRDS(paste0("../Data/Linear_models/Keep_models/Mod_",k,"_FALSE_aridity_",grazing_intensity,".rds"))
    
    resid_mod=visreg::visreg(fit = model_spa_stat,xvar="Grazing",plot=F) 
    
    mod_cov = boot(data=resid_mod$res, statistic=boot_function_lm, R=1000, formula=visregRes~Grazing)
    
    d_slope=rbind(d_slope,
                  tibble(pval=twoside_pvalue(mod_cov$t),
                    q2=median(mod_cov$t),
                    q1=quantile(mod_cov$t,.025),
                    q3=quantile(mod_cov$t,.975),
                    With_cover=F,Grazing_intensity=grazing_intensity,
                    Stat=k,Driver="Grazing"
                  ))
    
    
    resid_mod=visreg::visreg(fit = model_spa_stat,xvar="Aridity",plot=F) 
    
    mod_cov = boot(data=resid_mod$res, statistic=boot_function_lm, R=1000, formula=visregRes~Aridity)
    
    d_slope=rbind(d_slope,
                  tibble(pval=twoside_pvalue(mod_cov$t),
                         q2=median(mod_cov$t),
                         q1=quantile(mod_cov$t,.025),
                         q3=quantile(mod_cov$t,.975),
                         With_cover=F,Grazing_intensity=grazing_intensity,
                         Stat=k,Driver="Aridity"
                  ))
    
    
  }
}

write.table(d_slope,"../Data/Linear_models/Slope_partial_residuals_aridity.csv",sep=";")

## >> 3) Model with aridity and grazing but without interactions ----

d_data=read.table("../Data/Spatial_structure_grazing.csv",sep=";")%>%
  Closer_to_normality(.)
d_data[,c(1,9:29,36,38:42,45:ncol(d_data))] = 
  apply(d_data[,c(1,9:29,36,38:42,45:ncol(d_data))],2,z_tranform)
d_data=Perform_PCA_spatial_struc(d_data)
save=d_data
d_slope=tibble()

for (k in c("perim_area_scaling","fmax_psd","PL_expo","Spectral_ratio",
            "core_area_land","division","contig","core_area",
            "flow_length","PLR","KS_dist",
            "Struct1","Struct2","rho_p")){
  
  for (grazing_intensity in c("low","high","all")){
    
    #for each we plot the slope against the partial residuals with and without cover
    
    if (k !="rho_p"){
      
      d_data_out=read.table(paste0("../Data/Linear_models/Keep_data/Data_",k,
                                   "_TRUE_aridity_no_inter_",grazing_intensity,".csv"),sep=" ")
      model_spa_stat=readRDS(paste0("../Data/Linear_models/Keep_models/Mod_",k,
                                    "_TRUE_aridity_no_inter_",grazing_intensity,".rds"))
      
      resid_mod=visreg::visreg(fit = model_spa_stat,xvar="Grazing",plot=F) 
      
      mod_cov = boot(data=resid_mod$res, statistic=boot_function_lm, R=1000, formula=visregRes~Grazing)
      
      d_slope=rbind(d_slope,
                    tibble(pval=twoside_pvalue(mod_cov$t),
                           q2=median(mod_cov$t),
                           q1=quantile(mod_cov$t,.025),
                           q3=quantile(mod_cov$t,.975),
                           With_cover=T,Grazing_intensity=grazing_intensity,
                           Stat=k,Driver="Grazing"
                    ))
      
      resid_mod=visreg::visreg(fit = model_spa_stat,xvar="Aridity",plot=F) 
      
      mod_cov = boot(data=resid_mod$res, statistic=boot_function_lm, R=1000, formula=visregRes~Aridity)
      
      d_slope=rbind(d_slope,
                    tibble(pval=twoside_pvalue(mod_cov$t),
                           q2=median(mod_cov$t),
                           q1=quantile(mod_cov$t,.025),
                           q3=quantile(mod_cov$t,.975),
                           With_cover=T,Grazing_intensity=grazing_intensity,
                           Stat=k,Driver="Aridity"
                    ))
    }
    d_data_out=read.table(paste0("../Data/Linear_models/Keep_data/Data_",k,
                                 "_FALSE_aridity_no_inter_",grazing_intensity,".csv"),sep=" ")
    model_spa_stat=readRDS(paste0("../Data/Linear_models/Keep_models/Mod_",k,
                                  "_FALSE_aridity_no_inter_",grazing_intensity,".rds"))
    
    resid_mod=visreg::visreg(fit = model_spa_stat,xvar="Grazing",plot=F) 
    
    mod_cov = boot(data=resid_mod$res, statistic=boot_function_lm, R=1000, formula=visregRes~Grazing)
    
    d_slope=rbind(d_slope,
                  tibble(pval=twoside_pvalue(mod_cov$t),
                         q2=median(mod_cov$t),
                         q1=quantile(mod_cov$t,.025),
                         q3=quantile(mod_cov$t,.975),
                         With_cover=F,Grazing_intensity=grazing_intensity,
                         Stat=k,Driver="Grazing"
                  ))
    
    
    resid_mod=visreg::visreg(fit = model_spa_stat,xvar="Aridity",plot=F) 
    
    mod_cov = boot(data=resid_mod$res, statistic=boot_function_lm, R=1000, formula=visregRes~Aridity)
    
    d_slope=rbind(d_slope,
                  tibble(pval=twoside_pvalue(mod_cov$t),
                         q2=median(mod_cov$t),
                         q1=quantile(mod_cov$t,.025),
                         q3=quantile(mod_cov$t,.975),
                         With_cover=F,Grazing_intensity=grazing_intensity,
                         Stat=k,Driver="Aridity"
                  ))
    
    
  }
}

write.table(d_slope,"../Data/Linear_models/Slope_partial_residuals_aridity_no_inter.csv",sep=";")


# ---------------------- Step 4: Indirect effect of grazing  ----

dir.create("../Data/SEM/",showWarnings = F)


for (k in c("fmax_psd","PL_expo","perim_area_scaling","core_area_land","core_area")){
  
  for (graz in c("all","low","high")){
    
    save=Get_data_resid_SEM(k,graz)
    
    
    d_sem=save#%>%filter(., Grazing %in% c(0,1))
    SEM_distance=((psem(
      lmer(Resid_mod ~  (1|Site_ID) + Grazing + Org_C_v +Sand+Aridity+Total_N+Total_P, d_sem),
      lmer(Total_N ~  (1|Site_ID) + Grazing +Aridity+Sand, d_sem),
      lmer(Total_P ~  (1|Site_ID) + Grazing +Aridity+Sand, d_sem),
      lmer(Org_C_v ~  (1|Site_ID) + Grazing + Total_N +Total_P + Aridity+Sand, d_sem),
      Total_N%~~%Total_P
    )))
    
    SEM_distance_boot = bootEff(SEM_distance, R = 1000, seed = 13, parallel = "snow",ran.eff = "Site_ID",ncpus = 50)
    saveRDS(SEM_distance_boot,paste0("../Data/SEM/SEM_boot_TotN_",k,"_","graz_",graz,".rds"))
    
    
    d_sem=save#%>%filter(., Grazing %in% c(0,1))
    SEM_distance=((psem(
      lmer(Resid_mod ~  (1|Site_ID) + Grazing + Org_C_v +Sand+Aridity+Nitrate+Amonium, d_sem),
      lmer(Nitrate ~  (1|Site_ID) + Grazing +Aridity+Sand, d_sem),
      lmer(Amonium ~  (1|Site_ID) + Grazing +Aridity+Sand, d_sem),
      lmer(Org_C_v ~  (1|Site_ID) + Grazing + Nitrate +Amonium + Aridity+Sand, d_sem),
      Nitrate%~~%Amonium
    )))
    
    SEM_distance_boot = bootEff(SEM_distance, R = 1000, seed = 13, parallel = "snow",ran.eff = "Site_ID",ncpus = 50)
    saveRDS(SEM_distance_boot,paste0("../Data/SEM/SEM_boot_Nit_",k,"_","graz_",graz,".rds"))
    
  }
}


# ---------------------- Step 5: Effect on the resilience ----
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
post_param=read.table("../Data/Inferrence/Eby_1_neigh/param_inferred_1_neigh.csv",sep=";")

d_data=read.table("../Data/Spatial_structure_grazing.csv",sep=";")
d_data$bimod=sapply(1:nrow(d_data),function(x){
  if (dip.test(post_param[,x])$p.value<.05 | dip.test(post_param[,x+504])$p.value<.05){
    return("bimod")
  }else {
    return("unimod")
  }
})
write.table(which(d_data$bimod!="bimod"),"../Data/Inferrence/Eby_1_neigh/Keeping_sites_1_neigh.csv",sep=";")




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


## >> 5) Running mixed-effect models with aridity (ABC uncertainty) ----
n_neigh=1
d=read.table(paste0("../Data/Inferrence/Eby_",n_neigh,"_neigh/param_inferred_",n_neigh,"_neigh.csv"),sep=";",header=T)
d_data=read.table(paste0("../Data/Spatial_structure_grazing.csv"),sep=";")%>%
  Closer_to_normality(.)
n_sites=504


keep_sites=read.table(paste0("../Data/Inferrence/Eby_",n_neigh,"_neigh/Keeping_sites_",n_neigh,"_neigh.csv"),sep=";")$x

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

d2=read.table(paste0("../Data/Inferrence/Eby_",n_neigh,"_neigh/Resilience_metrics_",n_neigh,"_neigh.csv"),sep=";")%>%
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
  
  mod_predictors=gsub("\n     ","","Aridity + Grazing + Sand + Sp_richness + Org_C +
      Lat + Long_cos + Long_sin + Slope + Elevation + Grazing * Aridity + Grazing * Org_C + ( 1 | Plot_n)")
  
  for (Grazing_intensity in c("Low","High","Full")){
    
    if (Grazing_intensity=="Low"){
      Sub=d2%>%filter(., Grazing %in% c(0,1))
    }else if (Grazing_intensity=="High"){
      Sub=d2%>%filter(., Grazing %in% c(2,3))
    }else{
      Sub=d2
    }
    
    Sub$Grazing=scale(Sub$Grazing)[,1]
    
    model_q=lmer(formula = paste("q ~ ",mod_predictors), data = Sub,na.action = na.fail ,REML ="FALSE")
    model_abs=lmer(formula = paste("abs_dist ~ ",mod_predictors), data = Sub,na.action = na.fail ,REML ="FALSE")
    model_rela=lmer(formula = paste("rela_dist ~ ",mod_predictors), data = Sub,na.action = na.fail ,REML ="FALSE")
    
    
    #Getting partial prediction
    
    resid_mod_abs=visreg::visreg(fit = model_abs,xvar="Grazing",plot=F) 
    mod_abs=lm(visregRes~Grazing,resid_mod_abs$res)
    
    resid_mod_rela=visreg::visreg(fit = model_rela,xvar="Grazing",plot=F) 
    mod_rela=lm(visregRes~Grazing,resid_mod_rela$res)
    
    resid_mod_q=visreg::visreg(fit = model_q,xvar="Grazing",plot=F) 
    mod_q=lm(visregRes~Grazing,resid_mod_q$res)
    
    d_partial=rbind(d_partial,
                    tibble(slope=c(summary(mod_abs)$coefficient[2,1],
                                   summary(mod_rela)$coefficient[2,1],
                                   summary(mod_q)$coefficient[2,1]),
                           Stat=c("Distance to tipping (abs)","Distance to tipping (rela)","q (Spatial structure)"),
                           Type=Grazing_intensity,Metric="Grazing"))
    
    resid_mod_abs=visreg::visreg(fit = model_abs,xvar="Aridity",plot=F) 
    mod_abs=lm(visregRes~Aridity,resid_mod_abs$res)
    
    resid_mod_rela=visreg::visreg(fit = model_rela,xvar="Aridity",plot=F) 
    mod_rela=lm(visregRes~Aridity,resid_mod_rela$res)
    
    resid_mod_q=visreg::visreg(fit = model_q,xvar="Aridity",plot=F) 
    mod_q=lm(visregRes~Aridity,resid_mod_q$res)
    
    d_partial=rbind(d_partial,
                    tibble(slope=c(summary(mod_abs)$coefficient[2,1],
                                   summary(mod_rela)$coefficient[2,1],
                                   summary(mod_q)$coefficient[2,1]),
                           Stat=c("Distance to tipping (abs)","Distance to tipping (rela)","q (Spatial structure)"),
                           Type=Grazing_intensity,Metric="Aridity"))
    
    model_p=summary(model_q)
    model_q=summary(model_q)
    model_abs=summary(model_abs)
    model_rela=summary(model_rela)
    model_size=summary(model_q)
    
    
    d_mod=rbind(d_mod,tibble(Aridity=c(model_p$coefficients[2,1],model_q$coefficients[2,1],model_abs$coefficients[2,1],
                                       model_rela$coefficients[2,1],model_size$coefficients[2,1]),
                             Grazing=c(model_p$coefficients[3,1],model_q$coefficients[3,1],model_abs$coefficients[3,1],
                                       model_rela$coefficients[3,1],model_size$coefficients[3,1]),
                             Sand=c(model_p$coefficients[4,1],model_q$coefficients[4,1],model_abs$coefficients[4,1],
                                    model_rela$coefficients[4,1],model_size$coefficients[4,1]),
                             Sp_richness=c(model_p$coefficients[5,1],model_q$coefficients[5,1],model_abs$coefficients[5,1],
                                           model_rela$coefficients[5,1],model_size$coefficients[5,1]),
                             # Cover=c(model_p$coefficients[6,1],model_q$coefficients[6,1],model_abs$coefficients[6,1],
                             #         model_rela$coefficients[6,1],model_size$coefficients[6,1]),
                             Org_C=c(model_p$coefficients[6,1],model_q$coefficients[6,1],model_abs$coefficients[6,1],
                                     model_rela$coefficients[6,1],model_size$coefficients[6,1]),
                             Lat=c(model_p$coefficients[7,1],model_q$coefficients[7,1],model_abs$coefficients[7,1],
                                   model_rela$coefficients[7,1],model_size$coefficients[7,1]),
                             Long_cos=c(model_p$coefficients[8,1],model_q$coefficients[8,1],model_abs$coefficients[8,1],
                                        model_rela$coefficients[8,1],model_size$coefficients[8,1]),
                             Long_sin=c(model_p$coefficients[9,1],model_q$coefficients[9,1],model_abs$coefficients[9,1],
                                        model_rela$coefficients[9,1],model_size$coefficients[9,1]),
                             Slope=c(model_p$coefficients[10,1],model_q$coefficients[10,1],model_abs$coefficients[10,1],
                                     model_rela$coefficients[10,1],model_size$coefficients[10,1]),
                             Elevation=c(model_p$coefficients[11,1],model_q$coefficients[11,1],model_abs$coefficients[11,1],
                                         model_rela$coefficients[11,1],model_size$coefficients[11,1]),
                             Aridity_Grazing=c(model_p$coefficients[12,1],model_q$coefficients[12,1],model_abs$coefficients[12,1],
                                               model_rela$coefficients[12,1],model_size$coefficients[12,1]),
                             Grazing_Org_C=c(model_p$coefficients[13,1],model_q$coefficients[13,1],model_abs$coefficients[13,1],
                                             model_rela$coefficients[13,1],model_size$coefficients[13,1]),
                             Param=c("p","q","Absolute distance","Relative distance","Size tipping"),
                             Type_grazing=Grazing_intensity))
    
    
    
    
  }
  print(k)
  
}

write.table(d_mod,paste0("../Data/Inferrence/Eby_",n_neigh,"_neigh/Drivers_stability_metrics_no_cover_ABC_uncertainty.csv"),sep=";")
write.table(d_partial,paste0("../Data/Inferrence/Eby_",n_neigh,"_neigh/Partial_residuals_grazing_no_cover_ABC_uncertainty.csv"),sep=";")


## >> 6) Running mixed-effect models with aridity (Data uncertainty) ----

n_neigh=1
d=read.table(paste0("../Data/Inferrence/Eby_",n_neigh,"_neigh/param_inferred_",n_neigh,"_neigh.csv"),sep=";",header=T)
d_data=read.table(paste0("../Data/Spatial_structure_grazing.csv"),sep=";")%>%
  Closer_to_normality(.)
n_sites=504

keep_sites=read.table(paste0("../Data/Inferrence/Eby_",n_neigh,"_neigh/Keeping_sites_",n_neigh,"_neigh.csv"),sep=";")$x

d=tibble(Site=1:n_sites,mean_p=apply(d[,(1:n_sites)],2,mean),sd_p=apply(d[,(1:n_sites)],2,sd),
         mean_q=apply(d[,(n_sites+1):(2*n_sites)],2,mean),sd_q=apply(d[,(n_sites+1):(2*n_sites)],2,sd),
         median_q=apply(d[,(n_sites+1):(2*n_sites)],2,median),
         Plot_n=d_data$Site_ID,
         Aridity=d_data$Aridity,
         Clim1=d_data$Clim1,
         Clim2=d_data$Clim2,
         Clim3=d_data$Clim3,
         Clim4=d_data$Clim4,
         Sp_richness=d_data$Sp_richness,
         Type_veg=d_data$Type_veg,
         Org_C=d_data$Org_C,
         Org_C_v=d_data$Org_C_v,
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

d2=read.table(paste0("../Data/Inferrence/Eby_",n_neigh,"_neigh/Resilience_metrics_",n_neigh,"_neigh.csv"),sep=";")%>%
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

d2=tibble(p=scale(sapply(1:nrow(d),function(x){return(rnorm(1,d$mean_p[x],d$sd_p[x]))}))[,1],
          q=scale(logit(d$median_q))[,1],
          Size_tipping=scale(sapply(1:nrow(d),function(x){return(rnorm(1,d$Size_mean[x],d$Size_sd[x]))}))[,1],
          abs_dist=scale(d$abs_median)[,1],
          rela_dist=scale(d$relativ_median)[,1],
          
          #site related variables
          Sand=(d$Sand-mean(d$Sand,na.rm=T))/sd(d$Sand,na.rm = T),
          Sp_richness=(d$Sp_richness-mean(d$Sp_richness,na.rm=T))/sd(d$Sp_richness,na.rm = T),
          Grazing=d$Grazing,
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

d_partial=d_mod=tibble()
for (Grazing_intensity in c("Low","High","Full")){
  
  if (Grazing_intensity=="Low"){
    Sub=d2%>%filter(., Grazing %in% c(0,1))
  }else if (Grazing_intensity=="High"){
    Sub=d2%>%filter(., Grazing %in% c(2,3))
  }else{
    Sub=d2
  }
  
  mod_predictors=gsub("\n     ","","Aridity + Grazing + Sand + Sp_richness+ Org_C_v + Org_C_v * Grazing +
        Lat + Long_cos + Long_sin + Slope + Elevation + Grazing * Aridity + ( 1 | Plot_n)")
  
  model_q=lmer(formula = paste("q ~ ",mod_predictors), data = Sub,na.action = na.fail ,REML ="FALSE")
  model_abs=lmer(formula = paste("abs_dist ~ ",mod_predictors), data = Sub,na.action = na.fail ,REML ="FALSE")
  model_rela=lmer(formula = paste("rela_dist ~ ",mod_predictors), data = Sub,na.action = na.fail ,REML ="FALSE")
  model_cover=lmer(formula = paste("Cover ~ ",mod_predictors), data = Sub,na.action = na.fail ,REML ="FALSE")
  
  #Getting partial prediction
  boot_function = function(formula, data, sampled_id) {
    d = data[sampled_id,] 
    fit = lm(formula, data=d) 
    return(coef(fit)[2]) #return coefficient estimates of model
  }
  
  #GRAZING
  
  resid_mod_abs=visreg::visreg(fit = model_abs,xvar="Grazing",plot=F) 
  model_abs_boot=boot(data=resid_mod_abs$res,statistic=boot_function,R=500, formula=visregRes~Grazing)$t[,1]
  
  resid_mod_rela=visreg::visreg(fit = model_rela,xvar="Grazing",plot=F) 
  model_rela_boot=boot(data=resid_mod_rela$res,statistic=boot_function,R=500, formula=visregRes~Grazing)$t[,1]
  
  resid_mod_q=visreg::visreg(fit = model_q,xvar="Grazing",plot=F) 
  model_q_boot=boot(data=resid_mod_q$res,statistic=boot_function,R=500, formula=visregRes~Grazing)$t[,1]

  resid_mod_cover=visreg::visreg(fit = model_cover,xvar="Grazing",plot=F) 
  model_cover_boot=boot(data=resid_mod_cover$res,statistic=boot_function,R=500, formula=visregRes~Grazing)$t[,1]
  
  d_mod=rbind(d_mod,data.frame(Stat=model_abs_boot)%>%
                add_column(., Param="Absolute distance",Type=Grazing_intensity,Metric="Grazing"))
  d_mod=rbind(d_mod,data.frame(Stat=model_rela_boot)%>%
                add_column(., Param="Relative distance",Type=Grazing_intensity,Metric="Grazing"))
  d_mod=rbind(d_mod,data.frame(Stat=model_q_boot)%>%
                add_column(., Param="q",Type=Grazing_intensity,Metric="Grazing"))
  d_mod=rbind(d_mod,data.frame(Stat=model_cover_boot)%>%
                add_column(., Param="Cover",Type=Grazing_intensity,Metric="Grazing"))
  
  
  #ARIDITY
  resid_mod_abs=visreg::visreg(fit = model_abs,xvar="Aridity",plot=F) 
  model_abs_boot=boot(data=resid_mod_abs$res,statistic=boot_function,R=500, formula=visregRes~Aridity)$t[,1]
  
  resid_mod_rela=visreg::visreg(fit = model_rela,xvar="Aridity",plot=F) 
  model_rela_boot=boot(data=resid_mod_rela$res,statistic=boot_function,R=500, formula=visregRes~Aridity)$t[,1]
  
  resid_mod_q=visreg::visreg(fit = model_q,xvar="Aridity",plot=F) 
  model_q_boot=boot(data=resid_mod_q$res,statistic=boot_function,R=500, formula=visregRes~Aridity)$t[,1]
  
  resid_mod_cover=visreg::visreg(fit = model_cover,xvar="Aridity",plot=F) 
  model_cover_boot=boot(data=resid_mod_cover$res,statistic=boot_function,R=500, formula=visregRes~Aridity)$t[,1]
  
  d_mod=rbind(d_mod,data.frame(Stat=model_abs_boot)%>%
                add_column(., Param="Absolute distance",Type=Grazing_intensity,Metric="Aridity"))
  d_mod=rbind(d_mod,data.frame(Stat=model_rela_boot)%>%
                add_column(., Param="Relative distance",Type=Grazing_intensity,Metric="Aridity"))
  d_mod=rbind(d_mod,data.frame(Stat=model_q_boot)%>%
                add_column(., Param="q",Type=Grazing_intensity,Metric="Aridity"))
  d_mod=rbind(d_mod,data.frame(Stat=model_cover_boot)%>%
                add_column(., Param="Cover",Type=Grazing_intensity,Metric="Aridity"))
}

write.table(d_mod,paste0("../Data/Inferrence/Eby_",n_neigh,"_neigh/Partial_residuals_grazing_no_cover_data_uncertainty.csv"),sep=";")



## >> 7) SEM on the resilience of drylands ----

n_neigh=1
d=read.table(paste0("../Data/Inferrence/Eby_",n_neigh,"_neigh/param_inferred_",n_neigh,"_neigh.csv"),sep=";",header=T)
d_data=read.table(paste0("../Data/Spatial_structure_grazing.csv"),sep=";")%>%
  Closer_to_normality(.)
n_sites=504
keep_sites=read.table(paste0("../Data/Inferrence/Eby_",n_neigh,"_neigh/Keeping_sites_",n_neigh,"_neigh.csv"),sep=";")$x

d=tibble(Site=1:n_sites,mean_p=apply(d[,(1:n_sites)],2,mean),sd_p=apply(d[,(1:n_sites)],2,sd),
         mean_q=apply(d[,(n_sites+1):(2*n_sites)],2,mean),sd_q=apply(d[,(n_sites+1):(2*n_sites)],2,sd),
         median_q=apply(d[,(n_sites+1):(2*n_sites)],2,median),
         Site_ID=d_data$Site_ID,
         Aridity=d_data$Aridity,
         Clim1=d_data$Clim1,
         Clim2=d_data$Clim2,
         Clim3=d_data$Clim3,
         Clim4=d_data$Clim4,
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

d2=read.table(paste0("../Data/Inferrence/Eby_",n_neigh,"_neigh/Resilience_metrics_",n_neigh,"_neigh.csv"),sep=";")%>%
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

d2=tibble(p=scale(sapply(1:nrow(d),function(x){return(rnorm(1,d$mean_p[x],d$sd_p[x]))}))[,1],
          q=scale(logit(d$median_q))[,1],
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

for (k in 1:3){
  if (k==1){d_sem=save%>%filter(., Grazing %in% c(0,1))
  }else if (k==2){
    d_sem=save%>%filter(., Grazing %in% c(2,3))
  }else{
    d_sem=save
  }
  
  SEM_distance_nit=psem(
    lmer(Resid_mod ~ (1|Site_ID) + Cover + q, d_sem),
    lmer(Cover ~ (1|Site_ID) + Aridity + Sand  + Org_C_v + Grazing, d_sem),
    lmer(q ~ (1|Site_ID) + Aridity + Nitrate + Sand  + Org_C_v + Grazing, d_sem),
    lmer(Nitrate ~ (1|Site_ID) +  Aridity + Sand + Grazing, d_sem),
    lmer(Org_C_v ~ (1|Site_ID) +  Aridity + Sand + Grazing+Nitrate, d_sem),
    q%~~%Cover,Sand%~~%Aridity,
    Cover%~~%Nitrate
  )
  print(summary(SEM_distance_nit))
  
  SEM_distance_NP=psem(
    lmer(Resid_mod ~ (1|Site_ID) + Cover + q, d_sem),
    lmer(Cover ~ (1|Site_ID) + Aridity  + Sand  + Org_C_v + Grazing, d_sem),
    lmer(q ~ (1|Site_ID) + Aridity  + Sand  + Total_N + Org_C_v  + Grazing, d_sem),
    lmer(Total_N ~ (1|Site_ID) +  Aridity + Sand + Grazing, d_sem),
    lmer(Org_C_v ~ (1|Site_ID) + Aridity + Sand + Grazing +Total_N, d_sem),
    lm(Sand ~ Aridity, d_sem),
    q%~~%Cover,
    Cover%~~%Total_N
  )
  #print(summary(SEM_distance_NP))
}

# ---------------------- Step 6: Island of fertility ----

## >> 1) Running the models ----
Run_model_fertility_nutrients=function(id){
  
  list_mod=expand.grid(with_cover=c(T),
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
      + Aridity*Grazing 
      + rho_p + Type_veg 
      + Sand 
      + (1|Site_ID)")))
    
  } else{
    
    formula_mod=formula(formula_(paste("value ~ Long_cos + Long_sin + Lattitude + Slope + Elevation 
      + Aridity + Grazing 
      + Aridity*Grazing 
      + Type_veg 
      + Sand 
      + (1|Site_ID)")))
  }
  
  
  save=d_data_mod
  
  for(grazing_intensity in c("low","high","all")){
    
    if (grazing_intensity =="low"){
      d_data_mod=save%>%filter(., Grazing %in% c(0,1))
    } else if (grazing_intensity =="high"){
      d_data_mod=save%>%filter(., Grazing %in% c(2,3))
    } else{
      d_data_mod=save
    }
    
    d_data_mod$Grazing=scale(d_data_mod$Grazing)[,1]
    
    model_spa_stat  = lmer(formula_mod, d_data_mod,
                           na.action = na.fail ,REML ="FALSE")
    
    #we remove potential outliers
    rm.outliers = romr.fnc(model_spa_stat, d_data_mod, trim=2.5)
    d_data_out = rm.outliers$data
    
    #saving data
    write.table(d_data_out,paste0("../Data/Linear_models_nutrients/Keep_data/Data_",stat,"_",with_cover,"_aridity_",grazing_intensity,".csv"))
    
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
    
    saveRDS(model_spa_stat,paste0("../Data/Linear_models_nutrients/Keep_models/Mod_",stat,"_",with_cover,"_aridity_",grazing_intensity,".rds"))
    
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
  write.table(d_all,paste0("../Data/Linear_models_nutrients/Importance/Importance_",stat,"_",with_cover,"_aridity.csv"),sep=";")
  write.table(d_all2,paste0("../Data/Linear_models_nutrients/Estimator/Estimators_model_",stat,"_",with_cover,"_aridity.csv"),sep=";")
  
}

library(parallel)

mclapply(1:30,Run_model_fertility_nutrients,mc.cores = 1)

d_all=d_all2=tibble()
for (k in list.files("../Data/Linear_models_nutrients/Importance",paste0("aridity.csv"))){
  d=read.table(paste0("../Data/Linear_models_nutrients/Importance/",k),sep=";")
  if (nrow(d)!=0){
    d_all=rbind(d_all,d)
  }
}
for (k in list.files("../Data/Linear_models_nutrients/Estimator",paste0("aridity.csv"))){
  d2=read.table(paste0("../Data/Linear_models_nutrients/Estimator/",k),sep=";")
  if (nrow(d2)!=0){
    d_all2=rbind(d_all2,d2)
  }
}
write.table(d_all,paste0("../Data/Linear_models_nutrients/Importance_nutrients.csv"),sep=";")
write.table(d_all2,paste0("../Data/Linear_models_nutrients/Estimators_model_nutrients.csv"),sep=";")



## >> 2) Extracting residuals ----

d_data=read.table("../Data/Spatial_structure_grazing.csv",sep=";")%>%
  Closer_to_normality(.)
d_data[,c(1,9:29,36,38:42,45:ncol(d_data))] = 
  apply(d_data[,c(1,9:29,36,38:42,45:ncol(d_data))],2,z_tranform)
d_data=Perform_PCA_spatial_struc(d_data)
save=d_data
d_slope=tibble()

for (k in c("Org_C","Org_C_v","Total_N","lnTotal_N","Nitrate","lnNitrate")){
  
  for (grazing_intensity in c("low","high","all")){
    
    #for each we plot the slope against the partial residuals with and without cover
    
     d_data_out=read.table(paste0("../Data/Linear_models_nutrients/Keep_data/Data_",k,
                                 "_TRUE_aridity_",grazing_intensity,".csv"),sep=" ")
    model_spa_stat=readRDS(paste0("../Data/Linear_models_nutrients/Keep_models/Mod_",k,
                                  "_TRUE_aridity_",grazing_intensity,".rds"))
    
    resid_mod=visreg::visreg(fit = model_spa_stat,xvar="Grazing",plot=F) 
    
    mod_cov = boot(data=resid_mod$res, statistic=boot_function_lm, R=1000, formula=visregRes~Grazing)
    
    d_slope=rbind(d_slope,
                  tibble(pval=twoside_pvalue(mod_cov$t),
                         q2=median(mod_cov$t),
                         q1=quantile(mod_cov$t,.025),
                         q3=quantile(mod_cov$t,.975),
                         With_cover=T,Grazing_intensity=grazing_intensity,
                         Stat=k,Driver="Grazing"
                  ))
    
    resid_mod=visreg::visreg(fit = model_spa_stat,xvar="Aridity",plot=F) 
    
    mod_cov = boot(data=resid_mod$res, statistic=boot_function_lm, R=1000, formula=visregRes~Aridity)
    
    d_slope=rbind(d_slope,
                  tibble(pval=twoside_pvalue(mod_cov$t),
                         q2=median(mod_cov$t),
                         q1=quantile(mod_cov$t,.025),
                         q3=quantile(mod_cov$t,.975),
                         With_cover=T,Grazing_intensity=grazing_intensity,
                         Stat=k,Driver="Aridity"
                  ))

  }
}

write.table(d_slope,"../Data/Linear_models_nutrients/Slope_partial_residuals_aridity_grazing.csv",sep=";")


# ---------------------- Step 7: Grazing as a category ----
## >> 1) Run importance ----

Run_model_importance_aridity_no_inter=function(id){
  
  list_mod=rbind(expand.grid(with_cover=c(T,F),
                             Stats=c("perim_area_scaling","fmax_psd","PL_expo","Spectral_ratio",
                                     "core_area_land","division","fractal_dim","contig","core_area",
                                     "flow_length","PLR","KS_dist",
                                     "Struct1","Struct2")), #We also add the case where the stats correspond to the first two components of the PCA
                 tibble(with_cover=F,Stats="rho_p"))
  
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
  
  for(grazing_intensity in c("low","high","all")){
    
    if (grazing_intensity =="low"){
      d_data_mod=save%>%filter(., Grazing %in% c(0,1))
    } else if (grazing_intensity =="high"){
      d_data_mod=save%>%filter(., Grazing %in% c(2,3))
    } else{
      d_data_mod=save
    }
    
    d_data_mod$Grazing=scale(d_data_mod$Grazing)[,1]
    
    model_spa_stat  = lmer(formula_mod, d_data_mod,
                           na.action = na.fail ,REML ="FALSE")
    
    #we remove potential outliers
    rm.outliers = romr.fnc(model_spa_stat, d_data_mod, trim=2.5)
    d_data_out = rm.outliers$data
    
    #saving data
    write.table(d_data_out,paste0("../Data/Linear_models/Keep_data/Data_",stat,"_",with_cover,"_aridity_no_inter_",grazing_intensity,".csv"))
    
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
    
    if (dim(select_model%>%filter(., delta<2))[1]==1){ #one model
      
      result_select=model_spa_stat
      R2=select_model%>%filter(., AICc<min(AICc)+2)
      
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
    
    saveRDS(model_spa_stat,paste0("../Data/Linear_models/Keep_models/Mod_",stat,"_",with_cover,"_aridity_no_inter_",grazing_intensity,".rds"))
    
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
  write.table(d_all,paste0("../Data/Linear_models/Importance/Importance_",stat,"_",with_cover,"_aridity_no_inter.csv"),sep=";")
  write.table(d_all2,paste0("../Data/Linear_models/Estimator/Estimators_model_",stat,"_",with_cover,"_aridity_no_inter.csv"),sep=";")
  
}

library(parallel)

mclapply(1:30,Run_model_importance_aridity_no_inter,mc.cores = 1)

d_all=d_all2=tibble()
for (k in list.files("../Data/Linear_models/Importance",paste0("aridity_no_inter"))){
  d=read.table(paste0("../Data/Linear_models/Importance/",k),sep=";")
  if (nrow(d)!=0){
    d_all=rbind(d_all,d)
  }
}
for (k in list.files("../Data/Linear_models/Estimator",paste0("aridity_no_inter"))){
  d2=read.table(paste0("../Data/Linear_models/Estimator/",k),sep=";")
  if (nrow(d)!=0){
    d_all2=rbind(d_all2,d2)
  }
}
write.table(d_all,paste0("../Data/Linear_models/Importance_","aridity_no_inter.csv"),sep=";")
write.table(d_all2,paste0("../Data/Linear_models/Estimators_model_","aridity_no_inter.csv"),sep=";")



library(parallel)

mclapply(1:28,Run_model_importance,mc.cores = 28)


d_all=d_all2=tibble()
for (k in list.files("../Data/Linear_models/Grazing_factor/Importance","TRUE")){
  d_all=rbind(d_all,read.table(paste0("../Data/Linear_models/Grazing_factor/Importance/",k),sep=";"))
}
for (k in list.files("../Data/Linear_models/Grazing_factor/Estimator","TRUE")){
  d_all2=rbind(d_all2,read.table(paste0("../Data/Linear_models/Grazing_factor/Estimator/",k),sep=";"))
}
write.table(d_all,"../Data/Linear_models/Grazing_factor/Importance.csv",sep=";")
write.table(d_all2,"../Data/Linear_models/Grazing_factor/Estimators_model.csv",sep=";")



d_all=d_all2=tibble()
for (k in list.files("../Data/Linear_models/Grazing_factor/Importance","FALSE")){
  d_all=rbind(d_all,read.table(paste0("../Data/Linear_models/Grazing_factor/Importance/",k),sep=";"))
}
for (k in list.files("../Data/Linear_models/Grazing_factor/Estimator","FALSE")){
  d_all2=rbind(d_all2,read.table(paste0("../Data/Linear_models/Grazing_factor/Estimator/",k),sep=";"))
}
write.table(d_all,"../Data/Linear_models/Grazing_factor/Importance_no_cov.csv",sep=";")
write.table(d_all2,"../Data/Linear_models/Grazing_factor/Estimators_model_no_cov.csv",sep=";")


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
  
  d_data_out=read.table(paste0("../Data/Linear_models/Grazing_factor/Keep_data/Data_",k,"_TRUE.csv"),sep=" ")
  model_spa_stat=readRDS(paste0("../Data/Linear_models/Grazing_factor/Keep_models/Mod_",k,"_TRUE.rds"))
  d_data$Grazing=as.numeric(d_data$Grazing)
  d_data_out=read.table(paste0("../Data/Linear_models/Keep_data/Data_",k,"_TRUE.csv"),sep=" ")
  model_spa_stat2=readRDS(paste0("../Data/Linear_models/Keep_models/Mod_",k,"_TRUE.rds"))
  d_AIC=rbind(d_AIC,tibble(AIC_factor=AIC(model_spa_stat,model_spa_stat2)[1,2],AIC_num=AIC(model_spa_stat,model_spa_stat2)[2,2])%>%
                add_column(., Stat=k))
  
}
d_AIC

## >> 3) Analyse residuals ----

pdf("../Figures/Linear_models/Grazing_factor/Test_factor_with_cov.pdf",width = 6,height = 4)


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
  
  d_data_out=read.table(paste0("../Data/Linear_models/Grazing_factor/Keep_data/Data_",k,"_TRUE.csv"),sep=" ")
  d_data_out$Grazing=as.factor(d_data_out$Grazing)
  model_spa_stat=readRDS(paste0("../Data/Linear_models/Grazing_factor/Keep_models/Mod_",k,"_TRUE.rds"))
  
  
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


pdf("../Figures/Linear_models/Grazing_factor/Test_factor_no_cov.pdf",width = 6,height = 4)


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
  
  
  d_data_out=read.table(paste0("../Data/Linear_models/Grazing_factor/Keep_data/Data_",k,"_FALSE.csv"),sep=" ")
  d_data_out$Grazing=as.factor(d_data_out$Grazing)
  model_spa_stat=readRDS(paste0("../Data/Linear_models/Grazing_factor/Keep_models/Mod_",k,"_FALSE.rds"))
  
  
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
  
  d_data_out=read.table(paste0("../Data/Linear_models/Grazing_factor/Keep_data/Data_",k,"_TRUE.csv"),sep=" ")
  d_data_out$Grazing=as.factor(d_data_out$Grazing)
  model_spa_stat=readRDS(paste0("../Data/Linear_models/Grazing_factor/Keep_models/Mod_",k,"_TRUE.rds"))
  
  
  resid_mod=visreg::visreg(fit = model_spa_stat,xvar="Grazing",plot=F,ylab=as.character(k)) 
  
  mod_res=lm("visregRes ~ Grazing",resid_mod$res)
  
  d_resid=rbind(d_resid,tibble(Low=confint(mod_res)[-1,1],
                               High=confint(mod_res)[-1,2],
                               Stat=k,Name=c("Low","Medium","High"),Cover=T))
  
  
  d_data$Grazing=as.factor(d_data$Grazing)
  
  d_data_out=read.table(paste0("../Data/Linear_models/Grazing_factor/Keep_data/Data_",k,"_FALSE.csv"),sep=" ")
  d_data_out$Grazing=as.factor(d_data_out$Grazing)
  model_spa_stat=readRDS(paste0("../Data/Linear_models/Grazing_factor/Keep_models/Mod_",k,"_FALSE.rds"))
  
  
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

ggsave("../Figures/Linear_models/Grazing_factor/Residuals_with_factor_grazing.pdf",p,width = 7,height = 5)


