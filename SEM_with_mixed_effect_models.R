rm(list=ls())
source("./Structure_grazing_function.R")


d_data=read.table("../Data/Spatial_structure_grazing.csv",sep=";")%>%
  Closer_to_normality(.)

d_data[,c(1,9:29,35,37:41,44:ncol(d_data))] = apply(d_data[,c(1,9:29,35,37:41,44:ncol(d_data))],2,z_tranform)

#Getting the multifunctionality index
d_data$MF=z_tranform(rowSums(d_data[,55:(ncol(d_data))],na.rm = T))

#Adding the first two components from the PCA
d_data=Perform_PCA_spatial_struc(d_data)

d_indirect=tibble() #to save indirect effects of grazing on the spatial structure
dir.create("../Figures/Step1_Understanding_grazing/SEM_lmer/",showWarnings = F)

param_list=expand.grid(Cov=c("with_cover","without_cover"),
                       Stat=c("perim_area_scaling","PL_expo","Cond_H","fmax_psd",
                              "PLR","flow_length","mean_perim_area","core_area_land",
                              "fractal_dim","division","contig",
                              "Struct1","Struct2"),
                       MF=c(T,F)) #MF instead of organic carbon

for (ID in 1:nrow(param_list)){
  
  dir.create(paste0("../Figures/Step1_Understanding_grazing/SEM_lmer/",param_list$Stat[ID]),showWarnings = F)
  
  with_cover=param_list$Cov[ID];k=as.character(param_list$Stat[ID])
  
  model_lmer=lmer(paste("value ~ (1|Site_ID) + Long_cos + Long_sin + Lattitude + Slope + Elevation + Clim1 + Clim2 + Clim3 + Clim4 + Type_veg",
                      ifelse(with_cover=="without_cover","+ rho_p","")),REML = F,
                data = d_data%>%melt(., measure.vars=k)%>%filter(., !is.na(value)),
                na.action = na.omit)
  
  #we remove potential outliers
  rm.outliers = romr.fnc(model_lmer, d_data%>%melt(., measure.vars=k)%>%filter(., !is.na(value)),
                         trim=2.5)
  d_data_out = rm.outliers$data
  
  #We first control for all covariates and extract the residuals
  model_lmer=lmer(paste("value ~ (1|Site_ID) + Long_cos + Long_sin + Lattitude + Slope + Elevation + Clim1 + Clim2 + Clim3 + Clim4 + Type_veg",
                      ifelse(with_cover=="without_cover","+ rho_p","")),
                data = d_data_out,REML = F,
                na.action = na.fail)
  
  resid_model=residuals(model_lmer) #extract residuals
  
  save=d_data[as.numeric(names(resid_model)),]%>% #add it to the dataframe 
    add_column(., Resid_mod=resid_model)
  
  if (param_list$MF[ID]){ save$Org_C=save$MF}
  
  #DOING the SEMs
  
  if (with_cover=="with_cover"){
    
    #SEM with all grazing intensity
    
    d_sem=save
    all_d=summary(psem(
      lmer(Resid_mod ~ (1|Site_ID) + Grazing + Woody + rho_p , d_sem),
      lmer(Woody ~ (1|Site_ID) + Grazing + Org_C + Sand, d_sem),
      lmer(Sand ~ (1|Site_ID) + Grazing , d_sem),
      lmer(rho_p ~ (1|Site_ID) + Grazing + Org_C + Sand , d_sem),
      lmer(Org_C ~ (1|Site_ID) + Grazing, d_sem)
    ))
    
    #SEM with low grazing intensity
    
    d_sem=filter(save,Grazing %in% 0:1)
    low_graz=summary(psem(
      lmer(Resid_mod ~ (1|Site_ID) + Grazing + Woody + rho_p , d_sem),
      lmer(Woody ~ (1|Site_ID) + Grazing + Org_C + Sand, d_sem),
      lmer(Sand ~ (1|Site_ID) + Grazing , d_sem),
      lmer(rho_p ~ (1|Site_ID) + Grazing + Org_C + Sand , d_sem),
      lmer(Org_C ~ (1|Site_ID) + Grazing, d_sem)
    ))
    
    #SEM with high grazing intensity
    
    d_sem=filter(save,Grazing %in% 2:3)
    high_graz=summary(psem(
      lmer(Resid_mod ~ (1|Site_ID) + Grazing + Woody + rho_p , d_sem),
      lmer(Woody ~ (1|Site_ID) + Grazing + Org_C + Sand, d_sem),
      lmer(Sand ~ (1|Site_ID) + Grazing , d_sem),
      lmer(rho_p ~ (1|Site_ID) + Grazing + Org_C + Sand , d_sem),
      lmer(Org_C ~ (1|Site_ID) + Grazing, d_sem)
    ))
    
    
  }else {
    
    #SEM with all grazing intensity
    
    d_sem=save
    all_d=summary(psem(
      lmer(Resid_mod ~ (1|Site_ID) + Grazing + Woody , d_sem),
      lmer(Woody ~ (1|Site_ID) + Grazing + Org_C + Sand, d_sem),
      lmer(Sand ~ (1|Site_ID) + Grazing , d_sem),
      lmer(Org_C ~ (1|Site_ID) + Grazing, d_sem)
    ))
    
    #SEM with low grazing intensity
    
    d_sem=filter(save,Grazing %in% 0:1)
    low_graz=summary(psem(
      lmer(Resid_mod ~ (1|Site_ID) + Grazing + Woody , d_sem),
      lmer(Woody ~ (1|Site_ID) + Grazing + Org_C + Sand, d_sem),
      lmer(Sand ~ (1|Site_ID) + Grazing , d_sem),
      lmer(Org_C ~ (1|Site_ID) + Grazing, d_sem)
    ))
    
    #SEM with high grazing intensity
    
    d_sem=filter(save,Grazing %in% 2:3)
    high_graz=summary(psem(
      lmer(Resid_mod ~ (1|Site_ID) + Grazing + Woody , d_sem),
      lmer(Woody ~ (1|Site_ID) + Grazing + Org_C + Sand, d_sem),
      lmer(Sand ~ (1|Site_ID) + Grazing , d_sem),
      lmer(Org_C ~ (1|Site_ID) + Grazing, d_sem)
    ))
  }
  
  
  d_indirect=rbind(d_indirect,
                   Get_indirect_effects_grazing(all_d)%>%
                     add_column(., Type="All",Stat=k,With_cover=with_cover,MF=param_list$MF[ID]),
                   Get_indirect_effects_grazing(high_graz)%>%
                     add_column(., Type="High",Stat=k,With_cover=with_cover,MF=param_list$MF[ID]),
                   Get_indirect_effects_grazing(low_graz)%>%
                     add_column(., Type="Low",Stat=k,With_cover=with_cover,MF=param_list$MF[ID]))
  
  #ploting
  
  if(with_cover=="with_cover"){
    
    Plot_SEM_with_cover(summary_sem = all_d,pdf_ = T,
                        name = paste0("../Figures/Step1_Understanding_grazing/SEM_lmer/",param_list$Stat[ID],"/",
                                      "SEM_",param_list$Stat[ID],"_",param_list$Cov[ID],"_MF_",param_list$MF[ID],"_all"),
                        title_ = paste0("SEM_",k,"_all"),
                        name_var = k,MF = param_list$MF[ID])
    
    Plot_SEM_with_cover(summary_sem = low_graz,pdf_ = T,
                        name = paste0("../Figures/Step1_Understanding_grazing/SEM_lmer/",param_list$Stat[ID],"/",
                                      "SEM_",param_list$Stat[ID],"_",param_list$Cov[ID],"_MF_",param_list$MF[ID],"_low"),
                        title_ = paste0("SEM_",k,"_low"),
                        name_var = k,MF = param_list$MF[ID])
    
    Plot_SEM_with_cover(summary_sem = high_graz,pdf_ = T,
                        name = paste0("../Figures/Step1_Understanding_grazing/SEM_lmer/",param_list$Stat[ID],"/",
                                      "SEM_",param_list$Stat[ID],"_",param_list$Cov[ID],"_MF_",param_list$MF[ID],"_high"),
                        title_ = paste0("SEM_",k,"_high"),
                        name_var = k,MF = param_list$MF[ID])
  }else {
    
    Plot_SEM_without_cover(summary_sem = all_d,pdf_ = T,
                           name = paste0("../Figures/Step1_Understanding_grazing/SEM_lmer/",param_list$Stat[ID],"/",
                                         "SEM_",param_list$Stat[ID],"_",param_list$Cov[ID],"_MF_",param_list$MF[ID],"_all"),
                           title_ = paste0("SEM_",k,"_all"),
                           name_var = k,MF = param_list$MF[ID])
    
    Plot_SEM_without_cover(summary_sem = low_graz,pdf_ = T,
                           name = paste0("../Figures/Step1_Understanding_grazing/SEM_lmer/",param_list$Stat[ID],"/",
                                         "SEM_",param_list$Stat[ID],"_",param_list$Cov[ID],"_MF_",param_list$MF[ID],"_low"),
                           title_ = paste0("SEM_",k,"_low"),
                           name_var = k,MF = param_list$MF[ID])
    
    Plot_SEM_without_cover(summary_sem = high_graz,pdf_ = T,
                           name = paste0("../Figures/Step1_Understanding_grazing/SEM_lmer/",param_list$Stat[ID],"/",
                                         "SEM_",param_list$Stat[ID],"_",param_list$Cov[ID],"_MF_",param_list$MF[ID],"_high"),
                           title_ = paste0("SEM_",k,"_high"),
                           name_var = k,MF = param_list$MF[ID])
  }
  
}

