# ---------------------- Step 2: Running mixed-effect models  ----

Run_model_importance_aridity=function(id){
  
  list_mod=rbind(expand.grid(with_cover=c(T,F),
                             Stats=c("perim_area_scaling","fmax_psd","PL_expo","Spectral_ratio",
                                     "core_area_land","division","fractal_dim","contig","core_area",
                                     "flow_length","PLR","KS_dist","Shape_metric",
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
      + Sand + Org_C_v 
      + Grazing * Org_C_v
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
  
  for(grazing_intensity in c("low","high","all","grazed")){
    
    if (grazing_intensity =="low"){
      d_data_mod=save%>%filter(., Grazing %in% c(0,1))
    } else if (grazing_intensity =="high"){
      d_data_mod=save%>%filter(., Grazing %in% c(2,3))
    } else if (grazing_intensity =="grazed"){
      d_data_mod=save%>%filter(., Grazing %in% c(1,2,3))
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
                               "flow_length","PLR","KS_dist","Shape_metric",
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
  
  for(grazing_intensity in c("low","high","all","grazed")){
    
    if (grazing_intensity =="low"){
      d_data_mod=save%>%filter(., Grazing %in% c(0,1))
    } else if (grazing_intensity =="high"){
      d_data_mod=save%>%filter(., Grazing %in% c(2,3))
    } else if (grazing_intensity =="grazed"){
      d_data_mod=save%>%filter(., Grazing %in% c(1,2,3))
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
                          dc(Org_C_v & Grazing, Org_C_v : Grazing) &
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
## >> 1) Model with aridity and grazing ----

d_data=read.table("../Data/Spatial_structure_grazing.csv",sep=";")%>%
  Closer_to_normality(.)
d_data[,c(1,9:29,36,38:42,45:ncol(d_data))] = 
  apply(d_data[,c(1,9:29,36,38:42,45:ncol(d_data))],2,z_tranform)
d_data=Perform_PCA_spatial_struc(d_data)
save=d_data
d_slope=tibble()

for (k in c("perim_area_scaling","fmax_psd","PL_expo","Spectral_ratio",
            "core_area_land","division","contig","core_area",
            "flow_length","PLR","KS_dist","Shape_metric",
            "Struct1","Struct2","rho_p","mean_psd","cv_psd")){
  
  for (grazing_intensity in c("low","high","all","grazed")){
    
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

## >> 2) Model with aridity and grazing but without interactions ----

d_data=read.table("../Data/Spatial_structure_grazing.csv",sep=";")%>%
  Closer_to_normality(.)
d_data[,c(1,9:29,36,38:42,45:ncol(d_data))] = 
  apply(d_data[,c(1,9:29,36,38:42,45:ncol(d_data))],2,z_tranform)
d_data=Perform_PCA_spatial_struc(d_data)
save=d_data
d_slope=tibble()

for (k in c("perim_area_scaling","fmax_psd","PL_expo","Spectral_ratio",
            "core_area_land","division","contig","core_area",
            "flow_length","PLR","KS_dist","Shape_metric",
            "Struct1","Struct2","rho_p","mean_psd","cv_psd")){
  
  for (grazing_intensity in c("low","high","all","grazed")){
    
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
    
    d_sem=save#%>%filter(., Grazing %in% c(0,1))
    SEM_distance=((psem(
      lmer(Resid_mod ~  (1|Site_ID) + Grazing + Org_C_v +Sand+Aridity+Nitrate+Amonium, d_sem),
      lmer(Nitrate ~  (1|Site_ID) + Grazing +Aridity+Sand, d_sem),
      lmer(Amonium ~  (1|Site_ID) + Grazing +Aridity+Sand, d_sem),
      lmer(Org_C_v ~  (1|Site_ID) + Grazing + Nitrate +Amonium + Aridity+Sand, d_sem),
      Nitrate%~~%Amonium
    )))
    
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
  }else if (Grazing_intensity=="grazed"){
    Sub=d2%>%filter(., Grazing %in% c(1,2,3))
  }else{
    Sub=d2
  }
  
  mod_predictors=gsub("\n     ","","Aridity + Grazing + Sand + Sp_richness+ Org_C_v + 
        Lat + Long_cos + Long_sin + Slope + Elevation + ( 1 | Plot_n)")
  
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
                       Stats=c("Org_C","Org_C_v","Total_N","lnTotal_N","Total_P","lnTotal_P"))
  
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

mclapply(1:6,Run_model_fertility_nutrients,mc.cores = 1)

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

for (k in c("Org_C","Org_C_v","Total_N","lnTotal_N","Total_P","lnTotal_P")){
  
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


# ---------------------- Step 7: Moving average aridity & sand ----


Run_model_moving_average_aridity=function(id){
  
  list_mod=rbind(expand.grid(with_cover=c(T,F),
                             Stats=c("perim_area_scaling","fmax_psd","PL_expo","Spectral_ratio",
                                     "core_area_land","division","fractal_dim","contig","core_area",
                                     "flow_length","PLR","KS_dist","mean_psd",
                                     "Struct1","Struct2")), #We also add the case where the stats correspond to the first two components of the PCA
                 tibble(with_cover=F,Stats="rho_p"))
  
  d_slope=tibble()#for keeping all information
  
  d_data=read.table("../Data/Spatial_structure_grazing.csv",sep=";")%>%
    Closer_to_normality(.)
  
  d_data[,c(1,9:29,38:42,45:ncol(d_data))] = 
    apply(d_data[,c(1,9:29,38:42,45:ncol(d_data))],2,z_tranform)
  
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
    d_data_mod$Grazing=scale(d_data_mod$Grazing)[,1]
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

mclapply(1:30,Run_model_moving_average_aridity,mc.cores = 30)


d_slope=tibble()
for (k in list.files("../Data/Moving_average_aridity",paste0("Moving"))){
  d=read.table(paste0("../Data/Moving_average_aridity/",k),sep=";")
  d_slope=rbind(d_slope,d)
}
write.table(d_slope,paste0("../Data/Moving_average_aridity/Moving_average_aridity.csv"),sep=";")



Run_model_moving_average_sand=function(id){
  
  list_mod=rbind(expand.grid(with_cover=c(T,F),
                             Stats=c("perim_area_scaling","fmax_psd","PL_expo","Spectral_ratio",
                                     "core_area_land","division","fractal_dim","contig","core_area",
                                     "flow_length","PLR","KS_dist","mean_psd",
                                     "Struct1","Struct2")), #We also add the case where the stats correspond to the first two components of the PCA
                 tibble(with_cover=F,Stats="rho_p"))
  
  d_slope=tibble()#for keeping all information
  
  d_data=read.table("../Data/Spatial_structure_grazing.csv",sep=";")%>%
    Closer_to_normality(.)
  
  d_data[,c(1,9:29,38:42,45,47:ncol(d_data))] = 
    apply(d_data[,c(1,9:29,38:42,45,47:ncol(d_data))],2,z_tranform)
  
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
  
  for (sand_thresh in seq(52.5,max(save$Sand),length.out=100)){ #keeping at least 90 sites
    
    d_data_mod=save%>%filter(., Sand<sand_thresh)
    d_data_mod$Grazing=scale(d_data_mod$Grazing)[,1]
    d_data_mod$Sand=scale(d_data_mod$Sand)[,1]
    
    model_spa_stat  = lm(formula_mod, d_data_mod,
                         na.action = na.fail )
    
    resid_mod=visreg::visreg(fit = model_spa_stat,xvar="Sand",plot=F) 
    
    mod_cov = boot(data=resid_mod$res, statistic=boot_function_lm, R=1000, formula=visregRes~Sand)
    
    #Merge in a df
    d_slope=rbind(d_slope,
                  tibble(pval=twoside_pvalue(mod_cov$t),
                         q2=median(mod_cov$t),
                         q1=quantile(mod_cov$t,.025),
                         q3=quantile(mod_cov$t,.975),
                         With_cover=with_cover,Sand_thresh=sand_thresh,
                         Stat=stat,Driver="Sand"
                  ))
  }
  
  write.table(d_slope,paste0("../Data/Moving_average_sand/Moving_average_sand_",stat,"_",with_cover,".csv"),sep=";")
  print("--")
}

library(parallel)

mclapply(1:30,Run_model_moving_average_sand,mc.cores = 30)


d_slope=tibble()
for (k in list.files("../Data/Moving_average_sand",paste0("Moving"))){
  d=read.table(paste0("../Data/Moving_average_sand/",k),sep=";")
  d_slope=rbind(d_slope,d)
}
write.table(d_slope,paste0("../Data/Moving_average_sand/Moving_average_sand.csv"),sep=";")






