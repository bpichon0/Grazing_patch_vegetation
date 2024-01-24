rm(list=ls())
source("./Structure_grazing_function.R")

#State at 16/01: does not work well. Abort

list_null_metrics=list.files("../Data/Metrics_null",".csv")
d_data=read.table("../Data/Spatial_structure_grazing.csv",sep=";")

d_all_null=tibble()
for (file_id in list_null_metrics){
  
  name_file=gsub("Metric_","",gsub("_null.csv","",file_id))
  file_null=read.table(paste0("../Data/Metrics_null/",file_id),sep=";")
  d_data[which(d_data$Full_name==name_file),2:30]=(d_data[which(d_data$Full_name==name_file),2:30]-
                                                     apply(file_null[,-ncol(file_null)],2,mean,na.rm=T))/
    apply(file_null[,-ncol(file_null)],2,sd,na.rm=T)
}


write.table(d_data,"../Data/Spatial_structure_grazing_null_expectation.csv",sep=";")


rm(list=ls())
source("./Structure_grazing_function.R")

# ---------------------- Step 2: Running mixed-effect models  ----

Run_model_importance_aridity=function(id){
  
  list_mod=rbind(expand.grid(with_cover=c(T,F),
                             Stats=c("perim_area_scaling","fmax_psd","PL_expo","Spectral_ratio",
                                     "core_area_land","division","fractal_dim","contig","core_area",
                                     "flow_length","PLR","KS_dist","Shape_metric",
                                     "Struct1","Struct2")), #We also add the case where the stats correspond to the first two components of the PCA
                 tibble(with_cover=F,Stats="rho_p"))
  
  d_all=d_all2=tibble()#for keeping all information
  
  d_data=read.table("../Data/Spatial_structure_grazing_null_expectation.csv",sep=";")%>%
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
    write.table(d_data_out,paste0("../Data/Linear_models_null/Keep_data/Data_",stat,"_",with_cover,"_aridity_",grazing_intensity,".csv"))
    
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
    
    saveRDS(model_spa_stat,paste0("../Data/Linear_models_null/Keep_models/Mod_",stat,"_",with_cover,"_aridity_",grazing_intensity,".rds"))
    
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
  write.table(d_all,paste0("../Data/Linear_models_null/Importance/Importance_",stat,"_",with_cover,"_aridity.csv"),sep=";")
  write.table(d_all2,paste0("../Data/Linear_models_null/Estimator/Estimators_model_",stat,"_",with_cover,"_aridity.csv"),sep=";")
  
}

library(parallel)

mclapply(1:30,Run_model_importance_aridity,mc.cores = 30)

d_all=d_all2=tibble()
for (k in list.files("../Data/Linear_models_null/Importance",paste0("aridity.csv"))){
  d=read.table(paste0("../Data/Linear_models_null/Importance/",k),sep=";")
  if (nrow(d)!=0){
    d_all=rbind(d_all,d)
  }
}
for (k in list.files("../Data/Linear_models_null/Estimator",paste0("aridity.csv"))){
  d2=read.table(paste0("../Data/Linear_models_null/Estimator/",k),sep=";")
  if (nrow(d2)!=0){
    d_all2=rbind(d_all2,d2)
  }
}
write.table(d_all,paste0("../Data/Linear_models_null/Importance_","aridity.csv"),sep=";")
write.table(d_all2,paste0("../Data/Linear_models_null/Estimators_model_","aridity.csv"),sep=";")


Run_model_importance_aridity_no_inter=function(id){
  
  list_mod=rbind(expand.grid(with_cover=c(T,F),
                             Stats=c("perim_area_scaling","fmax_psd","PL_expo","Spectral_ratio",
                                     "core_area_land","division","fractal_dim","contig","core_area",
                                     "flow_length","PLR","KS_dist","Shape_metric",
                                     "Struct1","Struct2")), #We also add the case where the stats correspond to the first two components of the PCA
                 tibble(with_cover=F,Stats="rho_p"))
  
  d_all=d_all2=tibble()#for keeping all information
  
  d_data=read.table("../Data/Spatial_structure_grazing_null_expectation.csv",sep=";")%>%
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
    write.table(d_data_out,paste0("../Data/Linear_models_null/Keep_data/Data_",stat,"_",with_cover,"_aridity_no_inter_",grazing_intensity,".csv"))
    
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
    
    saveRDS(model_spa_stat,paste0("../Data/Linear_models_null/Keep_models/Mod_",stat,"_",with_cover,"_aridity_no_inter_",grazing_intensity,".rds"))
    
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
  write.table(d_all,paste0("../Data/Linear_models_null/Importance/Importance_",stat,"_",with_cover,"_aridity_no_inter.csv"),sep=";")
  write.table(d_all2,paste0("../Data/Linear_models_null/Estimator/Estimators_model_",stat,"_",with_cover,"_aridity_no_inter.csv"),sep=";")
  
}

library(parallel)

mclapply(1:30,Run_model_importance_aridity_no_inter,mc.cores = 30)

d_all=d_all2=tibble()
for (k in list.files("../Data/Linear_models_null/Importance",paste0("aridity_no_inter"))){
  d=read.table(paste0("../Data/Linear_models_null/Importance/",k),sep=";")
  if (nrow(d)!=0){
    d_all=rbind(d_all,d)
  }
}
for (k in list.files("../Data/Linear_models_null/Estimator",paste0("aridity_no_inter"))){
  d2=read.table(paste0("../Data/Linear_models_null/Estimator/",k),sep=";")
  if (nrow(d)!=0){
    d_all2=rbind(d_all2,d2)
  }
}
write.table(d_all,paste0("../Data/Linear_models_null/Importance_","aridity_no_inter.csv"),sep=";")
write.table(d_all2,paste0("../Data/Linear_models_null/Estimators_model_","aridity_no_inter.csv"),sep=";")



# ---------------------- Step 3: Analyzing the residuals  ----
## >> 2) Model with aridity and grazing ----

d_data=read.table("../Data/Spatial_structure_grazing_null_expectation.csv",sep=";")%>%
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
      
      d_data_out=read.table(paste0("../Data/Linear_models_null/Keep_data/Data_",k,"_TRUE_aridity_",grazing_intensity,".csv"),sep=" ")
      model_spa_stat=readRDS(paste0("../Data/Linear_models_null/Keep_models/Mod_",k,"_TRUE_aridity_",grazing_intensity,".rds"))
      
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
    
    d_data_out=read.table(paste0("../Data/Linear_models_null/Keep_data/Data_",k,"_FALSE_aridity_",grazing_intensity,".csv"),sep=" ")
    model_spa_stat=readRDS(paste0("../Data/Linear_models_null/Keep_models/Mod_",k,"_FALSE_aridity_",grazing_intensity,".rds"))
    
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

write.table(d_slope,"../Data/Linear_models_null/Slope_partial_residuals_aridity.csv",sep=";")

## >> 3) Model with aridity and grazing but without interactions ----

d_data=read.table("../Data/Spatial_structure_grazing_null_expectation.csv",sep=";")%>%
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
      
      d_data_out=read.table(paste0("../Data/Linear_models_null/Keep_data/Data_",k,
                                   "_TRUE_aridity_no_inter_",grazing_intensity,".csv"),sep=" ")
      model_spa_stat=readRDS(paste0("../Data/Linear_models_null/Keep_models/Mod_",k,
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
    d_data_out=read.table(paste0("../Data/Linear_models_null/Keep_data/Data_",k,
                                 "_FALSE_aridity_no_inter_",grazing_intensity,".csv"),sep=" ")
    model_spa_stat=readRDS(paste0("../Data/Linear_models_null/Keep_models/Mod_",k,
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

write.table(d_slope,"../Data/Linear_models_null/Slope_partial_residuals_aridity_no_inter.csv",sep=";")


