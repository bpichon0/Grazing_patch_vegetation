rm(list=ls())
source("./Structure_grazing_function.R")

# ---------------------- Step 1: Effects of grazing on the spatial structure  ----
## >> 1) Running mixed-effect models without controlling for veg. cover ---- 

Run_model_importance=function(id){
  
  dir.create("./Data/Linear_models_factor/Keep_data",showWarnings = F)
  dir.create("./Data/Linear_models_factor/Keep_models",showWarnings = F)
  
  
  list_mod=rbind(expand.grid(with_cover=c(F),
                             Stats=c("perim_area_scaling","fmax_psd","PL_expo",
                                     "core_area_land","core_area","Small_patches",
                                     "flow_length","mean_psd","moran_I")))
  
  
  d_data=read.table("./Data/Spatial_structure_grazing.csv",sep=";")%>%
    Closer_to_normality(.)
  
  
  
  d_data=Perform_PCA_spatial_struc(d_data)
  
  stat=list_mod$Stats[id]
  
  d_data_mod=d_data%>%melt(., measure.vars=stat)%>%
    filter(., !is.na(value),!is.na(Herbivores))
  
  with_cover=list_mod$with_cover[id]
  
  if (with_cover){
    
    formula_mod=formula(formula_(paste("value ~ Long_cos + Long_sin + Lattitude + Slope + Elevation 
      + Aridity + Grazing  + Herbivores
      + Aridity*Grazing + MAT
      + rho_p + Type_veg 
      + Sand + Org_C_v 
      + (1|Site_ID)")))
    
  } else{
    
    formula_mod=formula(formula_(paste("value ~ Long_cos + Long_sin + Lattitude + Slope + Elevation 
      + Aridity + Grazing  + Herbivores
      + Aridity*Grazing  + MAT
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
    
    d_data[,colnames(dplyr::select_if(d_data, is.numeric))]=
      apply(d_data[,colnames(dplyr::select_if(d_data, is.numeric))],2,
            function(x){return((x-mean(x,na.rm=T))/sd(x,na.rm = T))})
    
    model_spa_stat  = lmer(formula_mod, d_data_mod,
                           na.action = na.fail ,REML ="FALSE")
    
    #we remove potential outliers
    rm.outliers = romr.fnc(model_spa_stat, d_data_mod, trim=2.5)
    d_data_out = rm.outliers$data
    
    #saving data
    write.table(d_data_out,paste0("./Data/Linear_models_factor/Keep_data/Data_",
                                  stat,"_",with_cover,"_aridity_",grazing_intensity,".csv"))
    
    model_spa_stat = lmer(formula_mod, data = d_data_out, 
                          na.action = na.fail,REML ="FALSE")

    
    saveRDS(model_spa_stat,#MuMIn::get.models(select_model,subset = delta<2)[[1]],
            paste0("./Data/Linear_models_factor/Keep_models/Mod_",stat,"_",
                   with_cover,"_aridity_",grazing_intensity,".rds"))
    
    
  }  
  
  
}


lapply(1:9,Run_model_importance)


## >> 2) Analyzing the residuals ---- 

d_data=read.table("./Data/Spatial_structure_grazing.csv",sep=";")%>%
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
with_interaction=T
d_slope=tibble()

for (grazing_intensity in c("binary","all")){
  
  for (k in c("perim_area_scaling","fmax_psd","PL_expo",
              "core_area_land","core_area","Small_patches",
              "flow_length","mean_psd","moran_I","rho_p")){
    
    #for each we plot the slope against the partial residuals
    
    d_data_out=read.table(paste0("./Data/Linear_models_factor/Keep_data/Data_",k,"_FALSE_aridity_",grazing_intensity,".csv"),sep=" ")
    model_spa_stat=readRDS(paste0("./Data/Linear_models_factor/Keep_models/Mod_",k,"_FALSE_aridity_",grazing_intensity,".rds"))
    
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
    
    resid_mod=visreg::visreg(fit = model_spa_stat,xvar="MAT",plot=F) 
    mod_cov = boot(data=resid_mod$res, statistic=boot_function_lm_arid, R=1000, formula=visregRes~MAT)
    
    d_slope=rbind(d_slope,
                  tibble(pval=twoside_pvalue(mod_cov$t),
                         q2=median(mod_cov$t),
                         q1=quantile(mod_cov$t,.025),
                         q3=quantile(mod_cov$t,.975),
                         q1_90=quantile(mod_cov$t,.05),
                         q3_90=quantile(mod_cov$t,.95),
                         With_cover=F,Grazing_intensity=grazing_intensity,
                         ID_grazing=c("all"),
                         Stat=k,Driver="MAT"
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

write.table(d_slope,paste0("./Data/Linear_models_factor/Slope_partial_residuals_aridity.csv"),sep=";")


# ---------------------- Step 2: Effects of grazing on the spatial structure cover control residuals ----
## >> 1) Running mixed-effect model ----

Run_model_importance=function(id){
  
  dir.create("./Data/Linear_models_factor_cover_control/Keep_data",showWarnings = F)
  dir.create("./Data/Linear_models_factor_cover_control/Keep_models",showWarnings = F)
  
  list_mod=rbind(expand.grid(with_cover=c(F),
                             Stats=c("perim_area_scaling","fmax_psd","PL_expo",
                                     "core_area_land","core_area","Small_patches",
                                     "flow_length","mean_psd","moran_I","Struct1","Struct2")))
  
  d_all=d_all2=tibble()#for keeping all information
  
  d_data=read.table("./Data/Spatial_structure_grazing.csv",sep=";")%>%
    Closer_to_normality(.)%>%Perform_PCA_spatial_struc(.)
  
  d_data=Perform_PCA_spatial_struc(d_data)
  
  stat=list_mod$Stats[id]
  
  d_data_mod=d_data%>%melt(., measure.vars=stat)%>%
    dplyr::filter(., !is.na(value),!is.na(Herbivores))
  
  with_cover=list_mod$with_cover[id]
  
  if (with_cover){
    
    formula_mod=formula(formula_(paste("value ~ Long_cos + Long_sin + Lattitude + Slope + Elevation 
      + Aridity + Grazing + Herbivores
      + Aridity*Grazing + MAT
      + rho_p + Type_veg 
      + Sand + Org_C_v
      + (1|Site_ID)")))
    
  } else{
    
    formula_mod=formula(formula_(paste("value ~ Long_cos + Long_sin + Lattitude + Slope + Elevation 
      + Aridity + Grazing + Herbivores
      + Aridity*Grazing  + MAT
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
    
    d_data[,colnames(dplyr::select_if(d_data, is.numeric))]=
      apply(d_data[,colnames(dplyr::select_if(d_data, is.numeric))],2,
            function(x){return((x-mean(x,na.rm=T))/sd(x,na.rm = T))})
    
    model_spa_stat  = lmer(formula_mod, d_data_mod,
                           na.action = na.fail ,REML ="FALSE")
    
    #we remove potential outliers
    rm.outliers = romr.fnc(model_spa_stat, d_data_mod, trim=2.5)
    d_data_out = rm.outliers$data
    
    #saving data
    write.table(d_data_out,paste0("./Data/Linear_models_factor_cover_control/Keep_data/Data_",
                                  stat,"_",with_cover,"_aridity_",grazing_intensity,".csv"))
    
    model_spa_stat = lmer(formula_mod, data = d_data_out, 
                          na.action = na.fail,REML ="FALSE")
    
    
    saveRDS(model_spa_stat,
            paste0("./Data/Linear_models_factor_cover_control/Keep_models/Mod_",stat,"_",
                   with_cover,"_aridity_",grazing_intensity,".rds"))
    
    
  }  
  
}

lapply(1:11,Run_model_importance)


## >> 2) Analyzing the residuals ---- 

d_data=read.table("./Data/Spatial_structure_grazing.csv",sep=";")%>%
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

d_slope=tibble()

for (grazing_intensity in c("binary","all")){
  
  for (k in c("perim_area_scaling","fmax_psd","PL_expo",
              "core_area_land","core_area","Small_patches",
              "flow_length","mean_psd","moran_I","Struct1","Struct2")){
    
    #for each we plot the slope against the partial residuals with and without cover
    
    d_data_out=read.table(paste0("./Data/Linear_models_factor_cover_control/Keep_data/Data_",k,"_FALSE_aridity_",grazing_intensity,".csv"),sep=" ")
    model_spa_stat=readRDS(paste0("./Data/Linear_models_factor_cover_control/Keep_models/Mod_",k,"_FALSE_aridity_",grazing_intensity,".rds"))
    
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
    
    
    resid_mod=visreg::visreg(fit = model_spa_stat,xvar="MAT",plot=F) 
    
    
    mod_cov = boot(data=resid_mod$res, statistic=boot_function_lm_arid, R=1000, formula=visregRes~MAT)
    
    d_slope=rbind(d_slope,
                  tibble(pval=twoside_pvalue(mod_cov$t),
                         q2=median(mod_cov$t),
                         q1=quantile(mod_cov$t,.025),
                         q3=quantile(mod_cov$t,.975),
                         q1_90=quantile(mod_cov$t,.05),
                         q3_90=quantile(mod_cov$t,.95),
                         With_cover=F,Grazing_intensity=grazing_intensity,
                         ID_grazing=c("all"),
                         Stat=k,Driver="MAT"
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

write.table(d_slope,paste0("./Data/Linear_models_factor_cover_control/Slope_partial_residuals_aridity.csv"),sep=";")





# ---------------------- Step 3: Path analysis with traits ----

#Data of spatial structure

d_data=Perform_PCA_spatial_struc(
  Closer_to_normality_CWM(
    read.table("./Data/Spatial_structure_control_cover_with_traits.csv",sep=";")))%>%
  mutate(., Grazing=as.factor(Grazing),ID=as.character(ID))%>%
  mutate(., Dev_MaxLS=-Dev_MaxLS,#positive values = facilitation
         Grazing=scale(as.numeric(Grazing))[,1])

d_data[,colnames(dplyr::select_if(d_data, is.numeric))]=apply(d_data[,colnames(dplyr::select_if(d_data, is.numeric))],
                                                    2,
                                                    function(x){return((x-mean(x,na.rm=T))/sd(x,na.rm = T))})

d=read.table("./Data/Data_SEM.csv",sep=";")
d[,colnames(dplyr::select_if(d, is.numeric))]=
  apply(d[,colnames(dplyr::select_if(d, is.numeric))],2,
        function(x){return((x-mean(x,na.rm=T))/sd(x,na.rm = T))})

d_mod1=dplyr::filter(d,!is.na(CWM_MaxLS))%>%
  mutate(., Grazing=as.numeric(Grazing))

mod1=lmer(CWM_MaxLS~Aridity+I(Aridity^2)+(1|Site_ID)+Grazing+Woody+MAT+I(MAT^2)+
            Lattitude+Long_cos+Long_sin+Elevation,data=d_mod1,na.action = na.fail)

d_mod3=dplyr::filter(d,!is.na(Woody))%>%
  mutate(., Grazing=as.numeric(Grazing))

mod3=lmer(Woody~Aridity+I(Aridity^2)+Grazing+(1|Site_ID)+MAT+I(MAT^2)+
            Lattitude+Long_cos+Long_sin+Elevation,data=d_mod3,na.action = na.fail)


#Since woody and CWM correlated, we perform PLS
#First for DEV_MAXLS

d_mod2=d%>%
  mutate(., Grazing=as.numeric(Grazing))%>%
  add_column(., MAT2=.$MAT^2)

d_mod2=d_mod2[,c("Dev_MaxLS","CWM_MaxLS","Aridity","Woody","Grazing","Lattitude","Long_cos","Long_sin","Elevation","MAT","MAT2")]%>%
  dplyr::filter(., !is.na(Dev_MaxLS),!is.na(CWM_MaxLS))%>%
  dplyr::mutate(., Dev_MaxLS= -Dev_MaxLS) #positive = facilitation

cv.modpls = cv.plsR(Dev_MaxLS~.,data=d_mod2,K=10,nt=10,modele = "pls",
                    grouplist = createFolds(d_mod2[,1], k = 10, list = F, returnTrain = FALSE),NK=200,verbose = F)

#We check number of components
res_cv=cvtable(summary(cv.modpls,MClassed=T))
#5 or 3 components according to CV

weights_mod=res_cv$CVPress[c(4,5)] #117 and 47

mod_pls = plsR(Dev_MaxLS~.,data=d_mod2,nt=10,pvals.expli=TRUE,typeVC = "adaptative")

#2 or 3 components according to AIC 


#We perform the two PLS
mod_pls = plsR(Dev_MaxLS~.,data=d_mod2,nt=4,pvals.expli=TRUE,typeVC = "adaptative")
boot_pls_mod_4= bootpls(mod_pls,typeboot="fmodel_np",R=2000,verbose = T)
save_mod_pls=boot_pls_mod_4

mod_pls = plsR(Dev_MaxLS~.,data=d_mod2,nt=5,pvals.expli=TRUE,typeVC = "adaptative")
boot_pls_mod_5= bootpls(mod_pls,typeboot="fmodel_np",R=2000,verbose = T)
#R2 = 0.39

#And coefficient by the relative support of the two PLS
d_result_mod2 = data.frame(predictor=c("CWM_MaxLS","Aridity","Woody","Grazing","Lattitude","Long_cos","Long_sin","Elevation","MAT","MAT2"),
                          q2=cbind(boot_pls_mod_4$t0[-1,],boot_pls_mod_5$t0[-1,]) %*% prop.table(weights_mod), #weighted mean
                          mean=cbind(apply(boot_pls_mod_4$t[,-1],2,mean),
                                     apply(boot_pls_mod_5$t[,-1],2,mean))%*% prop.table(weights_mod),
                          q1=apply(cbind(apply(boot_pls_mod_4$t[,-1],2,quantile,.025),
                                         apply(boot_pls_mod_5$t[,-1],2,quantile,.025)),1,
                                   function(x){return(min((x)))}),#return largest incertitude
                          q3=apply(cbind(apply(boot_pls_mod_4$t[,-1],2,quantile,.975),
                                         apply(boot_pls_mod_5$t[,-1],2,quantile,.975)),1,
                                   function(x){return(max((x)))}),#return largest incertitude
                          sd=apply(cbind(apply(boot_pls_mod_4$t[,-1],2,sd),
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
         apply(save_mod_pls$t[,-1],2,twoside_pvalue),#apply(boot_pls_mod_8$t[,-1],2,twoside_pvalue) give the same
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


write.table(d_all,"./Data/CI_links_dsep.csv",sep=";")

# ---------------------- Step 4: Using the model with different spatial structures ----
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

saveRDS(object = list_landscape,"./Data/Minimal_examples_stats_landscapes.rds")

write.table(d,"./Data/Minimal_examples_stats.csv",sep=";")




## >> 2) Effect size parameter ----

parameters=matrix(c(0.001,.3,0.05,.14,0.1,0,.001,.319,0.05,.202,0.17,0,0.001,.325,0.05,.20802,0.18,0),9,2,byrow = T)
d=read.table("./Data/Minimal_examples_stats.csv",sep=";")%>%
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

write.table(d_effect,"./Data/Model_coefficients.csv",sep=";")




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
write.table(d_effect,"./Data/Model_coefficients.csv",sep=";")




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
    
  }
  
}

# ---------------------- Step 5: Sensitivity spatial resolution ----
## >> 1) Running mixed-effect model ----

Run_model_importance_aridity=function(id){
  
  dir.create("./Data/Linear_models_resolution/Keep_data",showWarnings = F)
  dir.create("./Data/Linear_models_resolution/Keep_models",showWarnings = F)
  
  list_mod=rbind(expand.grid(with_cover=c(F),
                             Stats=c("perim_area_scaling","fmax_psd","PL_expo",
                                     "core_area_land","core_area","Small_patches",
                                     "flow_length","mean_psd","moran_I")),
                 tibble(with_cover=F,
                        Stats=c("rho_p")))
  
  d_all=d_all2=tibble()#for keeping all information
  
  d_data=read.table("./Data/Spatial_structure_grazing.csv",sep=";")%>%
    Closer_to_normality(.)
  
  
  d_data=Perform_PCA_spatial_struc(d_data)
  
  stat=list_mod$Stats[id]
  
  d_data_mod=d_data%>%melt(., measure.vars=stat)%>%
    filter(., !is.na(value),!is.na(Herbivores))
  
  with_cover=list_mod$with_cover[id]
  
  if (with_cover){
    
    formula_mod=formula(formula_(paste("value ~ Long_cos + Long_sin + Lattitude + Slope + Elevation 
      + Aridity + Grazing + Resolution + Herbivores + MAT
      + Aridity*Grazing 
      + rho_p + Type_veg 
      + Sand + Org_C_v 
      + (1|Site_ID)")))
    
  } else{
    
    formula_mod=formula(formula_(paste("value ~ Long_cos + Long_sin + Lattitude + Slope + Elevation 
      + Aridity + Grazing + Resolution + Herbivores + MAT
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
    
    d_data[,colnames(dplyr::select_if(d_data, is.numeric))]=
      apply(d_data[,colnames(dplyr::select_if(d_data, is.numeric))],2,
            function(x){return((x-mean(x,na.rm=T))/sd(x,na.rm = T))})
    
    model_spa_stat  = lmer(formula_mod, d_data_mod,
                           na.action = na.fail ,REML ="FALSE")
    
    #we remove potential outliers
    rm.outliers = romr.fnc(model_spa_stat, d_data_mod, trim=2.5)
    d_data_out = rm.outliers$data
    
    #saving data
    write.table(d_data_out,paste0("./Data/Linear_models_resolution/Keep_data/Data_",
                                  stat,"_",with_cover,"_aridity_",grazing_intensity,".csv"))
    
    model_spa_stat = lmer(formula_mod, data = d_data_out, 
                          na.action = na.fail,REML ="FALSE")
    saveRDS(model_spa_stat,#MuMIn::get.models(select_model,subset = delta<2)[[1]],
            paste0("./Data/Linear_models_resolution/Keep_models/Mod_",stat,"_",
                   with_cover,"_aridity_",grazing_intensity,".rds"))
    
    
  }  
  
}

library(parallel)

mclapply(1:10,Run_model_importance_aridity,mc.cores = 10)


## >> 2) Analyzing the residuals ---- 

d_data=read.table("./Data/Spatial_structure_grazing.csv",sep=";")%>%
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

for (with_interaction in c(T)){
  d_slope=tibble()
  
  for (grazing_intensity in c("binary","all")){
    
    for (k in c("perim_area_scaling","fmax_psd","PL_expo",
                "core_area_land","core_area","Small_patches",
                "flow_length","mean_psd","moran_I")){
      
      #for each we plot the slope against the partial residuals with and without cover
      
      d_data_out=read.table(paste0("./Data/Linear_models_resolution/Keep_data/Data_",k,"_FALSE_aridity_",grazing_intensity,".csv"),sep=" ")
      model_spa_stat=readRDS(paste0("./Data/Linear_models_resolution/Keep_models/Mod_",k,"_FALSE_aridity_",grazing_intensity,".rds"))
      
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
  
  write.table(d_slope,paste0("./Data/Linear_models_resolution/Slope_partial_residuals_aridity.csv"),sep=";")
}







