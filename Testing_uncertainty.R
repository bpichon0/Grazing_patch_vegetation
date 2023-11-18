rm(list=ls())
source("./Structure_grazing_function.R")

Run_models=function(id_sim){
  
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
  
  nsim=50
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
    
    
    model_q=(lmer(formula = paste("q ~ ",mod_predictors), data = d2,na.action = na.fail ,REML ="FALSE"))
    model_rela=(lmer(formula = paste("rela_dist ~ ",mod_predictors), data = d2,na.action = na.fail ,REML ="FALSE"))
    model_size=(lmer(formula = paste("Size_tipping ~ ",mod_predictors), data = d2,na.action = na.fail ,REML ="FALSE"))
    
    #Getting partial prediction
    
    #Full grazing pressure
    resid_mod_size=visreg::visreg(fit = model_size,xvar="Grazing",plot=F) 
    mod_size=confint(lm(visregRes~Grazing,resid_mod_size$res))
    
    mod_size=ifelse(sign(mod_size[2,1] != mod_size[2,2]),"No_significant",
                   ifelse(mod_size[2,1]>0,"Positive","Negative"))
    
    
    resid_mod_rela=visreg::visreg(fit = model_rela,xvar="Grazing",plot=F) 
    mod_rela=confint(lm(visregRes~Grazing,resid_mod_rela$res))
    
    mod_rela=ifelse(sign(mod_rela[2,1] != mod_rela[2,2]),"No_significant",
                    ifelse(mod_rela[2,1]>0,"Positive","Negative"))
    
    resid_mod_q=visreg::visreg(fit = model_q,xvar="Grazing",plot=F) 
    mod_q=confint(lm(visregRes~Grazing,resid_mod_q$res))
    
    mod_q=ifelse(sign(mod_q[2,1] != mod_q[2,2]),"No_significant",
                    ifelse(mod_q[2,1]>0,"Positive","Negative"))
    
    d_partial=rbind(d_partial,
                    tibble(slope=c(mod_size,
                                   mod_rela,
                                   mod_q),
                           Stat=c("Distance to tipping","Size tipping","q (Spatial structure)"),
                           Type="Full"))
    
    
    model_q=confint(model_q)
    model_rela=confint(model_rela)
    model_size=confint(model_size)
    
    name_pred=rownames(model_q)[4:nrow(model_q)]
    signif_q=sapply(4:nrow(model_q),function(x){
      return(ifelse(sign(model_q[x,1] != model_q[x,2]),"No_significant",
                    ifelse(model_q[x,1]>0,"Positive","Negative")))
    })
    signif_rela=sapply(4:nrow(model_rela),function(x){
      return(ifelse(sign(model_rela[x,1] != model_rela[x,2]),"No_significant",
                    ifelse(model_rela[x,1]>0,"Positive","Negative")))
    })
    signif_size=sapply(4:nrow(model_size),function(x){
      return(ifelse(sign(model_size[x,1] != model_size[x,2]),"No_significant",
                    ifelse(model_size[x,1]>0,"Positive","Negative")))
    })
  
    
    d_mod=rbind(d_mod,tibble(Name=name_pred,
                             Signif=signif_q,
                             Param=c("q"),
                             Type_grazing="All"),
                tibble(Name=name_pred,
                       Signif=signif_rela,
                       Param=c("Relative distance"),
                       Type_grazing="All"),
                tibble(Name=name_pred,
                       Signif=signif_size,
                       Param=c("Size tipping"),
                       Type_grazing="All"))
    
    
    #High grazing pressure
    
    Sub=d2%>%filter(., Grazing %in% c(2,3))
    
    model_q=(lmer(formula = paste("q ~ ",mod_predictors), data = Sub,na.action = na.fail ,REML ="FALSE"))
    model_rela=(lmer(formula = paste("rela_dist ~ ",mod_predictors), data = Sub,na.action = na.fail ,REML ="FALSE"))
    model_size=(lmer(formula = paste("Size_tipping ~ ",mod_predictors), data = Sub,na.action = na.fail ,REML ="FALSE"))
    
    #Getting partial prediction
    
    #Full grazing pressure
    resid_mod_size=visreg::visreg(fit = model_size,xvar="Grazing",plot=F) 
    mod_size=confint(lm(visregRes~Grazing,resid_mod_size$res))
    
    mod_size=ifelse(sign(mod_size[2,1] != mod_size[2,2]),"No_significant",
                    ifelse(mod_size[2,1]>0,"Positive","Negative"))
    
    
    resid_mod_rela=visreg::visreg(fit = model_rela,xvar="Grazing",plot=F) 
    mod_rela=confint(lm(visregRes~Grazing,resid_mod_rela$res))
    
    mod_rela=ifelse(sign(mod_rela[2,1] != mod_rela[2,2]),"No_significant",
                    ifelse(mod_rela[2,1]>0,"Positive","Negative"))
    
    resid_mod_q=visreg::visreg(fit = model_q,xvar="Grazing",plot=F) 
    mod_q=confint(lm(visregRes~Grazing,resid_mod_q$res))
    
    mod_q=ifelse(sign(mod_q[2,1] != mod_q[2,2]),"No_significant",
                 ifelse(mod_q[2,1]>0,"Positive","Negative"))
    
    d_partial=rbind(d_partial,
                    tibble(slope=c(mod_size,
                                   mod_rela,
                                   mod_q),
                           Stat=c("Distance to tipping","Size tipping","q (Spatial structure)"),
                           Type="High"))
    
    model_q=confint(model_q)
    model_rela=confint(model_rela)
    model_size=confint(model_size)
    
    
    name_pred=rownames(model_q)[4:nrow(model_q)]
    signif_q=sapply(4:nrow(model_q),function(x){
      return(ifelse(sign(model_q[x,1] != model_q[x,2]),"No_significant",
                    ifelse(model_q[x,1]>0,"Positive","Negative")))
    })
    signif_rela=sapply(4:nrow(model_rela),function(x){
      return(ifelse(sign(model_rela[x,1] != model_rela[x,2]),"No_significant",
                    ifelse(model_rela[x,1]>0,"Positive","Negative")))
    })
    signif_size=sapply(4:nrow(model_size),function(x){
      return(ifelse(sign(model_size[x,1] != model_size[x,2]),"No_significant",
                    ifelse(model_size[x,1]>0,"Positive","Negative")))
    })
    
    
    d_mod=rbind(d_mod,tibble(Name=name_pred,
                             Signif=signif_q,
                             Param=c("q"),
                             Type_grazing="High"),
                tibble(Name=name_pred,
                       Signif=signif_rela,
                       Param=c("Relative distance"),
                       Type_grazing="High"),
                tibble(Name=name_pred,
                       Signif=signif_size,
                       Param=c("Size tipping"),
                       Type_grazing="High"))
    
    
    #Low grazing pressure
    
    Sub=d2%>%filter(., Grazing %in% c(0,1))
    
    model_q=(lmer(formula = paste("q ~ ",mod_predictors), data = Sub,na.action = na.fail ,REML ="FALSE"))
    model_rela=(lmer(formula = paste("rela_dist ~ ",mod_predictors), data = Sub,na.action = na.fail ,REML ="FALSE"))
    model_size=(lmer(formula = paste("Size_tipping ~ ",mod_predictors), data = Sub,na.action = na.fail ,REML ="FALSE"))
    
    #Getting partial prediction
    
    #Full grazing pressure
    resid_mod_size=visreg::visreg(fit = model_size,xvar="Grazing",plot=F) 
    mod_size=confint(lm(visregRes~Grazing,resid_mod_size$res))
    
    mod_size=ifelse(sign(mod_size[2,1] != mod_size[2,2]),"No_significant",
                    ifelse(mod_size[2,1]>0,"Positive","Negative"))
    
    
    resid_mod_rela=visreg::visreg(fit = model_rela,xvar="Grazing",plot=F) 
    mod_rela=confint(lm(visregRes~Grazing,resid_mod_rela$res))
    
    mod_rela=ifelse(sign(mod_rela[2,1] != mod_rela[2,2]),"No_significant",
                    ifelse(mod_rela[2,1]>0,"Positive","Negative"))
    
    resid_mod_q=visreg::visreg(fit = model_q,xvar="Grazing",plot=F) 
    mod_q=confint(lm(visregRes~Grazing,resid_mod_q$res))
    
    mod_q=ifelse(sign(mod_q[2,1] != mod_q[2,2]),"No_significant",
                 ifelse(mod_q[2,1]>0,"Positive","Negative"))
    
    d_partial=rbind(d_partial,
                    tibble(slope=c(mod_size,
                                   mod_rela,
                                   mod_q),
                           Stat=c("Distance to tipping","Size tipping","q (Spatial structure)"),
                           Type="Low"))
    
    
    model_q=confint(model_q)
    model_rela=confint(model_rela)
    model_size=confint(model_size)
    
    name_pred=rownames(model_q)[4:nrow(model_q)]
    signif_q=sapply(4:nrow(model_q),function(x){
      return(ifelse(sign(model_q[x,1] != model_q[x,2]),"No_significant",
                    ifelse(model_q[x,1]>0,"Positive","Negative")))
    })
    signif_rela=sapply(4:nrow(model_rela),function(x){
      return(ifelse(sign(model_rela[x,1] != model_rela[x,2]),"No_significant",
                    ifelse(model_rela[x,1]>0,"Positive","Negative")))
    })
    signif_size=sapply(4:nrow(model_size),function(x){
      return(ifelse(sign(model_size[x,1] != model_size[x,2]),"No_significant",
                    ifelse(model_size[x,1]>0,"Positive","Negative")))
    })
    
    
    d_mod=rbind(d_mod,tibble(Name=name_pred,
                             Signif=signif_q,
                             Param=c("q"),
                             Type_grazing="Low"),
                tibble(Name=name_pred,
                       Signif=signif_rela,
                       Param=c("Relative distance"),
                       Type_grazing="Low"),
                tibble(Name=name_pred,
                       Signif=signif_size,
                       Param=c("Size tipping"),
                       Type_grazing="Low"))
    
    print(k)
    
  }
  
  write.table(d_mod,paste0("../Data/Inferrence/Test/Drivers_stability_metrics_no_cover_test_",id_sim,".csv"),sep=";")
  write.table(d_partial,paste0("../Data/Inferrence/Test/Partial_residuals_grazing_no_cover_test_",id_sim,".csv"),sep=";")
  
}
  

library(parallel)
mclapply(1:20,Run_models,mc.cores = 20)

