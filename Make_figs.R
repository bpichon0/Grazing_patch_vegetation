rm(list=ls())
source("./Structure_grazing_function.R")



dir.create("../Figures/Final_figs",showWarnings = F)
dir.create("../Figures/Final_figs/SI",showWarnings = F)
# -------------------- Main figures -------------------------

## >> Figure 1: Importance and effect of grazing no cover ----

d_slope=read.table("../Data/Step1_Understanding_grazing/Slope_partial_residuals.csv",sep=";")%>%
  filter(., With_cover==1)%>%
  Filter_relevant_stats(.)%>%  
  Rename_spatial_statistics(.)%>%
  add_column(., Signif=sapply(1:nrow(.),function(x){
    return(is_signif(.$pval[x]))}))%>%
  arrange(., slope)%>%
  add_column(., Order_f=1:nrow(.))%>%
  mutate(.,Stat = fct_reorder(Stat, Order_f))%>%
  add_column(., order_signif=sapply(1:nrow(.),
                                    function(x){return(ifelse(.$slope[x]<0,.075,-.12))}))

p2=ggplot(d_slope)+
  geom_pointrange(aes(x=slope,y=Stat,xmin=Low_int,xmax=High_int,color=pval<.05),shape=17,size=.8)+
  the_theme+
  scale_color_manual(values=c("grey","black"))+
  #geom_text(aes(x=order_signif,y=1:8,label=Signif))+
  geom_vline(xintercept = 0,linetype=9)+
  labs(x=substitute(paste(beta," (partial res. Grazing)")),y="",color="")+
  theme(legend.position = "none")



#Importance based on R2
d_R2=read.table(paste0("../Data/Step1_Understanding_grazing/Importance.csv"),sep=";")%>%
  filter(., With_cover==1)%>%
  dplyr::select(., Sp_stat,R2m,R2C)%>%
  dplyr::rename(., Stat=Sp_stat)%>%
  Filter_relevant_stats(.)%>%  
  Rename_spatial_statistics(.)%>%
  add_column(., Order_f=sapply(1:nrow(.),function(x){
    which(d_slope$Stat==.$Stat[x])
  }))%>%
  arrange(., Order_f)%>%
  mutate(.,Stat = fct_reorder(Stat, Order_f))
  

d_all2=read.table(paste0("../Data/Step1_Understanding_grazing/Estimators_model.csv"),sep=";")%>%
  dplyr::rename(., observed=Median)%>%
  filter(., With_cover==1)%>%
  Filter_relevant_stats(.)%>%  
  Organize_df(., "bar")%>%
  group_by(., Type_pred,Stat)%>%
  dplyr::summarise(., sum_effect=sum(abs(observed)),.groups = "keep")%>%
  Rename_spatial_statistics(.)%>%group_by(., Stat)%>%
  dplyr::summarise(., .groups="keep",sum_effect=sum_effect/sum(sum_effect),Type_pred=Type_pred)%>%
  add_column(., Order_f=sapply(1:nrow(.),function(x){
    which(d_slope$Stat==.$Stat[x])
  }))%>%
  arrange(., Order_f)%>%
  mutate(.,Stat = fct_reorder(Stat, Order_f))
  


p1=ggplot(d_all2)+
  geom_bar(aes(x=100*sum_effect,y=Stat,fill=Type_pred),stat="identity",width = .75)+
  geom_text(data=d_R2,aes(y=Stat,x=20,label=paste0("R²m = ",round(R2m,2),", ","R²c = ",round(R2C,2))),size=3)+
  the_theme+
  labs(x="Importance (%)",fill="")+
  scale_fill_manual(values=c("Geography"="#FFB15B","Grazing"="#F3412A",
                             "Abiotic"="#4564CE","Vegetation"="#6DD275","Climatic"="#CC80C7"))+
  theme(axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank())

p_tot=ggarrange(
  ggarrange(p1+theme(legend.position = "none"),
            p2+theme(axis.text.y = element_blank(),axis.ticks.y = element_blank()),
            ncol=2,labels = letters[1:2],widths = c(1,.5),align = "h"),
  ggarrange(ggplot()+theme_void(),get_legend(p1),ggplot()+theme_void(),ncol=3,widths = c(.2,1,.6)),
  nrow=2,heights = c(1,.1))


ggsave("../Figures/Final_figs/Importance_effect_grazing_partial_res.pdf",
       p_tot,width = 8,height = 4)



## >> Figure 2: Interactive effect of grazing ----

d_all=read.table(paste0("../Data/Step1_Understanding_grazing/Importance.csv"),sep=";")%>%
  add_column(., With_cover=rep(c(T,F),each=nrow(.)/2))%>%
  filter(., With_cover==T)

d_all2=read.table(paste0("../Data/Step1_Understanding_grazing/Estimators_model.csv"),sep=";")%>%
  dplyr::rename(., observed=Median)


cover=T;index=1
for (k in c("perim_area_scaling","core_area")){
  
  title=ifelse(k=="perim_area_scaling","Fractal scaling area-perimeter","% of core vegetation")
  
  loc_pval=.05+max(d_all2%>% #localization of pvalues in the plot
                     filter(., Stat==k , With_cover==cover,term!="Type, both")%>%#not enougth data
                     Organize_df(., "predictor")%>%
                     dplyr::select(.,q3)%>%dplyr::pull(.))
  
  assign(paste0("p1_",index),ggplot(d_all2%>%
                filter(., Stat==k , With_cover==cover,term!="Type, both")%>%
                Organize_df(., "predictor")%>%
                add_column(., Is_signif_pval=sapply(1:nrow(.),function(x){return(is_signif(.$pvalue[x]))})))+
    geom_pointrange(aes(x=observed,y=term,xmin=q1,xmax=q3,color=Type_pred))+
    geom_text(aes(x=loc_pval,y=term,label=Is_signif_pval))+
    the_theme+  
    ggtitle(title)+
    labs(x="",y="",color="")+
    geom_vline(xintercept = 0,linetype=9)+
    scale_color_manual(values=c("Geography"="#FFB15B","Grazing"="#F3412A",
                                "Abiotic"="#4564CE","Vegetation"="#6DD275","Climatic"="#CC80C7"))+
    theme(legend.position = "none",title = element_text(size=10)))

  index=index+1 
}


list_resid=list(stat=c("Struct2", "Struct2",
                       "Struct1"),
                metric=c("Org_C","Clim2",
                         "Clim1"),
                title=c("Partial res. PC2 struct.","Partial res. PC2 struct.",
                        "Partial res. PC1 struct."),
                xwrite=c(-3,-1,-1),ywrite=c(-3,-1.5,-2))

for (ID in 1:length(list_resid$stat)){
  
  stat=list_resid$stat[ID]
  metric=list_resid$metric[ID]
  
  d_data_out=Get_data_without_outliers(stat,T)
  
  model_spa_stat=readRDS(paste0("../Data/Step1_Understanding_grazing/Keep_models/Mod_",stat,"_TRUE.rds"))
  
  resid_mod=visreg::visreg(fit = model_spa_stat,xvar = metric,by="Grazing",plot=F) #extracting residuals
  resid_mod$res$Grazing=rep(0:3,as.vector(table(d_data_out$Grazing))) #correcting vigreg output
  
  # Getting coeff across grazing gradient 
  beta_mean=round(coef(summary(lm(data=resid_mod$res%>%
                            melt(., measure.vars=metric,value.name = "value2"),visregRes~value2)))[2,1],2)
  beta_sd=round(coef(summary(lm(data=resid_mod$res%>%
                            melt(., measure.vars=metric,value.name = "value2"),visregRes~value2)))[2,2],2)
  pval=coef(summary(lm(data=resid_mod$res%>%
                            melt(., measure.vars=metric,value.name = "value2"),visregRes~value2)))[2,4]

  assign(paste0("p2_",ID),
         ggplot(resid_mod$res%>%
                  mutate(., Grazing=as.character(Grazing))%>%
                  dplyr::select(., -value)%>%
                  melt(., measure.vars=metric)%>%arrange(., Grazing,decreasing=T))+
           geom_point(aes(value,visregRes),color="gray30",alpha=.2)+
           geom_smooth(aes(value,visregRes,fill=Grazing,color=Grazing),method = "lm",level=.95)+
           # annotate("text",x=list_resid$xwrite[ID],y=list_resid$ywrite[ID],
           #               label=paste0("Est. = ",beta_mean," ± ",beta_sd,is_signif(pval)),size=4)+
           the_theme+
           scale_color_manual(values=c("0"="#FFC170",'1'="#43FFF4","2"="#5196F9","3"="#012E71"))+
           scale_fill_manual(values=c("0"="#FFC170",'1'="#43FFF4","2"="#5196F9","3"="#012E71"))+
           labs(x=ifelse(metric=="Org_C","Facilitation",metric),y=list_resid$title[ID]))
  
}

p=ggarrange(ggarrange(p1_1,p1_2,ncol=2,labels=letters[1:2],hjust = -12),
          ggarrange(p2_1,p2_2,p2_3,ncol=3,labels = c(letters[3],""),common.legend = T,legend = "bottom",hjust = -7),
          nrow=2,heights = c(1,.5))
ggsave("../Figures/Final_figs/Interaction_with_grazing.pdf",p,width = 10,height = 7)



## >> Figure 3: Indirect effects of grazing ----

d_data=read.table("../Data/Spatial_structure_grazing.csv",sep=";")%>%
  Closer_to_normality(.)

d_data[,c(1,9:29,35,37:41,44:ncol(d_data))] = apply(d_data[,c(1,9:29,35,37:41,44:ncol(d_data))],2,z_tranform)

#Getting the multifunctionality index
d_data$MF=z_tranform(rowSums(d_data[,55:(ncol(d_data))],na.rm = T))

#Adding the first two components from the PCA
d_data=Perform_PCA_spatial_struc(d_data)
d_indirect=tibble() #to save indirect effects of grazing on the spatial structure
dir.create("../Figures/Step1_Understanding_grazing/SEM",showWarnings = F)

param_list=expand.grid(Cov=c("with_cover"),
                       Stat=c("perim_area_scaling","PL_expo","fmax_psd",
                              "PLR","flow_length",
                              "contig","core_area"))

for (ID in 1:nrow(param_list)){
  
  dir.create(paste0("../Figures/Step1_Understanding_grazing/SEM/",param_list$Stat[ID]),showWarnings = F)
  
  with_cover=param_list$Cov[ID];k=as.character(param_list$Stat[ID])
  
  model_lmer=lmer(paste("value ~ (1|Site_ID) + Long_cos + Long_sin + Lattitude + Slope + Elevation + Clim1 + Clim2 + Clim3 + Clim4 + Type_veg",
                        ifelse(with_cover=="without_cover","+ rho_p","")),
                  data = d_data%>%melt(., measure.vars=k)%>%filter(., !is.na(value)),
                  na.action = na.omit,REML ="FALSE")
  
  #we remove potential outliers
  rm.outliers = romr.fnc(model_lmer, d_data%>%melt(., measure.vars=k)%>%filter(., !is.na(value)),
                         trim=2.5)
  d_data_out = rm.outliers$data
  
  #We first control for all covariates and extract the residuals
  model_lmer=lmer(paste("value ~ (1|Site_ID) + Long_cos + Long_sin + Lattitude + Slope + Elevation + Clim1 + Clim2 + Clim3 + Clim4 + Type_veg",
                        ifelse(with_cover=="without_cover","+ rho_p","")),
                  data = d_data_out,
                  na.action = na.fail,REML ="FALSE")
  
  resid_model=residuals(model_lmer) #extract residuals
  
  save=d_data[as.numeric(names(resid_model)),]%>% #add it to the dataframe 
    add_column(., Resid_mod=resid_model)

  #DOING the SEMs
  r.squaredGLMM(model_lmer)
  
  #SEM with all grazing intensity
  
  d_sem=save
  all_d=summary(psem(
    lm(Resid_mod ~ Grazing + Woody + rho_p , d_sem),
    lm(Woody ~ Grazing + Org_C + Sand, d_sem),
    lm(Sand ~ Grazing , d_sem),
    lm(rho_p ~ Grazing + Org_C + Sand , d_sem),
    lm(Org_C ~ Grazing, d_sem)
  ))
  
  #SEM with low grazing intensity
  
  d_sem=filter(save,Grazing %in% 0:1)
  low_graz=summary(psem(
    lm(Resid_mod ~ Grazing + Woody + rho_p , d_sem),
    lm(Woody ~ Grazing + Org_C + Sand, d_sem),
    lm(Sand ~ Grazing , d_sem),
    lm(rho_p ~ Grazing + Org_C + Sand , d_sem),
    lm(Org_C ~ Grazing, d_sem)
  ))
  
  #SEM with high grazing intensity
  
  d_sem=filter(save,Grazing %in% 2:3)
  high_graz=summary(psem(
    lm(Resid_mod ~ Grazing + Woody + rho_p , d_sem),
    lm(Woody ~ Grazing + Org_C + Sand, d_sem),
    lm(Sand ~ Grazing , d_sem),
    lm(rho_p ~ Grazing + Org_C + Sand , d_sem),
    lm(Org_C ~ Grazing, d_sem)
  ))
  
  d_indirect=rbind(d_indirect,
                   Get_indirect_effects_grazing(all_d)%>%
                     add_column(., Type="All",Stat=k),
                   Get_indirect_effects_grazing(high_graz)%>%
                     add_column(., Type="High",Stat=k),
                   Get_indirect_effects_grazing(low_graz)%>%
                     add_column(., Type="Low",Stat=k))
  
}


p=ggplot(d_indirect%>%
         filter(., Type !="All")%>%
         Rename_spatial_statistics(.)%>%
         mutate(., Stat=recode_factor(Stat,"mean_perim_area"="Mean ratio perimeter-area")))+
  geom_bar(aes(x=Indirect,y=Stat,fill=Type),position="dodge", stat="identity",width = .65)+
  facet_wrap(.~Path,scales = "free",nrow = 2,labeller = label_bquote(rows = "Path"==.(Path)))+
  scale_fill_manual(values=c("#3F6CE6","grey"))+
  labs(x="Indirect effect",y="",fill="")+
  the_theme+
  theme(axis.title.y = element_blank())+
  geom_vline(xintercept = 0)


ggsave("../Figures/Final_figs/Indirect_effects_SEM.pdf",p,
       width=4,height = 7)




# -------------------- SI figures ---------------------------


## >> Selection criteria for images ----
info_kmean=read.table("../Data/Landscapes/Kmean_clust_info.csv",sep=",",header = T)%>%
  filter(., Size!=200,Dataset=="biodesert")%>%
  add_column(., Quadratic_error=(.$field_cover-.$img_cover)**2)

trade_off=tibble(threshold=seq(0,.05,length.out=100))%>%
  add_column(., N_site_kept=sapply(1:nrow(.),function(x){
    return(length(which(info_kmean$Quadratic_error>.$threshold[x])))
  }))

p1=ggplot(trade_off)+
  geom_point(aes(x=N_site_kept,threshold))+
  the_theme+
  geom_hline(yintercept = 0.01,linetype=9,color="#FF9C4B",lwd=2)+
  labs(x="# of images removed",y="(Field - landscape cover)²")


p2=ggplot(info_kmean%>%
            mutate(., status=recode_factor(status,"kept"="Kept",
                                           "removed"="Removed (error > 0.01)",
                                           "removed_visual"="Removed after visual inspection")))+
  geom_point(aes(x=field_cover,y=img_cover,color=status))+
  geom_smooth(data=info_kmean%>%filter(., status=="kept"),
              aes(x=field_cover,y=img_cover),se = F,color="black",method = "lm")+
  geom_abline(slope=1,intercept = 0,linetype=9)+
  the_theme+
  labs(x="Field cover",y="Landscape cover",color="")+
  scale_color_manual(values=c("black","#FF9C4B","#D7A3FF"))

ggsave("../Figures/Final_figs/SI/Criteria_selection_images.pdf",
       ggarrange(p1,p2,ncol = 2,labels = letters[1:2],widths = c(.8,1)),
       width = 10,height = 4)

## >> Distribution of covariates across sites ----
d_data=read.table("../Data/Spatial_structure_grazing.csv",sep=";")%>%
  Closer_to_normality(.)

d_data[,c(1,9:29,35,37:41,44:ncol(d_data))] = apply(d_data[,c(1,9:29,35,37:41,44:ncol(d_data))],2,z_tranform)

#Getting the multifunctionality index
d_data$MF=z_tranform(rowSums(d_data[,55:(ncol(d_data))],na.rm = T))


p1=read.table("../Data/Spatial_structure_grazing.csv",sep=";")%>%
  Closer_to_normality(.)%>%
  dplyr::rename(., Cover=rho_p)%>%
  melt(., measure.vars=c("Cover","Woody"))%>%
  ggplot(.)+
  geom_density(aes(x=value),fill="#6DD275",alpha=.5)+
  facet_wrap(.~variable,scales = "free",ncol=2)+
  the_theme+
  theme(strip.background.x = element_blank())+
  labs(x="Value",y="")
  
  
p2=read.table("../Data/Spatial_structure_grazing.csv",sep=";")%>%
  Closer_to_normality(.)%>%
  dplyr::rename(., 'PC1 clim'=Clim1, 'PC3 clim'=Clim3,
                'PC2 clim'=Clim2, 'PC4 clim'=Clim4)%>%
  melt(., measure.vars=c("PC1 clim","PC2 clim","PC3 clim","PC4 clim"))%>%
  ggplot(.)+
  geom_density(aes(x=value),fill="#CC80C7",alpha=.5)+
  facet_wrap(.~variable,scales = "free",ncol=4)+
  the_theme+
  theme(strip.background.x = element_blank())+
  labs(x="Value",y="")

p3=read.table("../Data/Spatial_structure_grazing.csv",sep=";")%>%
  Closer_to_normality(.)%>%
  add_column(., MF=d_data$MF)%>%
  dplyr::rename(., Facilitation=Org_C,Multifunctionality=MF)%>%
  melt(., measure.vars=c("Facilitation","Sand","Slope","Multifunctionality"))%>%
  ggplot(.)+
  geom_density(aes(x=value),fill="#4564CE",alpha=.5)+
  facet_wrap(.~variable,scales = "free",ncol=4)+
  the_theme+
  theme(strip.background.x = element_blank())+
  labs(x="Value",y="")

p4=read.table("../Data/Spatial_structure_grazing.csv",sep=";")%>%
  Closer_to_normality(.)%>%
  dplyr::rename(., "Longitude (cos)"=Long_cos, "Longitude (sin)"=Long_sin)%>%
  melt(., measure.vars=c("Longitude (cos)","Longitude (sin)","Lattitude","Elevation"))%>%
  ggplot(.)+
  geom_density(aes(x=value),fill="#FFB15B",alpha=.5)+
  facet_wrap(.~variable,scales = "free",ncol=4)+
  the_theme+
  theme(strip.background.x = element_blank())+
  labs(x="Value",y="")


p_tot=ggarrange(ggarrange(ggplot()+theme_void(),p1,ggplot()+theme_void(),ncol=3,widths = c(.4,1,.4)),
                p2,
                p3,
                p4,nrow=4,labels=letters[1:4])
ggsave("../Figures/Final_figs/SI/Density_predictor_images.pdf",p_tot,width = 8,height = 10)
  






## >> Correlation between predictors ----

d_data=read.table("../Data/Spatial_structure_grazing.csv",sep=";")%>%
  Closer_to_normality(.)

d_data[,c(1,9:29,35,37:41,44:ncol(d_data))] = apply(d_data[,c(1,9:29,35,37:41,44:ncol(d_data))],2,z_tranform)
d_data=d_data%>%
  dplyr::rename(., Facilitation=Org_C,
                "Longitude (cos)"= Long_cos,
                "Longitude (sin)"= Long_sin,
                Cover=rho_p)

correlation_matrix=round(cor(d_data[,c(paste0("Clim",1:4),"Elevation","Slope",
                                 "Longitude (cos)",
                                 "Longitude (sin)",
                                 "Facilitation",
                                 "Cover","Sand","Woody")],use = "na.or.complete"),2)

correlation_matrix[lower.tri(correlation_matrix)]=NA

p=ggplot(correlation_matrix%>%melt(.))+
  geom_tile(aes(x=Var1,Var2,fill=value))+
  the_theme+
  geom_text(aes(x=Var1,Var2,label=value))+
  scale_fill_gradient2(low="red","white","blue",midpoint = 0,na.value = "white")+
  theme(axis.text.x = element_text(angle=60,hjust=1))+
  labs(x="",y="")

ggsave("../Figures/Final_figs/SI/Correlation_predictors.pdf",p,width = 6.5,height = 6.5)

## >> Grazing on vegetation cover ----

d_data=read.table("../Data/Spatial_structure_grazing.csv",sep=";")%>%
  Closer_to_normality(.)

d_data[,c(1,9:29,35,37:41,44:ncol(d_data))] = 
  apply(d_data[,c(1,9:29,35,37:41,44:ncol(d_data))],2,z_tranform)

d_data=Perform_PCA_spatial_struc(d_data)

formula_mod=formula(formula_(paste("rho_p ~ Long_cos + Long_sin + Lattitude + Slope + Elevation 
      + Clim1 + Clim2 + Clim3 + Clim4 + Grazing 
      + Woody + Woody*Grazing 
      + Clim1*Grazing + Clim2*Grazing + Clim3*Grazing + Clim4*Grazing  
      + Type_veg + Type_veg*Grazing
      + Sand + Sand*Grazing + Org_C + Org_C*Grazing
      + (1|Site_ID)")))



model_cover  = lmer(formula_mod, d_data,
                    na.action = na.fail ,REML ="FALSE")

#we remove potential outliers
rm.outliers = romr.fnc(model_cover, d_data, trim=2.5)
d_data_out = rm.outliers$data

model_cover = lmer(formula_mod, data = d_data_out, 
                   na.action = na.fail,REML ="FALSE")

# #do some model selection
# select_model=dredge(model_cover, subset = ~ Slope & Elevation &
#                       Long_cos & Long_sin & Lattitude &
#                       dc(Woody & Grazing, Woody : Grazing) &
#                       dc(Clim1 & Grazing, Clim1 : Grazing) &
#                       dc(Clim2 & Grazing, Clim2 : Grazing) &
#                       dc(Clim3 & Grazing, Clim3 : Grazing) &
#                       dc(Clim4 & Grazing, Clim4 : Grazing) &
#                       dc(rho_p & Grazing, rho_p : Grazing) &
#                       dc(Org_C & Grazing, Org_C : Grazing) &
#                       dc(Sand  & Grazing, Sand  : Grazing) &
#                       dc(Type_veg & Grazing, Type_veg : Grazing),
#                     options(na.action = "na.fail") )

#best model
model_cover=MuMIn::get.models(select_model, subset = delta == 0)[[1]]
# saveRDS(model_cover,paste0("../Data/Step1_Understanding_grazing/Keep_models/Mod_cover_TRUE.rds"))
model_cover=readRDS(paste0("../Data/Step1_Understanding_grazing/Keep_models/Mod_cover_TRUE.rds"))
R2_mod=r.squaredGLMM(model_cover)
boot_mod=bootstrap(model_cover,.f=fixef,type = "parametric",500)

#Merge in a df
d_pred=tibble(Median=apply(boot_mod$replicates[,-1],2,median),
              q1=apply(boot_mod$replicates[,-1],2,quantile,.025), 
              q3=apply(boot_mod$replicates[,-1],2,quantile,.975),
              pvalue=apply(boot_mod$replicates[,-1],2,twoside_pvalue),
              term=colnames(boot_mod$replicates)[-1],
              Stat="Cover",
              R2m=R2_mod[1])

loc_pval=.05+max(d_pred%>% #localization of pvalues in the plot
                   Organize_df(., "predictor")%>%
                   dplyr::select(.,q3)%>%dplyr::pull(.))

p1_1=ggplot(d_pred%>%
              filter(., term!="Type, both")%>%
              Organize_df(., "predictor")%>%
              add_column(., Is_signif_pval=sapply(1:nrow(.),function(x){return(is_signif(.$pvalue[x]))})))+
  geom_pointrange(aes(x=Median,y=term,xmin=q1,xmax=q3,color=Type_pred))+
  geom_text(aes(x=loc_pval,y=term,label=Is_signif_pval))+
  geom_hline(yintercept = c(5.5,17.5,21.5,25.5),linetype=9)+
  the_theme+
  labs(x="",y="",color="")+
  geom_vline(xintercept = 0,linetype=9)+
  scale_color_manual(values=c("Geography"="#FFB15B","Grazing"="#F3412A",
                              "Abiotic"="#4564CE","Vegetation"="#6DD275","Climatic"="#CC80C7"))+
  theme(legend.position = "none")

p1_2=ggplot(d_pred%>%
              add_column(., xplot=0)%>%
              Organize_df(., "bar")%>%
              group_by(., Type_pred,Stat,xplot)%>%
              dplyr::summarise(., sum_effect=sum(abs(Median)),.groups = "keep")%>%
              mutate(., sum_effect=sum_effect/sum(sum_effect)))+
  geom_bar(aes(x=xplot,y=sum_effect,fill=Type_pred),stat="identity",width = .05)+
  geom_text(aes(x=0,y=1.05,label=paste0("R² \n = \n ",unique(round(d_pred$R2m,2)))),size=3)+
  the_theme+
  labs(y="Relative effects of estimates",fill="")+
  scale_fill_manual(values=c("Geography"="#FFB15B","Grazing"="#F3412A",
                             "Abiotic"="#4564CE","Vegetation"="#6DD275","Climatic"="#CC80C7"))+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),axis.ticks.x = element_blank(),
        axis.line.x = element_blank())

p1=ggarrange(
  ggarrange(p1_1,p1_2+theme(legend.position = "none"),widths = c(1,.2),labels = letters[1:2]),
  get_legend(p1_2),nrow=2,heights = c(1,.15)
)

#extracting residuals of cover against grazing
resid_mod=visreg::visreg(fit = model_cover,xvar="Grazing",plot=F) 

p2=ggplot(NULL)+
  geom_jitter(data=resid_mod$res,
              aes(Grazing,visregRes),width = .1,color="gray",alpha=.5)+
  geom_pointrange(data=resid_mod$res%>%
                    dplyr::group_by(., Grazing)%>%
                    dplyr::summarise(., .groups = "keep",
                                     q1=quantile(visregRes,.25),q3=quantile(visregRes,.75),q2=median(visregRes)),
                  aes(x=Grazing,y=q2,ymin=q1,ymax=q3),color="black")+
  the_theme+
  scale_color_manual(values=c("0"="#FFC170",'1'="#43FFF4","2"="#5196F9","3"="#012E71"))+
  scale_fill_manual(values=c("0"="#FFC170",'1'="#43FFF4","2"="#5196F9","3"="#012E71"))+
  labs(x="Grazing",y="Residuals of cover")


#extracting residuals of cover against woody
resid_mod=visreg::visreg(fit = model_cover,xvar="Woody",plot=F) 

p3=ggplot(resid_mod$res)+
  geom_point(aes(Woody,visregRes),color="gray30",alpha=.2)+
  geom_smooth(aes(Woody,visregRes),color="black",fill="gray",method = "lm",level=.95)+
  the_theme+
  scale_color_manual(values=c("0"="#FFC170",'1'="#43FFF4","2"="#5196F9","3"="#012E71"))+
  scale_fill_manual(values=c("0"="#FFC170",'1'="#43FFF4","2"="#5196F9","3"="#012E71"))+
  labs(x="% woody cover",y="Residuals of cover")

#extracting residuals of cover against grazing
resid_mod=visreg::visreg(fit = model_cover,xvar="Org_C",by="Grazing",plot=F) 
resid_mod$res$Grazing=rep(0:3,as.vector(table(d_data_out$Grazing))) #correcting vigreg output

p4=ggplot(resid_mod$res%>%
            mutate(., Grazing=as.character(Grazing)))+
  geom_point(aes(Org_C,visregRes),color="gray30",alpha=.2)+
  geom_smooth(aes(Org_C,visregRes,fill=Grazing,color=Grazing),method = "lm",level=.95)+
  the_theme+
  scale_color_manual(values=c("0"="#FFC170",'1'="#43FFF4","2"="#5196F9","3"="#012E71"))+
  scale_fill_manual(values=c("0"="#FFC170",'1'="#43FFF4","2"="#5196F9","3"="#012E71"))+
  labs(x="Facilitation",y="Residuals on cover")

ggsave(paste0("../Figures/Final_figs/SI/Grazing_on_cover.pdf"),
       ggarrange(
         ggarrange(ggplot()+theme_void(),p1,ggplot()+theme_void(),ncol=3,widths = c(.25,1,.25)),
         ggarrange(
           ggarrange(p2,p3,p4+theme(legend.position="none"),ncol=3),
           get_legend(p4),nrow=2,heights = c(1,.3)),
         nrow=2,heights = c(1, .3)),
       width = 8,height = 9)



## >> Standardize predictors ----

cover=1
d_all=read.table(paste0("../Data/Step1_Understanding_grazing/Importance.csv"),sep=";")%>%
  filter(., With_cover==cover)%>%
  dplyr::rename(., Stat=Sp_stat)%>%
  Filter_relevant_stats(.)%>%
  Rename_spatial_statistics(.)

d_all2=read.table(paste0("../Data/Step1_Understanding_grazing/Estimators_model.csv"),sep=";")%>%
  Filter_relevant_stats(.)%>%
  Rename_spatial_statistics(.)

id=1
for (k in unique(d_all2$Stat)){
  
  loc_pval=.075+max(d_all2%>% #localization of pvalues in the plot
                     filter(., Stat==k , With_cover==cover)%>%
                     Organize_df(., "predictor")%>%
                     dplyr::select(.,q3)%>%dplyr::pull(.))
  
  assign(paste0("p_",id),ggplot(d_all2%>%
                filter(., Stat==k , With_cover==cover)%>%
                add_column(., loc_pval=loc_pval)%>%
                Organize_df(., "predictor")%>%
                add_column(., Is_signif_pval=sapply(1:nrow(.),
                                                    function(x){return(is_signif(.$pvalue[x]))})))+
    geom_pointrange(aes(x=Median,y=term,xmin=q1,xmax=q3,color=Type_pred))+
    geom_text(aes(x=loc_pval,y=term,label=Is_signif_pval))+
    the_theme+
    labs(x="",y="",color="")+
    ggtitle(k)+
    geom_vline(xintercept = 0,linetype=9)+
    scale_color_manual(values=c("Geography"="#FFB15B","Grazing"="#F3412A",
                                "Abiotic"="#4564CE","Vegetation"="#6DD275","Climatic"="#CC80C7"))+
    theme(legend.position = "none",title = element_text(size=9)))
  
  id=id+1
  
}
p_tot=ggarrange(p_1,p_2,p_3,p_4,ncol=2,nrow=2)
ggsave(paste0("../Figures/Final_figs/SI/Standardize_coef_a.pdf"),p_tot,width = 9,height = 7)
p_tot=ggarrange(p_5,p_6,p_7,p_8,ncol=2,nrow=2)
ggsave(paste0("../Figures/Final_figs/SI/Standardize_coef_b.pdf"),p_tot,width = 9,height = 7)



## >> Correlation metrics-cover ----
d_data=read.table("../Data/Spatial_structure_grazing.csv",sep=";")%>%
  Closer_to_normality(.)
d_data=Perform_PCA_spatial_struc(d_data)

p=ggplot(d_data%>%melt(., measure.vars=c("perim_area_scaling","fmax_psd","Cond_H","PL_expo","Spectral_ratio",
                                         "core_area_land","division","fractal_dim","contig","core_area",
                                         "flow_length","PLR",
                                         "Struct1","Struct2"))%>%
           dplyr::rename(., Stat=variable)%>%
           Filter_relevant_stats(.)%>%
           Rename_spatial_statistics(.))+
  geom_point(aes(x=rho_p,y=value),color="grey")+
  the_theme+
  facet_wrap(.~Stat,scales = "free",ncol = 4)+
  labs(x="Vegetation cover",y="Spatial metrics")

ggsave("../Figures/Final_figs/SI/Correlation_metrics_cover.pdf",p,width = 9,height = 4)

## >> Importance and effect of grazing without cover ----

d_slope=read.table("../Data/Step1_Understanding_grazing/Slope_partial_residuals.csv",sep=";")%>%
  filter(., With_cover==0)%>%
  Filter_relevant_stats(.)%>%  
  Rename_spatial_statistics(.)%>%
  add_column(., Signif=sapply(1:nrow(.),function(x){
    return(is_signif(.$pval[x]))}))%>%
  arrange(., slope)%>%
  add_column(., Order_f=1:nrow(.))%>%
  mutate(.,Stat = fct_reorder(Stat, Order_f))%>%
  add_column(., order_signif=sapply(1:nrow(.),
                                    function(x){return(ifelse(.$slope[x]<0,.075,-.12))}))

p2=ggplot(d_slope)+
  geom_pointrange(aes(x=slope,y=Stat,xmin=Low_int,xmax=High_int,color=pval<.05),shape=17,size=.8)+
  the_theme+
  scale_color_manual(values=c("grey","black"))+
  #geom_text(aes(x=order_signif,y=1:8,label=Signif))+
  geom_vline(xintercept = 0,linetype=9)+
  labs(x=substitute(paste(beta," (partial res. Grazing)")),y="",color="")+
  theme(legend.position = "none")



#Importance based on R2
d_R2=read.table(paste0("../Data/Step1_Understanding_grazing/Importance.csv"),sep=";")%>%
  filter(., With_cover==0)%>%
  dplyr::select(., Sp_stat,R2m,R2C)%>%
  dplyr::rename(., Stat=Sp_stat)%>%
  Filter_relevant_stats(.)%>%  
  Rename_spatial_statistics(.)%>%
  add_column(., Order_f=sapply(1:nrow(.),function(x){
    which(d_slope$Stat==.$Stat[x])
  }))%>%
  arrange(., Order_f)%>%
  mutate(.,Stat = fct_reorder(Stat, Order_f))


d_all2=read.table(paste0("../Data/Step1_Understanding_grazing/Estimators_model.csv"),sep=";")%>%
  dplyr::rename(., observed=Median)%>%
  filter(., With_cover==0)%>%
  Filter_relevant_stats(.)%>%  
  Organize_df(., "bar")%>%
  group_by(., Type_pred,Stat)%>%
  dplyr::summarise(., sum_effect=sum(abs(observed)),.groups = "keep")%>%
  Rename_spatial_statistics(.)%>%group_by(., Stat)%>%
  dplyr::summarise(., .groups="keep",sum_effect=sum_effect/sum(sum_effect),Type_pred=Type_pred)%>%
  add_column(., Order_f=sapply(1:nrow(.),function(x){
    which(d_slope$Stat==.$Stat[x])
  }))%>%
  arrange(., Order_f)%>%
  mutate(.,Stat = fct_reorder(Stat, Order_f))



p1=ggplot(d_all2)+
  geom_bar(aes(x=100*sum_effect,y=Stat,fill=Type_pred),stat="identity",width = .75)+
  geom_text(data=d_R2,aes(y=Stat,x=20,label=paste0("R²m = ",round(R2m,2),", ","R²c = ",round(R2C,2))),size=3)+
  the_theme+
  labs(x="Importance (%)",fill="")+
  scale_fill_manual(values=c("Geography"="#FFB15B","Grazing"="#F3412A",
                             "Abiotic"="#4564CE","Vegetation"="#6DD275","Climatic"="#CC80C7"))+
  theme(axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank())

p_tot=ggarrange(
  ggarrange(p1+theme(legend.position = "none"),
            p2+theme(axis.text.y = element_blank(),axis.ticks.y = element_blank()),
            ncol=2,labels = letters[1:2],widths = c(1,.5),align = "h"),
  ggarrange(ggplot()+theme_void(),get_legend(p1),ggplot()+theme_void(),ncol=3,widths = c(.2,1,.6)),
  nrow=2,heights = c(1,.1))



ggsave("../Figures/Final_figs/SI/Importance_effect_grazing_without_cover.pdf",
       p_tot,width = 8,height = 4)



## >> Importance ----

d_all=read.table(paste0("../Data/Step1_Understanding_grazing/Importance_nointer.csv"),sep=";")

with_cov=T

mean_imp=as_tibble(t(colMeans(d_all[which(d_all$With_cover==with_cov),c(-1,-ncol(d_all))])))

d_prep=rbind(cbind("Sp_stat"="Averaged across \n spatial statistics",
                   (mean_imp-mean_imp)%>% #to have white box
                     add_column(., With_cover=T)),
             cbind("Sp_stat"="Averaged across \n spatial statistics",
                   as_tibble(t(colMeans(d_all[which(d_all$With_cover==F),c(-1,-ncol(d_all))])))%>%
                     add_column(., With_cover=F)),
             d_all)%>%
  
  add_column(., "R2"=0)%>%
  dplyr::relocate(., R2, .before = With_cover)%>%
  filter(., With_cover==with_cov)%>%
  melt(., measure.vars=colnames(.)[5:(ncol(.)-1)])%>%
  mutate(., Sp_stat=recode_factor(Sp_stat,"fmax"="Max patch",
                                  "perim_area_scaling"="Perim-area scaling",
                                  "PL"="PSD exponent",
                                  "spectral_ratio"="Spectral ratio",
                                  "Cond_H"="Conditional entropy"))%>%
  mutate(., variable=recode_factor(variable,"Org_C"="Facilitation",
                                   "Interactions"="Interactions \n with grazing"))%>%
  #arranging order for predictors (xaxis)
  add_column(., Order_stat=rep(rev(1:(length(unique(d_all$Sp_stat))+1)),
                               nrow(.)/(length(unique(d_all$Sp_stat))+1)))%>%
  mutate(Sp_stat = fct_reorder(Sp_stat, Order_stat))%>%
  
  #arranging order for stats
  add_column(., Order_f=rep(c(sapply(c(4:ncol(mean_imp)),
                                     function(x){
                                       return(which(round(mean_imp[1,x]%>%pull(.),5) ==
                                                      round(sort(as.numeric(mean_imp[1,c(4:ncol(mean_imp))]),
                                                                 decreasing = T),5)))}),12),
                            each=(length(unique(d_all$Sp_stat))+1)))%>%
  mutate(variable = fct_reorder(variable, Order_f))

n_var=length(unique(d_prep$variable))
n_metric=length(unique(d_prep$Sp_stat))

p=ggplot(d_prep)+
  geom_tile(aes(x=variable,y=Sp_stat,fill=value))+
  geom_text(data=tibble(x=n_var,y=1:(n_metric-1),label=rev(round(unique(d_prep$R2m)[-1],2))),
            aes(x=x,y=y,label=paste0(label)),size=3)+
  geom_text(data=tibble(x=1:(n_var-1),y=(n_metric),
                        label=round(sort(as.numeric(
                          as_tibble(t(colMeans(d_all[d_all$With_cover==with_cov,
                                                     c(-1,-ncol(d_all))]))))[-c(1:3)],decreasing = T),2)),
            aes(x=x,y=y,label=paste0(label)),size=3)+
  scale_fill_gradientn(colours = colorRampPalette(c("white","#A129DA",alpha("#07299C",.3)))(100))+
  the_theme+
  theme(axis.text.x = element_text(angle=c(rep(60,(n_var-1)),0),hjust=c(rep(1,(n_var-1)),0)))+
  labs(x="Predictors",y="Spatial statistic",fill="Importance (%)")

ggsave(paste0("../Figures/Final_figs/SI/Importance_grazing_no_inter.pdf"),p,width = 6,height = 5)




