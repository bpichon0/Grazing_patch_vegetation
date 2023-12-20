rm(list=ls())
source("./Structure_grazing_function.R")

# -------------------- Main figures -------------------------

## >> Figure 1: Importance grazing with other variables ----

#Importance using multimodel selection 

with_cov=T
d_all=read.table(paste0("../Data/Linear_models/Importance_aridity_no_inter.csv"),sep=";")%>%
  dplyr::rename(., Stat=Sp_stat)%>%
  Filter_relevant_stats(.)%>%
  filter(., Grazing_intensity=="all",With_cover==with_cov)%>%
  dplyr::select(., -Interactions,-Grazing_intensity,-With_cover,-Woody)



mean_imp=as_tibble(t(colMeans(d_all[,-c(1,2)])))

d_prep=rbind(cbind("Stat"="Averaged across \n spatial statistics",
                   (mean_imp-mean_imp)), #to have white box
             d_all[,-2])%>%
  
  add_column(., "R2"=0)%>%
  dplyr::relocate(., R2, .after = Cover)%>%
  melt(., measure.vars=colnames(.)[4:(ncol(.))])%>%
  Rename_spatial_statistics(.)%>%
  mutate(., variable=recode_factor(variable,"Org_C"="Org. Carb.",
                                   "Interactions"="Interactions \n with grazing",
                                   "Type_veg"="Type vegetation"))%>%
  #arranging order for predictors (xaxis)
  add_column(., Order_stat=rep(rev(1:(length(unique(d_all$Stat))+1)),
                               nrow(.)/(length(unique(d_all$Stat))+1)))%>%
  mutate(Stat = fct_reorder(Stat, Order_stat))%>%
  Rename_spatial_statistics(.)%>%
  
  #arranging order for stats
  add_column(., Order_f=rep(c(sapply(c(3:ncol(mean_imp)),
                                     function(x){
                                       return(which(round(mean_imp[1,x]%>%pull(.),5) ==
                                                      round(sort(as.numeric(mean_imp[1,c(3:ncol(mean_imp))]),
                                                                 decreasing = T),5)))}),12),
                            each=(length(unique(d_all$Stat))+1)))%>%
  mutate(variable = fct_reorder(variable, Order_f))

n_var=length(unique(d_prep$variable))
n_metric=length(unique(d_prep$Stat))

p=ggplot(d_prep)+
  geom_tile(aes(x=variable,y=Stat,fill=value))+
  geom_text(data=tibble(x=n_var,y=1:(n_metric-1),label=rev(round(unique(d_prep$R2m)[-1],2))),
            aes(x=x,y=y,label=paste0(label)),size=3)+
  geom_text(data=tibble(x=1:(n_var-1),y=(n_metric),
                        label=round(sort(as.numeric(mean_imp)[-c(1,2)],decreasing = T),2)),
            aes(x=x,y=y,label=paste0(label)),size=3)+
  scale_fill_gradientn(colours = colorRampPalette(c("white","grey","grey30"))(100))+
  the_theme+
  theme(axis.text.x = element_text(angle=60,hjust=1))+
  labs(x="",y="Spatial statistic",fill="Akaike weight importance")

ggsave("../Figures/Importance_grazing_other_variables.pdf",
       p,width = 6,height = 5)


## >> Figure 2: Partial residuals high-low grazing intensity ----

d_slope=read.table("../Data/Linear_models/Slope_partial_residuals_aridity.csv",sep=";")%>%
  filter(., With_cover==0,Driver=="Grazing")%>%
  Filter_relevant_stats(.)%>%  
  Rename_spatial_statistics(.)%>%
  mutate(., Grazing_intensity=recode_factor(Grazing_intensity,
                                            "all"="All dataset",
                                            "low"="Low grazing",
                                            "high"="High grazing"))%>%
  arrange(., q2)%>%
  add_column(., Order_f=1:nrow(.))%>%
  mutate(.,Stat = fct_reorder(Stat, Order_f))%>%
  add_column(., Signif=.$pval<.05)%>%
  arrange(., Stat,q2)%>%
  add_column(., Order_stat=rep(1:3,nrow(.)/3))%>%
  mutate(.,Grazing_intensity = fct_reorder(Grazing_intensity, Order_stat))
  

p=ggplot(d_slope)+
  geom_pointrange(aes(x=q2,y=Stat,xmin=q1,xmax=q3,color=Grazing_intensity,shape=Signif),
                  size=.8,position=position_jitterdodge())+
  the_theme+
  geom_vline(xintercept = 0,linetype=9)+
  scale_color_manual(values=c("grey30","#B685CA","#92D0C5"))+
  scale_shape_manual(values=c(1,16))+
  labs(x=substitute(paste(beta," (partial res. Grazing)")),y="",color="")+
  theme(legend.box = "vertical")+guides(shape="none")

ggsave("../Figures/Grazing_intensity_partial_residuals.pdf",p,width = 6,height = 6)

## >> Figure 3: Indirect effects of grazing on PCs ----

Name_grazing=c("Low grazing intensity","High grazing intensity")
pdf(paste0("../Figures/SEM_consumption_recycling.pdf"),width = 12,height = 5)
stat="core_area_land"
par(mfrow=c(1,2))

id=1
for (graz in c("low","high")){
  
  save=Get_data_resid_SEM(stat,graz)
  
  d_sem=save#%>%filter(., Grazing %in% c(0,1))
  SEM_distance=summary((psem(
    lmer(Resid_mod ~  (1|Site_ID) + Org_C_v+Grazing+Sand, d_sem),
    lmer(Nitrate ~  (1|Site_ID) + Grazing +Aridity+Sand, d_sem),
    lm(Org_C_v ~  Nitrate  + Aridity+Sand+Grazing, d_sem),
    Sand%~~%Aridity
  )))
  
  Plot_SEM(SEM_distance,pdf_ = F,title_ = paste0(" ",Name_grazing[id],
                                                 "\n Goodness of fit: Fisher C stat = ",
                                                 round(SEM_distance$Cstat,2)[1],
                                                 ",df = ",SEM_distance$Cstat[2],
                                                 ", Pval = ",SEM_distance$Cstat[3]),
           name_var = "% landscape   \n covered by core    ",type_N = "Nitrate",
           label_cex = 1.1,edge_cex = 1.2)
  
  id=id+1
}
dev.off()

## >> Figure 4: Predictors of drylands stability ----

n_neigh=1
d_partial=read.table(paste0("../Data/Inferrence/Eby_",n_neigh,"_neigh","/Partial_residuals_grazing_no_cover_data_uncertainty.csv"),sep=";")%>%
  add_column(., Type_effect=paste0(.$Type,"_",.$Metric))%>%
  filter(., Type_effect %!in% c("Low_Aridity","High_Aridity"),Param!="Relative distance")%>%
  mutate(., Type_effect = recode_factor(Type_effect, "Full_Aridity"="Aridity","Low_Grazing"="Grazing (low intensity)",
                                        "High_Grazing"="Grazing (high intensity)","Full_Grazing"="Grazing (all intensities)"))%>%
  mutate(., Param = recode_factor(Param, "Absolute distance"="Distance to the tipping point"))%>%
  mutate(., Type_effect=as.character(Type_effect),Param=recode_factor(Param,
                                                                      "q"="Parameter q (spatial aggregation)"))%>%
  add_column(., Order=sapply(1:nrow(.),function(x){
    if (.$Param[x]=="Cover"){
      return(1)
    }else if (.$Param[x]=="Parameter q (spatial aggregation)"){
      return(2)
    }else{
      return(3)
    }
  }))

d_quantile=as.data.frame(d_partial%>%
  dplyr::group_by(., Param,Type_effect)%>%
  dplyr::summarise(.,q1=quantile(Stat,.025),
                   pval=twoside_pvalue(Stat),
                   q2=median(Stat),
                   q3=quantile(Stat,.975),.groups = "keep"))%>%
  add_column(., Signif=.$pval<.1)%>%
  add_column(., Order=sapply(1:nrow(.),function(x){
    if (.$Param[x]=="Cover"){
      return(1)
    }else if (.$Param[x]=="Parameter q (spatial aggregation)"){
      return(2)
    }else{
      return(3)
    }
  }))


p=ggplot(NULL)+
  geom_flat_violin(data=d_partial,
                   aes(y=Stat,x=Type_effect,fill=Type_effect),position = position_nudge(x = 0.1, y = 0), 
                   alpha = 0.8,scale = "width",width=.3)+
  geom_pointrange(data=d_quantile,
                    aes(y=q2,ymin=q1,ymax=q3,x=Type_effect,color=Type_effect,shape=Signif),
                  position = position_nudge(x = -0.1, y = 0))+
  coord_flip()+
  the_theme2+
  labs(x="",y="Standardized effect",fill="Fixed effect  ")+
  geom_hline(yintercept = 0,linetype=9)+
  guides(fill="none")+
  facet_wrap(.~fct_reorder(Param, Order))+
  scale_shape_manual(values=c(1,16))+
  theme(legend.position = "none")+
  scale_color_manual(values=c("Grazing (all intensities)"="grey30","Grazing (low intensity)"="#92D0C5",
                              "Grazing (high intensity)"="#B685CA","Aridity"="#E6BB6B"))+
  scale_fill_manual(values=c("Grazing (all intensities)"="grey50","Grazing (low intensity)"="#92D0C5",
                             "Grazing (high intensity)"="#B685CA","Aridity"="#E6BB6B"))

ggsave("../Figures/Stability_aridity_grazing.pdf",p,width=9,height = 4)


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
  geom_point(aes(threshold,N_site_kept))+
  the_theme+
  geom_vline(xintercept = 0.01,color="red",lwd=1)+
  labs(y="# of images removed",x="(Field - landscape cover)?")


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

ggsave("../Figures/SI/Criteria_selection_images.pdf",
       ggarrange(ggarrange(p1,p2+theme(legend.position = "none"),ncol = 2,labels = letters[1:2],widths = c(.8,1)),
                 ggarrange(ggplot()+theme_void(),get_legend(p2),ncol=3,widths = c(1,.7)),
                 nrow=2,heights = c(1,.1)),
       width = 7,height = 3.5)

## >> Distribution of covariates across sites ----
d_data=read.table("../Data/Spatial_structure_grazing.csv",sep=";")%>%
  Closer_to_normality(.)

d_data[,c(1,9:29,36,38:42,45:ncol(d_data))] = apply(d_data[,c(1,9:29,36,38:42,45:ncol(d_data))],2,z_tranform)

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
ggsave("../Figures/SI/Density_predictor_images.pdf",p_tot,width = 8,height = 10)
  


## >> Predictors and grazing intensity ----



d_slope=read.table("../Data/Linear_models/Slope_partial_residuals_aridity.csv",sep=";")%>%
  filter(., With_cover==0,Driver=="Grazing")%>%
  Filter_relevant_stats(.)%>%  
  Rename_spatial_statistics(.)%>%
  add_column(., Signif=sapply(1:nrow(.),function(x){
    return(is_signif(.$pval[x]))}))%>%
  arrange(., slope)%>%
  add_column(., Order_f=1:nrow(.))%>%
  mutate(.,Stat = fct_reorder(Stat, Order_f))%>%
  add_column(., order_signif=sapply(1:nrow(.),
                                    function(x){return(ifelse(.$slope[x]<0,.075,-.12))}))

p=ggplot(d_slope%>%mutate(., Grazing_intensity=recode_factor(Grazing_intensity,
                                                             "all"="All dataset",
                                                             "low"="Low grazing",
                                                             "high"="High grazing")))+
  geom_pointrange(aes(x=slope,y=Stat,xmin=Low_int,xmax=High_int,color=Grazing_intensity),shape=17,size=.8,position=position_jitter(height = .4))+
  the_theme+
  geom_vline(xintercept = 0,linetype=9)+
  scale_color_manual(values=c("grey","lightblue","lightgreen"))+
  labs(x=substitute(paste(beta," (partial res. Grazing)")),y="",color="")


ggsave("../Figures/SI/Grazing_intensity_partial_residuals_without_cover.pdf",p,width = 6,height = 7)


## >> Correlation between predictors ----

d_data=read.table("../Data/Spatial_structure_grazing.csv",sep=";")%>%
  Closer_to_normality(.)

d_data[,c(1,9:29,36,38:42,45:ncol(d_data))] = apply(d_data[,c(1,9:29,36,38:42,45:ncol(d_data))],2,z_tranform)

d_data=d_data%>%
  dplyr::rename(., "Organic C."=Org_C,
                "Longitude (cos)"= Long_cos,
                "Longitude (sin)"= Long_sin,
                "Total N"=Total_N,
                Cover=rho_p)

cor_mat_pred=psych::corr.test(d_data[,c("Aridity","Elevation","Slope",
                                        "Longitude (cos)",
                                        "Longitude (sin)",
                                        "Organic C.",
                                        "Cover","Sand")],adjust="none")
correlation_matrix=round(cor_mat_pred$r,2)

correlation_matrix[lower.tri(correlation_matrix)]=NA
diag(correlation_matrix)=NA

p=ggplot(correlation_matrix%>%
           melt(.)%>%
           add_column(., pval=sapply(1:nrow(.), function(x){
             if (is.na(.$value[x])){
               return(NA)
             }else{return(melt(cor_mat_pred$p)$value[x])}
           })))+
  geom_tile(aes(x=Var1,Var2,fill=value))+
  the_theme+
  geom_text(aes(x=Var1,Var2,label=ifelse(pval<.1,"","X")))+
  scale_fill_gradient2(low="red",mid="white",high = "blue",midpoint = 0,na.value = "white")+
  theme(axis.text.x = element_text(angle=60,hjust=1))+
  labs(x="",y="",fill="")

ggsave("../Figures/SI/Correlation_predictors.pdf",p,width = 5,height = 5)

## >> Box-plot spatial structure and grazing ----

d_data=read.table("../Data/Spatial_structure_grazing.csv",sep=";")%>%
  Closer_to_normality(.)

p=ggplot(d_data%>%
           melt(., measure.vars = c("PLR","PL_expo","fmax_psd","flow_length",
                                    "perim_area_scaling","core_area","contig",
                                    "core_area_land","KS_dist"))%>%
           mutate(., variable=recode_factor(variable,
                                            "contig"="Contiguity",
                                            "core_area_land"="% landscape covered by core",
                                            "core_area"="Mean % of core in patches",
                                            "flow_length"="Bare soil connectivity",
                                            "perim_area_scaling"="Fractal scaling area, perim.",
                                            "PLR"="PLR",
                                            "KS_dist"="Distance to null expect.",
                                            "PL_expo"="Power-law exp. of the PSD",
                                            "fmax_psd"="log (largest patch)")))+
  geom_boxplot(aes(x=Grazing,y=value,group=Grazing),outlier.shape = NA)+
  geom_jitter(aes(x=Grazing,y=value,color=Grazing>1),size=.5,alpha=.5)+
  scale_color_manual(values=c("#92D0C5","#B685CA"))+
  facet_wrap(.~variable,scales = "free")+
  the_theme2+
  labs(x="Grazing intensity",y="")+
  theme(legend.position = "none")

ggsave("../Figures/SI/Distribution_spatial_stat_grazing.pdf",p,width = 11,height = 7)

## >> Spatial metrics and aridity ----

d_data=read.table("../Data/Spatial_structure_grazing.csv",sep=";")%>%
  Closer_to_normality(.)

p=ggplot(d_data%>%
           melt(., measure.vars = c("PLR","PL_expo","fmax_psd","flow_length",
                                    "perim_area_scaling","core_area","contig",
                                    "core_area_land","KS_dist"))%>%
           mutate(., variable=recode_factor(variable,
                                            "contig"="Contiguity",
                                            "core_area_land"="% landscape covered by core",
                                            "core_area"="Mean % of core in patches",
                                            "flow_length"="Bare soil connectivity",
                                            "perim_area_scaling"="Fractal scaling area, perim.",
                                            "PLR"="PLR",
                                            "KS_dist"="Distance to null expect.",
                                            "PL_expo"="Power-law exp. of the PSD",
                                            "fmax_psd"="log (largest patch)")))+
  geom_point(aes(x=Aridity,y=value),size=.5,alpha=.5)+
  facet_wrap(.~variable,scales = "free")+
  the_theme2+
  labs(x="Grazing intensity",y="")+
  theme(legend.position = "none")

ggsave("../Figures/SI/Distribution_spatial_stat_aridity.pdf",p,width = 11,height = 7)

## >> Grazing on vegetation cover ----

d_effects=read.table(paste0("../Data/Linear_models/Estimators_model_aridity.csv"),sep=";")%>%
  filter(., With_cover==0,Stat=="rho_p")

loc_pval=.05+max(d_effects%>% #localization of p-values in the plot
                   Organize_df(., "predictor")%>%
                   dplyr::select(.,q3)%>%dplyr::pull(.))

#Standardize effects 

p1=ggplot(d_effects%>%
           filter(., term!="Type, both")%>%
           Organize_df(., "predictor")%>%
           add_column(., Is_signif_pval=sapply(1:nrow(.),function(x){
             if (sign(.$q1[x])==sign(.$q3[x])){return("*")
             }else{return("")}}))%>%
           mutate(., Grazing_intensity=recode_factor(Grazing_intensity,"all"="All dataset",
                                                     "low"="Low grazing intensity",
                                                     "high"="High grazing intensity")))+
  geom_pointrange(aes(x=Median,y=term,xmin=q1,xmax=q3,color=Type_pred))+
  geom_text(aes(x=loc_pval,y=term,label=Is_signif_pval))+
  the_theme2+
  labs(x="Standardized effect on vegetation cover",y="",color="")+
  geom_vline(xintercept = 0,linetype=9)+
  scale_color_manual(values=c("Geography"="#FFB15B","Grazing"="#F3412A",
                              "Abiotic"="#4564CE","Vegetation"="#6DD275","Climatic"="#CC80C7"))+
  theme(legend.position = "none")+facet_wrap(.~Grazing_intensity,scales = "free")


#Then, partial residuals

d_slope=read.table("../Data/Linear_models/Slope_partial_residuals_aridity.csv",sep=";")%>%
  filter(., With_cover==0)%>%
  filter(., Stat == "rho_p",Driver=="Grazing")%>%
  mutate(., Stat="Cover")%>%
  mutate(., Grazing_intensity=recode_factor(Grazing_intensity,
                                            "all"="Grazing (all intensities)",
                                            "low"="Grazing (low intensity)",
                                            "high"="Grazing (high intensity)"))%>%
  arrange(., q2)%>%
  add_column(., Signif=.$pval<.05)

p2=ggplot(d_slope)+
  geom_pointrange(aes(x=q2,y=Stat,xmin=q1,xmax=q3,color=Grazing_intensity,shape=Signif),
                  size=.8,position=position_jitterdodge())+
  the_theme+
  geom_vline(xintercept = 0,linetype=9)+
  scale_color_manual(values=c("Grazing (all intensities)"="grey30","Grazing (low intensity)"="#92D0C5",
                              "Grazing (high intensity)"="#B685CA"))+
  scale_shape_manual(values=c(1,16))+
  labs(x=substitute(paste(beta," (partial res. Grazing)")),y="",color="")+
  theme(legend.position = "right")+guides(shape="none",fill=guide_legend(ncol=3))

p_tot=ggarrange(p1,ggarrange(ggplot()+theme_void(),p2,ggplot()+theme_void(),ncol=3,widths = c(.3,1,.3)),
                nrow=2,labels=letters[1:2],heights = c(1,.2))

ggsave(paste0("../Figures/SI/Predictors_of_vegetation_cover.pdf"),p_tot, width = 14,height = 7)



## >> Importance of grazing and other variables with vegetation cover ----

#Importance using multimodel selection 
with_cov=F
d_all=read.table(paste0("../Data/Linear_models/Importance_aridity_no_inter.csv"),sep=";")%>%
  dplyr::rename(., Stat=Sp_stat)%>%
  Filter_relevant_stats(.)%>%
  filter(., Grazing_intensity=="all",With_cover==with_cov)%>%
  dplyr::select(., -Interactions,-Grazing_intensity,-With_cover,-Woody,-Cover)



mean_imp=as_tibble(t(colMeans(d_all[,-c(1,2)])))

d_prep=rbind(cbind("Stat"="Averaged across \n spatial statistics",
                   (mean_imp-mean_imp)), #to have white box
             d_all[,-2])%>%
  
  add_column(., "R2"=0)%>%
  dplyr::relocate(., R2, .after = Grazing)%>%
  melt(., measure.vars=colnames(.)[4:(ncol(.))])%>%
  Rename_spatial_statistics(.)%>%
  mutate(., variable=recode_factor(variable,"Org_C"="Organic C.",
                                   "Interactions"="Interactions \n with grazing"))%>%
  #arranging order for predictors (xaxis)
  add_column(., Order_stat=rep(rev(1:(length(unique(d_all$Stat))+1)),
                               nrow(.)/(length(unique(d_all$Stat))+1)))%>%
  mutate(Stat = fct_reorder(Stat, Order_stat))%>%
  Rename_spatial_statistics(.)%>%
  
  #arranging order for stats
  add_column(., Order_f=rep(c(sapply(c(3:ncol(mean_imp)),
                                     function(x){
                                       return(which(round(mean_imp[1,x]%>%pull(.),5) ==
                                                      round(sort(as.numeric(mean_imp[1,c(3:ncol(mean_imp))]),
                                                                 decreasing = T),5)))}),12),
                            each=(length(unique(d_all$Stat))+1)))%>%
  mutate(variable = fct_reorder(variable, Order_f))

n_var=length(unique(d_prep$variable))
n_metric=length(unique(d_prep$Stat))

p=ggplot(d_prep)+
  geom_tile(aes(x=variable,y=Stat,fill=value))+
  geom_text(data=tibble(x=n_var,y=1:(n_metric-1),label=rev(round(unique(d_prep$R2m)[-1],2))),
            aes(x=x,y=y,label=paste0(label)),size=3)+
  geom_text(data=tibble(x=1:(n_var-1),y=(n_metric),
                        label=round(sort(as.numeric(mean_imp)[-c(1,2)],decreasing = T),2)),
            aes(x=x,y=y,label=paste0(label)),size=3)+
  scale_fill_gradientn(colours = colorRampPalette(c("white","grey","grey30"))(100))+
  the_theme+
  theme(axis.text.x = element_text(angle=60,hjust=1))+
  labs(x="Predictors",y="Spatial statistic",fill="Frequency in \n best models")

ggsave("../Figures/SI/Importance_grazing_other_variables_without_cover.pdf",
       p,width = 6,height = 5)


## >> Partial residuals high-low grazing intensity (with vegetation cover) ----

d_slope=read.table("../Data/Linear_models/Slope_partial_residuals_aridity.csv",sep=";")%>%
  filter(., With_cover==1,Driver=="Grazing")%>%
  Filter_relevant_stats(.)%>%  
  Rename_spatial_statistics(.)%>%
  mutate(., Grazing_intensity=recode_factor(Grazing_intensity,
                                            "all"="All dataset",
                                            "low"="Low grazing",
                                            "high"="High grazing"))%>%
  arrange(., q2)%>%
  add_column(., Order_f=1:nrow(.))%>%
  mutate(.,Stat = fct_reorder(Stat, Order_f))%>%
  add_column(., Signif=.$pval<.05)%>%
  arrange(., Stat,q2)%>%
  add_column(., Order_stat=rep(1:3,nrow(.)/3))%>%
  mutate(.,Grazing_intensity = fct_reorder(Grazing_intensity, Order_stat))


p=ggplot(d_slope)+
  geom_pointrange(aes(x=q2,y=Stat,xmin=q1,xmax=q3,color=Grazing_intensity,shape=Signif),
                  size=.8,position=position_jitterdodge())+
  the_theme+
  geom_vline(xintercept = 0,linetype=9)+
  scale_color_manual(values=c("grey30","#B685CA","#92D0C5"))+
  scale_shape_manual(values=c(1,16))+
  labs(x=substitute(paste(beta," (partial res. Grazing)")),y="",color="")+
  theme(legend.box = "vertical")+guides(shape="none")


ggsave("../Figures/SI/Grazing_intensity_partial_residuals_with_cover.pdf",p,width = 6,height = 6)



## >> Standardize predictors ----

d_all2=read.table(paste0("../Data/Linear_models/Estimators_model_aridity.csv"),sep=";")%>%
  Filter_relevant_stats(.)%>%
  Rename_spatial_statistics(.)%>%
  filter(., Grazing_intensity=="all",With_cover==1)


id=1
for (k in unique(d_all2$Stat)){
  
  assign(paste0("p_",id),ggplot(d_all2%>%
                                  filter(., Stat==k)%>%
                                  Organize_df(., "predictor"))+
           geom_pointrange(aes(x=Median,y=term,xmin=q1,xmax=q3,color=Type_pred))+
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
ggsave(paste0("../Figures/SI/Standardize_coef_a_cov.pdf"),p_tot,width = 9,height = 9)
p_tot=ggarrange(p_5,p_6,p_7,p_8,ncol=2,nrow=2)
ggsave(paste0("../Figures/SI/Standardize_coef_b_cov.pdf"),p_tot,width = 9,height = 9)



## >> Standardize predictors no cover ----

d_all2=read.table(paste0("../Data/Linear_models/Estimators_model_aridity.csv"),sep=";")%>%
  Filter_relevant_stats(.)%>%
  Rename_spatial_statistics(.)%>%
  filter(., Grazing_intensity=="all",With_cover==0)

loc_pval=.05+max(d_all2%>% #localization of p-values in the plot
                   Organize_df(., "predictor")%>%
                   dplyr::select(.,q3)%>%dplyr::pull(.))

id=1
for (k in unique(d_all2$Stat)){
  
  assign(paste0("p_",id),ggplot(d_all2%>%
                                  filter(., Stat==k)%>%
                                  Organize_df(., "predictor")%>%
                                  add_column(., Is_signif_pval=sapply(1:nrow(.),function(x){
                                    if (sign(.$q1[x])==sign(.$q3[x])){return("*")
                                    }else{return("")}})))+
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
ggsave(paste0("../Figures/SI/Standardize_coef_a_no_cov.pdf"),p_tot,width = 9,height = 9)
p_tot=ggarrange(p_5,p_6,p_7,p_8,ncol=2,nrow=2)
ggsave(paste0("../Figures/SI/Standardize_coef_b_no_cov.pdf"),p_tot,width = 9,height = 9)



## >> Structural equation models (SEM) spatial stat, grazing intensity----

name_stat=c("log \n (largest patch)","Power-law exp. \n of the PSD","% landscape   \n covered by core    ",
            "Bare soil \n connectivity","Mean % of \n core in patches","Vegetation \n cover")
Name_grazing=c("Low grazing","High grazing")#,"All grazing")
ID=1

for (stat in c("fmax_psd","PL_expo","core_area_land","flow_length","core_area","rho_p")){
  
  pdf(paste0("../Figures/SI/SEM_cons_recycling_stat_",stat,".pdf"),width = 10,height = 7)#width = 14,height = 7)
  #layout(mat = matrix(c(1, 2, 3, 4,5,6),  nrow = 2, ncol = 3, byrow = F))
  layout(mat = matrix(c(1, 2, 3, 4),  nrow = 2, ncol = 2, byrow = F))
  
  id=1
  for (graz in c("low","high")){#,"all")){
    
    save=Get_data_resid_SEM(stat,graz)
    
    
    d_sem=save#%>%filter(., Grazing %in% c(0,1))
    SEM_distance=summary((psem(
      lmer(Resid_mod ~  (1|Site_ID) + Org_C_v+Grazing+Sand, d_sem),
      lmer(Nitrate ~  (1|Site_ID) + Grazing +Aridity+Sand, d_sem),
      lm(Org_C_v ~  Nitrate  + Aridity+Sand+Grazing, d_sem),
      Sand%~~%Aridity
    )))
    
    Plot_SEM(SEM_distance,pdf_ = F,title_ = paste0(Name_grazing[id],
                                                   "; Goodness of fit: Fisher C stat = ",
                                                   round(SEM_distance$Cstat,2)[1],
                                                   ",df = ",SEM_distance$Cstat[2],
                                                   ", Pval = ",SEM_distance$Cstat[3]),
             name_var = name_stat[ID],type_N = "Nitrate")
    
    SEM_distance=summary((psem(
      lmer(Resid_mod ~  (1|Site_ID) + Org_C_v+Grazing+Sand, d_sem),
      lmer(Total_N ~  (1|Site_ID) + Grazing +Aridity+Sand, d_sem),
      lm(Org_C_v ~  Total_N  + Aridity+Sand+Grazing, d_sem),
      Sand%~~%Aridity
    )))
    Plot_SEM(SEM_distance,pdf_ = F,title_ = paste0(Name_grazing[id],
                                                   "; Goodness of fit: Fisher C stat = ",
                                                   round(SEM_distance$Cstat,2)[1],
                                                   ",df = ",SEM_distance$Cstat[2],",
                                                   Pval = ",SEM_distance$Cstat[3]),
             name_var = name_stat[ID],type_N = "Total N.")
    id=id+1
  }
  dev.off()
  
  ID=ID+1
}



## >> Spatial resolution: correlation and AIC ----

#Correlation spatial resolution and spatial metrics

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
  geom_point(aes(x=Resolution,y=value),color="grey")+
  the_theme+
  facet_wrap(.~Stat,scales = "free",ncol = 4)+
  labs(x="Spatial resolution (m per pixel)",y="Spatial metrics")

ggsave("../Figures/SI/Correlation_metrics_Resolution.pdf",p,width = 9,height = 4)

#To evaluate the importance of spatial resolution, we compare the best models selected before
#with and without spatial resolution as a predictor

d_AIC=tibble()
list_mod=expand.grid(with_cover=c(T),
                     Stats=c("perim_area_scaling","fmax_psd","Cond_H","PL_expo","Spectral_ratio",
                             "core_area_land","division","fractal_dim","contig","core_area",
                             "flow_length","PLR",
                             "Struct1","Struct2")) #We also add the case where the stats correspond to the first two components of the PCA

for (id in 1:nrow(list_mod)){  
  
  stat=list_mod$Stats[id]
  
  d_data_out=  read.table(paste0("../Data/Linear_models/Keep_data/Data_",stat,"_TRUE.csv"),sep=";")
  if (ncol(d_data_out)==1){
    d_data_out=  read.table(paste0("../Data/Linear_models/Keep_data/Data_",stat,"_TRUE.csv"),sep=" ")
  }
  model_spa_stat=readRDS(paste0("../Data/Linear_models/Keep_models/Mod_",stat,"_TRUE.rds"))
  
  mod_with_res=lmer(formula = formula(
    paste(paste(formula(model_spa_stat))[2],paste(formula(model_spa_stat))[1],paste(paste(formula(model_spa_stat))[3],"+ Resolution")
    )),data = d_data_out,REML = F)
  
  d_AIC=rbind(d_AIC,tibble(AIC_no_res=AICc(mod_with_res,model_spa_stat)[2,2],
                           AIC_res=AICc(mod_with_res,model_spa_stat)[1,2],
                           Stat=stat))
}
d_AIC%>%Filter_relevant_stats(.)%>%Rename_spatial_statistics(.)%>%
  add_column(., Delta_AIC=.$AIC_no_res-.$AIC_res)

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

ggsave("../Figures/SI/Correlation_metrics_cover.pdf",p,width = 9,height = 4)


## >> Box-plot nutrients and grazing ----

d_data=read.table("../Data/Spatial_structure_grazing.csv",sep=";")%>%
  Closer_to_normality(.)

p=ggplot(d_data%>%
         melt(., measure.vars = c("Org_C","Org_C_v","Total_N","lnTotal_N","Nitrate","lnNitrate"))%>%
           mutate(., variable=recode_factor(variable,
                                            "lnTotal_N"="Total nitrogen (plants-bare)",
                                            "Total_N"="Total nitrogen",
                                            "Org_C"="Organic carbon (plants-bare)",
                                            "Org_C_v"="Organic carbon",
                                            "lnNitrate"="Nitrate (plants-bare)")))+
  geom_boxplot(aes(x=Grazing,y=value,group=Grazing),outlier.shape = NA)+
  geom_jitter(aes(x=Grazing,y=value,color=Grazing>1),size=.5,alpha=.5)+
  scale_color_manual(values=c("#92D0C5","#B685CA"))+
  facet_wrap(.~variable,scales = "free")+
  the_theme2+
  labs(x="Grazing intensity",y="")+
  theme(legend.position = "none")

ggsave("../Figures/SI/Distribution_N_C_cycling_grazing.pdf",p,width = 9,height = 5)


## >> Partial residuals high-low grazing/aridity N/C ----

d_slope=read.table("../Data/Linear_models_nutrients/Slope_partial_residuals_aridity_grazing.csv",sep=";")%>%
  add_column(., Driver_name=paste0(.$Driver,"_",.$Grazing_intensity))%>%
  mutate(., Driver_name=recode_factor(Driver_name,
                                      "Grazing_all"="All grazing",
                                      "Aridity_all"="Aridity",
                                      "Grazing_low"="Low grazing",
                                            "Grazing_high"="High grazing"))%>%
  filter(., Driver_name %in% c("All grazing","Aridity","Low grazing","High grazing"))%>%
  add_column(., Signif=.$pval<.05)%>%
  mutate(., Stat=recode_factor(Stat,
                               "lnTotal_N"="Total N (plants-bare)",
                               "Total_N"="Total N",
                               "Org_C"="Organic C. (plants-bare)",
                               "Org_C_v"="Organic C.",
                               "lnNitrate"="Nitrate (plants-bare)",
  ))



p=ggplot(d_slope)+
  geom_pointrange(aes(x=q2,y=Stat,xmin=q1,xmax=q3,color=Driver_name,shape=Signif),
                  size=.8,position=position_jitterdodge())+
  the_theme+
  geom_vline(xintercept = 0,linetype=9)+
  scale_color_manual(values=c("black","#EAB38C","#92D0C5","#B685CA"))+
  scale_shape_manual(values=c(1,16))+
  labs(x=substitute(paste(beta," (partial res. Grazing)")),y="",color="")+
  theme(legend.box = "vertical")+guides(shape="none")

ggsave("../Figures/SI/Nutrients_partial_res.pdf",p,width = 7,height = 5)



