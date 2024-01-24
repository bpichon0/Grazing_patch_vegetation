rm(list=ls())
source("./Structure_grazing_function.R")

# -------------------- Main figures -------------------------

## >> Figure 0: Chosen landscape ----
landscape=Get_empirical_site_all("biodesert",89,3)
psd_land=spatialwarnings::patchsizes(landscape>0)
psd_tibble=tibble(patch_size=psd_land)
psd_tibble$freq=sapply(1:nrow(psd_tibble),function(x){
  return(length(which(psd_tibble$patch_size>=psd_tibble$patch_size[x]))/nrow(psd_tibble))
})
p1=ggplot(psd_tibble)+
  geom_point(aes(x=patch_size,y=freq))+
  the_theme+
  scale_x_log10()+
  scale_y_log10()+
  labs(x="Patch size (x)",y=expression("Frequency (P">="x)"))

#identify largest patch

all_patches=get_patches(landscape,class = "1")$layer_1$class_1

ggplot(all_patches%>%melt(.)) +
  geom_raster(aes(x = Var1, y = Var2,
                  fill = (value==as.numeric(names(sort(table(all_patches)))[length(sort(table(all_patches)))])))) +
  coord_fixed() +
  theme_transparent() +theme(legend.position = "none")+
  scale_fill_manual(values=c("black","#1F4E79"),na.value = "white")


perim_patch=landscapemetrics::lsm_p_perim(raster(landscape))%>%filter(., class==1)%>%pull(., value)
area_patch=landscapemetrics::lsm_p_area(raster(landscape))%>%filter(., class==1)%>%pull(., value)
lsm_c_pafrac(raster(landscape), directions = 8, verbose = TRUE)%>%filter(., class==1)%>%pull(., value)

## >> Figure 1: Importance grazing with other variables ----

#Importance using multimodel selection 

with_cov=F
d_all=read.table(paste0("../Data/Linear_models_factor/Importance_aridity_no_inter_factor.csv"),sep=";")%>%
  dplyr::rename(., Stat=Sp_stat)%>%
  Filter_relevant_stats(.)%>%
  filter(., Grazing_intensity=="all",With_cover==with_cov)%>%
  dplyr::select(., -Interactions,-Grazing_intensity,-With_cover,-Woody)



mean_imp=as_tibble(t(colMeans(d_all[,-c(1,2)])))

d_prep=rbind(cbind("Stat"="Averaged across \n spatial statistics",
                   (mean_imp-mean_imp)), #to have white box
             d_all[,-2])%>%
  
  add_column(., "R2"=0)%>%
  dplyr::select(., -Cover)%>%
  dplyr::relocate(., R2, .after = Grazing)%>%
  melt(., measure.vars=colnames(.)[4:(ncol(.))])%>%
  add_column(., Type_stat=sapply(1:nrow(.),function(x){
    if (.$Stat[x] %in% c("PL_expo","fmax_psd")){
      return("Patch-size")
    }else if (.$Stat[x] =="flow_length"){
      return("Hydro.")
    }else if (.$Stat[x] =="Averaged across \n spatial statistics"){
      return("Hydro.")
    }else{
      return("Patch geometry")
    }
  }))%>%
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
                                                                 decreasing = T),5)))})),
                            each=(length(unique(d_all$Stat))+1)))%>%
  mutate(variable = fct_reorder(variable, Order_f))

n_var=length(unique(d_prep$variable))
n_metric=length(unique(d_prep$Stat))

p=ggplot(d_prep)+
  geom_tile(aes(x=variable,y=Stat,fill=value))+
  geom_text(data=tibble(x=n_var,y=1:(n_metric-1),label=(round(unique(d_prep$R2m)[-1],2))),
            aes(x=x,y=y,label=paste0(label)),size=3)+
  geom_text(data=tibble(x=1:(n_var-1),y=(n_metric),
                        label=round(sort(as.numeric(mean_imp)[-c(1,2,length(mean_imp))],decreasing = T),2)),
            aes(x=x,y=y,label=paste0(label)),size=3)+
  scale_fill_gradientn(colours = colorRampPalette(c("white","grey","grey30"))(100))+
  the_theme+
  theme(axis.text.x = element_text(angle=60,hjust=1))+
  labs(x="",y="Spatial statistic",fill="Akaike weight importance")

ggsave("../Figures/Importance_grazing_other_variables.pdf",
       p,width = 6,height = 5)


## >> Figure 2: Partial residuals grazing intensity spatial structure ----

d_slope=read.table("../Data/Linear_models_factor/Slope_partial_residuals_aridity_no_inter.csv",sep=";")%>%
  filter(., With_cover==0,Driver=="Grazing")%>%
  Filter_relevant_stats(.)%>%  
  Rename_spatial_statistics(.)%>%
  mutate(., ID_grazing=as.character(ID_grazing))%>%
  arrange(., q2)%>%
  add_column(., Order_f=1:nrow(.))%>%
  mutate(.,Stat = fct_reorder(Stat, Order_f))%>%
  add_column(., Signif=.$pval<.05)%>%
  arrange(., Stat,q2)%>%
  add_column(., Order_stat=rep(3:1,nrow(.)/3))%>%
  mutate(.,ID_grazing = fct_reorder(ID_grazing, Order_stat))%>%
  add_column(., Type_stat=sapply(1:nrow(.),function(x){
    if (.$Stat[x] %in% c("Power-law exp. of the PSD","log (largest patch)","Mean patch size")){
      return("Patch-size")
    }else if (.$Stat[x] =="Bare soil connectivity"){
      return("Hydro.")
    }else{
      return("Patch geometry")
    }
  }))


p=ggplot(d_slope)+
  geom_vline(xintercept = 0,linetype=9)+
  
  geom_linerange(aes(x=q2,y=Stat,xmin=q1_90,xmax=q3_90,color=as.factor(ID_grazing)),
                 lwd=2,position=position_jitterdodge(seed=123))+   
  geom_linerange(aes(x=q2,y=Stat,xmin=q1,xmax=q3,color=as.factor(ID_grazing)),
                 lwd=.8,position=position_jitterdodge(seed=123))+  
  geom_pointrange(aes(x=q2,y=Stat,xmin=q1,xmax=q3,fill=as.factor(ID_grazing)),color="transparent",
                 size=.3,position=position_jitterdodge(seed=123),shape=21)+  
  
  the_theme2+
  scale_color_manual(values=c("#FBD2A5","#FF8888","#C17F9D"),
                     labels=c("Low grazing (1)","Medium grazing (2)","High grazing (3)"))+
  scale_fill_manual(values=c("grey20","grey20","grey20"))+
  scale_shape_manual(values=c(21,16))+
  labs(x=substitute(paste(beta," (partial res. Grazing)")),y="",color="")+
  guides(shape="none",fill="none")+
  theme(legend.position="bottom")+
  facet_grid(Type_stat ~ ., space = "free_y", 
             scale = "free_y")+
  theme(panel.background = element_rect(fill="white"),legend.text = element_text(size=10))


ggsave("../Figures/Grazing_intensity_partial_residuals.pdf",p,width = 7,height = 5)


## >> Figure 3: Indirect effects of grazing on PCs ----

Name_grazing=c("Low grazing intensity","High grazing intensity","Grazed (1,2,3)","All data")
pdf(paste0("../Figures/SEM_consumption_recycling.pdf"),width = 12,height = 5)

par(mfrow=c(1,2))

id=1
for (stat in c("core_area_land","PL_expo","mean_psd")){
  for (graz in c("low","high","grazed","all")){
    
    save=Get_data_resid_SEM(stat,graz)
    
    d_sem=save#%>%filter(., Grazing %in% c(0,1))
    SEM_distance=summary((psem(
      lmer(Resid_mod ~  (1|Site_ID) + Org_C_v+Grazing+Sand, d_sem),
      lmer(Nitrate ~  (1|Site_ID) + Grazing +Aridity+Sand, d_sem),
      lmer(Org_C_v ~  (1|Site_ID) + Nitrate  + Aridity+Sand+Grazing, d_sem),
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
}
dev.off()


## >> Figure 3bis: Fertility island grazing ---- 

d_partial=read.table(paste0("../Data/Linear_models_fertility_factor/Slope_partial_residuals_aridity_grazing.csv"),sep=";")%>%
  filter(., With_cover==F)%>%
  mutate(., ID_grazing = recode_factor(ID_grazing, "none"="Aridity","1"="Grazing (low intensity)",
                                       "2"="Grazing (medium intensity)","3"="Grazing (high intensity)"))%>%
  add_column(., Signif=.$pval<.05)%>%
  filter(., Stat %in% c("lnTotal_N","Total_N","Org_C","Org_C_v","Nitrate","lnNitrate"))%>%
  mutate(., Stat=recode_factor(Stat,
                               "lnTotal_N"="Total N (plants-bare)",
                               "Total_N"="Total N",
                               "lnNitrate"="Nitrate (plants-bare)",
                               "Org_C"="Organic C. (plants-bare)",
                               "Org_C_v"="Organic C.",
  ))


p=ggplot(d_partial)+
  geom_linerange(aes(x=q2,y=Stat,xmin=q1_90,xmax=q3_90,color=as.factor(ID_grazing)),
                 lwd=2,position=position_jitterdodge(seed=123))+   
  geom_linerange(aes(x=q2,y=Stat,xmin=q1,xmax=q3,color=as.factor(ID_grazing)),
                 lwd=.8,position=position_jitterdodge(seed=123))+  
  geom_pointrange(aes(x=q2,y=Stat,xmin=q1,xmax=q3,fill=as.factor(ID_grazing)),color="transparent",
                  size=.3,position=position_jitterdodge(seed=123),shape=21)+  
  
  the_theme2+
  geom_vline(xintercept = 0,linetype=9)+
  scale_color_manual(values=c("lightgreen","#FBD2A5","#FF8888","#C17F9D"),
                     labels=c("Aridity","Low grazing (1)  ",
                              "Medium grazing (2)  ","High grazing (3)  "))+
  scale_fill_manual(values=c("grey20","grey20","grey20","grey20"))+
  scale_shape_manual(values=c(21,16))+
  labs(x=substitute(paste(beta," (partial res. Grazing)")),y="",color="")+
  theme(legend.box = "vertical")+guides(shape="none",fill="none")

ggsave("../Figures/Nutrients_partial_res.pdf",p,width = 7,height = 4)

## >> Figure 4: Predictors of drylands stability & SEM ----

# Figure 4a, grazing effects on the distance to the tipping point 

d_partial=read.table(paste0("../Data/Inferrence/Grazing_on_resilience_factor.csv"),sep=";")%>%
  filter(., With_inter==F)%>%
  mutate(., Grazing_id = recode_factor(Grazing_id, "none"="Aridity","1"="Grazing (low intensity)",
                                        "2"="Grazing (medium intensity)","3"="Grazing (high intensity)"))%>%
  mutate(., Param = recode_factor(Param, "Absolute distance"="Distance to the tipping point"))%>%
  mutate(., Param=recode_factor(Param,
                                                                      "q"="Aggregation parameter (q)"))%>%
  add_column(., Order=sapply(1:nrow(.),function(x){
    if (.$Param[x]=="Cover"){
      return(1)
    }else if (.$Param[x]=="Aggregation parameter (q)"){
      return(2)
    }else{
      return(3)
    }
  }))

d_quantile=as.data.frame(d_partial%>%
  dplyr::group_by(., Param,Grazing_id)%>%
  dplyr::summarise(.,q1=quantile(Stat,.025),
                   pval=twoside_pvalue(Stat),
                   q2=median(Stat),
                   q3=quantile(Stat,.975),
                   q1_90=quantile(Stat,.05),
                   q3_90=quantile(Stat,.95),
                   .groups = "keep"))%>%
  add_column(., Signif=.$pval<.1)%>%
  add_column(., Order=sapply(1:nrow(.),function(x){
    if (.$Param[x]=="Cover"){
      return(1)
    }else if (.$Param[x]=="Aggregation parameter (q)"){
      return(2)
    }else{
      return(3)
    }
}))


id=1 
for (k in unique(d_quantile$Param)[c(1,2,4)]){
  assign(paste0("p",id),
         ggplot(d_quantile%>%filter(.,Param==k,Grazing_id!="Aridity"))+
           geom_vline(xintercept = 0,linetype=9)+
           
           geom_linerange(aes(x=q2,y=Grazing_id,xmin=q1_90,xmax=q3_90,color=as.factor(Grazing_id)),
                          lwd=2,position=position_jitterdodge(seed=123))+   
           geom_linerange(aes(x=q2,y=Grazing_id,xmin=q1,xmax=q3,color=as.factor(Grazing_id)),
                          lwd=.8,position=position_jitterdodge(seed=123))+  
           geom_pointrange(aes(x=q2,y=Grazing_id,xmin=q1,xmax=q3,fill=as.factor(Grazing_id)),color="transparent",
                           size=.3,position=position_jitterdodge(seed=123),shape=21)+  
           
           the_theme2+
           scale_color_manual(values=c("#FBD2A5","#FF8888","#C17F9D"),
                              labels=c("Low grazing (1)  ",
                                       "Medium grazing (2)  ","High grazing (3)  "))+
           scale_fill_manual(values=c("grey20","grey20","grey20","grey20"))+
           scale_shape_manual(values=c(21,16))+
           labs(x=substitute(paste(beta," (partial res. Grazing)")),y="",color="")+
           guides(shape="none",fill="none")+
           theme(legend.position="bottom")+
           facet_wrap(.~Param,scales="free"))
  
  id=id+1
}

p_tot=ggarrange(p3+theme(axis.text.y = element_blank(),axis.ticks.y = element_blank()),
                p1+theme(axis.text.y = element_blank(),axis.ticks.y = element_blank()),
                p2+theme(axis.text.y = element_blank(),axis.ticks.y = element_blank())
                ,ncol=3,common.legend = T,legend="bottom")

p_tot=ggarrange(p_tot,ggplot()+theme_void(),nrow=2,labels = letters[1:2],heights = c(1,.75))
ggsave("../Figures/Stability_aridity_grazing.pdf",p_tot,width=7,height = 6)


# >> Figure 4b: SEM distance to the tipping point 

Total_effects=read.table("../Data/SEM/Total_effects_SEM_resilience.csv",sep=";")
Direct_effects=read.table("../Data/SEM/Direct_effects_SEM_resilience.csv",sep=";")


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

## >> Correlation between spatial stats ----

d_data=read.table("../Data/Spatial_structure_grazing.csv",sep=";")%>%
  Closer_to_normality(.)

d_data[,c(1,9:29,36,38:42,45:ncol(d_data))] = apply(d_data[,c(1,9:29,36,38:42,45:ncol(d_data))],2,z_tranform)


cor_mat_pred=psych::corr.test(d_data[,c("PLR","PL_expo","mean_psd","flow_length",
                                        "perim_area_scaling","core_area",
                                        "core_area_land","KS_dist","Shape_metric")],adjust="none")

rownames(cor_mat_pred$r)=colnames(cor_mat_pred$r)=c("PLR","Power-law exp. of the PSD","Mean patch size",
                                                    "Bare soil connectivity","Fractal scaling area, perim.",
                                                    "Mean % of core pixels in patches","% landscape covered by core pixels",
                                                    "Distance to null expect.","Patch shape complexity")

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
                                    "perim_area_scaling","core_area","Shape_metric",
                                    "core_area_land","KS_dist"))%>%
           mutate(., variable=recode_factor(variable,
                                            "Shape_metric"="Patch shape complexity",
                                            "core_area_land"="% landscape covered by core pixels",
                                            "core_area"="Mean % of core pixels in patches",
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
                                            "core_area_land"="% landscape covered by core pixels",
                                            "core_area"="Mean % of core pixels in patches",
                                            "flow_length"="Bare soil connectivity",
                                            "perim_area_scaling"="Fractal scaling area, perim.",
                                            "PLR"="PLR",
                                            "KS_dist"="Distance to null expect.",
                                            "PL_expo"="Power-law exp. of the PSD",
                                            "fmax_psd"="log (largest patch)")))+
  geom_point(aes(x=Aridity,y=value),size=.5,alpha=.5)+
  facet_wrap(.~variable,scales = "free")+
  the_theme2+
  labs(x="Aridity",y="")+
  theme(legend.position = "none")

ggsave("../Figures/SI/Distribution_spatial_stat_aridity.pdf",p,width = 11,height = 7)

## >> Importance of grazing and other variables with vegetation cover ----

#Importance using multimodel selection 
with_cov=T
d_all=read.table(paste0("../Data/Linear_models_factor/Importance_aridity_no_inter_factor.csv"),sep=";")%>%
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
  geom_text(data=tibble(x=n_var,y=1:(n_metric-1),label=(round(unique(d_prep$R2m)[-1],2))),
            aes(x=x,y=y,label=paste0(label)),size=3)+
  geom_text(data=tibble(x=1:(n_var-1),y=(n_metric),
                        label=round(sort(as.numeric(mean_imp)[-c(1,2)],decreasing = T),2)),
            aes(x=x,y=y,label=paste0(label)),size=3)+
  scale_fill_gradientn(colours = colorRampPalette(c("white","grey","grey30"))(100))+
  the_theme+
  theme(axis.text.x = element_text(angle=60,hjust=1))+
  labs(x="Predictors",y="Spatial statistic",fill="Frequency in \n best models")

ggsave("../Figures/SI/Importance_grazing_other_variables_with_cover.pdf",
       p,width = 6,height = 5)


## >> Standardize predictors ----

d_all=read.table(paste0("../Data/Linear_models_factor/Estimators_model_aridity_no_inter_factor.csv"),sep=";")%>%
  Filter_relevant_stats(.)%>%
  Rename_spatial_statistics(.)%>%
  filter(., Grazing_intensity=="all",With_cover==0)%>%
  add_column(., pval=sapply(1:nrow(.),function(x){
    if (sign(.$q1[x])==sign(.$q3[x])){
      return("*")
    }else{
      return("")
    }
  }))

loc_pval=.05+max(d_all%>% #localization of p-values in the plot
                   Organize_df(., "predictor")%>%
                   dplyr::select(.,q3)%>%dplyr::pull(.))


id=1
for (k in unique(d_all$Stat)){
  
  assign(paste0("p_",id),ggplot(d_all%>%
                                  filter(., Stat==k)%>%
                                  Organize_df(., "predictor")%>%
                                  add_column(., Is_signif_pval=sapply(1:nrow(.),function(x){
                                    if (sign(.$q1[x])==sign(.$q3[x])){return("*")
                                    }else{return("")}})))+
           geom_pointrange(aes(x=Median,y=term,xmin=q1,xmax=q3,color=Type_pred))+
           geom_text(aes(x=loc_pval,y=term,label=pval))+
           the_theme+
           labs(x="",y="",color="")+
           ggtitle(k)+
           geom_vline(xintercept = 0,linetype=9)+
           scale_color_manual(values=c("Geography"="#FFB15B","Grazing"="#F3412A",
                                       "Abiotic"="#4564CE","Vegetation"="#6DD275","Climatic"="#CC80C7"))+
           theme(legend.position = "none",title = element_text(size=9)))
  
  id=id+1
  
}
p_tot=ggarrange(p_1,p_2,p_3,ncol=3,nrow=1)
ggsave(paste0("../Figures/SI/Standardize_coef_a_no_cov.pdf"),p_tot,width = 12,height = 4)
p_tot=ggarrange(p_4,p_5,p_6,ncol=3,nrow=1)
ggsave(paste0("../Figures/SI/Standardize_coef_b_no_cov.pdf"),p_tot,width = 12,height = 4)



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
list_mod=expand.grid(with_cover=c(F,T),
                     Stats=c("PLR","PL_expo","fmax_psd","flow_length",
                             "perim_area_scaling","core_area",
                             "core_area_land","KS_dist","Shape_metric")) 

for (id in 1:nrow(list_mod)){  
  
  stat=list_mod$Stats[id]
  with_cover=list_mod$with_cover[id]
  
  d_data_out=  read.table(paste0("../Data/Linear_models_factor/Keep_data/Data_",stat,"_",with_cover,"_aridity_all.csv"),sep=";") 
  if (ncol(d_data_out)==1){
    d_data_out=  read.table(paste0("../Data/Linear_models_factor/Keep_data/Data_",stat,"_",with_cover,"_aridity_all.csv"),sep=" ")
  }
  model_spa_stat=readRDS(paste0("../Data/Linear_models_factor/Keep_models/Mod_",stat,"_",with_cover,"_aridity_all.rds"))
  
  mod_with_res=lmer(formula = formula(
    paste(paste(formula(model_spa_stat))[2],paste(formula(model_spa_stat))[1],paste(paste(formula(model_spa_stat))[3],"+ Resolution") #we add the spatial resolution to the model
    )),data = d_data_out,REML = F)
  
  d_AIC=rbind(d_AIC,tibble(AIC_no_res=AICc(mod_with_res,model_spa_stat)[2,2],
                           AIC_res=AICc(mod_with_res,model_spa_stat)[1,2],
                           Stat=stat,
                           With_cover=with_cover))
}
d_AIC%>%Filter_relevant_stats(.)%>%Rename_spatial_statistics(.)%>%
  add_column(., Delta_AIC=.$AIC_no_res-.$AIC_res)

## >> Correlation metrics-cover ----

d_data=read.table("../Data/Spatial_structure_grazing.csv",sep=";")%>%
  Closer_to_normality(.)
d_data=Perform_PCA_spatial_struc(d_data)

p=ggplot(d_data%>%melt(., measure.vars=c("perim_area_scaling","fmax_psd","PL_expo",
                                         "core_area_land","division","fractal_dim","core_area",
                                         "flow_length","PLR","Shape_metric",
                                         "KS_dist","mean_psd"))%>%
           dplyr::rename(., Stat=variable)%>%
           Filter_relevant_stats(.)%>%
           Rename_spatial_statistics(.))+
  geom_point(aes(x=rho_p,y=value),color="grey")+
  the_theme+
  facet_wrap(.~Stat,scales = "free",ncol = 3)+
  labs(x="Vegetation cover",y="Spatial metrics")

ggsave("../Figures/SI/Correlation_metrics_cover.pdf",p,width = 7,height = 6)


## >> Box-plot nutrients and grazing ----

d_data=read.table("../Data/Spatial_structure_grazing.csv",sep=";")%>%
  Closer_to_normality(.)%>%
  melt(., measure.vars = c("Org_C","Org_C_v","Total_N","lnTotal_N"))%>%
  mutate(., variable=recode_factor(variable,
                                   "lnTotal_N"="Total nitrogen (plants-bare)",
                                   "Total_N"="Total nitrogen",
                                   "Org_C"="Organic carbon (plants-bare)",
                                   "Org_C_v"="Organic carbon",
                                   "lnNitrate"="Nitrate (plants-bare)"))

p=ggplot(d_data)+
  geom_boxplot(aes(x=Grazing,y=value,group=Grazing),outlier.shape = NA)+
  geom_jitter(aes(x=Grazing,y=value,color=Grazing>1),size=.5,alpha=.5,height = 0)+
  scale_color_manual(values=c("#92D0C5","#B685CA"))+
  facet_wrap(.~variable,scales = "free")+
  the_theme2+
  labs(x="Grazing intensity",y="")+
  theme(legend.position = "none")

ggsave("../Figures/SI/Distribution_N_C_cycling_grazing.pdf",p,width = 6,height = 5)


## >> Moving average aridity ----

d_moving=read.table("../Data/Moving_average_aridity/Moving_average_aridity.csv",sep=";")%>%
  filter(., With_cover==0)%>%
  Filter_relevant_stats(.)%>%  
  Rename_spatial_statistics(.)

p=ggplot(d_moving)+
  geom_line(aes(x=Aridity_thresh,y=q2),lwd=1)+
  geom_ribbon(aes(x=Aridity_thresh,y=q2,ymin=q1,ymax=q3),alpha=.6,fill="grey")+
  the_theme2+
  geom_hline(yintercept = 0,color="#9CD891")+
  facet_wrap(.~Stat,scales = "free")+
  labs(x="Aridity",y=substitute(paste(beta," (partial res. aridity)")))


ggsave("../Figures/SI/Moving_average_aridity_sand.pdf",p,width = 7,height = 5)




