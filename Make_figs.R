rm(list=ls())
source("./Structure_grazing_function.R")

# -------------------- Main figures -------------------------

## >> Figure 1: Importance grazing with other variables ----

#Importance using multimodel selection 
with_cov=F
d_all=read.table(paste0("../Data/Linear_models_factor_cover_control/Importance_aridity_no_inter_factor.csv"),sep=";")%>%
  dplyr::rename(., Stat=Sp_stat)%>%
  Filter_relevant_stats(.)%>%
  filter(., Grazing_intensity=="all",With_cover==with_cov)%>%
  dplyr::select(., -Interactions,-Grazing_intensity,-With_cover,-Woody)



mean_imp=as_tibble(t(colMeans(d_all[,-c(1,2,ncol(d_all))])))

d_prep=rbind(cbind("Stat"="Averaged across \n spatial statistics",
                   (mean_imp-mean_imp)), #to have white box
             d_all[,-c(2,ncol(d_all))])%>%
  
  # add_column(., "R2"=0)%>%
  # dplyr::select(., -Cover)%>%
  # dplyr::relocate(., R2, .after = Grazing)%>%
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
                                                                 decreasing = T),5)))})),
                            each=(length(unique(d_all$Stat))+1)))%>%
  mutate(variable = fct_reorder(variable, Order_f))


n_var=length(unique(d_prep$variable))
n_metric=length(unique(d_prep$Stat))

p=ggplot(d_prep)+
  geom_tile(aes(x=variable,y=Stat,fill=value))+
  geom_text(data=tibble(x=1:(n_var),y=n_metric,
                        label=round(sort(as.numeric(mean_imp)[-c(1,2)],decreasing = T),2)),
            aes(x=x,y=y,label=paste0(label)),size=3)+
  scale_fill_gradientn(colours = colorRampPalette(c("white","grey","grey30"))(100))+
  the_theme+
  # facet_grid(Type_stat ~ ., space = "free_y", 
  #            scale = "free_y")+
  theme(panel.background = element_rect(fill="white"),panel.border = element_rect(colour = "transparent", fill=NA))+
  theme(axis.text.x = element_text(angle=60,hjust=1))+
  labs(x="",y="Spatial statistic",fill="Akaike weight importance")

ggsave("../Figures/Importance_grazing_other_variables.pdf",
       p,width = 6,height = 5)


## >> Figure 2: Partial residuals grazing intensity spatial structure ----

d_slope=read.table("../Data/Linear_models_factor_cover_control/Slope_partial_residuals_aridity_no_inter.csv",sep=";")%>%
  filter(., With_cover==0,Driver=="Grazing",Grazing_intensity=="all")%>%
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
    if (.$Stat[x] %in% c("Power-law exp. of the PSD","log (largest patch)",
                         "Mean patch size","Number of smallest patches")){
      return("Patch-size")
    }else if (.$Stat[x] =="Bare soil connectivity"){
      return("Hydro.")
    }else if (.$Stat[x] =="Spatial autocorrelation of vege."){
      return("Autocorr.")
    }else{
      return("Geom.")
    }
  }))


p=ggplot(d_slope)+
  geom_vline(xintercept = 0,linetype=9)+
  
  geom_linerange(aes(x=q2,y=Stat,xmin=q1_90,xmax=q3_90,color=as.factor(ID_grazing)),
                 lwd=1.5,position=position_jitterdodge(seed=123))+   
  geom_linerange(aes(x=q2,y=Stat,xmin=q1,xmax=q3,color=as.factor(ID_grazing)),
                 lwd=.5,position=position_jitterdodge(seed=123))+  
  geom_pointrange(aes(x=q2,y=Stat,xmin=q1,xmax=q3,fill=as.factor(ID_grazing)),color="transparent",
                  size=.15,position=position_jitterdodge(seed=123),shape=21)+  
  
  the_theme2+
  scale_color_manual(values=c("#FBD2A5","#FF8888","#C17F9D"),
                     labels=c("Low grazing (1)","Medium grazing (2)","High grazing (3)"))+
  scale_fill_manual(values=c("grey20","grey20","grey20"))+
  scale_shape_manual(values=c(21,16))+
  labs(x=substitute(paste(beta," (partial res. Grazing)")),y="",color="")+
  guides(shape="none",fill="none")+
  facet_grid(Type_stat ~ ., space = "free_y", 
             scale = "free_y")+
  theme(panel.background = element_rect(fill="white"),
        legend.text = element_text(size=8))+
  guides(color=guide_legend(ncol=1))+
  theme(legend.position="bottom")#c(.2,.65))

ggsave("../Figures/Grazing_intensity_partial_residuals.pdf",p,width = 4,height = 5)


## >> Figure 3: Illustrating the effects with simulations ----

list_landscape=readRDS("../Data/Simulations/Minimal_examples_stats_landscapes.rds")
d=read.table("../Data/Simulations/Minimal_examples_stats.csv",sep=";")

p_stat=ggplot(d%>%
                add_column(., ID=rep(c("High","Medium","Low","Random"),3),
                           Cover=rep(c("High","Medium","Low"),each=4))%>%
                melt(.,id.vars = c("ID","Cover"))%>%
                add_column(., Order=rep(4:1,nrow(.)/4))%>%
                mutate(.,ID = fct_reorder(ID, Order))%>%
                # filter(., ID!="Medium")%>%
                dplyr::rename(., Stat=variable)%>%
                Filter_relevant_stats(.)%>%
                filter(.,ID!="Random")%>%
                mutate(., Stat=recode_factor(Stat,
                                             "core_area"="Mean % of core \n pixels in patches",
                                             "mean_psd"="Mean patch size",
                                             "flow_length"="Bare soil connectivity",
                                             "fractal_dim"="Fractal dimension",
                                             "Small_patches"="# of smallest patches",
                                             "moran_I"="Spatial autocorr. of vege.",
                                             "PL_expo"="Power-law exp. \n of the PSD",
                                             "fmax_psd"="log (largest patch)")))+
  geom_line(aes(x=ID,y=value,group=interaction(Stat,Cover),color=Cover),lwd=2)+
  geom_point(aes(x=ID,y=value),shape=21,size=3,fill="white",color="black")+
  facet_wrap(.~Stat,scales = "free",nrow=2)+
  the_theme2+
  theme(axis.text.x =element_text(angle = 60,hjust = 1,size=12),
        legend.title = element_text(size=14),
        legend.text = element_text(size=14),
        strip.text.x.top = element_text(size=13),
        axis.title.x = element_text(size=14))+
  scale_fill_manual(values=c("#D27C3F","#BE91D0","#8BD882"))+
  scale_color_manual(values=c("#D27C3F","#BE91D0","#8BD882"))+
  labs(x="Aggregation level",y="",fill="Vegetation cover")+
  guides(color="none",fill="none")

theme_land=theme(
  panel.background = element_rect(fill = "transparent",
                                  colour = NA_character_), 
  panel.grid.major = element_blank(), 
  panel.grid.minor = element_blank(), 
  plot.background = element_rect(fill = "transparent",
                                 colour = NA_character_), 
  legend.background = element_rect(fill = "transparent"),
  legend.box.background = element_rect(fill = "transparent"),
  legend.key = element_rect(fill = "transparent"),
  axis.title = element_blank(),
  axis.ticks = element_blank(),
  axis.text = element_blank()
)

id=1
for (cover_id in 1:3){
  for (aggregation_id in 1:3){
    assign(paste0("p",aggregation_id,cover_id),
           Plot_lanscape(list_landscape[[id]])+theme_land+
             theme(plot.margin = unit(c(0,2,0,-2), 'lines')))
    id=id+1
  }
}


p_landscapes=ggarrange(p13+
                         theme(panel.border = element_rect(colour = "#D27C3F",fill=NA,size=2)),
                       p12+
                         theme(panel.border = element_rect(colour = "#D27C3F",fill=NA,size=2)),
                       p11+
                         theme(panel.border = element_rect(colour = "#D27C3F",fill=NA,size=2)),
                       p23+
                         theme(panel.border = element_rect(colour = "#8BD882",fill=NA,size=2)),
                       p22+
                         theme(panel.border = element_rect(colour = "#8BD882",fill=NA,size=2)),
                       p21+
                         theme(panel.border = element_rect(colour = "#8BD882",fill=NA,size=2)),
                       p33+
                         theme(panel.border = element_rect(colour = "#BE91D0",fill=NA,size=2)),
                       p32+
                         theme(panel.border = element_rect(colour = "#BE91D0",fill=NA,size=2)),
                       p31+
                         theme(panel.border = element_rect(colour = "#BE91D0",fill=NA,size=2)),
                       nrow=3,ncol=3)
p_arrow_top=ggplot(NULL) + 
  geom_line(data=tibble(x=c(0,1),y=c(1,1)),
            aes(x=x,y=y),
            arrow = arrow(length=unit(0.30,"cm"), type = "closed"))+
  theme_transparent()+
  ylim(.96,1.2)+
  xlim(0,1)+
  annotate("text",x=.5,y=1.1,label = "Increasing spatial aggregation",size=6,
           family = "NewCenturySchoolbook")+
  theme(plot.margin = unit(rep(.5,4), 'lines'))

p_arrow_left=ggplot(NULL) + 
  geom_line(data=tibble(x=c(1,1),y=c(0,1)),
            aes(x=x,y=y),
            arrow = arrow(length=unit(0.30,"cm"), type = "closed"))+
  annotate("text",x=.9,y=.5,label = "Increasing vegetation cover",angle = 90,size=6,
           family = "NewCenturySchoolbook")+
  theme_transparent()+
  ylim(0,1)+
  xlim(.8,1.1)+
  theme(plot.margin = unit(rep(.5,4), 'lines'))

p_landscapes=ggarrange(p_arrow_top,p_landscapes,nrow=2,heights = c(.1,1))
p_landscapes=ggarrange(ggarrange(ggplot()+theme_void(),p_arrow_left,nrow=2,heights = c(.1,1)),
                       p_landscapes,ncol=2,widths = c(.1,1))

p_tot=ggarrange(p_landscapes,p_stat,nrow=2,heights = c(1.15,1),labels = c("a","b"))
ggsave("../Figures/Examples_stat_toy_landscapes.pdf",p_tot,width = 10,height = 14)

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

cor_mat_pred=psych::corr.test(d_data[,c("PL_expo","fmax_psd","flow_length",
                                        "moran_I","core_area",
                                        "mean_psd")],adjust="none")

rownames(cor_mat_pred$r)=colnames(cor_mat_pred$r)=c("Power-law exp. of the PSD","log (largest patch)",
                                                    "Bare soil connectivity","Spatial autocorrelation of vege.",
                                                    "Mean % of core pixels in patches",
                                                    #"Number of smallest patches",
                                                    "Mean patch size")

correlation_matrix=round(cor_mat_pred$r,2)

correlation_matrix[lower.tri(correlation_matrix)]=NA
diag(correlation_matrix)=NA

p1=ggplot(correlation_matrix%>%
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



d_data=read.table("../Data/Spatial_structure_grazing_control_cover.csv",sep=";")

cor_mat_pred=psych::corr.test(d_data[,c("PL_expo","fmax_psd","flow_length",
                                        "moran_I","core_area",
                                        "mean_psd")],adjust="none")

rownames(cor_mat_pred$r)=colnames(cor_mat_pred$r)=c("Power-law exp. of the PSD","log (largest patch)",
                                                    "Bare soil connectivity","Spatial autocorrelation of vege.",
                                                    "Mean % of core pixels in patches",
                                                    #"Number of smallest patches",
                                                    "Mean patch size")

correlation_matrix=round(cor_mat_pred$r,2)

correlation_matrix[lower.tri(correlation_matrix)]=NA
diag(correlation_matrix)=NA

p2=ggplot(correlation_matrix%>%
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

ggsave("../Figures/SI/Correlation_spatial_stats.pdf",
       ggarrange(p1+ggtitle("Raw spatial statistics"),
                 p2+ggtitle("With cover correction"),ncol=2),
       width = 10,height = 5)

## >> Box-plot spatial structure and grazing ----

d_data=read.table("../Data/Spatial_structure_grazing.csv",sep=";")%>%
  Closer_to_normality(.)

p=ggplot(d_data%>%
           melt(., measure.vars = c("PL_expo","fmax_psd","flow_length",
                                    "moran_I","core_area",
                                    "mean_psd"))%>%
           mutate(., variable=recode_factor(variable,
                                            "core_area"="Mean % of core pixels in patches",
                                            "flow_length"="Bare soil connectivity",
                                            "mean_psd"="Mean patch size",
                                            "moran_I"="Spatial autocorrelation of vege.",
                                            "PL_expo"="Power-law exp. of the PSD",
                                            "fmax_psd"="log (largest patch)")))+
  geom_boxplot(aes(x=Grazing,y=value,group=as.factor(Grazing)),outlier.shape = NA)+
  geom_jitter(aes(x=Grazing,y=value,color=as.factor(Grazing)),size=.5,alpha=.5)+
  scale_color_manual(values=c("grey","#FBD2A5","#FF8888","#C17F9D"))+
  facet_wrap(.~variable,scales = "free")+
  the_theme2+
  labs(x="Grazing intensity",y="")+
  theme(legend.position = "none")

ggsave("../Figures/SI/Distribution_spatial_stat_grazing.pdf",p,width = 11,height = 7)



d_data=read.table("../Data/Spatial_structure_grazing_control_cover.csv",sep=";")

p=ggplot(d_data%>%
           melt(., measure.vars = c("PL_expo","fmax_psd","flow_length",
                                    "moran_I","core_area",
                                    "mean_psd"))%>%
           mutate(., variable=recode_factor(variable,
                                            "core_area"="Mean % of core pixels in patches",
                                            "flow_length"="Bare soil connectivity",
                                            "mean_psd"="Mean patch size",
                                            "moran_I"="Spatial autocorrelation of vege.",
                                            "PL_expo"="Power-law exp. of the PSD",
                                            "fmax_psd"="log (largest patch)")))+
  geom_boxplot(aes(x=Grazing,y=value,group=as.factor(Grazing)),outlier.shape = NA)+
  geom_jitter(aes(x=Grazing,y=value,color=as.factor(Grazing)),size=.5,alpha=.5)+
  scale_color_manual(values=c("grey","#FBD2A5","#FF8888","#C17F9D"))+
  facet_wrap(.~variable,scales = "free")+
  the_theme2+
  labs(x="Grazing intensity",y="")+
  theme(legend.position = "none")

ggsave("../Figures/SI/Distribution_spatial_stat_grazing_control_cover.pdf",p,width = 11,height = 7)


## >> Spatial metrics and aridity ----

d_data=read.table("../Data/Spatial_structure_grazing.csv",sep=";")%>%
  Closer_to_normality(.)

p=ggplot(d_data%>%
           melt(., measure.vars = c("PL_expo","fmax_psd","flow_length",
                                    "perim_area_scaling","core_area",
                                    "core_area_land","mean_psd"))%>%
           mutate(., variable=recode_factor(variable,
                                            "core_area_land"="% landscape covered by core pixels",
                                            "core_area"="Mean % of core pixels in patches",
                                            "flow_length"="Bare soil connectivity",
                                            "perim_area_scaling"="Fractal scaling area, perim.",
                                            "PL_expo"="Power-law exp. of the PSD",
                                            "fmax_psd"="log (largest patch)")))+
  geom_point(aes(x=Aridity,y=value),size=2,alpha=.5)+
  facet_wrap(.~variable,scales = "free")+
  the_theme2+
  labs(x="Aridity",y="")+
  theme(legend.position = "none")

ggsave("../Figures/SI/Distribution_spatial_stat_aridity.pdf",p,width = 11,height = 7)

## >> Importance of grazing and other variables without vegetation cover ----

#Importance using multimodel selection 

with_cov=F
d_all=read.table(paste0("../Data/Linear_models_factor/Importance_aridity_no_inter_factor.csv"),sep=";")%>%
  dplyr::rename(., Stat=Sp_stat)%>%
  Filter_relevant_stats(.)%>%
  filter(., Grazing_intensity=="all",With_cover==with_cov)%>%
  dplyr::select(., -Interactions,-Grazing_intensity,-With_cover,-Woody)



mean_imp=as_tibble(t(colMeans(d_all[,-c(1,2,ncol(d_all))])))

d_prep=rbind(cbind("Stat"="Averaged across \n spatial statistics",
                   (mean_imp-mean_imp)), #to have white box
             d_all[,-c(2,ncol(d_all))])%>%
  
  # add_column(., "R2"=0)%>%
  # dplyr::relocate(., R2, .after = Grazing)%>%
  melt(., measure.vars=colnames(.)[4:(ncol(.))])%>%
  add_column(., Type_stat=sapply(1:nrow(.),function(x){
    if (.$Stat[x] %in% c("Power-law exp. of the PSD","log (largest patch)",
                         "Mean patch size","Number of smallest patches")){
      return("Patch-size")
    }else if (.$Stat[x] =="Bare soil connectivity"){
      return("Hydro.")
    }else if (.$Stat[x] =="Spatial autocorrelation of vege."){
      return("Autocorr.")
    }else{
      return("Geom.")
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
  # geom_text(data=tibble(x=n_var,y=1:(n_metric-1),label=(round(unique(d_prep$R2m)[-1],2))),
  #           aes(x=x,y=y,label=paste0(label)),size=3)+
  geom_text(data=tibble(x=1:(n_var),y=(n_metric),
                        label=round(sort(as.numeric(mean_imp)[-c(1,2)],decreasing = T),2)),
            aes(x=x,y=y,label=paste0(label)),size=3)+
  scale_fill_gradientn(colours = colorRampPalette(c("white","grey","grey30"))(100))+
  the_theme+
  theme(axis.text.x = element_text(angle=60,hjust=1))+
  labs(x="Predictors",y="Spatial statistic",fill="Frequency in \n best models")

ggsave("../Figures/SI/Importance_grazing_other_variables_without_cover.pdf",
       p,width = 6,height = 5)


## >> Standardize predictors ----

d_all=read.table(paste0("../Data/Linear_models_factor_cover_control/Estimators_model_aridity_no_inter_factor.csv"),sep=";")%>%
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

loc_pval=1.5


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
ggsave(paste0("../Figures/SI/Standardize_coef_a.pdf"),p_tot,width = 12,height = 4)
p_tot=ggarrange(p_4,p_5,p_6,ncol=3,nrow=1)
ggsave(paste0("../Figures/SI/Standardize_coef_b.pdf"),p_tot,width = 12,height = 4)


#Same but without controlling for vegetation cover

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

loc_pval=1.5


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
ggsave(paste0("../Figures/SI/Standardize_coef_a_no_cover.pdf"),p_tot,width = 12,height = 4)
p_tot=ggarrange(p_4,p_5,p_6,ncol=3,nrow=1)
ggsave(paste0("../Figures/SI/Standardize_coef_b_no_cover.pdf"),p_tot,width = 12,height = 4)


## >> Spatial resolution: correlation, AIC and partial residuals ----

#Correlation spatial resolution and spatial metrics

d_data=read.table("../Data/Spatial_structure_grazing.csv",sep=";")

p=ggplot(d_data%>%melt(., measure.vars=c("PL_expo","fmax_psd","flow_length",
                                         "moran_I","core_area",
                                         "mean_psd"))%>%
            dplyr::rename(., Stat=variable)%>%
            Filter_relevant_stats(.)%>%
            Rename_spatial_statistics(.))+
  geom_point(aes(x=Resolution,y=value),color="grey")+
  the_theme2+
  facet_wrap(.~Stat,scales = "free",nrow = 2)+
  labs(x="Vegetation cover",y="Spatial metric")


ggsave("../Figures/SI/Correlation_metrics_Resolution.pdf",p,width = 9,height = 5)


#To evaluate the importance of spatial resolution, we compare the best models selected before
#with and without spatial resolution as a predictor

d_AIC=tibble()
list_mod=expand.grid(with_cover=c(F,T),
                     Stats=c("mean_psd","PL_expo","fmax_psd","flow_length",
                             "perim_area_scaling","core_area",
                             "core_area_land","KS_dist","Shape_metric")) 

for (id in 1:nrow(list_mod)){  
  
  stat=list_mod$Stats[id]
  with_cover=list_mod$with_cover[id]
  
  d_data_out=  read.table(paste0("../Data/Linear_models_factor_cover_control/Keep_data/Data_",stat,"_",with_cover,"_aridity_no_inter_all.csv"),sep=";") 
  if (ncol(d_data_out)==1){
    d_data_out=  read.table(paste0("../Data/Linear_models_factor_cover_control/Keep_data/Data_",stat,"_",with_cover,"_aridity_no_inter_all.csv"),sep=" ")
  }
  model_spa_stat=readRDS(paste0("../Data/Linear_models_factor_cover_control/Keep_models/Mod_",stat,"_",with_cover,"_aridity_no_inter_all.rds"))
  
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




d_slope=read.table("../Data/Linear_models_resolution/Slope_partial_residuals_aridity_no_inter.csv",sep=";")%>%
  filter(., With_cover==0,Driver=="Grazing",Grazing_intensity=="all")%>%
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
    if (.$Stat[x] %in% c("Power-law exp. of the PSD","log (largest patch)",
                         "Mean patch size","Number of smallest patches")){
      return("Patch-size")
    }else if (.$Stat[x] =="Bare soil connectivity"){
      return("Hydro.")
    }else if (.$Stat[x] =="Spatial autocorrelation of vege."){
      return("Autocorr.")
    }else{
      return("Geom.")
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
  facet_grid(Type_stat ~ ., space = "free_y", 
             scale = "free_y")+
  theme(panel.background = element_rect(fill="white"),
        legend.text = element_text(size=8))+
  guides(color=guide_legend(ncol=1))+
  theme(legend.position="bottom")#c(.2,.65))

ggsave("../Figures/SI/Partial_residuals_with_spatial_resolution.pdf",p,width = 5,height = 5)


## >> Correlation metrics-cover ----

d_data=read.table("../Data/Spatial_structure_grazing.csv",sep=";")

p1=ggplot(d_data%>%melt(., measure.vars=c("PL_expo","fmax_psd","flow_length",
                                          "moran_I","core_area",
                                          "Small_patches","mean_psd"))%>%
            dplyr::rename(., Stat=variable)%>%
            Filter_relevant_stats(.)%>%
            Rename_spatial_statistics(.))+
  geom_point(aes(x=rho_p,y=value),color="grey")+
  the_theme2+
  facet_wrap(.~Stat,scales = "free",nrow = 2)+
  labs(x="Vegetation cover",y="Spatial metric")


d_data=read.table("../Data/Spatial_structure_grazing_control_cover.csv",sep=";")

p2=ggplot(d_data%>%melt(., measure.vars=c("PL_expo","fmax_psd","flow_length",
                                          "moran_I","core_area",
                                          "Small_patches","mean_psd"))%>%
            dplyr::rename(., Stat=variable)%>%
            Filter_relevant_stats(.)%>%
            Rename_spatial_statistics(.))+
  geom_point(aes(x=rho_p,y=value),color="grey")+
  the_theme2+
  facet_wrap(.~Stat,scales = "free",nrow = 2)+
  labs(x="Vegetation cover",y="Spatial metric")


ggsave("../Figures/SI/Correlation_metrics_cover.pdf",
       ggarrange(p1+ggtitle("Not controlled for vegetation cover"),
                 p2+ggtitle("Controlled for vegetation cover"),
                 nrow = 2,
                 labels = letters[1:2]),width = 10,height = 12)


## >> Box-plot nutrients and grazing ----
# 
# d_data=read.table("../Data/Spatial_structure_grazing.csv",sep=";")%>%
#   Closer_to_normality(.)%>%
#   melt(., measure.vars = c("Org_C","Org_C_v","Total_N","lnTotal_N"))%>%
#   mutate(., variable=recode_factor(variable,
#                                    "lnTotal_N"="Total nitrogen (plants-bare)",
#                                    "Total_N"="Total nitrogen",
#                                    "Org_C"="Organic carbon (plants-bare)",
#                                    "Org_C_v"="Organic carbon",
#                                    "lnNitrate"="Nitrate (plants-bare)"))
# 
# p=ggplot(d_data)+
#   geom_boxplot(aes(x=Grazing,y=value,group=Grazing),outlier.shape = NA)+
#   geom_jitter(aes(x=Grazing,y=value,color=Grazing>1),size=.5,alpha=.5,height = 0)+
#   scale_color_manual(values=c("#92D0C5","#B685CA"))+
#   facet_wrap(.~variable,scales = "free")+
#   the_theme2+
#   labs(x="Grazing intensity",y="")+
#   theme(legend.position = "none")
# 
# ggsave("../Figures/SI/Distribution_N_C_cycling_grazing.pdf",p,width = 6,height = 5)


## >> Moving average aridity ----
 
# d_moving=read.table("../Data/Moving_average_aridity/Moving_average_aridity.csv",sep=";")%>%
#   filter(., With_cover==0)%>%
#   Filter_relevant_stats(.)%>%  
#   Rename_spatial_statistics(.)
# 
# p=ggplot(d_moving)+
#   geom_line(aes(x=Aridity_thresh,y=q2),lwd=1)+
#   geom_ribbon(aes(x=Aridity_thresh,y=q2,ymin=q1,ymax=q3),alpha=.6,fill="grey")+
#   the_theme2+
#   geom_hline(yintercept = 0,color="#9CD891")+
#   facet_wrap(.~Stat,scales = "free")+
#   labs(x="Aridity",y=substitute(paste(beta," (partial res. aridity)")))
# 
# 
# ggsave("../Figures/SI/Moving_average_aridity_sand.pdf",p,width = 7,height = 5)





## >> Partial residuals without accounting for vegetation cover ----

d_slope=read.table("../Data/Linear_models_factor/Slope_partial_residuals_aridity_no_inter.csv",sep=";")%>%
  filter(., With_cover==0,Driver=="Grazing",Grazing_intensity=="all")%>%
  filter(., Stat %in% c("PL_expo","fmax_psd","flow_length",
                        "moran_I","core_area","Small_patches",
                        "mean_psd"))%>%
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
    if (.$Stat[x] %in% c("Power-law exp. of the PSD","log (largest patch)",
                         "Mean patch size","Number of smallest patches")){
      return("Patch-size")
    }else if (.$Stat[x] =="Bare soil connectivity"){
      return("Hydro.")
    }else if (.$Stat[x] =="Spatial autocorrelation of vege."){
      return("Autocorr.")
    }else{
      return("Geom.")
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
  facet_grid(Type_stat ~ ., space = "free_y", 
             scale = "free_y")+
  theme(panel.background = element_rect(fill="white"),
        legend.text = element_text(size=8))+
  guides(color=guide_legend(ncol=1))+
  theme(legend.position="bottom")#c(.2,.65))

ggsave("../Figures/SI/Grazing_intensity_partial_residuals_without_cover.pdf",p,width = 5,height = 5.5)



## >> Aridity vs grazing spatial structure ----

d_slope=read.table("../Data/Linear_models_factor_cover_control/Slope_partial_residuals_aridity_no_inter.csv",sep=";")%>%
  filter(., With_cover==0,Grazing_intensity=="binary")%>%
  Filter_relevant_stats(.)%>%  
  Rename_spatial_statistics(.)%>%
  mutate(., ID_grazing=as.character(ID_grazing))%>%
  arrange(., q2)%>%
  add_column(., Order_f=1:nrow(.))%>%
  mutate(.,Stat = fct_reorder(Stat, Order_f))%>%
  add_column(., Signif=.$pval<.05)%>%
  arrange(., Stat,q2)%>%
  add_column(., Type_stat=sapply(1:nrow(.),function(x){
    if (.$Stat[x] %in% c("Power-law exp. of the PSD","log (largest patch)",
                         "Mean patch size","Number of smallest patches")){
      return("Patch-size")
    }else if (.$Stat[x] =="Bare soil connectivity"){
      return("Hydro.")
    }else if (.$Stat[x] =="Spatial autocorrelation of vege."){
      return("Autocorr.")
    }else{
      return("Geom.")
    }
  }))

p=ggplot(d_slope)+
  geom_vline(xintercept = 0,linetype=9)+
  
  geom_linerange(aes(x=q2,y=Stat,xmin=q1_90,xmax=q3_90,color=as.factor(Driver)),
                 lwd=2,position=position_jitterdodge(seed=123))+   
  geom_linerange(aes(x=q2,y=Stat,xmin=q1,xmax=q3,color=as.factor(Driver)),
                 lwd=.8,position=position_jitterdodge(seed=123))+  
  geom_pointrange(aes(x=q2,y=Stat,xmin=q1,xmax=q3,fill=as.factor(Driver)),color="transparent",
                  size=.3,position=position_jitterdodge(seed=123),shape=21)+  
  
  the_theme2+
  scale_color_manual(values=c("orange","grey50"))+
  scale_fill_manual(values=c("grey20","grey20","grey20"))+
  scale_shape_manual(values=c(21,16))+
  labs(x=substitute(paste(beta," (partial res. Grazing)")),y="",color="")+
  guides(shape="none",fill="none")+
  facet_grid(Type_stat ~ ., space = "free_y", 
             scale = "free_y")+
  theme(panel.background = element_rect(fill="white"),
        legend.text = element_text(size=8))+
  guides(color=guide_legend(ncol=1))+
  theme(legend.position="bottom")#c(.2,.65))

ggsave("../Figures/SI/Grazing_intensity_partial_residuals_aridity_grazing.pdf",p,width = 5,height = 5)



## >> Toy examples Schneider model PSD ----

list_landscape=readRDS("../Data/Simulations/Minimal_examples_stats_landscapes.rds")
d=read.table("../Data/Simulations/Minimal_examples_stats.csv",sep=";")


id=1
for (cover_id in 1:3){
  for (aggregation_id in 1:4){
    assign(paste0("p_psd",aggregation_id,cover_id),
           Plot_psd_raw(list_landscape[[id]])+
             scale_y_log10(limits=c(1,1000)))
    id=id+1
  }
}


p_psd=ggarrange(
  p_psd11+ggtitle("Cover = 0.65, \nHigh aggregation")+
    theme(panel.border = element_rect(colour = "#D27C3F",fill=NA,size=2)),
  p_psd12+ggtitle("Cover = 0.65, \nMedium aggregation")+
    theme(panel.border = element_rect(colour = "#D27C3F",fill=NA,size=2)),
  p_psd13+ggtitle("Cover = 0.65, \nLow aggregation")+
    theme(panel.border = element_rect(colour = "#D27C3F",fill=NA,size=2)),
  # p_psd41+ggtitle("Cover = 0.65, n random"),
  p_psd21+ggtitle("Cover = 0.35, \nHigh aggregation")+
    theme(panel.border = element_rect(colour = "#8BD882",fill=NA,size=2)),
  p_psd22+ggtitle("Cover = 0.35, \nMedium aggregation")+
    theme(panel.border = element_rect(colour = "#8BD882",fill=NA,size=2)),
  p_psd23+ggtitle("Cover = 0.35, \nLow aggregation")+
    theme(panel.border = element_rect(colour = "#8BD882",fill=NA,size=2)),
  # p_psd42+ggtitle("Cover = 0.35, n random"),
  p_psd31+ggtitle("Cover = 0.17, \nHigh aggregation")+
    theme(panel.border = element_rect(colour = "#BE91D0",fill=NA,size=2)),
  p_psd32+ggtitle("Cover = 0.17, \nMedium aggregation")+
    theme(panel.border = element_rect(colour = "#BE91D0",fill=NA,size=2)),
  p_psd33+ggtitle("Cover = 0.17, \nLow aggregation")+
    theme(panel.border = element_rect(colour = "#BE91D0",fill=NA,size=2)),
  # p_psd43+ggtitle("Cover = 0.17, \n random"),
  nrow=3,ncol=3)

ggsave("../Figures/SI/Examples_stat_toy_PSD.pdf",p_psd,width = 9,height = 8)

## >> Interaction aridity and grazing ----

d=tibble()

name_stats=c("Power-law exp. \n of the PSD","log (largest patch)",
             "Bare soil connectivity","Spatial autocorrelation \n of vege.",
             "Number of \n smallest patches",
             "Mean % of core \n pixels in patches",
             "Mean patch size")

color_grazing=c("grey","#FBD2A5","#FF8888","#C17F9D")

id=1

for (k in c("PL_expo","fmax_psd",
            "flow_length","moran_I","Small_patches",
            "core_area","mean_psd")){
  
  d_data_out=read.table(paste0("../Data/Linear_models_factor_cover_control/Keep_data/Data_",k,"_FALSE_aridity_all.csv"),sep=" ")
  
  model_spa_stat=readRDS(paste0("../Data/Linear_models_factor_cover_control/Keep_models/Mod_",k,"_FALSE_aridity_all.rds"))
  
  d_data_out$Grazing = as.factor(d_data_out$Grazing)
  d_data_out$Grazing = relevel(d_data_out$Grazing,ref = "0")
  
  partial_res=visreg::visreg(fit = model_spa_stat,xvar="Aridity",by="Grazing",plot=F)
  
  for (grazing_id in unique(partial_res$res$Grazing)){
    summary_lm=summary(lm(data = filter(partial_res$res,Grazing==grazing_id),visregRes~Aridity))
    
    d=rbind(d,tibble(Stat=k,
                     pval=summary_lm$coefficients[2,4],
                     Grazing_id=grazing_id))
    
    if (summary_lm$coefficients[2,4]<.05){
      
      assign(paste0("p_",grazing_id),
             ggplot(filter(partial_res$res,Grazing==grazing_id))+
               geom_point(aes(x=Aridity,y=visregRes),fill=color_grazing[as.numeric(grazing_id)+1],color="black",shape=21)+
               geom_smooth(aes(x=Aridity,y=visregRes),method = "lm",color="black",fill=color_grazing[as.numeric(grazing_id)+1])+
               the_theme2+
               labs(x="Aridity",y=name_stats[id]))
    } else{
      assign(paste0("p_",grazing_id),
             ggplot(filter(partial_res$res,Grazing==grazing_id))+
               geom_point(aes(x=Aridity,y=visregRes),fill=color_grazing[as.numeric(grazing_id)+1],color="black",shape=21)+
               the_theme2+
               labs(x="Aridity",y=name_stats[id]))
    }
    
  }
  
  if (k=="PL_expo"){
    p_0=p_0+ggtitle("Ungrazed")
    p_1=p_1+ggtitle("Low grazing")
    p_2=p_2+ggtitle("Medium grazing")
    p_3=p_3+ggtitle("High grazing")
  }
  if (k!="mean_psd"){
    p_0=p_0+theme(axis.title.x = element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank())
    p_1=p_1+theme(axis.title.x = element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank())
    p_2=p_2+theme(axis.title.x = element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank())
    p_3=p_3+theme(axis.title.x = element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank())
  }
  
  p_1=p_1+theme(axis.title.y = element_blank())
  p_2=p_2+theme(axis.title.y = element_blank())
  p_3=p_3+theme(axis.title.y = element_blank())
  assign(paste0("p_tot_",id),
         ggarrange(p_0,p_1,p_2,p_3,ncol=4))
  
  id=id+1
}

p_fig=ggarrange(p_tot_1,
                p_tot_2,
                p_tot_3,
                p_tot_4,
                p_tot_5,
                p_tot_6,
                p_tot_7,
                nrow=7,align = "hv",
                labels=letters[1:7])

ggsave("../Figures/SI/Interactive_effects.pdf",p_fig,width = 9,height = 14)

## >> Correlation between traits (pair plot and PCA) ----

d_traits_norm=d_traits%>%Closer_to_normality_traits(.)
corr_pred=corr.test(d_traits_norm[,7:(ncol(d_traits_norm)-1)],use = "na.or.complete",adjust = "none")
corr_pred$r=round(corr_pred$r,2)
corr_pred$r[lower.tri(corr_pred$r)]=NA
diag(corr_pred$r)=NA
p=ggplot(corr_pred$r%>%
           melt(.)%>%
           add_column(., pval=sapply(1:nrow(.), function(x){
             if (is.na(.$value[x])){
               return(NA)
             }else{return(melt(corr_pred$p)$value[x])}
           })))+
  geom_tile(aes(x=Var1,Var2,fill=value))+
  the_theme+
  geom_text(aes(x=Var1,Var2,label=ifelse(pval<.1,"","X")))+
  scale_fill_gradient2(low="red",mid="white",high = "blue",midpoint = 0,na.value = "white")+
  theme(axis.text.x = element_text(angle=60,hjust=1))+
  labs(x="",y="",fill="")
ggsave("../Figures/SI/Correlation_between_traits.pdf",p,width = 6,height = 5)





d_traits_norm=d_traits%>%Closer_to_normality_traits(.)

struct_variables=colnames(d_traits_norm)[7:(ncol(d_traits_norm)-1)]
res.comp=imputePCA(d_traits_norm[,which(colnames(d_traits_norm) %in% struct_variables)],ncp=4,scale = T) 

if ("completeObs" %in% names(res.comp)){
  res.pca=PCA(res.comp$completeObs, ncp = 4,  graph=F)
}else {
  res.pca=PCA(res.comp, ncp = 4,  graph=F)
}

p1=factoextra::fviz_pca_var(res.pca,col.var="#2746B1")+
  the_theme+
  ggtitle("")+
  labs(x=paste0("PC 1 (",round(res.pca$eig[1,2],1),")"),
       y=paste0("PC 2 (",round(res.pca$eig[2,2],1),")"))

p2=factoextra::fviz_pca_var(res.pca,col.var="#2746B1",axes = c(1,3))+
  the_theme+
  ggtitle("")+
  labs(x=paste0("PC 1 (",round(res.pca$eig[1,2],1),")"),
       y=paste0("PC 2 (",round(res.pca$eig[3,2],1),")"))

p3=factoextra::fviz_pca_var(res.pca,col.var="#2746B1",axes = c(2,3))+
  the_theme+
  ggtitle("")+
  labs(x=paste0("PC 1 (",round(res.pca$eig[2,2],1),")"),
       y=paste0("PC 2 (",round(res.pca$eig[3,2],1),")"))


ggsave("../Figures/SI/PCA_on_traits.pdf",ggarrange(p1,p2,p3,ncol=3,labels = letters[1:3]),width = 12,height = 5)
