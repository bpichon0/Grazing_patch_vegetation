rm(list=ls())
source("./Structure_grazing_function.R")



# -------------------- Main figures -------------------------

## >> Figure 1: Importance and effect of grazing with cover ----

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
  geom_text(data=d_R2,aes(y=Stat,x=20,label=paste0("R2m = ",round(R2m,2),", ","R2c = ",round(R2C,2))),size=3)+
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

list_resid=list(stat=c("perim_area_scaling", "core_area",
                       "PL_expo","core_area"),
                metric=c("Type_veg","Type_veg",
                         "Clim2","Clim1"),
                title=c("Fractal scaling area, perim.","Mean % of core in patches",
                        "Power-law exp. of the PSD","Mean % of core in patches"),
                x_axis=c("Grazing intensity","Grazing intensity",
                         "PC2 climatic variables", 'PC1 climatic variables'),
                xwrite=c(-3,-1,-1),ywrite=c(-3,-1.5,-2))

for (ID in 1:length(list_resid$stat)){
  
  stat=list_resid$stat[ID]
  metric=list_resid$metric[ID]
  
  d_data_out=  read.table(paste0("../Data/Step1_Understanding_grazing/Keep_data/Data_",stat,"_TRUE.csv"),sep=";")
  if (ncol(d_data_out)==1){
    d_data_out=  read.table(paste0("../Data/Step1_Understanding_grazing/Keep_data/Data_",stat,"_TRUE.csv"),sep=" ")
  }
  model_spa_stat=readRDS(paste0("../Data/Step1_Understanding_grazing/Keep_models/Mod_",stat,"_TRUE.rds"))
  
  
  if (metric=="Type_veg"){
    
    resid_mod=visreg::visreg(fit = model_spa_stat,xvar = "Grazing",by="Type_veg",plot=F) #extracting residuals
    
    summary(lm(visregRes~Grazing,resid_mod$res%>%filter(.,Type_veg=="Grass_Shrub")))
    summary(lm(visregRes~Grazing,resid_mod$res%>%filter(.,Type_veg=="Forest")))
    summary(lm(visregRes~Grazing,resid_mod$res%>%filter(.,Type_veg=="Grassland")))
    summary(lm(visregRes~Grazing,resid_mod$res%>%filter(.,Type_veg=="Shrubland")))
    
    
    if (ID==1){
      annot_mod=annotate("text", x = 2.5, y = ifelse(stat=="core_area",-1.1,-1.5), 
               label = "P slopes < 0.05")
        
    }else {
      annot_mod=annotate("text", x = 2.5, y = ifelse(stat=="core_area",-1.1,-1.5), 
                         label = "Shrubland : NS")
    }
      
    assign(paste0("p_",ID),ggplot(resid_mod$res%>%
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
      labs(x=list_resid$x_axis[ID],y=list_resid$title[ID],fill="",color=""))

  }else {
    
    resid_mod=visreg::visreg(fit = model_spa_stat,xvar = metric,by="Grazing",plot=F) #extracting residuals
    resid_mod$res$Grazing=rep(0:3,as.vector(table(d_data_out$Grazing))) #correcting vigreg output
    
    summary(lm(formula(paste("visregRes~Grazing+",metric,"+",metric,"*Grazing")),resid_mod$res))
    
    if (ID==3){
      annot_mod=annotate("text", x = 1, y = -1.45, parse = TRUE,
               label = "'Y = ' * alpha[1] * 'G' * + alpha[2] * 'C2***'* + alpha[3] *  'C2 x G***'")
    }else {
      annot_mod=annotate("text", x = 0.25, y = -1.65, parse = TRUE,
                         label = "'Y = ' * alpha[1] * 'G' * + alpha[2] * 'C1'* + alpha[3] *  'C1 x G***'")
    }
    
    assign(paste0("p_",ID),
           ggplot(resid_mod$res%>%
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
             labs(x=list_resid$x_axis[ID],y=list_resid$title[ID],fill="",color=""))
    }
}


p=ggarrange(ggarrange(p_1,p_2,ncol=2,labels=letters[1:2],
                      common.legend = T,legend = "bottom"),
          ggarrange(p_3,p_4,ncol=2,labels = letters[3:4],
                    common.legend = T,legend = "bottom"),
          nrow=2,heights = c(1,1))
ggsave("../Figures/Final_figs/Interaction_with_grazing.pdf",p,width = 8,height = 7)



## >> Figure 3: Indirect effects of grazing ----

d_data=read.table("../Data/Spatial_structure_grazing.csv",sep=";")%>%
  Closer_to_normality(.)

d_data[,c(1,9:29,35,37:41,44:ncol(d_data))] = apply(d_data[,c(1,9:29,35,37:41,44:ncol(d_data))],2,z_tranform)

#Getting the multifunctionality index
d_data$MF=z_tranform(rowSums(d_data[,55:(ncol(d_data))],na.rm = T))

#Adding the first two components from the PCA
d_data=Perform_PCA_spatial_struc(d_data)


param_list=tibble(Stat=c("Struct1","Struct2",
                         "perim_area_scaling","PL_expo","fmax_psd",
                         "PLR","flow_length","core_area_land",
                         "contig"),
                  Name_var=c("PC1 sp. struct.","PC2 sp. struct.",
                             "Fractal scaling area, perim.",
                             "Power-law exp. of the PSD",
                             "log (largest patch)",
                             "PLR","Bare soil connectivity",
                             "% landscape covered by core",
                             "Contiguity"))
d_indirect=tibble()

id_plot=1
pdf(paste0("../Figures/Final_figs/SI/SEM_grazing_",letters[id_plot],".pdf"),
    width = 6,height = 4.5)
par(mfrow=c(2,2))

for (ID in 1:nrow(param_list)){
  
  k=as.character(param_list$Stat[ID])
  
  model_lmer=lm("value ~ Long_cos + Long_sin + Lattitude + Slope + Elevation +
                        Clim1 + Clim2 + Clim3 + Clim4 + Sand + Type_veg",
                data = d_data%>%melt(., measure.vars=k)%>%filter(., !is.na(value)),
                na.action = na.omit)
  
  #we remove potential outliers
  rm.outliers = romr.fnc(model_lmer, d_data%>%melt(., measure.vars=k)%>%filter(., !is.na(value)),
                         trim=2.5)
  d_data_out = rm.outliers$data
  
  #We first control for all covariates and extract the residuals
  model_lmer=lm("value ~ Long_cos + Long_sin + Lattitude + Slope + Elevation +
                        Clim1 + Clim2 + Clim3 + Clim4 + Sand + Type_veg",
                data = d_data_out,
                na.action = na.fail)
  
  resid_model=residuals(model_lmer) #extract residuals
  
  save=d_data[as.numeric(names(resid_model)),]%>% #add it to the dataframe 
    add_column(., Resid_mod=resid_model)
  
  
  #DOING the SEMs
  
  d_sem=save
  # all_d=summary(psem(
  #   lmer(Resid_mod ~ (1|Site_ID) + Grazing + Woody + rho_p + Org_C, d_sem),
  #   lmer(Woody ~ (1|Site_ID) + Grazing + Org_C , d_sem),
  #   lmer(rho_p ~ (1|Site_ID) + Grazing + Org_C   , d_sem),
  #   lmer(Org_C ~ (1|Site_ID) + Grazing, d_sem)
  # ))
  
  #SEM with low grazing intensity
  
  d_sem=filter(save,Grazing %in% 0:1)
  # low_graz=summary(psem(
  #   lmer(Resid_mod ~ (1|Site_ID) + Grazing + Woody + rho_p + Org_C, d_sem),
  #   lmer(Woody ~ (1|Site_ID) + Grazing + Org_C , d_sem),
  #   lmer(rho_p ~ (1|Site_ID) + Grazing + Org_C   , d_sem),
  #   lmer(Org_C ~ (1|Site_ID) + Grazing, d_sem)
  # ))
  
  low_graz=summary(psem(
    lmer(Resid_mod ~ (1|Site_ID) + Grazing + Woody +
           rho_p + Org_C + Sand + lnTotal_N +
           lnTotal_P , d_sem),
    lmer(Woody ~ (1|Site_ID) + Grazing + Org_C , d_sem),
    lmer(rho_p ~ (1|Site_ID) + Grazing + Org_C  , d_sem),
    lmer(Sand ~ (1|Site_ID) + Grazing   , d_sem),
    lmer(lnTotal_P ~ (1|Site_ID) + Grazing   , d_sem),
    lmer(lnTotal_N ~ (1|Site_ID) + Grazing   , d_sem),
    lmer(Org_C ~ (1|Site_ID) + Grazing, d_sem)
  ))
  
  
  print(low_graz)
  
  
  #SEM with high grazing intensity
  
  d_sem=filter(save,Grazing %in% 2:3)
  # high_graz=summary(psem(
  #   lmer(Resid_mod ~ (1|Site_ID) + Grazing + Woody + rho_p + Org_C, d_sem),
  #   lmer(Woody ~ (1|Site_ID) + Grazing + Org_C , d_sem),
  #   lmer(rho_p ~ (1|Site_ID) + Grazing + Org_C  , d_sem),
  #   lmer(Org_C ~ (1|Site_ID) + Grazing, d_sem)
  # ))
  high_graz=summary(psem(
    lmer(Resid_mod ~ (1|Site_ID) + Grazing + Woody +
           rho_p + Org_C + Sand + lnTotal_N +
           lnTotal_P , d_sem),
    lmer(Woody ~ (1|Site_ID) + Grazing + Org_C , d_sem),
    lmer(rho_p ~ (1|Site_ID) + Grazing + Org_C  , d_sem),
    lmer(Sand ~ (1|Site_ID) + Grazing   , d_sem),
    lmer(lnTotal_P ~ (1|Site_ID) + Grazing   , d_sem),
    lmer(lnTotal_N ~ (1|Site_ID) + Grazing   , d_sem),
    lmer(Org_C ~ (1|Site_ID) + Grazing, d_sem)
  ))
  
  print(high_graz)
  
  # Plot_SEM_with_cover(summary_sem = low_graz,pdf_ = F,
  #                     name = "",
  #                     title_ = "Low grazing pressure",
  #                     name_var = param_list$Name_var[ID],MF = param_list$MF[ID])
  # 
  # Plot_SEM_with_cover(summary_sem = high_graz,pdf_ = F,
  #                     name = "",
  #                     title_ = "High grazing pressure",
  #                     name_var = param_list$Name_var[ID],MF = param_list$MF[ID])
  
  #Take all indirect effects 
  d_indirect=rbind(d_indirect,
                   Get_indirect_effects_grazing(all_d)%>%
                     add_column(., Stat=param_list$Name_var[ID],Intensity="All"),
                   Get_indirect_effects_grazing(low_graz)%>%
                     add_column(., Stat=param_list$Name_var[ID],Intensity="Low"),
                   Get_indirect_effects_grazing(high_graz)%>%
                     add_column(., Stat=param_list$Name_var[ID],Intensity="High"))
  
  
  #new pdf
  # if (ID %in% c(2,4,6,8)){
  #   id_plot=id_plot+1
  #   dev.off()
  #   if (letters[id_plot]=="e"){
  #     pdf(paste0("../Figures/Final_figs/SI/SEM_grazing_",letters[id_plot],".pdf"),
  #         width = 7,height = 2.5)
  #     
  #     par(mfrow=c(1,2))
  #     
  #   }else{
  #     pdf(paste0("../Figures/Final_figs/SI/SEM_grazing_",letters[id_plot],".pdf"),
  #         width = 6,height = 4.5)
  #     par(mfrow=c(2,2))
  #   }
  # }
  
}
dev.off()

p=ggplot(d_indirect%>%filter(., Intensity!="All")%>%
         dplyr::group_by(., Stat,Intensity)%>%
         dplyr::summarise(., .groups = "keep",
                          Effect=sum(Effect),Path=Path),
       aes(y=Stat,x=Effect,fill=Intensity))+
  geom_bar(position="dodge",  stat="identity",alpha=.6, width=.4) +
  the_theme+
  labs(x="Total effect on the spatial structure",y="",fill="")+
  geom_vline(xintercept = 0)+
  scale_fill_manual(values=c("#B7151C","#E6A0B8"))

ggsave("../Figures/Final_figs/Total_effect_sp_struct.pdf",p,width = 5,height = 4)

## >> Figure 4: Predictors of dryland stability ----

d_mod=read.table("../Data/Inferrence/Drivers_stability_metrics_no_cover.csv",sep=";")


for (k in 1:3){
  assign(paste0("p_",k),
         ggplot(d_mod%>%melt(., id.vars=c("Param","Type_grazing"))%>%
                  filter(., Param %!in% c("p","Absolute distance"))%>%
                  mutate(., Param=recode_factor(Param,"p"="Parameter p, (Cover)",
                                                "q" ="Parameter q, (spatial structure)",
                                                "Relative distance"="Distance to desertification point",
                                                "Size tipping" ="Size of the tipping point"))%>%
                  filter(., Param==unique(.$Param)[k])%>%
                  dplyr::group_by(., variable,Param,Type_grazing)%>%
                  dplyr::summarise(., .groups = "keep",q2=median(value),q1=quantile(value,.025),q3=quantile(value,.975))%>%
                  add_column(., color_pred=sapply(1:nrow(.),function(x){
                    if (.$variable[x] %in% c("Lat","Elevation","Slope","Long_cos","Long_sin")){return("COV")
                    }else{return("NO_COV")}
                  }))%>%filter(., color_pred=="NO_COV")%>%
                  mutate(., variable=recode_factor(variable,
                                                   "Sp_richness"="Species richness",
                                                   "Clim1"="PC1 climatic",
                                                   "Clim2"="PC2 climatic",
                                                   "Clim3"="PC3 climatic",
                                                   "Clim4"="PC4 climatic",
                                                   "Org_C"="Organic carbon",
                                                   "SR"="Species richness",
                                                   "Clim1_Grazing"="PC1 climatic X Grazing",
                                                   "Clim2_Grazing"="PC2 climatic X Grazing",
                                                   "Grazing_Org_C"="Organic carbon X Grazing")))+
           geom_pointrange(aes(y=variable,x=q2,xmin=q1,xmax=q3,color=color_pred))+
           the_theme+
           facet_grid(Type_grazing~Param,scales="free")+
           labs(y="",x="Standardized effect",fill="Fixed effect  ")+
           geom_vline(xintercept = 0,linetype=9)+
           scale_color_manual(values=c("grey","black"))+
           guides(fill="none")+
           theme(legend.position = "none"))
}

p_tot=ggarrange(p_1,
                p_2+theme(axis.text.y = element_blank(),
                          axis.title.y = element_blank(),
                          axis.ticks.y = element_blank()),
                p_3+theme(axis.text.y = element_blank(),
                          axis.title.y = element_blank(),
                          axis.ticks.y = element_blank()),
                ncol = 3,widths = c(1.5,1,1))

ggsave("../Figures/Final_figs/LME_stability_spatial_struct_no_cov.pdf",p_tot,width=15,height = 6)


d_partial=read.table("../Data/Inferrence/Partial_residuals_grazing_no_cover.csv",sep=";")

p=ggplot(NULL)+
  geom_flat_violin(data=d_partial,
                   aes(y=slope,x=Stat),fill="grey",position = position_nudge(x = 0.1, y = 0), 
                   alpha = 0.8,scale = "width",width=.3)+
  coord_flip()+
  geom_pointrange(data=d_partial%>%melt(., id.vars=c("Stat","Type"))%>%
                    dplyr::group_by(., Stat,Type)%>%
                    dplyr::summarise(., .groups = "keep",q2=median(value),q1=quantile(value,.05),q3=quantile(value,.95)),
                  aes(y=q2,ymin=q1,ymax=q3,x=Stat))+
  the_theme+
  labs(x="",y="Standardized effect",fill="Fixed effect  ")+
  geom_hline(yintercept = 0,linetype=9)+
  guides(fill="none")+
  facet_wrap(.~Type)+
  theme(legend.position = "none")

ggsave("../Figures/Final_figs/Stabilit_partial_effect_no_cover.pdf",p,width=6,height = 5)


  


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

ggsave("../Figures/Final_figs/SI/Criteria_selection_images.pdf",
       ggarrange(ggarrange(p1,p2+theme(legend.position = "none"),ncol = 2,labels = letters[1:2],widths = c(.8,1)),
                 ggarrange(ggplot()+theme_void(),get_legend(p2),ncol=3,widths = c(1,.7)),
                 nrow=2,heights = c(1,.1)),
       width = 7,height = 3.5)

## >> Distribution of co-variates across sites ----
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
#                     extra=c(R2m=function(x) r.squaredGLMM(x)[1],
#                             R2c=function(x) r.squaredGLMM(x)[2]),
#                     options(na.action = "na.fail") )

# saveRDS(select_model,"../Data/Step1_Understanding_grazing/Keep_models/Mod_cover.rds")
select_model=readRDS("../Data/Step1_Understanding_grazing/Keep_models/Mod_cover.rds")


#best models (delta AIC <2)
model_cover=model.avg(select_model, subset = delta < 2)



R2=select_model%>%filter(., AICc<min(AICc)+2)

summary_coef=confint(model_cover)

#Merge in a df
d_pred=tibble(Median=confint(model_cover,level = 1e-9)[-1,1],
              q1=summary_coef[-1,1],
              q3=summary_coef[-1,2],
              pvalue=summary(model_cover)$coefmat.full[-1,5],
              term=rownames(summary_coef)[-1],
              Stat="Cover",
              R2m=mean(R2$R2m),
              R2c=mean(R2$R2c))

loc_pval=.05+max(d_pred%>% #localization of pvalues in the plot
                   Organize_df(., "predictor")%>%
                   dplyr::select(.,q3)%>%dplyr::pull(.))

p1_1=ggplot(d_pred%>%
              filter(., term!="Type, both")%>%
              Organize_df(., "predictor")%>%
              add_column(., Is_signif_pval=sapply(1:nrow(.),function(x){return(is_signif(.$pvalue[x]))})))+
  geom_pointrange(aes(x=Median,y=term,xmin=q1,xmax=q3,color=Type_pred))+
  geom_text(aes(x=loc_pval,y=term,label=Is_signif_pval))+
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
  geom_text(aes(x=0,y=1.05,label=paste0("R? \n = \n ",unique(round(d_pred$R2m,2)))),size=3)+
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
model_cover=MuMIn::get.models(select_model, subset = delta == 0)[[1]]

resid_mod=visreg::visreg(fit = model_cover,xvar="Grazing",plot=F)

summary(lm(visregRes~Grazing,resid_mod$res))


p2=ggplot(NULL)+
  geom_jitter(data=resid_mod$res,
              aes(Grazing,visregRes),width = .1,color="gray",alpha=.5)+
  geom_pointrange(data=resid_mod$res%>%
                    dplyr::group_by(., Grazing)%>%
                    dplyr::summarise(., .groups = "keep",
                                     q1=quantile(visregRes,.25),q3=quantile(visregRes,.75),q2=median(visregRes)),
                  aes(x=Grazing,y=q2,ymin=q1,ymax=q3),color="black")+
  annotate("text", x = 1.5, y = 1, parse = TRUE,
           label = "'Y = ' * alpha[1] * ' G***'")+
  the_theme+
  scale_color_manual(values=c("0"="#FFC170",'1'="#43FFF4","2"="#5196F9","3"="#012E71"))+
  scale_fill_manual(values=c("0"="#FFC170",'1'="#43FFF4","2"="#5196F9","3"="#012E71"))+
  labs(x="Grazing",y="Residuals of cover")


#extracting residuals of cover against woody
resid_mod=visreg::visreg(fit = model_cover,xvar="Woody",plot=F)
summary(lm(visregRes~Woody,resid_mod$res))

p3=ggplot(resid_mod$res)+
  geom_point(aes(Woody,visregRes),color="gray30",alpha=.2)+
  geom_smooth(aes(Woody,visregRes),color="black",fill="gray",method = "lm",level=.95)+
  annotate("text", x = 0, y = 1, parse = TRUE,
           label = "'Y = ' * alpha[1] * ' W***'")+
  the_theme+
  scale_color_manual(values=c("0"="#FFC170",'1'="#43FFF4","2"="#5196F9","3"="#012E71"))+
  scale_fill_manual(values=c("0"="#FFC170",'1'="#43FFF4","2"="#5196F9","3"="#012E71"))+
  labs(x="% woody cover",y="Residuals of cover")

#extracting residuals of cover against grazing
resid_mod=visreg::visreg(fit = model_cover,xvar="Org_C",by="Grazing",plot=F)
resid_mod$res$Grazing=rep(0:3,as.vector(table(d_data_out$Grazing))) #correcting vigreg output

summary(lm(visregRes~Grazing+Org_C+Org_C*Grazing,resid_mod$res))


p4=ggplot(resid_mod$res%>%
            mutate(., Grazing=as.character(Grazing)))+
  geom_point(aes(Org_C,visregRes),color="gray30",alpha=.2)+
  geom_smooth(aes(Org_C,visregRes,fill=Grazing,color=Grazing),method = "lm",level=.95)+
  annotate("text", x = -1, y = 1.5, parse = TRUE,
           label = "'Y = ' * alpha[1] * ' W***' * + alpha[2] * 'F**'* + alpha[3] *  'F x G***'")+
  the_theme+
  scale_color_manual(values=c("0"="#FFC170",'1'="#43FFF4","2"="#5196F9","3"="#012E71"))+
  scale_fill_manual(values=c("0"="#FFC170",'1'="#43FFF4","2"="#5196F9","3"="#012E71"))+
  labs(x="Facilitation",y="Residuals on cover")

ggsave(paste0("../Figures/Final_figs/SI/Grazing_on_cover.pdf"),
       ggarrange(
         ggarrange(ggplot()+theme_void(),p1,ggplot()+theme_void(),ncol=3,widths = c(.25,1,.25)),
         ggarrange(
           ggarrange(p2,p3,p4+theme(legend.position="none"),ncol=3,labels = c(letters[3],"","")),
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
ggsave(paste0("../Figures/Final_figs/SI/Standardize_coef_a.pdf"),p_tot,width = 9,height = 9)
p_tot=ggarrange(p_5,p_6,p_7,p_8,ncol=2,nrow=2)
ggsave(paste0("../Figures/Final_figs/SI/Standardize_coef_b.pdf"),p_tot,width = 9,height = 9)



## >> Standardize predictors no cover ----

cover=0
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
ggsave(paste0("../Figures/Final_figs/SI/Standardize_coef_a_no_cov.pdf"),p_tot,width = 9,height = 9)
p_tot=ggarrange(p_5,p_6,p_7,p_8,ncol=2,nrow=2)
ggsave(paste0("../Figures/Final_figs/SI/Standardize_coef_b_no_cov.pdf"),p_tot,width = 9,height = 9)



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

## >> Correlation metrics-resolution ----
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

ggsave("../Figures/Final_figs/SI/Correlation_metrics_Resolution.pdf",p,width = 9,height = 4)

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
  geom_text(data=d_R2,aes(y=Stat,x=20,label=paste0("R?m = ",round(R2m,2),", ","R?c = ",round(R2C,2))),size=3)+
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

d_all=read.table(paste0("../Data/Step1_Understanding_grazing/Importance_no_inter.csv"),sep=";")%>%
  dplyr::rename(., Stat=Sp_stat)%>%
  Filter_relevant_stats(.)%>%
  dplyr::select(., -Interactions)


with_cov=T

mean_imp=as_tibble(t(colMeans(d_all[which(d_all$With_cover==with_cov),c(-1,-ncol(d_all))])))

d_prep=rbind(cbind("Stat"="Averaged across \n spatial statistics",
                   (mean_imp-mean_imp)%>% #to have white box
                     add_column(., With_cover=T)),
             cbind("Stat"="Averaged across \n spatial statistics",
                   as_tibble(t(colMeans(d_all[which(d_all$With_cover==F),c(-1,-ncol(d_all))])))%>%
                     add_column(., With_cover=F)),
             d_all)%>%
  
  add_column(., "R2"=0)%>%
  dplyr::relocate(., R2, .before = With_cover)%>%
  filter(., With_cover==with_cov)%>%
  melt(., measure.vars=colnames(.)[5:(ncol(.)-1)])%>%
  Rename_spatial_statistics(.)%>%
  mutate(., variable=recode_factor(variable,"Org_C"="Facilitation",
                                   "Interactions"="Interactions \n with grazing"))%>%
  #arranging order for predictors (xaxis)
  add_column(., Order_stat=rep(rev(1:(length(unique(d_all$Stat))+1)),
                               nrow(.)/(length(unique(d_all$Stat))+1)))%>%
  mutate(Stat = fct_reorder(Stat, Order_stat))%>%
  Rename_spatial_statistics(.)%>%
  
  #arranging order for stats
  add_column(., Order_f=rep(c(sapply(c(4:ncol(mean_imp)),
                                     function(x){
                                       return(which(round(mean_imp[1,x]%>%pull(.),5) ==
                                                      round(sort(as.numeric(mean_imp[1,c(4:ncol(mean_imp))]),
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
                        label=round(sort(as.numeric(
                          as_tibble(t(colMeans(d_all[d_all$With_cover==with_cov,
                                                     c(-1,-ncol(d_all))]))))[-c(1:3)],decreasing = T),2)),
            aes(x=x,y=y,label=paste0(label)),size=3)+
  scale_fill_gradientn(colours = colorRampPalette(c("white","#A129DA",alpha("#07299C",.3)))(100))+
  the_theme+
  theme(axis.text.x = element_text(angle=c(rep(60,(n_var-1)),0),hjust=c(rep(1,(n_var-1)),0)))+
  labs(x="Predictors",y="Spatial statistic",fill="Frequency in best models")

ggsave(paste0("../Figures/Final_figs/SI/Importance_grazing_no_inter.pdf"),p,width = 6,height = 5)






## >> All interactions ----


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
               "core_area_land","Struct1","Struct2")){
  
  title=list_title[ID_title]  
  
  ID_plot=1
  list_plot=list()
  ID_x=1
  
  for (metric in c("Type_veg","Org_C","Clim1","Clim2","Clim3","Clim4","Sand","Woody")){
    
    name_x=list_name_x[ID_x]
    
    d_data_out=  read.table(paste0("../Data/Step1_Understanding_grazing/Keep_data/Data_",stat,"_TRUE.csv"),sep=";")
    if (ncol(d_data_out)==1){
      d_data_out=  read.table(paste0("../Data/Step1_Understanding_grazing/Keep_data/Data_",stat,"_TRUE.csv"),sep=" ")
    }
    model_spa_stat=readRDS(paste0("../Data/Step1_Understanding_grazing/Keep_models/Mod_",stat,"_TRUE.rds"))
    
    
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


  ggsave(paste0("../Figures/Final_figs/SI/All_interactions_",stat,".pdf"),p,width = 4*length(list_plot),height = 4)
  
  ID_title=ID_title+1
  
}










## >> Inference x-y stats ----

keep_sites=read.table("../Data/Inferrence/Keeping_sites.csv",sep=";")
x_y_stat=read.table("../Data/Inferrence/x_y_stat_all_100.csv",sep=";")%>%
  filter(., Site_ID %in% keep_sites$x)

list_plots=list()
name_plot=c("Cover","# neighbors","Clustering","Skewness","Variance","Autocorrelation",
            "SDR","PLR","Exponent p.l.","CV PSD","Frac. max")

for (i in 1:11){
  d_fil=cbind(filter(x_y_stat%>%
                       melt(., id.vars=c("Site_ID","Method", "Type")),variable==colnames(x_y_stat)[i])%>%
                filter(., Type=="Sim")%>%dplyr::rename(., value_sim=value),
              filter(x_y_stat%>%
                       melt(., id.vars=c("Site_ID","Method", "Type")),variable==colnames(x_y_stat)[i])%>%
                filter(., Type=="Obs")%>%dplyr::rename(., value_obs=value)%>%dplyr::select(., value_obs))
  
  list_plots[[i]]=ggplot(d_fil)+
    geom_point(aes(x=value_obs,y=value_sim),color="#96C3DC",alpha=.75)+the_theme+
    labs(x="",y="")+
    geom_abline(slope=1,intercept = 0,color="black")+
    ggtitle(name_plot[i])+
    theme(title = element_text(size=10))
}

p=annotate_figure(ggarrange(plotlist=list_plots,ncol = 4,nrow = 3),
                  left=text_grob("Closest simulations",rot=90,color="black",size=15,face ="bold",vjust=1,family = "NewCenturySchoolbook"),
                  bottom = text_grob("Observed spatial statistic",color="black",size=15,face="bold",vjust=-1,family = "NewCenturySchoolbook"),
                  top = text_grob(paste0("Nkeep = ",unique(x_y_stat$Nkeep)),color="black",size=15,face="bold",vjust=1,family = "NewCenturySchoolbook"))


ggsave(paste0("../Figures/Final_figs/SI/X_Y_inferrence_stats.pdf"),p,width = 8,height = 7)


## >> Inference correlation with variables ----

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

ggplot(d)+
  geom_boxplot(aes(x=as.factor(Grazing),y=mean_q))+
  the_theme





## >> Boxplot grazing and spatial structure ----

pdf("../Figures/Final_figs/SI/Grazing_sp_struct.pdf",width = 5,height = 4)
for (k in c(1:23,70,71)){
  boxplot(d_data[,k]~d_data$Grazing,ylab=colnames(d_data)[k],xlab="Grazing")
}
dev.off()

pdf("../Figures/Final_figs/SI/Grazing_nutrients.pdf",width = 5,height = 4)
for (k in c(45,57:69)){
  boxplot(d_data[,k]~d_data$Grazing,ylab=colnames(d_data)[k],xlab="Grazing")
}
dev.off()
