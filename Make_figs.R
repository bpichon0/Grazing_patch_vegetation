rm(list=ls())
source("./Structure_grazing_function.R")

# -------------------- Main figures -------------------------

## >> Figure 2: Partial residuals grazing intensity spatial structure ----

d_slope=read.table("../Data/Linear_models_factor_cover_control/Slope_partial_residuals_aridity.csv",sep=";")%>%
  dplyr::filter(., With_cover==0,Grazing_intensity=="all",Driver!="Herbivores")%>%
  Filter_relevant_stats(.)%>%  
  Rename_spatial_statistics(.)%>%
  mutate(., ID_grazing=as.character(ID_grazing))%>%
  arrange(., q2)%>%
  add_column(., Order_f=1:nrow(.))%>%
  mutate(.,Stat = fct_reorder(Stat, Order_f))%>%
  add_column(., Signif=.$pval<.05)%>%
  arrange(., Stat,q2)%>%
  add_column(., Order_stat=rep(4:1,nrow(.)/4))%>%
  mutate(.,ID_grazing = fct_reorder(ID_grazing, Order_stat))%>%
  add_column(., Type_stat=sapply(1:nrow(.),function(x){
    if (.$Stat[x] %in% c("Power-law exp. of the PSD","log (largest patch)",
                         "Mean patch size","Number of smallest patches")){
      return("Patch-size")
    }else if (.$Stat[x] =="Bare soil connectivity"){
      return("Hydro.")
    }else if (.$Stat[x] =="Spatial autocorrelation of vege."){
      return("Autocor.")
    }else{
      return("Geom.")
    }
  }))%>%
  mutate(., Stat=recode_factor(Stat,
                               "Mean % of core pixels in patches"="Mean % of core \n pixels in patches",
                               "flow_length"="Bare soil connectivity",
                               "Spatial autocorrelation of vege."="Spatial autocorr. \n of vege.",
                               "Power-law exp. of the PSD"="Power-law exponent \n of the PSD"))


p1=ggplot(d_slope)+
  geom_vline(xintercept = 0,linetype=9)+
  
  geom_linerange(aes(x=q2,y=Stat,xmin=q1_90,xmax=q3_90,color=as.factor(ID_grazing)),
                 lwd=1.5,position=position_jitterdodge(seed=123))+   
  geom_linerange(aes(x=q2,y=Stat,xmin=q1,xmax=q3,color=as.factor(ID_grazing)),
                 lwd=.5,position=position_jitterdodge(seed=123))+  
  geom_pointrange(aes(x=q2,y=Stat,xmin=q1,xmax=q3,fill=as.factor(ID_grazing)),color="transparent",
                  size=.15,position=position_jitterdodge(seed=123),shape=21)+  
  
  the_theme2+
  scale_color_manual(values=c("#FBD2A5","#FF8888","#C17F9D","grey"),
                     labels=c("Low grazing (1)","Medium grazing (2)","High grazing (3)","Aridity"))+
  scale_fill_manual(values=c("grey20","grey20","grey20","grey20"))+
  scale_shape_manual(values=c(21,16))+
  labs(x=substitute(paste(beta," (partial residuals)")),y="",color="")+
  guides(shape="none",fill="none")+
  facet_grid(Type_stat ~ Driver, space = "free_y", 
             scale = "free_y")+
  theme(panel.background = element_rect(fill="white"),
        legend.text = element_text(size=8))+
  guides(color=guide_legend(ncol=4))+
  theme(legend.position="bottom")#c(.2,.65))

d_slope=read.table("../Data/Linear_models_factor_cover_control/Slope_partial_residuals_aridity.csv",sep=";")%>%
  dplyr::filter(., With_cover==0,Grazing_intensity=="all",Stat=="Struct1",Driver!="Herbivores")%>%
  Rename_spatial_statistics(.)%>%
  mutate(., ID_grazing=as.character(ID_grazing))%>%
  arrange(., q2)%>%
  add_column(., Order_f=1:nrow(.))%>%
  mutate(.,Stat = fct_reorder(Stat, Order_f))%>%
  add_column(., Signif=.$pval<.05)%>%
  arrange(., Stat,q2)%>%
  add_column(., Order_stat=rep(4:1,nrow(.)/4))%>%
  mutate(.,ID_grazing = fct_reorder(ID_grazing, Order_stat))%>%
  add_column(., Type_stat=sapply(1:nrow(.),function(x){
    if (.$Stat[x] %in% c("Power-law exp. of the PSD","log (largest patch)",
                         "Mean patch size","Number of smallest patches")){
      return("Patch-size")
    }else if (.$Stat[x] =="Bare soil connectivity"){
      return("Hydro.")
    }else if (.$Stat[x] =="Spatial autocorrelation of vege."){
      return("Autocor.")
    }else{
      return("Geom.")
    }
  }))


p2=ggplot(d_slope%>%mutate(., Type_stat=""))+
  geom_vline(xintercept = 0,linetype=9)+
  
  geom_linerange(aes(x=q2,y=Stat,xmin=q1_90,xmax=q3_90,color=as.factor(ID_grazing)),
                 lwd=1.5,position=position_jitterdodge(seed=123))+   
  geom_linerange(aes(x=q2,y=Stat,xmin=q1,xmax=q3,color=as.factor(ID_grazing)),
                 lwd=.5,position=position_jitterdodge(seed=123))+  
  geom_pointrange(aes(x=q2,y=Stat,xmin=q1,xmax=q3,fill=as.factor(ID_grazing)),color="transparent",
                  size=.15,position=position_jitterdodge(seed=123),shape=21)+  
  
  the_theme2+
  scale_color_manual(values=c("#FBD2A5","#FF8888","#C17F9D","grey"),
                     labels=c("Low grazing (1)","Medium grazing (2)","High grazing (3)","Aridity"))+
  scale_fill_manual(values=c("grey20","grey20","grey20","grey20"))+
  scale_shape_manual(values=c(21,16))+
  labs(x=substitute(paste(beta," (partial residuals)")),y="",color="")+
  guides(shape="none",fill="none")+
  facet_grid(Type_stat ~ ., space = "free_y", 
             scale = "free_y")+
  theme(panel.background = element_rect(fill="white"),
        legend.text = element_text(size=8))+
  guides(color=guide_legend(ncol=4))+
  theme(legend.position="none")#c(.2,.65))

p_tot=ggarrange(
  p1+theme(axis.text.y = element_text(size=10),legend.text = element_text(size=11))
    ,ggarrange(ggplot()+theme_void(),
                             p2+theme(axis.text.y = element_text(size=10))
                  ,ggplot()+theme_void(),ncol=3,widths = c(.7,1.2,.5)),nrow = 2,labels = letters[1:2],heights = c(1,.25))

ggsave("../Figures/Grazing_intensity_partial_residuals.pdf",p_tot,width = 8,height = 6)


## >> Figure 3: Illustrating the effects with simulations ----

d_slope=read.table("../Data/Simulations/Model_coefficients.csv",sep=";")%>%
  Rename_spatial_statistics(.)%>%
  add_column(., Type_stat=sapply(1:nrow(.),function(x){
    if (.$Stat[x] %in% c("Power-law exp. of the PSD","log (largest patch)",
                         "Mean patch size","Number of smallest patches")){
      return("Patch-size")
    }else if (.$Stat[x] =="Bare soil connectivity"){
      return("Hydro.")
    }else if (.$Stat[x] =="Spatial autocorrelation of vege."){
      return("Autocor.")
    }else{
      return("Geom.")
    }
  }))

d_slope$shape_="Coherent with obs."
d_slope$shape_[which(d_slope$Stat=="log (largest patch)" | (d_slope$Cover=="High" &  d_slope$Stat=="Mean patch size"))]="Different from obs."

p1=ggplot(d_slope)+
  geom_vline(xintercept = 0,linetype=9)+
  geom_point(aes(x=q2,y=Stat,fill=Cover,color=Cover,shape=shape_),size=5)+  
  the_theme2+
  labs(x=substitute(paste(beta," increasing grazing over aridity")),y="",color="")+
  scale_fill_manual(values=rev(c("lightgreen","forestgreen","#025a2a")),
                    labels=rev(c("Low cover","Medium over","High cover")))+
  scale_color_manual(values=rev(c("lightgreen","forestgreen","#025a2a")),
                     labels=rev(c("Low cover","Medium over","High cover")))+
  scale_shape_manual(values=rev(c(23,21)))+
  guides()+
  facet_grid(Type_stat ~ ., space = "free_y", 
             scale = "free_y")+
  theme(panel.background = element_rect(fill="white"),legend.box="vertical",
        legend.text = element_text(size=14),strip.text.y = element_text(size=16),
        axis.text.y = element_text(size=16),axis.title.x = element_text(size=16))+
  labs(fill="",color="",shape="")+
  geom_text(data=tibble(x=-.95,y=1.1,Type_stat="Patch-size",label="*"),aes(x=x,y=y,label=label),size=10)



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

list_landscape=readRDS("../Data/Simulations/Minimal_examples_stats_landscapes.rds")

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
                         theme(panel.border = element_rect(colour = "#025a2a",fill=NA,size=2)),
                       p12+
                         theme(panel.border = element_rect(colour = "#025a2a",fill=NA,size=2)),
                       p11+
                         theme(panel.border = element_rect(colour = "#025a2a",fill=NA,size=2)),
                       p23+
                         theme(panel.border = element_rect(colour = "forestgreen",fill=NA,size=2)),
                       p22+
                         theme(panel.border = element_rect(colour = "forestgreen",fill=NA,size=2)),
                       p21+
                         theme(panel.border = element_rect(colour = "forestgreen",fill=NA,size=2)),
                       p33+
                         theme(panel.border = element_rect(colour = "lightgreen",fill=NA,size=2)),
                       p32+
                         theme(panel.border = element_rect(colour = "lightgreen",fill=NA,size=2)),
                       p31+
                         theme(panel.border = element_rect(colour = "lightgreen",fill=NA,size=2)),
                       nrow=3,ncol=3)
p_arrow_top=ggplot(NULL) + 
  geom_line(data=tibble(x=c(0,1),y=c(1,1)),
            aes(x=x,y=y),
            arrow = arrow(length=unit(0.30,"cm"), type = "closed"))+
  theme_transparent()+
  ylim(.96,1.2)+
  xlim(0,1)+
  annotate("text",x=.5,y=1.1,label = "Increasing spatially heterogeneous mortality",size=7,
           family = "NewCenturySchoolbook")+
  theme(plot.margin = unit(rep(.5,4), 'lines'))

p_arrow_left=ggplot(NULL) + 
  geom_line(data=tibble(x=c(1,1),y=c(0,1)),
            aes(x=x,y=y),
            arrow = arrow(length=unit(0.30,"cm"), type = "closed"))+
  annotate("text",x=.9,y=.5,label = "Increasing vegetation cover",angle = 90,size=7,
           family = "NewCenturySchoolbook")+
  theme_transparent()+
  ylim(0,1)+
  xlim(.8,1.1)+
  theme(plot.margin = unit(rep(.5,4), 'lines'))

p_landscapes=ggarrange(p_arrow_top,p_landscapes,nrow=2,heights = c(.1,1))
p_landscapes=ggarrange(ggarrange(ggplot()+theme_void(),p_arrow_left,nrow=2,heights = c(.1,1)),
                       p_landscapes,ncol=2,widths = c(.1,1))

p_tot=ggarrange(p_landscapes,p1,ncol=2,widths = c(1.4,1),labels = c("a","b"),font.label = list(size=24))
ggsave("../Figures/Effect_size_simulations.pdf",p_tot,width = 18,height = 8)


## >> Figure 4: Traits ----
#See outputs of the SEMs in Step 4 of "Structure_grazing_main.R
d_all=read.table("../Data/Dsep/CI_links_dsep.csv",sep=";")%>%
  dplyr::filter(., Variable %!in% c("(Intercept)","Elevation","Long_cos","Long_sin","Lattitude"))

#Computing indirect effects

AI_woody=d_all$q2[which(d_all$Response=="Woody" & d_all$Variable=="Aridity")]

Grazing_woody=d_all$q2[which(d_all$Response=="Woody" & d_all$Variable=="Grazing")]

AI_size_D=d_all$q2[which(d_all$Response=="CWM_MaxLS" & d_all$Variable=="Aridity")]

Grazing_size_D=d_all$q2[which(d_all$Response=="CWM_MaxLS" & d_all$Variable=="Grazing")]

AI_size_I=d_all$q2[which(d_all$Response=="Woody" & d_all$Variable=="Aridity")] * 
  d_all$q2[which(d_all$Response=="CWM_MaxLS" & d_all$Variable=="Woody")]

Grazing_size_I=d_all$q2[which(d_all$Response=="Woody" & d_all$Variable=="Grazing")] * 
  d_all$q2[which(d_all$Response=="CWM_MaxLS" & d_all$Variable=="Woody")]

AI_size_T=sum(c(d_all$q2[which(d_all$Response=="Woody" & d_all$Variable=="Aridity")],
                d_all$q2[which(d_all$Response=="CWM_MaxLS" & d_all$Variable=="Aridity")])*
                c(d_all$q2[which(d_all$Response=="CWM_MaxLS" & d_all$Variable=="Woody")],1))

Grazing_size_T=sum(c(d_all$q2[which(d_all$Response=="Woody" & d_all$Variable=="Grazing")],
                     d_all$q2[which(d_all$Response=="CWM_MaxLS" & d_all$Variable=="Grazing")])*
                     c(d_all$q2[which(d_all$Response=="CWM_MaxLS" & d_all$Variable=="Woody")],1))

AI_DevLS_D=d_all$q2[which(d_all$Response=="Dev_MaxLS" & d_all$Variable=="Aridity")]

Grazing_DevLS_D=d_all$q2[which(d_all$Response=="Dev_MaxLS" & d_all$Variable=="Grazing")]

AI_DevLS_I=AI_size_T*
  d_all$q2[which(d_all$Response=="Dev_MaxLS" & d_all$Variable=="CWM_MaxLS")]+
  d_all$q2[which(d_all$Response=="Woody" & d_all$Variable=="Aridity")]*
  d_all$q2[which(d_all$Response=="Dev_MaxLS" & d_all$Variable=="Woody")]

Grazing_DevLS_I=Grazing_size_T*
  d_all$q2[which(d_all$Response=="Dev_MaxLS" & d_all$Variable=="CWM_MaxLS")]+
  d_all$q2[which(d_all$Response=="Woody" & d_all$Variable=="Grazing")]*
  d_all$q2[which(d_all$Response=="Dev_MaxLS" & d_all$Variable=="Woody")]

AI_DevLS_T=AI_size_T*
  d_all$q2[which(d_all$Response=="Dev_MaxLS" & d_all$Variable=="CWM_MaxLS")]+
  d_all$q2[which(d_all$Response=="Woody" & d_all$Variable=="Aridity")]*
  d_all$q2[which(d_all$Response=="Dev_MaxLS" & d_all$Variable=="Woody")]+
  d_all$q2[which(d_all$Response=="Dev_MaxLS" & d_all$Variable=="Aridity")]

Grazing_DevLS_T=Grazing_size_T*
  d_all$q2[which(d_all$Response=="Dev_MaxLS" & d_all$Variable=="CWM_MaxLS")]+
  d_all$q2[which(d_all$Response=="Woody" & d_all$Variable=="Grazing")]*
  d_all$q2[which(d_all$Response=="Dev_MaxLS" & d_all$Variable=="Woody")]+
  d_all$q2[which(d_all$Response=="Dev_MaxLS" & d_all$Variable=="Grazing")]

AI_PC1=AI_size_T*d_all$q2[which(d_all$Response=="Struct_PC1" & d_all$Variable=="CWM_MaxLS")]+ #Size
  
  AI_DevLS_T*d_all$q2[which(d_all$Response=="Struct_PC1" & d_all$Variable=="Dev_MaxLS")]+    #dev maxLS
  
  AI_woody*d_all$q2[which(d_all$Response=="Struct_PC1" & d_all$Variable=="Woody")]              #Woody

Grazing_PC1=Grazing_size_T*d_all$q2[which(d_all$Response=="Struct_PC1" & d_all$Variable=="CWM_MaxLS")]+ #Size
  
  Grazing_DevLS_T*d_all$q2[which(d_all$Response=="Struct_PC1" & d_all$Variable=="Dev_MaxLS")]+    #dev maxLS
  
  Grazing_woody*d_all$q2[which(d_all$Response=="Struct_PC1" & d_all$Variable=="Woody")]              #Woody



d_indirect=tibble(
  Indirect=c(0,0,AI_size_I,Grazing_size_I,AI_DevLS_I,Grazing_DevLS_I,AI_PC1,Grazing_PC1),
  Direct=c(AI_woody,Grazing_woody,AI_size_D,Grazing_size_D,AI_DevLS_D,Grazing_DevLS_D,0,0),
  Total=c(AI_woody,Grazing_woody,AI_size_T,Grazing_size_T,AI_DevLS_T,Grazing_DevLS_T,AI_PC1,Grazing_PC1),
  Effect=rep(c("Aridity","Grazing"),4),
  Response=c(rep("% of Woody",2),rep("Average plant size (cLS)",2),rep("Trait spatial aggregation (DevLS)",2),rep("PC 1 sp. struct.",2))
)%>%melt(., id.vars=c("Effect","Response"))


for (k in 1:4){
  assign(paste0("p2",k),
         ggplot(d_indirect[-which(d_indirect$Response=="PC 1 sp. struct." & d_indirect$variable %in% c("Direct","Indirect")|
                                    d_indirect$Response=="% of Woody" & d_indirect$variable %in% c("Total","Indirect")),]%>%
                  dplyr::filter(., Response==unique(d_indirect$Response)[k]))+
                             geom_bar(aes(x=value,y=(variable),
                                          group=interaction(variable,Response,Effect),fill=Effect),
                                      stat="identity", 
                                      width =.5, position=position_dodge())+
                             the_theme2+
                             geom_hline(yintercept = 0)+
                             facet_wrap(.~Response,scales = "free")+
                             coord_flip()+
                             labs(x="Effect size",y="",color="")+
                             guides(fill="none")+
                             scale_fill_manual(values=c("#FFE699","#F3A875"))
         
         )
}


p_tot=ggarrange(
  ggarrange(p22+theme(panel.border = element_rect(color="#A8D18F",size=1),
                      axis.line = element_line(color="#A8D18F",size=1)),
            p21+theme(panel.border = element_rect(color="#A8D18F",size=1),
                      axis.line = element_line(color="#A8D18F",size=1)),ncol=2,widths = c(1,.5)),
  ggarrange(p23+theme(panel.border = element_rect(color="#A8D18F",size=1),
                      axis.line = element_line(color="#A8D18F",size=1)),
            p24+theme(panel.border = element_rect(color="#C376EF",size=1),
                      axis.line = element_line(color="#C376EF",size=1)),ncol=2,widths = c(1,.5)),
  nrow=2
)
p_tot=ggarrange(ggplot()+theme_void(),p_tot,ncol=2,widths =  c(1,.8),labels = letters[1:2])
ggsave("../Figures/SEM_traits_spatial_structure.pdf",p_tot,width = 12,height = 6)

 

## >> Figure 5: Type herbivores ----

d_slope=read.table("../Data/Linear_models_factor_cover_control/Slope_partial_residuals_aridity.csv",sep=";")%>%
  dplyr::filter(., With_cover==0,Grazing_intensity=="all",Driver=="Herbivores")%>%
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
      return("Autocor.")
    }else{
      return("Geom.")
    }
  }))%>%
  mutate(., Stat=recode_factor(Stat,
                               "Mean % of core pixels in patches"="Mean % of core \n pixels in patches",
                               "flow_length"="Bare soil connectivity",
                               "Spatial autocorrelation of vege."="Spatial autocorr. \n of vege.",
                               "Power-law exp. of the PSD"="Power-law exponent \n of the PSD"))

p=ggplot(d_slope)+
  geom_vline(xintercept = 0,linetype=9)+
  
  geom_linerange(aes(x=q2,y=Stat,xmin=q1_90,xmax=q3_90,color=as.factor(ID_grazing)),
                 lwd=1.5,position=position_jitterdodge(seed=123))+   
  geom_linerange(aes(x=q2,y=Stat,xmin=q1,xmax=q3,color=as.factor(ID_grazing)),
                 lwd=.5,position=position_jitterdodge(seed=123))+  
  geom_pointrange(aes(x=q2,y=Stat,xmin=q1,xmax=q3,fill=as.factor(ID_grazing)),color="transparent",
                  size=.15,position=position_jitterdodge(seed=123),shape=21)+  
  
  the_theme2+
  scale_fill_manual(values=c("grey20","grey20","grey20","grey20"))+
  scale_shape_manual(values=c(21,16))+
  labs(x=substitute(paste(beta," (partial residuals, relative to cattle)")),y="",color="")+
  guides(shape="none",fill="none")+
  theme(panel.background = element_rect(fill="white"),
        legend.text = element_text(size=8))+
  guides(color=guide_legend(ncol=4))+
  scale_color_manual(values=c("orange","violet","lightblue"))+
  facet_grid(Type_stat ~ ., space = "free_y", 
             scale = "free_y")
  


# 
# d_data2=Perform_PCA_spatial_struc(
#   Closer_to_normality_CWM(
#     read.table("../Data/Spatial_structure_control_cover_with_traits.csv",sep=";")))%>%
#   mutate(., Grazing=as.factor(Grazing),ID=as.character(ID))%>%
#   mutate(., Dev_MaxLS=-Dev_MaxLS,#positive values = facilitation
#          Grazing=scale(as.numeric(Grazing))[,1])
# 
# 
# p2=ggplot(d_data2%>%dplyr::filter(., !is.na(Herbivores)))+
#   geom_boxplot(aes(x=Herbivores,y=CWM_MaxLS,fill=Herbivores),width=.3,outlier.shape = NA)+
#   scale_fill_manual(values=c("orange","violet","lightblue","forestgreen"))+
#   the_theme2+
#   labs(x="",y="Average plant size \n (cLS)")
# 
# 
# p3=ggplot(d_data2%>%dplyr::filter(., !is.na(Herbivores)))+
#   geom_boxplot(aes(x=Herbivores,y=Dev_MaxLS,fill=Herbivores),width=.3,outlier.shape = NA)+
#   scale_fill_manual(values=c("orange","violet","lightblue","forestgreen"))+
#   the_theme2+
#   labs(x="",y="Trait aggregation \n (DevLS)")
# 
# 
# p4=ggplot(d_data2%>%dplyr::filter(., !is.na(Herbivores)))+
#   geom_boxplot(aes(x=Herbivores,y=Woody,fill=Herbivores),width=.3,outlier.shape = NA)+
#   scale_fill_manual(values=c("orange","violet","lightblue","forestgreen"))+
#   the_theme2+
#   labs(x="",y="% of woody \n species")
# 
# p_tot=ggarrange(p,ggarrange(p2,p3,p4,nrow=3,labels = letters[2:4],common.legend = T,legend = "none"),
#           labels = c("a",""),ncol=2,widths = c(1.2,1))

p_tot=p
ggsave("../Figures/Type_herbivores.pdf",p_tot,width = 5,height = 5)

## >> Figure 6: Consequence resilience ----

# 
# d=read.table("./All_d.csv",sep=";")%>%
#   add_column(., Cover_binary=.$Cover>0.03)%>%
#   mutate(., 
#          g0=rep(seq(0,.5,length.out=150),each=300),
#          b=rep(rep(rev(seq(0,1,length.out=150)),each=2),150))
# d$Cover[is.na(d$Cover)]="Desert"
# d$Type_cover="Normal"
# d$Type_cover[which(d$Cover=="Desert")]="Desert"
# d$Cover[d$Cover %in% c("Vegetation","Desert")]=NA
# d$Cover=as.numeric(d$Cover)
# 
# d_restor=d%>%dplyr::filter(.,Branch=="Restoration")
# d_bistab=d%>%dplyr::filter(.,Branch=="Degradation")%>%
#   add_column(., Bistab=sapply(1:nrow(.),function(x){
#     if (.$Cover_binary[x]==d_restor$Cover_binary[x]){
#       return("No")
#     }else{
#       return("Yes")
#     }
#   }))
# 
# 
# d_bistab$Cover[d_bistab$Cover==0]=NA
# p=ggplot(NULL)+
#   geom_point(data=d_bistab,
#              aes(x=1-b,g0,color=Cover))+
#   scale_color_viridis_c(option = "G",na.value = "#D8C0A5")+
#   the_theme2+
#   labs(x="Decreasing recruitment probability (i.e., increasing stress)",y="Increasing spatially heterogeneous mortality",color="Vegetation \n   cover")+
#   theme(axis.ticks = element_blank(),axis.text = element_blank(),axis.title.x = element_text(size=16),
#         axis.title.y = element_text(size=16),legend.title = element_text(size = 16))+
#   annotate("text",x=.75,y=.4,label="No vegetation",size=6)
# 
# size_bistab=sapply(unique(d$g0),function(x){
#   d_fil=dplyr::filter(d_bistab,g0==x)
#   return(ifelse(any(d_fil$Bistab=="Yes"),table(d_fil$Bistab)[2]/150,0))
# })
# b_crit=sapply(unique(d$g0),function(x){
#   d_fil=dplyr::filter(d_bistab,g0==x)
#   return(max(d_fil$b[which(is.na(d_fil$Cover))]))
# })
# 
# p2=ggplot(tibble(b_crit=1-b_crit,g0=unique(d$g0)))+
#   geom_smooth(aes(x=g0,y=b_crit),shape=21,size=2,fill="lightblue",color="lightblue")+
#   labs(x="Increasing spatially heterogeneous mortality",y="Recruitment prob. at which \n the ecosystemt desertify")+
#   the_theme2+
#   theme(axis.ticks = element_blank(),axis.text = element_blank(),axis.title.x = element_text(size=16),
#         axis.title.y = element_text(size=16))
# 
# p3=ggplot(tibble(size_bistab=size_bistab,g0=unique(d$g0)))+
#   geom_smooth(aes(x=g0,y=size_bistab),shape=21,size=2,fill="lightblue",color="lightblue")+
#   labs(x="Increasing spatially heterogeneous mortality",y="Width of the bistability region")+
#   the_theme2+
#   theme(axis.ticks = element_blank(),axis.text = element_blank(),axis.title.x = element_text(size=16),
#         axis.title.y = element_text(size=16))

d=read.table("../Data/Simulations/All_simulations_PA.csv",sep=";")%>%
  add_column(., Cover_binary=.$Cover>0.03)%>%
  mutate(., 
         g0=rep(seq(0,.3,length.out=50),each=100),
         b=rep(rep(rev(seq(0,1,length.out=50)),each=2),50))
d$Cover[is.na(d$Cover)]="Desert"
d$Type_cover="Normal"
d$Type_cover[which(d$Cover=="Desert")]="Desert"
d$Cover[d$Cover %in% c("Vegetation","Desert")]=NA
d$Cover=as.numeric(d$Cover)

d_restor=d%>%dplyr::filter(.,Branch=="Restoration")
d_bistab=d%>%dplyr::filter(.,Branch=="Degradation")%>%
  add_column(., Bistab=sapply(1:nrow(.),function(x){
    if (.$Cover_binary[x]==d_restor$Cover_binary[x]){
      return("No")
    }else{
      return("Yes")
    }
  }))


d_bistab$Cover[d_bistab$Cover==0]=NA

size_bistab=sapply(unique(d$g0),function(x){
  d_fil=dplyr::filter(d_bistab,g0==x)
  return(ifelse(any(d_fil$Bistab=="Yes"),table(d_fil$Bistab)[2]/50,0))
})

cover_desert=sapply(unique(d$g0),function(x){
 d_fil=dplyr::filter(d_bistab,g0==x)
 return(d_fil$Cover[min(which(is.na(d_fil$Cover)))-1])
})


p1=ggplot(tibble(size_bistab=size_bistab,g0=unique(d$g0)))+
  geom_smooth(aes(x=g0,y=size_bistab),size=3,fill="lightblue",color="lightblue")+
  labs(x="Increasing spatially heterogeneous mortality",y="Width of the bistability region")+
  the_theme2+
  theme(axis.ticks = element_blank(),axis.text = element_blank(),axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14))

p2=ggplot(tibble(cover_desert=cover_desert,g0=unique(d$g0)))+
  geom_smooth(aes(x=g0,y=cover_desert),size=3,fill="lightblue",color="lightblue")+
  labs(x="Increasing spatially heterogeneous mortality",y="Cover at the desertification point")+
  the_theme2+
  theme(axis.ticks = element_blank(),axis.text = element_blank(),axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14))

p3=ggplot(NULL)+
  geom_point(data=d%>%dplyr::filter(g0==.0),
             aes(x=1-b,Cover),shape=21,size=2,fill="lightblue",color="black")+
  the_theme2+
  labs(x="Decreasing recruitment prob. (i.e., higher stress)",y="Vegetation cover")+
  theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(),axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14),legend.title = element_text(size = 14))+
  annotate("text",x=.75,y=.9,label="Spatially homogeneous \n stress",size=5)+
  ylim(0,1)
  
p4=ggplot(NULL)+
  geom_point(data=d%>%dplyr::filter(g0==.3),
             aes(x=1-b,Cover),shape=21,size=2,fill="lightblue",color="black")+
  the_theme2+
  labs(x="Decreasing recruitment prob. (i.e., higher stress)",y="Vegetation cover")+
  theme(axis.ticks.X = element_blank(),axis.text.X = element_blank(),axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14),legend.title = element_text(size = 14))+
  annotate("text",x=.75,y=.9,label="Spatially heterogeneous \n stress",size=5)+
  ylim(0,1)


p0=ggplot(NULL)+
  geom_tile(data=d_bistab,
             aes(x=1-b,g0,fill=Cover))+
  scale_fill_viridis_c(option = "G",na.value = "#D8C0A5")+
  the_theme2+
  labs(x="Decreasing recruitment probability (i.e., increasing stress)",y="Increasing spatially \n heterogeneous mortality",
       fill="Vegetation \n    cover")+
  theme(axis.ticks = element_blank(),axis.text = element_blank(),axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14),legend.title = element_text(size = 14))+
  annotate("text",x=.75,y=.15,label="No vegetation",size=5)+
  geom_hline(yintercept = c(0,.3))


p_tot=ggarrange(ggarrange(p0,
                          ggarrange(p3+theme(axis.title.x = element_blank()),p4,nrow=2,labels = c("b","c"),font.label = list(size=22)),
                          ncol=2,labels = c("a",""),font.label = list(size=22),widths = c(1,.7)),
                ggarrange(p2,p1,ncol=2,labels = c("d","e"),font.label = list(size=22)),
                nrow = 2,heights = c(1,.8))
ggsave("../Figures/Consequence_resilience.pdf",p_tot,width = 12,height = 12)

# -------------------- SI figures ---------------------------
## >> Selection criteria for images ----
info_kmean=read.table("../Data/Landscapes/Kmean_clust_info.csv",sep=",",header = T)%>%
  dplyr::filter(., Size!=200,Dataset=="biodesert")%>%
  add_column(., Quadratic_error=(.$field_cover-.$img_cover)**2)

trade_off=tibble(threshold=seq(0,.05,length.out=100))%>%
  add_column(., N_site_kept=sapply(1:nrow(.),function(x){
    return(length(which(info_kmean$Quadratic_error>.$threshold[x])))
  }))

p1=ggplot(trade_off)+
  geom_point(aes(threshold,N_site_kept))+
  the_theme+
  geom_vline(xintercept = 0.01,color="red",lwd=1)+
  labs(y="# of images kept",x="(Field cover - landscape cover)^2")


p2=ggplot(info_kmean%>%
            mutate(., status=recode_factor(status,"kept"="Kept",
                                           "removed"="Removed (error > 0.01)",
                                           "removed_visual"="Removed after visual inspection")))+
  geom_point(aes(x=field_cover,y=img_cover,color=status))+
  geom_smooth(data=info_kmean%>%dplyr::filter(., status=="kept"),
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
                                    "moran_I","core_area","Small_patches",
                                    "mean_psd"))%>%
           mutate(., variable=recode_factor(variable,
                                            "core_area"="Mean % of core pixels in patches",
                                            "flow_length"="Bare soil connectivity",
                                            "Small_patches"="Number of smallest patches",
                                            "mean_psd"="Mean patch size",
                                            "moran_I"="Spatial autocorrelation of vege.",
                                            "PL_expo"="Power-law exp. of the PSD",
                                            "fmax_psd"="log (largest patch)")))+
  geom_boxplot(aes(x=Grazing,y=value,group=as.factor(Grazing)),outlier.shape = NA)+
  geom_jitter(aes(x=Grazing,y=value,color=as.factor(Grazing)),size=.5,alpha=.5)+
  scale_color_manual(values=c("grey","#FBD2A5","#FF8888","#C17F9D"),
                     label=c("Ungrazed","Low","Medium","High"))+
  facet_wrap(.~variable,scales = "free",nrow = 3)+
  scale_x_continuous(breaks = c(0,1,2,3),labels=c("Ungrazed","Low","Medium","High"))+
  the_theme2+
  labs(x="Grazing intensity",y="")+
  theme(legend.position = "none")

ggsave("../Figures/SI/Distribution_spatial_stat_grazing.pdf",p,width = 11,height = 7)


d_data=read.table("../Data/Spatial_structure_grazing.csv",sep=";")%>%
  Closer_to_normality(.)

p=ggplot(d_data)+
  geom_boxplot(aes(x=Grazing,y=rho_p,group=as.factor(Grazing)),outlier.shape = NA)+
  geom_jitter(aes(x=Grazing,y=rho_p,color=as.factor(Grazing)),size=1,alpha=.5)+
  scale_color_manual(values=c("grey","#FBD2A5","#FF8888","#C17F9D"),
                     label=c("Ungrazed","Low","Medium","High"))+
  the_theme2+
  labs(x="Grazing intensity",y="Vegetation cover")+
  scale_x_continuous(breaks = c(0,1,2,3),labels=c("Ungrazed","Low","Medium","High"))+
  theme(legend.position = "none")

ggsave("../Figures/SI/Grazing_on_cover.pdf",p,width = 5,height = 4)

d_data=read.table("../Data/Spatial_structure_grazing_control_cover.csv",sep=";")

p=ggplot(d_data%>%
           melt(., measure.vars = c("PL_expo","fmax_psd","flow_length",
                                    "moran_I","core_area","Small_patches",
                                    "mean_psd"))%>%
           mutate(., variable=recode_factor(variable,
                                            "core_area"="Mean % of core pixels in patches",
                                            "flow_length"="Bare soil connectivity",
                                            "Small_patches"="Number of smallest patches",
                                            "mean_psd"="Mean patch size",
                                            "moran_I"="Spatial autocorrelation of vege.",
                                            "PL_expo"="Power-law exp. of the PSD",
                                            "fmax_psd"="log (largest patch)")))+
  geom_boxplot(aes(x=Grazing,y=value,group=as.factor(Grazing)),outlier.shape = NA)+
  geom_jitter(aes(x=Grazing,y=value,color=as.factor(Grazing)),size=.5,alpha=.5)+
  scale_color_manual(values=c("grey","#FBD2A5","#FF8888","#C17F9D"))+
  facet_wrap(.~variable,scales = "free",nrow = 3)+
  scale_x_continuous(breaks = c(0,1,2,3),labels=c("Ungrazed","Low","Medium","High"))+
  the_theme2+
  labs(x="Grazing intensity",y="")+
  theme(legend.position = "none")

ggsave("../Figures/SI/Distribution_spatial_stat_grazing_control_cover.pdf",p,width = 11,height = 7)




d_data=read.table("../Data/Spatial_structure_grazing.csv",sep=";")%>%
  Closer_to_normality(.)

p=ggplot(d_data%>%
           melt(., measure.vars = c("PL_expo","fmax_psd","flow_length",
                                    "moran_I","core_area","Small_patches",
                                    "mean_psd"))%>%
           mutate(., variable=recode_factor(variable,
                                            "core_area"="Mean % of core pixels in patches",
                                            "flow_length"="Bare soil connectivity",
                                            "Small_patches"="Number of smallest patches",
                                            "mean_psd"="Mean patch size",
                                            "moran_I"="Spatial autocorrelation of vege.",
                                            "PL_expo"="Power-law exp. of the PSD",
                                            "fmax_psd"="log (largest patch)")))+
  geom_boxplot(aes(x=Grazing,y=value,group=as.factor(Grazing)),outlier.shape = NA)+
  geom_jitter(aes(x=Grazing,y=value,color=as.factor(Grazing)),size=.5,alpha=.5)+
  scale_color_manual(values=c("orange","violet","lightblue","forestgreen"))+
  facet_wrap(.~variable,scales = "free",nrow = 3)+
  scale_x_continuous(breaks = c(0,1,2,3),labels=c("Ungrazed","Low","Medium","High"))+
  the_theme2+
  labs(x="Grazing intensity",y="")+
  theme(legend.position = "none")

ggsave("../Figures/SI/Distribution_spatial_stat_herbivores.pdf",p,width = 11,height = 7)

d_data=read.table("../Data/Spatial_structure_grazing_control_cover.csv",sep=";")

p=ggplot(d_data%>%
           melt(., measure.vars = c("PL_expo","fmax_psd","flow_length",
                                    "moran_I","core_area","Small_patches",
                                    "mean_psd"))%>%
           mutate(., variable=recode_factor(variable,
                                            "core_area"="Mean % of core pixels in patches",
                                            "flow_length"="Bare soil connectivity",
                                            "Small_patches"="Number of smallest patches",
                                            "mean_psd"="Mean patch size",
                                            "moran_I"="Spatial autocorrelation of vege.",
                                            "PL_expo"="Power-law exp. of the PSD",
                                            "fmax_psd"="log (largest patch)")))+
  geom_boxplot(aes(x=Grazing,y=value,group=as.factor(Grazing)),outlier.shape = NA)+
  geom_jitter(aes(x=Grazing,y=value,color=as.factor(Grazing)),size=.5,alpha=.5)+
  scale_color_manual(values=c("orange","violet","lightblue","forestgreen"))+
  facet_wrap(.~variable,scales = "free",nrow = 3)+
  scale_x_continuous(breaks = c(0,1,2,3),labels=c("Ungrazed","Low","Medium","High"))+
  the_theme2+
  labs(x="Grazing intensity",y="")+
  theme(legend.position = "none")

ggsave("../Figures/SI/Distribution_spatial_stat_herbivores_control_cover.pdf",p,width = 11,height = 7)



## >> Spatial metrics and aridity ----

d_data=read.table("../Data/Spatial_structure_grazing.csv",sep=";")%>%
  Closer_to_normality(.)

p=ggplot(d_data%>%
           melt(., measure.vars = c("PL_expo","fmax_psd","flow_length",
                                    "core_area","mean_psd","Small_patches",
                                    "moran_I"))%>%
           mutate(., variable=recode_factor(variable,
                                            "mean_psd"="Mean patch size",
                                            "core_area"="Mean % of core pixels in patches",
                                            "Small_patches"="Number of smallest patches",
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
  labs(x="Spatial resolution (pixel size in m²)",y="Spatial statistic")


ggsave("../Figures/SI/Correlation_metrics_Resolution.pdf",p,width = 9,height = 5)


d_slope=read.table("../Data/Linear_models_resolution/Slope_partial_residuals_aridity.csv",sep=";")%>%
  dplyr::filter(., With_cover==0,Driver=="Grazing",Grazing_intensity=="all")%>%
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
       ggarrange(p1+ggtitle("Not controlled for vegetation cover")+theme(title = element_text(size=14)),
                 p2+ggtitle("Controlled for vegetation cover")+theme(title = element_text(size=14)),
                 nrow = 2,
                 labels = letters[1:2]),width = 10,height = 12)


## >> Partial residuals without accounting for vegetation cover ----

d_slope=read.table("../Data/Linear_models_factor/Slope_partial_residuals_aridity.csv",sep=";")%>%
  dplyr::filter(., With_cover==0,Driver=="Grazing",Grazing_intensity=="all")%>%
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

ggsave("../Figures/SI/Grazing_intensity_partial_residuals_without_cover.pdf",p,width = 5,height = 5.5)



## >> Toy examples Schneider model PSD and dynamics with mean psd ----

list_landscape=readRDS("../Data/Simulations/Minimal_examples_stats_landscapes.rds")
d=read.table("../Data/Simulations/Minimal_examples_stats.csv",sep=";")


id=1
for (cover_id in 1:3){
  for (aggregation_id in 1:3){
    assign(paste0("p_psd",cover_id,aggregation_id),
           Plot_psd_raw(list_landscape[[id]])+
             scale_y_log10(limits=c(1,1000)))
    id=id+1
  }
}


p_psd=ggarrange(
  p_psd11+ggtitle("Cover = 0.65, \nHigh asso. protection")+
    theme(panel.border = element_rect(colour = "#D27C3F",fill=NA,size=2)),
  p_psd12+ggtitle("Cover = 0.65, \nMedium asso. protection")+
    theme(panel.border = element_rect(colour = "#D27C3F",fill=NA,size=2)),
  p_psd13+ggtitle("Cover = 0.65, \nNo asso. protection")+
    theme(panel.border = element_rect(colour = "#D27C3F",fill=NA,size=2)),
  # p_psd41+ggtitle("Cover = 0.65, n random"),
  p_psd21+ggtitle("Cover = 0.35, \nHigh asso. protection")+
    theme(panel.border = element_rect(colour = "#8BD882",fill=NA,size=2)),
  p_psd22+ggtitle("Cover = 0.35, \nMedium asso. protection")+
    theme(panel.border = element_rect(colour = "#8BD882",fill=NA,size=2)),
  p_psd23+ggtitle("Cover = 0.35, \nNo asso. protection")+
    theme(panel.border = element_rect(colour = "#8BD882",fill=NA,size=2)),
  # p_psd42+ggtitle("Cover = 0.35, n random"),
  p_psd31+ggtitle("Cover = 0.17, \nHigh asso. protection")+
    theme(panel.border = element_rect(colour = "#BE91D0",fill=NA,size=2)),
  p_psd32+ggtitle("Cover = 0.17, \nMedium asso. protection")+
    theme(panel.border = element_rect(colour = "#BE91D0",fill=NA,size=2)),
  p_psd33+ggtitle("Cover = 0.17, \nNo asso. protection")+
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
    summary_lm=summary(lm(data = dplyr::filter(partial_res$res,Grazing==grazing_id),visregRes~Aridity))
    
    d=rbind(d,tibble(Stat=k,
                     pval=summary_lm$coefficients[2,4],
                     Grazing_id=grazing_id))
    
    if (summary_lm$coefficients[2,4]<.05){
      
      assign(paste0("p_",grazing_id),
             ggplot(dplyr::filter(partial_res$res,Grazing==grazing_id))+
               geom_point(aes(x=Aridity,y=visregRes),fill=color_grazing[as.numeric(grazing_id)+1],color="black",shape=21)+
               geom_smooth(aes(x=Aridity,y=visregRes),method = "lm",color="black",fill=color_grazing[as.numeric(grazing_id)+1])+
               the_theme2+
               labs(x="Aridity",y=name_stats[id]))
    } else{
      assign(paste0("p_",grazing_id),
             ggplot(dplyr::filter(partial_res$res,Grazing==grazing_id))+
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
                p_tot_5,
                p_tot_2,
                p_tot_3,
                p_tot_4,
                p_tot_6,
                p_tot_7,
                nrow=7,align = "hv",
                labels=letters[1:7])

ggsave("../Figures/SI/Interactive_effects.pdf",p_fig,width = 9,height = 14)


# 
# 
# d=tibble()
# 
# name_stats=c("Power-law exp. \n of the PSD","log (largest patch)",
#              "Bare soil connectivity","Spatial autocorrelation \n of vege.",
#              "Number of \n smallest patches",
#              "Mean % of core \n pixels in patches",
#              "Mean patch size")
# 
# color_grazing=c("grey","#FBD2A5","#FF8888","#C17F9D")
# 
# for (k in c("PL_expo","fmax_psd",
#             "flow_length","moran_I","Small_patches",
#             "core_area","mean_psd")){
#   
#   d_data_out=read.table(paste0("../Data/Linear_models_factor_cover_control/Keep_data/Data_",k,"_FALSE_aridity_all.csv"),sep=" ")
#   
#   model_spa_stat=readRDS(paste0("../Data/Linear_models_factor_cover_control/Keep_models/Mod_",k,"_FALSE_aridity_all.rds"))
#   
#   d_data_out$Grazing = as.factor(d_data_out$Grazing)
#   d_data_out$Grazing = relevel(d_data_out$Grazing,ref = "0")
#   
#   partial_res=visreg::visreg(fit = model_spa_stat,xvar="Aridity",by="Grazing",plot=F)
#   
#   for (grazing_id in unique(partial_res$res$Grazing)){
#     summary_lm=(lm(data = dplyr::filter(partial_res$res,Grazing==grazing_id),visregRes~Aridity))
#     
#     d=rbind(d,tibble(Stat=k,
#                      q1=confint(summary_lm,level = .95)[2,1],
#                      q3=confint(summary_lm,level = .95)[2,2],
#                      q2=confint(summary_lm,level = 0)[2,2],
#                      q3_90=confint(summary_lm,level = .9)[2,2],
#                      q1_90=confint(summary_lm,level = .9)[2,1],
#                      ID_grazing=grazing_id))
#     
#     
#   }
# }
# 
# 
# ggplot(d%>%Rename_spatial_statistics(.)%>%
#          add_column(., Type_stat=sapply(1:nrow(.),function(x){
#            if (.$Stat[x] %in% c("Power-law exp. of the PSD","log (largest patch)",
#                                 "Mean patch size","Number of smallest patches")){
#              return("Patch-size")
#            }else if (.$Stat[x] =="Bare soil connectivity"){
#              return("Hydro.")
#            }else if (.$Stat[x] =="Spatial autocorrelation of vege."){
#              return("Autocor.")
#            }else{
#              return("Geom.")
#            }
#          })))+
#   geom_vline(xintercept = 0,linetype=9)+
#   
#   geom_linerange(aes(x=q2,y=Stat,xmin=q1_90,xmax=q3_90,color=as.factor(ID_grazing)),
#                  lwd=1.5,position=position_jitterdodge(seed=123))+   
#   geom_linerange(aes(x=q2,y=Stat,xmin=q1,xmax=q3,color=as.factor(ID_grazing)),
#                  lwd=.5,position=position_jitterdodge(seed=123))+  
#   geom_pointrange(aes(x=q2,y=Stat,xmin=q1,xmax=q3,fill=as.factor(ID_grazing)),color="transparent",
#                   size=.15,position=position_jitterdodge(seed=123),shape=21)+  
#   
#   the_theme2+
#   scale_color_manual(values=c("grey","#FBD2A5","#FF8888","#C17F9D"),
#                      labels=c("Ungrazed","Low grazing","Medium grazing","High grazing"))+
#   scale_fill_manual(values=c("grey20","grey20","grey20","grey20"))+
#   scale_shape_manual(values=c(21,16))+
#   labs(x=substitute(paste(beta," (partial residuals)")),y="",color="")+
#   guides(shape="none",fill="none")+
#   facet_grid(Type_stat ~ ., space = "free_y", 
#              scale = "free_y")+
#   theme(panel.background = element_rect(fill="white"),
#         legend.text = element_text(size=8))+
#   guides(color=guide_legend(ncol=4))+
#   theme(legend.position="bottom")

## >> Interaction herbivores and grazing ----

d=tibble()

for (k in c("PL_expo","fmax_psd",
            "flow_length","moran_I","Small_patches",
            "core_area","mean_psd")){
  
  d_data_out=read.table(paste0("../Data/Linear_models_factor_cover_control/Keep_data/Data_",k,"_FALSE_aridity_all.csv"),sep=" ")
  
  model_spa_stat=readRDS(paste0("../Data/Linear_models_factor_cover_control/Keep_models/Mod_",k,"_FALSE_aridity_all.rds"))
  
  d_data_out$Grazing = as.factor(d_data_out$Grazing)
  d_data_out$Grazing = relevel(d_data_out$Grazing,ref = "0")
  
  partial_res=visreg::visreg(fit = model_spa_stat,xvar="Grazing",by="Herbivores",plot=F)
  partial_res2=visreg::visreg(fit = model_spa_stat,xvar="Grazing",by="Herbivores",plot=F,alpha = .1)
  
  d=rbind(d,tibble(Stat=k,
                   q1=partial_res$fit$visregLwr,
                   q3=partial_res$fit$visregUpr,
                   q2=partial_res$fit$visregFit,
                   q1_90=partial_res2$fit$visregLwr,
                   q3_90=partial_res2$fit$visregUpr,
                   Grazing_id=partial_res$fit$Grazing,
                   Herbivores_id=partial_res$fit$Herbivores))
  
}


p=ggplot(d%>%Rename_spatial_statistics(.)%>%
           add_column(., Type_stat=sapply(1:nrow(.),function(x){
             if (.$Stat[x] %in% c("Power-law exp. of the PSD","log (largest patch)",
                                  "Mean patch size","Number of smallest patches")){
               return("Patch-size")
             }else if (.$Stat[x] =="Bare soil connectivity"){
               return("Hydro.")
             }else if (.$Stat[x] =="Spatial autocorrelation of vege."){
               return("Autocor.")
             }else{
               return("Geom.")
             }
           })))+
  geom_vline(xintercept = 0,linetype=9)+
  
  # geom_linerange(aes(y=q2,x=Grazing_id,ymin=q1_90,ymax=q3_90,color=as.factor(Herbivores_id)),
  #                lwd=1.5,position=position_jitterdodge(seed=123))+   
  geom_line(aes(y=q2,x=Grazing_id,color=as.factor(Herbivores_id),group=interaction(Stat,Herbivores_id)))+
  geom_pointrange(aes(y=q2,x=Grazing_id,ymin=q1,ymax=q3,fill=as.factor(Herbivores_id),color=as.factor(Herbivores_id)),
                  position=position_jitterdodge(seed=123))+  
  the_theme2+
  scale_fill_manual(values=c("grey20","grey20","grey20","grey20"))+
  scale_shape_manual(values=c(21,16))+
  labs(x=substitute(paste(beta," (partial residuals, relative to cattle)")),y="",color="")+
  guides(shape="none",fill="none")+
  geom_hline(yintercept = 0)+
  theme(panel.background = element_rect(fill="white"),
        legend.text = element_text(size=8))+
  guides(color=guide_legend(ncol=4))+
  scale_color_manual(values=c("lightgreen","orange","violet","lightblue"))+
  facet_wrap(. ~ Stat, scales = "free")


## >> Boxplot variables and grazing ----

save_data=read.table("../Data/Spatial_structure_with_traits.csv",sep=";")

pdf("../Figures/SI/Exploration_drivers_sp_struc.pdf",width = 10,height = 4)
for (spatial_metric in colnames(save_data)[c(6,8,9,11,14,17,19,96,97)]){
  
  
  #We first control for all covariates and extract the residuals
  model_lmer=lm("value ~ Long_cos + Long_sin + Lattitude + Elevation",
                data = save_data%>%melt(., measure.vars=spatial_metric)%>%dplyr::filter(., !is.na(value)),
                na.action = na.fail)
  
  resid_model=residuals(model_lmer) #extract residuals of the distance to the tipping point after controlling for all covariates
  
  d_sem=save_data[as.numeric(names(resid_model)),]%>% #add it to the dataframe 
    add_column(., Resid_mod=resid_model)%>%
    mutate(., Grazing==as.factor(Grazing))
  d_sem$Grazing=as.character(d_sem$Grazing)
  
  
  
  p1=ggplot(d_sem)+
    geom_point(aes(x=Woody,y=Resid_mod,fill=Grazing,group=Grazing),shape=21,size=2)+
    scale_fill_manual(values=c("0"="green","1"="#FBD2A5","2"="#FF8888","3"="#C17F9D"))+
    the_theme2
  
  p2=ggplot(d_sem)+
    geom_point(aes(x=Herb,y=Resid_mod,fill=Grazing,group=Grazing),shape=21,size=2)+
    scale_fill_manual(values=c("0"="green","1"="#FBD2A5","2"="#FF8888","3"="#C17F9D"))+
    the_theme2
  
  p3=ggplot(d_sem)+
    geom_point(aes(x=Grass,y=Resid_mod,fill=Grazing,group=Grazing),shape=21,size=2)+
    scale_fill_manual(values=c("0"="green","1"="#FBD2A5","2"="#FF8888","3"="#C17F9D"))+
    the_theme2
  print(ggarrange(p1+ggtitle(spatial_metric),p2+ggtitle(spatial_metric),p3+ggtitle(spatial_metric),ncol=3))
  
  
  p1=ggplot(d_sem)+
    geom_point(aes(x=Total_N,y=Resid_mod,fill=Grazing,group=Grazing),shape=21,size=2)+
    scale_fill_manual(values=c("0"="green","1"="#FBD2A5","2"="#FF8888","3"="#C17F9D"))+
    the_theme2
  
  p2=ggplot(d_sem)+
    geom_point(aes(x=Org_C_v,y=Resid_mod,fill=Grazing,group=Grazing),shape=21,size=2)+
    scale_fill_manual(values=c("0"="green","1"="#FBD2A5","2"="#FF8888","3"="#C17F9D"))+
    the_theme2
  
  p3=ggplot(d_sem)+
    geom_point(aes(x=Org_C_tot,y=Resid_mod,fill=Grazing,group=Grazing),shape=21,size=2)+
    scale_fill_manual(values=c("0"="green","1"="#FBD2A5","2"="#FF8888","3"="#C17F9D"))+
    the_theme2
  print(ggarrange(p1,p2,p3,ncol=3))
  
  
  p1=ggplot(d_sem)+
    geom_point(aes(x=Productivity,y=Resid_mod,fill=Grazing,group=Grazing),shape=21,size=2)+
    scale_fill_manual(values=c("0"="green","1"="#FBD2A5","2"="#FF8888","3"="#C17F9D"))+
    the_theme2
  
  p2=ggplot(d_sem)+
    geom_point(aes(x=Fertility,y=Resid_mod,fill=Grazing,group=Grazing),shape=21,size=2)+
    scale_fill_manual(values=c("0"="green","1"="#FBD2A5","2"="#FF8888","3"="#C17F9D"))+
    the_theme2
  
  p3=ggplot(d_sem)+
    geom_point(aes(x=Forage_Quality,y=Resid_mod,fill=Grazing,group=Grazing),shape=21,size=2)+
    scale_fill_manual(values=c("0"="green","1"="#FBD2A5","2"="#FF8888","3"="#C17F9D"))+
    the_theme2
  print(ggarrange(p1,p2,p3,ncol=3))
  
  
  p1=ggplot(d_sem)+
    geom_point(aes(x=Dev_LDMC,y=Resid_mod,fill=Grazing,group=Grazing),shape=21,size=2)+
    scale_fill_manual(values=c("0"="green","1"="#FBD2A5","2"="#FF8888","3"="#C17F9D"))+
    the_theme2
  
  p2=ggplot(d_sem)+
    geom_point(aes(x=Dev_LA,y=Resid_mod,fill=Grazing,group=Grazing),shape=21,size=2)+
    scale_fill_manual(values=c("0"="green","1"="#FBD2A5","2"="#FF8888","3"="#C17F9D"))+
    the_theme2
  
  p3=ggplot(d_sem)+
    geom_point(aes(x=Dev_SLA,y=Resid_mod,fill=Grazing,group=Grazing),shape=21,size=2)+
    scale_fill_manual(values=c("0"="green","1"="#FBD2A5","2"="#FF8888","3"="#C17F9D"))+
    the_theme2
  print(ggarrange(p1,p2,p3,ncol=3))
  
  
  p1=ggplot(d_sem)+
    geom_point(aes(x=Dev_MaxH,y=Resid_mod,fill=Grazing,group=Grazing),shape=21,size=2)+
    scale_fill_manual(values=c("0"="green","1"="#FBD2A5","2"="#FF8888","3"="#C17F9D"))+
    the_theme2
  
  p2=ggplot(d_sem)+
    geom_point(aes(x=Dev_LNC,y=Resid_mod,fill=Grazing,group=Grazing),shape=21,size=2)+
    scale_fill_manual(values=c("0"="green","1"="#FBD2A5","2"="#FF8888","3"="#C17F9D"))+
    the_theme2
  
  p3=ggplot(d_sem)+
    geom_point(aes(x=Dev_LCC,y=Resid_mod,fill=Grazing,group=Grazing),shape=21,size=2)+
    scale_fill_manual(values=c("0"="green","1"="#FBD2A5","2"="#FF8888","3"="#C17F9D"))+
    the_theme2
  print(ggarrange(p1,p2,p3,ncol=3))
  
  
  
  
  p1=ggplot(d_sem)+
    geom_point(aes(x=Dev_PC1,y=Resid_mod,fill=Grazing,group=Grazing),shape=21,size=2)+
    scale_fill_manual(values=c("0"="green","1"="#FBD2A5","2"="#FF8888","3"="#C17F9D"))+
    the_theme2
  
  p2=ggplot(d_sem)+
    geom_point(aes(x=Dev_PC2,y=Resid_mod,fill=Grazing,group=Grazing),shape=21,size=2)+
    scale_fill_manual(values=c("0"="green","1"="#FBD2A5","2"="#FF8888","3"="#C17F9D"))+
    the_theme2
  
  p3=ggplot(d_sem)+
    geom_point(aes(x=Dev_PC3,y=Resid_mod,fill=Grazing,group=Grazing),shape=21,size=2)+
    scale_fill_manual(values=c("0"="green","1"="#FBD2A5","2"="#FF8888","3"="#C17F9D"))+
    the_theme2
  print(ggarrange(p1,p2,p3,ncol=3))
}
dev.off()


pdf("../Figures/SI/EXploration_grazing_effects_mecanisms.pdf.",width = 5,height = 4)
{d_sem=Closer_to_normality_CWM(save_data)
  d_sem$Grazing=as.character(d_sem$Grazing)
  spatial_metric=""
  
  print(ggplot(d_sem)+
          geom_boxplot(aes(x=Grazing,y=Woody,fill=Grazing,group=Grazing))+
          scale_fill_manual(values=c("0"="green","1"="#FBD2A5","2"="#FF8888","3"="#C17F9D"))+
          the_theme2+ggtitle(spatial_metric))
  print(ggplot(d_sem)+
          geom_boxplot(aes(x=Grazing,y=Grass,fill=Grazing,group=Grazing))+
          scale_fill_manual(values=c("0"="green","1"="#FBD2A5","2"="#FF8888","3"="#C17F9D"))+
          the_theme2+ggtitle(spatial_metric))
  print(ggplot(d_sem)+
          geom_boxplot(aes(x=Grazing,y=Herb,fill=Grazing,group=Grazing))+
          scale_fill_manual(values=c("0"="green","1"="#FBD2A5","2"="#FF8888","3"="#C17F9D"))+
          the_theme2+ggtitle(spatial_metric))
  print(ggplot(d_sem)+
          geom_boxplot(aes(x=Grazing,y=Dev_LA,fill=Grazing,group=Grazing))+
          scale_fill_manual(values=c("0"="green","1"="#FBD2A5","2"="#FF8888","3"="#C17F9D"))+
          the_theme2+ggtitle(spatial_metric))
  print(ggplot(d_sem)+
          geom_boxplot(aes(x=Grazing,y=Dev_LL,fill=Grazing,group=Grazing))+
          scale_fill_manual(values=c("0"="green","1"="#FBD2A5","2"="#FF8888","3"="#C17F9D"))+
          the_theme2+ggtitle(spatial_metric))
  print(ggplot(d_sem)+
          geom_boxplot(aes(x=Grazing,y=Dev_SLA,fill=Grazing,group=Grazing))+
          scale_fill_manual(values=c("0"="green","1"="#FBD2A5","2"="#FF8888","3"="#C17F9D"))+
          the_theme2+ggtitle(spatial_metric))
  print(ggplot(d_sem)+
          geom_boxplot(aes(x=Grazing,y=Dev_LCC,fill=Grazing,group=Grazing))+
          scale_fill_manual(values=c("0"="green","1"="#FBD2A5","2"="#FF8888","3"="#C17F9D"))+
          the_theme2+ggtitle(spatial_metric))
  print(ggplot(d_sem)+
          geom_boxplot(aes(x=Grazing,y=Dev_LDMC,fill=Grazing,group=Grazing))+
          scale_fill_manual(values=c("0"="green","1"="#FBD2A5","2"="#FF8888","3"="#C17F9D"))+
          the_theme2+ggtitle(spatial_metric))
  print(ggplot(d_sem)+
          geom_boxplot(aes(x=Grazing,y=Dev_MaxH,fill=Grazing,group=Grazing))+
          scale_fill_manual(values=c("0"="green","1"="#FBD2A5","2"="#FF8888","3"="#C17F9D"))+
          the_theme2+ggtitle(spatial_metric))
  print(ggplot(d_sem)+
          geom_boxplot(aes(x=Grazing,y=Dev_Phenol,fill=Grazing,group=Grazing))+
          scale_fill_manual(values=c("0"="green","1"="#FBD2A5","2"="#FF8888","3"="#C17F9D"))+
          the_theme2+ggtitle(spatial_metric))
  print(ggplot(d_sem)+
          geom_boxplot(aes(x=Grazing,y=Dev_SLA,fill=Grazing,group=Grazing))+
          scale_fill_manual(values=c("0"="green","1"="#FBD2A5","2"="#FF8888","3"="#C17F9D"))+
          the_theme2+ggtitle(spatial_metric))
  print(ggplot(d_sem)+
          geom_boxplot(aes(x=Grazing,y=Dev_PC1,fill=Grazing,group=Grazing))+
          scale_fill_manual(values=c("0"="green","1"="#FBD2A5","2"="#FF8888","3"="#C17F9D"))+
          the_theme2+ggtitle(spatial_metric))
  
  print(ggplot(d_sem)+
          geom_boxplot(aes(x=Grazing,y=Dev_PC2,fill=Grazing,group=Grazing))+
          scale_fill_manual(values=c("0"="green","1"="#FBD2A5","2"="#FF8888","3"="#C17F9D"))+
          the_theme2+ggtitle(spatial_metric))
  print(ggplot(d_sem)+
          geom_boxplot(aes(x=Grazing,y=Dev_PC3,fill=Grazing,group=Grazing))+
          scale_fill_manual(values=c("0"="green","1"="#FBD2A5","2"="#FF8888","3"="#C17F9D"))+
          the_theme2+ggtitle(spatial_metric))
  
  
  
  
  print(ggplot(d_sem)+
          geom_boxplot(aes(x=Grazing,y=CWM_LA,fill=Grazing,group=Grazing))+
          scale_fill_manual(values=c("0"="green","1"="#FBD2A5","2"="#FF8888","3"="#C17F9D"))+
          the_theme2+ggtitle(spatial_metric))
  print(ggplot(d_sem)+
          geom_boxplot(aes(x=Grazing,y=CWM_LL,fill=Grazing,group=Grazing))+
          scale_fill_manual(values=c("0"="green","1"="#FBD2A5","2"="#FF8888","3"="#C17F9D"))+
          the_theme2+ggtitle(spatial_metric))
  print(ggplot(d_sem)+
          geom_boxplot(aes(x=Grazing,y=CWM_SLA,fill=Grazing,group=Grazing))+
          scale_fill_manual(values=c("0"="green","1"="#FBD2A5","2"="#FF8888","3"="#C17F9D"))+
          the_theme2+ggtitle(spatial_metric))
  print(ggplot(d_sem)+
          geom_boxplot(aes(x=Grazing,y=CWM_LCC,fill=Grazing,group=Grazing))+
          scale_fill_manual(values=c("0"="green","1"="#FBD2A5","2"="#FF8888","3"="#C17F9D"))+
          the_theme2+ggtitle(spatial_metric))
  print(ggplot(d_sem)+
          geom_boxplot(aes(x=Grazing,y=CWM_LDMC,fill=Grazing,group=Grazing))+
          scale_fill_manual(values=c("0"="green","1"="#FBD2A5","2"="#FF8888","3"="#C17F9D"))+
          the_theme2+ggtitle(spatial_metric))
  print(ggplot(d_sem)+
          geom_boxplot(aes(x=Grazing,y=CWM_MaxH,fill=Grazing,group=Grazing))+
          scale_fill_manual(values=c("0"="green","1"="#FBD2A5","2"="#FF8888","3"="#C17F9D"))+
          the_theme2+ggtitle(spatial_metric))
  print(ggplot(d_sem)+
          geom_boxplot(aes(x=Grazing,y=CWM_Phenol,fill=Grazing,group=Grazing))+
          scale_fill_manual(values=c("0"="green","1"="#FBD2A5","2"="#FF8888","3"="#C17F9D"))+
          the_theme2+ggtitle(spatial_metric))
  print(ggplot(d_sem)+
          geom_boxplot(aes(x=Grazing,y=CWM_SLA,fill=Grazing,group=Grazing))+
          scale_fill_manual(values=c("0"="green","1"="#FBD2A5","2"="#FF8888","3"="#C17F9D"))+
          the_theme2+ggtitle(spatial_metric))
  print(ggplot(d_sem)+
          geom_boxplot(aes(x=Grazing,y=CWM_PC1,fill=Grazing,group=Grazing))+
          scale_fill_manual(values=c("0"="green","1"="#FBD2A5","2"="#FF8888","3"="#C17F9D"))+
          the_theme2+ggtitle(spatial_metric))
  
  print(ggplot(d_sem)+
          geom_boxplot(aes(x=Grazing,y=CWM_PC2,fill=Grazing,group=Grazing))+
          scale_fill_manual(values=c("0"="green","1"="#FBD2A5","2"="#FF8888","3"="#C17F9D"))+
          the_theme2+ggtitle(spatial_metric))
  print(ggplot(d_sem)+
          geom_boxplot(aes(x=Grazing,y=CWM_PC3,fill=Grazing,group=Grazing))+
          scale_fill_manual(values=c("0"="green","1"="#FBD2A5","2"="#FF8888","3"="#C17F9D"))+
          the_theme2+ggtitle(spatial_metric))
  
  
  print(ggplot(d_sem)+
          geom_boxplot(aes(x=Grazing,y=Fertility,fill=Grazing,group=Grazing))+
          scale_fill_manual(values=c("0"="green","1"="#FBD2A5","2"="#FF8888","3"="#C17F9D"))+
          the_theme2+ggtitle(spatial_metric))
  print(ggplot(d_sem)+
          geom_boxplot(aes(x=Grazing,y=Forage_Quality,fill=Grazing,group=Grazing))+
          scale_fill_manual(values=c("0"="green","1"="#FBD2A5","2"="#FF8888","3"="#C17F9D"))+
          the_theme2+ggtitle(spatial_metric))
  
}
dev.off()
## >> Basic plots traits ----

d_biodesert_traits=read.table("../Data/Biodesert_all_sites_traits.csv",sep=";",quote = "")
colnames(d_biodesert_traits)=gsub("X.","",colnames(d_biodesert_traits))

p=ggplot(d_biodesert_traits)+
  geom_point(aes(x=CWM_MaxLS.,y=Woody.),shape=21,size=2,fill="lightblue",color="black")+
  geom_smooth(aes(x=CWM_MaxLS.,y=Woody.),fill="lightblue",color="black",method = "lm")+
  the_theme2+
  labs(x="CWM Lateral spread",y="% of Woody")+
  ylim(0,1)


ggsave("../Figures/SI/Woody_MaxLS_spatial_structure.pdf",p,width = 5,height = 4)


## >> Shrubland woody connectivity
d_data=Perform_PCA_spatial_struc(Closer_to_normality_CWM(read.table("../Data/Spatial_structure_control_cover_with_traits.csv",sep=";"))%>%
                                   mutate(., Grazing=recode_factor(Grazing,"0"="Ungrazed","1"="Low","2"="Medium","3"="High")),plot = F)

p1=ggplot(d_biodesert_traits)+
  geom_boxplot(aes(VEG.,Woody.))+
  the_theme2+
  labs(x="Type of vegetation",y="% of Woody")+
  scale_x_discrete(labels=c("Grassland","Open woodlands \n with shrubs",
                            "Shrubland","Forest"))

p2=ggplot(d_data)+
  geom_boxplot(aes(Type_veg,flow_length))+
  the_theme2+
  labs(x="Type of vegetation",y="Bare soil connectivity")+
  scale_x_discrete(labels=c("Forest","Open woodlands \n with shrubs","Grassland",
                            "Shrubland"))


p14=ggplot(d_data)+
  geom_point(aes(x=Woody,y=flow_length),shape=21,size=2,fill="lightblue",color="black")+
  geom_smooth(aes(x=Woody,y=flow_length),fill="lightblue",color="black",method = "lm")+
  the_theme2+
  labs(x="% of Woody",y="Bare soil connectivity")

p24=ggplot(d_data)+
  geom_point(aes(x=CWM_MaxLS,y=flow_length),shape=21,size=2,fill="lightblue",color="black")+
  geom_smooth(aes(x=CWM_MaxLS,y=flow_length),fill="lightblue",color="black",method = "lm")+
  the_theme2+
  labs(x="CWM Lateral spread",y="Bare soil connectivity")


p34=ggplot(d_data)+
  geom_point(aes(x=-Dev_MaxLS,y=flow_length),shape=21,size=2,fill="lightblue",color="black")+
  geom_smooth(aes(x=-Dev_MaxLS,y=flow_length),fill="lightblue",color="black",method = "lm")+
  the_theme2+
  labs(x="Trait aggregation: lateral spread",y="Bare soil connectivity")

p3=ggarrange(p14,p24,p34,ncol=3,labels = letters[2:4])

ggsave("../Figures/SI/Type_vegetation.pdf",ggarrange(ggarrange(p1+theme(axis.text.x = element_text(hjust = 1,angle = 60)),
                                                               p2+theme(axis.text.x = element_text(hjust = 1,angle = 60)),
                                                               ncol=2,labels = letters[1:2]),
                                                     p3,
                                                     nrow=2,heights = c(1,.7),labels = c("","c"))
       ,width = 10,height = 6)




## >> Simple plots traits and grazing/aridity ----


p1=ggplot(d_biodesert_traits%>%
            mutate(., Dev_MaxLS.=-Dev_MaxLS.)%>%
            melt(., measure.vars = c("Woody.","CWM_MaxLS.","Dev_MaxLS."))%>%
            mutate(., variable=recode_factor(variable,
                                             "Woody."="% of Woody",
                                             "Dev_MaxLS."="Trait spatial aggregation: \n Lateral spread",
                                             #"Dev_LDMC"="Dev. LDMC from rand. expec.",
                                             #"Dev_Phenol"="Dev. Phenolics from rand. expec.",
                                             "CWM_MaxLS."="CWM Lateral spread")))+
  geom_boxplot(aes(x=Grazing.,y=value,group=as.factor(Grazing.)),outlier.shape = NA)+
  geom_jitter(aes(x=Grazing.,y=value,color=as.factor(Grazing.)),size=.5,alpha=.5)+
  scale_color_manual(values=c("grey","#FBD2A5","#FF8888","#C17F9D"))+
  facet_wrap(.~variable,scales = "free",ncol = 5)+
  the_theme2+
  labs(x="Grazing intensity",y="")+
  theme(legend.position = "none")

p2=ggplot(d_biodesert_traits%>%
            mutate(., Dev_MaxLS.=-Dev_MaxLS.)%>%
            melt(., measure.vars = c("Woody.","CWM_MaxLS.","Dev_MaxLS."))%>%
            mutate(., variable=recode_factor(variable,
                                             "Woody."="% of Woody",
                                             "Dev_MaxLS."="Trait spatial aggregation: \n Lateral spread",
                                             #"Dev_LDMC"="Dev. LDMC from rand. expec.",
                                             #"Dev_Phenol"="Dev. Phenolics from rand. expec.",
                                             "CWM_MaxLS."="CWM Lateral spread")))+
  geom_point(aes(x=Aridity.,y=value),size=1,alpha=.5,color="grey")+
  facet_wrap(.~variable,scales = "free",ncol = 5)+
  the_theme2+
  labs(x="Aridity level",y="")+
  theme(legend.position = "none")

ggsave("../Figures/SI/Distribution_traits_grazing_aridity.pdf",
       ggarrange(p1,p2,nrow=2,labels = letters[1:2]),width = 8,height = 5)



