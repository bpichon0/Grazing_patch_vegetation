rm(list=ls())
source("./Structure_grazing_function.R")

# ---------------- Step 0: Preliminary analyses ----



## >> Distributions of metrics ----
d_data=read.table("../Data/Spatial_structure_grazing.csv",sep=";")

p=ggplot(d_data%>%
           melt(., id.vars=c("Site_ID","Sub_id","Full_name","Resolution","Type_veg","Nurse","Herbivores")))+
  geom_histogram(aes(x=value),fill=alpha("gray40",.75))+
  facet_wrap(.~variable,scales = "free")+
  the_theme
ggsave("../Figures/Step0_preliminary/Distrib_metrics.pdf",p,width = 8,height = 10)




## >> Correlation between metrics ----
d_data=read.table("../Data/Spatial_structure_grazing.csv",sep=";")

mat_cor=cor(d_data[,c(4:24)],use = "na.or.complete",method = "pearson")
mat_cor[upper.tri(mat_cor)]=NA
diag(mat_cor)=NA

p1=ggplot(mat_cor%>%melt(.))+
  geom_tile(aes(x=Var1,y=Var2,fill=value))+
  scale_fill_gradient2(low="red","white","blue")+
  the_theme+
  theme(axis.text.x = element_text(angle=60,hjust=1),axis.title = element_blank())

#then, correlation graph

mat_cor=abs(rcorr(as.matrix(d_data[,c(4:24)]),type = "pearson")$r)
mat_cor[rcorr(as.matrix(d_data[,c(4:24)]),type = "pearson")$P>.05]=0

network_EWS = igraph::simplify(
  graph_from_adjacency_matrix(mat_cor,
                              weighted = TRUE,
                              mode = c("undirected")))

E(network_EWS)$width = E(network_EWS)$weight*5

# Finding the best modularity partition

n_sim=500
all_modularity=sapply(1:n_sim,function(x){
  set.seed(x)
  clust = cluster_louvain(network_EWS) 
  V(network_EWS)$community=clust$membership
  return(modularity(network_EWS, clust$membership,weights = E(network_EWS)))
})

set.seed(sample(which(all_modularity==max(all_modularity))))
clust = cluster_louvain(network_EWS) 
V(network_EWS)$community=clust$membership
mod = modularity(network_EWS, clust$membership,weights = E(network_EWS))
layout_new=layout.fruchterman.reingold(network_EWS,weights=E(network_EWS)$weight)

color_modules=c("#A684DE","#72CCC3","#E6967F","#B9E4AC","#8E354F")

pdf("../Figures/Step0_preliminary/Correlations_EWS.pdf",width = 7,height = 7)
#print(p1)
plot(network_EWS,
     vertex.color = color_modules[V(network_EWS)$community],
     vertex.frame.width = 1,vertex.label.color="black",
     edge.curved = .3,edge.color = "grey",
     layout = layout_new,frame = TRUE
)
dev.off()



#correlation of metrics with cover
p=ggplot(d_data%>%melt(., measure.vars=colnames(d_data)[4:29]))+
  geom_point(aes(x=rho_p,y=value),color="grey")+
  the_theme+
  facet_wrap(.~variable,scales = "free")+
  labs(x="Vegetation cover",y="Spatial metrics")

ggsave("../Figures/Step0_preliminary/Correlation_metrics_cover.pdf",p,width = 12,height = 10)


## >> Slopes change in metric with grazing & aridity ----

dir.create("../Figures/Step0_preliminary/Trends/",showWarnings = F)

#Getting the trend in model
d_sim=read.table("../Data/Spatial_structure_simulations.csv",sep=";")
d_trend_mod=tibble()
for (k in colnames(d_sim)[4:29]){ d_trend_mod=rbind(d_trend_mod,Compute_trend_model(d_sim,k))}

#Getting the trend in data
d_data=read.table("../Data/Spatial_structure_grazing.csv",sep=";")
stat_id=c(4:9,11:29)
d_data[,c(stat_id,35)]=apply(d_data[,c(stat_id,35)],2,z_tranform) #standardize variables

d_trend_data=d_signif_grazing=d_signif_interaction=tibble()
for (k in colnames(d_data)[stat_id]){ 
  d_trend_data=rbind(d_trend_data,Compute_trend_data(d_data,k))
  d_signif_grazing=rbind(d_signif_grazing,Test_grazing_effect(d_data,k))
  d_signif_interaction=rbind(d_signif_interaction,Test_interaction_effect(d_data,k))
}


pA=ggplot(
  rbind(
    d_trend_data%>%add_column(., Type="Data")%>%dplyr::select(., -p_val_full,),
    d_trend_mod%>%add_column(., Type="Model"))%>%
    filter(., Stat %in% unique(.$Stat)[1:10],Driver=="Aridity"))+
  
  geom_pointrange(aes(x=q2,xmin=q1,xmax=q3,y=Stat,color=q2>0,alpha=p_val_indiv<0.05))+
  geom_vline(xintercept = 0,linetype=9)+
  geom_hline(yintercept = seq(.5,15.5,5),linetype=9)+
  
  scale_color_manual(values = rev(c("red","blue")))+
  scale_alpha_manual(values=c(.2,1))+
  facet_grid(.~Type,scales = "free")+
  the_theme+
  labs(x="Slope against aridity",y="Spatial statistics",color="",alpha="")+
  guides(color=F,alpha=F)

pB=ggplot(
  rbind(
    d_trend_data%>%add_column(., Type="Data")%>%dplyr::select(., -p_val_full,),
    d_trend_mod%>%add_column(., Type="Model"))%>%
    filter(., Stat %in% unique(.$Stat)[1:10],Driver=="Grazing"))+
  
  geom_pointrange(aes(x=q2,xmin=q1,xmax=q3,y=Stat,color=q2>0,alpha=p_val_indiv<0.05))+
  geom_vline(xintercept = 0,linetype=9)+
  geom_hline(yintercept = seq(.5,15.5,5),linetype=9)+
  geom_text(data=d_signif_grazing%>%
              add_column(., Type="Data")%>%
              filter(., p_val_full<.05,Stat %in% unique(.$Stat)[1:10]),
            aes(y=Stat,x=.5),label="*",size=6)+
  
  scale_color_manual(values = rev(c("red","blue")))+
  scale_alpha_manual(values=c(.2,1))+
  facet_grid(.~Type,scales = "free")+
  the_theme+
  labs(x="Slope against grazing",y="Spatial statistics",color="",alpha="")+
  guides(color=F,alpha=F)

ggsave("../Figures/Step0_preliminary/Trends/Trends_data_models.pdf",ggarrange(pA,pB,nrow=2),width=6,height = 6)






pA=ggplot(
  rbind(
    d_trend_data%>%add_column(., Type="Data")%>%dplyr::select(., -p_val_full,),
    d_trend_mod%>%add_column(., Type="Model"))%>%
    filter(., Driver=="Aridity",Stat %!in% unique(.$Stat)[1:10]))+
  
  geom_pointrange(aes(x=q2,xmin=q1,xmax=q3,y=Stat,color=q2>0,alpha=p_val_indiv<0.05))+
  geom_vline(xintercept = 0,linetype=9)+
  geom_hline(yintercept = seq(.5,15.5,5),linetype=9)+
  
  scale_color_manual(values = rev(c("red","blue")))+
  scale_alpha_manual(values=c(.2,1))+
  facet_grid(.~Type,scales = "free")+
  the_theme+
  labs(x="Slope against aridity",y="Spatial statistics",color="",alpha="")+
  guides(color=F,alpha=F)

pB=ggplot(
  rbind(
    d_trend_data%>%add_column(., Type="Data")%>%dplyr::select(., -p_val_full),
    d_trend_mod%>%add_column(., Type="Model"))%>%
    filter(., Driver=="Grazing",Stat %!in% unique(.$Stat)[1:10]))+
  
  geom_pointrange(aes(x=q2,xmin=q1,xmax=q3,y=Stat,color=q2>0,alpha=p_val_indiv<0.05))+
  geom_vline(xintercept = 0,linetype=9)+
  geom_hline(yintercept = seq(.5,15.5,5),linetype=9)+
  geom_text(data=d_signif_grazing%>%
              add_column(., Type="Data")%>%
              filter(., p_val_full<.05,Stat %!in% unique(.$Stat)[1:10]),
            aes(y=Stat,x=.5),label="*",size=6)+
  
  scale_color_manual(values = rev(c("red","blue")))+
  scale_alpha_manual(values=c(.2,1))+
  facet_grid(.~Type,scales = "free")+
  the_theme+
  labs(x="Slope against grazing",y="Spatial statistics",color="",alpha="")+
  guides(color=F,alpha=F)

ggsave("../Figures/Step0_preliminary/Trends/Trends_data_models2.pdf",ggarrange(pA,pB,nrow=2),width=6,height = 12)


## >> Metrics grazing and different herbivores  ----

dir.create("../Figures/Step0_preliminary/Grazing_herb/",showWarnings = F)

d_data=read.table("../Data/Spatial_structure_grazing.csv",sep=";")

p1=ggplot(NULL)+
  geom_jitter(data=d_data%>%
                melt(.,measure.vars=colnames(d_data)[4:14]),
              aes(x=Grazing,y=value),width=.1,color="gray",alpha=.4)+
  geom_pointrange(data = d_data%>%
                    melt(., measure.vars=colnames(d_data)[4:14])%>%
                    group_by(., Grazing,variable)%>%
                    dplyr::summarise(., med_q2=median(value,na.rm=T),
                                     q1=quantile(value,.25,na.rm=T),
                                     q3=quantile(value,.75,na.rm=T),.groups = "keep"),
                  aes(x=Grazing,ymin=q1,ymax=q3,y=med_q2))+
  facet_wrap(.~variable,scales = "free")+
  the_theme

ggsave("../Figures/Step0_preliminary/Grazing_herb/Stats_grazing.pdf",p1,width = 8,height = 6)

p2=ggplot(NULL)+
  geom_jitter(data=d_data%>%
                melt(.,measure.vars=colnames(d_data)[15:29]),
              aes(x=Grazing,y=value),width=.1,color="gray",alpha=.4)+
  geom_pointrange(data = d_data%>%
                    melt(., measure.vars=colnames(d_data)[15:29])%>%
                    group_by(., Grazing,variable)%>%
                    dplyr::summarise(., med_q2=median(value,na.rm=T),
                                     q1=quantile(value,.25,na.rm=T),
                                     q3=quantile(value,.75,na.rm=T),.groups = "keep"),
                  aes(x=Grazing,ymin=q1,ymax=q3,y=med_q2))+
  facet_wrap(.~variable,scales = "free")+
  the_theme

ggsave("../Figures/Step0_preliminary/Grazing_herb/Stats_grazing2.pdf",p2,width = 8,height = 9)


p3=ggplot(NULL)+
  geom_jitter(data=d_data%>%filter(., Herbivores %in% c("Cattle","Goat","Horse","Sheep"))%>%
                melt(.,measure.vars=colnames(d_data)[4:14]),
              aes(x=Herbivores,y=value),width=.1,color="gray",alpha=.4)+
  geom_pointrange(data = d_data%>%filter(., Herbivores %in% c("Cattle","Goat","Horse","Sheep"))%>%
                    melt(., measure.vars=colnames(d_data)[4:14])%>%
                    group_by(., Herbivores,variable)%>%
                    dplyr::summarise(., med_q2=median(value,na.rm=T),
                                     q1=quantile(value,.25,na.rm=T),
                                     q3=quantile(value,.75,na.rm=T),.groups = "keep"),
                  aes(x=Herbivores,ymin=q1,ymax=q3,y=med_q2))+
  facet_wrap(.~variable,scales = "free")+
  the_theme

ggsave("../Figures/Step0_preliminary/Grazing_herb/Stats_herbivores.pdf",p3,width = 8,height = 6)

p4=ggplot(NULL)+
  geom_jitter(data=d_data%>%filter(., Herbivores %in% c("Cattle","Goat","Horse","Sheep"))%>%
                melt(.,measure.vars=colnames(d_data)[15:29]),
              aes(x=Herbivores,y=value),width=.1,color="gray",alpha=.4)+
  geom_pointrange(data = d_data%>%filter(., Herbivores %in% c("Cattle","Goat","Horse","Sheep"))%>%
                    melt(., measure.vars=colnames(d_data)[15:29])%>%
                    group_by(., Herbivores,variable)%>%
                    dplyr::summarise(., med_q2=median(value,na.rm=T),
                                     q1=quantile(value,.25,na.rm=T),
                                     q3=quantile(value,.75,na.rm=T),.groups = "keep"),
                  aes(x=Herbivores,ymin=q1,ymax=q3,y=med_q2))+
  facet_wrap(.~variable,scales = "free")+
  the_theme

ggsave("../Figures/Step0_preliminary/Grazing_herb/Stats_herbivores2.pdf",p4,width = 8,height = 9)





## >> Interaction grazing-aridity ----

d_data=read.table("../Data/Spatial_structure_grazing.csv",sep=";")

d_trend_data=tibble()
for (k in colnames(d_data)[4:14]){
  d_trend_data=rbind(d_trend_data,do.call(rbind,lapply(unique(d_data$Grazing),function(x){
    d2=melt(d_data,measure.vars=k)%>%filter(., Grazing==x)
    model_mix=lmer("value ~ Aridity + (1|Site_ID)",data=d2%>%filter(., !is.na(value)))
    sumary_model=summary(model_mix)
    d=data.frame(Stat_value=sumary_model$coefficients[1,1] + d2$Aridity*sumary_model$coefficients[2,1],
                 Aridity=d2$Aridity,Grazing=x,Stat=k)
    
    return(d)
  }))
  )
}



p=ggplot(d_trend_data)+
  geom_point(data=d_data%>%
               melt(., measure.vars=colnames(d_data)[4:14])%>%
               dplyr::rename(., Stat=variable),
             aes(x=Aridity,y=value,color=as.factor(Grazing)),alpha=.5)+
  geom_line(aes(x=Aridity,y=Stat_value,color=as.factor(Grazing)),lwd=1.5)+
  theme_classic()+theme(legend.position = "bottom")+
  facet_wrap(.~Stat,scales = "free")+
  scale_color_manual(values=c("#F7F597","#EF943C","#C57348","#AD0B0B"))+
  the_theme+
  labs(color="Grazing intensity",y="Spatial statistics")

ggsave("../Figures/Step0_preliminary/Trends/Interaction_aridity_grazing.pdf",p,width = 8,height = 6)



d_data=read.table("../Data/Spatial_structure_grazing.csv",sep=";")

d_trend_data=tibble()
for (k in colnames(d_data)[15:29]){
  d_trend_data=rbind(d_trend_data,do.call(rbind,lapply(unique(d_data$Grazing),function(x){
    d2=melt(d_data,measure.vars=k)%>%filter(., Grazing==x)
    model_mix=lmer("value ~ Aridity + (1|Site_ID)",data=d2%>%filter(., !is.na(value)))
    sumary_model=summary(model_mix)
    d=data.frame(Stat_value=sumary_model$coefficients[1,1] + d2$Aridity*sumary_model$coefficients[2,1],
                 Aridity=d2$Aridity,Grazing=x,Stat=k)
    
    return(d)
  }))
  )
}



p=ggplot(d_trend_data)+
  geom_point(data=d_data%>%
               melt(., measure.vars=colnames(d_data)[15:29])%>%
               dplyr::rename(., Stat=variable),
             aes(x=Aridity,y=value,color=as.factor(Grazing)),alpha=.5)+
  geom_line(aes(x=Aridity,y=Stat_value,color=as.factor(Grazing)),lwd=1.5)+
  theme_classic()+theme(legend.position = "bottom")+
  facet_wrap(.~Stat,scales = "free")+
  scale_color_manual(values=c("#F7F597","#EF943C","#C57348","#AD0B0B"))+
  the_theme+
  labs(color="Grazing intensity",y="Spatial statistics")

ggsave("../Figures/Step0_preliminary/Trends/Interaction_aridity_grazing2.pdf",p,width = 8,height = 9)


## >> Interaction herbivores-aridity ----

d_data=read.table("../Data/Spatial_structure_grazing.csv",sep=";")

d_trend_data=tibble()
for (k in colnames(d_data)[4:14]){
  d_trend_data=rbind(d_trend_data,do.call(rbind,lapply(c("Cattle","Goat","Horse","Sheep"),function(x){
    d2=melt(d_data,measure.vars=k)%>%filter(., Herbivores==x)
    model_mix=lmer("value ~ Aridity + (1|Site_ID)",data=d2%>%filter(., !is.na(value)))
    sumary_model=summary(model_mix)
    d=data.frame(Stat_value=sumary_model$coefficients[1,1] + d2$Aridity*sumary_model$coefficients[2,1],
                 Aridity=d2$Aridity,Herbivores=x,Stat=k)
    
    return(d)
  }))
  )
}


p=ggplot(d_trend_data)+
  geom_point(data=d_data%>%
               filter(.,Herbivores %in% c("Cattle","Goat","Horse","Sheep"))%>%
               melt(., measure.vars=colnames(d_data)[4:14])%>%
               dplyr::rename(., Stat=variable),
             aes(x=Aridity,y=value,color=as.factor(Herbivores)),alpha=.5)+
  geom_line(aes(x=Aridity,y=Stat_value,color=as.factor(Herbivores)),lwd=1.5)+
  theme_classic()+theme(legend.position = "bottom")+
  facet_wrap(.~Stat,scales = "free")+
  scale_color_manual(values=c("#A55EAD","#93C18D","#E0AD7C","#7C87E0"))+
  the_theme+
  labs(color="Type of herbivores",y="Spatial statistics")
ggsave("../Figures/Step0_preliminary/Trends/Interaction_aridity_herbivores.pdf",p,width = 8,height = 6)



d_data=read.table("../Data/Spatial_structure_grazing.csv",sep=";")

d_trend_data=tibble()
for (k in colnames(d_data)[15:29]){
  d_trend_data=rbind(d_trend_data,do.call(rbind,lapply(c("Cattle","Goat","Horse","Sheep"),function(x){
    d2=melt(d_data,measure.vars=k)%>%filter(., Herbivores==x)
    model_mix=lmer("value ~ Aridity + (1|Site_ID)",data=d2%>%filter(., !is.na(value)))
    sumary_model=summary(model_mix)
    d=data.frame(Stat_value=sumary_model$coefficients[1,1] + d2$Aridity*sumary_model$coefficients[2,1],
                 Aridity=d2$Aridity,Herbivores=x,Stat=k)
    
    return(d)
  }))
  )
}


p=ggplot(d_trend_data)+
  geom_point(data=d_data%>%
               filter(.,Herbivores %in% c("Cattle","Goat","Horse","Sheep"))%>%
               melt(., measure.vars=colnames(d_data)[15:29])%>%
               dplyr::rename(., Stat=variable),
             aes(x=Aridity,y=value,color=as.factor(Herbivores)),alpha=.5)+
  geom_line(aes(x=Aridity,y=Stat_value,color=as.factor(Herbivores)),lwd=1.5)+
  theme_classic()+theme(legend.position = "bottom")+
  facet_wrap(.~Stat,scales = "free")+
  scale_color_manual(values=c("#A55EAD","#93C18D","#E0AD7C","#7C87E0"))+
  the_theme+
  labs(color="Type of herbivores",y="Spatial statistics")
ggsave("../Figures/Step0_preliminary/Trends/Interaction_aridity_herbivores2.pdf",p,width = 8,height = 9)


## >> Inference resilience ----
dir.create("../Figures/Step0_preliminary/Inferrence/",showWarnings = F)
#x-y stats

x_y_stat=read.table("../Data/Inferrence/x_y_stat_all_100.csv",sep=";")
list_plots=list()
name_plot=c("Cover","# neighbors","Clustering","Skewness","Variance","Autocorrelation",
            "SDR","PLR","Exponent p.l.","CV PSD","Frac. max")
index=1
for (i in c(1:11)){
  d_fil=cbind(filter(x_y_stat%>%
                       melt(., id.vars=c("Site_ID","Method", "Type")),variable==colnames(x_y_stat)[i])%>%
                filter(., Type=="Sim")%>%dplyr::rename(., value_sim=value),
              filter(x_y_stat%>%
                       melt(., id.vars=c("Site_ID","Method", "Type")),variable==colnames(x_y_stat)[i])%>%
                filter(., Type=="Obs")%>%dplyr::rename(., value_obs=value)%>%dplyr::select(., value_obs))
  
  list_plots[[index]]=ggplot(d_fil)+
    geom_point(aes(x=value_obs,y=value_sim),color="#96C3DC",alpha=.75)+the_theme+
    labs(x="",y="")+
    geom_abline(slope=1,intercept = 0,color="black")+
    ggtitle(name_plot[i])+
    theme(title = element_text(size=10))
  
  index=index+1
}

p=annotate_figure(ggarrange(plotlist=list_plots,ncol = 4,nrow = 3),
                  left=text_grob("Closest simulations",rot=90,color="black",size=15,face ="bold",vjust=1,family = "NewCenturySchoolbook"),
                  bottom = text_grob("Observed spatial statistic",color="black",size=15,face="bold",vjust=-1,family = "NewCenturySchoolbook"))
ggsave("../Figures/Step0_preliminary/Inferrence/x_y_inferrence_stats.pdf",p,width = 10,height = 8)

# Inferred parameters 

params=read.table("../Data/Inferrence/param_inferred.csv",sep=";")
d_data=read.table("../Data/Spatial_structure_grazing.csv",sep=";")

nsite=ncol(params)/3
param_infered=tibble(p2=apply(params[,1:nsite],2,median),
                     p1=apply(params[,1:nsite],2,quantile,.25),
                     p3=apply(params[,1:nsite],2,quantile,.75),
                     q2=apply(params[,(nsite+1):(2*nsite)],2,median),
                     q1=apply(params[,(nsite+1):(2*nsite)],2,quantile,.25),
                     q3=apply(params[,(nsite+1):(2*nsite)],2,quantile,.75),
)

d_data=cbind(d_data,param_infered)

p1_1=ggplot(d_data)+
  geom_pointrange(aes(x=Aridity,y=p2,ymax=p3,ymin=p1),color="gray",alpha=.5)+
  geom_smooth(aes(x=Aridity,y=p2),method = "lm",color="black")+
  the_theme+
  labs(y="Median of p posterior")

p1_2=ggplot(d_data)+
  geom_pointrange(aes(x=Aridity,y=q2,ymax=q3,ymin=q1),color="gray",alpha=.5)+
  geom_smooth(aes(x=Aridity,y=q2),method = "lm",color="black")+
  the_theme+
  labs(y="Median of q posterior")

p2_1=ggplot(d_data)+
  geom_jitter(aes(x=Grazing,y=p2),width=.15,color="gray",alpha=.5)+
  geom_pointrange(data=d_data%>%group_by(., Grazing)%>%
                    dplyr::summarise(., .groups = "keep",
                                     med_p=median(p2,na.rm=T),
                                     q1_p=quantile(p2,.25),
                                     q3_p=quantile(p2,.75)),
                  aes(x=Grazing,y=med_p,ymax=q3_p,ymin=q1_p))+
  the_theme+
  labs(y="Median of p posterior")

p2_2=ggplot(d_data)+
  geom_jitter(aes(x=Grazing,y=q2),width=.15,color="gray",alpha=.5)+
  geom_pointrange(data=d_data%>%group_by(., Grazing)%>%
                    dplyr::summarise(., .groups = "keep",
                                     med_q=median(q2,na.rm=T),
                                     q1_q=quantile(q2,.25),
                                     q3_q=quantile(q2,.75)),
                  aes(x=Grazing,y=med_q,ymax=q3_q,ymin=q1_q))+
  the_theme+
  labs(y="Median of q posterior")

p_tot=ggarrange(p1_1,p1_2,p2_1,p2_2,nrow=2,ncol=2,labels = c("",letters[1],"",letters[2]))

ggsave("../Figures/Step0_preliminary/Inferrence/Params_grazing_aridity.pdf",p_tot,width = 7,height = 7)


## >> PCA ----

d_pca=list(list(id=1,sumstat=c(1:14)),list(id=2,sumstat=c(1:3,15:20)),
           list(id=3,sumstat=c(1:3,21:29)))

for (k in 1:length(d_pca)){
  
  sumstat_name=colnames(d)[d_pca[[k]]$sumstat]
  res.comp=imputePCA(d[,which(colnames(d) %in% sumstat_name)],ncp=3,scale = T) 
  
  if ("completeObs" %in% names(res.comp)){
    res.pca=PCA(res.comp$completeObs, ncp = 3,  graph=F)
  }else {
    res.pca=PCA(res.comp, ncp = 3,  graph=F)
  }
  
  axes_for_plot=tibble(x=c(1,1,2),y=c(2,3,3))
  
  for (i in 1:3){
    assign(paste0("p",i),
           d%>%
             add_column(., PC1=res.pca$ind$coord[,axes_for_plot$x[i]],PC2=res.pca$ind$coord[,axes_for_plot$y[i]])%>%
             ggplot(.) +
             geom_hline(yintercept = 0, lty = 2) +
             geom_vline(xintercept = 0, lty = 2) +
             geom_point(aes(x = PC1, y = PC2, color = as.factor(Grazing),fill=as.factor(Grazing)))+
             scale_color_viridis_d()+
             scale_fill_viridis_d()+
             labs(x=paste0("PC ",axes_for_plot$x[i]," (",round(res.pca$eig[axes_for_plot$x[i],2], 1)," %)"),
                  y=paste0("PC ",axes_for_plot$y[i]," (",round(res.pca$eig[axes_for_plot$y[i],2], 1)," %)"),
                  color="Grazing intensity",fill="")+
             ggtitle("")+guides()+
             theme_classic()+theme(legend.position = "bottom")+
             guides(fill="none")+
             theme(legend.box = "vertical")
    )
  }
  
  p=ggarrange(ggarrange(p1+theme(legend.position = "none"),
                        p2+theme(legend.position = "none"),
                        p3+theme(legend.position = "none"),
                        ncol=3,align = "hv"),ggarrange(ggplot()+theme_void(),get_legend(p2),ggplot()+theme_void(),ncol=3,widths = c(.3,1,.3)),
              nrow=2,heights = c(1,.2))
  
  ggsave(paste0("../Figures/Step0_preliminary/PCA_grazing_spatial_structure_",d_pca[[k]]$id,".pdf"),p,width = 11,height = 5)
  
  #by type of vegetation
  for (i in 1:3){
    assign(paste0("p",i),
           d%>%
             add_column(., PC1=res.pca$ind$coord[,axes_for_plot$x[i]],PC2=res.pca$ind$coord[,axes_for_plot$y[i]])%>%
             ggplot(.) +
             geom_hline(yintercept = 0, lty = 2) +
             geom_vline(xintercept = 0, lty = 2) +
             geom_point(aes(x = PC1, y = PC2, color = as.factor(Type_veg),fill=as.factor(Type_veg)))+
             labs(x=paste0("PC ",axes_for_plot$x[i]," (",round(res.pca$eig[axes_for_plot$x[i],2], 1)," %)"),
                  y=paste0("PC ",axes_for_plot$y[i]," (",round(res.pca$eig[axes_for_plot$y[i],2], 1)," %)"),
                  color="Type vegetation",fill="")+
             ggtitle("")+guides()+
             theme_classic()+theme(legend.position = "bottom")+
             guides(fill="none")+
             theme(legend.box = "vertical")
    )
  }
  
  p=ggarrange(ggarrange(p1+theme(legend.position = "none"),
                        p2+theme(legend.position = "none"),
                        p3+theme(legend.position = "none"),
                        ncol=3,align = "hv"),ggarrange(ggplot()+theme_void(),get_legend(p2),ggplot()+theme_void(),ncol=3,widths = c(.3,1,.3)),
              nrow=2,heights = c(1,.2))
  
  ggsave(paste0("../Figures/Step0_preliminary/PCA_grazing_spatial_structure_",d_pca[[k]]$id,"_veg.pdf"),p,width = 11,height = 5)
  
  #by type of herbivores
  
  
  
  for (i in 1:3){
    assign(paste0("p",i),
           d%>%
             add_column(., PC1=res.pca$ind$coord[,axes_for_plot$x[i]],PC2=res.pca$ind$coord[,axes_for_plot$y[i]])%>%
             ggplot(.) +
             geom_hline(yintercept = 0, lty = 2) +
             geom_vline(xintercept = 0, lty = 2) +
             geom_point(aes(x = PC1, y = PC2, color = as.factor(Herbivores),fill=as.factor(Herbivores)))+
             labs(x=paste0("PC ",axes_for_plot$x[i]," (",round(res.pca$eig[axes_for_plot$x[i],2], 1)," %)"),
                  y=paste0("PC ",axes_for_plot$y[i]," (",round(res.pca$eig[axes_for_plot$y[i],2], 1)," %)"),
                  color="Type vegetation",fill="")+
             ggtitle("")+guides()+
             theme_classic()+theme(legend.position = "bottom")+
             guides(fill="none")+
             theme(legend.box = "vertical")
    )
  }
  
  p=ggarrange(ggarrange(p1+theme(legend.position = "none"),
                        p2+theme(legend.position = "none"),
                        p3+theme(legend.position = "none"),
                        ncol=3,align = "hv"),ggarrange(ggplot()+theme_void(),get_legend(p2),ggplot()+theme_void(),ncol=3,widths = c(.3,1,.3)),
              nrow=2,heights = c(1,.2))
  
  ggsave(paste0("../Figures/Step0_preliminary/PCA_grazing_spatial_structure_",d_pca[[k]]$id,"_herbivore.pdf"),p,width = 11,height = 5)
  
  
}





## >> Type of plants ----

d_data=d_data%>%
  add_column(.,
             Fern=sapply(1:nrow(.),function(x){return(d_biodesert$RFCQ[which(d_biodesert$ID==.$Site_ID[x])])}),
             Grass=sapply(1:nrow(.),function(x){return(d_biodesert$RGCQ[which(d_biodesert$ID==.$Site_ID[x])])}),
             Herb=sapply(1:nrow(.),function(x){return(d_biodesert$RHCQ[which(d_biodesert$ID==.$Site_ID[x])])}),
             Woody=sapply(1:nrow(.),function(x){return(d_biodesert$RWCQ[which(d_biodesert$ID==.$Site_ID[x])])})
  )

ggplot(d_data%>%melt(., measure.vars = colnames(d_data)[52:55])%>%
         dplyr::group_by(variable,Grazing)%>%
         dplyr::summarise(., .groups = "keep",mean_val=mean(value),sd_val=sd(value)))+
  geom_pointrange(aes(x=Grazing,y=mean_val,ymin=mean_val-sd_val,ymax=mean_val+sd_val))+
  the_theme+
  facet_wrap(.~variable,scales = "free")



# ---------------- Step 1: Understanding grazing on vegetation cover----

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
select_model=dredge(model_cover, subset = ~ Slope & Elevation &
                      Long_cos & Long_sin & Lattitude &
                      dc(Woody & Grazing, Woody : Grazing) &
                      dc(Clim1 & Grazing, Clim1 : Grazing) &
                      dc(Clim2 & Grazing, Clim2 : Grazing) &
                      dc(Clim3 & Grazing, Clim3 : Grazing) &
                      dc(Clim4 & Grazing, Clim4 : Grazing) &
                      dc(rho_p & Grazing, rho_p : Grazing) &
                      dc(Org_C & Grazing, Org_C : Grazing) &
                      dc(Sand  & Grazing, Sand  : Grazing) &
                      dc(Type_veg & Grazing, Type_veg : Grazing),
                    options(na.action = "na.fail") )

#best model
model_cover=MuMIn::get.models(select_model, subset = delta == 0)[[1]]
saveRDS(model_cover,paste0("../Data/Step1_Understanding_grazing/Keep_models/Mod_cover_TRUE.rds"))

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

ggsave(paste0("../Figures/Final_figs/SI/Predictors_COVER.pdf"),
       ggarrange(p1,ggarrange(ggarrange(p2,p3,p4+theme(legend.position="none"),ncol=3),
                              get_legend(p4),nrow=2,heights = c(1,.3)),
                 nrow=2,heights = c(1, .3)),
       width = 6,height = 9)





# ---------------- Step 2: Understanding grazing on spatial structure ----

## >> Importance ----
dir.create("../Figures/Step1_Understanding_grazing/Importance",showWarnings = F)

for (with_interactions in c("_nointer","","_grazing_square","_grazing_no_constrains")){
  
  d_all=read.table(paste0("../Data/Step1_Understanding_grazing/Importance",with_interactions,".csv"),sep=";")
  
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
  
  ggsave(paste0("../Figures/Step1_Understanding_grazing/Importance/Importance_grazing",with_interactions,".pdf"),p,width = 6,height = 5)
  
  
}

#Importance based on R2
d_R2=read.table(paste0("../Data/Step1_Understanding_grazing/Importance.csv"),sep=";")%>%
  filter(., With_cover==T)%>%
  dplyr::select(., Sp_stat,R2m)%>%
  dplyr::rename(., Stat=Sp_stat)%>%
  Rename_spatial_statistics(.)

d_all2=read.table(paste0("../Data/Step1_Understanding_grazing/Estimators_model.csv"),sep=";")%>%
  filter(., With_cover==T)%>%
  Organize_df(., "bar")%>%
  group_by(., Type_pred,Stat)%>%
  dplyr::summarise(., sum_effect=sum(abs(Median)),.groups = "keep")%>%
  Rename_spatial_statistics(.)%>%group_by(., Stat)%>%
  dplyr::summarise(., .groups="keep",sum_effect=sum_effect/sum(sum_effect),Type_pred=Type_pred)


p=ggplot(d_all2)+
  geom_bar(aes(x=100*sum_effect,y=Stat,fill=Type_pred),stat="identity",width = .75)+
  geom_text(data=d_R2,aes(y=Stat,x=8,label=paste0("r² = ",round(R2m,2))))+
  the_theme+
  labs(x="Importance (%)",fill="")+
  scale_fill_manual(values=c("Geography"="#FFB15B","Grazing"="#F3412A",
                             "Abiotic"="#4564CE","Vegetation"="#6DD275","Climatic"="#CC80C7"))+
  theme(axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank())

ggsave("../Figures/Step1_Understanding_grazing/Importance/Importance_grazing_effects.pdf",p,width = 6,height = 6)





## >> Direct & interactions: LME models ----

for (with_interactions in c("","_nointer","_grazing_square")){
  
  for (cover in c(1,0)){

    d_all=read.table(paste0("../Data/Step1_Understanding_grazing/Importance",with_interactions,".csv"),sep=";")%>%
      filter(., With_cover==cover)
    
    d_all2=read.table(paste0("../Data/Step1_Understanding_grazing/Estimators_model",with_interactions,".csv"),sep=";")
  
    for (k in unique(d_all2$Stat)){
      
      loc_pval=.05+max(d_all2%>% #localization of pvalues in the plot
                         filter(., Stat==k , With_cover==cover,term!="Type, both")%>%#not enougth data
                         Organize_df(., "predictor")%>%
                         dplyr::select(.,q3)%>%dplyr::pull(.))
      
      p1_1=ggplot(d_all2%>%
                    filter(., Stat==k , With_cover==cover,term!="Type, both")%>%
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
      
      p1_2=ggplot(d_all2%>%
                    filter(., Stat==k , With_cover==cover)%>%
                    add_column(., xplot=0)%>%
                    Organize_df(., "bar")%>%
                    group_by(., Type_pred,Stat,xplot)%>%
                    dplyr::summarise(., sum_effect=sum(abs(Median)),.groups = "keep")%>%
                    mutate(., sum_effect=sum_effect/sum(sum_effect)))+
        geom_bar(aes(x=xplot,y=sum_effect,fill=Type_pred),stat="identity",width = .05)+
        geom_text(aes(x=0,y=1.05,label=paste0("R² \n = \n ",round(d_all$R2m[which(d_all$Sp_stat==k)],2))),size=3)+
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
      
      ggsave(paste0("../Figures/Step1_Understanding_grazing/Pred",ifelse(cover,"_with_cover","_without_cover"),
                    "/Preditor_",k,with_interactions,".pdf"),p1,width = 6,height = 6)
    }
  }
}



# Saving all slopes against grazing and ploting it
d_all2=read.table("../Data/Step1_Understanding_grazing/Estimators_model.csv",sep=";")

p=ggplot(d_all2%>%
           filter(., term=="Grazing",With_cover==T)%>%
           Rename_spatial_statistics(.)%>%
           add_column(., Signif=sapply(1:nrow(.),function(x){
             return(is_signif(.$pvalue[x]))}))%>%
           arrange(., Median)%>%
           add_column(., Order_f=1:nrow(.))%>%
           mutate(.,Stat = fct_reorder(Stat, Order_f))%>%
           add_column(., order_signif=sapply(1:nrow(.),
                                             function(x){return(ifelse(.$Median[x]<0,.4,-.25))})))+
  geom_pointrange(aes(x=Median,y=Stat,xmin=q1,xmax=q3,color=Median>0,alpha=pvalue<.05),shape=17)+
  the_theme+
  scale_color_manual(values=c("#FB4343","#1C5ECC"))+
  scale_alpha_manual(values=c(.4,1))+
  geom_text(aes(x=order_signif,y=1:14,label=Signif))+
  geom_vline(xintercept = 0,linetype=9)+
  geom_hline(yintercept = c(3.5,6.5,9.5),linetype=9)+
  labs(x=substitute(paste(beta," (Grazing)")),y="",color="")+
  theme(legend.position = "none")

ggsave("../Figures/Step1_Understanding_grazing/Beta_coef_stats_grazing.pdf",p,width = 6,height = 4)


## >> Indirect effects: SEM ----

### SEM to link grazing to spatial patterns through the changes in facilitation, sand content 

d_data=read.table("../Data/Spatial_structure_grazing.csv",sep=";")%>%
  Closer_to_normality(.)

d_data[,c(1,9:29,35,37:41,44:ncol(d_data))] = apply(d_data[,c(1,9:29,35,37:41,44:ncol(d_data))],2,z_tranform)

#Getting the multifunctionality index
d_data$MF=z_tranform(rowSums(d_data[,55:(ncol(d_data))],na.rm = T))

#Adding the first two components from the PCA
d_data=Perform_PCA_spatial_struc(d_data)

d_indirect=tibble() #to save indirect effects of grazing on the spatial structure
dir.create("../Figures/Step1_Understanding_grazing/SEM",showWarnings = F)

param_list=expand.grid(Cov=c("with_cover","without_cover"),
                       Stat=c("perim_area_scaling","PL_expo","Cond_H","fmax_psd",
                              "PLR","flow_length","mean_perim_area","core_area_land",
                              "fractal_dim","division","contig","core_area",
                              "Struct1","Struct2"),
                       MF=c(T,F)) #MF instead of organic carbon

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
  
  if (param_list$MF[ID]){ save$Org_C=save$MF}
  
  #DOING the SEMs
  
  if (with_cover=="with_cover"){
    
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
    
    
  }else {
    
    #SEM with all grazing intensity
    
    d_sem=save
    all_d=summary(psem(
      lm(Resid_mod ~ Grazing + Woody , d_sem),
      lm(Woody ~ Grazing + Org_C + Sand, d_sem),
      lm(Sand ~ Grazing , d_sem),
      lm(Org_C ~ Grazing, d_sem)
    ))
    
    #SEM with low grazing intensity
    
    d_sem=filter(save,Grazing %in% 0:1)
    low_graz=summary(psem(
      lm(Resid_mod ~ Grazing + Woody , d_sem),
      lm(Woody ~ Grazing + Org_C + Sand, d_sem),
      lm(Sand ~ Grazing , d_sem),
      lm(Org_C ~ Grazing, d_sem)
    ))
    
    #SEM with high grazing intensity
    
    d_sem=filter(save,Grazing %in% 2:3)
    high_graz=summary(psem(
      lm(Resid_mod ~ Grazing + Woody , d_sem),
      lm(Woody ~ Grazing + Org_C + Sand, d_sem),
      lm(Sand ~ Grazing , d_sem),
      lm(Org_C ~ Grazing, d_sem)
    ))
  }
  
  
  d_indirect=rbind(d_indirect,
                   Get_indirect_effects_grazing(all_d)%>%
                     add_column(., Type="All",Stat=k,With_cover=with_cover,MF=param_list$MF[ID]),
                   Get_indirect_effects_grazing(high_graz)%>%
                     add_column(., Type="High",Stat=k,With_cover=with_cover,MF=param_list$MF[ID]),
                   Get_indirect_effects_grazing(low_graz)%>%
                     add_column(., Type="Low",Stat=k,With_cover=with_cover,MF=param_list$MF[ID]))
  
  #ploting
  
  if(with_cover=="with_cover"){
    
    Plot_SEM_with_cover(summary_sem = all_d,pdf_ = T,
                        name = paste0("../Figures/Step1_Understanding_grazing/SEM/",param_list$Stat[ID],"/",
                                      "SEM_",param_list$Stat[ID],"_",param_list$Cov[ID],"_MF_",param_list$MF[ID],"_all"),
                        title_ = paste0("SEM_",k,"_all"),
                        name_var = k,MF = param_list$MF[ID])
    
    Plot_SEM_with_cover(summary_sem = low_graz,pdf_ = T,
                        name = paste0("../Figures/Step1_Understanding_grazing/SEM/",param_list$Stat[ID],"/",
                                      "SEM_",param_list$Stat[ID],"_",param_list$Cov[ID],"_MF_",param_list$MF[ID],"_low"),
                        title_ = paste0("SEM_",k,"_low"),
                        name_var = k,MF = param_list$MF[ID])
    
    Plot_SEM_with_cover(summary_sem = high_graz,pdf_ = T,
                        name = paste0("../Figures/Step1_Understanding_grazing/SEM/",param_list$Stat[ID],"/",
                                      "SEM_",param_list$Stat[ID],"_",param_list$Cov[ID],"_MF_",param_list$MF[ID],"_high"),
                        title_ = paste0("SEM_",k,"_high"),
                        name_var = k,MF = param_list$MF[ID])
  }else {
    
    Plot_SEM_without_cover(summary_sem = all_d,pdf_ = T,
                           name = paste0("../Figures/Step1_Understanding_grazing/SEM/",param_list$Stat[ID],"/",
                                         "SEM_",param_list$Stat[ID],"_",param_list$Cov[ID],"_MF_",param_list$MF[ID],"_all"),
                           title_ = paste0("SEM_",k,"_all"),
                           name_var = k,MF = param_list$MF[ID])
    
    Plot_SEM_without_cover(summary_sem = low_graz,pdf_ = T,
                           name = paste0("../Figures/Step1_Understanding_grazing/SEM/",param_list$Stat[ID],"/",
                                         "SEM_",param_list$Stat[ID],"_",param_list$Cov[ID],"_MF_",param_list$MF[ID],"_low"),
                           title_ = paste0("SEM_",k,"_low"),
                           name_var = k,MF = param_list$MF[ID])
    
    Plot_SEM_without_cover(summary_sem = high_graz,pdf_ = T,
                           name = paste0("../Figures/Step1_Understanding_grazing/SEM/",param_list$Stat[ID],"/",
                                         "SEM_",param_list$Stat[ID],"_",param_list$Cov[ID],"_MF_",param_list$MF[ID],"_high"),
                           title_ = paste0("SEM_",k,"_high"),
                           name_var = k,MF = param_list$MF[ID])
  }
}

p=ggplot(d_indirect%>%
           filter(., With_cover=="with_cover"))+
  geom_bar(aes(x=Path,y=Indirect,fill=Type),position="dodge", stat="identity",width = .5)+
  facet_wrap(.~Stat,scales = "free",nrow = 5)+
  scale_fill_manual(values=c("#A1CA7A","#BB7CB8","#D29271"))+
  the_theme+
  geom_hline(yintercept = 0)

ggsave("../Figures/Step1_Understanding_grazing/Indirect_effects_grazing.pdf",p,width = 7,height = 9)


## >> Residuals analysis: interactions ----

dir.create("../Figures/Step1_Understanding_grazing/Residuals",showWarnings = F)
for (k in list.files("../Data/Step1_Understanding_grazing/Keep_models/",".rds")
     [-grep("cover",list.files("../Data/Step1_Understanding_grazing/Keep_models/",".rds"))]){
  
  stat=gsub("_TRUE.rds","",gsub("Mod_","",k))
  with_cover=T
  
  d_data=read.table("../Data/Spatial_structure_grazing.csv",sep=";")%>%
    Closer_to_normality(.)
  
  d_data[,c(1,9:29,35,37:41,44:ncol(d_data))] = 
    apply(d_data[,c(1,9:29,35,37:41,44:ncol(d_data))],2,z_tranform)
  
  d_data=Perform_PCA_spatial_struc(d_data)
  
  d_data_mod=d_data%>%melt(., measure.vars=stat)%>%
    filter(., !is.na(value))
  
  
  if (with_cover){
    
    formula_mod=formula(formula_(paste("value ~ Long_cos + Long_sin + Lattitude + Slope + Elevation 
      + Clim1 + Clim2 + Clim3 + Clim4 + Grazing 
      + Woody + Woody*Grazing
      + Clim1*Grazing + Clim2*Grazing + Clim3*Grazing + Clim4*Grazing  
      + rho_p + rho_p*Grazing + Type_veg + Type_veg*Grazing
      + Sand + Sand*Grazing + Org_C + Org_C*Grazing
      + (1|Site_ID)")))
    
  } else{
    
    formula_mod=formula(formula_(paste("value ~ Long_cos + Long_sin + Lattitude + Slope + Elevation 
      + Clim1 + Clim2 + Clim3 + Clim4 + Grazing 
      + Woody + Woody*Grazing
      + Clim1*Grazing + Clim2*Grazing + Clim3*Grazing + Clim4*Grazing  
      + Type_veg + Type_veg*Grazing
      + Sand + Sand*Grazing + Org_C + Org_C*Grazing
      + (1|Site_ID)")))
  }
  
  model_spa_stat  = lmer(formula_mod, d_data_mod,
                         na.action = na.fail ,REML ="FALSE")
  
  #we remove potential outliers
  rm.outliers = romr.fnc(model_spa_stat, d_data_mod, trim=2.5)
  d_data_out = rm.outliers$data
  
  model_spa_stat = lmer(formula_mod, data = d_data_out, 
                        na.action = na.fail,REML ="FALSE")
  
  model_spa_stat=readRDS(paste0("../Data/Step1_Understanding_grazing/Keep_models/",k))
  
  model_spa_stat = lmer(formula(model_spa_stat), data = d_data_out, 
                        na.action = na.fail,REML ="FALSE")
  
  for (metric in c("Clim2","Org_C","Woody","Clim1")){
    
    if (metric %in% rownames(summary(model_spa_stat)$coefficients)){
      
      resid_mod=visreg::visreg(fit = model_spa_stat,xvar = metric,by="Grazing",plot=F) #extracting residuals
      resid_mod$res$Grazing=rep(0:3,as.vector(table(d_data_out$Grazing))) #correcting vigreg output
      
      p=ggplot(resid_mod$res%>%
                 mutate(., Grazing=as.character(Grazing))%>%
                 dplyr::select(., -value)%>%
                 melt(., measure.vars=metric))+
        geom_point(aes(value,visregRes),color="black",alpha=.2)+
        geom_smooth(aes(value,visregRes,fill=Grazing,color=Grazing),method = "lm",level=.95)+
        the_theme+
        facet_wrap(.~Grazing,scales = "free",labeller = label_bquote(cols = "Grazing"==.(Grazing)),ncol = 4)+
        scale_color_manual(values=c("0"="#6EDA22",'1'="#FFF100","2"="#FBAC45","3"="#E82736"))+
        scale_fill_manual(values=c("0"="#6EDA22",'1'="#FFF100","2"="#FBAC45","3"="#E82736"))+
        labs(x=ifelse(metric=="Org_C","Facilitation",metric),y=paste0("Partial residuals of ",stat))+
        theme(legend.position = "none")
      
      ggsave(paste0("../Figures/Step1_Understanding_grazing/Residuals/",stat,"_",metric,".pdf"),p,width = 8,height = 3)
    }
    
  }
  
  if ("Grazing" %in% rownames(summary(model_spa_stat)$coefficients)){
    
    metric="Grazing"
    resid_mod=visreg::visreg(fit = model_spa_stat,xvar = "Grazing",plot=F) #extracting residuals
    
    p=ggplot(resid_mod$res%>%
               mutate(., Grazing=as.character(Grazing)))+
      geom_jitter(aes(Grazing,visregRes),color="black",alpha=.2,width = .1)+
      geom_pointrange(data=resid_mod$res%>%
                        group_by(., Grazing)%>%
                        mutate(., Grazing=as.character(Grazing))%>%
                        dplyr::summarise(., .groups = "keep",Mean_res=median(visregRes),q1=quantile(visregRes,.25),q3=quantile(visregRes,.75)),
                      aes(Grazing,y=Mean_res,ymin=q1,ymax=q3,fill=Grazing),color="black",shape=24,lwd=1)+
      the_theme+
      scale_color_manual(values=c("0"="#6EDA22",'1'="#FFF100","2"="#FBAC45","3"="#E82736"))+
      scale_fill_manual(values=c("0"="#6EDA22",'1'="#FFF100","2"="#FBAC45","3"="#E82736"))+
      labs(x=ifelse(metric=="Org_C","Facilitation",metric),y=paste0("Partial residuals of ",stat))+
      theme(legend.position = "none")
    
    ggsave(paste0("../Figures/Step1_Understanding_grazing/Residuals/",stat,"_",metric,".pdf"),p,width = 5,height = 3)
    
    
  }
  
}





## >> Residuals analysis: grazing effect ----

pdf("../Figures/Step1_Understanding_grazing/Residual_effect_grazing_cover.pdf",width = 6,height = 4)

d_slope=tibble()

d_data=read.table("../Data/Spatial_structure_grazing.csv",sep=";")%>%
  Closer_to_normality(.)
d_data[,c(1,9:29,35,37:41,44:ncol(d_data))] = 
  apply(d_data[,c(1,9:29,35,37:41,44:ncol(d_data))],2,z_tranform)
d_data=Perform_PCA_spatial_struc(d_data)

save=d_data

for (k in c("perim_area_scaling","fmax_psd","Cond_H","PL_expo","Spectral_ratio",
            "core_area_land","division","fractal_dim","contig","core_area",
            "flow_length","PLR",
            "Struct1","Struct2")){
  
  #for each we plot the slope against the partial residuals with and without cover
  
  par(mfrow=c(1,2))
  
  d_data_mod=save%>%melt(., measure.vars=k)%>%
    filter(., !is.na(value))
  
  formula_mod=formula(formula_(paste("value ~ Long_cos + Long_sin + Lattitude + Slope + Elevation 
      + Clim1 + Clim2 + Clim3 + Clim4 + Grazing 
      + Woody + Woody*Grazing
      + Clim1*Grazing + Clim2*Grazing + Clim3*Grazing + Clim4*Grazing  
      + rho_p + rho_p*Grazing + Type_veg + Type_veg*Grazing
      + Sand + Sand*Grazing + Org_C + Org_C*Grazing
      + (1|Site_ID)")))
  
  
  model_spa_stat  = lmer(formula_mod, d_data_mod,
                         na.action = na.fail ,REML ="FALSE")
  
  #we remove potential outliers
  rm.outliers = romr.fnc(model_spa_stat, d_data_mod, trim=2.5)
  d_data_out = rm.outliers$data
  
  model_spa_stat = lmer(formula_mod, data = d_data_out, 
                        na.action = na.fail,REML ="FALSE")
  
  # model_spa_stat=readRDS(paste0("../Data/Step1_Understanding_grazing/Keep_models/Mod_",k,"_TRUE.rds"))
  
  resid_mod=visreg::visreg(fit = model_spa_stat,xvar="Grazing",plot=T) 
  mtext(paste0(k,", with cover"))
  
  mod_cov=lm(visregRes~Grazing,resid_mod$res)
  
  d_slope=rbind(d_slope,
                tibble(pval=summary(mod_cov)$coefficient[2,4],
                       slope=summary(mod_cov)$coefficient[2,1],
                       Low_int=confint(mod_cov)[2,1],
                       High_int=confint(mod_cov)[2,2],
                       With_cover=T,
                       Stat=k
                ))
  
  
  formula_mod=formula(formula_(paste("value ~ Long_cos + Long_sin + Lattitude + Slope + Elevation 
      + Clim1 + Clim2 + Clim3 + Clim4 + Grazing 
      + Woody + Woody*Grazing
      + Clim1*Grazing + Clim2*Grazing + Clim3*Grazing + Clim4*Grazing  
      + Type_veg + Type_veg*Grazing
      + Sand + Sand*Grazing + Org_C + Org_C*Grazing
      + (1|Site_ID)")))
  
  
  model_spa_stat  = lmer(formula_mod, d_data_mod,
                         na.action = na.fail ,REML ="FALSE")
  
  #we remove potential outliers
  rm.outliers = romr.fnc(model_spa_stat, d_data_mod, trim=2.5)
  d_data_out = rm.outliers$data
  
  model_spa_stat = lmer(formula_mod, data = d_data_out, 
                        na.action = na.fail,REML ="FALSE")
  # model_spa_stat=readRDS(paste0("../Data/Step1_Understanding_grazing/Keep_models/Mod_",k,"_FALSE.rds"))
  
  resid_mod=visreg::visreg(fit = model_spa_stat,xvar="Grazing",plot=T) 
  mtext(paste0(k,", without cover"))
  
  mod_no_cov=lm(visregRes~Grazing,resid_mod$res)
  
  d_slope=rbind(d_slope,
                tibble(pval=summary(mod_no_cov)$coefficient[2,4],
                       slope=summary(mod_no_cov)$coefficient[2,1],
                       Low_int=confint(mod_no_cov)[2,1],
                       High_int=confint(mod_no_cov)[2,2],
                       With_cover=F,
                       Stat=k
                ))
}
dev.off()

write.table(d_slope,"../Data/Step1_Understanding_grazing/Slope_partial_residuals.csv",sep=";")




## >> RDA analysis ----

d_data=read.table("../Data/Spatial_structure_grazing.csv",sep=";")%>%
  Closer_to_normality(.)

d_data[,c(1,9:29,35,37:41,44:ncol(d_data))] = apply(d_data[,c(1,9:29,35,37:41,44:ncol(d_data))],2,z_tranform)

d_data$MF=rowSums(d_data[,55:(ncol(d_data))])

responses_var=c("fractal_dim","contig","perim_area_scaling","mean_perim_area",
                "core_area","core_area_land","fmax_psd","Cond_H","PLR","cv_psd",
                "Spectral_ratio","flow_length")

predictor_var=c("Aridity","Grazing")
controlled_var=c("Sand","Org_C","Long_cos",
                 "Long_sin","Lattitude","Slope","Elevation")

row_to_remove=unlist(sapply(1:nrow(d_data),function(x){
  if(any(is.na(d_data[x,responses_var]))){return(x)}}))

#Doing the RDA analysis
Y=as.matrix(d_data[-row_to_remove,responses_var])
X=as.matrix(d_data[-row_to_remove,predictor_var])
W=as.matrix(d_data[-row_to_remove,controlled_var])

RDA_spatial_struct=rda(Y ~ X + Condition(W))

data_rda1=as.data.frame(summary(RDA_spatial_struct)$species)

data_rda2=as.data.frame(scores(RDA_spatial_struct, display="bp", choices=c(1, 2), scaling=1))

#Here distance between 

p=ggplot(NULL)+ 
  geom_text(data=data_rda1%>%add_column(., labels=rownames(data_rda1)),
            aes(x=RDA1, y=RDA2,label=labels),size=4)+
  geom_segment(data=data_rda2, aes(x=0, xend=RDA1, y=0, yend=RDA2), 
               arrow=arrow(length=unit(0.05,"npc")),lwd=2,color="#23529E") +
  geom_text(data=data_rda2%>%add_column(., labels=gsub("X","",rownames(data_rda2))), 
            aes(x=RDA1,y=RDA2,label=labels,
                hjust=0.5*(1-sign(RDA1)),vjust=0.5*(1-sign(RDA2))), 
            color="#0C2F69", size=4)+
  geom_hline(yintercept=0, linetype=9)+
  geom_vline(xintercept=0, linetype=9)+
  the_theme+
  xlim(-.7,1)+
  labs(x=paste0("RDA 1 (",round(summary(RDA_spatial_struct)$cont$importance[1,1],2),")"),
       y=paste0("RDA 2 (",round(summary(RDA_spatial_struct)$cont$importance[1,2],2),")"))

ggsave("../Figures/Step1_Understanding_grazing/RDA_analysis.pdf",p,width = 5,height = 4)

#Testing significativity of the RDA
anova.cca(RDA_spatial_struct, permutations = 1000)
#the RDA is supported here, based on the permutation test


## >> SEM with aridity & grazing ----

d_data=read.table("../Data/Spatial_structure_grazing.csv",sep=";")%>%
  Closer_to_normality(.)
d_data[,c(1,9:29,35,37:41,44:ncol(d_data))] = apply(d_data[,c(1,9:29,35,37:41,44:ncol(d_data))],2,z_tranform)

#Getting the multifunctionality index
d_data$MF=z_tranform(rowSums(d_data[,55:(ncol(d_data))],na.rm = T))

#Adding the first two components from the PCA
d_data=Perform_PCA_spatial_struc(d_data)

dir.create("../Figures/Step1_Understanding_grazing/SEM_aritidy_grazing/",showWarnings = F)


for (k in c("perim_area_scaling","PL_expo","Cond_H","fmax_psd",
            "PLR","flow_length","mean_perim_area","core_area_land",
            "fractal_dim","division","contig","core_area",
            "Struct1","Struct2")){
  
  dir.create(paste0("../Figures/Step1_Understanding_grazing/SEM_aritidy_grazing/",k),showWarnings = F)
  
  model_lmer=lmer(paste("value ~ Long_cos + Long_sin + Lattitude + Slope + Elevation + Clim1 + Clim4 + Type_veg + (1|Site_ID)"),
                  data = d_data%>%melt(., measure.vars=k)%>%filter(., !is.na(value)),
                  na.action = na.omit,REML ="FALSE")
  
  #we remove potential outliers
  rm.outliers = romr.fnc(model_lmer, d_data%>%melt(., measure.vars=k)%>%filter(., !is.na(value)),
                         trim=2.5)
  d_data_out = rm.outliers$data
  
  #We first control for all covariates and extract the residuals
  model_lmer=lmer(paste("value ~ Long_cos + Long_sin + Lattitude + Slope + Elevation + Clim1 + Clim4 + Type_veg + (1|Site_ID)"),
                  data = d_data_out,
                  na.action = na.fail,REML ="FALSE")
  
  resid_model=residuals(model_lmer) #extract residuals
  
  save=d_data[as.numeric(names(resid_model)),]%>% #add it to the dataframe 
    add_column(., Resid_mod=resid_model)
  
  #DOING the SEMs
  
  d_sem=save
  all_d=summary(psem(
    lm(Resid_mod ~ MF + rho_p + Woody + Sp_richness + Sand + Grazing, d_sem),
    lm(Woody ~ Grazing + Aridity, d_sem),
    lm(rho_p ~ Grazing + Aridity, d_sem),
    lm(MF ~ Grazing + Aridity, d_sem),
    lm(Sp_richness ~ Grazing + Aridity, d_sem),
    lm(Sand ~ Grazing + Aridity, d_sem)
  ))
  
  #SEM with low grazing intensity
  
  d_sem=filter(save,Grazing %in% 0:1)%>%filter(., !is.na(N_transfo))
  low_graz=summary(psem(
    lm(Resid_mod ~ MF + rho_p + Woody + Sp_richness + Sand + Grazing, d_sem),
    lm(Woody ~ Grazing + Aridity, d_sem),
    lm(rho_p ~ Grazing + Aridity, d_sem),
    lm(MF ~ Grazing + Aridity, d_sem),
    lm(Sp_richness ~ Grazing + Aridity, d_sem),
    lm(Sand ~ Grazing + Aridity, d_sem)
  ))
  
  #SEM with high grazing intensity
  
  d_sem=filter(save,Grazing %in% 2:3)%>%filter(., !is.na(N_transfo))
  high_graz=summary(psem(
    lm(Resid_mod ~ MF + rho_p + Woody + Sp_richness + Sand + Grazing, d_sem),
    lm(Woody ~ Grazing + Aridity, d_sem),
    lm(rho_p ~ Grazing + Aridity, d_sem),
    lm(MF ~ Grazing + Aridity, d_sem),
    lm(Sp_richness ~ Grazing + Aridity, d_sem),
    lm(Sand ~ Grazing + Aridity, d_sem)
  ))
  
  
  #ploting
  Plot_SEM_Aridity_Grazing(summary_sem = all_d,pdf_ = T,
                      name = paste0("../Figures/Step1_Understanding_grazing/SEM_aritidy_grazing/",k,"/",
                                    "SEM_",k,"_all"),
                      title_ = paste0("SEM_",k,"_all"),
                      name_var = k)
  
  Plot_SEM_Aridity_Grazing(summary_sem = low_graz,pdf_ = T,
                      name = paste0("../Figures/Step1_Understanding_grazing/SEM_aritidy_grazing/",k,"/",
                                    "SEM_",k,"_low"),
                      title_ = paste0("SEM_",k,"_low"),
                      name_var = k)
  
  Plot_SEM_Aridity_Grazing(summary_sem = high_graz,pdf_ = T,
                      name = paste0("../Figures/Step1_Understanding_grazing/SEM_aritidy_grazing/",k,"/",
                                    "SEM_",k,"_high"),
                      title_ = paste0("SEM_",k,"_high"),
                      name_var = k)
  
}



## >> Lasso penalization for multivariate models ----

d_data=read.table("../Data/Spatial_structure_grazing.csv",sep=";")%>%
  Closer_to_normality(.)

d_data[,c(1,7:29,35,37:41,(44:(ncol(d_data))))] = apply(d_data[,c(1,7:29,35,37:41,(44:(ncol(d_data))))],2,z_tranform)

d_data$MF=rowSums(d_data[,55:(ncol(d_data))])

d_data=Perform_PCA_spatial_struc(d_data)

#response variables
responses_var=c("fractal_dim","contig","perim_area_scaling","mean_perim_area",
                "core_area","core_area_land","fmax_psd","Cond_H","PLR","cv_psd",
                "Spectral_ratio","flow_length","Struct1","Struct2")

predictor_var=c("Slope","Lattitude","Long_cos","Long_sin","Elevation",
                "Clim1","Clim2","Clim3","Clim4",
                "Grazing","rho_p","Woody",
                "Org_C","Sand")

row_to_remove=unlist(sapply(1:nrow(d_data),function(x){if(any(is.na(d_data[x,responses_var]))){return(x)}}))

#response matrix
Y=as.matrix(d_data[-row_to_remove,responses_var])

#preditors matrix
X=as.matrix(d_data[-row_to_remove,predictor_var])

# random_effect=as.vector(d_data$Site_ID)
# E_hat = residuals(lmer(Y ~ X + (1|random_effect),REML = "FALSE"))

#Getting residuals
E_hat = residuals(lm(Y ~ X,na.action = na.omit))

# Testing the type of error 
whitening_choice(E_hat, c( "nonparam","AR1","ARMA"), pAR = 1, qMA = 0)

square_root_inv_hat_Sigma = whitening(E_hat, "nonparam")

#Finding the most important variables for each predictor using Lasso penalization
#and reorganizing the df to sort predictors as decreasing importance/frequency

Frequencies=variable_selection(Y = Y,X = X,square_root_inv_hat_Sigma,nb_repli = 1000)%>%
  dplyr::rename(., term=Names_of_X,Stats=Names_of_Y)%>%
  dcast(., ... ~ Stats, value.var="frequency")%>%
  add_column(., Average_freq=rowMeans(.[,-1]))%>%
  arrange(., Average_freq)%>%
  add_column(., Order_f=rev(1:nrow(.)))%>%
  mutate(term = fct_reorder(term, Order_f))%>%
  add_column(., save_Average_freq=.$Average_freq)%>%
  mutate(., Average_freq=0)%>%
  melt(., id.vars=c("term","Order_f","save_Average_freq"))

p=ggplot(data = Frequencies%>%
           mutate(., value=sapply(1:nrow(.),function(x){if(.$value[x]<.5){return(NA)}else{return(.$value[x])}}))) +
  geom_tile(aes(x = term, y =variable, fill = value),size = 0.75) +
  geom_text(data=tibble(x=1:14,y=15,
                        label=rev(round(unique(Frequencies$save_Average_freq),2))),
            aes(x=x,y=y,label=label),size=3)+
  labs(x="",y="",fill="~ Importance")+
  scale_fill_gradientn(colours = colorRampPalette(c("white","#A129DA",alpha("#07299C",.3)))(100),na.value = "white")+
  scale_color_gradientn(colours = colorRampPalette(c("white","#A129DA",alpha("#07299C",.3)))(100),na.value = "white")+
  the_theme+
  theme(axis.text.x = element_text(angle = 60,hjust=1))

ggsave("../Figures/Step1_Understanding_grazing/Importance/Importance_lasso_penalized.pdf",p,width = 6,height = 5)


## >> OTHER 1/ SEM with many soil functioning variables (Not working well) ----

d_data=read.table("../Data/Spatial_structure_grazing.csv",sep=";")%>%
  Closer_to_normality(.)
d_data[,c(1,9:29,35,37:41,44:ncol(d_data))] = apply(d_data[,c(1,9:29,35,37:41,44:ncol(d_data))],2,z_tranform)

#Getting the multifunctionality index
d_data$MF=z_tranform(rowSums(d_data[,55:(ncol(d_data))],na.rm = T))

#Adding the first two components from the PCA
d_data=Perform_PCA_spatial_struc(d_data)

dir.create("../Figures/Step1_Understanding_grazing/SEM_all_soil/",showWarnings = F)


for (k in c("perim_area_scaling","PL_expo","Cond_H","fmax_psd",
            "PLR","flow_length","mean_perim_area","core_area_land",
            "fractal_dim","division","contig",
            "Struct1","Struct2")){
  
  dir.create(paste0("../Figures/Step1_Understanding_grazing/SEM_all_soil/",k),showWarnings = F)
  
  model_lmer=lmer(paste("value ~ Long_cos + Long_sin + Lattitude + Slope + Elevation + Clim1 + Clim2 + Clim3 + Clim4 + Type_veg + (1|Site_ID)",
                        ifelse(with_cover=="without_cover","+ rho_p","")),
                  data = d_data%>%melt(., measure.vars=k)%>%filter(., !is.na(value)),
                  na.action = na.omit,REML ="FALSE")
  
  #we remove potential outliers
  rm.outliers = romr.fnc(model_lmer, d_data%>%melt(., measure.vars=k)%>%filter(., !is.na(value)),
                         trim=2.5)
  d_data_out = rm.outliers$data
  
  #We first control for all covariates and extract the residuals
  model_lmer=lmer(paste("value ~ Long_cos + Long_sin + Lattitude + Slope + Elevation + Clim1 + Clim2 + Clim3 + Clim4 + Type_veg + (1|Site_ID)",
                        ifelse(with_cover=="without_cover","+ rho_p","")),
                  data = d_data_out,
                  na.action = na.fail,REML ="FALSE")
  
  resid_model=residuals(model_lmer) #extract residuals
  
  save=d_data[as.numeric(names(resid_model)),]%>% #add it to the dataframe 
    add_column(., Resid_mod=resid_model)
  
  #DOING the SEMs
  
  d_sem=save%>%filter(., !is.na(N_transfo))
  all_d=summary(psem(
    lm(Resid_mod ~ Grazing + rho_p + Woody , d_sem),
    lm(Woody ~ Grazing + Beta_gluco + Phosphatase + Amonium + Total_N + Nitrate + N_transfo +
         Hexose + N_minera + Total_P + Aromatic + Org_C + Sand, d_sem),
    lm(rho_p ~ Grazing + Beta_gluco + Phosphatase + Amonium + Total_N + Nitrate + N_transfo +
         Hexose + N_minera + Total_P + Aromatic + Org_C + Sand, d_sem),
    lm(Beta_gluco ~ Grazing , d_sem),
    lm(Phosphatase ~ Grazing , d_sem),
    lm(Amonium ~ Grazing , d_sem),
    lm(Total_N ~ Grazing , d_sem),
    lm(Nitrate ~ Grazing , d_sem),
    lm(N_transfo ~ Grazing , d_sem),
    lm(Hexose ~ Grazing , d_sem),
    lm(N_minera ~ Grazing , d_sem),
    lm(Aromatic ~ Grazing , d_sem),
    lm(Total_P ~ Grazing , d_sem),
    lm(Org_C ~ Grazing , d_sem),
    lm(Sand ~ Grazing , d_sem)
  ))
  
  #SEM with low grazing intensity
  
  d_sem=filter(save,Grazing %in% 0:1)%>%filter(., !is.na(N_transfo))
  low_graz=summary(psem(
    lm(Resid_mod ~ Grazing + rho_p + Woody , d_sem),
    lm(Woody ~ Grazing + Beta_gluco + Phosphatase + Amonium + Total_N + Nitrate + N_transfo +
         Hexose + N_minera + Total_P + Aromatic + Org_C + Sand, d_sem),
    lm(rho_p ~ Grazing + Beta_gluco + Phosphatase + Amonium + Total_N + Nitrate + N_transfo +
         Hexose + N_minera + Total_P + Aromatic + Org_C + Sand, d_sem),
    lm(Beta_gluco ~ Grazing , d_sem),
    lm(Phosphatase ~ Grazing , d_sem),
    lm(Amonium ~ Grazing , d_sem),
    lm(Total_N ~ Grazing , d_sem),
    lm(Nitrate ~ Grazing , d_sem),
    lm(N_transfo ~ Grazing , d_sem),
    lm(Hexose ~ Grazing , d_sem),
    lm(N_minera ~ Grazing , d_sem),
    lm(Aromatic ~ Grazing , d_sem),
    lm(Total_P ~ Grazing , d_sem),
    lm(Org_C ~ Grazing , d_sem),
    lm(Sand ~ Grazing , d_sem)
  ))
  
  #SEM with high grazing intensity
  
  d_sem=filter(save,Grazing %in% 2:3)%>%filter(., !is.na(N_transfo))
  high_graz=summary(psem(
    lm(Resid_mod ~ Grazing + rho_p + Woody , d_sem),
    lm(Woody ~ Grazing + Beta_gluco + Phosphatase + Amonium + Total_N + Nitrate + N_transfo +
         Hexose + N_minera + Total_P + Aromatic + Org_C + Sand, d_sem),
    lm(rho_p ~ Grazing + Beta_gluco + Phosphatase + Amonium + Total_N + Nitrate + N_transfo +
         Hexose + N_minera + Total_P + Aromatic + Org_C + Sand, d_sem),
    lm(Beta_gluco ~ Grazing , d_sem),
    lm(Phosphatase ~ Grazing , d_sem),
    lm(Amonium ~ Grazing , d_sem),
    lm(Total_N ~ Grazing , d_sem),
    lm(Nitrate ~ Grazing , d_sem),
    lm(N_transfo ~ Grazing , d_sem),
    lm(Hexose ~ Grazing , d_sem),
    lm(N_minera ~ Grazing , d_sem),
    lm(Aromatic ~ Grazing , d_sem),
    lm(Total_P ~ Grazing , d_sem),
    lm(Org_C ~ Grazing , d_sem),
    lm(Sand ~ Grazing , d_sem)
  ))
  
  
  #ploting
  
  if(with_cover=="with_cover"){
    
    Plot_SEM_with_cover(summary_sem = all_d,pdf_ = T,
                        name = paste0("../Figures/Step1_Understanding_grazing/SEM/",param_list$Stat[ID],"/",
                                      "SEM_",param_list$Stat[ID],"_",param_list$Cov[ID],"_MF_",param_list$MF[ID],"_all"),
                        title_ = paste0("SEM_",k,"_all"),
                        name_var = k,MF = param_list$MF[ID])
    
    Plot_SEM_with_cover(summary_sem = low_graz,pdf_ = T,
                        name = paste0("../Figures/Step1_Understanding_grazing/SEM/",param_list$Stat[ID],"/",
                                      "SEM_",param_list$Stat[ID],"_",param_list$Cov[ID],"_MF_",param_list$MF[ID],"_low"),
                        title_ = paste0("SEM_",k,"_low"),
                        name_var = k,MF = param_list$MF[ID])
    
    Plot_SEM_with_cover(summary_sem = high_graz,pdf_ = T,
                        name = paste0("../Figures/Step1_Understanding_grazing/SEM/",param_list$Stat[ID],"/",
                                      "SEM_",param_list$Stat[ID],"_",param_list$Cov[ID],"_MF_",param_list$MF[ID],"_high"),
                        title_ = paste0("SEM_",k,"_high"),
                        name_var = k,MF = param_list$MF[ID])
  }else {
    
    Plot_SEM_without_cover(summary_sem = all_d,pdf_ = T,
                           name = paste0("../Figures/Step1_Understanding_grazing/SEM/",param_list$Stat[ID],"/",
                                         "SEM_",param_list$Stat[ID],"_",param_list$Cov[ID],"_MF_",param_list$MF[ID],"_all"),
                           title_ = paste0("SEM_",k,"_all"),
                           name_var = k,MF = param_list$MF[ID])
    
    Plot_SEM_without_cover(summary_sem = low_graz,pdf_ = T,
                           name = paste0("../Figures/Step1_Understanding_grazing/SEM/",param_list$Stat[ID],"/",
                                         "SEM_",param_list$Stat[ID],"_",param_list$Cov[ID],"_MF_",param_list$MF[ID],"_low"),
                           title_ = paste0("SEM_",k,"_low"),
                           name_var = k,MF = param_list$MF[ID])
    
    Plot_SEM_without_cover(summary_sem = high_graz,pdf_ = T,
                           name = paste0("../Figures/Step1_Understanding_grazing/SEM/",param_list$Stat[ID],"/",
                                         "SEM_",param_list$Stat[ID],"_",param_list$Cov[ID],"_MF_",param_list$MF[ID],"_high"),
                           title_ = paste0("SEM_",k,"_high"),
                           name_var = k,MF = param_list$MF[ID])
  }
  
}



p=ggplot(d_indirect%>%
           filter(., With_cover=="with_cover"))+
  geom_bar(aes(x=Path,y=Indirect,fill=Type),position="dodge", stat="identity",width = .5)+
  facet_wrap(.~Stat,scales = "free",nrow = 5)+
  scale_fill_manual(values=c("#A1CA7A","#BB7CB8","#D29271"))+
  the_theme+
  geom_hline(yintercept = 0)

ggsave("../Figures/Step1_Understanding_grazing/Indirect_effects_grazing.pdf",p,width = 7,height = 9)



d_data=read.table("../Data/Spatial_structure_grazing.csv",sep=";")%>%
  Closer_to_normality(.)

d_data[,c(1,7:29,35,37:41,(44:(ncol(d_data))))] = apply(d_data[,c(1,7:29,35,37:41,(44:(ncol(d_data))))],2,z_tranform)

d_data$MF=z_tranform(rowSums(d_data[,55:(ncol(d_data))]))


#then, correlation graph
pdf("../Figures/Step1_Understanding_grazing/Correlation_functioning_metrics2.pdf",width = 12,height = 12)


list_grazing=list(0:3,c(0,1),c(2,3))
name_plots=c("All data","Low grazing pressure","High grazing pressure")

par(mfrow=c(3,2),mar=rep(4,4))
for (id in 1:length(list_grazing)){
  
  #correlation functioning metrics
  mat_cor=(rcorr(as.matrix(d_data[which(d_data$Grazing %in% list_grazing[[id]]),c(55:65)]),type = "pearson")$r)
  mat_cor[rcorr(as.matrix(d_data[,c(55:65)]),type = "pearson")$P>0.05]=0
  mat_cor[mat_cor<.15]=0
  
  
  
  network_EWS = igraph::simplify(
    graph_from_adjacency_matrix(mat_cor,
                                weighted = TRUE,
                                mode = c("undirected")))
  
  E(network_EWS)$width = E(network_EWS)$weight*3
  E(network_EWS)$lty   = sapply(E(network_EWS)$weight,function(x){ifelse(x>0,1,2)})
  
  layout_new=layout.circle(network_EWS)
  
  
  plot(network_EWS,
       vertex.color = c("#E6967F"),
       vertex.frame.width = 1,vertex.label.color="black",
       edge.curved = .3,edge.color = "grey",
       vertex.size=300*abs(colMeans(d_data[which(d_data$Grazing %in% list_grazing[[id]]),55:65],na.rm = T)),
       layout = layout_new,frame = TRUE,main=name_plots[id]
  )
  
  #correlation spatial structure 
  d_struc=d_data
  
  mat_cor=(rcorr(as.matrix(d_struc[which(d_struc$Grazing %in% list_grazing[[id]]),c(7:24)]),type = "pearson")$r)
  mat_cor[rcorr(as.matrix(d_struc[,c(7:24)]),type = "pearson")$P>0.05]=0
  mat_cor[mat_cor<.15]=0
  
  network_EWS = igraph::simplify(
    graph_from_adjacency_matrix(mat_cor,
                                weighted = TRUE,
                                mode = c("undirected")))
  
  E(network_EWS)$width = E(network_EWS)$weight*3
  
  layout_new=layout.circle(network_EWS)
  
  plot(network_EWS,
       vertex.color = c("#B9E4AC"),
       vertex.frame.width = 1,vertex.label.color="black",
       edge.curved = .3,edge.color = "grey",
       vertex.size=300*abs(colMeans(d_struc[which(d_struc$Grazing %in% list_grazing[[id]]),7:24],na.rm = T)),
       layout = layout_new,frame = TRUE,main=name_plots[id]
  )
}
dev.off()



#Same but this time the colors correspond to the slope of change along grazing gradient

saling_size_nodes=tibble(Function=c(500,500,200),Sp_stats=c(500,300,250))

formula_mod=formula(formula_(paste("value ~ Long_cos + Long_sin + Lattitude + Slope + Elevation 
      + Clim1 + Clim2 + Clim3 + Clim4
      + Woody + rho_p + Type_veg 
      + Sand + (1|Site_ID)")))



pdf("../Figures/Step1_Understanding_grazing/Correlation_functioning_metrics.pdf",width = 12,height = 12)
par(mfrow=c(3,2),mar=rep(4,4))

for (id in 1:length(list_grazing)){
  
  #correlation functioning metrics
  mat_cor=(rcorr(as.matrix(d_data[which(d_data$Grazing %in% list_grazing[[id]]),c(55:65)]),type = "pearson")$r)
  mat_cor[rcorr(as.matrix(d_data[,c(55:65)]),type = "pearson")$P>0.05]=0
  mat_cor[mat_cor<.15]=0
  
  network_EWS = igraph::simplify(
    graph_from_adjacency_matrix(mat_cor,
                                weighted = TRUE,
                                mode = c("undirected")))
  
  E(network_EWS)$width = abs(E(network_EWS)$weight)*3
  E(network_EWS)$lty   = sapply(E(network_EWS)$weight,function(x){ifelse(x>0,1,2)})
  layout_new=layout.circle(network_EWS)
  
  #We compute the slopes along Grazing gradient
  d_slope=tibble()
  d_struc=filter(d_data,Grazing %in% list_grazing[[id]])
  
  for (k in c(colnames(d_struc)[c(55:65)])){
    
    dat_lmer=d_struc%>%melt(., measure.vars=k)%>%filter(., !is.na(value))
    #first we fit the complete model and extract the residuals
    
    full_mod=lmer(formula_mod,
                  data=dat_lmer,
                  na.action = na.omit)  
    
    #residuals
    dat_lmer$resids=residuals(full_mod)
    
    mod_stat=lm(resids ~ Grazing,data=dat_lmer,na.action = na.omit)  
    d_slope=rbind(d_slope,tibble(Estimate=summary(mod_stat)$coefficients[2,1],
                                 pval=summary(mod_stat)$coefficients[2,4])%>%
                    add_column(., Stat=k)%>%
                    dplyr::rename(., Slope=Estimate))
  }
  
  plot(network_EWS,
       vertex.color = sapply(1:nrow(d_slope),function(x){
         if (d_slope$pval[x]<.05){ #significative slope
           if (d_slope$Slope[x]>0){
             return("#617CBB")
           }else {
             return("#CE5A5A")
           }
         }else {
           return("gray")
         }
       }),
       vertex.frame.width = 1,vertex.label.color="black",
       edge.curved = .3,edge.color = "grey",
       vertex.size=saling_size_nodes$Function[id]*abs(d_slope$Slope),
       layout = layout_new,frame = TRUE,main=name_plots[id]
  )
  
  #correlation spatial structure 
  d_struc=d_data
  
  mat_cor=(rcorr(as.matrix(d_struc[which(d_struc$Grazing %in% list_grazing[[id]]),c(7:24)]),type = "pearson")$r)
  mat_cor[rcorr(as.matrix(d_struc[,c(7:24)]),type = "pearson")$P>0.05]=0
  mat_cor[mat_cor<.15]=0
  
  
  network_EWS = igraph::simplify(
    graph_from_adjacency_matrix(mat_cor,
                                weighted = TRUE,
                                mode = c("undirected")))
  
  E(network_EWS)$width = abs(E(network_EWS)$weight)*3
  E(network_EWS)$lty   = sapply(E(network_EWS)$weight,function(x){ifelse(x>0,1,2)})
  
  layout_new=layout.circle(network_EWS)
  
  #We compute the slopes along Grazing gradient
  d_slope=tibble()
  d_struc=filter(d_data,Grazing %in% list_grazing[[id]])
  
  for (k in c(colnames(d_struc)[c(7:24)])){
    
    dat_lmer=d_struc%>%melt(., measure.vars=k)%>%filter(., !is.na(value))
    #first we fit the complete model and extract the residuals
    
    full_mod=lmer(formula_mod,
                  data=dat_lmer,
                  na.action = na.omit)  
    
    #residuals
    dat_lmer$resids=residuals(full_mod)
    
    mod_stat=lm(resids ~ Grazing,data=dat_lmer,na.action = na.omit)  
    d_slope=rbind(d_slope,tibble(Estimate=summary(mod_stat)$coefficients[2,1],
                                 pval=summary(mod_stat)$coefficients[2,4])%>%
                    add_column(., Stat=k)%>%
                    dplyr::rename(., Slope=Estimate))
  }
  
  plot(network_EWS,
       vertex.color = sapply(1:nrow(d_slope),function(x){
         if (d_slope$pval[x]<.05){
           if (d_slope$Slope[x]>0){
             return("#617CBB")
           }else {
             return("#CE5A5A")
           }
         }else {
           return("gray")
         }
       }),
       vertex.frame.width = 1,vertex.label.color="black",
       edge.curved = .3,edge.color = "grey",
       vertex.size=saling_size_nodes$Sp_stats[id]*abs(d_slope$Slope),
       layout = layout_new,frame = TRUE,main=name_plots[id]
  )
}

dev.off()



## >> OTHER 2/ SEM with mixed effect models ----

d_data=read.table("../Data/Spatial_structure_grazing.csv",sep=";")%>%
  Closer_to_normality(.)

d_data[,c(1,9:29,35,37:41,44:ncol(d_data))] = apply(d_data[,c(1,9:29,35,37:41,44:ncol(d_data))],2,z_tranform)

#Getting the multifunctionality index
d_data$MF=z_tranform(rowSums(d_data[,55:(ncol(d_data))],na.rm = T))

#Adding the first two components from the PCA
d_data=Perform_PCA_spatial_struc(d_data)

d_indirect=tibble() #to save indirect effects of grazing on the spatial structure
dir.create("../Figures/Step1_Understanding_grazing/SEM_lmer/",showWarnings = F)

param_list=expand.grid(Cov=c("with_cover","without_cover"),
                       Stat=c("perim_area_scaling","PL_expo","Cond_H","fmax_psd",
                              "PLR","flow_length","mean_perim_area","core_area_land",
                              "fractal_dim","division","contig",
                              "Struct1","Struct2"),
                       MF=c(T,F)) #MF instead of organic carbon

for (ID in 1:nrow(param_list)){
  
  dir.create(paste0("../Figures/Step1_Understanding_grazing/SEM_lmer/",param_list$Stat[ID]),showWarnings = F)
  
  with_cover=param_list$Cov[ID];k=as.character(param_list$Stat[ID])
  
  model_lmer=lmer(paste("value ~ (1|Site_ID) + Long_cos + Long_sin + Lattitude + Slope + Elevation + Clim1 + Clim2 + Clim3 + Clim4 + Type_veg",
                        ifelse(with_cover=="without_cover","+ rho_p","")),REML = F,
                  data = d_data%>%melt(., measure.vars=k)%>%filter(., !is.na(value)),
                  na.action = na.omit)
  
  #we remove potential outliers
  rm.outliers = romr.fnc(model_lmer, d_data%>%melt(., measure.vars=k)%>%filter(., !is.na(value)),
                         trim=2.5)
  d_data_out = rm.outliers$data
  
  #We first control for all covariates and extract the residuals
  model_lmer=lmer(paste("value ~ (1|Site_ID) + Long_cos + Long_sin + Lattitude + Slope + Elevation + Clim1 + Clim2 + Clim3 + Clim4 + Type_veg",
                        ifelse(with_cover=="without_cover","+ rho_p","")),
                  data = d_data_out,REML = F,
                  na.action = na.fail)
  
  resid_model=residuals(model_lmer) #extract residuals
  
  save=d_data[as.numeric(names(resid_model)),]%>% #add it to the dataframe 
    add_column(., Resid_mod=resid_model)
  
  if (param_list$MF[ID]){ save$Org_C=save$MF}
  
  #DOING the SEMs
  
  if (with_cover=="with_cover"){
    
    #SEM with all grazing intensity
    
    d_sem=save
    all_d=summary(psem(
      lmer(Resid_mod ~ (1|Site_ID) + Grazing + Woody + rho_p , d_sem),
      lmer(Woody ~ (1|Site_ID) + Grazing + Org_C + Sand, d_sem),
      lmer(Sand ~ (1|Site_ID) + Grazing , d_sem),
      lmer(rho_p ~ (1|Site_ID) + Grazing + Org_C + Sand , d_sem),
      lmer(Org_C ~ (1|Site_ID) + Grazing, d_sem)
    ))
    
    #SEM with low grazing intensity
    
    d_sem=filter(save,Grazing %in% 0:1)
    low_graz=summary(psem(
      lmer(Resid_mod ~ (1|Site_ID) + Grazing + Woody + rho_p , d_sem),
      lmer(Woody ~ (1|Site_ID) + Grazing + Org_C + Sand, d_sem),
      lmer(Sand ~ (1|Site_ID) + Grazing , d_sem),
      lmer(rho_p ~ (1|Site_ID) + Grazing + Org_C + Sand , d_sem),
      lmer(Org_C ~ (1|Site_ID) + Grazing, d_sem)
    ))
    
    #SEM with high grazing intensity
    
    d_sem=filter(save,Grazing %in% 2:3)
    high_graz=summary(psem(
      lmer(Resid_mod ~ (1|Site_ID) + Grazing + Woody + rho_p , d_sem),
      lmer(Woody ~ (1|Site_ID) + Grazing + Org_C + Sand, d_sem),
      lmer(Sand ~ (1|Site_ID) + Grazing , d_sem),
      lmer(rho_p ~ (1|Site_ID) + Grazing + Org_C + Sand , d_sem),
      lmer(Org_C ~ (1|Site_ID) + Grazing, d_sem)
    ))
    
    
  }else {
    
    #SEM with all grazing intensity
    
    d_sem=save
    all_d=summary(psem(
      lmer(Resid_mod ~ (1|Site_ID) + Grazing + Woody , d_sem),
      lmer(Woody ~ (1|Site_ID) + Grazing + Org_C + Sand, d_sem),
      lmer(Sand ~ (1|Site_ID) + Grazing , d_sem),
      lmer(Org_C ~ (1|Site_ID) + Grazing, d_sem)
    ))
    
    #SEM with low grazing intensity
    
    d_sem=filter(save,Grazing %in% 0:1)
    low_graz=summary(psem(
      lmer(Resid_mod ~ (1|Site_ID) + Grazing + Woody , d_sem),
      lmer(Woody ~ (1|Site_ID) + Grazing + Org_C + Sand, d_sem),
      lmer(Sand ~ (1|Site_ID) + Grazing , d_sem),
      lmer(Org_C ~ (1|Site_ID) + Grazing, d_sem)
    ))
    
    #SEM with high grazing intensity
    
    d_sem=filter(save,Grazing %in% 2:3)
    high_graz=summary(psem(
      lmer(Resid_mod ~ (1|Site_ID) + Grazing + Woody , d_sem),
      lmer(Woody ~ (1|Site_ID) + Grazing + Org_C + Sand, d_sem),
      lmer(Sand ~ (1|Site_ID) + Grazing , d_sem),
      lmer(Org_C ~ (1|Site_ID) + Grazing, d_sem)
    ))
  }
  
  
  d_indirect=rbind(d_indirect,
                   Get_indirect_effects_grazing(all_d)%>%
                     add_column(., Type="All",Stat=k,With_cover=with_cover,MF=param_list$MF[ID]),
                   Get_indirect_effects_grazing(high_graz)%>%
                     add_column(., Type="High",Stat=k,With_cover=with_cover,MF=param_list$MF[ID]),
                   Get_indirect_effects_grazing(low_graz)%>%
                     add_column(., Type="Low",Stat=k,With_cover=with_cover,MF=param_list$MF[ID]))
  
  #ploting
  
  if(with_cover=="with_cover"){
    
    Plot_SEM_with_cover(summary_sem = all_d,pdf_ = T,
                        name = paste0("../Figures/Step1_Understanding_grazing/SEM_lmer/",param_list$Stat[ID],"/",
                                      "SEM_",param_list$Stat[ID],"_",param_list$Cov[ID],"_MF_",param_list$MF[ID],"_all"),
                        title_ = paste0("SEM_",k,"_all"),
                        name_var = k,MF = param_list$MF[ID])
    
    Plot_SEM_with_cover(summary_sem = low_graz,pdf_ = T,
                        name = paste0("../Figures/Step1_Understanding_grazing/SEM_lmer/",param_list$Stat[ID],"/",
                                      "SEM_",param_list$Stat[ID],"_",param_list$Cov[ID],"_MF_",param_list$MF[ID],"_low"),
                        title_ = paste0("SEM_",k,"_low"),
                        name_var = k,MF = param_list$MF[ID])
    
    Plot_SEM_with_cover(summary_sem = high_graz,pdf_ = T,
                        name = paste0("../Figures/Step1_Understanding_grazing/SEM_lmer/",param_list$Stat[ID],"/",
                                      "SEM_",param_list$Stat[ID],"_",param_list$Cov[ID],"_MF_",param_list$MF[ID],"_high"),
                        title_ = paste0("SEM_",k,"_high"),
                        name_var = k,MF = param_list$MF[ID])
  }else {
    
    Plot_SEM_without_cover(summary_sem = all_d,pdf_ = T,
                           name = paste0("../Figures/Step1_Understanding_grazing/SEM_lmer/",param_list$Stat[ID],"/",
                                         "SEM_",param_list$Stat[ID],"_",param_list$Cov[ID],"_MF_",param_list$MF[ID],"_all"),
                           title_ = paste0("SEM_",k,"_all"),
                           name_var = k,MF = param_list$MF[ID])
    
    Plot_SEM_without_cover(summary_sem = low_graz,pdf_ = T,
                           name = paste0("../Figures/Step1_Understanding_grazing/SEM_lmer/",param_list$Stat[ID],"/",
                                         "SEM_",param_list$Stat[ID],"_",param_list$Cov[ID],"_MF_",param_list$MF[ID],"_low"),
                           title_ = paste0("SEM_",k,"_low"),
                           name_var = k,MF = param_list$MF[ID])
    
    Plot_SEM_without_cover(summary_sem = high_graz,pdf_ = T,
                           name = paste0("../Figures/Step1_Understanding_grazing/SEM_lmer/",param_list$Stat[ID],"/",
                                         "SEM_",param_list$Stat[ID],"_",param_list$Cov[ID],"_MF_",param_list$MF[ID],"_high"),
                           title_ = paste0("SEM_",k,"_high"),
                           name_var = k,MF = param_list$MF[ID])
  }
  
}


## >> OTHER 3/ SEM replace cover ----

Plot_SEM_replace_cover=function(summary_sem,pdf_=F,name="SEM",title_="",name_var="",MF=F){
  
  l=summary_sem$coefficients[,c("Predictor","Response","Estimate","P.Value")]
  l$color="#9ED471"
  l$color[which(l$Estimate<0)]="#EFC46A"
  l$color[which(l$P.Value>0.1)]="gray"
  
  g = graph.data.frame(l, directed=T)
  g= g %>% set_edge_attr("color", value =l$color)
  g= g %>% set_vertex_attr("name", value =c("Grazing","Woody","Facilitation","Sand",
                                            name_var))
  coord=data.frame(label=vertex_attr(g, "name"),
                   lab2=c("Grazing","Woody","Facilitation","Sand",
                          name_var),
                   x=c(-10,10,0,0,10),y=c(0,-15,30,-30,15))
  EL=as_edgelist(g)
  EL=cbind(EL,l[,3])
  
  name_edge=sapply(1:length(summary_sem$coefficients$P.Value),function(x){#adding pvalue
    return(paste0(round(l[x,3],3),is_signif(summary_sem$coefficients$P.Value[x])))
  })
  
  asi=abs(l[,3])/0.05
  asi[asi<5]=5
  
  if(pdf_){
    pdf(paste0("./",name,".pdf"),width = 7,height = 4)
    qgraph(EL,layout=as.matrix(coord[,c("x","y")]),edge.color=l$color,edge.labels=name_edge,
           border.color="white",label.cex=1.5,label.scale=F,title=title_,
           edge.label.cex = 1.5,edge.label.position=0.25,vsize2=4,vsize=30,
           shape="ellipse",edge.labels=T,fade=F,esize=5,asize=asi,
           mar=rep(3,4))
    dev.off()
    
  }else {
    qgraph(EL,layout=as.matrix(coord[,c("x","y")]),edge.color=l$color,edge.labels=name_edge,
           border.color="white",label.cex=1.5,label.scale=F,
           edge.label.cex = 1.5,edge.label.position=0.25,vsize2=4,vsize=30,
           shape="ellipse",edge.labels=T,fade=F,esize=5,asize=asi,
           mar=rep(3,4))
    
  }
}


d_data=read.table("../Data/Spatial_structure_grazing.csv",sep=";")%>%
  Closer_to_normality(.)

d_data[,c(1,9:29,35,37:41,44:ncol(d_data))] = apply(d_data[,c(1,9:29,35,37:41,44:ncol(d_data))],2,z_tranform)

#Getting the multifunctionality index
d_data$MF=z_tranform(rowSums(d_data[,55:(ncol(d_data))],na.rm = T))

#Adding the first two components from the PCA
d_data=Perform_PCA_spatial_struc(d_data)

d_indirect=tibble() #to save indirect effects of grazing on the spatial structure
dir.create("../Figures/Step1_Understanding_grazing/SEMtest",showWarnings = F)

param_list=expand.grid(Cov=c("with_cover"),
                       Stat=c("perim_area_scaling","flow_length",
                              "fractal_dim","core_area",
                              "Struct1","Struct2"),
                       MF=c(F)) #MF instead of organic carbon

for (ID in 1:nrow(param_list)){
  
  dir.create(paste0("../Figures/Step1_Understanding_grazing/SEMtest/",param_list$Stat[ID]),showWarnings = F)
  
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
  
  #SEM with all grazing intensity
  
  d_sem=save
  all_d=summary(psem(
    lm(Resid_mod ~ Grazing + Woody + Org_C + Sand , d_sem),
    lm(Woody ~ Grazing + Org_C + Sand, d_sem),
    lm(Sand ~ Grazing , d_sem),
    lm(Org_C ~ Grazing, d_sem)
  ))
  
  #SEM with low grazing intensity
  
  d_sem=filter(save,Grazing %in% 0:1)
  low_graz=summary(psem(
    lm(Resid_mod ~ Grazing + Woody + Org_C + Sand , d_sem),
    lm(Woody ~ Grazing + Org_C + Sand, d_sem),
    lm(Sand ~ Grazing , d_sem),
    lm(Org_C ~ Grazing, d_sem)
  ))
  
  #SEM with high grazing intensity
  
  d_sem=filter(save,Grazing %in% 2:3)
  high_graz=summary(psem(
    lm(Resid_mod ~ Grazing + Woody + Org_C + Sand, d_sem),
    lm(Woody ~ Grazing + Org_C + Sand, d_sem),
    lm(Sand ~ Grazing , d_sem),
    lm(Org_C ~ Grazing, d_sem)
  ))
  
  
  #ploting
  
  
  Plot_SEM_replace_cover(summary_sem = all_d,pdf_ = T,
                         name = paste0("../Figures/Step1_Understanding_grazing/SEMtest/",param_list$Stat[ID],"/",
                                       "SEM_",param_list$Stat[ID],"_",param_list$Cov[ID],"_MF_",param_list$MF[ID],"_all"),
                         title_ = paste0("SEM_",k,"_all"),
                         name_var = k,MF = param_list$MF[ID])
  
  Plot_SEM_replace_cover(summary_sem = low_graz,pdf_ = T,
                         name = paste0("../Figures/Step1_Understanding_grazing/SEMtest/",param_list$Stat[ID],"/",
                                       "SEM_",param_list$Stat[ID],"_",param_list$Cov[ID],"_MF_",param_list$MF[ID],"_low"),
                         title_ = paste0("SEM_",k,"_low"),
                         name_var = k,MF = param_list$MF[ID])
  
  Plot_SEM_replace_cover(summary_sem = high_graz,pdf_ = T,
                         name = paste0("../Figures/Step1_Understanding_grazing/SEMtest/",param_list$Stat[ID],"/",
                                       "SEM_",param_list$Stat[ID],"_",param_list$Cov[ID],"_MF_",param_list$MF[ID],"_high"),
                         title_ = paste0("SEM_",k,"_high"),
                         name_var = k,MF = param_list$MF[ID])
}
