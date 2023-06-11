rm(list=ls())
source("./Structure_grazing_function.R")

# ---------------- Step 0: Preliminary analyses ----

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

mat_cor=abs(cor(d_data[,c(4:24)],use = "na.or.complete",method = "pearson"))
mat_cor[mat_cor<.1]=0
library(igraph)
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



# ---------------- Step 1: Going big ----

## >> Importance ----

d_all=read.table("../Data/Step1_Going_big/Importance.csv",sep=";")

p=ggplot(rbind(cbind("Sp_stat"="Averaged across \n spatial statistics",as_tibble(t(colMeans(d_all[,-1])))),d_all)%>%
           add_column(., "R2"=0)%>%
           melt(., measure.vars=colnames(.)[5:ncol(.)])%>%
           mutate(., Sp_stat=recode_factor(Sp_stat,"fmax"="Max patch","perim_area_scaling"="Perim-area scaling",
                                            "PL"="PSD exponent","spectral_ratio"="Spectral ratio",
                                           "Cond_H"="Conditional entropy"))%>%
           mutate(., variable=recode_factor(variable,"Org_C"="Facilitation",
                                            "Interactions"="Interactions \n with grazing"))%>%
           add_column(., Order_stat=rep(rev(1:5),12))%>%
           mutate(Sp_stat = fct_reorder(Sp_stat, Order_stat))%>%
           add_column(., Order_f=rep(c(10,7,9,11,4,2,3,5,8,1,6,12),each=5))%>%
           mutate(variable = fct_reorder(variable, Order_f)))+
  geom_tile(aes(x=variable,y=Sp_stat,fill=value))+
  geom_text(data=tibble(x=12,y=1:4,label=round(d_all$R2m,2)),aes(x=x,y=y,label=paste0(label)))+
  geom_text(data=tibble(x=1:11,y=5,
                        label=round(sort(as.numeric(as_tibble(t(colMeans(d_all[,-1]))))[-c(1:3)],decreasing = T),2)),
            aes(x=x,y=y,label=paste0(label)))+
  scale_fill_gradient2(low = "white",mid="#F39565",high = "#C12121",midpoint = .5)+
  the_theme+
  theme(axis.text.x = element_text(angle=c(rep(60,11),0),hjust=c(rep(1,11),0)))+
  labs(x="Predictors",y="Spatial statistic",fill="Importance (%)")

ggsave("../Figures/Step1_Going_big/Importance_grazing.pdf",p,width = 7,height = 5)


## >> Direct & interactions: LME models ----

d_all2=read.table("../Data/Step1_Going_big/Estimators_model.csv",sep=";")


for (k in unique(d_all2$Stat)){
  
  p1_1=ggplot(d_all2%>%
                filter(., Stat==k)%>%
                filter(., term !="(Intercept)")%>%
                Organize_df(., "predictor"))+
    geom_pointrange(aes(x=observed,y=term,xmin=observed-se,xmax=observed+se,color=Type_pred))+
    the_theme+
    labs(x="",y="",color="")+
    geom_vline(xintercept = 0,linetype=9)+
    scale_color_manual(values=c("Geo"="#FFB15B","Grazing"="#F3412A",
                                "Abiotic"="#4564CE","Vegetation"="#6DD275","Climatic"="#CC80C7"))+
    theme(legend.position = "none")
  
  p1_2=ggplot(d_all2%>%
                filter(., Stat==k)%>%
                add_column(., xplot=0)%>%
                filter(., term !="(Intercept)")%>%
                Organize_df(., "bar")%>%
                group_by(., Type_pred,Stat,xplot)%>%
                dplyr::summarise(., sum_effect=sum(abs(observed)),.groups = "keep")%>%
                mutate(., sum_effect=sum_effect/sum(sum_effect)))+
    geom_bar(aes(x=xplot,y=sum_effect,fill=Type_pred),stat="identity",width = .05)+
    geom_text(aes(x=0,y=1.05,label=paste0("R² \n = \n ",round(d_all$R2m[which(d_all$Sp_stat==k)],2))),size=3)+
    the_theme+
    labs(y="Relative effects of estimates",fill="")+
    scale_fill_manual(values=c("Geo"="#FFB15B","Grazing"="#F3412A",
                               "Abiotic"="#4564CE","Vegetation"="#6DD275","Climatic"="#CC80C7"))+
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),axis.ticks.x = element_blank(),
          axis.line.x = element_blank())
  
  p1=ggarrange(
    ggarrange(p1_1,p1_2+theme(legend.position = "none"),widths = c(1,.2),labels = letters[1:2]),
    get_legend(p1_2),nrow=2,heights = c(1,.15)
  )
  
  
  ggsave(paste0("../Figures/Step1_Going_big/Preditor_",k,".pdf"),p1,width = 6,height = 6)
  
  
}



## >> Indirect effects: SEM ----

d_data=read.table("../Data/Spatial_structure_grazing.csv",sep=";")%>%
  add_column(., Long_sin=sin(.$Longitude),Long_cos=cos(.$Longitude))

d_data[,c(1,9:29,35,37:41,44:ncol(d_data))] = apply(d_data[,c(1,9:29,35,37:41,44:ncol(d_data))],2,z_tranform)


#We first control for all covariates and extract the residuals
for (with_cover in c("with_cover","without_cover")){
  
  dir.create(paste0("../Figures/Step1_Going_big/SEM_",with_cover),showWarnings = F)
  
  for (k in c("perim_area_scaling","PL_expo","Cond_H","fmax_psd","PLR","flow_length","mean_perim_area","core_area_land","fractal_dim","division",
              "contig")){
  
    
    
    model_lmer=lmer(paste("value ~ Long_cos + Long_sin + Lattitude + Slope + Elevation + Clim1 + Clim2 + Clim3 + Clim4 + Type_veg + (1|Site_ID)",
                  ifelse(with_cover=="without_cover","+ rho_p","")),
                    data = d_data%>%melt(., measure.vars=k)%>%filter(., !is.na(value)),
                    na.action = na.fail,REML ="FALSE")
    
    #potential outliers
    mcp.fnc(model_lmer)
    
    #we remove it
    rm.outliers = romr.fnc(model_lmer, d_data%>%melt(., measure.vars=k)%>%filter(., !is.na(value)),
                           trim=2.5)
    d_data_out = rm.outliers$data
    
    model_lmer=lmer(paste("value ~ Long_cos + Long_sin + Lattitude + Slope + Elevation + Clim1 + Clim2 + Clim3 + Clim4 + Type_veg + (1|Site_ID)",
                          ifelse(with_cover=="without_cover","+ rho_p","")),
                    data = d_data_out,
                    na.action = na.fail,REML ="FALSE")
    
    resid_model=residuals(model_lmer) #extract residuals
    
    save=d_data[as.numeric(names(resid_model)),]%>% #add it to the dataframe 
      add_column(., Resid_mod=resid_model)
    
    
    #DOING the SEMs
    
    if (with_cover=="with_cover"){
      
      #SEM with all grazing intensity
      
      d_sem=save
      all_d=summary(psem(
        lm(Resid_mod ~ Grazing + Woody + Sand + Org_C + rho_p, d_sem),
        lm(Woody ~ Grazing, d_sem),
        lm(Sand ~ Grazing , d_sem),
        lm(rho_p ~ Grazing , d_sem),
        lm(Org_C ~ Grazing, d_sem)
      ))
      
      #SEM with low grazing intensity
      
      d_sem=filter(save,Grazing %in% 0:1)
      low_graz=summary(psem(
        lm(Resid_mod ~ Grazing + Woody + Sand + Org_C + rho_p, d_sem),
        lm(Woody ~ Grazing, d_sem),
        lm(Sand ~ Grazing , d_sem),
        lm(rho_p ~ Grazing , d_sem),
        lm(Org_C ~ Grazing, d_sem)
      ))
      
      #SEM with high grazing intensity
      
      d_sem=filter(save,Grazing %in% 2:3)
      high_graz=summary(psem(
        lm(Resid_mod ~ Grazing + Woody + Sand + Org_C + rho_p, d_sem),
        lm(Woody ~ Grazing, d_sem),
        lm(Sand ~ Grazing , d_sem),
        lm(rho_p ~ Grazing , d_sem),
        lm(Org_C ~ Grazing, d_sem)
      ))
      
    }else {
      
      #SEM with all grazing intensity
      
      d_sem=save
      all_d=summary(psem(
        lm(Resid_mod ~ Grazing + Woody + Sand + Org_C , d_sem),
        lm(Woody ~ Grazing, d_sem),
        lm(Sand ~ Grazing , d_sem),
        lm(Org_C ~ Grazing, d_sem)
      ))
      
      #SEM with low grazing intensity
      
      d_sem=filter(save,Grazing %in% 0:1)
      low_graz=summary(psem(
        lm(Resid_mod ~ Grazing + Woody + Sand + Org_C, d_sem),
        lm(Woody ~ Grazing, d_sem),
        lm(Sand ~ Grazing , d_sem),
        lm(Org_C ~ Grazing, d_sem)
      ))
      
      #SEM with high grazing intensity
      
      d_sem=filter(save,Grazing %in% 2:3)
      high_graz=summary(psem(
        lm(Resid_mod ~ Grazing + Woody + Sand + Org_C, d_sem),
        lm(Woody ~ Grazing, d_sem),
        lm(Sand ~ Grazing , d_sem),
        lm(Org_C ~ Grazing, d_sem)
      ))
      
    }
    

    
    
    #ploting
    
    if(with_cover=="with_cover"){
      Plot_SEM_with_cover(all_d,T,name = paste0("../Figures/Step1_Going_big/",
                                                paste0("SEM_",with_cover,"/"),"SEM_",k,"_all"),
                          title_ = paste0("SEM_",k,"_all"),
                          name_var = k)
      Plot_SEM_with_cover(low_graz,T,name = paste0("../Figures/Step1_Going_big/",
                                                   paste0("SEM_",with_cover,"/"),"SEM_",k,"_low"),
                          title_ = paste0("SEM_",k,"_low"),
                          name_var = k)
      Plot_SEM_with_cover(high_graz,T,name = paste0("../Figures/Step1_Going_big/",
                                                    paste0("SEM_",with_cover,"/"),"SEM_",k,"_high"),
                          title_ = paste0("SEM_",k,"_high"),
                          name_var = k)
      
    }else {
      Plot_SEM_without_cover(all_d,T,name = paste0("../Figures/Step1_Going_big/",
                                                paste0("SEM_",with_cover,"/"),"SEM_",k,"_all"),
                          title_ = paste0("SEM_",k,"_all"),
                          name_var = k)
      Plot_SEM_without_cover(low_graz,T,name = paste0("../Figures/Step1_Going_big/",
                                                   paste0("SEM_",with_cover,"/"),"SEM_",k,"_low"),
                          title_ = paste0("SEM_",k,"_low"),
                          name_var = k)
      Plot_SEM_without_cover(high_graz,T,name = paste0("../Figures/Step1_Going_big/",
                                                    paste0("SEM_",with_cover,"/"),"SEM_",k,"_high"),
                          title_ = paste0("SEM_",k,"_high"),
                          name_var = k)
      
    }
    
  }

}







