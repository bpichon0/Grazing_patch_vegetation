## >> RF analysis  ----

#We want to understand the drivers of spatial structure
d_all=d_all2=tibble()#for keeping all information

d_data=read.table("../Data/Spatial_structure_grazing.csv",sep=";")%>%
  Closer_to_normality(.)

d_data[,c(1,9:29,35,37:41,44:ncol(d_data))] = apply(d_data[,c(1,9:29,35,37:41,44:ncol(d_data))],2,z_tranform)
pdf("./RF.pdf",p1,width = 13,height = 7)


for (stat in c("fractal_dim","fmax_psd","Spectral_ratio","PL_expo",
               "flow_length","perim_area_scaling","mean_perim_area",
               "core_area","core_area_land","contig")){
  
  subset_data=d_data[,c(which(colnames(d_data)==stat),34:ncol(d_data))]%>%
    melt(., measure.vars = stat)%>%
    filter(., !is.na(value))
  
  labels=subset_data[,"value"]
  preds=data.matrix(subset_data[,-which(colnames(subset_data)=="value")])
  
  set.seed(200)
  trainID = sample(1:nrow(subset_data),size = floor(nrow(subset_data)*0.7))
  validid = (1:nrow(subset_data))
  validid = validid[!(validid %in% trainID)]
  
  dtrain = xgb.DMatrix(data = preds[trainID,], label= labels[trainID])
  dtest = xgb.DMatrix(data = preds[validid,], label= labels[validid])
  dall = xgb.DMatrix(data = preds, label= labels)
  
  params_booster=list(booster = 'gbtree', 
                      eta = 0.2, 
                      gamma = 0, 
                      max.depth = 6, 
                      subsample = 1, 
                      colsample_bytree = 1, 
                      eval_metric="error",
                      colsample_bytree=1)
  
  bst.cv = xgb.cv(data = preds[trainID,],
                  label = as.numeric(labels)[trainID], 
                  params = params_booster,
                  nrounds = 200, 
                  nfold = 5,
                  print_every_n = 20,
                  verbose = 2)
  
  nrounds_best=which.min(bst.cv$evaluation_log$test_error_mean)
  
  set.seed(300)
  model_tuned = xgboost(data = preds[trainID,],
                        label = as.numeric(labels)[trainID],
                        booster = "gbtree",
                        eval_metric="error",
                        max.depth = 6, 
                        nround = nrounds_best, 
                        eta=.2,
                        gamma=0,
                        subsample=1,
                        nfold=5,
                        colsample_bytree=1)
  
  importance_matrix = xgb.importance(colnames(preds), model = model_tuned)
  
  
  
  #The SHAP
  
  library(fastshap)
  pfun = function(object, newdata) {
    predict(object, data = newdata)$predictions
  }
  predsi=preds
  EXPVAL=explain(model_tuned,X=predsi,pred_wrapper = pfun,exact = TRUE)
  predsi=as.data.frame(predsi)
  
  p1=importance_matrix%>%
    mutate(name = fct_reorder(Feature, Gain)) %>%
    ggplot( aes(x=name, y=Gain)) +
    geom_bar(stat="identity", fill="black", width=.4) +
    coord_flip() +
    labs(x="",y="Importance") +
    ggtitle(stat)+
    the_theme+
    theme(panel.grid.major = element_line(color = "gray90"))
  
  print(p1)
  
  shap_long_NEG=tibble()
  for(i in 1:ncol(EXPVAL)){
    var=colnames(predsi)[i]
    shap_long_NEG=rbind(shap_long_NEG,tibble(variable=var,value=EXPVAL[,i],rfvalue=predsi[,i]))
  }
  
  
  
  shap_long_NEG$Importance=importance_matrix$Gain[match(shap_long_NEG$variable,importance_matrix$Feature)]
  shap_long_NEG=shap_long_NEG %>%
    mutate(name=fct_reorder(variable, Importance))
  
  p1=ggplot(shap_long_NEG,aes(x=rfvalue,y=value))+
    geom_point(size=0.2,alpha=0.5,color="gray")+
    geom_smooth(color="palegreen")+
    xlab("\nFeature value")+
    ylab("SHAP value \n[impact on model output]\n")+
    ggtitle(stat)+
    facet_wrap(~variable,scales = "free",ncol = 7)+
    the_theme
  
  print(p1)
}

dev.off()










## >> Replicate Wang paper ----

d_data=read.table("../Data/Spatial_structure_grazing.csv",sep=";")%>%
  add_column(., Long_sin=sin(.$Longitude),Long_cos=cos(.$Longitude))

d_data=d_data%>%  
  add_column(.,PPQW=as.numeric(sapply(1:nrow(.),
                                      function(x){return(d_biodesert$RAWAQ[which(d_biodesert$ID==.$Site_ID[x])])})))

d_data[,c(1,9:29,35,37:41,44:ncol(d_data))] = apply(d_data[,c(1,9:29,35,37:41,44:ncol(d_data))],2,z_tranform)

d_all2=tibble()


stat="fmax_psd"
print(stat)

d_data_mod=d_data%>%melt(., measure.vars=stat)%>%
  filter(., !is.na(value),!is.na(PPQW))%>%
  filter(., Grazing !=0)

formula_mod=formula(formula_(paste("value ~ Long_cos + Long_sin + Lattitude + Slope + Elevation 
      + Aridity*rho_p + Aridity*Org_C
      + Aridity*Grazing
      + PPQW*Grazing + Org_C*Grazing
      + Sand*Grazing + rho_p*Grazing
      + (1|Site_ID) + (1|Sub_id)")))

model_spa_stat  = lmer(formula_mod, d_data_mod,
                       na.action = "na.fail" ,REML ="FALSE")

#we remove potential outliers
rm.outliers = romr.fnc(model_spa_stat, d_data_mod, trim=2.5)
d_data_out = rm.outliers$data

model_spa_stat = lmer(formula_mod, data = d_data_out, 
                      na.action = na.fail,REML ="FALSE")


#Get R? of the full model
model_spa_stat=lmer(formula_mod, data = d_data_out,
                    na.action = na.fail,REML ="FALSE")




test=visreg::visreg(fit = model_spa_stat,xvar = "Aridity",by = "Grazing",plot=T)
plot(ggpredict(model_spa_stat, c("Org_C", "Grazing")), residuals = TRUE, grid = TRUE)

ggplot(test$res)+
  geom_point(aes(x=Aridity,y=visregRes, color=as.factor(Grazing)))+the_theme+
  scale_color_manual(values=c("green","yellow","red"))+
  geom_smooth(aes(x=Aridity,y=visregRes, color=as.factor(Grazing)),method = "lm")

## >> Correlation multi-functionality with low and high grazing pressures ----

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





## >> Importance and effect of grazing standardize effects ----

#Importance based on R2
d_R2=read.table(paste0("../Data/Step1_Understanding_grazing/Importance.csv"),sep=";")%>%
  filter(., With_cover==1)%>%
  dplyr::select(., Sp_stat,R2m,R2C)%>%
  dplyr::rename(., Stat=Sp_stat)%>%
  Rename_spatial_statistics(.)

d_all2=read.table(paste0("../Data/Step1_Understanding_grazing/Estimators_model.csv"),sep=";")%>%
  dplyr::rename(., observed=Median)%>%
  filter(., With_cover==1)%>%
  Organize_df(., "bar")%>%
  group_by(., Type_pred,Stat)%>%
  dplyr::summarise(., sum_effect=sum(abs(observed)),.groups = "keep")%>%
  Rename_spatial_statistics(.)%>%group_by(., Stat)%>%
  dplyr::summarise(., .groups="keep",sum_effect=sum_effect/sum(sum_effect),Type_pred=Type_pred)


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


d_all2=read.table("../Data/Step1_Understanding_grazing/Estimators_model.csv",sep=";")

p2=ggplot(d_all2%>%
            filter(., term=="Grazing",
                   With_cover==1)%>%
            Rename_spatial_statistics(.)%>%
            add_column(., Signif=sapply(1:nrow(.),function(x){
              return(is_signif(.$pvalue[x]))}))%>%
            arrange(., Median)%>%
            add_column(., Order_f=1:nrow(.))%>%
            mutate(.,Stat = fct_reorder(Stat, Order_f))%>%
            add_column(., order_signif=sapply(1:nrow(.),
                                              function(x){return(ifelse(.$Median[x]<0,.4,-.25))})))+
  geom_pointrange(aes(x=Median,y=Stat,xmin=q1,xmax=q3,color=pvalue<.05),shape=17,size=.8)+
  the_theme+
  scale_color_manual(values=c("grey","black"))+
  geom_text(aes(x=order_signif,y=1:14,label=Signif))+
  geom_vline(xintercept = 0,linetype=9)+
  geom_hline(yintercept = c(3.5,6.5,9.5),linetype=9)+
  labs(x=substitute(paste(beta," (Grazing)")),y="",color="")+
  theme(legend.position = "none")

p_tot=ggarrange(
  ggarrange(p1+theme(legend.position = "none"),
            p2,
            ncol=2,labels = letters[1:2],widths = c(1,.75)),
  ggarrange(ggplot()+theme_void(),get_legend(p1),ggplot()+theme_void(),ncol=3,widths = c(.3,1,.3)),
  nrow=2,heights = c(1,.1))


ggsave("../Figures/Final_figs/Importance_effect_grazing.pdf",
       p_tot,width = 10,height = 6)




#Importance based on R2
d_R2=read.table(paste0("../Data/Step1_Understanding_grazing/Importance.csv"),sep=";")%>%
  filter(., With_cover==0)%>%
  dplyr::select(., Sp_stat,R2m,R2C)%>%
  dplyr::rename(., Stat=Sp_stat)%>%
  Rename_spatial_statistics(.)

d_all2=read.table(paste0("../Data/Step1_Understanding_grazing/Estimators_model.csv"),sep=";")%>%
  dplyr::rename(., observed=Median)%>%
  filter(., With_cover==0)%>%
  Organize_df(., "bar")%>%
  group_by(., Type_pred,Stat)%>%
  dplyr::summarise(., sum_effect=sum(abs(observed)),.groups = "keep")%>%
  Rename_spatial_statistics(.)%>%group_by(., Stat)%>%
  dplyr::summarise(., .groups="keep",sum_effect=sum_effect/sum(sum_effect),Type_pred=Type_pred)


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


d_all2=read.table("../Data/Step1_Understanding_grazing/SI/Estimators_model.csv",sep=";")

p2=ggplot(d_all2%>%
            filter(., term=="Grazing",
                   With_cover==0)%>%
            Rename_spatial_statistics(.)%>%
            add_column(., Signif=sapply(1:nrow(.),function(x){
              return(is_signif(.$pvalue[x]))}))%>%
            arrange(., Median)%>%
            add_column(., Order_f=1:nrow(.))%>%
            mutate(.,Stat = fct_reorder(Stat, Order_f))%>%
            add_column(., order_signif=sapply(1:nrow(.),
                                              function(x){return(ifelse(.$Median[x]<0,.4,-.25))})))+
  geom_pointrange(aes(x=Median,y=Stat,xmin=q1,xmax=q3,color=pvalue<.05),shape=17,size=.8)+
  the_theme+
  scale_color_manual(values=c("grey","black"))+
  geom_text(aes(x=order_signif,y=1:14,label=Signif))+
  geom_vline(xintercept = 0,linetype=9)+
  geom_hline(yintercept = c(3.5,6.5,9.5),linetype=9)+
  labs(x=substitute(paste(beta," (Grazing)")),y="",color="")+
  theme(legend.position = "none")

p_tot=ggarrange(
  ggarrange(p1+theme(legend.position = "none"),
            p2,
            ncol=2,labels = letters[1:2],widths = c(1,.75)),
  ggarrange(ggplot()+theme_void(),get_legend(p1),ggplot()+theme_void(),ncol=3,widths = c(.3,1,.3)),
  nrow=2,heights = c(1,.1))


ggsave("../Figures/Final_figs/SI/Importance_effect_grazing_without_cover.pdf",
       p_tot,width = 10,height = 6)
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
  geom_text(data=d_R2,aes(y=Stat,x=8,label=paste0("r? = ",round(R2m,2))))+
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
        geom_text(aes(x=0,y=1.05,label=paste0("R? \n = \n ",round(d_all$R2m[which(d_all$Sp_stat==k)],2))),size=3)+
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
  geom_pointrange(aes(x=Org_C,y=p2,ymax=p3,ymin=p1),color="gray",alpha=.5)+
  geom_smooth(aes(x=Org_C,y=p2),method = "lm",color="black")+
  the_theme+
  labs(y="Median of p posterior")

p2_2=ggplot(d_data)+
  geom_pointrange(aes(x=Org_C,y=q2,ymax=q3,ymin=q1),color="gray",alpha=.5)+
  geom_smooth(aes(x=Org_C,y=q2),method = "lm",color="black")+
  the_theme+
  labs(y="Median of q posterior")


p3_1=ggplot(d_data)+
  geom_pointrange(aes(x=Aridity,y=p2,ymax=p3,ymin=p1),color="gray",alpha=.5)+
  geom_smooth(aes(x=Aridity,y=p2),method = "lm",color="black")+
  the_theme+
  labs(y="Median of p posterior")

p3_2=ggplot(d_data)+
  geom_pointrange(aes(x=Aridity,y=q2,ymax=q3,ymin=q1),color="gray",alpha=.5)+
  geom_smooth(aes(x=Aridity,y=q2),method = "lm",color="black")+
  the_theme+
  labs(y="Median of q posterior")

p_tot=ggarrange(p1_1,p1_2,p2_1,p2_2,p3_1,p3_2,nrow=3,ncol=2,labels = c("",letters[1],"",letters[2],"",letters[3]))

ggsave("../Figures/Step0_preliminary/Inferrence/Params_grazing_aridity.pdf",p_tot,width = 7,height = 10)


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
                              "fractal_dim","core_area","core_area_land",
                              "PL_expo",
                              "Struct1","Struct2"),
                       MF=c(T)) #MF instead of organic carbon

for (ID in 1:nrow(param_list)){
  
  dir.create(paste0("../Figures/Step1_Understanding_grazing/SEMtest/",param_list$Stat[ID]),showWarnings = F)
  
  with_cover=param_list$Cov[ID];k=as.character(param_list$Stat[ID])
  
  model_lmer=lmer(paste("value ~ (1|Site_ID) + Long_cos + Long_sin + Lattitude + Slope + Elevation + Clim1 + Clim2 + Clim3 + Clim4 + rho_p",
                        ifelse(with_cover=="without_cover","+ rho_p","")),
                  data = d_data%>%melt(., measure.vars=k)%>%filter(., !is.na(value)),
                  na.action = na.omit,REML ="FALSE")
  
  #we remove potential outliers
  rm.outliers = romr.fnc(model_lmer, d_data%>%melt(., measure.vars=k)%>%filter(., !is.na(value)),
                         trim=2.5)
  d_data_out = rm.outliers$data
  
  #We first control for all covariates and extract the residuals
  model_lmer=lmer(paste("value ~ (1|Site_ID) + Long_cos + Long_sin + Lattitude + Slope + Elevation + Clim1 + Clim2 + Clim3 + Clim4 + rho_p",
                        ifelse(with_cover=="without_cover","+ rho_p","")),
                  data = d_data_out,
                  na.action = na.fail,REML ="FALSE")
  
  resid_model=residuals(model_lmer) #extract residuals
  
  d_data$Org_C=d_data$MF
  
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

## >> Lasso penalization for multivariate models ----

d_data=read.table("../Data/Spatial_structure_grazing.csv",sep=";")%>%
  Closer_to_normality(.)

d_data[,c(1,7:29,35,37:41,(44:(ncol(d_data))))] = apply(d_data[,c(1,7:29,35,37:41,(44:(ncol(d_data))))],2,z_tranform)

d_data$MF=rowSums(d_data[,55:(ncol(d_data))])

d_data=Perform_PCA_spatial_struc(d_data)

#response variables
responses_var=c("contig","perim_area_scaling",
                "core_area","core_area_land","fmax_psd","PLR","flow_length",
                "Struct1","Struct2")

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

ggsave("../Figures/Final_figs/SI/Importance_lasso_penalized.pdf",p,width = 6,height = 5)



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



## >> Sensitivity cover inferrence ----

## >> 1) How do spatial statistics change with uncertainty of cover ----

d=tibble();thresh=0.01
n_site=100

for (id in 1:n_site){ # Id of the site
  
  for (n_clust in c(2,3,4)){ # number of kmeans clusters
    
    for (binary_thresh in 1:(n_clust-1)){ # binary thresholds
      
      #getting the images and the kept matrix = "true" one 
      cover_site=d_biocom_old$Cover[id]
      img1=Get_png_empirical_site_biocom(id)
      kept_land=Get_empirical_site_biocom(id)
      
      #performing the kmean
      kmean_land=k_means_RGB(img1,n_clust)
      
      #defining the thresholds to binarize the matrix
      if (n_clust>2){
        val_scale=get_cut_grayscale_values(n_clust)[[binary_thresh]]
        site_land=binarize(kmean_land,val_scale[[1]],val_scale[[2]])
      }else {
        if (kept_land[1,1]!=kmean_land[1,1]){ #we order K-means classification
          kmean_land[kmean_land==0]=2
          kmean_land[kmean_land==1]=0
          kmean_land[kmean_land==2]=1
          site_land=kmean_land
        }
      }
      cover_img_binary=sum(site_land==1)/(nrow(site_land)**2) #cover of the binarized image
      
      if ((cover_img_binary-d_biocom_old$Cover[id])**2<thresh){ #we only keep the transformations close in term of cover
        d=rbind(d,Get_sumstat(site_land==1)%>%
                  add_column(., Nclust=n_clust,binary_thresh=binary_thresh,ID_site=id,Type="Changed"))
      } 
    }
  }
}
write.table(d,"../Data/Sensitivity_cover/Varying_threshold.csv",sep=";")

d_tot=rbind(d,d_biocom_old[1:100,14:24]%>%add_column(Nclust=0,binary_thresh=0,
                                                     ID_site=1:100,Type="Kept"))%>%
  arrange(., ID_site,Nclust)%>%
  add_column(., Error_cover=sapply(1:nrow(.),function(x){
    return((.$rho_p[x]-d_biocom_old$rho_p[.$ID_site[x]])**2)
  }))


#Doing a PCA

sumstat_name=colnames(d_tot)[1:11]
res.comp=imputePCA(d_tot[,which(colnames(d_tot) %in% sumstat_name)],ncp=3,scale = T) 

if ("completeObs" %in% names(res.comp)){
  res.pca=PCA(res.comp$completeObs, ncp = 3,  graph=F)
}else {
  res.pca=PCA(res.comp, ncp = 3,  graph=F)
}

axes_for_plot=tibble(x=c(1,1,2),y=c(2,3,3))

for (i in 1:3){
  assign(paste0("p",i),
         d_tot%>%
           mutate(., Error_cover=unlist(sapply(1:nrow(.),function(x){
             if (.$Error_cover[x]==0){
               return(NA)
             }else { 
               return(.$Error_cover[x])}
           })))%>%
           add_column(., PC1=res.pca$ind$coord[,axes_for_plot$x[i]],PC2=res.pca$ind$coord[,axes_for_plot$y[i]])%>%
           ggplot(.) +
           geom_hline(yintercept = 0, lty = 2) +
           geom_vline(xintercept = 0, lty = 2) +
           geom_point(aes(x = PC1, y = PC2, color = Error_cover,fill=Error_cover,shape=Type))+
           geom_line(aes(x = PC1, y = PC2, group=ID_site),color="gray",alpha=.5)+
           scale_color_viridis_c(na.value = "red")+
           scale_fill_viridis_c(na.value = "red")+
           scale_shape_manual(values = c("Kept"=25,"Changed"=21))+
           labs(x=paste0("PC ",axes_for_plot$x[i]," (",round(res.pca$eig[axes_for_plot$x[i],2], 1)," %)"),
                y=paste0("PC ",axes_for_plot$y[i]," (",round(res.pca$eig[axes_for_plot$y[i],2], 1)," %)"),
                color="Quadratic error on cover",fill="")+
           ggtitle("")+guides()+
           theme_classic()+theme(legend.position = "bottom")+
           guides(fill="none",shape="none")+
           theme(legend.box = "vertical")
  )
}

p=ggarrange(ggarrange(p1+theme(legend.position = "none"),
                      p2+theme(legend.position = "none"),
                      p3+theme(legend.position = "none"),
                      ncol=3,align = "hv"),ggarrange(ggplot()+theme_void(),get_legend(p2),ggplot()+theme_void(),ncol=3,widths = c(.3,1,.3)),
            nrow=2,heights = c(1,.2))

ggsave("../Figures/PCA_error.pdf",p,width = 11,height = 5)





## >> 2) Change in parameter inferred ----

NA_kept=100
id=1
d_biocom_old=read.table("../Data/Sensitivity_cover/Varying_threshold.csv",sep=";")
d_sim=read.table("../Data_new/All_new_sim2.csv",sep=";")%>%
  dplyr::relocate(., Pooling,.after =q )%>%
  #filter(., PL_expo>0,!is.na(PLR),ID %in% c(1:15))%>%
  filter(., PL_expo>0,!is.na(PLR))%>%
  dplyr::select(., -ID)

rownames(d_sim)=1:nrow(d_sim)


for (id_plot in id){
  `%!in%` = Negate(`%in%`)
  n_param=3
  
  
  d_param_infer_NN=array(0,c(NA_kept,nrow(d_biocom_old),3))
  d_param_infer_rej=array(0,c(NA_kept,nrow(d_biocom_old),3))
  
  d_NRMSE_sumstat=x_y_stat=tibble()
  
  name_plot=c("all","no_cv","no_cv_fmax","no_PL_PLR","not_the_4","only_fmax","no_PL","no_PLR")
  
  sumstat_kept=4:14
  
  
  
  for (empirical_id in 1:nrow(d_biocom_old)){
    print(empirical_id)
    target=d_biocom_old[empirical_id,sumstat_kept-2]
    matrix_param=d_sim[,1:3]
    
    mat_sumstat=rbind(d_sim[,sumstat_kept],target)
    
    #1) Boxcox
    
    for (x in 1:ncol(mat_sumstat)){
      if (any(is.na(target)) & x %!in% which(is.na(target))){
        if (colnames(mat_sumstat)[x] %in% c("skewness","moran_I","fmax_psd")){
          
          
          
          b=boxcox(lm(mat_sumstat[,x]+abs(min(mat_sumstat[,x]))+.5 ~ 1),plotit = F,eps = .05)     #Working with positive values
          lambda_x=b$x[which.max(b$y)]
          if (lambda_x !=0){ #to avoid errors
            mat_sumstat[,x] = (exp(mat_sumstat[,x]*(lambda_x)) -1)/(lambda_x)
          }
          
        }else {
          b=boxcox(lm(mat_sumstat[,x]+.5 ~ 1),plotit = F,eps = .05)
          lambda_x=b$x[which.max(b$y)]
          if (lambda_x !=0){ #to avoid errors
            mat_sumstat[,x] = (mat_sumstat[,x]^(lambda_x) -1)/(lambda_x)
          }
        }
      }
    }
    
    
    #2) Scaling
    
    for (x in 1:ncol(mat_sumstat)) mat_sumstat[,x] = (mat_sumstat[,x]-mean(mat_sumstat[,x],na.rm = T))/sd(mat_sumstat[,x],na.rm = T)
    
    if (any(is.na(mat_sumstat[nrow(mat_sumstat),]))){
      
      which_na=which(is.na(mat_sumstat[nrow(mat_sumstat),]))
      
      cross_valid=abc(target = mat_sumstat[nrow(mat_sumstat),-which_na],
                      param = matrix_param[-nrow(mat_sumstat),],sumstat = mat_sumstat[-nrow(mat_sumstat),-which_na], #removing the target data
                      tol = 1000/nrow(matrix_param),method = "rejection") #we keep the 1000 closest simulations for the first step
      
    }else {
      cross_valid=abc(target = mat_sumstat[nrow(mat_sumstat),],
                      param = matrix_param[-nrow(mat_sumstat),],sumstat = mat_sumstat[-nrow(mat_sumstat),], #removing the target data
                      tol = 1000/nrow(matrix_param),method = "rejection") #we keep the 1000 closest simulations for the first step
    }
    
    
    #Keeping 1000 simulations and doing the same steps again: normality, scaling and PLS
    
    mat_sumstat_step1=d_sim[as.numeric(rownames(cross_valid$ss)),sumstat_kept] #we keep information with the true values
    mat_sumstat_step1=rbind(mat_sumstat_step1,target)
    
    #again, first box cox
    which_na=which(is.na(mat_sumstat[nrow(mat_sumstat),]))
    
    for (x in 1:ncol(mat_sumstat_step1)){
      
      if (any(which_na) & x %!in% which_na){
        
        if (colnames(mat_sumstat_step1)[x] %in% c("skewness","moran_I","fmax_psd")){
          
          b=boxcox(lm(mat_sumstat_step1[,x]+abs(min(mat_sumstat_step1[,x]))+.5 ~ 1),plotit = F,eps = .05)     #Working with positive values
          lambda_x=b$x[which.max(b$y)]
          
          if (lambda_x !=0){ #to avoid errors
            mat_sumstat_step1[,x] = (exp(mat_sumstat_step1[,x]*(lambda_x)) -1)/(lambda_x)
          }
          
        }else {
          b=boxcox(lm(mat_sumstat_step1[,x] ~ 1),plotit = F,eps = .05)
          lambda_x=b$x[which.max(b$y)]
          if (lambda_x !=0){ #to avoid errors
            mat_sumstat_step1[,x] = (mat_sumstat_step1[,x]^(lambda_x) -1)/(lambda_x)
          }
        }
      }
    }
    
    #and normalization
    for (x in 1:ncol(mat_sumstat_step1)) {
      if (length(unique(mat_sumstat_step1[,x])) != 1){
        mat_sumstat_step1[,x] = (mat_sumstat_step1[,x]-mean(mat_sumstat_step1[,x],na.rm = T))/sd(mat_sumstat_step1[,x],na.rm = T)
      }
    }
    if (any(is.na(mat_sumstat_step1[nrow(mat_sumstat_step1),]))){
      
      which_na=which(is.na(mat_sumstat_step1[nrow(mat_sumstat_step1),]))
      
      cross_valid=abc(target = mat_sumstat_step1[nrow(mat_sumstat_step1),-which_na],
                      param = cross_valid$unadj.values,
                      sumstat = mat_sumstat_step1[-nrow(mat_sumstat_step1),-which_na], #removing the target data
                      tol = NA_kept/nrow(mat_sumstat_step1),method = "rejection",transf = c(rep("logit",2),"none"), #as parameters are proba, we perform logit regression
                      logit.bounds = matrix(c(0,1),3,2,byrow = T),
                      numnet = 10,sizenet = 15)
      
      cross_valid$ss=d_sim[as.numeric(rownames(cross_valid$ss)),sumstat_kept] #we keep information with the true values
      
      if(name_plot[id_plot]=="no_PLR"){
        
        x_y_stat=rbind(x_y_stat,target[-which_na]%>%
                         add_column(., PL_expo=NA, Site_ID=empirical_id,Method="rejection",Type="Obs")%>%
                         relocate(., PL_expo,.after =Spectral_ratio ))
        
        x_y_stat=rbind(x_y_stat,as_tibble(t(colMeans(cross_valid$ss)[-which_na]))%>%
                         add_column(.,PL_expo=NA, Site_ID=empirical_id,Method="rejection",Type="Sim")%>%
                         relocate(., PL_expo,.after =Spectral_ratio ))
      }else{
        
        x_y_stat=rbind(x_y_stat,target[-which_na]%>%
                         add_column(., PL_expo=NA, Site_ID=empirical_id,Method="rejection",Type="Obs")%>%
                         relocate(., PL_expo,.after =PLR ))
        
        x_y_stat=rbind(x_y_stat,as_tibble(t(colMeans(cross_valid$ss)[-which_na]))%>%
                         add_column(.,PL_expo=NA, Site_ID=empirical_id,Method="rejection",Type="Sim")%>%
                         relocate(., PL_expo,.after =PLR ))
        
        
      }
      
      
    }else {
      cross_valid=abc(target = mat_sumstat_step1[nrow(mat_sumstat_step1),],
                      param = cross_valid$unadj.values,
                      sumstat = mat_sumstat_step1[-nrow(mat_sumstat_step1),], #removing the target data
                      tol = NA_kept/nrow(mat_sumstat_step1),method = "rejection",transf = c(rep("logit",2),"none"), #as parameters are proba, we perform logit regression
                      logit.bounds = matrix(c(0,1),3,2,byrow = T),
                      numnet = 10,sizenet = 15)
      
      cross_valid$ss=d_sim[as.numeric(rownames(cross_valid$ss)),sumstat_kept] #we keep information with the true values
      
      x_y_stat=rbind(x_y_stat,target%>%
                       add_column(., Site_ID=empirical_id,Method="rejection",Type="Obs"))
      
      x_y_stat=rbind(x_y_stat,as_tibble(t(colMeans(cross_valid$ss)))%>%
                       add_column(.,Site_ID=empirical_id,Method="rejection",Type="Sim"))
      
    }
    
    
    cross_valid$ss=d_sim[as.numeric(rownames(cross_valid$ss)),sumstat_kept] #we keep information with the true values
    
    mat_sumstat=d_sim[,sumstat_kept]
    
    cross_valid$adj.values=cross_valid$adj.values
    
    #NRMSE for the summary statistics observed
    RMSE = sapply(1:ncol(cross_valid$ss),function(x){
      sqrt(sum((cross_valid$ss[,x]-target[,x])**2,na.rm = T)/nrow(cross_valid$ss) )
    }
    )
    
    RMSE_prior=sapply(1:ncol(mat_sumstat),function(x){
      sqrt(sum((mat_sumstat[,x]-target[,x])**2,na.rm = T)/nrow(mat_sumstat) )
    }
    )
    NRMSE = RMSE/RMSE_prior
    
    d_NRMSE_sumstat=rbind(d_NRMSE_sumstat,as_tibble(t(NRMSE)))
    
    d_param_infer_rej[,empirical_id,1]=cross_valid$unadj.values[,1] # we keep the whole distribution for p
    d_param_infer_rej[,empirical_id,2]=cross_valid$unadj.values[,2] # for q
    d_param_infer_rej[,empirical_id,3]=cross_valid$unadj.values[,3] # for the scale of observation
    
    
  }
  
  write.table(d_NRMSE_sumstat,paste0("../Data/Sensitivity_cover/NRMSE_sumstat.csv"),sep=";")
  
  write.table(x_y_stat,paste0("../Data/Sensitivity_cover/x_y_stat.csv"),sep=";")
  
  write.table(d_param_infer_rej,paste0("../Data/Sensitivity_cover/param_inferred.csv"),sep=";")
}


nsite=100
param_infered=read.table("../Data/Sensitivity_cover/param_inferred.csv",sep=";")
param_infered_biocom=read.table("../Data/Sensitivity_cover/posterior_param.csv",sep=";")
d_threshold=read.table("../Data/Sensitivity_cover/Varying_threshold.csv",sep=";")

d_inferrence=tibble(p=c(colMeans(param_infered_biocom[1:nsite]),colMeans(param_infered)[1:nrow(d_threshold)]),
                    q=c(colMeans(param_infered_biocom[c(1:nsite)+345]),colMeans(param_infered)[(nrow(d_threshold)+1):(2*nrow(d_threshold))]),
                    Site=c(1:nsite,d_threshold$ID_site),
                    Type=c(rep(1,nsite),d_threshold$Nclust))


p=ggplot(d_inferrence%>%melt(.,measure.vars=c("p","q")))+
  geom_line(aes(Type,value,group=Site),color="gray",lwd=.5,alpha=.5)+
  facet_wrap(.~variable)+
  scale_x_continuous(labels = c("True","Modif_1","Modif_2","Modif_3"))+
  the_theme


ggsave("../Figures/Change_parameters_no_info_cover.pdf",p,width = 8,height = 4)




## >> Importance with climatic variables ----


#Importance using climatic variables 

d_all=read.table(paste0("../Data/Linear_models/Importance_no_inter.csv"),sep=";")%>%
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
  labs(x="Predictors",y="Spatial statistic",fill="Frequency in best models")

ggsave(paste0("../Figures/SI/Importance_grazing_no_inter.pdf"),p,width = 6,height = 5)





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
    
    d_data_out=  read.table(paste0("../Data/Linear_models/Keep_data/Data_",stat,"_TRUE.csv"),sep=";")
    if (ncol(d_data_out)==1){
      d_data_out=  read.table(paste0("../Data/Linear_models/Keep_data/Data_",stat,"_TRUE.csv"),sep=" ")
    }
    model_spa_stat=readRDS(paste0("../Data/Linear_models/Keep_models/Mod_",stat,"_TRUE.rds"))
    
    
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
  
  
  ggsave(paste0("../Figures/SI/All_interactions_",stat,".pdf"),p,width = 4*length(list_plot),height = 4)
  
  ID_title=ID_title+1
  
}















## >> SEM_no_woody_2_soil_vars ----

rm(list=ls())
source("./Structure_grazing_function.R")

Plot_SEM_soil_var=function(summary_sem,pdf_=F,name="SEM",title_="",name_var="",MF=F,pval=.1){
  
  l=summary_sem$coefficients[,c("Predictor","Response","Estimate","P.Value")]
  l$color="#9ED471"
  l$color[which(l$Estimate<0)]="#EFC46A"
  l$color[which(l$P.Value>pval)]="gray"
  g = graph.data.frame(l, directed=T)
  g= g %>% set_edge_attr("color", value =l$color)
  
  coord=data.frame(label=c("Grazing","Amonium","Organic_C","Nitrate","Cover",paste0(name_var)),
                   lab2=c("Grazing","Amonium","Organic_C","Nitrate","Cover",paste0(name_var)),
                   x=c(-40,-20,-20,-20,-20,0),y=c(0,20,40,60,-40,0))
  
  EL=as_edgelist(g)
  EL=cbind(EL,l[,3])
  EL=as.data.frame(EL)%>%
    mutate(., V1=recode_factor(V1,"rho_p"="Cover","Org_C"="Organic_C","Resid_mod"=paste0(name_var)))%>%
    mutate(., V2=recode_factor(V2,"rho_p"="Cover","Org_C"="Organic_C","Resid_mod"=paste0(name_var)))
  EL=as.matrix(EL)
  
  name_edge=sapply(1:length(summary_sem$coefficients$P.Value),function(x){#adding pvalue
    if (is_signif(summary_sem$coefficients$P.Value[x])==""){
      return("")
    }else{
      return(paste0(round(l[x,3],2))) #,is_signif(summary_sem$coefficients$P.Value[x])))
    }
  })
  
  asi=abs(l[,3])/0.05
  asi[asi<5]=5
  
  if(pdf_){
    pdf(paste0("./",name,".pdf"),width = 7,height = 4)
    qgraph(EL,layout=as.matrix(coord[,c("x","y")]),edge.color=l$color,edge.labels=name_edge,
           border.color="white",label.cex=1.5,label.scale=F,title=title_,
           edge.label.cex = 1.2,edge.label.position=0.5,vsize2=4,vsize=30,
           shape="ellipse",edge.labels=T,fade=F,esize=5,asize=asi,
           mar=rep(3,4))
    dev.off()
    
  }else {
    qgraph(EL,layout=as.matrix(coord[,c("x","y")]),edge.color=l$color,edge.labels=name_edge,
           border.color="white",label.cex=1,label.scale=F,title=title_,
           edge.label.cex = 3,edge.label.position=0.5,vsize2=10,vsize=30,
           shape="ellipse",edge.labels=T,fade=F,esize=5,asize=asi,
           mar=rep(3,4))
    
  }
}

dir.create("../Figures/Linear_models/SEM4",showWarnings = F)

d_data=read.table("../Data/Spatial_structure_grazing.csv",sep=";")%>%
  Closer_to_normality(.)

d_data[,c(1,9:29,36,38:42,45:ncol(d_data))] = apply(d_data[,c(1,9:29,36,38:42,45:ncol(d_data))],2,z_tranform)

#Getting the multifunctionality index
d_data$MF=z_tranform(rowSums(d_data[,56:(ncol(d_data))],na.rm = T))

d_data=Perform_PCA_spatial_struc(d_data) #Adding the first two components from the PCA

param_list=expand.grid(Stat=c("Struct1","Struct2"),
                       MF=c(F)) #MF instead of organic carbon

d_indirect=tibble()

for (ID in 1:nrow(param_list)){
  
  
  k=as.character(param_list$Stat[ID])
  
  model_lmer=lm("value ~ Long_cos + Long_sin + Lattitude + Slope + Elevation +
                        Aridity + Sand + Type_veg",
                data = d_data%>%melt(., measure.vars=k)%>%filter(., !is.na(value)),
                na.action = na.omit)
  
  #we remove potential outliers
  rm.outliers = romr.fnc(model_lmer, d_data%>%melt(., measure.vars=k)%>%filter(., !is.na(value)),
                         trim=2.5)
  d_data_out = rm.outliers$data
  
  #We first control for all covariates and extract the residuals
  model_lmer=lm("value ~ Long_cos + Long_sin + Lattitude + Slope + Elevation +
                        Aridity + Sand + Type_veg",
                data = d_data_out,
                na.action = na.fail)
  
  resid_model=residuals(model_lmer) #extract residuals
  
  save=d_data[as.numeric(names(resid_model)),]%>% #add it to the dataframe 
    add_column(., Resid_mod=resid_model)
  
  #DOING the SEMs
  
  
  
  #SEM with all grazing intensity
  
  d_sem=save
  all_d=summary(psem(
    lm(Resid_mod ~  Grazing + Nitrate + Org_C + Amonium + rho_p, d_sem),
    lm(Nitrate ~  Grazing, d_sem),
    lm(Org_C ~  Grazing, d_sem),
    lm(rho_p ~  Grazing, d_sem),
    lm(Amonium ~  Grazing, d_sem)
  ))
  
  #SEM with low grazing intensity
  
  d_sem=filter(save,Grazing %in% 0:1)
  low_graz=summary(psem(
    lm(Resid_mod ~  Grazing + Nitrate + Org_C + Amonium + rho_p, d_sem),
    lm(Nitrate ~  Grazing, d_sem),
    lm(Org_C ~  Grazing, d_sem),
    lm(rho_p ~  Grazing, d_sem),
    lm(Amonium ~  Grazing, d_sem)
  ))
  
  #SEM with high grazing intensity
  
  d_sem=filter(save,Grazing %in% 2:3)
  high_graz=summary(psem(
    lm(Resid_mod ~  Grazing + Nitrate + Org_C + Amonium + rho_p, d_sem),
    lm(Nitrate ~  Grazing, d_sem),
    lm(Org_C ~  Grazing, d_sem),
    lm(rho_p ~  Grazing, d_sem),
    lm(Amonium ~  Grazing, d_sem)
  ))
  
  
  Plot_SEM_soil_var(summary_sem = low_graz,pdf_ = T,
                    name = paste0("../Figures/Linear_models/SEM4/",
                                  "SEM_type_veg_",param_list$Stat[ID],"_",param_list$Cov[ID],"_low"),
                    title_ = "Low grazing pressure",
                    name_var = param_list$Stat[ID],MF = param_list$MF[ID])
  
  Plot_SEM_soil_var(summary_sem = high_graz,pdf_ = T,
                    name = paste0("../Figures/Linear_models/SEM4/",
                                  "SEM_type_veg_",param_list$Stat[ID],"_",param_list$Cov[ID],"_high"),
                    title_ = "High grazing pressure",
                    name_var = param_list$Stat[ID],MF = param_list$MF[ID])
  
  Plot_SEM_soil_var(summary_sem = all_d,pdf_ = T,
                    name = paste0("../Figures/Linear_models/SEM4/",
                                  "SEM_type_veg_",param_list$Stat[ID],"_",param_list$Cov[ID],"_all"),
                    title_ = "All grazing pressure",
                    name_var = param_list$Stat[ID],MF = param_list$MF[ID])
  
}
## >> SEM_woody_1_soil_vars----
rm(list=ls())
source("./Structure_grazing_function.R")


d_data=read.table("../Data/Spatial_structure_grazing.csv",sep=";")%>%
  Closer_to_normality(.)

d_data[,c(1,9:29,36,38:42,45:ncol(d_data))] = apply(d_data[,c(1,9:29,36,38:42,45:ncol(d_data))],2,z_tranform)

#Getting the multifunctionality index
d_data$MF=z_tranform(rowSums(d_data[,56:(ncol(d_data))],na.rm = T))

d_data=Perform_PCA_spatial_struc(d_data) #Adding the first two components from the PCA

param_list=expand.grid(Stat=c("Struct1","Struct2",
                              "perim_area_scaling"),
                       MF=c(F),
                       soil_var=c("Nitrate","Amonium")) #MF instead of organic carbon

d_indirect=tibble()

for (ID in 1:nrow(param_list)){
  
  
  k=as.character(param_list$Stat[ID])
  
  model_lmer=lm("value ~ Long_cos + Long_sin + Lattitude + Slope + Elevation +
                        Aridity + Sand + Type_veg",
                data = d_data%>%melt(., measure.vars=k)%>%filter(., !is.na(value)),
                na.action = na.omit)
  
  #we remove potential outliers
  rm.outliers = romr.fnc(model_lmer, d_data%>%melt(., measure.vars=k)%>%filter(., !is.na(value)),
                         trim=2.5)
  d_data_out = rm.outliers$data
  
  #We first control for all covariates and extract the residuals
  model_lmer=lm("value ~ Long_cos + Long_sin + Lattitude + Slope + Elevation +
                        Aridity + Sand + Type_veg",
                data = d_data_out,
                na.action = na.fail)
  
  resid_model=residuals(model_lmer) #extract residuals
  
  save=d_data[as.numeric(names(resid_model)),]%>% #add it to the dataframe 
    add_column(., Resid_mod=resid_model)
  
  #DOING the SEMs
  if (param_list$soil_var[ID]=="Nitrate"){
    save$soil_var=save$Nitrate
  }else{
    save$soil_var=save$Amonium
  }
  
  
  #SEM with all grazing intensity
  
  d_sem=save
  all_d=summary(psem(
    lm(Resid_mod ~  Grazing + Woody + Org_C + soil_var + rho_p, d_sem),
    lm(Woody ~  Grazing , d_sem),
    lm(Org_C ~  Grazing, d_sem),
    lm(rho_p ~  Grazing, d_sem),
    lm(soil_var ~  Grazing, d_sem)
  ))
  
  #SEM with low grazing intensity
  
  d_sem=filter(save,Grazing %in% 0:1)
  low_graz=summary(psem(
    lm(Resid_mod ~  Grazing + Woody + Org_C + soil_var + rho_p, d_sem),
    lm(Woody ~  Grazing , d_sem),
    lm(Org_C ~  Grazing, d_sem),
    lm(rho_p ~  Grazing, d_sem),
    lm(soil_var ~  Grazing, d_sem)
  ))
  
  #SEM with high grazing intensity
  
  d_sem=filter(save,Grazing %in% 2:3)
  high_graz=summary(psem(
    lm(Resid_mod ~  Grazing + Woody + Org_C + soil_var + rho_p, d_sem),
    lm(Woody ~  Grazing , d_sem),
    lm(Org_C ~  Grazing, d_sem),
    lm(rho_p ~  Grazing, d_sem),
    lm(soil_var ~  Grazing, d_sem)
  ))
  
  
  Plot_SEM_1_soil_var(summary_sem = low_graz,pdf_ = T,
                      name = paste0("../Figures/Linear_models/SEM2/",
                                    "SEM_type_veg_",param_list$Stat[ID],"_",param_list$Cov[ID],"_low_",param_list$soil_var[ID]),
                      title_ = "Low grazing pressure",
                      name_var = param_list$Stat[ID],MF = param_list$MF[ID])
  
  Plot_SEM_1_soil_var(summary_sem = high_graz,pdf_ = T,
                      name = paste0("../Figures/Linear_models/SEM2/",
                                    "SEM_type_veg_",param_list$Stat[ID],"_",param_list$Cov[ID],"_high_",param_list$soil_var[ID]),
                      title_ = "High grazing pressure",
                      name_var = param_list$Stat[ID],MF = param_list$MF[ID])
  
  Plot_SEM_1_soil_var(summary_sem = all_d,pdf_ = T,
                      name = paste0("../Figures/Linear_models/SEM2/",
                                    "SEM_type_veg_",param_list$Stat[ID],"_",param_list$Cov[ID],"_all_",param_list$soil_var[ID]),
                      title_ = "All grazing pressure",
                      name_var = param_list$Stat[ID],MF = param_list$MF[ID])
  
}
## >> SEM_no_woody_1_soil_var----
rm(list=ls())
source("./Structure_grazing_function.R")

dir.create("../Figures/Linear_models/SEM3",showWarnings = F)

Plot_SEM_1_soil_var=function(summary_sem,pdf_=F,name="SEM",title_="",name_var="",MF=F,pval=.1,name_soil_var="Nitrate"){
  
  l=summary_sem$coefficients[,c("Predictor","Response","Estimate","P.Value")]
  l$color="#9ED471"
  l$color[which(l$Estimate<0)]="#EFC46A"
  l$color[which(l$P.Value>pval)]="gray"
  g = graph.data.frame(l, directed=T)
  g= g %>% set_edge_attr("color", value =l$color)
  
  coord=data.frame(label=c("Grazing","rho_p","Organic_C",name_soil_var,paste0(name_var)),
                   lab2=c("Grazing","rho_p","Organic_C",name_soil_var,paste0(name_var)),
                   x=c(-40,-20,-20,-20,0),y=c(0,-20,20,40,0))
  
  EL=as_edgelist(g)
  EL=cbind(EL,l[,3])
  EL=as.data.frame(EL)%>%
    mutate(., V1=recode_factor(V1,"Org_C"="Organic_C","Resid_mod"=paste0(name_var)))%>%
    mutate(., V2=recode_factor(V2,"Org_C"="Organic_C","Resid_mod"=paste0(name_var)))
  EL=as.matrix(EL)
  
  name_edge=sapply(1:length(summary_sem$coefficients$P.Value),function(x){#adding pvalue
    if (is_signif(summary_sem$coefficients$P.Value[x])==""){
      return("")
    }else{
      return(paste0(round(l[x,3],2))) #,is_signif(summary_sem$coefficients$P.Value[x])))
    }
  })
  
  asi=abs(l[,3])/0.05
  asi[asi<5]=5
  
  if(pdf_){
    pdf(paste0("./",name,".pdf"),width = 7,height = 4)
    qgraph(EL,layout=as.matrix(coord[,c("x","y")]),edge.color=l$color,edge.labels=name_edge,
           border.color="white",label.cex=1.5,label.scale=F,title=title_,
           edge.label.cex = 1.2,edge.label.position=0.5,vsize2=4,vsize=30,
           shape="ellipse",edge.labels=T,fade=F,esize=5,asize=asi,
           mar=rep(3,4))
    dev.off()
    
  }else {
    qgraph(EL,layout=as.matrix(coord[,c("x","y")]),edge.color=l$color,edge.labels=name_edge,
           border.color="white",label.cex=1,label.scale=F,title=title_,
           edge.label.cex = 3,edge.label.position=0.5,vsize2=10,vsize=30,
           shape="ellipse",edge.labels=T,fade=F,esize=5,asize=asi,
           mar=rep(3,4))
    
  }
}


d_data=read.table("../Data/Spatial_structure_grazing.csv",sep=";")%>%
  Closer_to_normality(.)

d_data[,c(1,9:29,36,38:42,45:ncol(d_data))] = apply(d_data[,c(1,9:29,36,38:42,45:ncol(d_data))],2,z_tranform)

#Getting the multifunctionality index
d_data$MF=z_tranform(rowSums(d_data[,56:(ncol(d_data))],na.rm = T))

d_data=Perform_PCA_spatial_struc(d_data) #Adding the first two components from the PCA

param_list=expand.grid(Stat=c("Struct1","Struct2"),
                       MF=c(F),
                       soil_var=c("Nitrate","Amonium")) #MF instead of organic carbon

d_indirect=tibble()

for (ID in 1:nrow(param_list)){
  
  
  k=as.character(param_list$Stat[ID])
  
  model_lmer=lm("value ~ Long_cos + Long_sin + Lattitude + Slope + Elevation +
                        Aridity + Sand + Type_veg",
                data = d_data%>%melt(., measure.vars=k)%>%filter(., !is.na(value)),
                na.action = na.omit)
  
  #we remove potential outliers
  rm.outliers = romr.fnc(model_lmer, d_data%>%melt(., measure.vars=k)%>%filter(., !is.na(value)),
                         trim=2.5)
  d_data_out = rm.outliers$data
  
  #We first control for all covariates and extract the residuals
  model_lmer=lm("value ~ Long_cos + Long_sin + Lattitude + Slope + Elevation +
                        Aridity + Sand + Type_veg",
                data = d_data_out,
                na.action = na.fail)
  
  resid_model=residuals(model_lmer) #extract residuals
  
  save=d_data[as.numeric(names(resid_model)),]%>% #add it to the dataframe 
    add_column(., Resid_mod=resid_model)
  
  #DOING the SEMs
  if (param_list$soil_var[ID]=="Nitrate"){
    save$soil_var=save$Nitrate
  }else{
    save$soil_var=save$Amonium
  }
  
  
  #SEM with all grazing intensity
  
  d_sem=save
  all_d=summary(psem(
    lm(Resid_mod ~  Grazing + rho_p + Org_C + soil_var,d_sem),
    lm(rho_p ~  Grazing , d_sem),
    lm(Org_C ~  Grazing, d_sem),
    lm(soil_var ~  Grazing, d_sem)
  ))
  
  #SEM with low grazing intensity
  
  d_sem=filter(save,Grazing %in% 0:1)
  low_graz=summary(psem(
    lm(Resid_mod ~  Grazing + rho_p + Org_C + soil_var, d_sem),
    lm(rho_p ~  Grazing , d_sem),
    lm(Org_C ~  Grazing, d_sem),
    lm(soil_var ~  Grazing, d_sem)
  ))
  
  #SEM with high grazing intensity
  
  d_sem=filter(save,Grazing %in% 2:3)
  high_graz=summary(psem(
    lm(Resid_mod ~  Grazing + rho_p + Org_C + soil_var, d_sem),
    lm(rho_p ~  Grazing , d_sem),
    lm(Org_C ~  Grazing, d_sem),
    lm(soil_var ~  Grazing, d_sem)
  ))
  
  
  Plot_SEM_1_soil_var(summary_sem = low_graz,pdf_ = T,
                      name = paste0("../Figures/Linear_models/SEM3/",
                                    "SEM_type_veg_",param_list$Stat[ID],"_",param_list$Cov[ID],"_low_",param_list$soil_var[ID]),
                      title_ = "Low grazing pressure",
                      name_var = param_list$Stat[ID],MF = param_list$MF[ID])
  
  Plot_SEM_1_soil_var(summary_sem = high_graz,pdf_ = T,
                      name = paste0("../Figures/Linear_models/SEM3/",
                                    "SEM_type_veg_",param_list$Stat[ID],"_",param_list$Cov[ID],"_high_",param_list$soil_var[ID]),
                      title_ = "High grazing pressure",
                      name_var = param_list$Stat[ID],MF = param_list$MF[ID])
  
  Plot_SEM_1_soil_var(summary_sem = all_d,pdf_ = T,
                      name = paste0("../Figures/Linear_models/SEM3/",
                                    "SEM_type_veg_",param_list$Stat[ID],"_",param_list$Cov[ID],"_all_",param_list$soil_var[ID]),
                      title_ = "All grazing pressure",
                      name_var = param_list$Stat[ID],MF = param_list$MF[ID])
  
}
## >> 3) Spatial metrics on simulations ----

#simulating and computing spatial statistics
param_space=rbind(tibble(b=seq(0,.8,.03),g0=0,Driver="Aridity"),
                  tibble(b=.8,g0=seq(0,0.5,.02),Driver="Grazing"))

d_sim=tibble()
param=Get_classical_param(g0=0)

for (k in 1:nrow(param_space)){
  print(k)
  
  param$g0=param_space$g0[k]
  param$b=param_space$b[k]
  
  ini_land=Get_initial_lattice()
  out_model=Run_spatial_model(params = param,ini=ini_land)
  
  d_sim=rbind(d_sim,Get_sumstat(out_model$State==1)%>%
                add_column(., g0=param_space$g0[k],b=param_space$b[k],Driver=param_space$Driver[k]))
}
write.table(d_sim,"../Data/Spatial_structure_simulations.csv",sep=";")


d_sim=read.table("../Data/Spatial_structure_simulations.csv",sep=";")
d_slope=tibble()
for (k in c("perim_area_scaling","fmax_psd","PL_expo","Spectral_ratio",
            "core_area_land","division","contig","core_area",
            "flow_length","PLR")){
  
  data_mod_aridity=d_sim%>%filter(., Driver=="Aridity",rho_p>0.01)%>%melt(., measure.vars=k)%>%
    filter(., !is.na(value))%>%dplyr::select(., value,b,rho_p)
  data_mod_grazing=d_sim%>%filter(., Driver=="Grazing",rho_p>0.01)%>%melt(., measure.vars=k)%>%
    filter(., !is.na(value))%>%dplyr::select(., value,g0,rho_p)
  
  data_mod_grazing=as.data.frame(apply(data_mod_grazing,2,scale))
  data_mod_aridity=as.data.frame(apply(data_mod_aridity,2,scale))
  
  model_spa_stat_aridity=lm("value ~ b + rho_p",data=data_mod_aridity)
  model_spa_stat_grazing=lm("value ~ g0 + rho_p",data=data_mod_grazing)
  
  resid_mod=visreg::visreg(fit = model_spa_stat_grazing,xvar="g0",plot=F) 
  mod_cov=lm(visregRes~g0,resid_mod$res)
  
  d_slope=rbind(d_slope,
                tibble(slope=summary(mod_cov)$coefficient[2,1],
                       Low_int=confint(mod_cov)[2,1],
                       High_int=confint(mod_cov)[2,2],
                       With_cover=T,
                       Stat=k,Driver="Grazing (g0)"
                ))
  
  resid_mod=visreg::visreg(fit = model_spa_stat_aridity,xvar="b",plot=F) 
  mod_cov=lm(visregRes~b,resid_mod$res)
  
  #as increase in aridity leads to the decrease of b in the model, we take the opposite sign of the effects
  
  d_slope=rbind(d_slope,
                tibble(slope=-summary(mod_cov)$coefficient[2,1],
                       Low_int=-confint(mod_cov)[2,1],
                       High_int=-confint(mod_cov)[2,2],
                       With_cover=T,
                       Stat=k,Driver="Aridity (b)"
                ))
  
  
  model_spa_stat_aridity=lm("value ~ b",data=data_mod_aridity)
  model_spa_stat_grazing=lm("value ~ g0",data=data_mod_grazing)
  
  resid_mod=visreg::visreg(fit = model_spa_stat_grazing,xvar="g0",plot=F) 
  mod_cov=lm(visregRes~g0,resid_mod$res)
  
  d_slope=rbind(d_slope,
                tibble(slope=summary(mod_cov)$coefficient[2,1],
                       Low_int=confint(mod_cov)[2,1],
                       High_int=confint(mod_cov)[2,2],
                       With_cover=F,
                       Stat=k,Driver="Grazing (g0)"
                ))
  
  resid_mod=visreg::visreg(fit = model_spa_stat_aridity,xvar="b",plot=F) 
  mod_cov=lm(visregRes~b,resid_mod$res)
  
  d_slope=rbind(d_slope,
                tibble(slope=-summary(mod_cov)$coefficient[2,1],
                       Low_int=-confint(mod_cov)[2,1],
                       High_int=-confint(mod_cov)[2,2],
                       With_cover=F,
                       Stat=k,Driver="Aridity (b)"
                ))
  
}

d_slope2=d_slope%>%
  filter(., With_cover==1,Driver=="Aridity (b)")%>%
  Filter_relevant_stats(.)%>%  
  Rename_spatial_statistics(.)%>%
  arrange(., slope)%>%
  add_column(., Order_f=1:nrow(.))%>%
  mutate(.,Stat = fct_reorder(Stat, Order_f))

p1=ggplot(d_slope2)+
  geom_pointrange(aes(x=slope,y=Stat,xmin=Low_int,xmax=High_int),shape=17,size=.8)+
  the_theme+
  scale_color_manual(values=c("grey","black"))+
  #geom_text(aes(x=order_signif,y=1:8,label=Signif))+
  geom_vline(xintercept = 0,linetype=9)+
  labs(x=substitute(paste(beta," (partial res. Grazing)")),y="",color="")+
  theme(legend.position = "none")

d_slope2=d_slope%>%
  filter(., With_cover==1,Driver=="Grazing (g0)")%>%
  Filter_relevant_stats(.)%>%  
  Rename_spatial_statistics(.)%>%
  arrange(., slope)%>%
  add_column(., Order_f=1:nrow(.))%>%
  mutate(.,Stat = fct_reorder(Stat, Order_f))

p2=ggplot(d_slope2)+
  geom_pointrange(aes(x=slope,y=Stat,xmin=Low_int,xmax=High_int),shape=17,size=.8)+
  the_theme+
  scale_color_manual(values=c("grey","black"))+
  #geom_text(aes(x=order_signif,y=1:8,label=Signif))+
  geom_vline(xintercept = 0,linetype=9)+
  labs(x=substitute(paste(beta," (partial res. Aridity)")),y="",color="")+
  theme(legend.position = "none")

ggarrange(p1,p2)


## >> Plot relative importance pred using standardized effects ----
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


## >> 1) With all climatic variables ----

dir.create("../Data/Linear_models/VIF",showWarnings = F)
dir.create("../Data/Linear_models/Auto_corr",showWarnings = F)
dir.create("../Data/Linear_models/Importance",showWarnings = F)
dir.create("../Data/Linear_models/Estimator",showWarnings = F)
dir.create("../Data/Linear_models/Keep_models",showWarnings = F)
dir.create("../Data/Linear_models/Keep_data",showWarnings = F)


Run_model_importance=function(id){
  
  list_mod=expand.grid(with_cover=c(T,F),
                       Stats=c("perim_area_scaling","fmax_psd","Cond_H","PL_expo","Spectral_ratio",
                               "core_area_land","division","fractal_dim","contig","core_area",
                               "flow_length","PLR","KS_dist",
                               "Struct1","Struct2")) #We also add the case where the stats correspond to the first two components of the PCA
  
  d_all=d_all2=tibble()#for keeping all information
  
  d_data=read.table("../Data/Spatial_structure_grazing.csv",sep=";")%>%
    Closer_to_normality(.)
  
  d_data[,c(1,9:29,36,38:42,45:ncol(d_data))] = 
    apply(d_data[,c(1,9:29,36,38:42,45:ncol(d_data))],2,z_tranform)
  
  d_data=Perform_PCA_spatial_struc(d_data)
  
  stat=list_mod$Stats[id]
  
  d_data_mod=d_data%>%melt(., measure.vars=stat)%>%
    filter(., !is.na(value))
  
  with_cover=list_mod$with_cover[id]
  
  if (with_cover){
    
    formula_mod=formula(formula_(paste("value ~ Long_cos + Long_sin + Lattitude + Slope + Elevation 
      + Clim1 + Clim2 + Clim3 + Clim4 + Grazing 
      + Woody + Woody*Grazing
      + Clim1*Grazing + Clim2*Grazing + Clim3*Grazing + Clim4*Grazing  
      + rho_p + rho_p*Grazing + Type_veg + Type_veg*Grazing
      + Sand + Sand*Grazing + Org_C_v + Org_C_v*Grazing
      + (1|Site_ID)")))
    
  } else{
    
    formula_mod=formula(formula_(paste("value ~ Long_cos + Long_sin + Lattitude + Slope + Elevation 
      + Clim1 + Clim2 + Clim3 + Clim4 + Grazing 
      + Woody + Woody*Grazing
      + Clim1*Grazing + Clim2*Grazing + Clim3*Grazing + Clim4*Grazing  
      + Type_veg + Type_veg*Grazing
      + Sand + Sand*Grazing + Org_C_v + Org_C_v*Grazing
      + (1|Site_ID)")))
  }
  
  model_spa_stat  = lmer(formula_mod, d_data_mod,
                         na.action = na.fail ,REML ="FALSE")
  
  #we remove potential outliers
  rm.outliers = romr.fnc(model_spa_stat, d_data_mod, trim=2.5)
  d_data_out = rm.outliers$data
  
  #saving data
  write.table(d_data_out,paste0("../Data/Linear_models/Keep_data/Data_",stat,"_",with_cover,".csv"))
  
  model_spa_stat = lmer(formula_mod, data = d_data_out, 
                        na.action = na.fail,REML ="FALSE")
  
  # #do some model selection
  select_model=dredge(model_spa_stat, subset = ~ Slope & Elevation &
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
                      extra=c(R2m=function(x) r.squaredGLMM(x)[1],
                              R2c=function(x) r.squaredGLMM(x)[2]),
                      options(na.action = "na.fail") )
  
  #extract the result of model selection
  result_select=model.avg(select_model, subset = delta < 2)
  
  #Get the importance of each metric
  importance_mod=sw(result_select)
  
  saveRDS(model_spa_stat,paste0("../Data/Linear_models/Keep_models/Mod_",stat,"_",with_cover,".rds"))
  
  R2=select_model%>%filter(., AICc<min(AICc)+2)
  
  d_all=rbind(d_all,cbind(tibble(Sp_stat=stat,
                                 N_outliers=rm.outliers$n.removed,
                                 R2m=mean(R2$R2m),R2C=mean(R2$R2c)),
                          Aggregate_importance(importance_mod))%>%
                add_column(., With_cover=with_cover))
  
  
  
  summary_coef=confint(result_select)
  
  #Merge in a df
  d_all2=tibble(Median=confint(result_select,level = 1e-9)[-1,1],
                q1=summary_coef[-1,1], 
                q3=summary_coef[-1,2],
                pvalue=summary(result_select)$coefmat.full[-1,5],
                term=rownames(summary_coef)[-1],
                Stat=stat,
                R2m=mean(R2$R2m),
                R2c=mean(R2$R2c))%>%
    add_column(., With_cover=with_cover)
  
  
  
  #Checking multicolinearity using VIF 
  
  vif_model=vif(model_spa_stat)
  write.table(vif_model,paste0("../Data/Linear_models/VIF/VIF_",stat,"_cover_",with_cover,".csv"),sep=";")
  
  
  #Checking spatial auto-correlation at different spatial scales (10, 30, 50)
  
  d_auto_corr=tibble()
  for (spatial_scales in c(10,20,50)){
    
    #Aggregating at different scales
    spatial_coords = cbind(X=as.numeric(d_data_out$Longitude), 
                           Y=as.numeric(d_data_out$Lattitude))
    
    KNN = knearneigh(spatial_coords, k=spatial_scales) 
    nb_KNN=knn2nb(KNN)
    
    #getting residuals
    resids_mod=residuals(model_spa_stat)
    
    #Testing for auto-correlation
    Moran_test=moran.test(resids_mod, nb2listw(nb_KNN))
    
    d_auto_corr=rbind(d_auto_corr,tibble(K=spatial_scales,
                                         Moran_stat=Moran_test$statistic,
                                         p_val=Moran_test$p.value))
    
  }
  
  # Saving autocorrelation checks
  write.table(d_auto_corr,paste0("../Data/Linear_models/Auto_corr/Auto_corr_",stat,"_cover_",with_cover,".csv"),sep=";")
  
  #Saving each DF
  write.table(d_all,paste0("../Data/Linear_models/Importance/Importance_",stat,"_",with_cover,".csv"),sep=";")
  write.table(d_all2,paste0("../Data/Linear_models/Estimator/Estimators_model_",stat,"_",with_cover,".csv"),sep=";")
  
}

library(parallel)

mclapply(1:30,Run_model_importance,mc.cores = 30)


d_all=d_all2=tibble()
for (k in list.files("../Data/Linear_models/Importance","E.csv")){
  d_all=rbind(d_all,read.table(paste0("../Data/Linear_models/Importance/",k),sep=";"))
}
for (k in list.files("../Data/Linear_models/Estimator","E.csv")){
  d_all2=rbind(d_all2,read.table(paste0("../Data/Linear_models/Estimator/",k),sep=";"))
}
write.table(d_all,"../Data/Linear_models/Importance.csv",sep=";")
write.table(d_all2,"../Data/Linear_models/Estimators_model.csv",sep=";")



Run_model_importance_no_inter=function(id){
  
  list_mod=expand.grid(with_cover=c(F,T),
                       Stats=c("perim_area_scaling","fmax_psd","Cond_H","PL_expo","Spectral_ratio",
                               "core_area_land","division","fractal_dim","contig","core_area",
                               "flow_length","PLR",
                               "Struct1","Struct2")) #We also add the case where the stats correspond to the first two components of the PCA
  
  d_all=d_all2=tibble()#for keeping all information
  
  d_data=read.table("../Data/Spatial_structure_grazing.csv",sep=";")%>%
    Closer_to_normality(.)
  
  d_data[,c(1,9:29,36,38:42,45:ncol(d_data))] = 
    apply(d_data[,c(1,9:29,36,38:42,45:ncol(d_data))],2,z_tranform)
  
  d_data=Perform_PCA_spatial_struc(d_data)
  
  stat=list_mod$Stats[id]
  
  d_data_mod=d_data%>%melt(., measure.vars=stat)%>%
    filter(., !is.na(value))
  
  with_cover=list_mod$with_cover[id]
  
  if (with_cover){
    
    formula_mod=formula(formula_(paste("value ~ Long_cos + Long_sin + Lattitude + Slope + Elevation 
      + Clim1 + Clim2 + Clim3 + Clim4 + Grazing 
      + Woody + Woody*Grazing
      + Type_veg + rho_p       + Sand + Org_C_v       + (1|Site_ID)")))
    
  } else{
    
    formula_mod=formula(formula_(paste("value ~ Long_cos + Long_sin + Lattitude + Slope + Elevation 
      + Clim1 + Clim2 + Clim3 + Clim4 + Grazing 
      + Woody + Woody*Grazing
      + Type_veg  + Sand + Org_C_v
      + (1|Site_ID)")))
  }
  
  model_spa_stat  = lmer(formula_mod, d_data_mod,
                         na.action = na.fail ,REML ="FALSE")
  
  #we remove potential outliers
  rm.outliers = romr.fnc(model_spa_stat, d_data_mod, trim=2.5)
  d_data_out = rm.outliers$data
  
  
  
  model_spa_stat = lmer(formula_mod, data = d_data_out, 
                        na.action = na.fail,REML ="FALSE")
  
  # #do some model selection
  select_model=dredge(model_spa_stat, subset = ~ Slope & Elevation &
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
                      extra=c(R2m=function(x) r.squaredGLMM(x)[1],
                              R2c=function(x) r.squaredGLMM(x)[2]),
                      options(na.action = "na.fail") )
  
  
  
  #extract the result of model selection
  result_select=summary(model.avg(select_model, subset = delta < 2))
  
  
  #Get the importance of each metric
  importance_mod=sw(result_select)
  R2=select_model%>%filter(., AICc<min(AICc)+2)
  
  d_all=rbind(d_all,cbind(tibble(Sp_stat=stat,
                                 N_outliers=rm.outliers$n.removed,
                                 R2m=mean(R2$R2m),R2C=mean(R2$R2c)),
                          Aggregate_importance(importance_mod))%>%
                add_column(., With_cover=with_cover))
  
  write.table(d_all,paste0("../Data/Linear_models/Importance/Importance_",stat,"_",with_cover,"nointer.csv"),sep=";")
  
}

library(parallel)

mclapply(1:30,Run_model_importance_no_inter,mc.cores = 30)


d_all=d_all2=tibble()
for (k in list.files("../Data/Linear_models/Importance","nointer")){
  d_all=rbind(d_all,read.table(paste0("../Data/Linear_models/Importance/",k),sep=";"))
}
for (k in list.files("../Data/Linear_models/Estimator","nointer")){
  d_all2=rbind(d_all2,read.table(paste0("../Data/Linear_models/Estimator/",k),sep=";"))
}
write.table(d_all,"../Data/Linear_models/Importance_no_inter.csv",sep=";")
write.table(d_all2,"../Data/Linear_models/Estimators_model_no_inter.csv",sep=";")


