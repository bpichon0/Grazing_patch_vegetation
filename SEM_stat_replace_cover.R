rm(list=ls())
source("./Structure_grazing_function.R")

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
