x = c("tidyverse", "ggpubr", "latex2exp", "reshape2", "simecol",
      "abc", "spatialwarnings", "FME","phaseR",
      "ggquiver", "scales","boot","RColorBrewer","ggnewscale","cluster","pls",
      "factoextra","FactoMineR","missMDA","GGally","diptest","raster","ape","abctools","viridis",
      "png","jpeg","landscapemetrics","lme4","lmeresampler","GGally","MuMIn",
      "LMERConvenienceFunctions","piecewiseSEM","qgraph")

#install pacakges if not installed already
install.packages(setdiff(x, rownames(installed.packages())))

#loading the packages
lapply(x, require, character.only = TRUE)



d_biocom_old=read.table("../Data/biocom_data.csv",sep=";")
d_Meta=read.table("../Data/Meta_data_sites.csv",sep=",",header = T)
d_biocom=readxl::read_xlsx("../Data/Final_biocom.xlsx")
d_biodesert=readxl::read_xlsx("../Data/Final_biodesert.xlsx")

the_theme=theme_classic()+theme(legend.position = "bottom",
                                strip.background = element_rect(fill=alpha("#CBB7D8",.6),color="black"),
                                strip.text.y = element_text(size = 10, angle = -90, face = "italic"),
                                strip.text.x = element_text(size = 10, face = "italic"),
                                legend.text = element_text(size = 10),text = element_text(family = "NewCenturySchoolbook"))
`%!in%` = Negate(`%in%`)

# 1) Getting matrices and plotting functions ----

rotate=function(x) t(apply(x, 2, rev))

Get_empirical_site_biocom=function(id){
  d_biocom=read.table("../Data/biocom_data.csv",sep=";")
  return(as.matrix(read.table(paste0("../../Linking_data_model_ABC/Data/Data_Biocom/landscapes/",d_biocom$File_ID[id],".txt"))))
}

Get_empirical_site_all=function(ERC="biocom",id,sub_id){
  return(as.matrix(read.table(paste0("../Data/Landscapes/Binary_landscapes/",
                                     ERC,"_50_",id,"_",sub_id,".csv"),sep=",")))
}

Get_png_empirical_site_biocom=function(id){
  d_biocom=read.table("../Data/biocom_data.csv",sep=";")
  img=readPNG(paste0("../../Linking_data_model_ABC/Data/Data_Biocom/png_img/",
                     ifelse(as.numeric(strsplit(d_biocom$File_ID[id],"-")[[1]][1])<100,
                            ifelse(as.numeric(strsplit(d_biocom$File_ID[id],"-")[[1]][1])<9,
                                   paste0("00",d_biocom$File_ID[id]),paste0("0",d_biocom$File_ID[id])),
                            d_biocom$File_ID[id]
                     ),
                     ".png"))
  return(img)
}

Get_png_empirical_site_biocom=function(ERC="biodesert",Size,ID,sub_ID){
  img=readJPEG(paste0("../Data/Landscapes/",ERC,"/",Size,"/",ID,"_",sub_ID,".jpeg"))
  return(img)
}


Plot_empirical=function(id,true_landscape=F){
  
  par(mfrow=c(1,1))
  
  if (true_landscape & is.character(id)){
    
    img1=readPNG(paste0("../../Linking_data_model_ABC/Data/Data_Biocom/png_img/",id))
    
  }else {
    if (is.numeric(id)){
      d_biocom=read.table("../Data/biocom_data.csv",sep=";")
      mat=as.matrix(read.table(paste0("../../Linking_data_model_ABC/Data/Data_Biocom/landscapes/",d_biocom$File_ID[id],".txt")))
      print(ggplot(mat%>%melt(.)) +
              geom_raster(aes(x = Var1, y = Var2,
                              fill = as.factor(value))) +
              coord_fixed() +
              theme_transparent() +theme(legend.position = "none")+
              scale_fill_manual(name = "Vegetation",
                                values = rev(c('black', 'white'))))
    }else {
      d_biocom=read.table("../Data/biocom_data.csv",sep=";")
      mat=as.matrix(read.table(paste0("../../Linking_data_model_ABC/Data/Data_Biocom/landscapes/",d_biocom$File_ID[which(d_biocom$File_ID==gsub(".png","",id))],".txt")))
      
      print(ggplot(mat%>%melt(.)) +
              geom_raster(aes(x = Var1, y = Var2,
                              fill = as.factor(value))) +
              coord_fixed() +
              theme_transparent() +theme(legend.position = "none")+
              scale_fill_manual(name = "Vegetation",
                                values = rev(c('black', 'white'))))
    }
  }
} 

Plot_lanscape=function(mat){
  
  if (mat[1,1] %in% c(F,T)){
    mat[mat==F]=0
    mat[mat==T]=1
  }
  
  if (length(unique(as.numeric(mat)))>2){
    print(ggplot(mat%>%melt(.)) +
      geom_raster(aes(x = Var1, y = Var2,
                      fill = (value))) +
      coord_fixed() +
      theme_transparent() +theme(legend.position = "none")+
      scale_fill_gradient2(low = "white",mid="gray",high="black"))
  } else {
    print(ggplot(mat%>%melt(.)) +
            geom_raster(aes(x = Var1, y = Var2,
                            fill = as.factor(value))) +
            coord_fixed() +
            theme_transparent() +theme(legend.position = "none")+
            scale_fill_manual(values=c("white","black")))

  }
}

myTryCatch=function(expr) {
  #' Catches errors and warnings in evaluation of expr
  #' 
  #' Particularly useful to detect optim problems in PL fit
  #' 
  #' @export
  
  warn=err=NULL
  
  value=withCallingHandlers(
    tryCatch(expr, error = function(e) {
      err <<- e
      NULL
    }), warning = function(w) {
      warn <<- w
      invokeRestart("muffleWarning")
    })
  
  list(value = value, warning = warn, error = err)
}


# 2) Image analysis ----
k_means_RGB = function(img, k) {
  
  imgDm = dim(img)
  
  # separate R,G,B
  
  imgRGB = data.frame(
    x = rep(1:imgDm[2], each = imgDm[1]),
    y = rep(imgDm[1]:1, imgDm[2]),
    R = as.vector(img[, , 1]),
    G = as.vector(img[, , 2]),
    B = as.vector(img[, , 3])
  )
  
  newBW = k_means(imgRGB, k, imgDm)
  
  return(newBW)
}

k_means = function(imgBands, k, dm) {
  
  k = as.integer(k)
  imgDm = dm
  kMeans = kmeans(imgBands[, c("R", "G", "B")], centers = k, nstart = 1)
  
  colors2 = c(1, 0)
  colors3 = c(1, 0.5, 0)
  colors4 = c(1, 0.7, 0.3, 0)
  
  colors = list(colors2, colors3, colors4)
  
  kMeansGray = kMeans$centers %>%
    as.data.frame() %>%
    mutate(num = 1:k, gray = 0.2126 * R + 0.7152 * G + 0.0722 * B) %>% # HDTV gray
    arrange(desc(gray))
  
  imgBands.update = cbind(imgBands, clust = kMeans$cluster) # %>% # add cluster info to img
  # ungroup
  
  imgBW = matrix(imgBands.update$clust,
                 nrow = imgDm[1],
                 ncol = imgDm[2]
  )
  
  # convert clusters numbers to gray scale
  
  newBW = imgBW # separate matrix
  for (l in 1:k) { # for each cluster category
    # replace clusters by ascending color (in gray : 1 -> 0)
    newBW[imgBW == kMeansGray$num[l]] = colors[[k - 1]][l]
  }
  
  return(newBW)
}

binarize = function(mat, cat0, cat1) {
  
  mat2 = mat
  mat2[mat %in% cat1] = 1 # full : black so gray is 0
  mat2[mat %in% cat0] = 0 # empty : white so gray is 1
  
  return(mat2)
}

get_cut_grayscale_values = function(nclust) {
  
  if (nclust == 2) {
    scale= c(1, 0)
  } else if (nclust == 3) {
    scale= c(1, 0.5, 0)
  } else if (nclust == 4) {
    scale=c(1, 0.7, 0.3, 0)
  }
  
  
  if (nclust == 2) {
    cut = list(scale)
  } else if (nclust == 3) {
    cut = list(
      list(scale[1], scale[2:3]),
      list(scale[1:2], scale[3])
    )
  } else if (nclust == 4) {
    cut = list(
      list(scale[1], scale[2:4]),
      list(scale[1:2], scale[3:4]),
      list(scale[1:3], scale[4])
    )
  }
  
  return(cut)
}


# 3) Spatial metrics functions ----

Get_sumstat=function(landscape,log_=T,slope=0){
  
  if (any(landscape==1 |landscape==T)){
    
    cover = sum(landscape) / (dim(landscape)[1]**2)
    
    # number of neighbors
    #vegetation clustering
    neighbors_mat = simecol::neighbors(x =landscape,state = 1, wdist =  matrix( c(0, 1, 0,1, 0, 1, 0, 1, 0), nrow = 3),bounds = 1)
    mean_nb_neigh = mean(neighbors_mat[which(landscape == 1)]) #mean number of plant neighbors
    mean_clustering = mean_nb_neigh / cover
    spatial_ews = generic_sews(landscape>0,4,moranI_coarse_grain = T)$value
    
    spectral_ratio = as.data.frame(spectral_sews(landscape>0,quiet=T))$value
    
    psd=spatialwarnings::patchdistr_sews(landscape>0)
    max_patchsize=max(psd$psd_obs)
    cv_patch=sd(psd$psd_obs)/mean(psd$psd_obs)
    PLR=spatialwarnings::raw_plrange(landscape>0)
    
  
    fit_psd=safe_psd_types(psd$psd_obs)
    
    if (log_){
      mean_clustering=log(mean_clustering)
      spectral_ratio=log(spectral_ratio)
      max_patchsize=log(max_patchsize/length(landscape))
    }
    
    
    
    #flow length
    flow_length=flowlength_sews(landscape>0,slope = slope,
                                cell_size = 50/sqrt(dim(landscape)[1]*dim(landscape)[2]))
  
    #power relationship between area and perimeter
    beta=lsm_c_pafrac(raster(landscape), directions = 4, verbose = TRUE)%>%filter(., class==1)%>%pull(., value)
    
    #mean perimeter-area ratio 
    mean_perim_area=lsm_c_para_mn(raster(landscape), directions = 4)%>%filter(., class==1)%>%pull(., value)
    
    #median, mean, sd of patch size (not in pixel unit but in m² to account for differences in resolutions)
    psd_scaled=lsm_p_area(raster(landscape),directions = 4)%>%filter(., class==1)%>%pull(., value)
    
    mean_psd=mean(psd_scaled)*1e4 #1e4 is to convert from ha to m²
    median_psd=median(psd_scaled)*1e4 #1e4 is to convert from ha to m²
    sd_psd=sd(psd_scaled)*1e4 #1e4 is to convert from ha to m²
  
    #core area metric
    core_area=lsm_c_cai_mn(raster(landscape),directions = 4)%>%filter(., class==1)%>%pull(., value)
    core_area_land=lsm_c_cpland(raster(landscape),directions = 4)%>%filter(., class==1)%>%pull(., value)
    
    #division of patches
    division=lsm_c_division(raster(landscape), directions = 4)%>%filter(., class==1)%>%pull(., value)
  
    #fractal dimension 
    fractal_dim=lsm_c_frac_mn(raster(landscape), directions = 4)%>%filter(., class==1)%>%pull(., value)
    
    #contig
    contig=lsm_c_contig_mn(raster(landscape),directions = 4)%>%filter(., class==1)%>%pull(., value)
    
    #All complexity measures at the landscape scale
    
    complex_land=calculate_lsm(raster(landscape), 
                               what = c("lsm_l_ent", "lsm_l_condent", "lsm_l_joinent","lsm_l_mutinf","lsm_l_relmutinf"),
                               full_name = TRUE,directions = 4,neighbourhood = 4)%>%
      dplyr::select(., value,name)
    
    Hx=lsm_l_ent(raster(landscape),neighbourhood = 4)%>%pull(., value)
    
    
    d=tibble(rho_p=cover,
             nb_neigh=mean_nb_neigh,clustering=mean_clustering,
             skewness=spatial_ews[2],variance=spatial_ews[1],moran_I=spatial_ews[3],
             Spectral_ratio=spectral_ratio,PLR=PLR,PL_expo=fit_psd["slope_best"],cv_psd=cv_patch,
             fmax_psd=max_patchsize,cutoff=fit_psd["cutoff_tpl"],slope_tpl=fit_psd["slope_tpl"],
             flow_length=flow_length$value, #flow length
             perim_area_scaling=beta, #scaling power relationship between area and perimeter 
             mean_perim_area=mean_perim_area, #mean perim/area 
             mean_psd=mean_psd, #mean patch size
             median_psd=median_psd, #median patch size
             sd_psd=sd_psd, #sd patch size
             core_area=core_area, #mean % of core area
             core_area_land=core_area_land,  #same but at the landscape scale
             division=division, #how much patches are divided or constitute big patches
             fractal_dim=fractal_dim, #fractal dimension
             contig=contig, #mean connectedness of cells in patches
             Cond_H=complex_land$value[1], #conditional entropy
             Shannon_H=complex_land$value[2], #shannon entropy
             Joint_H=complex_land$value[3], #joint entropy
             mutual_inf=complex_land$value[4], #mutual information
             relat_mutual_inf=complex_land$value[5] #relative mutual information
    )
  }else{
    d=tibble(rho_p=0,
             nb_neigh=0,
             clustering=0,
             skewness=0,variance=0,moran_I=0,
             Spectral_ratio=0,PLR=0,PL_expo=0,cv_psd=0,
             fmax_psd=0,cutoff=0,slope_tpl=0,
             flow_length=0,
             perim_area_scaling=0,
             mean_perim_area=0,
             mean_psd=0,
             median_psd=0,
             sd_psd=0,
             core_area=0,
             core_area_land=0,
             division=0,
             fractal_dim=0,
             contig=0,
             Cond_H=0,
             Shannon_H=0,
             Joint_H=0,
             mutual_inf=0,
             relat_mutual_inf=0
    )
    
    
  }
  
  return(d)
}

Get_spatial_resolution=function(landscape,n_meter=50){
  return(  Resolution = n_meter/sqrt(dim(landscape)[1]*dim(landscape)[2])
  )
}

safe_psd_types = function(psd) {
  
  if (length(unique(psd)) > 2) {
    
    # Trying fits for xmin = 1
    pl = safe_fit(pl_fit, psd, 1, F, c("plexpo"))
    tpl = safe_fit(tpl_fit, psd, 1, F, c("plexpo", "cutoff"))
    exp = safe_fit(exp_fit, psd, 1, F, c("cutoff"))
    lnorm = safe_fit(lnorm_fit, psd, 1, F, c("meanlog", "sdlog"))
    
    # detect warnings
    warnings = c(
      pl$warn,
      tpl$warn,
      exp$warn,
      lnorm$warn
    )
    
    bics = c(
      pl$bic,
      tpl$bic,
      exp$bic,
      lnorm$bic
    )
    
    comparison = data.frame(
      cbind(
        warn = warnings,
        bic = bics,
        type = 1:4
      )
    )
    # we want the type of the line without warn and with a normal bic
    best = comparison %>%
      filter(!warn & bic != -Inf) %>%
      filter(bic == min(bic)) %>%
      dplyr::select(type) %>%
      pull
    
    best_2 = comparison %>%
      slice_head(n = 3) %>% # exclude lnorm from comparison
      filter(!warn & bic != -Inf) %>%
      filter(bic == min(bic)) %>%
      dplyr::select(type) %>%
      pull
    
    # PSD shape parameters
    # we ensure that they are set to NA if they have warnings:
    # we need that every call returns vectors of the same lengths
    
    if (!tpl$warn) {
      slope_tpl  = tpl$plexpo
      cutoff_tpl  = tpl$cutoff
    } else {
      slope_tpl  = NA
      cutoff_tpl  = NA
    }
    
    if (!pl$warn) {
      slope_pl  = pl$plexpo
    } else {
      slope_pl  = NA
    }
    
    if (!is.na(slope_pl) & !is.na(slope_tpl)) {
      if (pl$bic < tpl$bic) {
        slope_best  = slope_pl
      } else {
        slope_best  = slope_tpl
      }
    } else {
      slope_best  = NA
    }
  } else {
    # if psd is too short, fitting anything does not make any sense:
    # return a NA vector
    best = NA
    best_2 = NA
    slope_pl  = NA
    slope_tpl  = NA
    cutoff_tpl = NA
    slope_best = NA
  }
  
  return(c(
    best = best,
    best_2 = best_2,
    slope_pl = slope_pl,
    slope_tpl = slope_tpl,
    cutoff_tpl = cutoff_tpl,
    slope_best = slope_best
  ))
}

safe_fit = function(fit_fun,psd,xmin = 1,bic_only = FALSE,force_output = NULL) {
  
  # define a compute bic function in internal scope
  # (needed for avoiding namespace issues when this is called by any function 
  # wrapped by spw::compute_indicator())
  
  compute_bic_local = function(fit, psd, xmin) {
    
    psd = psd[psd >= xmin] # fitting was done on x >= xmin
    
    return(fit$npars * log(length(psd)) - 2 * fit$ll)
  }
  
  try = myTryCatch(expr = {
    fit_fun(psd, xmin)
  })
  
  
  if (!length(try$error)) {
    fit = try$value
    warn = ifelse(length(try$warning) > 0, 1, 0) # did optim work well ?
    estimates = fit[names(fit) %in% c("plexpo", "cutoff", "meanlog", "sdlog")]
    
    if (any(sapply(estimates, is.nan))) { 
      # in some cases, some NaNs can be produced without warnings, catch them
      # estimates is a list and we want to keep it like that, thus the sapply()
      
      warn = 1
      estimates[sapply(estimates, is.nan)] = NA
    }
    
    bic = ifelse(
      warn,
      -Inf,
      compute_bic_local(fit, psd, xmin)
    )
  } else {
    warn = 1
    bic = -Inf
    estimates = NULL
  }
  
  out = list(
    bic = bic,
    warn = warn
  )
  
  if (!bic_only) {
    
    if (!length(estimates) & length(force_output)) {
      
      estimates = vector(mode = "list", length = length(force_output))
      names(estimates) = force_output
      
    }
    
    out = c(out, estimates)
  }
  
  return(out)
}


# 4) Statistical analysis functions ----

z_tranform=function(x){
  return((x-mean(x,na.rm=T))/sd(x,na.rm=T))
}



Compute_trend_data = function(d, spatial_stat,
                         both=T, # T = Grazing + Aridity
                         n_boot = 100, # number of boostrap
                         driver="Aridity" #if both = F, which driver ?
                         ){
  
  if (any(is.na(d[,spatial_stat]))){d=d[-which(is.na(d[,spatial_stat])),]} #removing the NA's
  
  if (both){ #both aridity and grazing in the model
    
    #full model
    model=d%>% lmer(formula = paste(spatial_stat, "~ Aridity + Grazing + ( 1 | Site_ID)"), data = .)
    
    #only random effect
    null=d%>% lmer(formula = paste(spatial_stat, "~ ( 1 | Site_ID)"), data = .)
    
    #null model for grazing (i.e. only aridity)
    null_grazing=d%>% lmer(formula = paste(spatial_stat, "~ Aridity + ( 1 | Site_ID)"), data = .)
    
    #null model for aridity (i.e. only grazing)
    null_aridity=d%>% lmer(formula = paste(spatial_stat, "~ Grazing + ( 1 | Site_ID)"), data = .)
    
    n_driver=2;drivers=c("Aridity","Grazing")  
    
  } else {
    
    model=d%>% lmer(formula = paste(spatial_stat, "~ ",driver," + ( 1 | Site_ID)"), data = .)
    null=d%>% lmer(formula = paste(spatial_stat, "~ ( 1 | Site_ID)"), data = .)
    n_driver=1
    drivers=driver
  } 
  
  #Type II test for model (compared to only random site effect)
  typeII_test= as_tibble(anova(model,null)) 
  p_val_II = typeII_test$`Pr(>Chisq)`[2]
  
  #boostraping (n = n_boot)
  lme_boot=bootstrap(model, .f = fixef, type = "parametric", B = n_boot)
  
  #getting slopes median and quantiles
  median_stat = apply(lme_boot$replicates[,2:(n_driver+1)],2,median)
  quantiles_stat = apply(lme_boot$replicates[,2:(n_driver+1)],2,quantile,c(.025,.975))

  #if joint model, testing effect of aridity and grazing in isolation + aggregating in df  
  if (both){
    
    typeII_grazing= as_tibble(anova(model,null_grazing)) 
    p_val_g = typeII_grazing$`Pr(>Chisq)`[2]

    typeII_aridity= as_tibble(anova(model,null_aridity)) 
    p_val_a = typeII_aridity$`Pr(>Chisq)`[2]
    
    d2=tibble(Stat=spatial_stat,
              p_val_full=p_val_II, #pvalue associated to the model test
              p_val_indiv=c(p_val_a,p_val_g), #pvalue associated to each predictor effect
              q2=median_stat, #median
              q1=quantiles_stat[1,], #quantile 2.5%
              q3=quantiles_stat[2,], #quantile 97.5%
              Driver=drivers
    )
    
  }else{
    d2=tibble(Stat=spatial_stat,
              p_val_full=p_val_II,#pvalue associated to the model test
              q2=median_stat, #median
              q1=quantiles_stat[1,], #quantile 2.5%
              q3=quantiles_stat[2,], #quantile 97.5%
              Driver=drivers
    )
  }
  
  
  return(d2)
  
}




Test_grazing_effect = function(d, spatial_stat){
  
  if (any(is.na(d[,spatial_stat]))){d=d[-which(is.na(d[,spatial_stat])),]} #removing the NA's
  
  
  #full model
  model=d%>% lmer(formula = paste(spatial_stat, "~ Aridity + Grazing + ( 1 | Site_ID)"), data = .)
  
  #null model for grazing (i.e. only aridity)
  null_grazing=d%>% lmer(formula = paste(spatial_stat, "~ Aridity + ( 1 | Site_ID)"), data = .)
  

  #Type II test for the effect of grazing
  typeII_test= as_tibble(anova(model,null_grazing)) 
  p_val_II = typeII_test$`Pr(>Chisq)`[2]
  
  d2=tibble(Stat=spatial_stat,
            p_val_full=p_val_II #pvalue associated to the model test
  )
  
  
  return(d2)
  
}

Test_interaction_effect = function(d, spatial_stat){
  
  if (any(is.na(d[,spatial_stat]))){d=d[-which(is.na(d[,spatial_stat])),]} #removing the NA's
  
  
  #full model
  model=d%>% lmer(formula = paste(spatial_stat, "~ Aridity + Grazing + Grazing*Aridity + ( 1 | Site_ID)"), data = .)
  
  #null model for grazing (i.e. only aridity)
  null_grazing=d%>% lmer(formula = paste(spatial_stat, "~ Aridity + Grazing + ( 1 | Site_ID)"), data = .)
  
  
  #Type II test for the effect of grazing
  typeII_test= as_tibble(anova(model,null_grazing)) 
  p_val_II = typeII_test$`Pr(>Chisq)`[2]
  
  d2=tibble(Stat=spatial_stat,
            p_val_full=p_val_II #pvalue associated to the model test
  )
  
  
  return(d2)
  
}

Compute_trend_model = function(d, spatial_stat,treshold=.05){ #below treshold, we estimate that there is no more cover
  
  if (any(is.na(d[,spatial_stat]))){d=d[-which(is.na(d[,spatial_stat])),]} #removing the NA's
  
  d_graz=d%>%filter(., Driver=="Grazing",rho_p>.05)
  d_arid=d%>%filter(., Driver=="Aridity",rho_p>.05)%>%mutate(., b=1-b) #to make similar interpretation as aridity
  
  d_arid[,spatial_stat]=z_tranform(d_arid[,spatial_stat])
  d_graz[,spatial_stat]=z_tranform(d_graz[,spatial_stat])
  
  
  #model grazing (g0)
  model_grazing=d_graz%>% lm(formula = paste(spatial_stat, "~ g0 "), data = .)
  
  #model aridity (b)
  model_aridity=d_arid%>% lm(formula = paste(spatial_stat, "~ b"), data = .)

  d2=tibble(Stat=spatial_stat,
            p_val_indiv=c(0,0), #no p val in models
            q2=c(model_aridity$coefficients[2],model_grazing$coefficients[2]), 
            q1=c(model_aridity$coefficients[2],model_grazing$coefficients[2]), #as no bootstrap, we set the same
            q3=c(model_aridity$coefficients[2],model_grazing$coefficients[2]), #as no bootstrap, we set the same 
            Driver=c("Aridity","Grazing")
  )
  
  return(d2)
  
}

Aggregate_importance=function(importance_df){
  
  importance_df=t(as.data.frame(importance_df))
  
  if (any(grep(":Grazing",colnames(importance_df)))){
    if (any(grep("Grazing:",colnames(importance_df)))){
      Interactions=mean(importance_df[,c(grep("Grazing:",colnames(importance_df)),
                                         grep(":Grazing",colnames(importance_df)))])
    }else {
      Interactions=mean(importance_df[,grep(":Grazing",colnames(importance_df))])
    }
  }else {
    if (any(grep("Grazing:",colnames(importance_df)))){
      Interactions=mean(importance_df[,grep("Grazing:",colnames(importance_df))])
    }
  }
  
  d=tibble(Interactions=ifelse(exists("Interactions"),Interactions,0),
         Sand =ifelse("Sand" %in% colnames(importance_df),importance_df[,"Sand"],0),
         Clim1=ifelse("Clim1" %in% colnames(importance_df),importance_df[,"Clim1"],0),
         Clim2=ifelse("Clim2" %in% colnames(importance_df),importance_df[,"Clim2"],0),
         Clim3=ifelse("Clim3" %in% colnames(importance_df),importance_df[,"Clim3"],0),
         Clim4=ifelse("Clim4" %in% colnames(importance_df),importance_df[,"Clim4"],0),
         Org_C=ifelse("Org_C" %in% colnames(importance_df),importance_df[,"Org_C"],0),
         Type_veg=ifelse("Type_veg" %in% colnames(importance_df),importance_df[,"Type_veg"],0),
         Grazing=ifelse("Grazing" %in% colnames(importance_df),importance_df[,"Grazing"],0),
         Cover=ifelse("rho_p" %in% colnames(importance_df),importance_df[,"rho_p"],0),
         Woody=ifelse("Woody" %in% colnames(importance_df),importance_df[,"Woody"],0),
  )
  
  return(d)  
  
}

Organize_df=function(df,type="predictor"){
  
  if (type=="predictor"){ #i.e. plots with standardize effects 
    
      df=df%>%
        mutate(., Stat=recode_factor(Stat,"fmax"="Max patch","perim_area_scaling"="Perim-area scaling",
                                 "PL"="PSD exponent","core_area"="Core-area" ))%>%
      mutate(.,Type_pred=unlist(sapply(1:nrow(.),function(x){
        if (.$term[x] %in% c("rho_p","Type_vegGrass_Shrub","Type_vegGrassland","Type_vegShrubland","Woody")) return("Vegetation")
        if (.$term[x] %in% c("Clim1","Clim2","Clim3","Clim4")) return("Climatic")
        if (.$term[x] %in% c("Org_C","Sand","Slope")) return("Abiotic")
        if (.$term[x] %in% c("Grazing:rho_p","Clim1:Grazing","Clim2:Grazing",
                             "Clim3:Grazing","Clim4:Grazing","Grazing:Type_vegShrubland",
                             "Grazing:Type_vegGrassland","Grazing:Type_vegGrass_Shrub",
                             "Grazing:Org_C","Grazing:Sand","Grazing","Grazing:Woody")) return("Grazing")
        if (.$term[x] %in% c("Long_cos","Long_sin","Lattitude","Elevation")) return("Geo")
      })))%>%
      mutate(., term=recode_factor(term,
                                   "rho_p"="Cover",
                                   "Type_vegGrassland"="Type, grassland",
                                   "Type_vegShrubland"="Type, shrubland",
                                   "Type_vegGrass_Shrub"="Type, both",
                                   "Org_C"="Facilitation",
                                   "Long_sin"="Longitude (sin)",
                                   "Long_cos"="Longitude (cos)",
                                   "Grazing:Type_vegShrubland"="Grazing * Type, shrubland",
                                   "Grazing:Type_vegGrassland"="Grazing * Type, grassland",
                                   "Grazing:Type_vegGrass_Shrub"="Grazing * Type, both",
                                   "Grazing:Org_C"="Grazing * Facilitation",
                                   "Clim1:Grazing"="Grazing * Clim1",
                                   "Clim2:Grazing"="Grazing * Clim2",
                                   "Clim3:Grazing"="Grazing * Clim3",
                                   "Clim4:Grazing"="Grazing * Clim4",
                                   "Grazing:Sand"="Grazing * Sand",
                                   "Grazing:Woody"="Grazing * Woody",
                                   "Grazing:rho_p"="Grazing * Cover"))%>%
        arrange(., Type_pred,observed)%>%
        add_column(., Order_f=1:nrow(.))%>%
        mutate(.,term = fct_reorder(term, Order_f))
    
  } else{
   
    df=df%>%
      mutate(.,Type_pred=unlist(sapply(1:nrow(.),function(x){
        if (.$term[x] %in% c("rho_p","Type_vegGrass_Shrub","Type_vegGrassland","Type_vegShrubland","Woody")) return("Vegetation")
        if (.$term[x] %in% c("Clim1","Clim2","Clim3","Clim4")) return("Climatic")
        if (.$term[x] %in% c("Org_C","Sand","Slope")) return("Abiotic")
        if (.$term[x] %in% c("Grazing:rho_p","Clim1:Grazing","Clim2:Grazing",
                             "Clim3:Grazing","Clim4:Grazing","Grazing:Type_vegShrubland",
                             "Grazing:Type_vegGrassland","Grazing:Type_vegGrass_Shrub",
                             "Grazing:Org_C","Grazing:Sand","Grazing","Grazing:Woody")) return("Grazing")
        if (.$term[x] %in% c("Long_cos","Long_sin","Lattitude","Elevation")) return("Geo")
      })))
    
  }
  return(df)
}

Plot_SEM_without_cover=function(summary_sem,pdf_=F,name="SEM",title_="",name_var=""){
  
  l=summary_sem$coefficients[,c("Predictor","Response","Estimate")]
  l$colo="#9ED471"
  l$colo[which(l$Estimate<0)]="#EFC46A"
  library(igraph)
  g <- graph.data.frame(l, directed=T)
  g= g %>% set_edge_attr("color", value =l$colo)
  g= g %>% set_vertex_attr("name", value =c("Grazing","Woody","Sand","Facilitation",
                                            name_var))
  coord=data.frame(label=vertex_attr(g, "name"),
                   lab2=c("Grazing","Woody","Sand","Facilitation",
                          name_var),
                   x=c(-10,0,0,0,10),y=c(0,-30,10,20,0))
  EL=as_edgelist(g)
  EL=cbind(EL,l[,3])
  name_edge=sapply(1:length(summary_sem$coefficients$P.Value),function(x){#adding pvalue
    return(paste0(round(l[x,3],3),ifelse(summary_sem$coefficients$P.Value[x]<.1,
                                         ifelse(summary_sem$coefficients$P.Value[x]<.05,
                                                ifelse(summary_sem$coefficients$P.Value[x]<.01,"***","**"),"*"),"")))
  })
  asi=abs(l[,3])/0.05
  asi[asi<5]=5
  
  if(pdf_){
    pdf(paste0("./",name,".pdf"),width = 7,height = 4)
    qgraph(EL,layout=as.matrix(coord[,c("x","y")]),edge.color=l$colo,edge.labels=name_edge,
           border.color="white",label.cex=1.5,label.scale=F,title=title_,
           edge.label.cex = 1.5,edge.label.position=0.5,vsize2=4,vsize=30,
           shape="ellipse",edge.labels=T,fade=F,esize=5,asize=asi,
           mar=rep(3,4))
    dev.off()
    
  }else {
    qgraph(EL,layout=as.matrix(coord[,c("x","y")]),edge.color=l$colo,edge.labels=name_edge,
           border.color="white",label.cex=1.5,label.scale=F,
           edge.label.cex = 1.5,edge.label.position=0.5,vsize2=12,vsize=25,
           shape="ellipse",edge.labels=T,fade=F,esize=5,asize=asi,
           mar=rep(3,4))
    
  }
}

Plot_SEM_with_cover=function(summary_sem,pdf_=F,name="SEM",title_="",name_var=""){
  
  
  l=summary_sem$coefficients[,c("Predictor","Response","Estimate")]
  l$colo="#9ED471"
  l$colo[which(l$Estimate<0)]="#EFC46A"
  library(igraph)
  g <- graph.data.frame(l, directed=T)
  g= g %>% set_edge_attr("color", value =l$colo)
  g= g %>% set_vertex_attr("name", value =c("Grazing","Woody","Sand","Facilitation","Cover",
                                            name_var))
  coord=data.frame(label=vertex_attr(g, "name"),
                   lab2=c("Grazing","Woody","Sand","Facilitation","Cover",
                          name_var),
                   x=c(-10,0,0,0,0,10),y=c(0,-30,-10,10,30,0))
  EL=as_edgelist(g)
  EL=cbind(EL,l[,3])
  name_edge=sapply(1:length(summary_sem$coefficients$P.Value),function(x){#adding pvalue
    return(paste0(round(l[x,3],3),ifelse(summary_sem$coefficients$P.Value[x]<.1,
                                         ifelse(summary_sem$coefficients$P.Value[x]<.05,
                                                ifelse(summary_sem$coefficients$P.Value[x]<.01,"***","**"),"*"),"")))
  })
  asi=abs(l[,3])/0.05
  asi[asi<5]=5
  
  if(pdf_){
    pdf(paste0("./",name,".pdf"),width = 7,height = 4)
    qgraph(EL,layout=as.matrix(coord[,c("x","y")]),edge.color=l$colo,edge.labels=name_edge,
           border.color="white",label.cex=1.5,label.scale=F,title=title_,
           edge.label.cex = 1.5,edge.label.position=0.5,vsize2=4,vsize=30,
           shape="ellipse",edge.labels=T,fade=F,esize=5,asize=asi,
           mar=rep(3,4))
    dev.off()
    
  }else {
    qgraph(EL,layout=as.matrix(coord[,c("x","y")]),edge.color=l$colo,edge.labels=name_edge,
           border.color="white",label.cex=1.5,label.scale=F,
           edge.label.cex = 1.5,edge.label.position=0.5,vsize2=12,vsize=25,
           shape="ellipse",edge.labels=T,fade=F,esize=5,asize=asi,
           mar=rep(3,4))
    
  }
}



# 5) Model functions ----



Get_params_model=function(delta,b,c,m,d,r,f,g0,z){
  
  return(list(
    
    delta  =  delta ,
    b      =  b     ,
    c      =  c     ,
    m      =  m     ,
    d      =  d     ,
    r      =  r     ,
    f      =  f     ,
    z      =  z     ,
    g0     =  g0
  ))
}

Get_classical_param=function(delta= 0.2,b=0.6,c= 0.3,m=.05,d= 0.2,r=0.0001,f=0.9,z=4,g0=0.2){
  return(Get_params_model(     delta= delta,b=b,c= c,m=m,d= d,r=r,f=f,z=z,g0=g0))
}

Get_initial_lattice=function(frac=c(.1,.1,.8),size=100){
  
  # Return a matrix of states with different initial fraction of vegetation, fertile and desert sites
  
  return(matrix(sample(-1:1,replace = T,size=size*size,prob = frac),ncol=size,nrow=size))
}


fourneighbors = function(landscape, state = 1, bounds = 1) {
  
  #Compute the number of neighbors in a given state for all the landscape 
  
  neighborhood = matrix(
    c(0, 1, 0,
      1, 0, 1,
      0, 1, 0),
    nrow = 3)
  nb_neighbors = simecol::neighbors(
    x = landscape,
    state = state,
    wdist = neighborhood,
    bounds = bounds)
  
  return(nb_neighbors)
  
}

#for (k in 1:length(param))assign(names(param)[k],param[[k]])

CA_spatial_model = function(init, params) {
  
  # Variables : 1 = vegetation, 0 = fertile, -1 = degraded   
  landscape = init
  rho_v = sum(landscape == 1) / length(landscape)
  
  # Neighbors :
  neigh_v = fourneighbors(landscape, state = 1, bounds = 1)

  
  with(params, {
    
    delta_t=.5 #more precision on the spatial model
    
    colonization = (delta * rho_v + (1 - delta) * neigh_v / z) *(b - c * rho_v )*delta_t
    
    # calculate regeneration, degradation & mortality rate
    death = (m + g0*(1 - neigh_v/z))*delta_t
    regeneration = (r + f * neigh_v / z)*delta_t
    degradation = d*delta_t 
    
    # Apply rules
    rnum = runif(length(landscape)) # one random number between 0 and 1 for each cell
    landscape_update = landscape
    
    ## New vegetation
    landscape_update[which(landscape == 0 & rnum <= colonization)] = 1
    
    ## New fertile
    landscape_update[which(landscape == 1 & rnum <= death)] = 0
    landscape_update[which(landscape == -1 & rnum <= regeneration)] = 0
    
    ## New degraded 
    landscape_update[which(landscape == 0 & rnum > colonization & rnum <= colonization + degradation)] = -1
    
    
    return(list(Rho_v = sum(landscape_update == 1) / length(landscape_update),
                Rho_f = sum(landscape_update == 0) / length(landscape_update),
                Rho_D = 1-(sum(landscape_update == 1) / length(landscape_update))-(sum(landscape_update == 0) / length(landscape_update)),
                Landscape=landscape_update))
  })
  
}


Run_spatial_model=function(time=seq(1,2000,1),params,ini,plot=F){
  
  
  d=tibble(Time=1,
           Rho_V=sum(ini == 1) / length(ini),
           Rho_F=sum(ini == 0) / length(ini),
           Rho_D=sum(ini == -1) / length(ini))
  state=list(Landscape=ini,Rho_v=d$Rho_V,Rho_f=d$Rho_F,Rho_D=d$Rho_D)
  
  for (k in 2:length(time)){
    
    params$dt=time[k]-time[k-1]
    
    state=CA_spatial_model(state$Landscape,params = params)
    d=rbind(d,tibble(Time=k,Rho_V=state$Rho_v,Rho_F=state$Rho_f,Rho_D=state$Rho_D))
  }
  
  return(list(d=d,State=state$Landscape)) #record time evolution and last landscape
  
}

Plot_landscape=function(landscape){
  landscape[landscape<1]=0
  image(landscape,xaxt = "n",yaxt ="n",col=rev(c("black","white")) )
}
