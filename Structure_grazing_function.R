x = c("tidyverse", "ggpubr", "latex2exp", "reshape2", "simecol","ggeffects","nlme",
      "abc", "spatialwarnings", "FME","phaseR","xgboost","igraph","MultiVarSel",
      "ggquiver", "scales","boot","RColorBrewer","ggnewscale","cluster","vegan",
      "factoextra","FactoMineR","missMDA","GGally","diptest","raster","ape","abctools","viridis",
      "png","jpeg","landscapemetrics","lme4","lmeresampler","lmerTest","GGally","MuMIn","multcompView",
      "LMERConvenienceFunctions","semEff","piecewiseSEM","qgraph","car","spdep","ParBayesianOptimization",
      "ggpattern","ggpmisc","ggpp","gradientForest","extendedForest","rfPermute","A3",
      "psych","emmeans","vegan","adiv","ade4","FD","ggh4x","deSolve","plsRglm","caret")

# source("https://gist.githubusercontent.com/benmarwick/2a1bb0133ff568cbe28d/raw/fb53bd97121f7f9ce947837ef1a4c65a73bffb3f/geom_flat_violin.R")

#loading the packages
lapply(x, require, character.only = TRUE)
#If not the package: install them by uncommenting the following line
#lapply(x, install.packages, character.only = TRUE)


the_theme=theme_classic()+theme(legend.position = "bottom",
                                strip.background = element_rect(fill="grey",color="black"),
                                strip.text.y = element_text(size = 10, angle = -90, face = "italic"),
                                strip.text.x = element_text(size = 10, face = "italic"),
                                legend.text = element_text(size = 10),text = element_text(family = "NewCenturySchoolbook"))


the_theme2 = theme_classic() + theme(
  legend.position = "bottom",
  panel.border = element_rect(colour = "black", fill=NA),
  strip.background = element_rect(fill = "transparent",color="transparent"),
  strip.text.y = element_text(size = 10, angle = -90),
  strip.text.x = element_text(size = 10),title = element_text(size=8),
  axis.title.y=element_text(size = 10),
  axis.title.x=element_text(size = 10),
  #legend.box="vertical",
  legend.text = element_text(size = 10), text = element_text(family = "NewCenturySchoolbook")
)


# scale_color_manual(values=c("black","green","yellow","red"))+
# cor(d_data[,"contig"],d_data[,c("rho_p","Org_C","Sand","Clim1","Clim2","Clim3","Clim4","Woody","Longitude","Slope","Elevation","Long_cos","Long_sin","Grazing")])

dir.create("./Figures/",showWarnings = F)
dir.create("./Figures/SI",showWarnings = F)

scale_values=function(x){return((x-mean(x,na.rm=T))/sd(x, na.rm = T))}

`%!in%` = Negate(`%in%`)

is_signif=function(x){
  if (x<=.001){
    return("*")
  }else if (x>.001 & x<=.01){
    return("*")
  }else if (x>.01 & x<.05){
    return("*")
  # }else if (x>.05 & x<.1){
  #   return("°")
  }else {
    return("")
  }
}

formula_=function(x){
  return(gsub("\n","",x))
}

twoside_pvalue=function(data){
  p1=sum(data>0)/length(data)
  p2=sum(data<0)/length(data)
  p=min(p1,p2)*2
  return(p)
}

isEmptyNumeric = function(x) {
  return(identical(x, numeric(0)))
}

get_bootstrapped_pval=function(x){
  return(ifelse(length(which(x>0))/length(x)>.5,length(which(x<0))/length(x),length(which(x>0))/length(x)))
}

Add_temperature_precip=function(d){
  return(d%>%
           add_column(., 
                      MAT=sapply(1:nrow(.),function(x){return(d_biodesert$AMT[which(d_biodesert$ID==.$ID[x])])}),
                      MAP=sapply(1:nrow(.),function(x){return(d_biodesert$RAI[which(d_biodesert$ID==.$ID[x])])}),
                      Max_Temp=sapply(1:nrow(.),function(x){return(d_biodesert$MAWM[which(d_biodesert$ID==.$ID[x])])}),
                      Precip_driest=sapply(1:nrow(.),function(x){return(d_biodesert$RADM[which(d_biodesert$ID==.$ID[x])])})))
}

# 1) Getting matrices and plotting functions ----

rotate=function(x) t(apply(x, 2, rev))

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
            theme(legend.position = "none")+
            scale_fill_gradient2(low = "white",mid="gray",high="black"))
  } else {
    print(ggplot(mat%>%melt(.)) +
            geom_raster(aes(x = Var1, y = Var2,
                            fill = as.factor(value))) +
            coord_fixed() +
            theme(legend.position = "none")+
            scale_fill_manual(values=c("white","black")))
    
  }
}

Plot_psd_raw=function(landscape){
  psd_id=spatialwarnings::patchsizes(landscape>0)
  psd_tibble=tibble(patch_size=psd_id)
  
  psd_tibble$freq=sapply(1:nrow(psd_tibble),function(x){
    return(length(which(psd_tibble$patch_size>=psd_tibble$patch_size[x])))
  })
  return(ggplot(psd_tibble)+
          geom_point(aes(x=patch_size,y=freq))+
          the_theme2+
          labs(x="Patch size",y="Number of patches")+
          scale_x_log10(limits=c(1,1000))+
          scale_y_log10(limits=c(1,1000)))
  
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
      err <= e
      NULL
    }), warning = function(w) {
      warn <= w
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

Get_KS_distance = function(mat,n_shuffle) {
  obs_psd = spatialwarnings::patchsizes(mat)
  # Compute KS distance with null N_SHUFFLE times
  all_ks = unlist(lapply(seq.int(n_shuffle), function(i) {
    null_mat = matrix(sample(mat), nrow = nrow(mat), ncol = ncol(mat))
    null_psd = spatialwarnings::patchsizes(null_mat)
    ks.test(obs_psd, null_psd)[["statistic"]]
  }))
  
  return( mean(all_ks) )
}

Get_sumstat=function(landscape,log_=T,slope=0,compute_KS=T){

  
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
    
    Small_patches=table(patchsizes(landscape>0))[1]
    Small_patches_fraction=table(patchsizes(landscape>0))[1]/sum(table(patchsizes(landscape>0)))
    
    fit_psd=safe_psd_types(psd$psd_obs)
    
    if (log_){
      mean_clustering=log(mean_clustering)
      spectral_ratio=log(spectral_ratio)
      max_patchsize=log(max_patchsize/length(landscape))
    }
    
    if (compute_KS){
      if(length(psd$psd_obs>0)){
        ks_dist  = Get_KS_distance(landscape>0,n_shuffle = 199)
      }else{
        ks_dist  = NA
      }
    }else{
      ks_dist=NA
    }

    #flow length
    flow_length=flowlength_sews(landscape>0,slope = slope,
                                cell_size = 50/sqrt(dim(landscape)[1]*dim(landscape)[2]))
    
    #power relationship between area and perimeter
    beta=lsm_c_pafrac(raster(landscape), directions = 8, verbose = TRUE)%>%dplyr::filter(., class==1)%>%pull(., value)
    
    #mean perimeter-area ratio 
    mean_perim_area=lsm_c_para_mn(raster(landscape), directions = 8)%>%dplyr::filter(., class==1)%>%pull(., value)
    
    #median, mean, sd of patch size (not in pixel unit but in m? to account for differences in resolutions)
    psd_scaled=lsm_p_area(raster(landscape),directions = 8)%>%dplyr::filter(., class==1)%>%pull(., value)
    
    mean_psd=mean(psd_scaled)*1e4 #1e4 is to convert from ha to m?
    median_psd=median(psd_scaled)*1e4 #1e4 is to convert from ha to m?
    sd_psd=sd(psd_scaled)*1e4 #1e4 is to convert from ha to m?

    #core area metric
    core_area=lsm_c_cai_mn(raster(landscape),directions = 8)%>%dplyr::filter(., class==1)%>%pull(., value)
    core_area_land=lsm_c_cpland(raster(landscape),directions = 8)%>%dplyr::filter(., class==1)%>%pull(., value)
    
    #division of patches
    division=lsm_c_division(raster(landscape), directions = 8)%>%dplyr::filter(., class==1)%>%pull(., value)
    
    #fractal dimension 
    fractal_dim=lsm_c_frac_mn(raster(landscape), directions = 8)%>%dplyr::filter(., class==1)%>%pull(., value)
    
    #contig
    contig=lsm_c_contig_mn(raster(landscape),directions = 8)%>%dplyr::filter(., class==1)%>%pull(., value)
    
    #shape
    Shape_metric=lsm_c_shape_mn(raster(landscape),directions = 8)%>%dplyr::filter(., class==1)%>%pull(., value)
    
    #All complexity measures at the landscape scale
    
    complex_land=calculate_lsm(raster(landscape), 
                               what = c("lsm_l_ent", "lsm_l_condent", "lsm_l_joinent","lsm_l_mutinf","lsm_l_relmutinf"),
                               full_name = TRUE,directions = 8,neighbourhood = 4)%>%
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
             Small_patches=Small_patches, #number of smaller patches
             Small_patches_fraction=Small_patches_fraction, #number of smaller patches
             sd_psd=sd_psd, #sd patch size
             KS_dist=ks_dist, #Kolmogorov distance with null expectation
             core_area=core_area, #mean % of core area
             core_area_land=core_area_land,  #same but at the landscape scale
             division=division, #how much patches are divided or constitute big patches
             fractal_dim=fractal_dim, #fractal dimension
             contig=contig, #mean connectedness of cells in patches
             Shape_metric=Shape_metric,#shape of vegetation patches
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
             Small_patches=0,
             Small_patches_fraction=0, 
             sd_psd=0,
             KS_dist=0,
             core_area=0,
             core_area_land=0,
             division=0,
             fractal_dim=0,
             contig=0,
             Shape_metric=0,
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
      dplyr::filter(!warn & bic != -Inf) %>%
      dplyr::filter(bic == min(bic)) %>%
      dplyr::select(type) %>%
      pull
    
    best_2 = comparison %>%
      slice_head(n = 3) %>% # exclude lnorm from comparison
      dplyr::filter(!warn & bic != -Inf) %>%
      dplyr::filter(bic == min(bic)) %>%
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

Closer_to_normality=function(df,controlled_by_cover=F){
  
  if (controlled_by_cover){
    df=df%>%
      mutate(., 
             Org_C_v=log(Org_C_v), 
             Slope=log(Slope+1))
    
  }else{
    df=df%>%
      mutate(., 
             flow_length=log(flow_length),
             cv_psd=log(cv_psd),
             mean_psd=log(mean_psd),
             sd_psd=log(sd_psd),
             median_psd=log(median_psd),
             Small_patches=log(Small_patches),
             core_area=log(core_area),
             Shape_metric=log(Shape_metric),
             Org_C_v=log(Org_C_v), 
             Slope=log(Slope+1))
    
  }
  return(df)
}




z_tranform=function(x){
  x[is.infinite(x)]=NA
  return((x-mean(x,na.rm=T))/sd(x,na.rm=T))
}



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

Closer_to_normality_CWM=function(df){
  
  
  return(df%>%
           mutate(., 
                  CWM_MaxLS=log(CWM_MaxLS+1),
           ))
}


Perform_PCA_spatial_struc=function(df,plot=F){
  save=df
  df=df%>%
    dplyr::rename("% landscape \n covered by core"="core_area_land",
                  "Mean % of \n core in patches"="core_area",
                  "Division"="division",
                  "Bare soil \n connectivity"="flow_length",
                  "Fractal dimension"="fractal_dim",
                  "Number of \n smallest patches"="Small_patches",
                  "Spatial autocorrelation \n of vege."="moran_I",
                  "Fractal scaling \n area, perim."="perim_area_scaling",
                  "PLR"="PLR",
                  "Mean patch \n size"="mean_psd",
                  "Distance to \n null expect."="KS_dist",
                  "Power-law exp. \n of the PSD"="PL_expo",
                  "Spectral ratio"="Spectral_ratio",
                  "log (largest patch)"="fmax_psd",
                  "Patch shape complexity"="Shape_metric")
  
  
  #Performing PCA on the spatial structure metrics 
  struct_variables=colnames(df)[c(6,9,11,14,17,19,22)]
  res.comp=imputePCA(df[,which(colnames(df) %in% struct_variables)],ncp=4,scale = T) 
  
  if ("completeObs" %in% names(res.comp)){
    res.pca=PCA(res.comp$completeObs, ncp = 4,  graph=F)
  }else {
    res.pca=PCA(res.comp, ncp = 4,  graph=F)
  }
  
  #ploting the PCA
  if (plot){
    print(factoextra::fviz_pca_var(res.pca,col.var="#2746B1",
                                   label.var=c("Spatial autocorr.","Power-law exp. \n of the PSD",
                                               "log (largest patch)","Bare soil \n connectivity",
                                               "Mean patch \n size","Number of the \n smallest patches",
                                               "% landscape covered \n by core pixels"))+
            the_theme+
            ggtitle("")+
            labs(x=paste0("PC 1 (",round(res.pca$eig[1,2],1),")"),
                 y=paste0("PC 2 (",round(res.pca$eig[2,2],1),")")))
  }
  
  #We extract the first 2 ones (2/3 of the variance) and add them to the data-frame
  save=save%>%
    add_column(.,Struct1=res.pca$ind$coord[,1],Struct2=res.pca$ind$coord[,2])
  return(save)
  
}


Filter_relevant_stats=function(df){
  return(df%>%dplyr::filter(., Stat %in% c("PL_expo","fmax_psd","flow_length",
                                    "moran_I","core_area","Small_patches",
                                    "mean_psd")))
}

Rename_spatial_statistics=function(df){
  return(df%>%
           mutate(., Stat=recode_factor(Stat,
                                        "core_area"="Circularity of patches",
                                        "mean_psd"="Mean patch size",
                                        "flow_length"="Bare soil connectivity",
                                        "PLR"="PLR",
                                        "Small_patches"="Number of smallest patches",
                                        "moran_I"="Spatial autocorrelation of vege.",
                                        "rho_p"="Vegetation cover",
                                        "PL_expo"="Power-law exp. of the PSD",
                                        "Struct1"="PC 1 spa. struc.",
                                        "Struct2" = "PC 2 spa. struc.",
                                        "Spectral_ratio"="Spectral ratio",
                                        "fmax_psd"="log (largest patch)")))
}

boot_function_lm = function(formula, data, indices) {
  d = data[indices,] 
  fit = lm(formula, data=d) 
  return(summary(fit)$coefficient[2,1])
}


# 5) Model functions ----



Get_params_model=function(delta=.1,b=.49,c=.2,m=.001,d=.1,r=0.001,f=.9,g0=0,z=4){
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

PA_version_model=function(t,state,param){
  
  names(state)=c("rho_pp","rho_pm","rho_mm","rho_p","rho_m")
  # # print(state)
  # if (any(state<0)){
  #   state[which(state<0)]=0
  # }
  # state[state<1e-5]=0
  with (as.list(c(state, param)),{
    
    
    rho_0  = 1     - rho_p  - rho_m
    rho_p0 = rho_p - rho_pp - rho_pm
    rho_0m = rho_m - rho_mm - rho_pm
    
    drho_pp=2 * rho_p0 * ((delta * rho_p)+(1-delta)/z + ((z-1)/z)*(1-delta)*(rho_p0/rho_0))* (b-c*rho_p) -
      2 * rho_pp * (m + g0 * (1-(rho_pp/rho_p)) * ((z-1)/z) + g0/z)
    
    drho_pm=rho_0m * ((delta * rho_p)+ ((z-1)/z)*(1-delta) *(rho_p0/rho_0))* (b-c*rho_p) +
      rho_p0 *  d -
      rho_pm * (r + f * (rho_pm / rho_m) * ((z-1)/z)  + f/z) -
      rho_pm * (m + g0 * (1-(rho_pp/rho_p))*((z-1)/z))
    
    
    drho_mm=2*rho_0m*d - 
      2 * rho_mm * (r + f * (rho_pm / rho_m) * ((z-1)/z) )
    
    drho_p=rho_0*((delta * rho_p)+(1-delta) *(rho_p0/rho_0))* (b-c*rho_p) -
      rho_p* (m + g0 * (1-(rho_pp/rho_p)))
    
    drho_m=rho_0 * d -
      rho_m * (r + f * (rho_pm / rho_m))     
    
    list(c(drho_pp,drho_pm,drho_mm,drho_p,drho_m))
  })
  
}

Compute_ode=function(param,method_ode="lsoda",ini_cover=c(0.8,.1,.1),time_max=500){
  state=c(ini_cover[1]*ini_cover[1],ini_cover[2]*ini_cover[1],ini_cover[2]*ini_cover[2],ini_cover[1],ini_cover[2])
  
  
  time=seq(0,time_max,1)
  dynamics = as.data.frame(ode(state,time,func=PA_version_model,parms=param ,method = method_ode))
  dynamics=dynamics[,c("time","4")]
  colnames(dynamics)=c("Time","Rho_p")
    
  return(as_tibble(dynamics))
}

Plot_landscape=function(landscape){
  landscape[landscape<1]=0
  image(landscape,xaxt = "n",yaxt ="n",col=rev(c("black","white")) )
}
