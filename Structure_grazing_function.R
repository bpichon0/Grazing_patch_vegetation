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


d_biocom_old=read.table("../Data/biocom_data.csv",sep=";")
d_Meta=read.table("../Data/Meta_data_sites.csv",sep=",",header = T)
d_biocom=readxl::read_xlsx("../Data/Final_biocom.xlsx")
d_biodesert=readxl::read_xlsx("../Data/Final_biodesert.xlsx")
d_traits=readxl::read_xlsx("../Data/Traits/Drypop_biodesert.xlsx")
save_traits=d_traits
meta_traits=readxl::read_xlsx("../Data/Traits/Drypop_biodesert.xlsx",sheet = 2)[,1:2]

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

dir.create("../Figures/",showWarnings = F)
dir.create("../Figures/SI",showWarnings = F)

scale_values=function(x){return((x-mean(x,na.rm=T))/sd(x, na.rm = T))}

`%!in%` = Negate(`%in%`)

is_signif=function(x){
  if (x<.001){
    return("***")
  }else if (x>.001 & x<.01){
    return("**")
  }else if (x>.01 & x<.05){
    return("*")
  }else if (x>.05 & x<.1){
    return(".")
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

Get_empirical_site_biocom=function(id){
  d_biocom=read.table("../Data/biocom_data.csv",sep=";")
  return(as.matrix(read.table(paste0("../../Linking_data_model_ABC/Data/Data_Biocom/landscapes/",d_biocom$File_ID[id],".txt"))))
}

Get_empirical_site_all=function(ERC="biocom",Size,id,sub_id,cropped="cropped"){
  
  if (cropped!=""){
    return(as.matrix(read.table(paste0("../Data/Landscapes/Binary_landscapes/",
                                       ERC,"_",Size,"_",cropped,"_",id,"_",sub_id,".csv"),sep=",")))
  }else{
    return(as.matrix(read.table(paste0("../Data/Landscapes/Binary_landscapes/",
                                       ERC,"_",Size,"_",id,"_",sub_id,".csv"),sep=",")))
  }
}

Get_empirical_site_biodesert=function(ID){
  
  info_kmean=read.table("../Data/Landscapes/Kmean_clust_info.csv",sep=",",header = T)%>%
    dplyr::filter(., Size!=200,Dataset=="biodesert",status=="kept") #keeping the kept sites
  
  return(as.matrix(read.table(paste0("../Data/Landscapes/Binary_landscapes/biodesert_",
                                     info_kmean$Size[ID],"_",info_kmean$Site[ID],"_",info_kmean$Image[ID],".csv"),
                              sep=",")))
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

Get_png_empirical_site_biodesert=function(ID){
  info_kmean=read.table("../Data/Landscapes/Kmean_clust_info.csv",sep=",",header = T)%>%
    dplyr::filter(., Size!=200,Dataset=="biodesert",status=="kept") #keeping the kept sites
  
  img=readJPEG(paste0("../Data/Landscapes/biodesert/All/",ID,"_",info_kmean$Image[ID],".jpeg"))
  return( grid::grid.raster(img))
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
          the_theme+
          labs(x="Patch size",y="Number of patches")+
          scale_x_log10()+
          scale_y_log10())
  
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
             Nitrate=log(Nitrate),
             Org_C=log(Org_C+5), # 5 is the rounded minimum of the function
             Org_C_v=log(Org_C_v), 
             Total_N=log(Total_N),
             lnTotal_N=log(lnTotal_N+2), 
             Total_P=log(Total_P),
             lnTotal_P=log(lnTotal_P+100), 
             lnNitrate=log(lnNitrate+100),
             Grass=log(Grass+0.01),
             Fertility=sqrt(Fertility),
             Herb=log(Herb+0.01),
             C_stock=log(C_stock),
             Forage_Quality=sqrt(Forage_Quality+0.01),
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
             Nitrate=log(Nitrate),
             Org_C=log(Org_C+5), # 5 is the rounded minimum of the function
             Org_C_v=log(Org_C_v), 
             Total_N=log(Total_N),
             lnTotal_N=log(lnTotal_N+2), 
             Total_P=log(Total_P),
             lnTotal_P=log(lnTotal_P+100), 
             lnNitrate=log(lnNitrate+100),
             Grass=log(Grass+0.01),
             Fertility=sqrt(Fertility),
             Herb=log(Herb+0.01),
             C_stock=log(C_stock),
             Forage_Quality=sqrt(Forage_Quality+0.01),
             Slope=log(Slope+1))
    
  }
  df$Dung_herbivores=log(df$Dung_herbivores+.1)
  return(df)
}

Closer_to_normality_traits=function(df){
  
  df=df%>%
    mutate(., 
           MaxH=log(MaxH),
           LL=log(LL+.1),
           SLA=sqrt(SLA),
           MaxLS=log(MaxLS),
           # Maxvolume=log(Maxvolume),
           LA=log(LA),
           Phenolics=sqrt(Phenolics),
           LNC=sqrt(LNC),
           
           
    )
  return(df)
}

Closer_to_normality_CWM=function(df){
  
  
  return(df%>%
           mutate(., 
                  CWM_MaxLS=log(CWM_MaxLS+1),
           ))
}


z_tranform=function(x){
  x[is.infinite(x)]=NA
  return((x-mean(x,na.rm=T))/sd(x,na.rm=T))
}

Aggregate_importance=function(importance_df,aridity=F){
  
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
  
  if (aridity){

    d=tibble(Interactions=ifelse(exists("Interactions"),Interactions,0),
             Sand =ifelse("Sand" %in% colnames(importance_df),importance_df[,"Sand"],0),
             Aridity=ifelse("Aridity" %in% colnames(importance_df),importance_df[,"Aridity"],0),
             Org_C=ifelse("Org_C_v" %in% colnames(importance_df),importance_df[,"Org_C_v"],0),
             Type_veg=ifelse("Type_veg" %in% colnames(importance_df),importance_df[,"Type_veg"],0),
             Grazing=ifelse("Grazing" %in% colnames(importance_df),importance_df[,"Grazing"],0),
             Cover=ifelse("rho_p" %in% colnames(importance_df),importance_df[,"rho_p"],0),
             Woody=ifelse("Woody" %in% colnames(importance_df),importance_df[,"Woody"],0),
    )
    
  }else {
    
    
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
    
    
  }
  
  return(d)  
  
}

Aggregate_models=function(subset_models){
  msList = get.models(subset_models, subset = delta < 2, method="REML")
  if (length(msList)>1){
    output = list(ms=subset_models, avg = summary(model.avg(msList)))
  }else{
    output = subset_models
  }
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


Organize_df=function(df,type="predictor"){
  
  if (any(df$term %in% c(".sigma","(Intercept)"))){
    df=dplyr::filter(df, term %!in% c(".sigma","(Intercept)"))
  }
  
  if (type=="predictor"){ #i.e. plots with standardize effects 
    
    df=df%>%
      mutate(., Stat=recode_factor(Stat,"fmax"="Max patch","perim_area_scaling"="Perim-area scaling",
                                   "PL"="PSD exponent","core_area"="Core-area" ))%>%
      mutate(.,Type_pred=unlist(sapply(1:nrow(.),function(x){
        if (.$term[x] %in% c("rho_p","Type_vegGrass_Shrub","Type_vegGrassland","Type_vegShrubland","Woody")) return("Vegetation")
        if (.$term[x] %in% c("Clim1","Clim2","Clim3","Clim4","Aridity")) return("Climatic")
        if (.$term[x] %in% c("Org_C","Sand","Slope","Org_C_v")) return("Abiotic")
        if (.$term[x] %in% unique(df$term)[grep("Grazing",unique(df$term))]) return("Grazing")
        if (.$term[x] %in% c("Long_cos","Long_sin","Lattitude","Elevation")) return("Geography")
      })))%>%
      mutate(., term=recode_factor(term,
                                   "rho_p"="Cover",
                                   "Aridity"="Aridity",
                                   "Type_vegGrassland"="Type, grassland",
                                   "Type_vegShrubland"="Type, shrubland",
                                   "Type_vegGrass_Shrub"="Type, both",
                                   "Org_C"="Organic C.",
                                   "Org_C_v"="Organic C.",
                                   "Long_sin"="Longitude (sin)",
                                   "Long_cos"="Longitude (cos)",
                                   "Grazing:Type_vegShrubland"="Grazing * Type, shrubland",
                                   "Grazing:Type_vegGrassland"="Grazing * Type, grassland",
                                   "Grazing:Type_vegGrass_Shrub"="Grazing * Type, both",
                                   "Grazing:Org_C"="Grazing * Organic C.",
                                   "Grazing:Org_C_v"="Grazing * Organic C.",
                                   "Clim1:Grazing"="Grazing * Clim1",
                                   "Clim2:Grazing"="Grazing * Clim2",
                                   "Clim3:Grazing"="Grazing * Clim3",
                                   "Clim4:Grazing"="Grazing * Clim4",
                                   "Aridity:Grazing"="Aridity * Grazing",
                                   "Grazing:Sand"="Grazing * Sand",
                                   "Grazing:Woody"="Grazing * Woody",
                                   "Grazing:rho_p"="Grazing * Cover"))%>%
      dplyr::arrange(., Type_pred,Median)%>%
      dplyr::arrange(.,-dplyr::row_number())%>%
      add_column(., Order_f=1:nrow(.))%>%
      mutate(.,term = fct_reorder(term, Order_f))
    
  } else{
    
    df=df%>%
      mutate(.,Type_pred=unlist(sapply(1:nrow(.),function(x){
        if (.$term[x] %in% c("rho_p","Type_vegGrass_Shrub","Type_vegGrassland","Type_vegShrubland","Woody")) return("Vegetation")
        if (.$term[x] %in% c("Clim1","Clim2","Clim3","Clim4","Aridity")) return("Climatic")
        if (.$term[x] %in% c("Org_C","Sand","Slope","Org_C_v")) return("Abiotic")
        if (.$term[x] %in% unique(df$term)[grep("Grazing",unique(df$term))]) return("Grazing")
        if (.$term[x] %in% c("Long_cos","Long_sin","Lattitude","Elevation")) return("Geography")
      })))
    
  }
  return(df)
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
                                        "Cond_H"="Conditional entropy",
                                        "contig"="Contiguity",
                                        "core_area_land"="% landscape covered by core pixels",
                                        "core_area"="Circularity of patches",
                                        "division"="Division",
                                        "mean_psd"="Mean patch size",
                                        "mean_perim_area"="Mean fraction",
                                        "flow_length"="Bare soil connectivity",
                                        "fractal_dim"="Fractal dimension",
                                        "perim_area_scaling"="Fractal scaling area, perim.",
                                        "PLR"="PLR",
                                        "Small_patches"="Number of smallest patches",
                                        "moran_I"="Spatial autocorrelation of vege.",
                                        "KS_dist"="Distance to null expect.",
                                        "rho_p"="Vegetation cover",
                                        "PL_expo"="Power-law exp. of the PSD",
                                        "Struct1"="PC 1 spa. struc.",
                                        "Struct2" = "PC 2 spa. struc.",
                                        "Shape_metric"="Patch shape complexity",
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

# 6) Trait functions ----

Get_CWT_traits=function(abundance,traits_site){
  
  if (sum(abundance)>1){abundance=abundance/sum(abundance)} #relative abundance
  
  if (any(abundance==0)){ #remove non-present species
    abundance=abundance[-which(abundance==0)]
  }
  traits_site=traits_site%>%
    add_column(., 
               Name_sp=paste0(.$Genus," ",.$Species))%>%
    dplyr::filter(., Name_sp %in% names(abundance))%>%
    dplyr::arrange(., Name_sp) #sorting the tibble to match abundances in the abundance vector
  
  if (length(abundance) != nrow(traits_site)){ #meaning that a species has no recorded trait
    abundance=abundance[-which(names(abundance) %!in% traits_site$Name_sp)] #we remove such species
    if (sum(abundance)){ #in case we remove a species that was alone
      abundance=sapply(abundance,function(x){return(NA)})
    }else{
      abundance=abundance/sum(abundance) # and rescale the abundances
    }
  }
  
  traits_site$Cover[which(traits_site$Name_sp %in% names(abundance))]=abundance #change species relative abundance
  
  #then compute CWM at the quadrat level
  if ("PC1_traits" %in% colnames(traits_site)){
    CWM_quadrat=melt(traits_site, 
                     measure.vars=c("LL","SLA","LDMC","LA","MaxH","MaxLS","d15N","d13C",
                                    "Phenolics","LNC","LCC",#"Maxvolume",
                                    "PC1_traits","PC2_traits","PC3_traits","PC4_traits"))%>%#for each trait 
      dplyr::group_by(.,Country,Site,Plot,variable)%>% 
      dplyr::summarise(., .groups = "keep",
                       CWM = sum(Cover*value,na.rm = T), # community weighted mean
                       CWM_var = sum((Cover/sum(Cover))*((value-sum((Cover/sum(Cover))*value,na.rm = T))**2),na.rm = T)  
    )
  }else{
    CWM_quadrat=melt(traits_site, 
                     measure.vars=c("LL","SLA","LDMC","LA","MaxH","MaxLS","d15N","d13C",
                                    #"Maxvolume",
                                    "Phenolics","LNC","LCC",))%>%#for each trait 
      dplyr::group_by(.,Country,Site,Plot,variable)%>% 
      dplyr::summarise(., .groups = "keep",
                       CWM = sum(Cover*value,na.rm = T), # community weighted mean
                       CWM_var = sum((Cover/sum(Cover))*((value-sum((Cover/sum(Cover))*value,na.rm = T))**2),na.rm = T)
      )
  }
  return(CWM_quadrat)
}

Get_CWT_life_forms=function(abundance,traits_site){
  
  if (sum(abundance)>1){abundance=abundance/sum(abundance)} #relative abundance
  
  if (any(abundance==0)){ #remove non-present species
    abundance=abundance[-which(abundance==0)]
  }
  traits_site=traits_site%>%
    add_column(., 
               Name_sp=paste0(.$Genus," ",.$Species))%>%
    dplyr::filter(., Name_sp %in% names(abundance))%>%
    dplyr::arrange(., Name_sp) #sorting the tibble to match abundances in the abundance vector
  
  if (length(abundance) != nrow(traits_site)){ #meaning that a species has no recorded trait
    abundance=abundance[-which(names(abundance) %!in% traits_site$Name_sp)] #we remove such species
    if (sum(abundance)){ #in case we remove a species that was alone
      abundance=sapply(abundance,function(x){return(NA)})
    }else{
      abundance=abundance/sum(abundance) # and rescale the abundances
    }
  }
  
  traits_site$Cover[which(traits_site$Name_sp %in% names(abundance))]=abundance #change species relative abundance

  return(traits_site%>% 
           dplyr::group_by(.,Country,Site,Plot)%>% 
           dplyr::do(., Diversity_community_levelLife_form_quadrat(.$Life_form,.$Cover))
  )
}

Get_diversity_indices=function(abundance_mat,traits_site,id_traits=c(8:18)){
  
  null_tibble=tibble(alpha_mean=1e9,alpha_var=1e9,alpha_mean_star=1e9,alpha_var_star=1e9,
                     beta_plot_composition=1e9,
                     Functional_dissimilarity_bar=1e9,Species_similarity_bar=1e9,Dissimilarity_gap_bar=1e9,beta_uniqueness_bar=1e9,
                     beta_redundancy_bar=1e9,Functional_dissimilarity_var=1e9,Species_similarity_var=1e9,Dissimilarity_gap_var=1e9,
                     beta_uniqueness_var=1e9,beta_redundancy_var=1e9,
                     Rao_beta=1e9,Rao_alpha=1e9,
                     Rao_Q_weighted_mean=1e9,Rao_Q_non_weighted_mean=1e9,
                     F_Richness_mean=1e9,
                     F_eveness_mean_weighted=1e9,F_Diversity_mean_weighted=1e9,F_Dispersion_mean_weighted=1e9,
                     F_eveness_mean_non_weighted=1e9,F_Diversity_mean_non_weighted=1e9,F_Dispersion_mean_non_weighted=1e9,
                     Rao_Q_weighted_var=1e9,Rao_Q_non_weighted_var=1e9,
                     F_Richness_var=1e9,F_eveness_var_weighted=1e9,
                     F_Diversity_var_weighted=1e9,F_Dispersion_var_weighted=1e9,
                     F_eveness_var_non_weighted=1e9,F_Diversity_var_non_weighted=1e9,
                     F_Dispersion_var_non_weighted=1e9
  )
  
  
  
  abundance_mat=apply(abundance_mat,2,function(x){return(x/sum(x,na.rm = T))}) #scaling into proportion
  
  if (any(is.na(colSums(abundance_mat)))){ 
    abundance_mat=abundance_mat[,-which(is.na(colSums(abundance_mat)))] 
  }
  
  if (any(is.na(rowSums(abundance_mat)))){ 
    abundance_mat=abundance_mat[-which(is.na(rowSums(abundance_mat))),] 
  }
  if (any(rowSums(abundance_mat)==0)){ 
    abundance_mat=abundance_mat[-which(rowSums(abundance_mat)==0),] 
  }
  if (any(colSums(abundance_mat)==0)){ 
    abundance_mat=abundance_mat[-which(colSums(abundance_mat)==0),] 
  }
  
  traits_site=traits_site%>%
    add_column(., 
               Name_sp=paste0(.$Genus," ",.$Species))%>%
    dplyr::filter(., Name_sp %in% rownames(abundance_mat))%>%
    dplyr::arrange(., Name_sp)%>% #sorting the tibble to match abundances in the abundance vector
    as.data.frame(.)
  
  #order matrix abundance
  abundance_mat=abundance_mat[order(rownames(abundance_mat)),]
  
  if (!is.null(dim(abundance_mat))){#only one species
    
    if (nrow(abundance_mat) > nrow(traits_site)){ #meaning that a species has no recorded trait
      abundance_mat=abundance_mat[-which(rownames(abundance_mat) %!in% traits_site$Name_sp),] #we remove such species
      if (is.null(dim(abundance_mat))){#only one species
        abundance_mat=(abundance_mat>0)+0
      }else{
        if (any(colSums(abundance_mat)==0)){ #in case we remove a species that was alone
          abundance_mat=abundance_mat[,-which(colSums(abundance_mat)==0)]
        }
        abundance_mat=apply(abundance_mat,2,function(x){return(x/sum(x,na.rm = T))}) # and rescale the abundances
      }
    }
    
    if (nrow(abundance_mat) < nrow(traits_site)){ #meaning that a species has no recorded abundance
      traits_site=traits_site[-which(traits_site$Name_sp  %!in% rownames(abundance_mat)),] #we remove such species
    }
    
    #Removing species without any measured traits
    if (any(rowSums(traits_site[,id_traits],na.rm = T)==0)){
      if (sum(rowSums(traits_site[,id_traits],na.rm = T)==0)!=nrow(traits_site)){
        which_sp=which(rowSums(traits_site[,id_traits],na.rm = T)==0)
        traits_site=traits_site[-which_sp,]
        abundance_mat=abundance_mat[-which_sp,]
        
        if (class(abundance_mat)[1]!="numeric"){
          abundance_mat=apply(abundance_mat,2,function(x){return(x/sum(x,na.rm = T))}) #rescaling
          
          if (any(apply(abundance_mat,2,function(x){return(sum(x,na.rm = T))})==0)){ #meaninig that we have removed species that were alone in a quadrat
            which_quadrats=which(apply(abundance_mat,2,function(x){return(sum(x,na.rm = T))})==0)
            abundance_mat=abundance_mat[,-which_quadrats]
          }
        }
      }
    }
    

    #Building a distance matrix from plant traits
    rownames(traits_site)=traits_site$Name_sp
    scale_traits=as.data.frame(scale(traits_site[,id_traits]))
    Dist_func = dist(scale_traits)
    Dist_func = Dist_func/max(Dist_func,na.rm = T)
    
    if (any(is.na(Dist_func))){
      #this means that at least two species have non common traits, 
      #we remove one by one the one having the least number of traits until there is not NA issue
      #this is generally the least abundant species
      
      
      while(any(is.na(Dist_func))){
        n_NA_species=sort(apply(as.matrix(Dist_func),2,function(x){return(sum(is.na(x)))}),decreasing = T)
        abundance_mat=abundance_mat[-which(traits_site$Complete_name==names(n_NA_species)[1]),]
        traits_site=traits_site[-which(traits_site$Complete_name==names(n_NA_species)[1]),]
        scale_traits=as.data.frame(scale(traits_site[,id_traits]))
        Dist_func = dist(scale_traits)
        Dist_func = Dist_func/max(Dist_func,na.rm = T)
      }
    }
    
    if (!is.null(dim(abundance_mat))){#only one species
      abundance_mat=apply(abundance_mat,2,function(x){return(x/sum(x,na.rm = T))}) # and rescale the abundances
      
      if (any(apply(abundance_mat,2,function(x){return(sum(x,na.rm = T))})==0)){ #Double checking
        which_quadrats=which(apply(abundance_mat,2,function(x){return(sum(x,na.rm = T))})==0)
        abundance_mat=abundance_mat[,-which_quadrats]
      }
      
      #Then, alpha functional diversity for each quadrat
      alpha_div=(FPdivparam(t(abundance_mat), Dist_func)) #Using the \alpha K index of Pavoine and Ricotta 2021
      alpha_div_bar=mean(alpha_div[,1],na.rm=T)
      alpha_div_var=var(alpha_div[,1],na.rm=T)
      
      alpha_div_star=(FPdivparam(t(abundance_mat), Dist_func)) #Using the\alpha K* index of Pavoine and Ricotta 2021
      alpha_div_star_bar=mean(alpha_div[,1],na.rm=T)
      alpha_div_star_var=var(alpha_div[,1],na.rm=T)
      
      #Then beta diversity in composition
      beta_plot_composition=dsimcom(t(abundance_mat),method=2)
      diag(beta_plot_composition)=NA
      beta_plot_composition=mean(beta_plot_composition,na.rm=T)
      
      #Then functional beta diversity 
      # beta_plot=betaUniqueness(t(abundance_mat), Dist_func)
      # 
      # D_KG=beta_plot$DKG  # Pairwise functional dissimilarities between plots 
      # diag(D_KG)=NA
      # 
      # S_BC=1-beta_plot$DR # Pairwise species similarity between plots of 
      # diag(S_BC)=NA
      # 
      # R_beta = beta_plot$DR- beta_plot$DKG # Pairwise dissimilarity gap between plots
      # diag(R_beta)=NA
      # 
      # beta_uniqueness = beta_plot$betaUniqueness # Pairwise beta uniqueness between plots
      # diag(beta_uniqueness)=NA
      # 
      # beta_redundancy = beta_plot$betaRedundancy # Pairwise beta redundancy between plots
      # diag(beta_redundancy)=NA
      
      D_KG_bar = NA  # mean(D_KG,na.rm=T)  #averaging
      S_BC_bar = NA  # mean(S_BC,na.rm=T) #averaging
      R_beta_bar = NA  # mean(R_beta,na.rm=T) #averaging
      beta_uniqueness_bar = NA  # mean(beta_uniqueness,na.rm=T) #averaging
      beta_redundancy_bar = NA  # mean(beta_redundancy,na.rm=T) #averaging
      
      D_KG_var = NA  # var(as.numeric(D_KG),na.rm = T) #variance
      S_BC_var = NA  # var(as.numeric(S_BC),na.rm = T) #variance
      R_beta_var = NA  # var(as.numeric(R_beta),na.rm = T) #variance
      beta_uniqueness_var = NA  # var(as.numeric(beta_uniqueness),na.rm = T) #variance
      beta_redundancy_var = NA  # var(as.numeric(beta_redundancy),na.rm = T) #variance
      
      
      Rao_number=try(EqRao(t(abundance_mat), Dist_func, formula = "QE",option = "eq"),silent = T)
      if (class(Rao_number)=="try-error"){
        Rao_beta  = NA
        Rao_alpha = NA
      }else{
        Rao_number=EqRao(t(abundance_mat), Dist_func, formula = "QE",option = "eq")
        Rao_beta  = Rao_number[which(rownames(Rao_number)=="beta"),1]
        Rao_alpha = Rao_number[which(rownames(Rao_number)=="alpha"),1]
      }
      
      # Order species and check for names in the rows
      rownames(traits_site)=paste0(traits_site$Genus," ",traits_site$Species)
      if (any(colSums(traits_site[,id_traits],na.rm = T)==0)){ #meaning only NA values
        id_traits=id_traits[-which(colSums(traits_site[,id_traits],na.rm = T)==0)]
      }
      traits_site=traits_site%>%dplyr::arrange(., Name_sp)
      abundance_mat=abundance_mat[order(rownames(abundance_mat)),]
      
      
      ## FUNCTIONAL DIV  
      
      #Functional diversity with FD, weighted by abundances
      funct_div=try(dbFD(traits_site[,id_traits],t(abundance_mat),w.abun = T,messages = F),silent = T)
      
      if (class(funct_div)=="try-error"){
        funct_div=try(dbFD(traits_site[,id_traits],t(abundance_mat),w.abun = T,messages = F,corr="cailliez"),silent = T)
      }
      if (class(funct_div)=="try-error"){
        funct_div=list(
          RaoQ=NA,
          FRic=NA,
          FEve=NA,
          FDiv=NA,
          FDis=NA
        )
      }
      
      if ("FDiv" %!in% names(funct_div)){
        funct_div$FDiv=NA
      }
      
      #Functional diversity with FD, nonweighted by abundances
      
      funct_div_non_weighted=try(dbFD(traits_site[,id_traits],t(abundance_mat),w.abun = F,messages = F),silent = T)
      
      if (class(funct_div_non_weighted)=="try-error"){
        funct_div_non_weighted=try(dbFD(traits_site[,id_traits],t(abundance_mat),w.abun = F,messages = F,corr="cailliez"),silent = T)
      }
      if (class(funct_div_non_weighted)=="try-error"){
        funct_div_non_weighted=list(
          RaoQ=NA,
          FRic=NA,
          FEve=NA,
          FDiv=NA,
          FDis=NA
        )
      }
      
      if ("FDiv" %!in% names(funct_div_non_weighted)){
        funct_div_non_weighted$FDiv=NA
      }
      
      return(tibble(alpha_mean=alpha_div_bar,
                    alpha_var=alpha_div_var,
                    alpha_mean_star=alpha_div_star_bar,
                    alpha_var_star=alpha_div_star_var,
                    #
                    beta_plot_composition=beta_plot_composition,
                    #
                    Functional_dissimilarity_bar=D_KG_bar,
                    Species_similarity_bar=S_BC_bar,
                    Dissimilarity_gap_bar=R_beta_bar,
                    beta_uniqueness_bar=beta_uniqueness_bar,
                    beta_redundancy_bar=beta_redundancy_bar,
                    Functional_dissimilarity_var=D_KG_var,
                    Species_similarity_var=S_BC_var,
                    Dissimilarity_gap_var=R_beta_var,
                    beta_uniqueness_var=beta_uniqueness_var,
                    beta_redundancy_var=beta_redundancy_var,
                    #
                    Rao_beta=Rao_beta,
                    Rao_alpha=Rao_alpha,
                    #
                    Rao_Q_weighted_mean=mean(funct_div$RaoQ,na.rm = T),
                    Rao_Q_non_weighted_mean=mean(funct_div_non_weighted$RaoQ,na.rm = T),
                    #
                    F_Richness_mean=mean(funct_div$FRic,na.rm = T),
                    #
                    F_eveness_mean_weighted=mean(funct_div$FEve,na.rm = T),
                    F_Diversity_mean_weighted=mean(funct_div$FDiv,na.rm = T),
                    F_Dispersion_mean_weighted=mean(funct_div$FDis,na.rm = T),
                    # non weighted
                    F_eveness_mean_non_weighted=mean(funct_div_non_weighted$FEve,na.rm = T),
                    F_Diversity_mean_non_weighted=mean(funct_div_non_weighted$FDiv,na.rm = T),
                    F_Dispersion_mean_non_weighted=mean(funct_div_non_weighted$FDis,na.rm = T),
                    #
                    #
                    #Same but with the variance
                    Rao_Q_weighted_var=var(funct_div$RaoQ,na.rm = T),
                    Rao_Q_non_weighted_var=var(funct_div_non_weighted$RaoQ,na.rm = T),
                    #
                    F_Richness_var=var(funct_div$FRic,na.rm = T),
                    #
                    F_eveness_var_weighted=var(funct_div$FEve,na.rm = T),
                    F_Diversity_var_weighted=var(funct_div$FDiv,na.rm = T),
                    F_Dispersion_var_weighted=var(funct_div$FDis,na.rm = T),
                    # non weighted
                    F_eveness_var_non_weighted=var(funct_div_non_weighted$FEve,na.rm = T),
                    F_Diversity_var_non_weighted=var(funct_div_non_weighted$FDiv,na.rm = T),
                    F_Dispersion_var_non_weighted=var(funct_div_non_weighted$FDis,na.rm = T)
                    
  
      ))
    }else{ #only one species
      return(null_tibble)
    }
  }else{ #only one species
    return(null_tibble)
  }
}


Diversity_community_level=function(traits_site,id_traits=c(8:18)){
  
  
  null_tibble=tibble(F_Richness_mean=1e9,
                     #
                     F_eveness_mean_weighted=1e9,
                     F_Diversity_mean_weighted=1e9,
                     F_Dispersion_mean_weighted=1e9,
                     F_RaoQ_mean_weighted=1e9,
                     # non weighted
                     F_eveness_mean_non_weighted=1e9,
                     F_Diversity_mean_non_weighted=1e9,
                     F_Dispersion_mean_non_weighted=1e9,
                     F_RaoQ_mean_non_weighted=1e9
                     
  )
  
  if (any(is.na(traits_site$Cover))){
    traits_site$Cover[which(is.na(traits_site$Cover))]=0
  }
  
  abundance=traits_site$Cover
  names(abundance)=traits_site$Complete_name
  
  if (sum(abundance)>1){abundance=abundance/sum(abundance)} #relative abundance
  
  if (any(abundance==0)){ #remove non-present species
    abundance=abundance[-which(abundance==0)]
  }
  
  abundance=abundance[order(names(abundance))]
  traits_site=traits_site%>%
    add_column(., 
               Name_sp=paste0(.$Genus," ",.$Species))%>%
    dplyr::filter(., Name_sp %in% names(abundance))%>%
    dplyr::arrange(., Name_sp) #sorting the tibble to match abundances in the abundance vector
  
  if (length(abundance) != nrow(traits_site)){ #meaning that a species has no recorded trait
    abundance=abundance[-which(names(abundance) %!in% traits_site$Name_sp)] #we remove such species
    if (sum(abundance)){ #in case we remove a species that was alone
      abundance=sapply(abundance,function(x){return(NA)})
    }else{
      abundance=abundance/sum(abundance) # and rescale the abundances
    }
  }
  traits_site$Cover[which(traits_site$Name_sp %in% names(abundance))]=abundance #change species relative abundance
  traits_site=as.data.frame(traits_site)
  
  if (length(abundance)>1){#only one species
    
    # Order species and check for names in the rows
    rownames(traits_site)=paste0(traits_site$Genus," ",traits_site$Species)
    if (any(colSums(traits_site[,id_traits],na.rm = T)==0)){ #meaning only NA values
      id_traits=id_traits[-which(colSums(traits_site[,id_traits],na.rm = T)==0)]
    }
    traits_site=traits_site%>%dplyr::arrange(., Name_sp)
    abundance=abundance[order(names(abundance))]
    
    abundance=matrix(abundance,ncol=1)
    rownames(abundance)=traits_site$Name_sp
    colnames(abundance)="Site"
    
    
    if (length(id_traits)>1){
      #Removing species without traits
      if (any(rowSums(traits_site[,id_traits],na.rm = T)==0)){
        which_sp=which(rowSums(traits_site[,id_traits],na.rm = T)==0)
        traits_site=traits_site[-which_sp,]
        abundance=abundance[-which_sp,]
        abundance=abundance/sum(abundance) #rescaling
      }
      
      ## FUNCTIONAL DIV  
      
      #Functional diversity with FD, weighted by abundances
      funct_div=try(dbFD(traits_site[,id_traits],t(abundance),w.abun = T,messages = F),silent = T)
      
      if (class(funct_div)=="try-error"){
        funct_div=try(dbFD(traits_site[,id_traits],t(abundance),w.abun = T,messages = F,corr="cailliez"),silent = T)
      }
      if (class(funct_div)=="try-error"){
        funct_div=list(
          RaoQ=NA,
          FRic=NA,
          FEve=NA,
          FDiv=NA,
          FDis=NA
        )
      }
      
      if ("FDiv" %!in% names(funct_div)){
        funct_div$FDiv=NA
      }
      
      
      #Functional diversity with FD, nonweighted by abundances
      
      funct_div_non_weighted=try(dbFD(traits_site[,id_traits],t(abundance),w.abun = F,messages = F),silent = T)
      
      if (class(funct_div_non_weighted)=="try-error"){
        funct_div_non_weighted=try(dbFD(traits_site[,id_traits],t(abundance),w.abun = F,messages = F,corr="cailliez"),silent = T)
      }
      if (class(funct_div_non_weighted)=="try-error"){
        funct_div_non_weighted=list(
          RaoQ=NA,
          FRic=NA,
          FEve=NA,
          FDiv=NA,
          FDis=NA
        )
      }
      
      if ("FDiv" %!in% names(funct_div_non_weighted)){
        funct_div_non_weighted$FDiv=NA
      }
      
      return(tibble(F_Richness_mean=mean(funct_div$FRic,na.rm = T),
                    #
                    F_eveness_mean_weighted=mean(funct_div$FEve,na.rm = T),
                    F_Diversity_mean_weighted=mean(funct_div$FDiv,na.rm = T),
                    F_Dispersion_mean_weighted=mean(funct_div$FDis,na.rm = T),
                    F_RaoQ_mean_weighted=mean(funct_div$RaoQ,na.rm = T),
                    # non weighted
                    F_eveness_mean_non_weighted=mean(funct_div_non_weighted$FEve,na.rm = T),
                    F_Diversity_mean_non_weighted=mean(funct_div_non_weighted$FDiv,na.rm = T),
                    F_Dispersion_mean_non_weighted=mean(funct_div_non_weighted$FDis,na.rm = T),
                    F_RaoQ_mean_non_weighted=mean(funct_div_non_weighted$RaoQ,na.rm = T)
      ))
    }else{
      return(null_tibble)
    }
  }else{ #only one species
    return(null_tibble)
  }
}


Life_form_quadrat=function(x,y){return(tibble(Woody=sum(y[which(x=="Wood")])/sum(y),
                                              Grass_woody=sum(y[which(x%in%c("Wood","Grass"))])/sum(y),
                                              Grass_herb=sum(y[which(x%in%c("Herb","herb","Fern","Grass"))])/sum(y),
                                              Herb_woody=sum(y[which(x%in%c("Wood","Herb","herb","Fern"))])/sum(y),
                                              Grass=sum(y[which(x=="Grass")])/sum(y),
                                              Herb=sum(y[which(x%in%c("Herb","herb","Fern"))])/sum(y)))}

Compute_pairwise_trait_distance=function(trait_vector){
  return(tibble(PwD=sum(sapply(1:length(trait_vector),function(x){
    sum(abs(trait_vector[x]-trait_vector[-x]))
  }))/2)) #since we do all pairwise, we divide by two between each pair is repeated twice
}

Compute_mean_var_qwt=function(mean_qwt,var_qwt,qwt_var){
  return(tibble(
    
    #mean qwt between quadrats
    mean_QWT_random=mean(mean_qwt),
    q025_QWT_random=quantile(mean_qwt,.025),
    q975_QWT_random=quantile(mean_qwt,.975),
    q05_QWT_random=quantile(mean_qwt,.05),
    q95_QWT_random=quantile(mean_qwt,.95),
    
    #var qwt between quadrats
    mean_QWT_var_random=mean(var_qwt),
    q025_QWT_var_random=quantile(var_qwt,.025),
    q975_QWT_Var_random=quantile(var_qwt,.975),
    q05_QWT_Var_random=quantile(var_qwt,.05),
    q95_QWT_Var_random=quantile(var_qwt,.95),
    
    #mean qwt-var between quadrats
    mean_QWVa_randomr=mean(qwt_var),
    q025_QWVar_random=quantile(qwt_var,.025),
    q975_QWVar_random=quantile(qwt_var,.975),
    q05_QWVar_random=quantile(qwt_var,.05),
    q95_QWVar_random=quantile(qwt_var,.95)
    )
  )
}

Compute_trait_variance=function(trait_vector){
  return(tibble(Variance=var(trait_vector))) #since we do all pairwise, we divide by two between each pair is repeated twice
}

Add_PCA_traits=function(d){
  #add the first 4 principal components of a PCA on all traits
  
  d_traits_norm=d%>%Closer_to_normality_traits(.)
  
  struct_variables=colnames(d_traits_norm)[8:18]
  res.comp=imputePCA(d_traits_norm[,which(colnames(d_traits_norm) %in% struct_variables)],ncp=4,scale = T) 
  
  if ("completeObs" %in% names(res.comp)){
    res.pca=PCA(res.comp$completeObs, ncp = 4,  graph=F)
  }else {
    res.pca=PCA(res.comp, ncp = 4,  graph=F)
  }
  
  return(d%>%
           add_column(., 
                      PC1_traits=res.pca$ind$coord[,1],
                      PC2_traits=res.pca$ind$coord[,2],
                      PC3_traits=res.pca$ind$coord[,3],
                      PC4_traits=res.pca$ind$coord[,4])%>%
           dplyr::relocate(., Life_form,.after =PC4_traits )) 
}

Add_palatability_traits=function(d){
  
  
  return(d%>%
           add_column(., 
                      Palatability)%>%
           dplyr::relocate(., Life_form,.after =PC4_traits )) 
}

Homogeneize_sites_name=function(Site_name){
  #Homogeneize site names across databases
  
  if (Site_name =="LaCampana"){
    Site_name="la Campana"
  }else if(Site_name == "Nyngan"){
    Site_name="Nyngam"
  }else if(Site_name == "Wuxin"){
    Site_name="Wuxing"
  }else if(Site_name == "Moztaza"){
    Site_name="Mostaza"
  }else if(Site_name == "Escourt"){
    Site_name="Escourt"
  }else if(Site_name == "Natab"){
    Site_name="Natab"
  }else if(Site_name == "Sandvel"){
    Site_name="Sandveld"
  }else if(Site_name == "Tsamsvlei"){
    Site_name="Tsamsvlei"
  }else if(Site_name == "Kruger"){
    Site_name="Kruger Park"
  }else if(Site_name == "Mara"){
    Site_name="MaraExperimentalFarm"
  }else if(Site_name == "ZaragozaSemiarido"){
    Site_name="Zaragoza semiarido"
  }
  return(Site_name)
}

Homogeneize_country_name=function(Country_name){
  #Homogeneize Country names across databases
  
  if (Country_name =="Kenia"){
    Country_name="Kenya"
  }else if(Country_name == "Kazajstan"){
    Country_name="Kazakhstan"
  }else if(Country_name == "Hungria"){
    Country_name="Hungary"
  }else if(Country_name == "Palestina"){
    Country_name="Palestine"
  }else if(Country_name == "Tunez"){
    Country_name="Tunisia"
  }
  return(Country_name)
}

Plot_distribution_CWM_obs_random=function(random,obs,type="CWM",type_plot="boxplot"){
  
  if (type_plot=="boxplot"){
    if (type=="CWM"){
      return(ggplot(rbind(random%>%mutate(., Randomization_id=as.character(Randomization_id)),obs))+
               geom_boxplot(aes(x=variable,y=CWM,group=interaction(variable,Randomization_id),fill=Randomization_id))+
               the_theme+
               scale_fill_manual(values=c("grey","lightblue"))+
               facet_wrap(.~variable,scales = "free"))
    }else{
      return(ggplot(NULL)+
               geom_boxplot(data=random,aes(x=variable,y=PwD,group=variable),fill="grey")+
               the_theme+
               geom_point(data=obs,aes(x=variable,PwD),color="lightblue",size=3)+
               facet_wrap(.~variable,scales = "free"))
    }
    
  }else{
    return(ggplot(rbind(ungroup(random),ungroup(obs)))+
             geom_histogram(aes(x=CWM,fill=as.character(Randomization_id),group=as.character(Randomization_id)))+
             scale_fill_manual(values=c("grey","lightblue"))+
             facet_wrap(.~variable,scales = "free")+
             the_theme)
  }
}
