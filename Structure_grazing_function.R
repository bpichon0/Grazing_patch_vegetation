x = c("tidyverse", "ggpubr", "latex2exp", "reshape2", "simecol",
      "abc", "spatialwarnings", "FME","phaseR",
      "ggquiver", "scales","boot","RColorBrewer","ggnewscale","cluster","pls",
      "factoextra","FactoMineR","missMDA","GGally","diptest","raster","ape","abctools","viridis",
      "png","jpeg","landscapemetrics")

#install pacakges if not installed already
install.packages(setdiff(x, rownames(installed.packages())))

#loading the packages
lapply(x, require, character.only = TRUE)



d_biocom_old=read.table("../Data/biocom_data.csv",sep=";")
d_Meta=read.table("../Data/Meta_data_sites.csv",sep=",",header = T)
d_biocom=readxl::read_xlsx("../Data/Final_biocom.xlsx")
d_biodesert=readxl::read_xlsx("../Data/Final_biodesert.xlsx")

the_theme=theme_classic()+theme(legend.position = "bottom",
                                strip.background = element_rect(fill="white",color="white"),
                                strip.text.y = element_text(size = 10, angle = -90, face = "italic"),
                                strip.text.x = element_text(size = 10, face = "italic"),
                                legend.text = element_text(size = 10),text = element_text(family = "NewCenturySchoolbook"))

# Getting matrices and plotting functions ----

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


# Image analysis ----
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


# Spatial metrics functions ----

Get_sumstat=function(landscape,log_=T,all=F,slope=0){
  
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
  factal_dim=lsm_c_frac_mn(raster(landscape), directions = 4)%>%filter(., class==1)%>%pull(., value)
  
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
           factal_dim=factal_dim, #fractal dimension
           contig=contig, #mean connectedness of cells in patches
           Cond_H=complex_land$value[1], #conditional entropy
           Shannon_H=complex_land$value[2], #shannon entropy
           Joint_H=complex_land$value[3], #joint entropy
           mutual_inf=complex_land$value[4], #mutual information
           relat_mutual_inf=complex_land$value[5] #relative mutual information
  )
  
  
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

