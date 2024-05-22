x = c("tidyverse", "ggpubr", "latex2exp", "reshape2", "simecol","ggeffects","nlme",
      "abc", "spatialwarnings", "FME","phaseR","xgboost","igraph","MultiVarSel",
      "ggquiver", "scales","boot","RColorBrewer","ggnewscale","cluster","vegan",
      "factoextra","FactoMineR","missMDA","GGally","diptest","raster","ape","abctools","viridis",
      "png","jpeg","landscapemetrics","lme4","lmeresampler","GGally","MuMIn","multcompView",
      "LMERConvenienceFunctions","semEff","piecewiseSEM","qgraph","car","spdep","ParBayesianOptimization",
      "ggpattern","ggpmisc","ggpp","gradientForest","extendedForest","rfPermute","A3",
      "psych","emmeans","vegan")

source("https://gist.githubusercontent.com/benmarwick/2a1bb0133ff568cbe28d/raw/fb53bd97121f7f9ce947837ef1a4c65a73bffb3f/geom_flat_violin.R")

#loading the packages
lapply(x, require, character.only = TRUE)



d_biocom_old=read.table("../Data/biocom_data.csv",sep=";")
d_Meta=read.table("../Data/Meta_data_sites.csv",sep=",",header = T)
d_biocom=readxl::read_xlsx("../Data/Final_biocom.xlsx")
d_biodesert=readxl::read_xlsx("../Data/Final_biodesert.xlsx")
d_traits=readxl::read_xlsx("../Data/Traits/Drypop_biodesert.xlsx")

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
    filter(., Size!=200,Dataset=="biodesert",status=="kept") #keeping the kept sites
  
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
    filter(., Size!=200,Dataset=="biodesert",status=="kept") #keeping the kept sites
  
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
    beta=lsm_c_pafrac(raster(landscape), directions = 8, verbose = TRUE)%>%filter(., class==1)%>%pull(., value)
    
    #mean perimeter-area ratio 
    mean_perim_area=lsm_c_para_mn(raster(landscape), directions = 8)%>%filter(., class==1)%>%pull(., value)
    
    #median, mean, sd of patch size (not in pixel unit but in m? to account for differences in resolutions)
    psd_scaled=lsm_p_area(raster(landscape),directions = 8)%>%filter(., class==1)%>%pull(., value)
    
    mean_psd=mean(psd_scaled)*1e4 #1e4 is to convert from ha to m?
    median_psd=median(psd_scaled)*1e4 #1e4 is to convert from ha to m?
    sd_psd=sd(psd_scaled)*1e4 #1e4 is to convert from ha to m?

    #core area metric
    core_area=lsm_c_cai_mn(raster(landscape),directions = 8)%>%filter(., class==1)%>%pull(., value)
    core_area_land=lsm_c_cpland(raster(landscape),directions = 8)%>%filter(., class==1)%>%pull(., value)
    
    #division of patches
    division=lsm_c_division(raster(landscape), directions = 8)%>%filter(., class==1)%>%pull(., value)
    
    #fractal dimension 
    fractal_dim=lsm_c_frac_mn(raster(landscape), directions = 8)%>%filter(., class==1)%>%pull(., value)
    
    #contig
    contig=lsm_c_contig_mn(raster(landscape),directions = 8)%>%filter(., class==1)%>%pull(., value)
    
    #shape
    Shape_metric=lsm_c_shape_mn(raster(landscape),directions = 8)%>%filter(., class==1)%>%pull(., value)
    
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
  
  return(df)
}
Closer_to_normality_traits=function(df){
  
  df=df%>%
    mutate(., 
           MaxH=log(MaxH),
           LL=log(LL+.1),
           SLA=sqrt(SLA),
           LA=log(LA),
           MaxLS=log(MaxLS),
           Maxvolume=log(Maxvolume),
           Phenolics=sqrt(Phenolics)
           )
  return(df)
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

Organize_df=function(df,type="predictor"){
  
  if (any(df$term %in% c(".sigma","(Intercept)"))){
    df=filter(df, term %!in% c(".sigma","(Intercept)"))
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
                  "Mean % of core in patches"="core_area",
                  "Division"="division",
                  "Bare soil connectivity"="flow_length",
                  "Fractal dimension"="fractal_dim",
                  "Number of smallest patches"="Small_patches",
                  "Spatial autocorrelation of vege."="moran_I",
                  "Fractal scaling \n area, perim."="perim_area_scaling",
                  "PLR"="PLR",
                  "Mean % of core pixels in patches"="core_area",
                  "Mean patch size"="mean_psd",
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
                                               "log (largest patch)","Bare soil connectivity",
                                               "Mean patch size","Number of the smallest patches",
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

Get_data_without_outliers=function(stat,with_cover=T){
  
  d_data=read.table("../Data/Spatial_structure_grazing.csv",sep=";")%>%
    Closer_to_normality(.)
  
  d_data[,c(1,9:29,38,40:44,47:ncol(d_data))] = 
    apply(d_data[,c(1,9:29,38,40:44,47:ncol(d_data))],2,z_tranform)
  
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
  return(d_data_out)
}

Filter_relevant_stats=function(df){
  return(df%>%filter(., Stat %in% c("PL_expo","fmax_psd","flow_length",
                                    "moran_I","core_area","Small_patches",
                                    "mean_psd")))
}

Rename_spatial_statistics=function(df){
  return(df%>%
           mutate(., Stat=recode_factor(Stat,
                                        "Cond_H"="Conditional entropy",
                                        "contig"="Contiguity",
                                        "core_area_land"="% landscape covered by core pixels",
                                        "core_area"="Mean % of core pixels in patches",
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

Get_data_resid_SEM=function(stat,grazing_intensity,control_cover=T){
  
  dir.create("../Data/SEM/",showWarnings = F)
  d_data=read.table("../Data/Spatial_structure_grazing.csv",sep=";")%>%
    Closer_to_normality(.)
  d_data[,c(1,9:29,36,38:42,45:ncol(d_data))] = apply(d_data[,c(1,9:29,36,38:42,45:ncol(d_data))],2,z_tranform)
  d_data=Perform_PCA_spatial_struc(d_data) #Adding the first two components from the PCA
  
  
  if (grazing_intensity=="low"){
    grazing_values=0:1
  } else if (grazing_intensity=="high"){
    grazing_values=2:3
  } else if (grazing_intensity=="all"){
    grazing_values=0:3
  } else if (grazing_intensity=="grazed"){
    grazing_values=1:3
  }
  
  k=stat
  if (control_cover){
    formula_model="value ~ Long_cos + Long_sin + Lattitude + Slope + Elevation + Type_veg + rho_p"
  }else{
    formula_model="value ~ Long_cos + Long_sin + Lattitude + Slope + Elevation + Type_veg"
  }
  
  model_lmer=lm(formula_model,
                data = d_data%>%melt(., measure.vars=k)%>%
                  filter(., !is.na(value),Grazing %in% grazing_values),
                na.action = na.omit)
  
  #we remove potential outliers
  rm.outliers = romr.fnc(model_lmer, d_data%>%melt(., measure.vars=k)%>%
                           filter(., !is.na(value),Grazing %in% grazing_values),
                         trim=2.5)
  d_data_out = rm.outliers$data
  
  #We first control for all covariates and extract the residuals
  model_lmer=lm(formula_model,
                data = d_data_out,
                na.action = na.fail)
  
  resid_model=residuals(model_lmer) #extract residuals
  
  save=d_data_out%>% #add it to the dataframe
    add_column(., Resid_mod=resid_model)%>%
    filter(., !is.na(Resid_mod))
  
  return(save)
}

# Plot_SEM=function(summary_sem,pdf_=F,name="SEM",title_=""
#                   ,name_var="",type_N="Nitrate",
#                   label_cex=1.5,edge_cex=2){
#   
#   l=as.data.frame(summary_sem$coefficients[,c("Predictor","Response","Estimate","P.Value")]%>%
#                     mutate(., 
#                            Predictor=recode_factor(Predictor,"Org_C_v"="Organic C.","Total_N"="Total N."),
#                            Response=recode_factor(Response,"Org_C_v"="Organic C.","Total_N"="Total N.")
#                     ))%>%
#     mutate(., Predictor=as.character(Predictor),Response=as.character(Response))
#   
#   l[nrow(l),1:2]=c("Aridity","Sand")
#   corr_arrow=l[nrow(l),c(2,1,3:ncol(l))];names(corr_arrow)=c("Predictor","Response","Estimate","P.Value")
#   l=rbind(l,corr_arrow)
#   
#   l$color="#9ED471"
#   l$color[which(l$Estimate<0)]="#EFC46A"
#   l$color[which(l$P.Value>.1)]="gray"
#   l$color[which(l$P.Value>.05 & l$P.Value<.1)]="lightblue"
#   g = graph.data.frame(l, directed=T)
#   g= g %>% set_edge_attr("color", value =l$color)
#   
#   coord=data.frame(label=c("Grazing","Organic C.",type_N,"Aridity","Sand",paste0(name_var)),
#                    lab2=c("Grazing","Organic C.",type_N,"Aridity","Sand",paste0(name_var)),
#                    x=c(10,-2,-2,-10,10,30),y=c(-10,25,-25,0,10,0))
#   
#   EL=as_edgelist(g)
#   EL=cbind(EL,l[,3])
#   EL=as.data.frame(EL)%>%
#     mutate(., V1=recode_factor(V1,"rho_p"="Cover","Org_C"="Organic C.","Resid_mod"=paste0(name_var)))%>%
#     mutate(., V2=recode_factor(V2,"rho_p"="Cover","Org_C"="Organic C.","Resid_mod"=paste0(name_var)))
#   EL=as.matrix(EL)
#   
#   name_edge=sapply(1:length(summary_sem$coefficients$P.Value),function(x){#adding pvalue
#     if (is_signif(summary_sem$coefficients$P.Value[x])==""){
#       return(paste0(round(l[x,3],2)))
#     }else{
#       return(paste0(round(l[x,3],2))) #,is_signif(summary_sem$coefficients$P.Value[x])))
#     }
#   })
#   name_edge=c(name_edge,name_edge[length(name_edge)])
#   
#   asi=abs(l[,3])/0.05
#   asi[asi<5]=5
#   
#   if(pdf_){
#     pdf(paste0("./",name,".pdf"),width = 9,height = 6)
#     qgraph(EL,layout=as.matrix(coord[,c("x","y")]),edge.color=l$color,edge.labels=name_edge,
#            border.color="white",label.cex=label_cex,label.scale=F,title=title_,
#            edge.label.cex = edge_cex,edge.label.position=0.5,vsize2=4,vsize=30,
#            shape="ellipse",edge.labels=T,fade=F,esize=5,asize=asi,
#            mar=rep(3,4))
#     dev.off()
#     
#   }else {
#     qgraph(EL,layout=as.matrix(coord[,c("x","y")]),edge.color=l$color,edge.labels=name_edge,
#            border.color="white",label.cex=label_cex,label.scale=F,title=title_,
#            edge.label.cex = edge_cex,edge.label.position=0.5,vsize2=4,vsize=30,
#            shape="ellipse",edge.labels=T,fade=F,esize=5,asize=asi,
#            mar=rep(3,4))
#     
#   }
# }
# 
# 
# Plot_SEM2=function(summary_sem,pdf_=F,name="SEM",title_=""
#                   ,name_var="",type_N="Nitrate",
#                   label_cex=1,edge_cex=1){
#   
#   l=as.data.frame(summary_sem$coefficients[,c("Predictor","Response","Estimate","P.Value")]%>%
#                     mutate(., 
#                            Predictor=recode_factor(Predictor,"Org_C_v"="Organic C."),
#                            Response=recode_factor(Response,"Org_C_v"="Organic C.")
#                     ))%>%
#     mutate(., Predictor=as.character(Predictor),Response=as.character(Response))
#   
#   l[nrow(l),1:2]=c("Cover","q")
#   corr_arrow=l[nrow(l),c(2,1,3:ncol(l))];names(corr_arrow)=c("Predictor","Response","Estimate","P.Value")
#   l=rbind(l,corr_arrow)
#   
#   l$color="#9ED471"
#   l$color[which(l$Estimate<0)]="#EFC46A"
#   l$color[which(l$P.Value>.1)]="gray"
#   l$color[which(l$P.Value>.05 & l$P.Value<.1)]="lightblue"
#   #l=l[-which(l$P.Value>.1),]
#   g = graph.data.frame(l, directed=T)
#   g= g %>% set_edge_attr("color", value =l$color)
#   
#   coord=data.frame(x=c(10,  #cover
#                        10, #q
#                        -30, #aridity
#                        -15, #organicC
#                        -30, #GRazing
#                        -15, #nitrate
#                        30), #resilience
#                    y=c(-10,
#                        10,
#                        -30,
#                        -7.5,
#                        30,
#                        7.5,
#                        0))
#   
#   EL=as_edgelist(g)
#   EL=cbind(EL,l[,3])
#   EL=as.data.frame(EL)%>%
#     mutate(., V1=recode_factor(V1,"rho_p"="Cover","Org_C"="Organic C.","Resid_mod"=paste0(name_var)))%>%
#     mutate(., V2=recode_factor(V2,"rho_p"="Cover","Org_C"="Organic C.","Resid_mod"=paste0(name_var)))
#   EL=as.matrix(EL)
#   
#   name_edge=sapply(1:length(summary_sem$coefficients$P.Value),function(x){#adding pvalue
#     if (is_signif(summary_sem$coefficients$P.Value[x])==""){
#       return(paste0(round(l[x,3],2)))
#     }else{
#       return(paste0(round(l[x,3],2))) #,is_signif(summary_sem$coefficients$P.Value[x])))
#     }
#   })
#   name_edge=c(name_edge,name_edge[length(name_edge)])
#   
#   asi=abs(l[,3])/0.05
#   asi[asi<5]=5
#   
#   
#   if(pdf_==F){
#     qgraph(EL,layout=as.matrix(coord[,c("x","y")]),edge.color=l$color,edge.labels=name_edge,
#            border.color="white",label.cex=label_cex,label.scale=F,title=title_,
#            edge.label.cex = edge_cex,edge.label.position=0.5,vsize2=4,vsize=30,
#            shape="ellipse",edge.labels=T,fade=F,esize=5,asize=asi,
#            mar=rep(3,4))
#   }else {
#     pdf(paste0("./",name,".pdf"),width = 9,height = 6)
#     
#     qgraph(EL,layout=as.matrix(coord[,c("x","y")]),edge.color=l$color,edge.labels=name_edge,
#            border.color="white",label.cex=label_cex,label.scale=F,title=title_,
#            edge.label.cex = edge_cex,edge.label.position=0.5,vsize2=4,vsize=30,
#            shape="ellipse",edge.labels=T,fade=F,esize=5,asize=asi,
#            mar=rep(3,4))
#     dev.off()
#     
#   }
# }
Plot_SEM3=function(summary_sem,pdf_=F,name="SEM",title_=""
                  ,name_var="",type_N="Nitrate",
                  label_cex=1,edge_cex=1){
  
  l=as.data.frame(summary_sem$coefficients[,c("Predictor","Response","Estimate","P.Value")]%>%
                    mutate(., 
                           Predictor=recode_factor(Predictor,"Org_C_v"="Organic C.","Total_N"="Total N."),
                           Response=recode_factor(Response,"Org_C_v"="Organic C.","Total_N"="Total N.")
                    ))%>%
    mutate(., Predictor=as.character(Predictor),Response=as.character(Response))
  
  corr_arrow=l[nrow(l),c(2,1,3:ncol(l))];names(corr_arrow)=c("Predictor","Response","Estimate","P.Value")
  l=rbind(l,corr_arrow)
  
  l$color="#9ED471"
  l$color[which(l$Estimate<0)]="#EFC46A"
  l$color[which(l$P.Value>.1)]="gray"
  l$color[which(l$P.Value>.05 & l$P.Value<.1)]="lightblue"
  g = graph.data.frame(l, directed=T)
  g= g %>% set_edge_attr("color", value =l$color)
  
  coord=data.frame(label=c("Grazing","Organic C.",type_N,"Aridity",paste0(name_var)),
                   lab2=c("Grazing","Organic C.",type_N,"Aridity",paste0(name_var)),
                   x=c(0,#orgC
                       -10,#grazing
                       -10,#aridity
                       0, #nitrate
                       10),#resilience
                   y=c(-5,
                       5,
                       -10,
                       10,
                       0))
  
  EL=as_edgelist(g)
  EL=cbind(EL,l[,3])
  EL=as.data.frame(EL)%>%
    mutate(., V1=recode_factor(V1,"rho_p"="Cover","Org_C"="Organic C.","Resid_mod"=paste0(name_var)))%>%
    mutate(., V2=recode_factor(V2,"rho_p"="Cover","Org_C"="Organic C.","Resid_mod"=paste0(name_var)))
  EL=as.matrix(EL)
  
  name_edge=sapply(1:length(summary_sem$coefficients$P.Value),function(x){#adding pvalue
    if (is_signif(summary_sem$coefficients$P.Value[x])==""){
      return(paste0(round(l[x,3],2)))
    }else{
      return(paste0(round(l[x,3],2))) #,is_signif(summary_sem$coefficients$P.Value[x])))
    }
  })
  name_edge=c(name_edge,name_edge[length(name_edge)])
  
  asi=abs(l[,3])/0.05
  asi[asi<5]=5
  
  if(pdf_){
    pdf(paste0("./",name,".pdf"),width = 9,height = 6)
    qgraph(EL,layout=as.matrix(coord[,c("x","y")]),edge.color=l$color,edge.labels=name_edge,
           border.color="white",label.cex=label_cex,label.scale=F,title=title_,
           edge.label.cex = edge_cex,edge.label.position=0.5,vsize2=4,vsize=30,
           shape="ellipse",edge.labels=T,fade=F,esize=5,asize=asi,
           mar=rep(3,4))
    dev.off()
    
  }else {
    qgraph(EL,layout=as.matrix(coord[,c("x","y")]),edge.color=l$color,edge.labels=name_edge,
           border.color="white",label.cex=label_cex,label.scale=F,title=title_,
           edge.label.cex = edge_cex,edge.label.position=0.5,vsize2=4,vsize=30,
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

Run_spatial_model=function(time=seq(1,2000,1),params,ini,plot=F,Fixed_cover=F,Cover_value=0){
  
  d=tibble(Time=1,
           Rho_V=sum(ini == 1) / length(ini),
           Rho_F=sum(ini == 0) / length(ini),
           Rho_D=sum(ini == -1) / length(ini))
  state=list(Landscape=ini,Rho_v=d$Rho_V,Rho_f=d$Rho_F,Rho_D=d$Rho_D)
  
  if (Fixed_cover){
    for (k in 2:300){
      
      params$dt=time[k]-time[k-1]
      
      state=CA_spatial_model(state$Landscape,params = params)
      d=rbind(d,tibble(Time=k,Rho_V=state$Rho_v,Rho_F=state$Rho_f,Rho_D=state$Rho_D))
      
    }
    
    while(state$Rho_v!=Cover_value){
      params$dt=time[k]-time[k-1]
      
      state=CA_spatial_model(state$Landscape,params = params)
      d=rbind(d,tibble(Time=k,Rho_V=state$Rho_v,Rho_F=state$Rho_f,Rho_D=state$Rho_D))
      k=k+1
    }
    
    
  }else{
    for (k in 2:length(time)){
      
      params$dt=time[k]-time[k-1]
      
      state=CA_spatial_model(state$Landscape,params = params)
      d=rbind(d,tibble(Time=k,Rho_V=state$Rho_v,Rho_F=state$Rho_f,Rho_D=state$Rho_D))
      
    }
    
  }
  
  return(list(d=d,State=state$Landscape)) #record time evolution and last landscape
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
    filter(., Name_sp %in% names(abundance))%>%
    dplyr::arrange(., Name_sp) #sorting the tibble to match abundances in the abundance vector
  
  if (length(abundance) != nrow(traits_site)){ #meaning that a species has no recorded trait
    abundance=abundance[-which(names(abundance) %!in% traits_site$Name_sp)] #we remove such species
    abundance=abundance/sum(abundance) # and rescale the abundances
  }

  traits_site$Cover[which(traits_site$Name_sp %in% names(abundance))]=abundance #change species relative abundance
  
  #then compute CWM at the quadrat level
  CWM_quadrat=melt(traits_site, 
                   measure.vars=c("LL","SLA","LDMC","LA","MaxH","MaxLS",
                                  "Maxvolume","Phenolics","LNC","LCC"))%>%#for each trait 
    dplyr::group_by(.,Country,Site,Plot,variable)%>% 
    dplyr::summarise(., .groups = "keep",
                     CWM = sum(Cover*value,na.rm = T) # community weighted mean
    )
  
  return(CWM_quadrat)
}


Compute_pairwise_trait_distance=function(trait_vector){
  return(tibble(PwD=sum(sapply(1:length(trait_vector),function(x){
    sum(abs(trait_vector[x]-trait_vector[-x]))
  }))/2 #since we do all pairwise, we divide by two between each pair is repeated twice
  ))
}


Add_PCA_traits=function(d){
  #add the first 3 principal components of a PCA on all traits
  
  d_traits_norm=d%>%Closer_to_normality_traits(.)
  
  struct_variables=colnames(d_traits_norm)[7:16]
  res.comp=imputePCA(d_traits_norm[,which(colnames(d_traits_norm) %in% struct_variables)],ncp=4,scale = T) 
  
  if ("completeObs" %in% names(res.comp)){
    res.pca=PCA(res.comp$completeObs, ncp = 4,  graph=F)
  }else {
    res.pca=PCA(res.comp, ncp = 4,  graph=F)
  }
  
  d_traits_norm=d_traits_norm%>%
    add_column(., 
               PC1=res.pca$ind$coord[,1],
               PC2=res.pca$ind$coord[,2],
               PC3=res.pca$ind$coord[,3])
  
 return(d_traits_norm) 
}



