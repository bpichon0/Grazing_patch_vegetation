x = c("tidyverse", "ggpubr", "latex2exp", "reshape2", "simecol","ggeffects","nlme",
      "abc", "spatialwarnings", "FME","phaseR","xgboost","igraph","MultiVarSel",
      "ggquiver", "scales","boot","RColorBrewer","ggnewscale","cluster","vegan",
      "factoextra","FactoMineR","missMDA","GGally","diptest","raster","ape","abctools","viridis",
      "png","jpeg","landscapemetrics","lme4","lmeresampler","GGally","MuMIn","multcompView",
      "LMERConvenienceFunctions","piecewiseSEM","qgraph","car","spdep","ParBayesianOptimization",
      "ggpattern")

source("https://gist.githubusercontent.com/benmarwick/2a1bb0133ff568cbe28d/raw/fb53bd97121f7f9ce947837ef1a4c65a73bffb3f/geom_flat_violin.R")

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


the_theme2 = theme_classic() + theme(
  legend.position = "bottom",
  panel.border = element_rect(colour = "black", fill=NA),
  strip.background = element_rect(fill = "transparent",color="transparent"),
  strip.text.y = element_text(size = 10, angle = -90),
  strip.text.x = element_text(size = 10),title = element_text(size=8),
  axis.title.y=element_text(size = 10),
  axis.title.x=element_text(size = 10),
  #legend.box="vertical",
  legend.text = element_text(size = 10), text = element_text(family = "Helvetica Light")
)


# scale_color_manual(values=c("black","green","yellow","red"))+
# cor(d_data[,"contig"],d_data[,c("rho_p","Org_C","Sand","Clim1","Clim2","Clim3","Clim4","Woody","Longitude","Slope","Elevation","Long_cos","Long_sin","Grazing")])

dir.create("../Figures/",showWarnings = F)
dir.create("../Figures/Final_figs",showWarnings = F)
dir.create("../Figures/Final_figs/SI",showWarnings = F)


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
    
    if(length(psd$psd_obs>0)){
      ks_dist  = Get_KS_distance(landscape>0,n_shuffle = 199)
    }else{
      ks_dist  = NA
    }
    

    #flow length
    flow_length=flowlength_sews(landscape>0,slope = slope,
                                cell_size = 50/sqrt(dim(landscape)[1]*dim(landscape)[2]))
    
    #power relationship between area and perimeter
    beta=lsm_c_pafrac(raster(landscape), directions = 4, verbose = TRUE)%>%filter(., class==1)%>%pull(., value)
    
    #mean perimeter-area ratio 
    mean_perim_area=lsm_c_para_mn(raster(landscape), directions = 4)%>%filter(., class==1)%>%pull(., value)
    
    #median, mean, sd of patch size (not in pixel unit but in m? to account for differences in resolutions)
    psd_scaled=lsm_p_area(raster(landscape),directions = 4)%>%filter(., class==1)%>%pull(., value)
    
    mean_psd=mean(psd_scaled)*1e4 #1e4 is to convert from ha to m?
    median_psd=median(psd_scaled)*1e4 #1e4 is to convert from ha to m?
    sd_psd=sd(psd_scaled)*1e4 #1e4 is to convert from ha to m?

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
             KS_dist=ks_dist, #Kolmogorov distance with null expectation
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
             KS_dist=0,
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

Closer_to_normality=function(df){
  
  df=df%>%
    mutate(., 
           flow_length=log(flow_length)+1,
           cv_psd=log(cv_psd)+1,
           mean_psd=log(mean_psd)+1,
           sd_psd=log(sd_psd)+1,
           median_psd=log(median_psd)+1,
           core_area=log(core_area)+1,
           Org_C=log(Org_C+5)+1,
           Slope=log(Slope+1)+1)
  
  return(df)
}

# 4) Statistical analysis functions ----

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
             Org_C=ifelse("Org_C" %in% colnames(importance_df),importance_df[,"Org_C"],0),
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
  
  if (type=="predictor"){ #i.e. plots with standardize effects 
    
    df=df%>%
      mutate(., Stat=recode_factor(Stat,"fmax"="Max patch","perim_area_scaling"="Perim-area scaling",
                                   "PL"="PSD exponent","core_area"="Core-area" ))%>%
      mutate(.,Type_pred=unlist(sapply(1:nrow(.),function(x){
        if (.$term[x] %in% c("rho_p","Type_vegGrass_Shrub","Type_vegGrassland","Type_vegShrubland","Woody")) return("Vegetation")
        if (.$term[x] %in% c("Clim1","Clim2","Clim3","Clim4")) return("Climatic")
        if (.$term[x] %in% c("Org_C","Sand","Slope")) return("Abiotic")
        if (.$term[x] %in% unique(df$term)[grep("Grazing",unique(df$term))]) return("Grazing")
        if (.$term[x] %in% c("Long_cos","Long_sin","Lattitude","Elevation")) return("Geography")
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
      dplyr::arrange(., Type_pred,Median)%>%
      dplyr::arrange(.,-dplyr::row_number())%>%
      add_column(., Order_f=1:nrow(.))%>%
      mutate(.,term = fct_reorder(term, Order_f))
    
  } else{
    
    df=df%>%
      mutate(.,Type_pred=unlist(sapply(1:nrow(.),function(x){
        if (.$term[x] %in% c("rho_p","Type_vegGrass_Shrub","Type_vegGrassland","Type_vegShrubland","Woody")) return("Vegetation")
        if (.$term[x] %in% c("Clim1","Clim2","Clim3","Clim4")) return("Climatic")
        if (.$term[x] %in% c("Org_C","Sand","Slope")) return("Abiotic")
        if (.$term[x] %in% unique(df$term)[grep("Grazing",unique(df$term))]) return("Grazing")
        if (.$term[x] %in% c("Long_cos","Long_sin","Lattitude","Elevation")) return("Geography")
      })))
    
  }
  return(df)
}

Plot_SEM_with_cover=function(summary_sem,pdf_=F,name="SEM",title_="",name_var="",MF=F,pval=.05){
  
  l=summary_sem$coefficients[,c("Predictor","Response","Estimate","P.Value")]
  l$color="#9ED471"
  l$color[which(l$Estimate<0)]="#EFC46A"
  l$color[which(l$P.Value>pval)]="gray"
  g = graph.data.frame(l, directed=T)
  g= g %>% set_edge_attr("color", value =l$color)
  
  coord=data.frame(label=c("Grazing","Woody","Cover","Facilitation",paste0(name_var)),
                   lab2=c("Grazing","Woody","Cover","Facilitation",paste0(name_var)),
                   x=c(-20,40,-40,20,0),y=c(0,-20,-20,0,-40))
  
  EL=as_edgelist(g)
  EL=cbind(EL,l[,3])
  EL=as.data.frame(EL)%>%
    mutate(., V1=recode_factor(V1,"rho_p"="Cover","Org_C"="Facilitation","Resid_mod"=paste0(name_var)))%>%
    mutate(., V2=recode_factor(V2,"rho_p"="Cover","Org_C"="Facilitation","Resid_mod"=paste0(name_var)))
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

Get_indirect_effects_grazing=function(summary_sem){
  
  coef=summary_sem$coefficients
  
  total_indirect=0
  if (all(coef$P.Value[c(7,3)]<.05)){
    total_indirect=total_indirect+prod(coef$Estimate[c(7,3)])
  }
  if (all(coef$P.Value[c(9,8,3)]<.05)){
    total_indirect=total_indirect+prod(coef$Estimate[c(9,8,3)])
  }
  if (all(coef$P.Value[c(5,2)]<.05)){
    total_indirect=total_indirect+prod(coef$Estimate[c(2,5)])
  }
  if (all(coef$P.Value[c(9,6,2)]<.05)){
    total_indirect=total_indirect+prod(coef$Estimate[c(9,6,2)])
  }
  if (all(coef$P.Value[c(9,4)]<.05)){
    total_indirect=total_indirect+prod(coef$Estimate[c(9,4)])
  }
  
  Direct_grazing=ifelse(coef$P.Value[c(1)]<.05,coef$Estimate[c(1)],0)

  return(tibble(Path=c("Indirect","Direct"),
                Effect=c(total_indirect,Direct_grazing)) )
  
}

Perform_PCA_spatial_struc=function(df){
  save=df
  df=df%>%
    dplyr::rename("Conditional entropy"="Cond_H",
                  "Contiguity"="contig",
                  "% landscape covered by core"="core_area_land",
                  "Mean % of core in patches"="core_area",
                  "Division"="division",
                  "Bare soil connectivity"="flow_length",
                  "Fractal dimension"="fractal_dim",
                  "Fractal scaling area, perim."="perim_area_scaling",
                  "PLR"="PLR",
                  "Distance to null expect."="KS_dist",
                  "Power-law exp. of the PSD"="PL_expo",
                  "Spectral ratio"="Spectral_ratio",
                  "log (largest patch)"="fmax_psd")
  
  
  #Performing PCA on the spatial structure metrics 
  struct_variables=colnames(df)[c(8:9,11,14:15,21,22,25)]
  res.comp=imputePCA(df[,which(colnames(df) %in% struct_variables)],ncp=4,scale = T) 
  
  if ("completeObs" %in% names(res.comp)){
    res.pca=PCA(res.comp$completeObs, ncp = 4,  graph=F)
  }else {
    res.pca=PCA(res.comp, ncp = 4,  graph=F)
  }
  
  #ploting the PCA
  p=factoextra::fviz_pca_var(res.pca,col.var="#2746B1",
                             label.var=c("Contiguity","% landscape covered by core",
                                         "Mean % of core in patches","Bare soil connectivity",
                                         "Fractal scaling area, perim.","PLR","Power-law exp. of the PSD",
                                         "log (largest patch)"))+
    the_theme+
    ggtitle("")+
    labs(x=paste0("PC 1 (",round(res.pca$eig[1,2],1),")"),
         y=paste0("PC 2 (",round(res.pca$eig[2,2],1),")"))
  ggsave("../Figures/Final_figs/SI/PCA_spatial_statistics.pdf",p,width = 7,height = 5)
  
  #We extract the first 2 ones (2/3 of the variance) and add them to the data-frame
  save=save%>%
    add_column(.,Struct1=res.pca$ind$coord[,1],Struct2=res.pca$ind$coord[,2])
  return(save)
  
}

Rename_spatial_statistics=function(df){
  return(df%>%
           mutate(., Stat=recode_factor(Stat,
                                        "Cond_H"="Conditional entropy",
                                        "contig"="Contiguity",
                                        "core_area_land"="% landscape covered by core",
                                        "core_area"="Mean % of core in patches",
                                        "division"="Division",
                                        "flow_length"="Bare soil connectivity",
                                        "fractal_dim"="Fractal dimension",
                                        "perim_area_scaling"="Fractal scaling area, perim.",
                                        "PLR"="PLR",
                                        "KS_dist"="Distance to null expect.",
                                        "PL_expo"="Power-law exp. of the PSD",
                                        "Struct1"="PC 1 spa. struc.",
                                        "Struct2" = "PC 2 spa. struc.",
                                        "Spectral_ratio"="Spectral ratio",
                                        "fmax_psd"="log (largest patch)")))
}

Get_data_without_outliers=function(stat,with_cover=T){
  
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
  return(d_data_out)
}

Filter_relevant_stats=function(df){
  return(df%>%filter(., Stat %in% c("PLR","PL_expo","fmax_psd","flow_length",
                                    "perim_area_scaling","core_area","contig",
                                    "core_area_land","KS_dist")))
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
