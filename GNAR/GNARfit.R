# adapt GNARfit function to include different weighting schemes 

GNARfit_weighting <- function(vts = covid_cases, 
                              net = covid_net_delaunay_gnar, 
                              numeric_vertices = FALSE,
                              alphaOrder=2, 
                              betaOrder=c(1,1), 
                              fact.var=NULL,
                              globalalpha=TRUE, 
                              tvnets=NULL, 
                              netsstart=NULL, 
                              ErrorIfNoNei=TRUE, 
                              inverse_distance = TRUE, 
                              # Great circle distance matrix with county names as row- / colnames 
                              distance_matrix = dist_urbanisation %>% as.matrix(), 
                              # data frame with column CountyName and column weight
                              weight_index = population_weight, 
                              # data frame with column CountyName and its numerical encoding 
                              county_index = NULL
                              ){
  #input checks
  stopifnot(is.GNARnet(net))
  stopifnot(ncol(vts) == length(net$edges))
  stopifnot(alphaOrder > 0)
  stopifnot(floor(alphaOrder) == alphaOrder)
  stopifnot(length(betaOrder) == alphaOrder)
  stopifnot(floor(betaOrder) == betaOrder)
  if(!is.null(fact.var)){
    stopifnot(length(fact.var) == length(net$edges))
    # if(sum(fact.var %in% c(0,1))!=length(fact.var)){
    #   cat("More than two (0/1) factor variables not yet supported")
    # }
    # stopifnot(sum(fact.var %in% c(0,1))==length(fact.var))
  }
  stopifnot(is.matrix(vts))
  # if(!globalalpha){
  #   cat("Individual alphas not yet supported")
  # }
  # stopifnot(globalalpha)
  stopifnot(is.logical(globalalpha))
  if(!is.null(tvnets)){
    cat("Time-varying networks not yet supported")
  }
  stopifnot(is.null(tvnets))
  
  
  
  useNofNei <- 1
  #cat("Note: input net should contain distances (not weights)")
  #end of input checks
  
  # construct list of all input 
  frbic <- list(nnodes=length(net$edges),
                alphas.in=alphaOrder,
                betas.in=betaOrder,
                fact.var=fact.var,
                globalalpha=globalalpha, 
                xtsp=tsp(vts), 
                time.in=nrow(vts), 
                net.in=net, 
                final.in=vts[(nrow(vts)-alphaOrder+1):nrow(vts),])
  
  # create design matrix necessary to fit GNAR model
  # dimensions x number of parameters (|alphaOrder| + |betaOrder|) 
  # values are covid numbers multiplied by inverse distance weight 
  dmat <- GNARdesign_weighting(vts = vts, 
                               net = net, 
                               numeric_vertices = numeric_vertices,
                               alphaOrder = alphaOrder, 
                               betaOrder = betaOrder, 
                               fact.var = fact.var,
                               globalalpha = globalalpha, 
                               tvnets = tvnets, 
                               netsstart = netsstart, 
                               inverse_distance = inverse_distance,
                               distance_matrix = distance_matrix,
                               weight_index = weight_index, 
                               county_index = county_index
                               )
  # dimensions 
  if(ErrorIfNoNei){
    if(any(apply(dmat==0, 2, all))){
      parname <- strsplit(names(which(apply(dmat==0, 2, all)))[1], split=NULL)[[1]]
      betastage <- parname[(which(parname==".")+1) :(length(parname))]
      stop("beta order too large for network, use max neighbour set smaller than ", betastage)
    }
  }
  
  # start for window to fit
  predt <- nrow(vts)-alphaOrder
  
  yvec <- NULL
  # for every vertex, filter out the y values for the output variable
  for(ii in 1:length(net$edges)){
    yvec <- c(yvec, vts[((alphaOrder+1):(predt+alphaOrder)),ii])
  }
  
  # can handle missing data
  if(sum(is.na(yvec))>0){
    yvec2 <- yvec[!is.na(yvec)]
    dmat2 <- dmat[!is.na(yvec),]
    modNoIntercept <- lm(yvec2~dmat2+0)
    
  }else{
    modNoIntercept <- lm(yvec~dmat+0)
  }
  
  out <- list(mod=modNoIntercept, 
              sd=coef(summary(modNoIntercept))[, "Std. Error"],
              y=yvec, 
              dd=dmat, 
              frbic=frbic)
  class(out) <- "GNARfit"
  return(out)
}
