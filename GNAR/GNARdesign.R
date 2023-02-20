# adapt GNARdesign function to include different weighting schemes 

GNARdesign_weighting <- function (vts = covid_cases, 
                                  net = covid_net_delaunay_gnar,
                                  numeric_vertices = FALSE, 
                                  alphaOrder = 2,
                                  betaOrder = c(1,1), 
                                  fact.var = NULL, 
                                  globalalpha = TRUE,
                                  tvnets = NULL, 
                                  netsstart = NULL, 
                                  inverse_distance = TRUE,
                                  # Great circle distance matrix with county names as row- / colnames 
                                  distance_matrix = dist_urbanisation %>% as.matrix(), 
                                  # data frame with column CountyName and column weight
                                  weight_index = population_weight, 
                                  # data frame with column CountyName and its numerical encoding 
                                  county_index = NULL
)
{
  stopifnot(is.GNARnet(net))
  stopifnot(ncol(vts) == length(net$edges))
  stopifnot(alphaOrder > 0)
  stopifnot(floor(alphaOrder) == alphaOrder)
  stopifnot(length(betaOrder) == alphaOrder)
  stopifnot(floor(betaOrder) == betaOrder)
  if(!is.null(fact.var)){
    stopifnot(length(fact.var) == length(net$edges))
    if(!globalalpha){
      stop("Use factors OR individual alphas")
    }
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
  
  # inverse_distance must be logical 
  stopifnot(is.logical(inverse_distance))
  
  
  #cat("Note: input net should contain distances (not weights)")
  #flip network so that NofNeighbours gives into node information
  netmat <- as.matrix(net, normalise=FALSE)
  if(!isSymmetric(netmat)){
    net <- as.GNARnet(t(netmat))
  }
  parNames <- parLoc <- NULL
  for (jj in 1:alphaOrder) {
    if(globalalpha){
      parNames <- c(parNames, paste("alpha", jj, sep = ""))
      parLoc <- c(parLoc, "a")
    }else{
      for(kk in 1:ncol(vts)){
        parNames <- c(parNames, paste("alpha", jj, "node", kk, sep=""))
        parLoc <- c(parLoc, "a")
      }
    }
    if (betaOrder[jj] > 0) {
      for (kk in 1:betaOrder[jj]) {
        parNames <- c(parNames, paste("beta", jj, ".",
                                      kk, sep = ""))
        parLoc <- c(parLoc, "b")
      }
    }
  }
  maxOrder <- alphaOrder
  predt <- nrow(vts) - maxOrder
  nnodes <- ncol(vts)
  if(globalalpha){
    dmat <- matrix(0, nrow = predt * nnodes, ncol = sum(c(alphaOrder,
                                                          betaOrder)), 
                   dimnames = list(NULL, parNames))
  }else{
    dmat <- matrix(0, nrow=predt*nnodes, ncol=sum(c(alphaOrder*nnodes,
                                                    betaOrder)), 
                   dimnames=list(NULL,parNames))
  }
  
  for (ii in 1:nnodes) {
    for (aa in 1:alphaOrder) {
      if(globalalpha){
        alphaLoc <- which(parLoc == "a")[aa]
      }else{
        alphaLoc <- which(parLoc=="a")[nnodes*(aa-1)+ii]
      }
      dmat[((predt * (ii - 1) + 1):(predt * ii)), alphaLoc] <- vts[((maxOrder + 1 - aa):(predt + (maxOrder - aa))), ii]
    }
  }
  if (sum(betaOrder) > 0) {
    betaN <- NULL
    betaTimes <- rep(1:alphaOrder, betaOrder)
    for (jj in 1:alphaOrder) {
      #betaTimes <- c(betaTimes, rep(jj, betaOrder[jj]))
      if (betaOrder[jj] > 0) {
        betaN <- c(betaN, 1:betaOrder[jj])
      }
    }
    
    # create county index 
    # if (!numeric_vertices) {
    #   county_index <- data.frame("index" = net$edges %>% unlist(), 
    #                              "CountyName" = net$edges %>% unlist() %>% names()) %>% 
    #     distinct()
    # }
    
    # adapt to incorporate Great Circle distances 
    for (ii in 1:nnodes) {
      NofNei <- NofNeighbours_named(node = ii, 
                                    stage = max(betaOrder),
                                    net = net, 
                                    numeric_vertices = numeric_vertices)
      # old: NofNei <- NofNeighbours(node = ii, 
      #                          stage = max(betaOrder), 
      #                          net = net)
      
      Nei <- NofNei$edges
      # require names for further neighbourhoods 
      # change in NofNei necessary 
  
      # name of current node
      vertex_name <- county_index[county_index$index == ii, ]$CountyName
      
      # must be done for all
      Wei <- list()
      
      # circle through all 
      for (degree_neighbourhood in seq(1, Nei %>% length())) {
        
        # filter distances between vertex and its neighbors
        if (numeric_vertices) {
          # name of neighbors 
          neighbour_names <- county_index[county_index$index %in% NofNei$edges[[degree_neighbourhood]], ]$CountyName
        } 
        if (! numeric_vertices) {
          neighbour_names <- names(NofNei$edges[[degree_neighbourhood]])
        }
        # not neighbour to itself 
        neighbour_names <- neighbour_names[neighbour_names != vertex_name]
        
        # distance between counties in km and assign names 
        Wei_v <- distance_matrix[vertex_name, neighbour_names] / 1000
        names(Wei_v) <- neighbour_names
        
        if (Wei_v %>% length() == 0) {
          Wei_v <- NA
        }
        
        if (inverse_distance) {
          weight_v <- 1
        }
        
        if (!inverse_distance) {
          weight <- weight_index[weight_index$CountyName %in% neighbour_names, ]
          weight_v <- weight[match(names(Wei_v), weight$CountyName), ]$weight
        }
        
        # inverse multiplication with population number for each county
        Wei[[degree_neighbourhood]] <- Wei_v * 1/weight_v
      }
      
      # old: Wei <- NofNei$dist
      
      if ((!is.null(Nei)) & (length(Nei) > 0)) {
        if (!is.null(Nei[[1]])&!is.na(Nei[[1]][1])) {
          
          # inverse distance weighting * other desired weighting
          Wei_w <- lapply(Wei, function(x){1/(x * sum(1/x))})
          
          for (bb in 1:sum(betaOrder)) {
            betaLoc <- which(parLoc == "b")[bb]
            if (length(Nei[[betaN[bb]]]) > 1) {
              # print(paste("node", ii, "betaN[bb]", betaN[bb]))
              # print("In length(Nei[[betaN[bb]]]) > 1")
              vts.cut <- vts[((maxOrder + 1 - betaTimes[bb]):(predt +
                                                                (maxOrder - betaTimes[bb]))), Nei[[betaN[bb]]]]
              for (kk in 1:nrow(vts.cut)) {
                if (any(is.na(vts.cut[kk, ]))) {
                  if (all(is.na(vts.cut[kk, ]))) {
                    #if there are no neighbours left at any time point, set to zero
                    vts.cut[kk, ] <- 0
                  }
                  else {
                    new.wei <- Wei[[betaN[bb]]][which(!is.na(vts.cut[kk,
                    ]))]
                    new.wei <- new.wei/sum(new.wei)
                    sub.val <- vts.cut[kk, which(!is.na(vts.cut[kk,
                    ]))] %*% new.wei
                    vts.cut[kk, which(is.na(vts.cut[kk,
                    ]))] <- sub.val
                  }
                }
              }
              dmat[((predt * (ii - 1) + 1):(predt *
                                              ii)), betaLoc] <- vts.cut %*% Wei[[betaN[bb]]]
            }
            else {
              # print(paste("node", ii, "betaN[bb]", betaN[bb]))
              # print("In length(Nei[[betaN[bb]]]) > 1 else")
              if ((length(Nei[[betaN[bb]]]) == 1) &
                  (!is.na(Nei[[betaN[bb]]]))) {
                # print("In (length(Nei[[betaN[bb]]]) == 1) & (!is.na(Nei[[betaN[bb]]]))")
                vts.cut <- vts[((maxOrder +
                                   1 - betaTimes[bb]):(predt + (maxOrder -
                                                                  betaTimes[bb]))), Nei[[betaN[bb]]]]
                #and if this is missing at any time point, set to zero
                vts.cut[is.na(vts.cut)] <- 0
                dmat[((predt * (ii - 1) + 1):(predt *
                                                ii)), betaLoc] <- vts.cut * Wei[[betaN[bb]]]
              }
              else {
                dmat[((predt * (ii - 1) + 1):(predt *
                                                ii)), betaLoc] <- 0
              }
            }
          }
        }
        else {
          for (bb in 1:sum(betaOrder)) {
            betaLoc <- which(parLoc == "b")[bb]
            dmat[((predt * (ii - 1) + 1):(predt * ii)),
                 betaLoc] <- 0
          }
        }
      }
      else {
        for (bb in 1:sum(betaOrder)) {
          betaLoc <- which(parLoc == "b")[bb]
          dmat[((predt * (ii - 1) + 1):(predt * ii)),
               betaLoc] <- 0
        }
      }
    }
  }
  if (is.null(fact.var)) {
    return(dmat)
  } else {
    #allow more than just two factors
    facun <- unique(fact.var)
    if(length(facun)==1){
      return(dmat)
    }else{
      dmcol <- ncol(dmat)
      dmatex <- dmat
      exnames <- paste(colnames(dmat), " '",facun[1],"'", sep="")
      for(ii in 2:length(facun)){ #duplicate matrix columns
        dmatex <- cbind(dmatex, dmat)
        #change names to reflect factors
        exnames <- c(exnames, paste(colnames(dmat), " '",facun[ii], "'", sep=""))
      }
      
      #for each unique factor, set other entries to 0
      
      for(ii in 1:length(facun)){
        dmatex[fact.var != facun[ii], ((ii-1)* dmcol + (1:dmcol))] <- 0
      }
      
      
      colnames(dmatex) <- exnames
      
      
      
      return(dmatex)
    }
    
  }
}
