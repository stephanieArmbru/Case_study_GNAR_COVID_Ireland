# adapt NofNeighbours to always include names of vertices



NofNeighbours_named <- function(node=1, 
                                stage=2, 
                                net=GNAR::fiveNet, 
                                numeric_vertices = FALSE){
  stopifnot(is.GNARnet(net))
  stopifnot(floor(stage)==stage)
  stopifnot(stage>0)
  stopifnot(node %in% 1:length(net$edges))
  
  
  tot.nei <- vector(mode="list", length=stage)
  tot.wei <- vector(mode="list", length=stage)
  #first get neighbours and weights
  if(!is.null(net$edges[[node]])){
    tot.nei[[1]] <- net$edges[[node]]
    tot.wei[[1]] <- net$dist[[node]]
  }
  
  #stopifnot(length(tot.nei[[1]])==length(tot.wei[[1]]))
  
  if(stage>1){
    for(ii in 2:stage){
      # only consider further stages if first stage not empty
      if(!is.null(tot.nei[[ii-1]][1])){
        # and if first stage is not NA
        if(!is.na(tot.nei[[ii-1]][1])){
          tmp.nei <- NULL
          tmp.wei <- NULL
          
          # run through all first neighbours 
          for(jj in 1:length(tot.nei[[ii-1]])){
            #get neighbours and weights
            tmp.nei <- c(tmp.nei, net$edges[[tot.nei[[ii-1]][jj]]])
            tmp.wei <- c(tmp.wei, tot.wei[[ii-1]][jj]+net$dist[[tot.nei[[ii-1]][jj]]])
            #so tmp.wei adds the distance
          }
          
          # both must have same length
          stopifnot(length(tmp.wei)==length(tmp.nei))
          
          #get rid of original node from list
          if(node %in% tmp.nei){
            posrem <- which(tmp.nei == node)
            tmp.nei <- tmp.nei[-posrem]
            tmp.wei <- tmp.wei[-posrem]
          }
          
          #get rid of nodes in previous lists
          if(sum(tmp.nei %in% unlist(tot.nei)) > 0){
            posrem <- which(tmp.nei%in%unlist(tot.nei))
            tmp.nei <- tmp.nei[-posrem]
            tmp.wei <- tmp.wei[-posrem]
          }
          
          #find minimum where have different paths to new node
          if(is.null(unique(tmp.nei))){
            tot.nei[[ii]] <- NA
          }else {
            
            if (numeric_vertices) {
              neighbourhood_names <- tmp.nei 
            }
            if (!numeric_vertices) {
              neighbourhood_names <- tmp.nei %>% names()
            }
            
            tmp.nei.df <- data.frame("CountyName" = neighbourhood_names, 
                                     "index" = tmp.nei) %>% 
              group_by(CountyName) %>% 
              summarise(index = min(index))
            
            rownames(tmp.nei.df) <- tmp.nei.df$CountyName
            
            tot.nei[[ii]] <- tmp.nei.df$index
            names(tot.nei[[ii]]) <- tmp.nei.df$CountyName
            
            # old: tot.nei[[ii]] <- unique(tmp.nei)
          }
          
          tot.wei[[ii]] <- rep(NA, length=length(tot.nei[[ii]]))
          
          if((length(tot.nei[[ii]])>0)&&(!is.na(tot.nei[[ii]][1]))){
            for(jj in 1:length(tot.nei[[ii]])){
              tot.wei[[ii]][jj] <- min(tmp.wei[tmp.nei==tot.nei[[ii]][jj]])
            }
          }else{
            tot.nei[[ii]] <- tot.wei[[ii]] <- NA
          }
        }else{
          tot.nei[[ii]] <- tot.wei[[ii]] <- NA
        }
      }else{
        tot.nei[[ii]] <- tot.wei[[ii]] <- NA
      }
    }
  }
  if(length(tot.nei[[1]])==0){
    tot.nei[[1]] <- tot.wei[[1]] <- NA
  }
  return(list(edges=tot.nei, dist=tot.wei))
}
