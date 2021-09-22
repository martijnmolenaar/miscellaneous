options(stringsAsFactors = FALSE)

library(dplyr)
library(pbmcapply)

## fils NaNs with zeros
fill.na <- function(x, fill = 0) {
  x2 <- x
  x2[is.nan(x2)] <- fill
  x2
}

## fils Inf with zeros
fill.inf <- function(x, fill = NA){
  x2 <- x
  x2[is.infinite(x2)] <- fill
  x2
}

## number of neighbours of vector within range
n_neighbours <- function(x, range = .1){
  max_dist <- diff(range(x)) * range
  sapply(x, function(x_n){
    sum(abs(x_n - x) < max_dist)
  })
}

## omit both vector indices of two vectors when at least one of them are NA 
omitNAs <- function(x,y){
  keep <- !is.na(x) & !is.na(y)
  list(x=x[keep], y=y[keep])
}

## simple model to simulate ion suppression
simulate_ion_suppression <- function(signal, max_suppression = .5){
  ion_s <- (-(signal / max(signal) * (1-max_suppression)) +1)
  signal * ion_s
}


### calculate_proportion and specificity matrix

# for instance, sample_location = '2021_03_03/W3/'

calculate_proportion_matrix <- function(sample_location){
  library(tidyr)
  
  overlap <- read.csv(file = paste0(sample_location, 'analysis/overlap_analysis2/overlap.regions.csv'))
  overlap$cell_id <- paste0("cell#",overlap$cell_id)
  overlap$am_id <- paste0("AM#",overlap$am_id)
  
  AM <- read.csv(file = paste0(sample_location, 'analysis/overlap_analysis2/ablation_mark.regions.csv'))
  AM$am_id <- paste0("AM#",AM$am_id)
  
  cell <- read.csv(file = paste0(sample_location, 'analysis/overlap_analysis2/cell.regions.csv'))
  cell$cell_id <- paste0("cell#",cell$cell_id)
  
  ### cast long overlap data.frame into wide data.frame
  
  overlap$cell_id <- factor(overlap$cell_id, levels = cell$cell_id)      
  overlap$am_id <- factor(overlap$am_id, levels = paste0("AM#",1:10000))
  
  overlap_matrix <- overlap[,colnames(overlap) %in% c("cell_id","am_id","area")] %>% 
    spread(key = am_id, value = area, drop = F, fill = 0)
  
  rownames(overlap_matrix) <- overlap_matrix$cell_id
  overlap_matrix <- overlap_matrix[,-1]    ## remove cell IDs to get numerical matrix
  
  ### calculate sampling proportions (= sampling ratios) from areas
  
  overlap_proportion_matrix <- overlap_matrix
  overlap_proportion_matrix <- rbind(overlap_proportion_matrix, total = AM$area)  ## add background
  
  overlap_proportion_matrix[] <-                                                       ## calculate proportions
    apply(overlap_proportion_matrix, 2, function(cell_i){
      cell_i / cell_i[length(cell_i)]
    })
  
  
  overlap_proportion_matrix <- overlap_proportion_matrix[rownames(overlap_proportion_matrix) != "total",]  ## remove background
  
  
  
  ### calculate sampling specificities
  overlap_specificity_matrix <-
    as.data.frame(t(apply(overlap_proportion_matrix, 1, function(cell_i) {
      specificities <- cell_i / sum(cell_i)
      specificities[is.nan(specificities)] <- 0
      return(specificities)
    })))
  
  list(overlap_proportion_matrix = overlap_proportion_matrix[rowSums(overlap_proportion_matrix) != 0,],
       overlap_specificity_matrix = overlap_specificity_matrix[rowSums(overlap_proportion_matrix) != 0,])
}


#### cell normalization definition by Rappez et al, Nature Methods
## faster version w/ matrices, is not completely the same as #1, cannot find the difference
## works now with feature matrix input

cell_normalization_Rappez_et_al <-
  function(overlap_proportion_matrix,
           overlap_specificity_matrix,
           AM_ion_intensity,
           AM_proportion_treshold_global = .3,
           AM_specificity_treshold_cell = 0, 
           AM_area_percentile = c(0,1),             ## not implemented yet
           AM_area = NULL,                          ## not implemented yet
           skip_AM_zeros = FALSE) {
    
    #browser()
    AM_proportions_global <- colSums(overlap_proportion_matrix)
    AM_proportions_filter <-
      AM_proportions_global > AM_proportion_treshold_global
    AM_specificity_filter <-
      overlap_specificity_matrix > AM_specificity_treshold_cell
    
    AM_proportions_global_filered <-      AM_proportions_global * AM_proportions_filter
    AM_specificity_filtered <-    as.matrix(overlap_specificity_matrix) * AM_specificity_filter
    
    norm_ion_int_cell <- 
      apply(as.matrix(AM_ion_intensity), 2, function(AM_ion_intensity_i){               ## iterate over features
        
        if(skip_AM_zeros) {
          AM_ion_intensity_j <- AM_ion_intensity_i
          AM_ion_intensity_j[is.na(AM_ion_intensity_j)] <- 0
          non_zeros <- AM_ion_intensity_j != 0
        } else {
          non_zeros <- T
        }
        
        #fill.na(
        cell_intensities_i <- rowSums(
            t( fill.inf(AM_ion_intensity_i[non_zeros] / AM_proportions_global_filered[non_zeros]) *   t(AM_specificity_filtered[,non_zeros]))
            , na.rm = T)   /     
            rowSums(AM_specificity_filtered[,non_zeros])
          #)
          
          has_noAM_cells <- rowSums(overlap_proportion_matrix[,AM_proportions_filter & non_zeros] != 0) == 0
          cell_intensities_i[has_noAM_cells] <- NA
          
          return(cell_intensities_i)
      })
    
    
    return(list(cell_intensities = norm_ion_int_cell))
    
    
  }


#### cell normalization definition by NNLS (Non-negative Least Squares), lineair inverse modeling
#### somewhat slow

cell_normalization_NNLS <-
  function(overlap_proportion_matrix,
           AM_ion_intensity,
           AM_proportion_treshold_per_cell = .3,
           AM_area_percentile = c(0,1),             ## not implemented yet
           AM_area = NULL,                          ## not implemented yet
           skip_AM_zeros = FALSE) {
    
    library(nnls)
    #browser()
    B <- as.matrix(AM_ion_intensity)
    
    A <- overlap_proportion_matrix
    A <- t(as.matrix(rbind(A, 1 - colSums(A))))   ## last row is background
    colnames(A)[dim(A)[2]] <- "background"
    
    AM_filter <- apply(A[,-dim(A)[2]], 1, function(Ai){any(Ai > AM_proportion_treshold_per_cell)})     ### at least .. overlap
    
    cell_intensities <-
      apply(B,2,function(B_i){
        
        if(skip_AM_zeros) {
          B_j <- B_i
          B_j[is.na(B_j)] <- 0
          non_zeros <- B_j != 0
          
          
        } else {
          non_zeros <- T
        }
        
        has_noAM_cells <- colSums(A[AM_filter & non_zeros,] != 0) == 0
        
        cell_intensities <- nnls(A = A[AM_filter & non_zeros, ],
             b = B_i[AM_filter & non_zeros])$x
        
        
        
        cell_intensities[has_noAM_cells] <- NA   ## return NA when 
        
        return(cell_intensities)
      })
    
    
    
    cell_intensities <-                               ## last entry is no cell but background
      cell_intensities[-dim(cell_intensities)[1],]
    
    names(cell_intensities) <- rownames(overlap_proportion_matrix)
    
    return(list(cell_intensities = cell_intensities))
    
  }

#### feature matrix for normalization quality assessment

calculate_normalization_features <-
  function(overlap_proportion_matrix,
           AM_ion_intensity,
           AM_proportion_treshold_global = .3,
           norm_ion_int_cell,
           AM_area_percentile = c(0,1),             ## not implemented yet
           AM_area = NULL,                          ## not implemented yet
           skip_AM_zeros = FALSE){
    
    library(igraph)
    
    AM_proportions_global <- colSums(overlap_proportion_matrix)
    AM_proportions_filter <-
      AM_proportions_global > AM_proportion_treshold_global
    
    AM_proportions_global_filered <-      AM_proportions_global * AM_proportions_filter
    
    if(skip_AM_zeros) {
      AM_ion_intensity_j <- AM_ion_intensity_i
      AM_ion_intensity_j[is.na(AM_ion_intensity_j)] <- 0
      non_zeros <- AM_ion_intensity_j != 0
    } else {
      non_zeros <- T
    }
    
    ### feature matrix
    
    AM_residuals <- as.matrix(t(overlap_proportion_matrix[,AM_proportions_filter & non_zeros])) %*% as.matrix(norm_ion_int_cell) - as.matrix(AM_ion_intensity[AM_proportions_filter & non_zeros])
    AM_residuals <- AM_residuals / as.matrix(AM_ion_intensity[AM_proportions_filter & non_zeros])
    
    cell_residuals <-
      sapply(seq_along(overlap_proportion_matrix[,AM_proportions_filter & non_zeros][,1]), function(cell_i){
        AM_boolean <- overlap_proportion_matrix[cell_i,AM_proportions_filter & non_zeros] > 0
        sum(AM_residuals[AM_boolean])
      })
    
    cell_co_AM_level <-  
      sapply(seq_along(overlap_proportion_matrix[,AM_proportions_filter & non_zeros][,1]), function(cell_i){
        AM_boolean <- overlap_proportion_matrix[cell_i,AM_proportions_filter & non_zeros] > 0
        fill.na(mean(colSums(
          as.matrix(overlap_proportion_matrix[,AM_proportions_filter & non_zeros][,AM_boolean])
          >0)))
      })
    
    cell_nr_of_co_AM <-  
      sapply(seq_along(overlap_proportion_matrix[,AM_proportions_filter & non_zeros][,1]), function(cell_i){
        AM_boolean <- overlap_proportion_matrix[cell_i,AM_proportions_filter & non_zeros] > 0
        sum(colSums(
          as.matrix(overlap_proportion_matrix[,AM_proportions_filter & non_zeros][,AM_boolean])
          >0) > 1)
      })
    
    cell_nr_of_non_co_AM <-  
      sapply(seq_along(overlap_proportion_matrix[,AM_proportions_filter & non_zeros][,1]), function(cell_i){
        AM_boolean <- overlap_proportion_matrix[cell_i,AM_proportions_filter & non_zeros] > 0
        sum(colSums(
          as.matrix(overlap_proportion_matrix[,AM_proportions_filter & non_zeros][,AM_boolean])
          >0) < 2)
      })
    
    cell_nr_of_AM <-  
      sapply(seq_along(overlap_proportion_matrix[,AM_proportions_filter & non_zeros][,1]), function(cell_i){
        sum(overlap_proportion_matrix[cell_i,AM_proportions_filter & non_zeros] > 0)
        
      })
    
    compensated_AM_ion_intensity <-  AM_ion_intensity / AM_proportions_global_filered
    
    cell_CV <-  
      sapply(seq_along(overlap_proportion_matrix[,AM_proportions_filter & non_zeros][,1]), function(cell_i){
        compensated_AM_ion_intensity_per_cell <- compensated_AM_ion_intensity[AM_proportions_filter & non_zeros][    overlap_proportion_matrix[cell_i,AM_proportions_filter & non_zeros] > 0 ]
        if(length(compensated_AM_ion_intensity_per_cell) > 1){
          sd(compensated_AM_ion_intensity_per_cell, na.rm = T) / mean(compensated_AM_ion_intensity_per_cell, na.rm = T)
        } else {
          0
        }
        
      })
    
    
    AM_links <- apply(overlap_proportion_matrix[,AM_proportions_filter & non_zeros],2,function(AM_i){
      output <- which(AM_i > 0)
      if(length(output > 0)){
        paste0("cell#",output)
      } else {
        output
      }
    })
    names(AM_links) <- paste0("AM#",seq_along(AM_links))
    vertices <- names(AM_links)
    AM_links <- AM_links[sapply(AM_links, length) > 0]
    
    vertices <- c(vertices, paste0("cell#",1:dim(overlap_proportion_matrix[,AM_proportions_filter & non_zeros])[1]))
    
    #cell_links <- cell_links[sapply(cell_links, length) > 0]
    
    edges <- 
      sapply(names(AM_links), function(i) {
        data.frame(from = i ,
                   to = AM_links[[i]])
      }, simplify = F) %>% bind_rows()
    
    v <- data.frame(id = vertices)
    
    v$shape <- ifelse(grepl("AM",v$id),"circle","square")
    v$color <- ifelse(grepl("AM",v$id),"gray30","lightpink")
    v$size <- ifelse(grepl("AM",v$id),3,10)
    
    G <- graph_from_data_frame(d=edges,  directed=F, vertices = v) 
    
    a <-
      data.frame(node_ID = V(G)$name, 
                 community = components(G)$membership) %>% group_by(community) %>%
      mutate(AMs_within_comm = sum(grepl("AM", node_ID)),
             cells_within_comm = sum(grepl("cell", node_ID)),
             abs_dermination = AMs_within_comm-cells_within_comm,
             relative_determination = AMs_within_comm / cells_within_comm) %>% filter(grepl("cell", node_ID))
    
    b <-
      data.frame(id = paste0("cell#",seq_along(overlap_proportion_matrix[,AM_proportions_filter & non_zeros][,1])),
                 cell_residuals, cell_nr_of_co_AM, cell_nr_of_non_co_AM, cell_co_AM_level, cell_nr_of_AM, cell_CV,
                 a)
    
    
    
    return(list(feature_matrix = b,
                graph = G))
    
  }


################# some retired functions
#### Weighted average as interpreted by Martijn
#### super slow because of iteration over cell (but direct matrix application is harder to understand)

cell_normalization_weighted_average <-
  function(overlap_proportion_matrix,
           AM_ion_intensity,
           AM_proportion_treshold_per_cell = .3) {
    
    cell_intensities <- unlist(pbmclapply(1:dim(overlap_proportion_matrix)[1], function(cell_i) {
      AM_proportions_by_cell <- t(overlap_proportion_matrix[cell_i, ])[, 1]
      
      AM_probs_filtered <-
        AM_proportions_by_cell[AM_proportions_by_cell > AM_proportion_treshold_per_cell]
      AM_ion_intensity_filtered <-
        AM_ion_intensity[AM_proportions_by_cell > AM_proportion_treshold_per_cell]
      
      norm_ion_int_cell <-
        weighted.mean(x = AM_ion_intensity_filtered / AM_probs_filtered, w = AM_probs_filtered)
      return(norm_ion_int_cell)
    }))
    
    cell_intensities[is.nan(cell_intensities)] <- 0
    
    return(cell_intensities)
  }

#### homogeneous_am_and_cell_by_sampling_ratio as interpreted by Martijn
#### super slow because of iteration over cell (but direct matrix application is harder to understand)

cell_normalization_homogeneous_am_and_cell_by_sampling_ratio <-
  function(overlap_proportion_matrix,
           AM_ion_intensity,
           AM_proportion_treshold_global = .3) {
    
    AM_proportions_global <- colSums(overlap_proportion_matrix)
    
    AM_filter <- (AM_proportions_global > AM_proportion_treshold_global) 
    
    AM_proportions_global_filered <- AM_proportions_global[AM_filter]
    
    AM_ion_intensity_filered <-   AM_ion_intensity[AM_filter]
    
    overlap_proportion_matrix_filtered <- overlap_proportion_matrix[,AM_filter]
    
    cell_intensities <- unlist(pbmclapply(1:dim(overlap_proportion_matrix_filtered)[1], function(cell_i) {
      norm_ion_int_cell <-
        sum((AM_ion_intensity_filered) * t(overlap_proportion_matrix_filtered[cell_i, ])) /   
        sum(t(overlap_proportion_matrix_filtered[cell_i, ]))
      
      return(norm_ion_int_cell)
    }))
    
    cell_intensities[is.nan(cell_intensities)] <- 0
    
    return(cell_intensities)
    
  }



#### cell normalization definition by Rappez et al, Nature Methods
#### super slow because of iteration over cell (but direct matrix application is harder to understand)

cell_normalization_Rappez_et_al_old <-
  function(overlap_proportion_matrix,
           overlap_specificity_matrix,
           AM_ion_intensity,
           AM_proportion_treshold_global = .3,
           AM_specificity_treshold_cell = 0) {
    
    cell_intensities <- unlist(pbmclapply(1:dim(overlap_proportion_matrix)[1], function(cell_i) {

      AM_proportions_global <- colSums(overlap_proportion_matrix)
      AM_filter <- (AM_proportions_global > AM_proportion_treshold_global) &
        (overlap_specificity_matrix[cell_i,] > AM_specificity_treshold_cell)
      
      AM_proportions_global_filered <-
        AM_proportions_global[AM_filter]
      
      AM_specificity_by_cell <-
        overlap_specificity_matrix[cell_i, ]
      
      AM_specificity_by_cell_filtered <-
        AM_specificity_by_cell[AM_filter]
      
      AM_ion_intensity_filered <-
        AM_ion_intensity[AM_filter]
      
      norm_ion_int_cell <-
        sum((AM_ion_intensity_filered / AM_proportions_global_filered) * AM_specificity_by_cell_filtered
        ) /     sum(AM_specificity_by_cell_filtered)
      return(norm_ion_int_cell)
    }))
    
    cell_intensities[is.nan(cell_intensities)] <- 0
    
    return(cell_intensities)
    
  }


### weighted_by_overlap_and_sampling_area, MM implementation

cell_normalization_SpaceM_weighted_by_overlap_and_sampling_area <-
  function(overlap_proportion_matrix,
           AM_areas,
           AM_ion_intensity,
           AM_proportion_treshold_global = .3) {
    
    library(geometry)
    
    am_intracellular_sampling_areas <- colSums(overlap_proportion_matrix)  #np.sum(overlap_matrix, axis=0)
    am_extracellular_penalty <- am_intracellular_sampling_areas / AM_areas
    AM_filter <- am_extracellular_penalty > AM_proportion_treshold_global
    
    sample_of_interest_weight <- as.matrix(overlap_proportion_matrix) / am_intracellular_sampling_areas
    
    coefficient_matrix <- t(
      t(am_extracellular_penalty[AM_filter] * sample_of_interest_weight[, AM_filter])
      / rowSums(sample_of_interest_weight, na.rm = T)
    )

    coefficient_matrix[is.nan(coefficient_matrix)] <- 0
    
    cell_intensities <-      t(t(AM_ion_intensity)[AM_filter] %*% t(coefficient_matrix))
    
    return(cell_intensities)
    
  }




