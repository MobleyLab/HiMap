###################################################################
## This file:
##      - runs the optimization via Fedorov.exchange
##      - computes the optimal designs with and without weighting,
##      - and outputs the optimal designs.
###################################################################

###################################################################
## Load packages, set ggplot theme, and import data
require( igraph)
require( tidyverse)
library(glue)

# Set ggplot2 theme
theme_set( theme_gray( base_size = 20))
theme_update( plot.title = element_text( hjust = 0.5), legend.position = 'top')
options( scipen = 999, stringsAsFactors = FALSE)
################################################################
# Optimization statistical engine function.
################################################################
Fedorov.exchange <- function( A.use, W.use, starting.rows, criterion, maxiter = 1e4, messages = TRUE){
  
  # Define the functions which compute the optimality criteria
  if( criterion == 'D'){ # D-optimal (minimizes determinant of the information matrix)
    crit.func <- function( A, W){
      return( det( t(A) %*% W %*% A))
    }
  }else if( criterion == 'A'){ # A-optimal (minimizes the trace of the inverse information matrix, equiavalent to minimizing the average variability of the estimates)
    crit.func <- function( A, W){
      inf.mat <- t(A) %*% W %*% A
      if( class( try( solve( inf.mat), silent = TRUE))[1] == "try-error"){
        return( -Inf)
      }else{
        return( -sum( diag( solve( inf.mat))))
      }
    }
  }else if( criterion == 'P'){ # Minimizes the average variability of all possible pairwise differences
    crit.func <- function( A, W){
      inf.mat <- t(A) %*% W %*% A
      if( class( try( solve( inf.mat), silent = TRUE))[1] == "try-error"){
        return( -Inf)
      }else{
        inv.inf <- solve( inf.mat)
        all.pairs.var <- NULL
        for( i in 1:(nrow( inv.inf) - 1)){
          for( j in (i + 1):nrow( inv.inf)){
            all.pairs.var <- c( all.pairs.var, ( inv.inf[i,i] + inv.inf[j,j] - 2*inv.inf[i,j]))
          }
        }
        return( -sum( all.pairs.var))
      }
    }
  }else if( criterion == 'mA'){ # Minimizes the maximum variability of the estimates
    crit.func <- function( A, W){
      inf.mat <- t(A) %*% W %*% A
      if( class( try( solve( inf.mat), silent = TRUE))[1] == "try-error"){
        return( -Inf)
      }else{
        return( -max( diag( solve( inf.mat))))
      }
    }
  }else if( criterion == 'mP'){ # Minimizes the maximum variability of the pairwise differences
    crit.func <- function( A, W){
      inf.mat <- t(A) %*% W %*% A
      if( class( try( solve( inf.mat), silent = TRUE))[1] == "try-error"){
        return( -Inf)
      }else{
        inv.inf <- solve( inf.mat)
        all.pairs.var <- NULL
        for( i in 1:(nrow( inv.inf) - 1)){
          for( j in (i + 1):nrow( inv.inf)){
            all.pairs.var <- c( all.pairs.var, (inv.inf[i,i] + inv.inf[j,j] - 2*inv.inf[i,j]))
          }
        }
        return( -max( all.pairs.var))
      }
    }
  }else if( criterion == 'negD'){ # Finds the worst design for the D-optimal criteria
    crit.func <- function( A, W){
      inf.mat <- t(A) %*% W %*% A
      if( class( try( solve( inf.mat), silent = TRUE))[1] == "try-error"){
        return( -Inf)
      }else{
        return( -det( t(A) %*% W %*% A))
      }
    }
  }else if( criterion == 'negA'){ # Finds the worst design for the A-optimal criteria
    crit.func <- function( A, W){
      inf.mat <- t(A) %*% W %*% A
      if( class( try( solve( inf.mat), silent = TRUE))[1] == "try-error"){
        return( -Inf)
      }else{
        return( sum( diag( solve( inf.mat))))
      }
    }
  }else if( criterion == 'random'){ # Finds a random design that's also identifiable
    n.pairs <- length( starting.rows[ -which( starting.rows %in% which( rowSums( abs( A.use)) == 1))])
    rand.iter <- 1
    repeat{
      current.rows <- c( sample( (1:nrow( A.use))[-which( rowSums( abs( A.use)) == 1)], n.pairs), which( rowSums( abs( A.use)) == 1))
      inf.mat <- t( A.use[ current.rows,]) %*% W.use[ current.rows, current.rows] %*% A.use[ current.rows,]
      print(paste0( 'Iteration ', rand.iter, ' to find seed design.'))
      rand.iter <- rand.iter + 1
      # Edits made here 04/11, test log(det(inf.mat)) instead of det(inf.mat)
      crit_i <- log(det( inf.mat))
      #print(det( inf.mat))
      if( !isTRUE( all.equal(crit_i, -Inf))){
        print(paste0("Finite log(det(inf.mat)) found, ", crit_i))
        if( class( try( solve( inf.mat), silent = TRUE))[1] != "try-error"){
          print("Invertible solution found.")
          break
        }
      }
    }
    return(
      list(
        A = A.use[ current.rows,],
        W = W.use[ current.rows, current.rows],
        criterion = NA,
        rows = current.rows
      )
    )
  }else{
    warning( 'Unknown Criterion')
    return(NULL)
  }
  
  # Controls the optimization iterations, below.
  # Initialize the design matrix
  current.rows <- starting.rows
  # Keep track of iterations
  current.iter <- 1
  # Run the algorithm
  repeat{
    # Compute the criterion for the current design
    current.D <- crit.func( A.use[ current.rows,], W.use[ current.rows, current.rows])
    
    if( messages){
      message( paste0( 'Iteration ', current.iter, ' (criterion = ', signif( current.D, 4), ')'))
    }
    # Is something going wrong with the messages?
    print(paste0( 'Iteration ', current.iter, ' (criterion = ', signif( current.D, 4), ')'))
    
    # Find all possible couples of current rows and potential rows
    candidate.rows <- setdiff( 1:nrow( A.use), current.rows)
    couples <- expand.grid( current = current.rows, candidate = candidate.rows, D = NA)
    
    # Compute the difference in criteria between all possible exchanges
    for( current.couple in 1:nrow( couples)){
      current.design <- A.use[ c( current.rows[ -which( current.rows == couples$current[ current.couple])], couples$candidate[ current.couple]),]
      current.weights <- W.use[ c( current.rows[ -which( current.rows == couples$current[ current.couple])], couples$candidate[ current.couple]), c(current.rows[-which(current.rows == couples$current[current.couple])], couples$candidate[current.couple])]
      couples$D[ current.couple] <- crit.func( current.design, current.weights)
    }
    couples$D.diff <- couples$D - current.D
    
    # If no exchanges improve the design, stop.  Otherwise, make the best swap and continue the algorithm
    if( sum( couples$D.diff > 0) == 0){
      break
    } else if( current.iter == maxiter) {
      # This might be a bug. What I want it to do is give up. Not return garbage below.
      print("The maximum number of iterations occured.")
      break
    }else{
      current.rows <- c( current.rows[ -which( current.rows == couples$current[ which.max( couples$D.diff)[1]])], couples$candidate[ which.max( couples$D.diff)[1]])
      current.iter <- current.iter + 1
    }
  }
  
  # Return the optimal design, weight matrix, criteria, and rows
  return(
    list(
      A = A.use[ current.rows,],
      W = W.use[ current.rows,current.rows],
      criterion = current.D,
      rows = current.rows
    )
  )
}
################################################################
# End optimization statistical engine.
################################################################

################################################################
# Import scores from python
################################################################
run_optimization <- function(ref_lig, dataframe, r_optim_types, num_edges){
    # ref_lig = reference ligand to use, currenlty must be manually selected
    # dataframe  = r dataframe generate from similarity scores
    # optim_type1 and optim_type2 = optimization type such as 'A' and 'D'.
    #                               currently two types are required.
    
    # Import atom pair scores
    optim_type1 <- r_optim_types[1]
    optim_type2 <- r_optim_types[2]
    k_numb <- num_edges
    print(paste("Preparing", optim_type1, "optimization"))
    print(paste("Preparing", optim_type2, "optimization"))

    sim.scores <- dataframe %>%
      as.matrix() %>%
      reshape2::melt( as.is = TRUE) %>%
      as_tibble() %>%
      setNames( c( 'LIGAND_1', 'LIGAND_2', 'SIM')) %>%
      filter( LIGAND_1 != LIGAND_2) %>%
      # Normalize weights.
      mutate( WEIGHT = (SIM^2 - min(SIM)^2)/(max(SIM)^2 - min(SIM)^2))
    
    ## Find the optimal designs
    # Define the list of ligands
    ligand.list <- sim.scores$LIGAND_1 %>% unique() %>% sort()
    reference.ligand <- c(ref_lig)
 
    # Generate full design
    A.full <- matrix( 0, ncol = length( ligand.list), nrow = choose( length( ligand.list), 2))
    for( i in 1:nrow( A.full)){
      A.full[i, t( combn( 1:length( ligand.list), 2))[i,]] <- c( 1, -1)
    }
    print( paste("Number of ligands:", length( ligand.list)))
    print(paste('The number of chosen edges is', k_numb))
    ###############

    # Define the A and W matrices for a full design with weights
    # Output array with a 1 at column of reference ligands. Rows: number of ref ligands
    a <- rbind(diag( ncol( A.full))[ which( ligand.list %in% reference.ligand),])
    A.use <- rbind( A.full, a)

    # Generate array of 2's for length of reference ligands
    ref_weights <- rep(2, length(reference.ligand))
    W.use <- diag( c(
      apply( A.full, 1, function(x){
        sim.scores %>%
          filter(
            LIGAND_1 == ligand.list[ which( x == 1)] &
              LIGAND_2 == ligand.list[ which( x == -1)]
          ) %>%
          pull( WEIGHT)
      }), ref_weights
    ))

    # Generate seed design to optimize.
    starting.rows <- Fedorov.exchange( A.use, W.use, c( 1:k_numb, nrow(A.use)), 'random')$rows
    # Begin optimization.
    a.opt.new <- Fedorov.exchange( A.use, W.use, starting.rows, optim_type1)$A %>%
      {.[which( rowSums( .) == 0),]} %>%
      apply( 1, function(x){
        tibble(
          LIGAND_1 = ligand.list[ which( x == 1)],
          LIGAND_2 = ligand.list[ which( x == -1)]
        )
      }) %>% {do.call( "rbind", .)}
    d.opt.new <- Fedorov.exchange( A.use, W.use, starting.rows, optim_type2)$A %>%
    {.[which( rowSums( .) == 0),]} %>%
      apply( 1, function(x){
        tibble(
          LIGAND_1 = ligand.list[ which( x == 1)],
          LIGAND_2 = ligand.list[ which( x == -1)]
        )
      }) %>% {do.call( "rbind", .)}

    ###################################################################
    ## Combine all designs and weights together
    string_A = paste(optim_type1, "-optimal")
    string_D = paste(optim_type2, "-optimal")
        
    analysis.dat <- bind_rows(
      a.opt.new %>% mutate( DESIGN = string_A),
      d.opt.new %>% mutate( DESIGN = string_D)
    ) %>%
      left_join( sim.scores %>% select( LIGAND_1, LIGAND_2, 'SIM_WEIGHT' = WEIGHT))
      
    ###################################################################
    ## Visualize designs

    # Compute the optimality criteria for each design.
    crit.dat <- analysis.dat %>%
      tidyr::crossing( as_tibble( t( setNames( rep( 0, length( ligand.list)), ligand.list)))) %>%
      rowwise() %>%
      do( mutate_if( as_tibble(.), grepl( .$LIGAND_1, names(.)), function(x){-1})) %>%
      do( mutate_if( as_tibble(.), grepl( .$LIGAND_2, names(.)), function(x){1})) %>%
      ungroup() %>%
      group_by( DESIGN) %>%
      do({
        design.mat <- {.} %>%
          select_at( ligand.list) %>%
          as.matrix() %>%
          rbind( diag( length( ligand.list))[ which( ligand.list == reference.ligand),])
        
        null.weight <- c( rep( 1, nrow(.)), 2)
        sim.weight <- c( .$SIM_WEIGHT, 2)
        
        # Removed solve from null.inner because det(inverse A) = 1/det(A)
        null.inner =  t( design.mat) %*% diag( null.weight) %*% design.mat
        weighted.inner =  t( design.mat) %*% diag( sim.weight) %*% design.mat
        
        # This is the hat matrix values, in development.
        # H diags
        #print("H diags")
        #print(diag(((design.mat %*% solve(t(design.mat) %*% diag( sim.weight) %*% design.mat) %*% t(design.mat) %*% diag( sim.weight)))))
        
        # Output critical data dataframe
        data.frame(
          A = sum( diag( solve( t( design.mat) %*% diag( null.weight) %*% design.mat))),
          D = det( solve( t( design.mat) %*% diag( null.weight) %*% design.mat))^(1/length( ligand.list)),
          A.ap = sum( diag( solve( t( design.mat) %*% diag( sim.weight) %*% design.mat))),
          D.ap = det( solve( t( design.mat) %*% diag( sim.weight) %*% design.mat))^(1/length( ligand.list))
          #H =  ((design.mat %*% solve(t(design.mat) %*% diag( sim.weight) %*% design.mat) %*% t(design.mat) %*% diag( sim.weight))),
          #Diag.Hi =  diag(((design.mat %*% solve(t(design.mat) %*% diag( sim.weight) %*% design.mat) %*% t(design.mat) %*% diag( sim.weight)))),
          #Diag.nullHi =  diag(((design.mat %*% solve(t(design.mat) %*% design.mat) %*% t(design.mat))))
          )
    
      }) %>%
      ungroup()
    
    # For in development.
    #print( "crit.dat, leverages")
    #print( crit.dat[, c("Diag.Hi", "Diag.nullHi")])
    
    print( "Critical Data:")
    print( crit.dat)
  
    # Plot the designs.
    # Collect vertices or node info for designs.
    vertex.dat <- analysis.dat %>%
      group_by( DESIGN) %>%
      do({
        net <- graph_from_data_frame(., directed = FALSE)
        bind_cols(
          as.data.frame( vertex.attributes( net)),
          as.data.frame( layout_with_gem( net))
        ) %>% setNames( c( 'LIGAND', 'X', 'Y'))
      }) %>% ungroup() %>%
      mutate( REF = ifelse( LIGAND %in% reference.ligand, 'Reference Ligand', 'Other')) %>%
      mutate( REF = factor( REF, levels = c( 'Reference Ligand', 'Other')))
    
    print("Vertex Data:")
    vertex.dat %>% as_tibble() %>% print(n=60)
    
    # Collect edge data for the designs.
    edge.dat <- analysis.dat %>%
      left_join( vertex.dat %>% select( 'LIGAND_1' = LIGAND, DESIGN, 'X1' = X, 'Y1' = Y)) %>%
      left_join( vertex.dat %>% select( 'LIGAND_2' = LIGAND, DESIGN, 'X2' = X, 'Y2' = Y))

    # Change name based on metric type. Might need to deprecate.
    vertex.dat <- vertex.dat %>%
      mutate( DESIGN = factor( DESIGN, levels = c( string_A,  string_D)))

    edge.dat <- edge.dat %>%
      mutate( DESIGN = factor( DESIGN, levels = c( string_A,  string_D)))
    
   # For in development.
    #leverage.dat <- crit.dat %>%
      #mutate( DESIGN = factor( DESIGN, levels = c( "A-optimal (lomap)",  "D-optimal (lomap)")))
        
    # Plotting of the weights colored.
    bottom.plot <- ggplot( vertex.dat, aes( x = X, y = Y)) +
      geom_segment( size = 1, data = edge.dat, aes( x = X1, y = Y1, xend = X2, yend = Y2, color = SIM_WEIGHT)) +
      geom_point( size = 10, pch = 21, aes( fill = REF)) +
      geom_text( size = 2, aes( x = X, y = Y, label = LIGAND)) +
      facet_wrap( ~ DESIGN, scales = "free", ncol = 2) +
      scale_color_gradientn( colors = c( "red1", "orange", "gold", "green", "springgreen"), limits = c( 0, 1)) +
      theme_bw( base_size = 20) + labs( fill = NULL, color = NULL) +
      theme( panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_blank()) +
      theme( axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank()) +
      theme( plot.title = element_text( hjust = 0.5), legend.position = 'top', legend.box = 'vertical', legend.key.width = unit(4, 'line'), legend.title=element_text(vjust = 0.9)) +
      labs( color = glue('Weights:'))
    now <- Sys.time()
    ggsave(paste0("Rgraph_colored", format(now, "_%Y%m%d_%H%M%S"), ".pdf"))
    
    # Write csv file
    csv_ofile = paste0("edge_data", format(now, "_%Y%m%d_%H%M%S"), ".csv")
    write.csv(edge.dat,csv_ofile, row.names = FALSE)
    
    #return(edge.dat)
}
