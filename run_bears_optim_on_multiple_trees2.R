run_bears_optim_on_multiple_trees2 <- function(BioGeoBEARS_run_object, newick_fns, model_name="", geog_fns=NULL, resfns=NULL, run_or_collect="collect", start_treenum=1, end_treenum=length(newick_fns), runslow=TRUE, plot_params=FALSE)
{
  defaults='
  BioGeoBEARS_run_object=BioGeoBEARS_run_object
  newick_fns=newick_fns
  model_name=model_name
  geog_fns=NULL
  resfns=NULL
  start_treenum=1
  end_treenum=2
  run_or_collect="both"
  runslow=runslow
  plot_params=FALSE
  '
  
    
  require(stringr)
  
  # Check the input newick_fns to make sure they end in .newick
  num_newick_strings = str_count(string=newick_fns, pattern="\\.newick")
  num_newick_strings_equals_1_TF = num_newick_strings == 1
  if (sum(num_newick_strings_equals_1_TF) != length(num_newick_strings_equals_1_TF))
  {
    error_txt = paste("\n\nERROR in run_bears_optim_on_multiple_trees(): All filenames in 'newick_fns' must have one and only one '.newick'.\nViolators in your 'newick_fns':\n\n", sep="")
    cat(error_txt)
    cat(newick_fns[num_newick_strings_equals_1_TF=FALSE], sep="\n")
    
    stop("Stopping on error in run_bears_optim_on_multiple_trees()")
  }
  
  # File names to store the state probabilities at nodes and corners
  if (model_name == "")
  {
    suffix_txt = ""
  } else {
    suffix_txt = paste("_", model_name, sep="")
  }
  stateprobs_at_nodes_fn = paste("state_probs_at_nodes_across_all_trees", suffix_txt, ".Rdata", sep="")
  stateprobs_at_corners_fn = paste("state_probs_at_corners_across_all_trees", suffix_txt, ".Rdata", sep="")
  
  # Get names of:
  # - geography filenames (geog_fns)
  # - results filenames (resfns)
  #newick_fns = output_trfns
  if (is.null(geog_fns))
  {
    geog_fns = gsub(pattern="newick", replacement="geog", x=newick_fns)
    cat(geog_fns, sep=", ")
  }
  
  # Error check for repeating a single geog file
  if (length(geog_fns) == 1)
  {
    geog_fns = rep(geog_fns, times=length(newick_fns))
  }
  
  if (is.null(resfns))
  {
    replacement_txt = paste(suffix_txt, ".Rdata", sep="")
    resfns = gsub(pattern="\\.newick", replacement=replacement_txt, x=newick_fns)
  }
  #print(runslow)
  if (runslow == TRUE)
  {
    row_start = 1
    row_end = 0
    for (i in start_treenum:end_treenum)
      #for (i in 1:length(newick_fns))
      #for (i in 1:2)
    {
      # If you just want to run the inferences and save the results
      if ((run_or_collect == "run") || (run_or_collect == "both"))
      {
        txt = paste("\n\nRunning inference under '", model_name, "' for tree #", i, " of ", start_treenum, "-", end_treenum, "...\nTree file: ", newick_fns[i], "\nGeog file: ", geog_fns[i], "\n\n", sep="")
        cat(txt)
        
        BioGeoBEARS_run_object$geogfn = geog_fns[i]
        BioGeoBEARS_run_object$trfn = newick_fns[i]
        resfn = resfns[i]
        
        
        # Make a tree table containing the node ages
        tree_table = prt(read.tree(newick_fns[[i]]))
        # The root age is the maximum node age:
        root_age = max(tree_table$time_bp)
        # Let's say the last time bin starts at 10 mya
        if (root_age > 5.35)
        {
          BioGeoBEARS_run_object$timesfn = "timeperiods2.txt"
          BioGeoBEARS_run_object$dispersal_multipliers_fn = "manual_dispersal_multipliers2.txt"
        } else {
          if (root_age > 4.85)
          {
            BioGeoBEARS_run_object$timesfn = "timeperiods3.txt"
            BioGeoBEARS_run_object$dispersal_multipliers_fn = "manual_dispersal_multipliers3.txt"
          } else {
            BioGeoBEARS_run_object$timesfn = "timeperiods4.txt"
            BioGeoBEARS_run_object$dispersal_multipliers_fn = "manual_dispersal_multipliers4.txt"
          }
        }
        
        
        
        
        # Check the max number of areas:
        #tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=BioGeoBEARS_run_object$geogfn)
        # max(rowSums(dfnums_to_numeric(tipranges@df), na.rm=TRUE))
        
        # You have to re-read the files, now that you've changed them!
        BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
        
        # And re-stratify the tree, if this is stratified
        if (length(BioGeoBEARS_run_object$timeperiods) > 1)
        {
          BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE, cut_fossils=FALSE)
        } # END if (length(BioGeoBEARS_run_object$timeperiods) > 1)
        
        # And, re-run the check
        check_BioGeoBEARS_run(BioGeoBEARS_run_object)
        
        res = bears_optim_run(BioGeoBEARS_run_object)
        res    
        
        save(res, file=resfn)
      } # END if ((run_or_collect == "run") || (run_or_collect == "both"))
      
      # Initialize the matrices if it's the first iteration
      if (i == start_treenum)
      {
        # If collecting results, create big empty tables to store the state probabilities at nodes and corners
        if ((run_or_collect == "collect") || (run_or_collect == "both"))
        {
          # For summarizing over the tree collection:
          # Set up an empty table to hold the state
          # probabilities of each run
          # numrows = numtrees * numnodes per tree
          # numcols = numstates
          example_trfn = newick_fns[[start_treenum]]
          example_tr = read.tree(example_trfn)
          numrows = length(example_tr$tip.label) + example_tr$Nnode
          numrows
          
          
          # Load results to res
          if (file.exists(resfns[[1]]) == FALSE)
          {
            stoptxt = paste0("STOP ERROR in run_bears_optim_on_multiple_trees(). Trying to load ", resfns[[1]], ", but the file is not found. Probably your 'run_or_collect' setting is 'run_or_collect='collect'', and you forgot to first run it on 'run_or_collect='run'' or 'run_or_collect='both''.")
            cat("\n\n")
            cat(stoptxt)
            cat("\n\n")
            stop(stoptxt)
          }
          load(resfns[[1]])
          numstates = ncol(res$ML_marginal_prob_each_state_at_branch_top_AT_node)
          
          # We will make these matrices double-size, just in case trees vary in size
          # (e.g., fossil inclusion or not), then cut the matrices down by excluding 
          # all NAs in the last column.
          
          # Create a big empty matrix (+1 to numcols for the OTUnames)
          # for the node states
          state_probs_at_nodes_across_all_trees = matrix(data=NA, nrow=length(resfns)*numrows*2, ncol=numstates+1)
          dim(state_probs_at_nodes_across_all_trees)
          
          # Create another matrix for the states at the corners (each corner is below
          # a node; thus probably nothing below the root)
          numrows = nrow(res$ML_marginal_prob_each_state_at_branch_bottom_below_node)
          numstates = ncol(res$ML_marginal_prob_each_state_at_branch_bottom_below_node)
          state_probs_at_corners_across_all_trees = matrix(data=NA, nrow=length(resfns)*numrows*2, ncol=numstates+1)
          dim(state_probs_at_corners_across_all_trees)
        } # end if ((run_or_collect == "collect") || (run_or_collect == "both"))
      } # end if (i == start_treenum)
      
      
      
      
      # If you want to run the summary after each tree, or just the summaries
      if ((run_or_collect == "collect") || (run_or_collect == "both"))
      {
        # Processing previously done inferences
        txt = paste("\nProcessing previous inferences under '", model_name, "' for tree #", i, " of ", start_treenum, "-", end_treenum, "...", sep="")
        cat(txt)
        
        
        # Do processing if desired
        # The goal is:
        # For each node and corner, get:
        # 1. The number of times that node appears in the sample
        #    (this should be approximately the PP of the bipartition)
        # 2. Whenever the node does appear, get the state probabilities
        #    at that node.
        # 3. Sum over all state probabilities and divide by #1
        # 4. Result is ancestral state probabilities averaged over the tree
        
        # For each tree, make an sorted text list of the OTUs descending from each node
        # Add to the state probabilities for each node
        trfn = newick_fns[i]
        tr = read.tree(trfn)
        trtable = prt(tr, printflag=FALSE, get_tipnames=TRUE)
        head(trtable)
        
        # Get the BioGeoBEARS results for this tree
        # Loads to "res"
        resfn = resfns[i]
        
        if (file.exists(resfns[i]) == FALSE)
        {
          stoptxt = paste0("STOP ERROR in run_bears_optim_on_multiple_trees(). Trying to load ", resfns[i], ", but the file is not found. Probably your 'run_or_collect' setting is 'run_or_collect='collect'', and you forgot to first run it on 'run_or_collect='run'' or 'run_or_collect='both''.")
          cat("\n\n")
          cat(stoptxt)
          cat("\n\n")
          stop(stoptxt)
        }
        load(resfn)
        
        names(res)
        state_probs_at_nodes = res$ML_marginal_prob_each_state_at_branch_top_AT_node
        state_probs_at_corners = res$ML_marginal_prob_each_state_at_branch_bottom_below_node
        dim(state_probs_at_corners)
        dim(trtable)
        
        # Store the results
        row_end = row_start-1 + nrow(state_probs_at_nodes)
        rownums = row_start:row_end
        # first columns are the state probabilities
        colnums = 1:numstates
        names_col = numstates + 1
        
        # Store the state probabilities at nodes
        state_probs_at_nodes_across_all_trees[rownums, colnums] = state_probs_at_nodes
        state_probs_at_corners_across_all_trees[rownums, colnums] = state_probs_at_corners
        
        # Store the OTUs descending from each node
        OTUnames = trtable$tipnames
        state_probs_at_nodes_across_all_trees[rownums, names_col] = OTUnames
        state_probs_at_corners_across_all_trees[rownums, names_col] = OTUnames
        
        # update row_start
        row_start = row_end + 1
      } # end if ((run_or_collect == FALSE) or (run_or_collect == "both"))
    } # end for (i in 1:length(newick_fns))
    
    # Delete rows that are all NA
    allNA <- function(tmprow)
    {
      row_is_allNA_TF = all(is.na(tmprow))
      return(row_is_allNA_TF)
    }
    
    if ((run_or_collect == "collect") || (run_or_collect == "both"))
    {
      rows_allNA_TF = apply(X=state_probs_at_nodes_across_all_trees, MARGIN=1, FUN=allNA)
      state_probs_at_nodes_across_all_trees = state_probs_at_nodes_across_all_trees[rows_allNA_TF==FALSE,]
      state_probs_at_corners_across_all_trees = state_probs_at_corners_across_all_trees[rows_allNA_TF==FALSE,]
      
      # Save the states		
      txt = paste("\n\nSaving state probabilities at nodes and corners across your tree sample to:\nworking directory: ", getwd(), "\n'", stateprobs_at_nodes_fn, "'\n'", stateprobs_at_corners_fn, "'\n(may be slow, ~1 minute)\n", sep="")
      cat(txt)
      
      # After the for-loop, save the ancestral states matrices if you like
      save(state_probs_at_nodes_across_all_trees, file=stateprobs_at_nodes_fn)
      save(state_probs_at_corners_across_all_trees, file=stateprobs_at_corners_fn)
    } # END if ((run_or_collect == "collect") || (run_or_collect == "both"))
    
  } else {
    
    # Load the states
    txt = paste("\n\nLoading state probabilities at nodes and corners across your tree sample from:\nworking directory: ", getwd(), "\n'", stateprobs_at_nodes_fn, "'\n'", stateprobs_at_corners_fn, "'\n(may be slow, ~1 minute)\n", sep="")
    cat(txt)
    # Load to: state_probs_at_nodes_across_all_trees
    load(file=stateprobs_at_nodes_fn)
    # Load to: state_probs_at_corners_across_all_trees
    load(file=stateprobs_at_corners_fn)
  } # end if runslow==TRUE
  
  
  # Also store parameter estimates
  optim_results_table = NULL
  get_optim_results = TRUE
  if (get_optim_results == TRUE)
  {
    cat("\n\nGetting ML parameter estimates for trees #", start_treenum, "-", end_treenum, ":\n", sep="")
    for (i in start_treenum:end_treenum)
      #for (i in 1:length(newick_fns))
      #for (i in 1:84)
    {
      cat(i, " ", sep="")
      
      # Get the BioGeoBEARS results for this tree
      # Loads to "res"
      if (file.exists(resfns[i]) == FALSE)
      {
        stoptxt = paste0("STOP ERROR in run_bears_optim_on_multiple_trees(). Trying to load ", resfns[i], ", but the file is not found. Probably your 'run_or_collect' setting is 'run_or_collect='collect'', and you forgot to first run it on 'run_or_collect='run'' or 'run_or_collect='both''.")
        cat("\n\n")
        cat(stoptxt)
        cat("\n\n")
        stop(stoptxt)
      }
      
      resfn = resfns[i]
      load(resfn)
      
      # Store the parameter estimates
      optim_results_table = rbind(optim_results_table, res$optim_result)
    } # end for (i in 1:length(newick_fns))
  } # end if (get_optim_results == TRUE)
  
  optim_results_mean = colMeans(optim_results_table)
  optim_results_sd = apply(X=optim_results_table, MARGIN=2, FUN=sd)
  
  # These results are atomic vectors, convert to data.frames as in optimx
  tmp_colnames = names(optim_results_mean)
  optim_results_mean = data.frame(matrix(data=optim_results_mean, nrow=1))
  names(optim_results_mean) = tmp_colnames
  optim_results_sd = data.frame(matrix(data=optim_results_sd, nrow=1))
  names(optim_results_sd) = tmp_colnames
  
  # Return results
  results_on_multiple_trees = list()
  
  
  # Add the state probabilities, if you get those...
  if ((run_or_collect == "collect") || (run_or_collect == "both"))
  {
    results_on_multiple_trees$state_probs_at_nodes_across_all_trees = state_probs_at_nodes_across_all_trees
    results_on_multiple_trees$state_probs_at_corners_across_all_trees = state_probs_at_corners_across_all_trees
  }
  
  # Filenames
  results_on_multiple_trees$newick_fns = newick_fns
  results_on_multiple_trees$geog_fns = geog_fns
  results_on_multiple_trees$resfns = resfns
  
  # Optim results (fast)
  results_on_multiple_trees$optim_results_table = optim_results_table
  results_on_multiple_trees$optim_results_mean = optim_results_mean
  results_on_multiple_trees$optim_results_sd = optim_results_sd
  
  
  if (plot_params == TRUE)
  {
    plot_params_from_multiple_trees(optim_results_table, BioGeoBEARS_run_object, optimx2012=FALSE)
  }
  
  return(results_on_multiple_trees)
}