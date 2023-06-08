CarnivalClass <- R6::R6Class(
  "CarnivalClass",
  public = list(
    scores_input = NULL, # a matrix or data frame (or a file path) consisting of genes on the row names and a single column
    organism = NULL,
    cbc_solver_path = NULL,
    datasets = NULL, # see OmnipathR::import_omnipath_interactions for options (ex. omnipath, dorothea)
    initiators = NULL, # initial nodes to use as the root of the network
                      # a vector of gene names 
    scores = NULL,
    progeny_scores = NULL,
    tf_activity_scores = NULL,
    network_sif = NULL,
    carnival_result = NULL,
    
    initialize = function(scores_input, organism, cbc_solver_path, datasets = c('omnipath'), initiators = NULL) {
      if (missing(scores_input) || missing(organism) || missing(cbc_solver_path)) {
        stop("All three arguments must be provided: scores_input, organism, cbc_solver_path")
      }
      
      self$scores_input <- scores_input
      self$organism <- organism
      self$cbc_solver_path <- cbc_solver_path
      self$datasets <- datasets 
      self$initiators <- initiators 
      if(!organism %in% c('Human', 'Mouse')){
        stop("Can't create a class for ",organism,". Allowed options are Human and Mouse")
      }
      if(is.data.frame(scores_input) || is.matrix(scores_input)) {
        self$scores <- as.matrix(scores_input)
      } else if(is.character(scores_input) && file.exists(scores_input)) {
        self$scores <- self$import_scores()
      } else {
          stop("Input must be a data frame, a matrix, or a valid file path.")
      }
      self$progeny_scores <- self$compute_progeny()
      self$tf_activity_scores <- self$compute_tf_activity()
      self$network_sif <- self$get_network_sif()
      # define initial nodes 
      if(is.null(self$initiators)) {
        iniMTX = base::setdiff(self$network_sif$source, self$network_sif$target)
        initiators = base::data.frame(base::matrix(data = NaN, nrow = 1, 
                                                   ncol = length(iniMTX)), 
                                      stringsAsFactors = F)
        colnames(initiators) = iniMTX
        self$initiators <- initiators
      } else {
        df <- base::data.frame(base::matrix(data = 1, nrow = 1, 
                                            ncol = length(self$initiators)), 
                               stringsAsFactors = F) 
        colnames(df) <- self$initiators
        self$initiators <- df
      }
      message("CarnivalClass object initialized.\n")
    },
    
    # import input file, perturbation readout  
    import_scores = function() {
      message(date(), " => importing perturbation data from ",self$scores_input)
      df <- read.csv(self$scores_input, row.names = 1)
      colnames(df)[1] <- 'score'
      return(as.matrix(df))
    },
    
    # compute pathway activity scores using progeny 
    compute_progeny = function() {
      message(date(), " => Calculating progeny pathway scores")
      progeny::progeny(self$scores, scale=TRUE, organism=self$organism, 
                       top = 100, perm = 10000, z_scores = FALSE)
    },
    # compute TF activities using dorothea/viper 
    compute_tf_activity = function() {
      message(date(), " => Computing TF activity scores")
      if(self$organism == 'Human') {
        data(dorothea_hs, package = "dorothea")
        dorothea_ann <- dorothea_hs
      } else if (self$organism == 'Mouse') {
        data(dorothea_mm, package = "dorothea")
        dorothea_ann <- dorothea_mm
      }
      regulons <- dorothea_ann[dorothea_ann$confidence %in% c('A', 'B', 'C'),]
      # compute TF activity scores using VIPER
      dorothea::run_viper(self$scores, regulons,
                          options =  list(minsize = 5, eset.filter = FALSE, 
                                          cores = 4, verbose = FALSE, nes = TRUE))
      
    },
    # import omnipath annotations and cleanup to get a SIF format
    # signed-directed graph 
    get_network_sif = function(datasets = c('omnipath')) {
      message(date(), " => importing OmniPath annotations in SIF format")
      
      # run carnival using progeny and TF activity scores
      organism_id <- NULL
      if(self$organism == 'Human') {
        organism_id <- 9606
      } else if(self$organism == 'Mouse') {
        organism_id <- 10090
      }
      # url
      omniR <- data.table::as.data.table(OmnipathR::import_omnipath_interactions(
        organism = organism_id, datasets = self$datasets))
      # convert to signed and directed 
      omnipath_sd <- omniR[consensus_direction == 1][consensus_stimulation == 1 | consensus_inhibition == 1]
      omnipath_sd[consensus_stimulation == 0]$consensus_stimulation <- -1
      omnipath_sd[consensus_inhibition == 1]$consensus_inhibition <- -1
      omnipath_sd[consensus_inhibition == 0]$consensus_inhibition <- 1
      
      # check consistency on consensus sign and select only those in a SIF format
      fields <- c('source_genesymbol', 'consensus_stimulation', 
                  'consensus_inhibition', 'target_genesymbol')
      sif <- subset(omnipath_sd, select = fields)[consensus_stimulation == consensus_inhibition]
      sif$consensus_stimulation <- NULL
      colnames(sif) <- c('source', 'interaction', 'target')
      # remove complexes
      sif <- unique(sif[!grepl(':|_', source)][!grepl(':|_', target)])
      message(date(), " => returning a network of ",nrow(sif), " interactions (complexes excluded)")
      return(sif)
    },
    runCarnival = function(threads = 20, absTFactivityThreshold = 2) {
      # pathway scores in list
      progenylist <- self$assignPROGENyScores(progeny = self$progeny_scores, 
                                         id = "gene", 
                                         access_idx = 1)
      # pick top TFs by inferred activity scores
      tfList <- self$tf_activity_scores[abs(self$tf_activity_scores) > absTFactivityThreshold,]
      message(date(), " => Using following TFs for causal signaling analysis: ",
              paste(names(tfList), collapse = " "))

      
      carnival_result = CARNIVAL::runCARNIVAL(inputObj = self$initiators,
                                     measObj = tfList, 
                                     netObj = self$network_sif, 
                                     weightObj = progenylist$score, 
                                     solverPath = self$cbc_solver_path, 
                                     solver =  "cbc",
                                     timelimit= 7200,
                                     threads = threads,
                                     mipGAP=0,
                                     poolrelGAP=0, 
                                     dir_name = file.path(getwd(), 'tmp_carnival'))
      self$carnival_result <- carnival_result
    }, 
    # detect communities in a given network (in sif format) 
    detect_communities = function(network) {
      df <- data.frame(network)
      g <- igraph::graph_from_data_frame(df[,c(1,3)])
      g <- igraph::set.edge.attribute(g, "interaction", value = df$Sign)
      # detect communities in the network
      comms <- igraph::communities(igraph::cluster_walktrap(g))
      #comms <- igraph::communities(igraph::cluster_label_prop(g))
      # sort by community size
      comms <- comms[names(sort(lengths(comms), decreasing = T))] 
      return(comms)
    },
    # assignPROGENyScores: Function taken from https://github.com/saezlab/transcriptutorial
    assignPROGENyScores = function (progeny = progeny, 
                                    id = "gene", access_idx = 1) {
      load(file = system.file("progenyMembers.RData",package="CARNIVAL"))
      if (id == "uniprot") {
        idx <- which(names(progenyMembers) == "uniprot")
        progenyMembers <- progenyMembers[[idx]]
      }
      else {
        idx <- which(names(progenyMembers) == "gene")
        progenyMembers <- progenyMembers[[idx]]
      }
      members <- matrix(data = , nrow = 1, ncol = 2)
      pathways <- colnames(progeny)
      ctrl <- intersect(x = access_idx, y = 1:nrow(progeny))
      if (length(ctrl) == 0) {
        stop("The indices you inserted do not correspond to \n              
             the number of rows/samples")
      }
      for (ii in 1:length(pathways)) {
        mm <- progenyMembers[[which(names(progenyMembers) ==
                                      pathways[ii])]]
        for (jj in 1:length(mm)) {
          members <- rbind(members, c(pathways[ii], mm[jj]))
        }
      }
      members <- members[-1, ]
      scores <- matrix(data = , nrow = nrow(progeny), ncol = nrow(members))
      colnames(scores) <- members[, 2]
      rownames(scores) <- rownames(progeny)
      members <- unique(members)
      for (i in 1:ncol(scores)) {
        for (j in 1:nrow(scores)) {
          scores[j, i] <- as.numeric(progeny[j, members[which(members[,
                                                                      2] == colnames(scores)[i]), 1]])
        }
      }
      pxList <- list()
      for (ii in 1:length(access_idx)) {
        pxList[[length(pxList) + 1]] <- as.data.frame(t(as.matrix(scores[access_idx[ii],
        ])))
      }
      names(pxList) <- rownames(progeny)[ctrl]
      return(pxList)
    }))



