CarnivalClass <- R6::R6Class(
  "CarnivalClass",
  public = list(
    scores_file_path = NULL,
    organism = NULL,
    cbc_solver_path = NULL,
    datasets = NULL, # see OmnipathR::import_omnipath_interactions for options (ex. omnipath, dorothea)
    scores = NULL,
    progeny_scores = NULL,
    tf_activity_scores = NULL,
    omnipath_sif = NULL,
    carnival_result = NULL,
    
    initialize = function(scores_file_path, organism, cbc_solver_path, datasets = c('omnipath')) {
      if (missing(scores_file_path) || missing(organism) || missing(cbc_solver_path)) {
        stop("All three arguments must be provided: scores_file_path, organism, cbc_solver_path")
      }
      
      self$scores_file_path <- scores_file_path
      self$organism <- organism
      self$cbc_solver_path <- cbc_solver_path
      self$datasets <- datasets 
      if(!organism %in% c('Human', 'Mouse')){
        stop("Can't create a class for ",organism,". Allowed options are Human and Mouse")
      }
      self$scores <- self$import_scores()
      self$progeny_scores <- self$compute_progeny()
      self$tf_activity_scores <- self$compute_tf_activity()
      self$omnipath_sif <- self$get_omnipath_sif()
      cat("CarnivalClass object initialized.\n")
    },
    
    # import input file, perturbation readout  
    import_scores = function() {
      message(date(), " => importing perturbation data from ",self$scores_file_path)
      df <- read.csv(self$scores_file_path, row.names = 1)
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
    get_omnipath_sif = function(datasets = c('omnipath')) {
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
      sif <- sif[!grepl(':', source) | grepl(':', target)]
      return(sif)
    },
    runCarnival = function() {
      # pathway scores in list
      progenylist <- self$assignPROGENyScores(progeny = self$progeny_scores, 
                                         id = "gene", 
                                         access_idx = 1)
      # tf scores in list
      tfList <- self$generateTFList(self$tf_activity_scores, top='all', access_idx = 1)
      
      iniMTX = base::setdiff(self$omnipath_sif$source, self$omnipath_sif$target)
      initiators = base::data.frame(base::matrix(data = NaN, nrow = 1, 
                                                 ncol = length(iniMTX)), 
                                    stringsAsFactors = F)
      colnames(initiators) = iniMTX
      print(head(initiators))
      
      carnival_result = CARNIVAL::runCARNIVAL( inputObj = initiators,
                                     measObj = tfList$score, 
                                     netObj = self$omnipath_sif, 
                                     weightObj = progenylist$score, 
                                     solverPath = self$cbc_solver_path, 
                                     solver =  "cbc",
                                     timelimit= 7200,
                                     threads = 20,
                                     mipGAP=0,
                                     poolrelGAP=0)
      self$carnival_result <- carnival_result
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
    },
    # generateTFList: Function taken from https://github.com/saezlab/transcriptutorial
    generateTFList = function (df = df, top = 50, access_idx = 1) {
      if (top == "all") {
        top <- nrow(df)
      }
      if (top > nrow(df)) {
        warning("Number of to TF's inserted exceeds the number of actual TF's in the\n            data frame. All the TF's will be considered.")
        top <- nrow(df)
      }
      ctrl <- intersect(x = access_idx, y = 1:ncol(df))
      if (length(ctrl) == 0) {
        stop("The indeces you inserted do not correspond to \n              
             the number of columns/samples")
      }
      returnList <- list()
      for (ii in 1:length(ctrl)) {
        tfThresh <- sort(x = abs(df[, ctrl[ii]]), decreasing = TRUE)[top]
        temp <- which(abs(df[, ctrl[ii]]) >= tfThresh)
        currDF <- matrix(data = , nrow = 1, ncol = top)
        colnames(currDF) <- rownames(df)[temp[1:top]]
        currDF[1, ] <- df[temp[1:top], ctrl[ii]]
        currDF <- as.data.frame(currDF)
        returnList[[length(returnList) + 1]] <- currDF
      }
      names(returnList) <- colnames(df)[ctrl]
      return(returnList)
    }))



