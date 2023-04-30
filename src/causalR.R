# causalR class definition

CausalrClass <- R6::R6Class(
  "CausalrClass",
  public = list(
    scores_file_path = NULL,
    upregulationThreshold = NULL, # score threshold to consider which genes are upregulated
    downregulationThreshold = NULL, # score threshold to consider which genes are downregulated
    organism = NULL,
    scores = NULL,
    omnipath_sif = NULL,
    datasets = NULL, # see OmnipathR::import_omnipath_interactions for options (ex. omnipath, dorothea)
    sif_path = NULL,
    ccg = NULL, 
    regulatedGenes = NULL, # computed from scores based on upregulationThreshold/downregulationThreshold
    result = NULL,
    
    initialize = function(scores_file_path, 
                          organism, 
                          upregulationThreshold, 
                          downregulationThreshold, 
                          datasets = 'omnipath') {
      if (missing(scores_file_path) || missing(organism) || 
          missing(upregulationThreshold) || missing(downregulationThreshold)) {
        stop("All arguments must be provided: scores_file_path, organism, 
             upregulationThreshold, downregulationThreshold")
      }
      
      self$scores_file_path <- scores_file_path
      self$organism <- organism
      self$upregulationThreshold <- upregulationThreshold
      self$downregulationThreshold <- downregulationThreshold
      self$datasets <- datasets 
      self$scores <- self$import_scores()
      self$omnipath_sif <- self$get_omnipath_sif()
      # harmonize data
      self$sif_path <- self$harmonize()
      # create a Computational Causal Graph (CCG)
      self$ccg <- self$compute_CCG()
      # filter up/down regulated genes 
      self$regulatedGenes <- self$get_regulated_genes()
      cat("CausalrClass object initialized.\n")
    },
    
    # import input file, perturbation readout  
    import_scores = function() {
      message(date(), " => importing perturbation data from ",self$scores_file_path)
      df <- read.csv(self$scores_file_path, row.names = 1)
      colnames(df)[1] <- 'score'
      return(as.matrix(df))
    },
    # import omnipath annotations and cleanup to get a SIF format
    # signed-directed graph 
    get_omnipath_sif = function(datasets = c('dorothea')) {
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
    
    # harmonize scores and interaction data 
    # remove genes that don't exist in both
    harmonize = function() {
      message(date(), " => harmonizing scores and network data")
      sif <- self$omnipath_sif
      genes <- intersect(rownames(self$scores), unique(c(sif$source, sif$target)))
      scores_subset <- self$scores[genes,,drop=F]
      # also update the network
      sif_subset <- data.frame(unique(sif[source %in% genes | target %in% genes]))
      sif_subset[sif_subset$interaction == 1,]$interaction <- 'Activates'
      sif_subset[sif_subset$interaction == -1,]$interaction <- 'Inhibits'
      # update 
      self$scores <- scores_subset
      self$omnipath_sif <- sif_subset
      sif_path <- file.path(getwd(), 'omnipath.sif')
      write.table(sif_subset, file = sif_path, 
                  sep = '\t', quote = F, 
                  row.names = F, col.names = F)
      return(sif_path)
    }, 
    
   
    compute_CCG = function() {
      # create a Computational Causal Graph
      message(date(), " => Creating the Causal Graph")
      ccg <- CausalR::CreateCCG(filename = self$sif_path)
      return(ccg)
    },
    
    get_regulated_genes = function() {
      message(date(), " => Filtering up/down regulated genes")
      # select genes that are up/down regulated (based on the experiment)
      df <- data.frame('nodes' = rownames(self$scores))
      df$status <- 0
      df[self$scores[,1] > self$upregulationThreshold,]$status <- 1
      df[self$scores[,1] < self$downregulationThreshold,]$status <- -1
      return(df)
    }, 

    # scan the network for each regulated gene to test the hypotheses that the genes 
    # are causal 
    # see Sequential Causal Analysis of Networks
    # 
    run_scanr = function(numCores = 32, numberOfDeltaToScan = 5, 
                         topNumGenes = 100) {
      r <- CausalR::runSCANR(network = self$ccg, experimentalData = self$regulatedGenes,
                             numberOfDeltaToScan = numberOfDeltaToScan, 
                             topNumGenes = topNumGenes, doParallel = TRUE, 
                             numCores = numCores, writeResultFiles = F, writeNetworkFiles = 'none')
      self$result <- r
      return(r)
    }
  ))
      






