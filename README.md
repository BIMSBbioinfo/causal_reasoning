<p align="left">
  <img alt="chocolate" src="img/chocolate.jpeg" width="25%" height="25%">
  <figcaption> Fig.1 - The amount of chocolate consumption is correlated with the prospect of winning a Nobel Prize </figcaption>
</p>

# Causal Reasoning Analysis 

A collection of tools, which, when given a readout of a perturbation experiment, runs tools such as CARNIVAL/causalR to infer the causal factors that leads to the observed outcome

# Dependencies


## CARNIVAL
R libraries
```
library(CARNIVAL)
library(progeny)
library(dorothea)
library(OmnipathR)
library(data.table)
library(igraph) 
```


CBC solver (for integer linear programming) optimisation
```
conda install -c conda-forge coincbc
```

## CausalR 

```
library(OmnipathR)
library(CausalR)
library(data.table)
```

## To Run Carnival

For more info on CARNIVAL, see https://www.bioconductor.org/packages/release/bioc/html/CARNIVAL.html

```
# import CarnivalClass
source('./src/carnival.R')

# instantiate the class; imports perturbation data
# calculates pathway and TF activity scores; processes omnipath data
carnival <- CarnivalClass$new(scores_file_path = 'data/lincs_troglitazone_PC3.csv',
                              organism = 'Human',
                              datasets = 'dorothea',
                              cbc_solver_path = '~/.conda/envs/causal_reasoning/bin/cbc')
# run carnival
carnival$runCarnival(threads = 32, absTFactivityThreshold = 2)

# detect sub modules in the resulting subnetworks
comms <- carnival$detect_communities(carnival$carnival_result$sifAll[[1]])

# plot subgraphs (plot the subgraph corresponding to the largest community)
carnival$plot_subgraph(carnival$get_subgraph(carnival$network_sif, #carnival$carnival_result$sifAll[[1]],
                                             comms[[1]]),
                       max_size = 50, margin = -0.5)


# access the results
carnival$carnival_result

# access other scores (pathway scores and TF activity scores)
carnival$progeny_scores
carnival$tf_activity_scores
```

## To Run CausalR 

For more info, see: https://www.bioconductor.org/packages/release/bioc/html/CausalR.html

```
# import CausalrClass
source('./src/causalR.R')

# instantiate the class
# it is important to define what is a good threshold for up/down regulation 
# based on the scores found in scores_file_path 
# datasets: which network to import (use dorothea network for testing, 
# causalR can take long to run on the whole omnipath network)
caus <- CausalrClass$new(scores_file_path = 'data/lincs_troglitazone_PC3.csv', 
                         organism = 'Human', 
                         upregulationThreshold = 1,
                         downregulationThreshold = -1,
                         datasets = 'dorothea') 

# Run "Sequential Causal Analysis of Networks" 
caus$run_scanr(numberOfDeltaToScan = 5, topNumGenes = 100)

# Access the results
caus$result

# Top genes who rank among the top genes most consistently:
caus$result$commonNodes

```





