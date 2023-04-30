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
carnival <- CarnivalClass$new(scores_file_path = './data/lincs_troglitazone_PC3.csv',
                              organism = 'Human', 
                              cbc_solver_path = '~/.conda/envs/causal_reasoning/bin/cbc')
# run CARNIVAL 
carnival$runCarnival()

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





