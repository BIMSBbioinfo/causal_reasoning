# causal_reasoning_omics
A collection of tools, which, when given a readout of a perturbation experiment, runs tools such as CARNIVAL/causalR to infer the causal factors that leads to the observed outcome

# Dependencies

R libraries
```
library(CARNIVAL)
library(progeny)
library(dorothea)
library(OmnipathR)
library(data.table)
```

## To Run Carnival

```
# import CarnivalClass
source('./src/carnival.R')

# instantiate the class; imports perturbation data, calculates pathway and TF activity scores; processes omnipath data
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


