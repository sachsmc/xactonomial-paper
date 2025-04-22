
# Improved small-sample inference for functions of parameters in the $k$-sample multinomial problem
## Reproducibility materials

### What is this?

This contains a series of R scripts to reproduce the simulation results in the above-titled paper. 

### How do I run this code? 

#### Prerequisites
You will need a recent version of R, along with the following packages `c("xactonomial", "data.table", "xtable")`. 

In the /example subdirectory, there are R scripts for several examples. No special hardward is needed for these, but they may take about an hour to run in total. 

The simulation study uses the `targets` package. To run the simulation study you will need the following packages in addition to the ones above `c("targets", "tarchetypes", "crew", "crew.cluster")`, and it is set up to run on a Slurm cluster. Read more about how to use it here: https://books.ropensci.org/targets/

#### Steps

1. Download the repository
2. Source the `_targets.R`
3. Run `tar_make()` to run the simulation. Warning, it takes a long time and it is set up to run on a SLURM cluster. If you have a similar cluster you can run `./run.sh &` at the console. 
4. Inspect the results using `tar_read` on one of the target names, with, e.g., `tar_read(xtable_list)`

```{r}
source("_targets.R")
tar_make()

tar_read(xtable_list)
```
