
# Exact confidence intervals for functions of parameters in the $k$-sample multinomial problem
## Reproducibility materials

### What is this?

This contains a series of R scripts to reproduce the simulation results in the above-titled paper. 

### How do I run this code? 

#### Prerequisites
You will need a recent version of R, along with the following packages `c("xactonomial", "data.table", "xtable")`. 

The `xactonomial` package is available here: https://github.com/sachsmc/xactonomial

This uses the `targets` package for reproducibility. Read more about how to use it here: https://books.ropensci.org/targets/

#### Steps

1. Download the repository
2. Source the `_targets.R`
3. Run `tar_make()` or `tar_make_clustermq()` to run the simulation. Warning, it takes a long time and it is set up to run on a SLURM cluster. If you have a similar cluster you can run `./run.sh &` at the console. 
4. Inspect the results using `tar_read` on one of the target names, with, e.g., `tar_read(xtable_list)`

```{r}
source("_targets.R")
tar_make()

tar_read(xtable_list)
```
