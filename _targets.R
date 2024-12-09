# Created by use_targets().
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline

# Load packages required to define the pipeline:
library(targets)
library(tarchetypes) # Load other packages as needed.

# Set target options:
tar_option_set(
  packages = c("xactonomial", "data.table", "xtable") # packages that your targets need to run
  # format = "qs", # Optionally set the default storage format. qs is fast.
  #
  # For distributed computing in tar_make(), supply a {crew} controller
  # as discussed at https://books.ropensci.org/targets/crew.html.
  # Choose a controller that suits your needs. For example, the following
  # sets a controller with 2 workers which will run as local R processes:
  #
  #   controller = crew::crew_controller_local(workers = 2)
  #
  # Alternatively, if you want workers to run on a high-performance computing
  # cluster, select a controller from the {crew.cluster} package. The following
  # example is a controller for Sun Grid Engine (SGE).
  #
  #   controller = crew.cluster::crew_controller_sge(
  #     workers = 50,
  #     # Many clusters install R as an environment module, and you can load it
  #     # with the script_lines argument. To select a specific verison of R,
  #     # you may need to include a version string, e.g. "module load R/4.3.0".
  #     # Check with your system administrator if you are unsure.
  #     script_lines = "module load R"
  #   )
  #
  # Set other options as needed.
)

# tar_make_clustermq() is an older (pre-{crew}) way to do distributed computing
# in {targets}, and its configuration for your machine is below.
options(clustermq.scheduler = "slurm")
options(clustermq.template = "clustermq.tmpl")
#options(clustermq.scheduler = "multicore")
# tar_make_future() is an older (pre-{crew}) way to do distributed computing
# in {targets}, and its configuration for your machine is below.
# Install packages {{future}}, {{future.callr}}, and {{future.batchtools}} to allow use_targets() to configure tar_make_future() options.

# Run the R scripts in the R/ folder with your custom functions:
tar_source()
# source("other_functions.R") # Source other scripts as needed.

simulation_settings <- list(
  list(n = 10,
  d = 2,
  theta = c(.4, .3, .3, .3, .3, .4)),
  list(n = 10,
  d = 2,
  theta = c(.7, .2, .1, .1, .2, .7)),
  list(n = 20,
  d = 2,
  theta = c(.4, .3, .3, .3, .3, .4)),
  list(n = 20,
  d = 2,
  theta = c(.7, .2, .1, .1, .2, .7)),
  list(n = 10,
  d = 2,
  theta = c(.3, .3, .2, .2, .2, .2, .3, .3)),
  list(n = 10,
  d = 2,
  theta = c(.6, .2, .1, .1, .1, .1, .2, .6)),
  list(n = 20,
  d = 2,
  theta = c(.3, .3, .2, .2, .2, .2, .3, .3)),
  list(n = 20,
  d = 2,
  theta = c(.6, .2, .1, .1, .1, .1, .2, .6)),
  list(n = 10,
  d = 3,
  theta = c(.7, .2, .1, .3, .3, .4, .1, .2, .7)),
  list(n = 5,
  d = 4,
  theta = c(.7, .2, .1, .3, .3, .4, .1, .2, .7, .4, .3, .3)),
  list(n = 5,
  d = 5,
  theta = c(.5, .2, .3, .3, .3, .4, .3, .2, .5, .4, .3, .3, .2, .3, .5))
)


bounds_settings <- list(
  list(n = 5,
  d = 2,
  theta = c(.45, .15, .3, .1, .05, .15, .4, .4)),
  list(n = 5,
  d = 2,
  theta = c(.75, .1, .1, .05, .05, .15, .15, .65)),
  list(n = 5,
  d = 2,
  theta = c(.85, .05, .05, .05, .1, .15, .05, .7)),
  list(n = 5,
  d = 2,
  theta = c(.85, .05, .05, .05, .9, .05, .04, .01)),
  list(n = 10,
  d = 2,
  theta = c(.45, .15, .3, .1, .05, .15, .4, .4)),
  list(n = 10,
  d = 2,
  theta = c(.75, .1, .1, .05, .05, .15, .15, .65)),
  list(n = 10,
  d = 2,
  theta = c(.85, .05, .05, .05, .1, .15, .05, .7)),
  list(n = 10,
  d = 2,
  theta = c(.85, .05, .05, .05, .9, .05, .04, .01)),
  list(n = 15,
  d = 2,
  theta = c(.45, .15, .3, .1, .05, .15, .4, .4)),
  list(n = 15,
  d = 2,
  theta = c(.75, .1, .1, .05, .05, .15, .15, .65)),
  list(n = 15,
  d = 2,
  theta = c(.85, .05, .05, .05, .1, .15, .05, .7)),
  list(n = 15,
  d = 2,
  theta = c(.85, .05, .05, .05, .9, .05, .04, .01)),
  list(n = 20,
  d = 2,
  theta = c(.45, .15, .3, .1, .05, .15, .4, .4)),
  list(n = 20,
  d = 2,
  theta = c(.75, .1, .1, .05, .05, .15, .15, .65)),
  list(n = 20,
  d = 2,
  theta = c(.85, .05, .05, .05, .1, .15, .05, .7)),
  list(n = 20,
  d = 2,
  theta = c(.85, .05, .05, .05, .9, .05, .04, .01))
)

max_settings <- list(
  list(n = 10, theta = c(2/10, 5/10, 3/10), maxit = 300),
  list(n = 10, theta = c(1/3, 1/3, 1/3), maxit = 1000),
  list(n = 20, theta = c(1/3, 1/3, 1/3), maxit = 1000),
  list(n = 30, theta = c(1/3, 1/3, 1/3), maxit = 1000),
  list(n = 40, theta = c(1/3, 1/3, 1/3), maxit = 1000),
  list(n = 20, theta = c(.4, .4, .2), maxit = 1000),
  list(n = 30, theta = c(.3, .3, .25, .15), maxit = 5000),
  list(n = 40, theta = c(.3, .3, .25, .15), maxit = 5000)

)

# Replace the target list below with your own:
list(
  tar_target(i, seq(1:2000)),
  tar_target(setting, simulation_settings),
  tar_target(setting_lb, bounds_settings),
  tar_target(setting_max, max_settings),
  tar_target(simulation_bc,
             run_one(i, setting[[1]]$n, d = setting[[1]]$d, setting[[1]]$theta, psi = psi_bc, limits = c(0, 1)),
             pattern = cross(i, setting),
             memory = "transient"
             ),
  tar_target(bc_summary,
             data.table(simulation_bc)[, .(
               cover_boot = mean(cover_boot),
               cover_xact = mean(cover_xact),
               true_psi = true_psi[1]), by = .(n, d, theta)][,
                    k := sapply(strsplit(theta, ", "), length) / d
                    ]
             ),
  tar_target(simulation_lb,
             run_one(i, setting_lb[[1]]$n, d = setting_lb[[1]]$d, setting_lb[[1]]$theta, psi = psi_lb, limits = c(-1, 1)),
             pattern = cross(i, setting_lb),
             memory = "transient"
             ),
  tar_target(lb_summary,
             data.table(simulation_lb)[, .(cover_boot = mean(cover_boot),
                                           cover_xact = mean(cover_xact),
                                           true_psi = true_psi[1]), by = .(n, d, theta)]
             ),
  tar_target(simulation_max,
            run_one(i, setting_max[[1]]$n, d = 1, setting_max[[1]]$theta, psi = psi_max, limits = c(0, 1),
                    maxit = setting_max[[1]]$maxit, csize = 200),
            pattern = cross(i, setting_max),
            memory = "transient"
            ),
  tar_target(max_summary,
             data.table(simulation_max)[, .(cover_boot = mean(cover_boot),
                                            cover_xact = mean(cover_xact),
                                            true_psi = true_psi[1], by = .(n, d, theta))]
  ),
  tar_target(finish,
             send_mail(bc_summary, lb_summary, max_summary)
             ),
  tar_target(xtable_list,
    list(
      lower_bound = print(xtable(lb_summary[order(n, true_psi), .(n, true_value = true_psi, xactonomial_coverage = cover_xact* 100,
      bootstrap_coverage = cover_boot * 100)], digits = c(-1, 0, 2, 1, 1)), include.rownames = FALSE),
      bhatty = print(xtable(bc_summary[order(n, d, true_psi), .(n, d, true_value = true_psi, xactonomial_coverage = cover_xact* 100,
      bootstrap_coverage = cover_boot * 100)], digits = c(-1, 0, 0, 2, 1, 1)), include.rownames = FALSE),
      max = print(xtable(max_summary[order(n, d, true_psi), .(n, d, true_value = true_psi, xactonomial_coverage = cover_xact* 100,
      bootstrap_coverage = cover_boot * 100)], digits = c(-1, 0, 0, 2, 1, 1)), include.rownames = FALSE),
    )
  )
)
