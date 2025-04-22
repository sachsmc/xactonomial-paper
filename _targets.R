# Created by use_targets().
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline

# Load packages required to define the pipeline:
library(targets)
library(tarchetypes) # Load other packages as needed.
library(crew)
library(crew.cluster)
# Set target options:
bigmem_controller <- crew_controller_slurm(
  workers = 50, 
  seconds_wall = 8000, 
  name = "bigmem_controller",
  options_cluster = crew_options_slurm(
    script_lines =c("module load gcc/13.2.0",
                    "module load openjdk/20.0.0",
                    "module load R/4.4.2"),
    memory_gigabytes_per_cpu = 10)
)

smallmem_controller <- crew_controller_slurm(
  workers = 250, 
  seconds_wall = 8000,
  name = "smallmem_controller",
  options_cluster = crew_options_slurm(
    script_lines =c("module load gcc/13.2.0",
                    "module load openjdk/20.0.0",
                    "module load R/4.4.2"),
    memory_gigabytes_per_cpu = 2)
)

tar_option_set(
  packages = c("xactonomial", "data.table", "xtable"), # packages that your targets need to run
  # format = "qs", # Optionally set the default storage format. qs is fast.
  #
    controller = crew_controller_group(smallmem_controller, bigmem_controller),
    resources = tar_resources(
    crew = tar_resources_crew(controller = "smallmem_controller")
  )
)

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
  list(n = 10, theta = c(2/10, 5/10, 3/10), maxit = 100),
  list(n = 10, theta = c(1/3, 1/3, 1/3), maxit = 100),
  list(n = 20, theta = c(1/3, 1/3, 1/3), maxit = 100),
  list(n = 30, theta = c(1/3, 1/3, 1/3), maxit = 100),
  list(n = 40, theta = c(1/3, 1/3, 1/3), maxit = 100),
  list(n = 20, theta = c(.4, .4, .2), maxit = 100),
  list(n = 30, theta = c(.3, .3, .25, .15), maxit = 100),
  list(n = 40, theta = c(.3, .3, .25, .15), maxit = 100)

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
             memory = "transient",
             resources = tar_resources(
               crew = tar_resources_crew(controller = "bigmem_controller")
             )
             ),
  tar_target(bc_summary,
             data.table(simulation_bc)[, .(
               cover_boot = mean(cover_boot),
               cover_xact = mean(cover_xact),
               width_xact = mean(width_xact),
               type1_err = mean(p_xact < 0.05),
               true_psi = true_psi[1]), by = .(n, d, theta)][,
                    k := sapply(strsplit(theta, ", "), length) / d
                    ]
             ),
  # tar_target(simulation_lb,
  #            run_one(i, setting_lb[[1]]$n, d = setting_lb[[1]]$d, setting_lb[[1]]$theta, psi = psi_lb, limits = c(-1, 1)),
  #            pattern = cross(i, setting_lb),
  #            memory = "transient"
  #            ),
  # tar_target(lb_summary,
  #            data.table(simulation_lb)[, .(cover_boot = mean(cover_boot),
  #                                          cover_xact = mean(cover_xact),
  #                                          width_xact = mean(width_xact),
  #                                          type1_err = mean(p_xact < 0.05),
  #                                          true_psi = true_psi[1]), by = .(n, d, theta)]
  #           ),
  tar_target(simulation_max,
            run_one(i, setting_max[[1]]$n, d = 1, setting_max[[1]]$theta, psi = psi_max,
                    limits = c(1/length(setting_max[[1]]$theta), 1),
                    maxit = setting_max[[1]]$maxit, csize = 200, psi_max = TRUE),
            pattern = cross(i, setting_max),
            memory = "transient"
            ),
  tar_target(max_summary,
             data.table(simulation_max)[, .(cover_boot = mean(cover_boot),
                                            cover_xact = mean(cover_xact),
                                            width_xact = mean(width_xact),
                                            type1_err = mean(p_xact < 0.05),
                                            true_psi = true_psi[1]), by = .(n, d, theta)]
  ),
  tar_target(finish,
             send_mail(bc_summary, max_summary)
             ),
  tar_target(xtable_list,
    list(
      bhatty = print(xtable(bc_summary[order(n, d, true_psi), .(n, d, true_value = true_psi, xactonomial_coverage = cover_xact* 100,
      bootstrap_coverage = cover_boot * 100)], digits = c(-1, 0, 0, 2, 1, 1)), include.rownames = FALSE),
      max = print(xtable(max_summary[order(n, d, true_psi), .(n, d, true_value = true_psi, xactonomial_coverage = cover_xact* 100,
      bootstrap_coverage = cover_boot * 100)], digits = c(-1, 0, 0, 2, 1, 1)), include.rownames = FALSE)
    )
  )
)
