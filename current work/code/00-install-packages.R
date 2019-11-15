# LOAD LIBRARIES FOR CASE-CROSSOVER DATA ANALYSIS ###
# ALSO INCLUDE HELPER FUNCTIONS

# Packages required for computation
suppressWarnings(suppressMessages(require(ipoptr)))
install.packages("parallel")
install.packages("coxed")
install.packages("Matrix")
install.packages("trustOptim")
install.packages("matrixStats") # for LogSumExp
install.packages("purrr")
install.packages("rstanarm")
install.packages("INLA")

# tidyverse packages- required for analysis
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(lubridate)

# other
suppressWarnings(suppressMessages(require(rstanarm)))

### HELPER FUNCTIONS ###
