#!/usr/bin/Rscript

args <- commandArgs(TRUE)
condition1 <- args[1]
condition2 <- args[2]

source("diff_exp_sleuth.R")

diff_exp_sleuth(condition1, condition2)
