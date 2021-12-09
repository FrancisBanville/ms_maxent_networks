## Activate project environment
import Pkg; Pkg.activate(".")

Pkg.instantiate()

## Load required packages

# Manipulating variables, data frames and files
import CSV 
using DataFrames 
using DelimitedFiles
using JLD 

# Solving systems of equations numerically 
import Ipopt
using JuMP 

# Doing statistics and models
using Distributions 
using LinearAlgebra 
using ProgressMeter
using SparseArrays
using Statistics
using StatsBase
using Turing

# Making plots
using Plots 
using StatsPlots 

# Analyzing ecological networks
using EcologicalNetworks 
using Mangal 
using RandomBooleanMatrices 

## Load custom functions
include(joinpath("code", "functions", "predict_links.jl")) # predict the number of links from the number of species
include(joinpath("code", "functions", "joint_degree_dist_maxent.jl")) # gives the joint degree distribution of maximum entropy given numbers of species and links
include(joinpath("code", "functions", "degree_dist_maxent.jl")) # gives the degree distribution of maximum entropy given numbers of species and links
include(joinpath("code", "functions", "joint_degree_seq.jl")) # gives the joint degree sequence of a unipartite network
include(joinpath("code", "functions", "joint_degree_dist.jl")) # gives the joint degree distribution of a unipartite network
include(joinpath("code", "functions", "network_maxent.jl")) # gives the joint degree distribution of a unipartite network

## Load scripts
include(joinpath("code", "01_import_mangal_metadata.jl"))
include(joinpath("code", "02_predict_networks.jl"))
include(joinpath("code", "03_make_figures.jl"))

