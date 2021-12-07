# Activate project environment
import Pkg; Pkg.activate(".")

Pkg.instantiate()

# Load required packages
import CSV
using DataFrames
using DelimitedFiles
using JLD

import Ipopt # solving a system of nonlinear equations numerically
using JuMP # solving a system of nonlinear equations numerically

using LinearAlgebra
using Plots
using ProgressMeter
using SparseArrays
using Statistics
using StatsBase
using StatsPlots
using Turing

using EcologicalNetworks
using Mangal
using RandomBooleanMatrices

# Load custom functions
include(joinpath("code", "functions", "predict_links.jl")) # predict the number of links from the number of species
include(joinpath("code", "functions", "joint_degree_dist_maxent.jl")) # gives the joint degree distribution of maximum entropy given numbers of species and links
include(joinpath("code", "functions", "degree_dist_maxent.jl")) # gives the degree distribution of maximum entropy given numbers of species and links

# Load scripts
include(joinpath("code", "01_import_mangal_metadata.jl"))
include(joinpath("code", "02_predict_networks.jl"))
include(joinpath("code", "03_make_figures.jl"))

