# Activate project environment
import Pkg; Pkg.activate(".")

Pkg.instantiate()

# Load required packages
import CSV
using DataFrames
using DelimitedFiles
using JLD

import Ipopt
using JuMP
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

using SimpleSDMLayers
using ArchGDAL
using GDAL

# Load custom functions
include(joinpath("src", "functions.jl")) 

# Load scripts
include(joinpath("src", "01_import_mangal_metadata.jl"))
include(joinpath("src", "02_predict_networks.jl"))
include(joinpath("src", "03_make_figures.jl"))
include(joinpath("src", "04_make_maps.jl"))


