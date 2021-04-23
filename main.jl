# Activate project environment
import Pkg; Pkg.activate(".")

Pkg.instantiate()

# Load required packages
import CSV
using DataFrames

using Plots
using ProgressMeter
using Statistics
using StatsBase
using StatsPlots
using Turing

using EcologicalNetworks
using Mangal
using RandomBooleanMatrices

# Load scripts
include("src/01_import_mangal_metadata.jl")
include("src/02_predict_networks.jl")
