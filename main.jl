# Activate project environment
import Pkg; Pkg.activate(".")

Pkg.instantiate()

# Load required packages
import CSV
using DataFrames

using Plots
using Statistics
using StatsBase
using StatsPlots
using Turing

using EcologicalNetworks
using Mangal

# Load scripts
include("src/01_import_mangal_metadata.jl")
