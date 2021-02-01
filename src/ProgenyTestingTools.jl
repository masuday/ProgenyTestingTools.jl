module ProgenyTestingTools

using Random
using Statistics
using SparseArrays
using LinearAlgebra

using CSV
using DataFrames
using OffsetArrays
using SparseMatrixDicts
using IterativeSolvers

# genetic parameters
"""
    par = PTParameters(var_p, h2_poly, h2_qtl, rep)

It defines the genetic parameters used in data generation supporting a single trait.
The parameters include `var_p` for the phenotypic variance, `h2_poly` for
the polygenic heritability, `h2_qtl` for the qtl heritability, and `rep` for 
the repeatability (for repeated records).
"""
struct PTParameters
   mu::Float64
   var_p::Float64
   h2_poly::Float64
   h2_qtl::Float64     # >0 if using genotypes
   rep::Float64
   var_g_sim::Float64  # true genetic variance in QMSim
end

# animal information
# multiple records but a single BV
mutable struct PTAnimal
   cg::Vector{Int}
   y::Vector{Float64}
   pe::Vector{Float64}
   e::Vector{Float64}
end

# central database for all animals
# will be updated when a new animal is born.
#
# maxCG: group x generation: increments when phenotype is assigned
# sg: group of sire
# dg: group of dam
mutable struct PTPopulation
   par::PTParameters
   hp::Bool
   maxAnimal::Int
   maxGroup::Int
   maxCG::Int
   df::DataFrame
#   male::Vector{Bool}
#   year::Vector{Int}
#   alive::Vector{Bool}
#   sampled::Vector{Bool}
#   pregnant::Vector{Bool}
#   genotyped::Vector{Bool}
#   pedigree::Vector{Tuple{Int,Int}}
#   inb::Vector{Float64}
#   pbv::Vector{Float64}
#   qbv::Vector{Float64}
#   tbv::Vector{Float64}
#   ebv::Vector{Float64}
#   gebv::Vector{Union{Missing,Float64}}
   animal::Vector{PTAnimal}
#   map::Union{Nothing,QMSimMap}
#   gfile::String
end

# group: a subset of the population
# sires, dams, id: linked to IDs in a population
mutable struct PTGroup
   pop::PTPopulation
   groupid::Int
   n::Int
   maxSire::Int
   maxDam::Int
   sires::Vector{Int}
   dams::Vector{Int}
   id::Vector{Int}
   generation::Vector{Int}
#   parentgroup::Vector{Tuple{Int,Int}}
end

export PTParameters, PTAnimal, PTPopulation, PTGroup
export check_parameters, get_var_poly, get_var_qtl, get_var_pe, get_var_error
export generate_population, migrate_from_hp!, generate_group, copy_group, add_new_animal!, 
       selectid, add_sires!, add_dams!, vacancy_for_sires, vacancy_for_dams, cull!

include("parameter.jl")
include("population.jl")
include("pedigree.jl")

end
