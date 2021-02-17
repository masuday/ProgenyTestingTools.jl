module ProgenyTestingTools

using Random
using Printf
using Statistics
using SparseArrays
using LinearAlgebra

using CSV
using QMSimFiles
using DataFrames
using OffsetArrays
using SparseMatrixDicts
using IterativeSolvers
using StatsBase: sample, Weights

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
"""
It defines a set of animals in a population.
The individual information will be stored in a dataframe `df` with the following columns.

- `id::Int` = sequential code (from 1)
- `male::Bool` = `true` for male, `false` for female
- `year::Int` = birth year
- `alive::Bool` = `true` for living individual
- `pregnant::Bool` = `true` for pregnant individual
- `genotyped::Bool` = `true` for genotyped individual
- `candidate::Bool` = `true` for candidate for selection
- `status::Int` = user-defined status code (default = 0)
- `sire::Int` = sire code
- `dam::Int` = dam code
- `siregroup::Int` = group ID of sire
- `damgroup::Int` = group ID of dam
- `nprog::Int` = number of progeny
- `nrecprog::Int` = number of progeny with record(s)
- `inb::Float64` = inbreeding coefficient (from 0.0 to 1.0)
- `pbv::Float64` = polygenic breeding value (PBV)
- `qbv::Float64` = QTL breeding value (QBV)
- `tbv::Float64` = true (total) breeding value: TBV = PBV + QBV
- `ebv::Float64` = estimated breeding value by pedigree BLUP
- `gebv::Float64` = genomic EBV by "pseudo" genomic prediction
- `rel::Float64` = user-supplied reliability of GEBV
- `nrec::Int` = number of record(s)
- `firsty::Float64` = the first observation
- `lasty::Float64` = the last observation
- `avgy::Float64` = average of observations
"""
mutable struct PTPopulation
   par::PTParameters
   hp::Bool
   maxAnimal::Int
   maxGroup::Int
   maxCG::Int
   df::DataFrame
   animal::Vector{PTAnimal}
   map::Union{Nothing,QMSimMap}
   gfile::String
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
end

import Base: show, write
export PTParameters, PTAnimal, PTPopulation, PTGroup
export check_parameters, get_var_poly, get_var_qtl, get_var_pe, get_var_error
export assign_phenotype!
export generate_population, migrate_from_hp!, generate_group, copy_group, add_new_animal!, 
       selectid, random_sampling, add_sires!, add_dams!, vacancy_for_sires, vacancy_for_dams, 
       mark_candidate!, change_status!, cull!
export update_inbreeding!
export mating!, calving!
export genetic_evaluation!, update_approximated_reliability!, write_files_for_blup

export test_mating!, regular_mating!, cull_old_bulls!, cull_old_cows!,
       mean_inbreeding, number_of_active_females, phenotyped_cows, expected_frequency_of_mating

include("parameter.jl")
include("population.jl")
include("phenotype.jl")
include("pedigree.jl")
include("mating.jl")
include("io.jl")
include("evaluation.jl")
include("batch.jl")

end
