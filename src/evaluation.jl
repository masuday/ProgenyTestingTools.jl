# functions for evaluation

"""
    genetic_evaluation!(pop::PTPopulation; method="blup", repeated=true, rel=0.5, verbose=false, updategebv=false, cg=false)

Perform genetic/genomic evaluation using `model`.
Put `repeated=true` if PE is considered in the model
This function supports the following models.

- `model="blup"`: traditional pedigree-based BLUP; stored as `ebv`
- `model="tbv"`: true breeding value (TBV) with reliability `rel` (defined as the squared correlation between TBV and predictions); stored as `gebv`

For "blup", the model does not have contemporary groups (CG); instead uses the general mean.
With `cg=true`, the group effect will be considered.

For "tbv", the prediction is `TBV + randn()*sqrt([(1-rel)/rel]*var_g)`, where `TBV` is the true breeding value, `randn()` is a standard-normal deviation, `rel` is the reliability, and `var_g` is the genetic variance defined as the initial parameter.
This formula does not account for a change of genetic variance due to selection or drift.
If you provide an array for `rel`, this function will take it as individual reliability to generate GEBV.
Note that, once GEBV is calculated, it will not be updated unless you have `updategebv=true`.
"""
function genetic_evaluation!(pop::PTPopulation; method::String="blup", repeated::Bool=true, rel::Union{Vector{Float64},Float64}=0.5, verbose::Bool=false, 
   keep_file::Bool=false, updategebv::Bool=false, cg::Bool=false)
   if method=="blup"
      write_files_for_blup(pop, "phenotype.txt", "pedigree.txt", repeated=repeated)
      genetic_evaluation_blup!(pop, "phenotype.txt", "pedigree.txt", repeated, verbose=verbose, cg=cg)
      if !keep_file
         rm("phenotype.txt")
         rm("pedigree.txt")
      end
   elseif method=="tbv"
      genetic_evaluation_tbv!(pop,rel,updategebv)
   else
      error("unknown method: $(method)")
   end
end

@inbounds function genetic_evaluation_tbv!(pop::PTPopulation, rel::Float64, updategebv::Bool)
   par = pop.par
   var_g = get_var_poly(par) + get_var_qtl(par)
   if rel>0.0
      var_error = ((1-max(rel,1.0))/max(rel,1.0)) * var_g
   else
      var_error = 100*var_g
   end
   sd_error = sqrt(var_error)
   if updategebv
      if sd_error>0
         pop.df[:,:gebv] .= pop.df[:,:tbv] .+ (randn(pop.maxAnimal).*sd_error)
      else
         pop.df[:,:gebv] .= pop.df[:,:tbv]
      end
   else
      @inbounds for i=1:pop.maxAnimal
         if ismissing(pop.gebv[i])
            pop.df[i,:gebv] = pop.df[i,:tbv] + (randn() * sd_error)
         end
      end
   end
end

function genetic_evaluation_tbv!(pop::PTPopulation, rel::Vector{Float64}, updategebv::Bool)
   par = pop.par
   var_g = get_var_poly(par) + get_var_qtl(par)
   @inbounds for i=1:pop.maxAnimal
      if !ismissing(rel[i]) && rel>0.0
         var_error = ((1-max(rel[i],1.0))/max(rel[i],1.0)) * var_g
      else
         var_error = 100*var_g
      end
      sd_error = sqrt(var_error)
      if updategebv || ismissing(pop.df[i,:gebv])
         pop.df[i,:gebv] = pop.df[i,:tbv] + (randn() * sd_error)
      end
   end
end

function genetic_evaluation_tbv!(pop::PTPopulation, updategebv::Bool)
   rel = pop.df.rel
   genetic_evaluation_tbv!(pop, rel, updategebv)
end

"""
    rel = approximated_reliability(pop::PTPopulation; de_extra::Float64=0.0)

Calculate (roughly) approximated reliability based on the number of recorded progeny and own records.
You can supply an extra daugher equivalent (DE) `de_extra` to the reliability if needed.
"""
function approximated_reliability(pop::PTPopulation; de_extra::Float64=0.0)
   par = pop.par
   h2 = par.h2_poly + par.h2_qtl
   rep = par.rep
   k = (4-h2)/h2
   ndau = pop.df[:,nrecprogeny]
   nrec = pop.df[:,nrec]
   rel = zeros(pop.maxAnimal)
   @inbounds for i=1:pop.maxAnimal
      # very rough approximation; see VanRaden and Wiggans (1991) for strict methods
      # PA
      s = pop.df[i,:sire]
      d = pop.df[i,:dam]
      if s>0; rel_s=rel[s]; else; rel_s=0.0; end
      if d>0; rel_d=rel[d]; else; rel_d=0.0; end
      rel_pa = (rel_s + rel_d)/4
      de_pa = k*rel_pa/(1-rel_pa)
      # own records
      if nrec[i]>0
         rel_y = nrec[i]*h2/(1 + (nrec[i]-1)*rep)
      else
         rel_y = 0.0
      end
      de_y = k*rel_y/(1-rel_y)
      # progeny contribution
      de_dau = ndau[i]
      # total DE
      de_animal = de_pa + de_y + de_dau + de_extra
      rel_animal = de_animal/(de_animal + k)
      rel[i] = rel_animal
   end
   return rel
end

function update_approximated_reliability!(pop::PTPopulation; de_extra::Float64=0.0)
   rel = approximated_reliability(pop,de_extra=de_extra)
   pop.df[:,:rel] .= rel
   return nothing
end

"""
    write_files_for_blup(pop::PTPopulation, datafile, pedfile; repeated=true)

Write the data and pedigree to files.
With `repeated=false`, this function writes only the first record per animal.
"""
function write_files_for_blup(pop::PTPopulation, datafile::String, pedfile::String; repeated::Bool=true)
   io = open(datafile, "w");
   write(io, pop, header=false, missing="0", repeated=repeated, skipmissing=true)
   close(io)

   io = open(pedfile, "w");
   write_pedigree(io, pop, header=false)
   close(io)
end

function genetic_evaluation_blup!(pop::PTPopulation, datafile::String, pedfile::String, repeated::Bool=true; algorithm::String="pcg", verbose::Bool=false, cg::Bool=false)
   (lhs,rhs) = build_mme(pop, datafile, pedfile, repeated, cg)
   if algorithm=="pcg"
      sol = similar(rhs)
      sol .= 0.0
      dlhs = diag(lhs)
      @inbounds for i=1:length(dlhs)
         if dlhs[i] ≈ 0.0
            dlhs[i] = 1.0
         end
      end
      dlhs = Diagonal(dlhs)
      #cg!(sol, lhs, rhs, Pl=dlhs, verbose=true, log=true)
      cg!(sol, lhs, rhs, Pl=dlhs, verbose=verbose)
   else
      sol = lhs \ rhs
   end
   pop.ebv .= sol[1:pop.maxAnimal]
end

function build_mme(pop::PTPopulation, datafile::String, pedfile::String, repeated::Bool, cg::Bool)
   par = pop.par
   var_a = get_var_poly(par)
   var_g = get_var_qtl(par)
   var_u = var_a + var_g
   var_pe = get_var_pe(par)
   var_e = get_var_error(par)
   var_p = par.var_p

   lambda_u = var_e/var_u
   lambda_pe = var_e/var_pe
   if repeated
      if var_pe<=0.0 && repeated
         @warn "no use of repeatability model because of repeatability <= 0"
         repeatability_model = false
      else
         repeatability_model = true
      end
   else
      repeatability_model = false
   end

   # details
   if repeatability_model
      # y = mu + u + pe + e
      pos_y = 13
      neff = 3
      random_type = ["u","pe","fixed"]
      offset = [0, pop.maxAnimal, pop.maxAnimal+pop.maxAnimal]
      if cg
         pos_eff = [1, 1, 12]
         nlev = [pop.maxAnimal, pop.maxAnimal, pop.maxCG]
      else
         pos_eff = [1, 1, 11]
         nlev = [pop.maxAnimal, pop.maxAnimal, 1]
      end
   else
      # y = mu + u + e
      pos_y = 13
      neff = 2
      random_type = ["u","fixed"]
      offset = [0, pop.maxAnimal]
      if cg
         pos_eff = [1, 12]
         nlev = [pop.maxAnimal, pop.maxCG]
      else
         pos_eff = [1, 11]
         nlev = [pop.maxAnimal, 1]
      end
   end
   
   return build_mme_lhs_rhs(datafile,pedfile,neff,pos_eff,offset,pos_y,nlev,random_type,lambda_u,lambda_pe)
end

function build_mme_lhs_rhs(datafile::String,pedfile::String,neff::Int,pos_eff::Vector{Int},offset::Vector{Int},pos_y::Int,nlev::Vector{Int},random_type::Vector{String},lambda_u::Float64,lambda_pe::Float64)
   neq = sum(nlev)
   lhs = SparseMatrixDict(neq,neq)
   rhs = zeros(neq)
   build_mme_lhs_rhs_fixed!(lhs,rhs,datafile,neff,pos_eff,offset,pos_y)
   for i in 1:neff
      if random_type[i]=="u"
         build_mme_lhs_rhs_ainverse!(lhs,pedfile,offset[i],lambda_u)
      elseif random_type[i]=="pe"
         build_mme_lhs_rhs_diagonal!(lhs,offset[i],nlev[i],lambda_pe)
      end
   end
   return sparse(lhs),rhs
end

function build_mme_lhs_rhs_fixed!(lhs::SparseMatrixDict,rhs::Vector{Float64},datafile::String,neff::Int,pos_eff::Vector{Int},offset::Vector{Int},pos_y::Int)
   # design matrix
   open(datafile) do io
      while !eof(io)
         line = readline(io)
         items = split(line)
         y = tryparse(Float64,items[pos_y])
         if isnothing(y) || y==0.0
            continue
         end
         level = tryparse.(Int,items[pos_eff])
         @inbounds for i in 1:neff
            addri = offset[i] + level[i]
            for j in 1:neff
               addrj = offset[j] + level[j]
               lhs[addri,addrj] = lhs[addri,addrj] + 1.0
            end
            rhs[addri] = rhs[addri] + y
         end
      end
   end
end

function build_mme_lhs_rhs_ainverse!(lhs::SparseMatrixDict,pedfile::String,offset::Int,lambda::Float64)  
   # A-inverse
   w = (1.0, -0.5, -0.5)
   open(pedfile) do io
      while !eof(io)
         line = readline(io)
         items = split(line)
         ped = tryparse.(Int,items[1:3])
         code = tryparse(Float64,items[4])/1000.0
         @inbounds for i in 1:3
            addri = offset + ped[i]
            for j in 1:3
               addrj = offset + ped[j]
               if ped[i]>0 && ped[j]>0
                  val = w[i]*w[j]*code*lambda
                  lhs[addri,addrj] = lhs[addri,addrj] + val
               end
            end
         end
      end
   end
end

function build_mme_lhs_rhs_diagonal!(lhs::SparseMatrixDict,offset::Int,nlev::Int,lambda::Float64)
   # identity matrix
   @inbounds for i in 1:nlev
      addri = offset + i
      lhs[addri,addri] = lhs[addri,addri] + lambda
   end
end
