# functions for phenotype

"""
    assign_phenotype!(group::PTGroup, gen=nothing; male=false, female=true, year=nothing, repeated=true, limit=nothing, idlist=1:group.pop.maxAnimal)

Put a phenotype to animals in generation `gen` in `group`, born in `year` (`nothing` to ignore the birth year).
If `gen` is `nothing` or a negative integer, this function gives all animals in this group.
By default, only females have phenotypes (`male=false` and `female=true`).
Also, an animal can have multiple (repeated records) due to `repeated=true`.
To limit the number of repeated reacord, use `limit` (default = `nothing`, no limitation).
Given `idlist`, only the animals in this list will have the phenotype; by default, all animals will have the phenotype.
The list has to have been sorted.
"""
function assign_phenotype!(group::PTGroup; gen::Union{Nothing,Int}=nothing, male::Bool=false, female::Bool=true, year::Union{Nothing,Int}=nothing,
   repeated::Bool=true, verbose::Bool=false, limit::Union{Nothing,Int}=nothing, idlist::Union{UnitRange{Int},Vector{Int}}=1:group.pop.maxAnimal)
   pop = group.pop
   par = group.pop.par
   if isnothing(gen)
      gencode = -1
   else
      gencode = gen
   end
   if gencode<0
      seqlist = collect(1:group.n)
      mingen = minimum(group.generation)
      maxgen = maximum(group.generation)
      cglist = Dict{Int,Int}()
      for i in mingen:maxgen
         pop.maxCG = pop.maxCG + 1
         cglist[i] = pop.maxCG
      end
   else
      seqlist = (1:group.n)[group.generation .== gencode]
      if length(seqlist)<1
         error("no such generation in this group")
      end
      pop.maxCG = pop.maxCG + 1
      cglist = Dict{Int,Int}(gen => pop.maxCG)
   end

   if !isnothing(year)
      yearlist = pop.df[group.id[seqlist],:year]
      seqlist = seqlist[yearlist .== year]
   end
   if verbose
      if length(seqlist)<1
         @warn "no animals will be phenotyped this time."
      end
   end

   _assign_phenotype!(group, seqlist, cglist, male, female, repeated, limit, idlist=idlist)

   return nothing
end

function assign_phenotype!(groups::Vector{PTGroup}; gen::Union{Nothing,Int}=nothing, male::Bool=false, female::Bool=true, year::Union{Nothing,Int}=nothing, 
   repeated::Bool=true, verbose::Bool=false, limit::Union{Nothing,Int}=nothing, idlist::Union{UnitRange{Int},Vector{Int}}=1:group.pop.maxAnimal)
   for group in groups
      assign_phenotype!(group, gen=gen, male=male, female=female, year=year, repeated=repeated, verbose=verbose, limit=limit, idlist=idlist)
   end
   return nothing
end

function _assign_phenotype!(group::PTGroup, seqlist::Vector{Int}, cglist::Dict{Int,Int}, male::Bool, female::Bool, 
   repeated::Bool, limit::Union{Nothing,Int}; idlist::Union{UnitRange{Int},Vector{Int}}=1:group.pop.maxAnimal)
   pop = group.pop
   par = group.pop.par

   # some CG would not be assigned if repeated=false
   var_a = get_var_poly(par)
   var_g = get_var_qtl(par)
   var_pe = get_var_pe(par)
   var_e = get_var_error(par)
   sd_pe = sqrt(var_pe)
   sd_e = sqrt(var_e)
   @inbounds for seq in seqlist
      id = group.id[seq]
      if !(male && pop.df[id,:male]) && !(female && !pop.df[id,:male]) && group.pop.df[id,:alive]
         # out of condition
         continue
      end
      if length(searchsorted(idlist,id))<1
         # not on the list
         continue
      end
      gen = group.generation[seq]
      cg = cglist[gen]
      animal = pop.animal[id]
      nrec = length(animal.y)
      if isnothing(limit)
         nlim = nrec+1
      else
         nlim = limit
      end
      if (nrec==0 || (nrec>=1 && repeated)) && nrec<=nlim
         tbv = pop.df[id,:tbv]
         if nrec==0
            this_pe = sd_pe*randn()
         else
            this_pe = animal.pe[1]
         end
         this_e = sd_e*randn()
         this_y = par.mu + tbv + this_pe + this_e
         push!(animal.cg, cg)
         push!(animal.y, this_y)
         push!(animal.pe, this_pe)
         push!(animal.e, this_e)
         s = pop.df[id,:sire]
         d = pop.df[id,:dam]
         if s>0
            pop.df[s,:nrecprog] = pop.df[s,:nrecprog] + 1
         end
         if d>0
            pop.df[d,:nrecprog] = pop.df[d,:nrecprog] + 1
         end
         pop.df[id,:nrec] = length(animal.y)
         pop.df[id,:firsty] = animal.y[1]
         pop.df[id,:lasty] = animal.y[end]
         pop.df[id,:avgy] = mean(animal.y)
      end
   end

   return nothing
end
