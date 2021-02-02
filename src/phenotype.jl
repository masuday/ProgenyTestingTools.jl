# functions for phenotype

"""
    assign_phenotype!(group::PTGroup, gen="all"; male=false, female=true, year=-1, repeated=true)

Put a phenotype to animals in generation `gen` in `group`, born in `year` (-1 to ignore the birth year).
If `gen` is `"all"` or a negative integer, this function gives all animals in this group.
By default, only females have phenotypes (`male=false` and `female=true`).
Also, an animal can have multiple (repeated records) due to `repeated=true`.
"""
function assign_phenotype!(group::PTGroup; gen::Union{Int,String}="all", male::Bool=false, female::Bool=true, year::Int=-1, repeated::Bool=true, verbose::Bool=false, limit::Int=100)
   pop = group.pop
   par = group.pop.par
   if typeof(gen) == String
      if gen=="all"
         gencode = -1
      else
         error("gen = `all` or integer number")
      end
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

   if year>=0
      yearlist = pop.df[group.id[seqlist],:year]
      seqlist = seqlist[yearlist .== year]
   end
   if verbose
      if length(seqlist)<1
         @warn "no animals will be phenotyped this time."
      end
   end

   _assign_phenotype!(group, seqlist, cglist, male, female, repeated)

   return nothing
end

function assign_phenotype!(groups::Vector{PTGroup}; gen::Union{Int,String}="all", male::Bool=false, female::Bool=true, year::Int=-1, repeated::Bool=true, verbose::Bool=false, limit::Int=100)
   for group in groups
      assign_phenotype!(group, gen=gen, male=false, female=true, year=-1, repeated=repeated, verbose=false)
   end
   return nothing
end

function _assign_phenotype!(group::PTGroup, seqlist::Vector{Int}, cglist::Dict{Int,Int}, male::Bool, female::Bool, repeated::Bool, limit::Int=100)
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
      gen = group.generation[seq]
      cg = cglist[gen]
      animal = pop.animal[id]
      nrec = length(animal.y)
      if (nrec==0 || (nrec>=1 && repeated)) && nrec<=limit
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
