# functions for mating

# mating between two groups
"""
    pregnant_dams = mating!(mgroup::PTGroup, fgroup::PTGroup, progeny="fgroup"; ...)

Perform mating between sires from `mgroup` and dams from `fgroup`.
The progeny will be in a group defined by `progeny`; `"mgroup"` for the male (sire) group or `"fgroup"` for the female (dam) group.
The pregnant dams (stored in the population database) will be excluded from the mating plan.
The generation of progeny is defined as the largest generation number of parents (in `progeny`) + 1.
It returns a list of dams that have been pregnant in this mating opportinity.
The following options are available.

- `n`: `"all"` for mating as many females as possible in particular `method` and `plan`. Or, you can limit the number of females for mating.
- `calving`: `true` for mating and calving at the same time in this function (but still returning `pregnant_dams`). By default, it is `false` to keep pregnancy for dams (use `calving!` to clear the pregnancy status).
- `method`: mating strategy. The default is `"dairy_standard_ai"`, in which a cow has only one chance to mate one of the sires.
- `plan`: `"once_per_female"` for using each dam once (default), or `"once_per_male"` for using each sire once.
- `pmp`: proportion of progeny male (default=0.5).
- `year`: birth year of progeny.
"""
function mating!(mgroup::PTGroup, fgroup::PTGroup, progeny="fgroup"; 
   n::Union{Int,String}="all", year::Int=0, method::String="dairy_standard_ai", plan::String="once_per_female", 
   pmp::Float64=0.5, calving::Bool=false, updateinb::Bool=true, verbose::Bool=false)

   if length(mgroup.sires)<1 || length(fgroup.dams)<1
      error("No sires or no dams")
   end

   if progeny=="mgroup"
      group = mgroup
      maxGen = maximum( group.generation[findall(x->in(x,group.sires),group.id)] )
   elseif progeny=="fgroup"
      group = fgroup
      maxGen = maximum( group.generation[findall(x->in(x,group.dams),group.id)] )
   else
      error("invalid argument progeny=$(progeny)")
   end

   sgroup = mgroup.groupid
   dgroup = fgroup.groupid
   non_pregnant_dams = selectid([:pregnant] => x->x==false, fgroup, idlist=fgroup.dams)
   if typeof(n) == String
      if n=="all"
         ndams = length(non_pregnant_dams)
      else
         error("n: to be `all`, or integer value")
      end
   else
      ndams = min(n,length(non_pregnant_dams))
   end
   if verbose
      if length(non_pregnant_dams)<1
         @warn "There is no non-pregnant dams in this group (gr=$(group.gid))."
      else
         println("Mating gr$(mgroup.gid) x gr$(fgroup.gid); non-pregnant dams:$(length(non_pregnant_dams)) out of $(length(fgroup.dams))")
      end
   end
   if method=="dairy_standard_ai"
      pregnant_dams = _mating_dairy_standard_ai!(group, mgroup.sires, non_pregnant_dams, sgroup, dgroup, maxGen+1, pmp, method="random", plan=plan, year=year, n=ndams)
   else
      error("unknown method $(method) for mating")
   end

   # keep the pregnancy status
   if !calving
      fgroup.pop.df[pregnant_dams,:pregnant] .= true
   end

   # update inbreeding
   if updateinb
      update_inbreeding!(group.pop.df.sire, group.pop.df.dam, group.pop.df.inb)
   end

   return pregnant_dams
end

# standard dairy AI mating
#  - litter size = 1
#  - one cow receives one service; producing n progeny from n dams.
#  - no action if the cow has already been serviced in the same year.
function _mating_dairy_standard_ai!(group::PTGroup, sires::Vector{Int}, dams::Vector{Int}, sgroup::Int, dgroup::Int, gen::Int, pmp::Float64; 
   method::String="random", plan::String="once_per_female", year::Int=0, n::Int=length(dams))
   var_a = get_var_poly(group.pop.par)
   var_g = get_var_qtl(group.pop.par)
   var_g_sim = group.pop.par.var_g_sim
   nai = 0
   ndams = length(dams)
   # non-pregnant dam list
   damlist = copy(dams)
   # pregnant dams
   pregnant_dams_id = Dict{Int,Bool}()
   if method=="random"
      shuffle!(damlist)
      if plan=="once_per_male"
         @inbounds for s in sires
            if s>0
               if isempty(damlist)
                  break
               end
               d = pop!(damlist)
               pregnant_dams_id[d] = true
               male = rand_male_calf(pmp)
               generate_and_add_animal!(group, s, d, sgroup, dgroup, male, gen, year, var_a, var_g, var_g_sim)
               nai = nai + 1
            end
            if nai>=n
               break
            end
         end
      elseif plan=="once_per_female"
         @inbounds for d in damlist
            if d>0
               pregnant_dams_id[d] = true
               s = rand(sires)
               male = rand_male_calf(pmp)
               generate_and_add_animal!(group, s, d, sgroup, dgroup, male, gen, year, var_a, var_g, var_g_sim)
               nai = nai + 1
            end
            if nai>=n
               break
            end
         end
      else
         error("unknown mating plan: $(plan)")
      end
   else
      error("unknown mating method: $(method)")
   end
   return Int.(keys(pregnant_dams_id))
end

function rand_male_calf(pmp::Float64)
   if rand()<=pmp
      return true
   else
      return false
   end
end

function generate_and_add_animal!(group::PTGroup, s::Int, d::Int, sgroup::Int, dgroup::Int, male::Bool, gen::Int, year::Int, var_a::Float64, var_g::Float64, var_g_sim::Float64)
   pop = group.pop
 
   # pedigree-based (polygenic) breeding valur
   pa = (pop.df[s,:pbv] + pop.df[d,:pbv])/2
   ms = get_ms_deviation(pop.df.inb[s], pop.df.inb[d], var_a)
   pbv = pa + ms

   # QTL-based breeding value
   if pop.df[s,:genotyped] && pop.df[d,:genotyped]
      gs = read_qmsim_individual_hdf5(pop.map,pop.gfile,s)
      gd = read_qmsim_individual_hdf5(pop.map,pop.gfile,d)
      ga = mating(pop.map,gs,gd)
      qbv = sqrt(var_g/var_g_sim)*ga.tbv
      hassnp = true
   else
      qbv = sqrt(var_g)*randn()  # just a random number
      hassnp = false
   end

   # total breeding values
   tbv = pbv + qbv

   # for future use
   cg = zeros(Int,0)
   y = zeros(Float64,0)
   pe = zeros(Float64,0)
   e = zeros(Float64,0)
   animal = PTAnimal(cg,y,pe,e)

   # add to population
   add_new_animal!(pop, male, animal, year=year, sire=s, dam=d, siregroup=sgroup, damgroup=dgroup, inb=0.0, pbv=pbv, qbv=qbv, tbv=tbv, ebv=0.0, gebv=missing, genotyped=hassnp)
   if hassnp
      add_qmsim_individual_hdf5(pop.map,pop.gfile,ga)
   end

   # add to group
   group.n = group.n + 1
   push!(group.id, pop.maxAnimal)
   push!(group.generation, gen)

   return nothing
end

function get_ms_deviation(inbs::Float64, inbd::Float64, var_a::Float64, k::Float64=0.0)
   sd_ms = sqrt( (1-k) * ( 0.5*(1 - 0.5*(inbs+inbd)) ) * var_a )
   return sd_ms*randn()
end

function update_polygenic_ms!(pop::PTPopulation, id::Int, k::Float64=0.0)
   s = pop.df[i,:sire]
   d = pop.df[i,:dam]
   if s==0
      us = 0.0
      inbs = 0.0
   else
      us = pop.df[s,:pbv]
      inbs = pop.df[s,:inb]
   end
   if d==0
      ud = 0.0
      inbd = 0.0
   else
      ud = pop.df[d,:pbv]
      inbd = pop.df[d,:inb]
   end
   pa = (us + ud)/2
   ms = get_ms_deviation(inbs, inbd, var_a, k)
   pop.df[id,:pbv] = pa + ms
end

function update_polygenic_ms!(pop::PTPopulation, idlist::Vector{Int}, k::Float64=0.0)
   for id in idlist
      update_polygenic_ms!(pop, id, k)
   end
end

# calving
"""
    calving(pop::PTPopulation)
    calving(group::PTGroup)

Change the pregnancy status for females from "pregnant" to "non-pregnant".
"""
function calving!(pop::PTPopulation)
   pop.df[:,:pregnant] .= false
   return nothing
end
function calving!(group::PTGroup)
   group.pop.df[group.dams,:pregnant] .= false
   return nothing
end
function calving!(groups::Vector{PTGroup})
   for group in groups
      calving!(group)
   end
   return nothing
end
