# functions for populations
"""
   pop = generate_population(par::PTParameters, nm::Int, nf::Int)
   pop = generate_population(par::PTParameters, qmsimfile, nm::Int, nf::Int; gfile)

Generates a new population with `nm` males and `nf` females.
For a historical population, you should provide the number of males and females.
For a new population, you do not need the numbers because the animals will be
generated in a course of simulation.

When using genotypes, you have to provide external genotypes, generated
by QMSim. The data should be converted to an HDF5 file by the QMSimData
package. The QMSim file must be needed for any populations because the file
defines the genetic architecture. Also, you should provide `gfile` to store
gennotypes that going to be generated in simulation.
```
# hostorical population with 5 males and 10 females
hp = generate_population(par, nm=5, nf=10)

# hostorical population with genotypes
hp = generate_population(par, "qmsimdata.h5", nm=5, nf=10)

# new population with no animals
pop = generate_population(par)

# new population with genotypes
pop = generate_population(par, "qmsimdata.h5", nm=0, nf=0, gfile="pop.h5")
```
"""
function generate_population(par::PTParameters; nm::Int=0, nf::Int=0, map::Union{Nothing,QMSimMap}=nothing, gfile="")
   # check parameters
   check_parameters(par)

   # initialization
   maxAnimal = nm+nf
   hp = ifelse(maxAnimal>0, true, false)
   maxGroup = 0
   maxCG = 0

   # generate polygenic breeding values
   pbv = zeros(Float64,maxAnimal)
   if par.h2_poly > 0.0
      var_poly = par.var_p * par.h2_poly
      pbv .= randn(maxAnimal)*sqrt(var_poly)
   end

   # breeding value due to QTL
   qbv = zeros(Float64,maxAnimal)

   # total breeding value
   tbv = zeros(Float64,maxAnimal)
   tbv .= pbv .+ qbv

   # genotyping status
   if isnothing(map)
      genotyped = fill(false, maxAnimal)
   else
      genotyped = fill(true, maxAnimal)
   end

   df = DataFrame(
      id = collect(1:maxAnimal),
      male = vcat([true for i=1:nm],[false for i=1:nf]),
      year = zeros(Int, maxAnimal),
      alive = fill(true, maxAnimal),
      pregnant = fill(false, maxAnimal),
      genotyped = genotyped,
      candidate = fill(false, maxAnimal),
      status = zeros(Int, maxAnimal),
      sire = zeros(Int, maxAnimal),
      dam = zeros(Int, maxAnimal),
      siregroup = zeros(Int, maxAnimal),
      damgroup = zeros(Int, maxAnimal),
      nprog = zeros(Int, maxAnimal),
      nrecprog = zeros(Int, maxAnimal),
      inb = zeros(Float64,maxAnimal),
      pbv = pbv,
      qbv = qbv,
      tbv = tbv,
      ebv = zeros(Float64,maxAnimal),
      gebv = Vector{Union{Missing,Float64}}(missing,maxAnimal),
      rel = Vector{Union{Missing,Float64}}(missing,maxAnimal),
      nrec = zeros(Int, maxAnimal),
      firsty = Vector{Union{Missing,Float64}}(missing,maxAnimal),
      lasty = Vector{Union{Missing,Float64}}(missing,maxAnimal),
      avgy = Vector{Union{Missing,Float64}}(missing,maxAnimal))

   # generate empty entry
   animal = Vector{PTAnimal}()
   @inbounds for id in 1:maxAnimal
      cg = zeros(Int,0)
      y = zeros(Float64,0)
      pe = zeros(Float64,0)
      e = zeros(Float64,0)
      push!(animal, PTAnimal(cg,y,pe,e))
   end

   return PTPopulation(par,hp,maxAnimal,maxGroup,maxCG,df,animal,map,gfile)
end

function generate_population(par::PTParameters, qmsimfile::String; nm::Int=0, nf::Int=0, map=nothing, gfile="")
   map = read_qmsim_map_hdf5(qmsimfile)
   if nm>0 || nf>0
      if gfile==qmsimfile || gfile==""
         @info "hostorical population: qmsimfile is the same as gfile, which will be read-only."
         pop = generate_population(par, nm=nm, nf=nf, map=map, gfile=qmsimfile)
      else
         throw(ArgumentError("For a historical population, qmsimfile must be the same as gfile."))
      end
   else
      if gfile==qmsimfile
         throw(ArgumentError("gfile must be provided as a different file as qmsimfile."))
      end
      # genomic data file created by QMSim/QMSimFiles
      if gfile==""
         h5file = "snp_"*randstring(4)*".h5"
      else
         h5file = gfile
      end
      pop = generate_population(par, nm=nm, nf=nf, map=map, gfile=h5file)
      create_qmsim_map_hdf5(map, pop.gfile)
   end
   return pop
end

"""
    add_new_animal!(pop::PTPopulation)

Add an animal to population `pop`.
The added animals lose the pedigree information.
"""
function add_new_animal!(pop::PTPopulation, male::Bool, animal::PTAnimal; 
                     genotyped=false, candidate=false, year=0, alive=true, pregnant=false, status=0, sire=0, dam=0, siregroup=0, damgroup=0, inb=0.0,
                     pbv=0.0, qbv=0.0, tbv=0.0, ebv=0.0, rel=missing, gebv=missing)
   pop.maxAnimal = pop.maxAnimal + 1
   (nrec,firsty,lasty,avgy) = get_phenotypic_information(animal.y)
   newdata = (pop.maxAnimal, male, year, alive, pregnant, genotyped, candidate, status, sire, dam, siregroup, damgroup, 0, 0, inb, pbv, qbv, tbv, ebv, gebv, rel, nrec, firsty, lasty, avgy)
   push!(pop.df, newdata)
   push!(pop.animal, animal)
   return nothing
end

function add_new_animal!(pop::PTPopulation, row::DataFrameRow{DataFrames.DataFrame,DataFrames.Index}, animal::PTAnimal)
   pop.maxAnimal = pop.maxAnimal + 1
   (nrec,firsty,lasty,avgy) = get_phenotypic_information(animal.y)
   newdata = (id=pop.maxAnimal, male=row.male, year=row.year, alive=row.alive, pregnant=row.pregnant, genotyped=row.genotyped, 
              candidate=row.candidate, status=row.status, sire=row.sire, dam=row.dam, siregroup=row.siregroup, damgroup=row.damgroup, nprog=0, nrecprog=0, 
              inb=row.inb, pbv=row.pbv, qbv=row.qbv, tbv=row.tbv, ebv=row.ebv, gebv=row.gebv, rel=row.rel,
              nrec=nrec, firsty=firsty, lasty=lasty, avgy=avgy)
   push!(pop.df, newdata)
   push!(pop.animal, animal)
   return nothing
end

function get_phenotypic_information(y::Vector{Float64})
   nrec = length(y)
   if nrec>0
      firsty = y[1]
      lasty = y[end]
      avgy = mean(y)
   else
      firsty = missing
      lasty = missing
      avgy = missing
   end
   return (nrec,firsty,lasty,avgy)
end

"""
    newidlist = migrate_from_hp!(hp::PTPopulation, pop::PTPopulation, idlist::Vector{Int}; year=0, halt=true)

Migrate some animals from a historical population `hp` to another population `pop`.
The migrated animals have new IDs in the new population.
The animals will have "dead" status in the historical population.
It returns a new ID list `newidlist` in the new population.
With `halt=true` (default), this function throws when there are no enough individuals in the historical population.
"""
function migrate_from_hp!(hp::PTPopulation, pop::PTPopulation, idlist::Vector{Int}; year::Union{Int,Vector{Int}}=0, halt::Bool=true)
   unique_idlist = unique(idlist)
   n = length(unique_idlist)
   if n==0 && halt
      error("empty ID list; maybe there is no more individual in the historical population.")
   end
   newidlist = Vector{Int}()
   @inbounds for id in unique_idlist
      # check id in the historiacl population
      if id > hp.maxAnimal
         error("No such ID $(id) in the first population")
      end
      if hp.df[id,:alive]
         # new animal for the new population
         add_new_animal!(pop, hp.df[id,:], hp.animal[id])
         pop.df[pop.maxAnimal,:id] = pop.maxAnimal
         # moved
         hp.df[id,:alive] = false
         # assigned new id
         push!(newidlist, pop.maxAnimal)
      else
         if halt
            error("id $(id) in hp is not alive; maybe there is no more individual in the historical population.")
         else
            @warn "id $(id) in hp is not alive, and has not been added to the population."
         end
      end

      # (hp) id to (pop) pop.maxAnimal
      # assigned new genotypes
      if hp.df[id,:genotyped]
         # copy genotypes from the database
         g = read_qmsim_individual_hdf5(hp.map,hp.gfile,id)
         add_qmsim_individual_hdf5(pop.map,pop.gfile,g)
         pop.df[pop.maxAnimal,:genotyped] = true
         pop.df[pop.maxAnimal,:qbv] = g.tbv
         pop.df[pop.maxAnimal,:tbv] = pop.df[pop.maxAnimal,:pbv] + pop.df[pop.maxAnimal,:qbv]
      end
   end
   # assign year
   assign_year!(pop, newidlist, year)
   return newidlist
end

"""
    assign_year!(pop::PTPopulation, idlist::Vector{Int}, year::Vector{Int})
    assign_year!(pop::PTPopulation, idlist::Vector{Int}, year::Int)

Assign birth year to animals defied by `idlist`.
By default, the "years" will be evenly assigned to the animals.
"""
function assign_year!(pop::PTPopulation, idlist::Vector{Int}, year::Vector{Int})
   m = length(year)
   j = 1
   @inbounds for i in idlist
      if 0<i && i<=pop.maxAnimal
         pop.df[i,:year] = year[j]
         j = j + 1
         if j>m
            j = 1
         end
      end
   end
   return nothing
end
function assign_year!(pop::PTPopulation, idlist::Vector{Int}, year::Int)
   assign_year!(pop, idlist, [year])
   return nothing
end

# new group with founder males and females
function generate_group(pop::PTPopulation; sires::Vector{Int}=Int[], dams::Vector{Int}=Int[], aliveonly::Bool=true, empty::Bool=false)
   if empty
      pop.maxGroup = pop.maxGroup + 1
      return PTGroup(pop, pop.maxGroup, 0, 0, 0, Int[], Int[], Int[], Int[])
   else
      return _generate_group(pop, sires=sires, dams=dams, aliveonly=aliveonly)
   end
end

function _generate_group(pop::PTPopulation; sires::Vector{Int}=Int[], dams::Vector{Int}=Int[], aliveonly::Bool=true)
   check_idlist(sires)
   check_idlist(dams)

   if aliveonly
      sirelist = sires[map(x->pop.df[x,:alive],sires)]
      damlist = dams[map(x->pop.df[x,:alive],dams)]
   else
      sirelist = sires
      damlist = dams
   end

   nsires = length(sirelist)
   ndams = length(damlist)
   hassire = ifelse(nsires>0, true, false)
   hasdam = ifelse(ndams>0, true, false)
   if !hassire && !hasdam
      throw(ArgumentError("no sire and dam lists; use generate_group(pop,empty=true) for an empty group"))
   end
  
   n = nsires + ndams
   pop.maxGroup = pop.maxGroup + 1
   idlist = [sirelist; damlist]
   sort!(idlist)
   generation = zeros(Int,n)
   return PTGroup(pop, pop.maxGroup, n, nsires, ndams, sirelist, damlist, idlist, generation)
end

function check_idlist(idlist::Vector{Int})
   n = length(idlist)
   if n>0 && n!=length(unique(idlist))
      throw(ArgumentError("possible duplicated animals in the list"))
   end
   if sum(idlist .< 1)>0
      throw(ArgumentError("zero or negative ID(s) in the list"))
   end
   return nothing
end

"""
    group = copy_group(orig::PTGroup)

Make a deep copy of the griginal group `orig`.
A new group code (gid) will be assigned.
"""
function copy_group(orig::PTGroup)
   pop = orig.pop
   pop.maxGroup = pop.maxGroup + 1
   n = length(orig.sires) + length(orig.dams)
   sires = copy(orig.sires)
   dams = copy(orig.dams)
   id = copy(orig.id)
   generation = copy(orig.generation)
   return PTGroup(pop, pop.maxGroup, n, orig.maxSire, orig.maxDam, sires, dams, id, generation)
end

# a general function to select IDs for a filter function
"""
    id = selectid(fun, pop, ...)
    id = selectid(fun, group, ...)
    id = selectid(fun, groups, ...)

Returns an array of IDs that meet the condition defined by a function `f` in a population, a group, or an array of groups.
The query function `f` should be reasonable for the dataframe stored in `pop`.
Some options are available.

- `aliveonly`: select only individuals alive (default=`true`).
- `limit`: maximum length of output (default=`nothing` i.e., all IDs)
- `idlist`: select IDs only in this list (default=`Int[]` i.e., no constrant).
- `allowempty`: allow an empty ID list (default=`true`); for `false`, the function throws an error if `length(id)==0`.
- `sortby`: sort IDs by items in the dataframe or particular condition
    - `:id`: no sorting
    - `:random`: randomized array
    - any comparable quantities (e.g., `:year`, `:tbv`, `:pbv`, `:qbv`, `:ebv`, `:gebv`, etc.)
- `rev`: revserse sort (default=`true`, i.e., descending "high to low" order)

When an array of groups is supplied, `aliveonly`, `idlist`, and `allowempty` will be aplied to each group; `limit`, `sortby` and `rev` will be applied only to the final output.

```julia
# males in the population
selectid(:male => x -> x==true, pop)

# culled males with unknown sire (coded by 0)
selectid([:male,:sire] => (x,y) -> x==true && y==0, group, aliveonly=false)

# top males with highest EBV in the group
selectid([:male] => x -> x==true, group, sortby=:ebv)
```
"""
function selectid(fun, pop::PTPopulation; limit::Union{Nothing,Int}=nothing, idlist::Vector{Int}=Int[], aliveonly::Bool=true,
   sortby::Symbol=:id, rev::Bool=true, allowempty::Bool=true)

   # symbols
   symb = get_symbol_array(fun)
   if !in(:alive,symb); push!(symb,:alive); end
   if !in(:id,symb);    push!(symb,:id);    end

   selected_id = _selectid(fun, symb, pop.df, idlist=idlist, aliveonly=aliveonly, sortby=sortby, rev=rev, allowempty=allowempty)

   leng = valid_length(selected_id, limit)
   if leng==length(selected_id)
      return selected_id
   else
      return selected_id[1:leng]
   end
end

# for groups
function selectid(fun, group::PTGroup; limit::Union{Nothing,Int}=nothing, idlist::Vector{Int}=Int[], aliveonly::Bool=true,
   sortby::Symbol=:id, rev::Bool=true, allowempty::Bool=true)

   # symbols
   symb = get_symbol_array(fun)
   if !in(:alive,symb); push!(symb,:alive); end
   if !in(:id,symb);    push!(symb,:id);    end

   if isempty(idlist)
      selected_id = _selectid(fun, symb, group.pop.df, idlist=group.id, aliveonly=aliveonly, sortby=sortby, rev=rev, allowempty=allowempty)
   else
      unique_id = sort(unique(intersect(idlist,group.id)))
      selected_id = _selectid(fun, symb, group.pop.df, idlist=unique_id, aliveonly=aliveonly, sortby=sortby, rev=rev, allowempty=allowempty)
   end

   leng = valid_length(selected_id, limit)
   if leng==length(selected_id)
      return selected_id
   else
      return selected_id[1:leng]
   end
end

function selectid(fun, groups::Vector{PTGroup}; limit::Union{Nothing,Int}=nothing, idlist::Vector{Int}=Int[], aliveonly::Bool=true, 
   sortby::Symbol=:id, rev::Bool=true, allowempty::Bool=true)

   if length(groups)==0
      return Int[]
   end

   # symbols
   symb = get_symbol_array(fun)
   if !in(:alive,symb); push!(symb,:alive); end
   if !in(:id,symb);    push!(symb,:id);    end

   selected_id = Int[]
   for group in groups
      if isempty(idlist)
         selected_id = [selected_id; _selectid(fun, symb, group.pop.df, idlist=group.id, aliveonly=aliveonly, allowempty=allowempty)]
      else
         unique_id = sort(unique(intersect(idlist,group.id)))
         selected_id = [selected_id; _selectid(fun, symb, group.pop.df, idlist=unique_id, aliveonly=aliveonly, allowempty=allowempty)]
      end   
   end

   selected_id = unique(selected_id)

   # sort
   if sortby==:id
      ret_id = selected_id
   elseif sortby==:random
      ret_id = shuffle(selected_id)
   else
      df = groups[1].pop.df
      perm = sortperm(df[selected_id,sortby],rev=rev)
      ret_id = selected_id[perm]
   end
   # check
   if !allowempty && isempty(ret_id)
      error("empty ID list generated by selectid")
   end
   # length check
   leng = valid_length(ret_id, limit)
   if leng==length(ret_id)
      return ret_id
   else
      return ret_id[1:leng]
   end
end

function valid_length(idlist::Vector{Int},limit::Union{Nothing,Int})
   listlen = length(idlist)
   if isnothing(limit)
      return listlen
   else
      if limit<0
         throw(ArgumentError("invalid limit < 0"))
      else
         return min(limit,listlen)
      end
   end
end

function selectid_general(fun, group::PTGroup; idlist::Vector{Int}=Int[], aliveonly::Bool=true, sortby::Symbol=:id, rev::Bool=true, allowempty::Bool=true)
   # symbols
   symb = get_symbol_array(fun)
   if !in(:alive,symb); push!(symb,:alive); end
   if !in(:id,symb);    push!(symb,:id);    end
   
   # including :generation
   if in(:generation,symb)
      if length(idlist)>0
         unique_id = sort(unique(intersect(idlist,group.id)))
         df = pop.df[unique_id,symb]
         df[:,:generation] = group.generation[unique_id]
      else
         df = pop.df[group.id,symb]
         df[:,:generation] = group.generation
      end
   else
      if length(idlist)>0
         unique_id = sort(unique(intersect(idlist,group.id)))
         df = group.pop.df[unique_id,symb]
      else
         df = group.pop.df[group.id,symb]
      end
   end

   return _selectid(fun, symb, df, aliveonly=aliveonly, sortby=sortby, rev=rev, allowempty=allowempty)
end

function _selectid(fun, symb::Vector{Symbol}, df::DataFrame; idlist::Vector{Int}=Int[], aliveonly::Bool=true, sortby::Symbol=:id, rev::Bool=true, allowempty::Bool=true)
   # filtering
   if aliveonly
      if length(idlist)>0
         alive_id = unique(idlist[findall(df[idlist,:alive])])
      else
         alive_id = df[df[!,:alive],:id]
      end
      selected_id = collect( filter(fun, df[alive_id,symb],view=true)[!,:id] )
   else
      if length(idlist)>0
         unique_id = unique(idlist)
         selected_id = collect( filter(fun, df[unique_id,symb],view=true)[!,:id] )
      else
         selected_id = collect( filter(fun, df[:,symb],view=true)[!,:id] )
      end
   end
   # sort
   if sortby==:id
      ret_id = selected_id
   elseif sortby==:random
      ret_id = shuffle(selected_id)
   else
      perm = sortperm(df[selected_id,sortby],rev=rev)
      ret_id = selected_id[perm]
   end
   # check
   if !allowempty && isempty(ret_id)
      error("empty ID list generated by selectid")
   end
   return ret_id
end

function get_symbol_array(fun)
   if typeof(fun[1])==Symbol
      symb = [fun[1]]
   else
      symb = copy(fun[1])
   end
   return symb
end

# random sampling
function random_sampling(pop::PTPopulation, n::Int; male::Bool=false, female::Bool=false, aliveonly::Bool=true, allowempty::Bool=true)
   if male && female
      idlist = selectid([:id] => x->x>0, pop, aliveonly=aliveonly, sortby=:random, allowempty=allowempty)
   elseif male
      idlist = selectid([:male] => x->x==true, pop, aliveonly=aliveonly, sortby=:random, allowempty=allowempty)
   elseif female
      idlist = selectid([:male] => x->x==!true, pop, aliveonly=aliveonly, sortby=:random, allowempty=allowempty)
   else
      throw(ArgumentError("needs male=true and/or female=true"))
   end
   if length(idlist)>=n
      return idlist[1:n]
   elseif allowempty
      return idlist
   else
      error("too few sampled individuals")
   end
end

function random_sampling(pop::PTPopulation; nm::Int=0, nf::Int=0, aliveonly::Bool=true, allowempty::Bool=true)
   if nm>0
      mlist = selectid([:male] => x->x==true, pop, aliveonly=aliveonly, sortby=:random, allowempty=allowempty)
      if length(mlist)>=nm
         return mlist[1:nm]
      elseif allowempty
         return mlist
      else
         error("too few sampled individuals")
      end
   else
      mlist = Int[]
   end
   if nf>0
      flist = selectid([:male] => x->x==!true, pop, aliveonly=aliveonly, sortby=:random, allowempty=allowempty)
      if length(flist)>=nf
         return flist[1:nf]
      elseif allowempty
         return flist
      else
         error("too few sampled individuals")
      end
   else
      flist = Int[]
   end
   idlist = [mlist; flist]
   return idlist
end

"""
    n = add_sires!(group::PTGroup, idlist::Vector{Int})

Add animals to a sire list in `group`.
The function reads the candidate list from the top, then adds the animal if it is not in the list.
If the size of list reaches the maximum (defined by `maxSire` in the group structure), this function stops adding any more.
It returns how many sires have been put into the sire list in the group.
"""
function add_sires!(group::PTGroup, idlist::Vector{Int}; verbose::Bool=false)
   n = length(group.sires)
   ntaken = 0
   @inbounds for i=1:length(idlist)
      candidate = idlist[i]
      if n>=group.maxSire
         break
      end
      if candidate>0
         if !group.pop.df[candidate,:alive]
            if verbose
               @warn "culled sire $(candidate); skip."
            end
            continue
         end
         if !in(candidate,group.sires) && group.pop.df[candidate,:male]
            push!(group.sires,candidate)
            n = n + 1
            ntaken = ntaken + 1
            if !in(candidate,group.id)
               push!(group.id,candidate)
               push!(group.generation,0)
            end
         end
      end
   end
   group.n = length(group.id)
   return ntaken
end

function add_sires!(group::PTGroup, id::Int; verbose::Bool=false)
   return add_sires!(group, [id], verbose=verbose)
end

"""
    n = add_dams!(group::PTGroup, idlist::Vector{Int})

Add animals to a sire list in `group`.
The function reads the candidate list from the top, then adds the animal if it is not in the list.
If the size of list reaches the maximum (defined by `maxDam` in the group structure), this function stops adding any more.
It returns how many dams have been taken.
"""
function add_dams!(group::PTGroup, idlist::Vector{Int}; verbose::Bool=false)
   n = length(group.dams)
   ntaken = 0
   @inbounds for i=1:length(idlist)
      candidate = idlist[i]
      if n>=group.maxDam
         break
      end
      if candidate>0
         if !group.pop.df[candidate,:alive]
            if verbose
               @warn "culled dam $(candidate); skip."
            end
            continue
         end
         if !in(candidate,group.dams) && !group.pop.df[candidate,:male]
            push!(group.dams,candidate)
            n = n + 1
            ntaken = ntaken + 1
            if !in(candidate,group.id)
               push!(group.id,candidate)
               push!(group.generation,0)
            end
         end
      end
   end
   group.n = length(group.id)
   return ntaken
end

function add_dams!(group::PTGroup, id::Int; verbose::Bool=false)
   return add_dams!(group, [id], verbose=verbose)
end

"""
    n = mark_candidate!(pop::PTPopulation, idlist::Vector{Int})

Mark individuals in `idlist` as candidates for selection.
It just puts a flag in the dataframe, and does nothing more.
"""
function mark_candidate!(pop::PTPopulation,idlist::Vector{Int}; check::Bool=true)
   for id in idlist
      if check
         if id<=0 || id>pop.maxAnimal || !pop.df[id,:alive]
            error("id $(id) out of range, dead, or inappropriate")
         end
      end
      pop.df[id,:candidate] = true
   end
   return nothing
end

"""
    n = change_status!(pop::PTPopulation, idlist::Vector{Int}, status=0)

Change status for individuals in `idlist`.
It just puts a number in the dataframe, and does nothing more.
The status code is arbitrary; the user can define it (the default code is 0).
This function does not check if the individual is alive or not.
"""
function change_status!(pop::PTPopulation,idlist::Vector{Int},status::Int=0; check::Bool=true)
   for id in idlist
      if check
         if id<=0 || id>pop.maxAnimal
            error("id $(id) out of range")
         end
      end
      pop.df[id,:status] = status
   end
   return nothing
end

"""
    n = vacancy_for_sires(group)

The number of vacant sires to be filled.
"""
function vacancy_for_sires(group)
   return group.maxSire - length(group.sires)
end

"""
    n = vacancy_for_dams(group)

The number of vacant dams to be filled.
"""
function vacancy_for_dams(group)
   return group.maxDam - length(group.dams)
end

"""
    n = cull!(pop::PTPopulation, idlist::Vector{Int})
    n = cull!(group::PTGroup, idlist::Vector{Int})

Remove the animal IDs in `idlist` from population `pop` or group `group`.
Returns the number of culled animals, `n`.
With `group`, the animals will also be removed from the population.
The status of the animals becomes "dead", i.e., `alive=false`.
"""
function cull!(pop::PTPopulation,idlist::Vector{Int})
   nculled = 0
   @inbounds for i in idlist
      if 0<i && i<=pop.maxAnimal && pop.df[i,:alive]
         pop.df[i,:alive] = false
         nculled = nculled + 1
      end
   end
   return nculled
end

function cull!(group::PTGroup,idlist::Vector{Int})
   pop = group.pop
   nculled = cull!(pop,idlist)
   setdiff!(group.sires,idlist)
   setdiff!(group.dams,idlist)
   return nculled
end

function cull!(groups::Vector{PTGroup},idlist::Vector{Int})
   nculled = 0
   for group in groups
      nculled = nculled + cull!(group,idlist)
   end
   return nculled
end
