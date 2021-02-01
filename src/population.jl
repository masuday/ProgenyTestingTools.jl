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
hp = generate_population(par, 5, 10)

# hostorical population with genotypes
hp = generate_population(par, "qmsimdata.h5", 5, 10)

# new population with no animals
pop = generate_population(par)

# new population with genotypes
pop = generate_population(par, "qmsimdata.h5", 0, 0, gfile="pop.h5")
```
"""
function generate_population(par::PTParameters; nm::Int=0, nf::Int=0)
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

   df = DataFrame(
      id = collect(1:maxAnimal),
      male = vcat([true for i=1:nm],[false for i=1:nf]),
      year = zeros(Int, maxAnimal),
      alive = fill(true, maxAnimal),
      pregnant = fill(false, maxAnimal),
      genotyped = fill(false, maxAnimal),
      sire = zeros(Int, maxAnimal),
      dam = zeros(Int, maxAnimal),
      siregroup = zeros(Int, maxAnimal),
      damgroup = zeros(Int, maxAnimal),
      nprog = zeros(Int, maxAnimal),
      nrecprog = zeros(Int, maxAnimal),
      inb = zeros(Float64,maxAnimal),
      pbv = zeros(Float64,maxAnimal),
      qbv = zeros(Float64,maxAnimal),
      tbv = zeros(Float64,maxAnimal),
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

   return PTPopulation(par,hp,maxAnimal,maxGroup,maxCG,df,animal) #,map,h5file)
end

"""
    add_animal!(pop::PTPopulation)

Add an animal to population `pop`.
The added animals lose the pedigree information.
"""
function add_animal!(pop::PTPopulation, male::Bool, animal::PTAnimal; 
                     genotyped=false, year=0, alive=true, pregnant=false, sire=0, dam=0, siregroup=0, damgroup=0, inb=0.0, pbv=0.0, qbv=0.0, tbv=0.0, ebv=0.0, rel=missing, gebv=missing)
   pop.maxAnimal = pop.maxAnimal + 1
   (nrec,firsty,lasty,avgy) = get_phenotypic_information(animal.y)
   newdata = (pop.maxAnimal, male, year, alive, pregnant, genotyped, sire, dam, siregroup, damgroup, 0, 0, inb, pbv, qbv, tbv, ebv, gebv, rel, nrec, firsty, lasty, avgy)
   push!(pop.df, newdata)
   push!(pop.animal, animal)
   return nothing
end

function add_animal!(pop::PTPopulation, row::DataFrameRow{DataFrames.DataFrame,DataFrames.Index}, animal::PTAnimal)
   pop.maxAnimal = pop.maxAnimal + 1
   (nrec,firsty,lasty,avgy) = get_phenotypic_information(animal.y)
   push!(pop.df, row)
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
    newidlist = migrate_from_hp!(hp::PTPopulation, pop::PTPopulation, idlist::Vector{Int}; year=0)

Migrate some animals from a historical population `hp` to another population `pop`.
The migrated animals have new IDs in the new population.
The animals will have "dead" status in the historical population.
It returns a new ID list `newidlist` in the new population.
"""
function migrate_from_hp!(hp::PTPopulation, pop::PTPopulation, idlist::Vector{Int}; year::Union{Int,Vector{Int}}=0)
   unique_idlist = unique(idlist)
   n = length(unique_idlist)
   newidlist = Vector{Int}()
   @inbounds for id in unique_idlist
      # check id in the historiacl population
      if id > hp.maxAnimal
         error("No such ID $(id) in the first population")
      end
      if hp.df[id,:alive]
         # new animal for the new population
         add_animal!(pop, hp.df[id,:], hp.animal[id])
         pop.df[pop.maxAnimal,:id] = pop.maxAnimal
         # moved
         hp.df[id,:alive] = false
         # assigned new id
         push!(newidlist, pop.maxAnimal)
      else
         @warn "id $(id) in hp is not alive, and has not been added to the populatiopn."
      end

      # (hp) id to (pop) pop.maxAnimal
      # assigned new genotypes
      #if hp.genotyped[id]
      #   # copy genotypes from the database
      #   g = read_qmsim_individual_hdf5(hp.map,hp.gfile,id)
      #   add_qmsim_individual_hdf5(pop.map,pop.gfile,g)
      #   pop.genotyped[pop.maxAnimal] = true
      #   pop.qbv[pop.maxAnimal] = g.tbv
      #   pop.tbv[pop.maxAnimal] = pop.pbv[pop.maxAnimal] + pop.qbv[pop.maxAnimal]
      #end
   end
   # assign year
   assign_year!(pop, unique_idlist, year)
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
function generate_group(pop::PTPopulation; sires::Vector{Int}=Int[], dams::Vector{Int}=Int[])
   nsires = length(sires)
   ndams = length(dams)
   hassire = ifelse(nsires>0, true, false)
   hasdam = ifelse(length(dams)>0, true, false)
   if !hassire && !hasdam
      throw(ArgumentError("no sire and dam lists"))
   end
   check_idlist(sires)
   check_idlist(dams)
   n = nsires + ndams
   pop.maxGroup = pop.maxGroup + 1
   idlist = [sires; dams]
   sort!(idlist)
   generation = zeros(Int,n)
   return PTGroup(pop, pop.maxGroup, n, nsires, ndams, sires, dams, idlist, generation)
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
# example 1: selectid(:male => x -> x==true, hp, aliveonly=true)
# example 2: selectid([:male,:sire] => (x,y) -> x==true && y==0, hp, aliveonly=true)
function selectid(fun, pop::PTPopulation; idlist::Vector{Int}=Int[], aliveonly::Bool=true, sortby::Symbol=:id, rev::Bool=true, allowempty::Bool=true)
   # symbols
   symb = get_symbol_array(fun)
   if !in(:alive,symb); push!(symb,:alive); end
   if !in(:id,symb);    push!(symb,:id);    end

   return _selectid(fun, symb, pop.df, idlist=idlist, aliveonly=aliveonly, sortby=sortby, rev=rev, allowempty=allowempty)
end

# for groups
function selectid(fun, group::PTGroup; aliveonly::Bool=true, sortby::Symbol=:id, rev::Bool=true, allowempty::Bool=true)
   # symbols
   symb = get_symbol_array(fun)
   if !in(:alive,symb); push!(symb,:alive); end
   if !in(:id,symb);    push!(symb,:id);    end

   return _selectid(fun, symb, group.pop.df, idlist=group.id, aliveonly=aliveonly, sortby=sortby, rev=rev, allowempty=allowempty)
end

function selectid_general(fun, group::PTGroup; idlist::Vector{Int}=Int[], aliveonly::Bool=true, sortby::Symbol=:id, rev::Bool=true, allowempty::Bool=true)
   # symbols
   symb = get_symbol_array(fun)
   if !in(:alive,symb); push!(symb,:alive); end
   if !in(:id,symb);    push!(symb,:id);    end
   
   # including :generation
   if in(:generation,symb)
      if length(idlist)>0
         unique_id = sort(unique(intersection(idlist,group.id)))
         df = pop.df[unique_id,symb]
         df[:,:generation] = group.generation[unique_id]
      else
         df = pop.df[group.id,symb]
         df[:,:generation] = group.generation
      end
   else
      if length(idlist)>0
         unique_id = sort(unique(intersection(idlist,group.id)))
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
      ret_id = randperm(selected_id)
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
