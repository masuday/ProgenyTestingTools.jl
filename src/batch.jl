# utility functions to run a batch of jobs

# assign the young bulls to random herd
function test_mating!(young_stud, herds, birthyear, n_test_mating_per_bull; pmp=0.5, method="random", updateinb=false, prob::Union{Nothing,Vector{Float64}}=nothing)
   n_herds = length(herds)
   n_pregnant = 0
   if method=="random"
      for bull in young_stud.sires
         # n mating per young bull
         for i in 1:n_test_mating_per_bull
            h = rand(1:n_herds)
            preg_dams = mating!(young_stud, herds[h], "fgroup", year=birthyear, method="dairy_standard_ai", n=1, pmp=pmp, updateinb=updateinb, prob=prob);
            n_pregnant = n_pregnant + length(preg_dams)
         end
      end
   else
      throw(ArgumentError("use `random` for mating method"))
   end
   return n_pregnant
end

function regular_mating!(bull_stud, herds, birthyear; method="random", updateinb=false, prob::Union{Nothing,Vector{Float64}}=nothing)
   n_herds = length(herds)
   n_pregnant = 0
   if method=="random"
      for h in 1:n_herds
         preg_dams = mating!(bull_stud, herds[h], "fgroup", year=birthyear, method="dairy_standard_ai", updateinb=updateinb, prob=prob) 
         n_pregnant = n_pregnant + length(preg_dams)
      end
   else
      throw(ArgumentError("use `random` for mating method"))
   end
   return n_pregnant   
end

function cull_old_bulls!(bull_stud,birth_year_of_culled_bulls)
   return cull!(bull_stud, selectid([:male,:year] => (x,y) -> x==true && y==birth_year_of_culled_bulls, bull_stud))
end

function cull_old_cows!(herds,birth_year_of_culled_cows)
   n_herds = length(herds)
   nculled = 0
   for h=1:n_herds
      nculled = nculled + cull!(herds[h], selectid([:male,:year] => (x,y) -> x==false && y==birth_year_of_culled_cows, herds[h]))
   end
   return nculled
end

function mean_inbreeding(pop; year::Int=-1)
   if year>=0
      return mean(pop.df[pop.df[:,:year].==year, :inb])
   else
      return mean(pop.df[:, :inb])
   end
end

function number_of_active_females(herd::PTGroup)
   females = selectid([:male]=>x->x==false, herd, aliveonly=true)
   nheifers = length(setdiff(females, herd.dams))
   ndams = length(herd.dams)
   return ndams,nheifers
end

function number_of_active_females(herds::Vector{PTGroup})
   n_herds = length(herds)
   ndams = 0
   nheifers = 0
   for h=1:n_herds
      (nd,nh) = number_of_active_females(herds[h])
      nheifers = nheifers + nh
      ndams = ndams + nd
   end
   return ndams,nheifers
end

function phenotyped_cows(herds::Vector{PTGroup}; year=-1, test_status=0, plan="limited", n_breeders=0, n_daughters=-1, verbose=false)
   n_herds = length(herds)
   pop = herds[1].pop
   if plan=="limited"
      # rule: All tested cows have phenotypes.
      # rule: All cows in breeder herds have phenotypes.
      # rule: Each bull has the upper limit on the number of daughters each year.
      tested_id = zeros(Int,0)
      breeder_cow_id = zeros(Int,0)
      standard_cow_id = zeros(Int,0)
      phenotyped_cow_id = zeros(Int,0)
      if year>0
         # tested daughters
         tested_id = selectid([:male,:year,:status] => (x,y,z) -> x==false && y==year && z==test_status, pop)
         # regular daughters
         for i=1:n_breeders
            breeder_cow_id = [breeder_cow_id; selectid([:male,:year] => (x,y) -> x==false && y==year, herds[i])]
         end
         for i=n_breeders+1:n_herds
            standard_cow_id = [standard_cow_id; selectid([:male,:year] => (x,y) -> x==false && y==year, herds[i])]
         end
         n_tested = length(tested_id)
         n_breeder = length(breeder_cow_id)
         n_standard = length(standard_cow_id)
         if n_daughters<0
            n_additional = n_standard
         else
            n_additional = min(max(0,n_daughters-n_breeder),n_standard)
         end

         if verbose
            println("  phenotyped cows:")
            println("    tested:$(n_tested) breeder:$(n_breeder) standard:$(n_additional); total:$(n_daughters)")
         end

         # tested daughters + breeder cows + additional cows
         phenotyped_cow_id = [phenotyped_cow_id; tested_id; breeder_cow_id]
         if n_additional>0
            shuffle!(standard_cow_id)
            phenotyped_cow_id = [phenotyped_cow_id; standard_cow_id[1:n_additional]]
         end
         sort!(phenotyped_cow_id)
         return phenotyped_cow_id
      else
         throw(ArgumentError("year not given"))
      end
   else
      throw(ArgumentError("plan not supported: $(plan)"))
   end
   return zeros(Int,0)
end

# To get N daughters on average with SD
# possible use:
#   Gamma(0.8,1.0); scale=65, shift=50; mean=100 std=55 upper=200 (lower=50)
#   Gamma(1.0,1.0); scale=40, shift=50; mean=90  std=40 upper=250 (lower=50)
#   or Gamma(1.0,1.0) upper=4.0 (98% percentile)
function expected_frequency_of_mating(mgroup::PTGroup,dist; by=:ebv, rev=true, scale=1.0, shift=0.0, upper=0.0)
   pop = mgroup.pop
   n = length(mgroup.sires)
   benchmark = pop.df[mgroup.sires,by]
   perm = sortperm(benchmark,rev=rev) # large to small by default
   iperm = invperm(perm)
   x = sort((rand(dist,n)*scale .+ shift),rev=true)  # large to small
   if upper>0
      x[x .> upper] .= upper
   end
   @show x
   x .= x[iperm]/sum(x)
   return x
end
