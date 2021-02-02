# utility functions to run a batch of jobs

# assign the young bulls to random herd
function test_mating!(young_stud, herds, birthyear, n_test_mating_per_bull; pmp=0.5, method="random", updateinb=false)
   n_herds = length(herds)
   n_pregnant = 0
   if method=="random"
      for bull in young_stud.sires
         # n mating per young bull
         for i in 1:n_test_mating_per_bull
            h = rand(1:n_herds)
            preg_dams = mating!(young_stud, herds[h], "fgroup", year=birthyear, method="dairy_standard_ai", n=1, pmp=pmp, updateinb=updateinb);
            n_pregnant = n_pregnant + length(preg_dams)
         end
      end
   else
      throw(ArgumentError("use `random` for mating method"))
   end
   return n_pregnant
end

function regular_mating!(bull_stud, herds, birthyear; method="random", updateinb=false)
   n_herds = length(herds)
   n_pregnant = 0
   if method=="random"
      for h in 1:n_herds
         preg_dams = mating!(bull_stud, herds[h], "fgroup", year=birthyear, method="dairy_standard_ai", updateinb=updateinb) 
         n_pregnant = n_pregnant + length(preg_dams)
      end
   else
      throw(ArgumentError("use `random` for mating method"))
   end
   return n_pregnant   
end
