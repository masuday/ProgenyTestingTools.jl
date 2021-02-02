# BLUP test
using ProgenyTestingTools
using DelimitedFiles
using SparseArrays

# mu, var_p, h2_poly, h2_qtl, rep
par = PTParameters(50, 100, 0.36, 0, 0, 1.0)
hp = generate_population(par,nm=20,nf=20)
pop = generate_population(par)
bulls1 = migrate_from_hp!(hp,pop,random_sampling(hp,2,male=true))
cows1 = migrate_from_hp!(hp,pop,random_sampling(hp,3,female=true))
mgroup1 = generate_group(pop, sires=bulls1)
fgroup1 = generate_group(pop, dams=cows1)
mating!(mgroup1, fgroup1, "fgroup", method="dairy_standard_ai")
assign_phenotype!(fgroup1, gen=1, repeated=true)
calving!(fgroup1)

fgroup1.dams = fgroup1.id[fgroup1.generation .== 1]
mating!(mgroup1, fgroup1, "fgroup", method="dairy_standard_ai")
assign_phenotype!(fgroup1, gen=2, repeated=true)
calving!(fgroup1)

fgroup1.dams = fgroup1.id[fgroup1.generation .== 2]
mating!(mgroup1, fgroup1, "fgroup", method="dairy_standard_ai")
assign_phenotype!(fgroup1, gen=3, repeated=true)
calving!(fgroup1)

fgroup1.dams = fgroup1.id[fgroup1.generation .== 3]
mating!(mgroup1, fgroup1, "fgroup", method="dairy_standard_ai")
assign_phenotype!(fgroup1, gen=4, repeated=true)
calving!(fgroup1)

write_files_for_blup(pop,"phenotype.txt","pedigree.txt", repeated=false)

(lhs,rhs) = ProgenyTestingTools.build_mme(pop,"phenotype.txt","pedigree.txt",false,false)
sol = lhs \ rhs

open("solutions.txt", "w") do io
   writedlm(io, sol)
end

# repeatability model
# mu, var_p, h2_poly, h2_qtl, rep
par = PTParameters(50, 100, 0.36, 0, 0.45, 1.0)
hp = generate_population(par,nm=20,nf=20)
pop = generate_population(par)
bulls1 = migrate_from_hp!(hp,pop,random_sampling(hp,2,male=true))
cows1 = migrate_from_hp!(hp,pop,random_sampling(hp,3,female=true))
mgroup1 = generate_group(pop, sires=bulls1)
fgroup1 = generate_group(pop, dams=cows1)
mating!(mgroup1, fgroup1, "fgroup", method="dairy_standard_ai")
assign_phenotype!(fgroup1, gen=1, repeated=true)
assign_phenotype!(fgroup1, gen=1, repeated=true)
calving!(fgroup1)

fgroup1.dams = fgroup1.id[fgroup1.generation .== 1]
mating!(mgroup1, fgroup1, "fgroup", method="dairy_standard_ai")
assign_phenotype!(fgroup1, gen=2, repeated=true)
assign_phenotype!(fgroup1, gen=2, repeated=true)
calving!(fgroup1)

fgroup1.dams = fgroup1.id[fgroup1.generation .== 2]
mating!(mgroup1, fgroup1, "fgroup", method="dairy_standard_ai")
assign_phenotype!(fgroup1, gen=3, repeated=true)
assign_phenotype!(fgroup1, gen=3, repeated=true)
calving!(fgroup1)

fgroup1.dams = fgroup1.id[fgroup1.generation .== 3]
mating!(mgroup1, fgroup1, "fgroup", method="dairy_standard_ai")
assign_phenotype!(fgroup1, gen=4, repeated=true)
assign_phenotype!(fgroup1, gen=4, repeated=true)
calving!(fgroup1)

write_files_for_blup(pop,"phenotype2.txt","pedigree2.txt", repeated=true)

(lhs,rhs) = ProgenyTestingTools.build_mme(pop,"phenotype2.txt","pedigree2.txt",true,false, verbose=true)
sol = lhs \ rhs

open("solutions2.txt", "w") do io
   writedlm(io, sol)
end
