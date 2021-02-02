using ProgenyTestingTools
using Test
using Random

@testset "parameters" begin
   # no repeatability
   var_a = 6.0
   var_g = 30.0
   var_u = var_a + var_g
   var_pe = 0.0
   var_e = 64.0

   # pissibility 1: rep = 0
   var_p = var_u + var_pe + var_e
   h2_poly = var_a/var_p
   h2_qtl = var_g/var_p
   rep = (var_u+var_pe)/var_p
   par = PTParameters(50, var_p, h2_poly, h2_qtl, rep, 1.0)
   @test isnothing( check_parameters(par) )
   @test get_var_poly(par) ≈ var_a
   @test get_var_qtl(par) ≈ var_g
   @test get_var_pe(par) ≈ var_pe
   @test get_var_error(par) ≈ var_e

   # pissibility 2: rep = heritability
   rep = h2_poly + h2_qtl
   par = PTParameters(50, var_p, h2_poly, h2_qtl, rep, 1.0)
   @test isnothing( check_parameters(par) )
   @test get_var_poly(par) ≈ var_a
   @test get_var_qtl(par) ≈ var_g
   @test get_var_pe(par) ≈ var_pe
   @test get_var_error(par) ≈ var_e

   # repeated records
   var_a = 6.0
   var_g = 30.0
   var_u = var_a + var_g
   var_pe = 9.0
   var_e = 55.0
   var_p = var_u + var_pe + var_e
   h2_poly = var_a/var_p
   h2_qtl = var_g/var_p
   rep = (var_u+var_pe)/var_p
   par = PTParameters(50, var_p, h2_poly, h2_qtl, rep, 1.0)
   par = PTParameters(50, var_p, h2_poly, h2_qtl, rep, 1.0)
   @test isnothing( check_parameters(par) )
   @test get_var_poly(par) ≈ var_a
   @test get_var_qtl(par) ≈ var_g
   @test get_var_pe(par) ≈ var_pe
   @test get_var_error(par) ≈ var_e

   # out of range
   @test_throws ErrorException check_parameters(PTParameters(50, -100, 0.1, 0.1, 0.5, 1.0))
   @test_throws ErrorException check_parameters(PTParameters(50, 100, -0.1, 0.1, 0.5, 1.0))
   @test_throws ErrorException check_parameters(PTParameters(50, 100, 1.1, 0.1, 0.5, 1.0))
   @test_throws ErrorException check_parameters(PTParameters(50, -100, 0.1, -0.1, 0.5, 1.0))
   @test_throws ErrorException check_parameters(PTParameters(50, -100, 0.1, 1.1, 0.5, 1.0))
   @test_throws ErrorException check_parameters(PTParameters(50, -100, 0.1, 0.1, -0.1, 1.0))
   @test_throws ErrorException check_parameters(PTParameters(50, -100, 0.1, 0.1, 0.1, 1.0))
   @test_throws ErrorException check_parameters(PTParameters(50, -100, 0.1, 0.1, 1.1, 1.0))
end

@testset "generate_population" begin
   par = PTParameters(50, 100, 0.5, 0, 0, 1.0)
   hp = generate_population(par,nm=5,nf=5)
   @test hp.maxAnimal == 10
   @test sum(hp.df[!,:male] .== true) == 5
   @test sum(hp.df[!,:male] .== false) == 5
end

@testset "migrate_from_hp!" begin
   par = PTParameters(50, 100, 0.5, 0, 0, 1.0)
   hp = generate_population(par,nm=5,nf=5)
   pop = generate_population(par)
   idlist = [3,8]
   newidlist = migrate_from_hp!(hp,pop,idlist)
   @test pop.df[1,:male]
   @test !pop.df[2,:male]
   @test hp.df[3,:tbv] ≈ pop.df[1,:tbv]
   @test hp.df[8,:tbv] ≈ pop.df[2,:tbv]
   @test pop.df[!,:id] ≈ [1,2]
   @test pop.maxAnimal == 2

   # duplicated ID
   hp = generate_population(par,nm=5,nf=5)
   pop = generate_population(par)
   idlist = [3,8,8,3,8]
   newidlist = migrate_from_hp!(hp,pop,idlist)
   @test pop.df[1,:male]
   @test !pop.df[2,:male]
   @test hp.df[3,:tbv] ≈ pop.df[1,:tbv]
   @test hp.df[8,:tbv] ≈ pop.df[2,:tbv]
   @test pop.df[!,:id] ≈ [1,2]
   @test pop.maxAnimal == 2
end

@testset "random_sampling" begin
   par = PTParameters(50, 100, 0.5, 0, 0, 1.0)
   hp = generate_population(par,nm=10,nf=10)
   cull!(hp,[10,20])

   # male only
   idlist = random_sampling(hp, 10, male=true, aliveonly=true, allowempty=true)
   @test all( sort(idlist) .== [1,2,3,4,5,6,7,8,9] )
   @test_throws ErrorException random_sampling(hp, 10, male=true, aliveonly=true, allowempty=false)
   idlist = random_sampling(hp, 10, male=true, aliveonly=false, allowempty=true)
   @test all( sort(idlist) .== [1,2,3,4,5,6,7,8,9,10] )
   ntests = 100
   ntested = 0
   for i in 1:ntests
      idlist = random_sampling(hp, 5, male=true, aliveonly=true, allowempty=true)
      if sum(map(x->in(x,[1,2,3,4,5,6,7,8,9]),idlist))==5
         ntested = ntested + 1
      end
   end
   @test ntests == ntested

   # female only
   idlist = random_sampling(hp, 10, female=true, aliveonly=true, allowempty=true)
   @test all( sort(idlist) .== [11,12,13,14,15,16,17,18,19] )
   @test_throws ErrorException random_sampling(hp, 10, female=true, aliveonly=true, allowempty=false)
   idlist = random_sampling(hp, 10, female=true, aliveonly=false, allowempty=true)
   @test all( sort(idlist) .== [11,12,13,14,15,16,17,18,19,20] )
   ntests = 100
   ntested = 0
   for i in 1:ntests
      idlist = random_sampling(hp, 5, female=true, aliveonly=true, allowempty=true)
      if sum(map(x->in(x,[11,12,13,14,15,16,17,18,19]),idlist))==5
         ntested = ntested + 1
      end
   end
   @test ntests == ntested

   # male and female
   idlist = random_sampling(hp, 20, male=true, female=true, aliveonly=true, allowempty=true)
   @test all( sort(idlist) .== [1,2,3,4,5,6,7,8,9,11,12,13,14,15,16,17,18,19] )
   @test_throws ErrorException random_sampling(hp, 20, male=true, female=true, aliveonly=true, allowempty=false)
   idlist = random_sampling(hp, 20, male=true, female=true, aliveonly=false, allowempty=true)
   @test all( sort(idlist) .== [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20] )
   ntests = 100
   ntested = 0
   for i in 1:ntests
      idlist = random_sampling(hp, 15, male=true, female=true, aliveonly=true, allowempty=true)
      if sum(map(x->in(x,[1,2,3,4,5,6,7,8,9,11,12,13,14,15,16,17,18,19]),idlist))==15
         ntested = ntested + 1
      end
   end
   @test ntests == ntested
end

@testset "assign_year" begin
   par = PTParameters(50, 100, 0.5, 0, 0, 1.0)
   # test 1
   hp = generate_population(par,nm=5,nf=5)
   pop = generate_population(par)
   idlist = [3,5,8,9]
   newidlist = migrate_from_hp!(hp,pop,idlist)
   ProgenyTestingTools.assign_year!(pop,newidlist,[10,20])
   @test pop.df[1,:year] == 10
   @test pop.df[2,:year] == 20
   @test pop.df[3,:year] == 10
   @test pop.df[4,:year] == 20
   # test 2
   hp = generate_population(par,nm=5,nf=5)
   pop = generate_population(par)
   idlist = [3,5,8,9]
   newidlist = migrate_from_hp!(hp,pop,idlist)
   ProgenyTestingTools.assign_year!(pop,newidlist,10)
   @test pop.df[1,:year] == 10
   @test pop.df[2,:year] == 10
   @test pop.df[3,:year] == 10
   @test pop.df[4,:year] == 10
end

@testset "generate_group (from population)" begin
   par = PTParameters(50, 100, 0.5, 0, 0, 1.0)
   hp = generate_population(par,nm=5,nf=5)
   pop = generate_population(par)
   hp_males = [2,4,5]
   hp_females = [7,8,10]
   pop_bulls = migrate_from_hp!(hp,pop,hp_males)
   pop_dams = migrate_from_hp!(hp,pop,hp_females)
   group = generate_group(pop, sires=pop_bulls[1:2], dams=pop_dams[1:2])
   @test group.n == 4 && group.maxSire==2 && group.maxDam == 2
   @test all(group.sires .== [1,2]) && all(group.dams .== [4,5])
end

@testset "copy_group" begin
   par = PTParameters(50, 100, 0.5, 0, 0, 1.0)
   hp = generate_population(par,nm=5,nf=5)
   pop = generate_population(par)
   hp_males = [2,4,5]
   hp_females = [7,8,10]
   pop_bulls = migrate_from_hp!(hp,pop,hp_males)
   pop_dams = migrate_from_hp!(hp,pop,hp_females)
   group = generate_group(pop, sires=pop_bulls[1:2], dams=pop_dams[1:2])
   group.generation .= [1,2,3,4]
   c = copy_group(group)
   @test c.n==group.n && c.maxSire==group.maxSire && c.maxDam==group.maxDam && c.groupid==group.groupid+1 && 
        all(c.id .== group.id) && all(c.generation .== group.generation) &&
        all(c.sires .== group.sires) && all(c.sires .== group.sires)
end

@testset "selectid from population" begin
   par = PTParameters(50, 100, 0.5, 0, 0, 1.0)
   hp = generate_population(par,nm=5,nf=5)
   pop = generate_population(par)
   hp_males = [3,2,4,5]
   hp_females = [8,7,10,9]
   pop_bulls = migrate_from_hp!(hp,pop,hp_males)
   pop_dams = migrate_from_hp!(hp,pop,hp_females)
   pop.df[2,:alive] = false
   pop.df[7,:alive] = false
   pop.df[!,:tbv] = [2,4,3,1, 2,1,3,4] * 1.0
   pop.df[!,:gebv] = [2,4,3,missing, 2,1,missing,4] * 1.0
   # simple selection: default aliveonly=true
   id = selectid(:male => x -> x==true, pop)
   @test all( id .== [1,3,4] )
   id = selectid(:male => x -> x==true, pop, aliveonly=false)
   @test all( id .== [1,2,3,4] )
   id = selectid(:male => x -> x==false, pop)
   @test all( id .== [5,6,8] )
   id = selectid(:male => x -> x==false, pop, aliveonly=false)
   @test all( id .== [5,6,7,8] )
   # multiple conditions
   id = selectid([:male,:tbv] => (x,y) -> x==true && y>=3.0, pop)
   @test all( id .== [3] )
   id = selectid([:male,:tbv] => (x,y) -> x==true && y>=3.0, pop, aliveonly=false)
   @test all( id .== [2,3] )
   # missing values
   id = selectid([:male,:gebv] => (x,y) -> x==true && !ismissing(y) && y>=2.0, pop)
   @test all( id .== [1,3] )
   id = selectid([:male,:gebv] => (x,y) -> x==true && !ismissing(y) && y>=2.0, pop, aliveonly=false)
   @test all( id .== [1,2,3] )
   # allowempty
   id = selectid([:male,:gebv] => (x,y) -> x==true && !ismissing(y) && y>=5.0, pop, allowempty=true)
   @test all( id .== [] )
   @test_throws ErrorException selectid([:male,:gebv] => (x,y) -> x==true && !ismissing(y) && y>=5.0, pop, allowempty=false)
   # particular id
   id = selectid(:male => x -> x==true, pop, idlist=[1,4,6,8])
   @test all( id .== [1,4] )
   id = selectid(:male => x -> x==true, pop, idlist=[1,1,6,4,1,8,4])
   @test all( id .== [1,4] )
   # sort
   id = selectid(:male => x -> x==true, pop, sortby=:tbv)
   @test all( id .== [3,1,4] )
   id = selectid(:male => x -> x==true, pop, sortby=:tbv, rev=false)
   @test all( id .== [4,1,3] )
   id = selectid(:male => x -> x==true, pop, sortby=:tbv, aliveonly=false)
   @test all( id .== [2,3,1,4] )
   id = selectid(:male => x -> x==true, pop, sortby=:gebv)
   @test all( id .== [4,3,1] )
   id = selectid([:male,:gebv] => (x,y) -> x==true && !ismissing(y), pop, sortby=:gebv)
   @test all( id .== [3,1] )
   id = selectid([:male,:gebv] => (x,y) -> x==true && !ismissing(y) && y>=5.0, pop, sortby=:gebv, allowempty=true)
   @test all( id .== [] )
end

@testset "selectid from group" begin
   par = PTParameters(50, 100, 0.5, 0, 0, 1.0)
   hp = generate_population(par,nm=5,nf=5)
   pop = generate_population(par)
   hp_males = [3,2,4,5]
   hp_females = [8,7,10,9]
   pop_bulls = migrate_from_hp!(hp,pop,hp_males)
   pop_dams = migrate_from_hp!(hp,pop,hp_females)
   pop.df[2,:alive] = false
   pop.df[7,:alive] = false
   pop.df[!,:tbv] = [2,4,3,1, 2,1,3,4] * 1.0
   pop.df[!,:gebv] = [2,4,3,missing, 2,1,missing,4] * 1.0
   group = generate_group(pop, sires=pop_bulls[1:3], dams=pop_dams[1:3], aliveonly=false)
   group.generation = collect(1:length(group.generation))*10
   # simple selection: default aliveonly=true
   id = selectid(:male => x -> x==true, group)
   @test all( id .== [1,3] )
   id = selectid(:male => x -> x==true, group, aliveonly=false)
   @test all( id .== [1,2,3] )
   id = selectid(:male => x -> x==false, group)
   @test all( id .== [5,6] )
   id = selectid(:male => x -> x==false, group, aliveonly=false)
   @test all( id .== [5,6,7] )
   # multiple conditions
   id = selectid([:male,:tbv] => (x,y) -> x==true && y>=3.0, group)
   @test all( id .== [3] )
   id = selectid([:male,:tbv] => (x,y) -> x==true && y>=3.0, group, aliveonly=false)
   @test all( id .== [2,3] )
   # missing values
   id = selectid([:male,:gebv] => (x,y) -> x==true && !ismissing(y) && y>=2.0, group)
   @test all( id .== [1,3] )
   id = selectid([:male,:gebv] => (x,y) -> x==true && !ismissing(y) && y>=2.0, group, aliveonly=false)
   @test all( id .== [1,2,3] )
   # allowempty
   id = selectid([:male,:gebv] => (x,y) -> x==true && !ismissing(y) && y>=5.0, group, allowempty=true)
   @test all( id .== [] )
   @test_throws ErrorException selectid([:male,:gebv] => (x,y) -> x==true && !ismissing(y) && y>=5.0, group, allowempty=false)
   # sort
   id = selectid(:male => x -> x==true, group, sortby=:tbv)
   @test all( id .== [3,1] )
   id = selectid(:male => x -> x==true, group, sortby=:tbv, rev=false)
   @test all( id .== [1,3] )
   id = selectid(:male => x -> x==true, group, sortby=:tbv, aliveonly=false)
   @test all( id .== [2,3,1] )
   id = selectid(:male => x -> x==true, group, sortby=:gebv)
   @test all( id .== [3,1] )
   id = selectid([:male,:gebv] => (x,y) -> x==true && !ismissing(y), group, sortby=:gebv)
   @test all( id .== [3,1] )
   id = selectid([:male,:gebv] => (x,y) -> x==true && !ismissing(y) && y>=5.0, group, sortby=:gebv, allowempty=true)
   @test all( id .== [] )

   # removed dead individuals
   group = generate_group(pop, sires=pop_bulls[1:3], dams=pop_dams[1:3], aliveonly=true)
   # simple selection: default aliveonly=true
   id = selectid(:male => x -> x==true, group)
   @test all( id .== [1,3] )
   id = selectid(:male => x -> x==true, group, aliveonly=false)
   @test all( id .== [1,3] )
   id = selectid(:male => x -> x==false, group)
   @test all( id .== [5,6] )
   id = selectid(:male => x -> x==false, group, aliveonly=false)
   @test all( id .== [5,6] )
   # multiple conditions
   id = selectid([:male,:tbv] => (x,y) -> x==true && y>=3.0, group)
   @test all( id .== [3] )
   id = selectid([:male,:tbv] => (x,y) -> x==true && y>=3.0, group, aliveonly=false)
   @test all( id .== [3] )
   # missing values
   id = selectid([:male,:gebv] => (x,y) -> x==true && !ismissing(y) && y>=2.0, group)
   @test all( id .== [1,3] )
   id = selectid([:male,:gebv] => (x,y) -> x==true && !ismissing(y) && y>=2.0, group, aliveonly=false)
   @test all( id .== [1,3] )
   # allowempty
   id = selectid([:male,:gebv] => (x,y) -> x==true && !ismissing(y) && y>=5.0, group, allowempty=true)
   @test all( id .== [] )
   @test_throws ErrorException selectid([:male,:gebv] => (x,y) -> x==true && !ismissing(y) && y>=5.0, group, allowempty=false)
   # sort
   id = selectid(:male => x -> x==true, group, sortby=:tbv)
   @test all( id .== [3,1] )
   id = selectid(:male => x -> x==true, group, sortby=:tbv, rev=false)
   @test all( id .== [1,3] )
   id = selectid(:male => x -> x==true, group, sortby=:tbv, aliveonly=false)
   @test all( id .== [3,1] )
   id = selectid(:male => x -> x==true, group, sortby=:gebv)
   @test all( id .== [3,1] )
   id = selectid([:male,:gebv] => (x,y) -> x==true && !ismissing(y), group, sortby=:gebv)
   @test all( id .== [3,1] )
   id = selectid([:male,:gebv] => (x,y) -> x==true && !ismissing(y) && y>=5.0, group, sortby=:gebv, allowempty=true)
   @test all( id .== [] )
end


@testset "add_sires! and add_dams!" begin
   par = PTParameters(50, 100, 0.5, 0, 0, 1.0)
   hp = generate_population(par,nm=5,nf=5)
   pop = generate_population(par)
   hp_males = [3,2,4,5]
   hp_females = [8,7,10,9]
   pop_bulls = migrate_from_hp!(hp,pop,hp_males)
   pop_dams = migrate_from_hp!(hp,pop,hp_females)
   pop.df[2,:alive] = false
   pop.df[7,:alive] = false
   group = generate_group(pop, sires=pop_bulls[3:3], dams=pop_dams[2:2])

   # capacity full
   @test 0 == add_sires!(group,[1,2])
   # extended
   group.maxSire = 2
   @test 1 == add_sires!(group,[1,2])
   group.maxSire = 3
   @test 0 == add_sires!(group,[1,3,1,2,3])
   @test 0 == add_sires!(group,[5,6,7,8])
   @test all( sort(group.sires) .==  [1,3])
   @test vacancy_for_sires(group) == 1

   # capacity full
   @test 0 == add_dams!(group,[5,6,7])
   # extended
   group.maxDam = 2
   @test 1 == add_dams!(group,[5,6,7])
   group.maxDam = 3
   @test 0 == add_dams!(group,[5,6,7])
   group.maxDam = 4
   @test 0 == add_dams!(group,[5,6,5,7,6,7])
   @test 0 == add_dams!(group,[1,2,3,4])
   @test all( sort(group.dams) .==  [5,6])
   @test vacancy_for_dams(group) == 2
end

@testset "cull!" begin
   par = PTParameters(50, 100, 0.5, 0, 0, 1.0)
   hp = generate_population(par,nm=5,nf=5)
   pop = generate_population(par)
   hp_males = [3,2,4,5]
   hp_females = [8,7,10,9]
   pop_bulls = migrate_from_hp!(hp,pop,hp_males)
   pop_dams = migrate_from_hp!(hp,pop,hp_females)
   cull!(pop,[2,7])
   @test !pop.df[2,:alive] && !pop.df[7,:alive]
   group = generate_group(pop, sires=pop_bulls[3:3], dams=pop_dams[2:2])

   group.maxSire = 3
   @test 2 == add_sires!(group,[1,2,3,4])
   @test all( sort(group.sires) .==  [1,3,4])
   @test 1 == cull!(group,[1])
   @test all( sort(group.sires) .==  [3,4])
   @test vacancy_for_sires(group) == 1
   @test 1 == cull!(group,[2,3])
   @test all( sort(group.sires) .==  [4])
   @test vacancy_for_sires(group) == 2
   @test 1 == cull!(group,[1,2,3,4,3,2,1])
   @test vacancy_for_sires(group) == 3

   group.maxDam = 4
   @test 2 == add_dams!(group,[5,6,7,8])
   @test all( sort(group.dams) .==  [5,6,8])
   @test 0 == cull!(group,[1])
   @test 0 == cull!(group,[1,2])
   @test 2 == cull!(group,[1,2,5,6])
   @test all( sort(group.dams) .==  [8])
   @test vacancy_for_dams(group) == 3
end

@testset "pedigree list and inbreeding" begin
   sires = [0,0,1,1,3,1,5]
   dams = [0,0,0,2,2,4,6]
   pedlist = ProgenyTestingTools.get_pedigree_list(sires,dams)
   f = ProgenyTestingTools.kernel_meuwissen_and_luo!(Float64, pedlist)
   ref_f = [ 0.0, 0.0, 0.0, 0.0, 0.0, 0.25, 0.15625]
   @test f ≈ ref_f
   g = ProgenyTestingTools.get_inbreeding(sires,dams)
   @test g ≈ ref_f

   par = PTParameters(50, 100, 0.5, 0, 0, 1.0)
   hp = generate_population(par,nm=3,nf=4)
   hp.df[:,:male] = [true,false,true,false,true,false,true]
   hp.df[:,:sire] = sires
   hp.df[:,:dam] = dams
   update_inbreeding!(hp)
   @test hp.df[:,:inb] ≈ ref_f
end

@testset "mating between groups! (dairy_standard_ai)" begin
   par = PTParameters(50, 100, 0.5, 0, 0, 1.0)
   hp = generate_population(par,nm=10,nf=10)
   pop = generate_population(par)
   bulls1 = migrate_from_hp!(hp,pop,random_sampling(hp,4,male=true))
   bulls2 = migrate_from_hp!(hp,pop,random_sampling(hp,4,male=true))
   cows = migrate_from_hp!(hp,pop,random_sampling(hp,8,female=true))
   mgroup1 = generate_group(pop, sires=bulls1)
   mgroup2 = generate_group(pop, sires=bulls2)
   fgroup = generate_group(pop, dams=cows)
   for i=1:4
      mating!(mgroup1, fgroup, "fgroup", n=1, method="dairy_standard_ai", plan="once_per_female", calving=false)
   end
   mating!(mgroup2, fgroup, "fgroup", n="all", method="dairy_standard_ai", plan="once_per_female", calving=false)

   # find progeny from sires in mgroup
   # extract sire group 1 by x->x[1]==mgroup1.gid
   @test all(fgroup.generation[9:16] .== 1)
   prog_id = selectid([:siregroup] => x->x==mgroup1.groupid, fgroup)
   @test all( (17 .<= prog_id) .& (prog_id .<= 20) )
   prog_id = selectid([:siregroup] => x->x==mgroup2.groupid, fgroup)
   @test all( (21 .<= prog_id) .& (prog_id .<= 24) )

   @test length(selectid([:pregnant] => x->x==true, pop))==8
   calving!(fgroup)
   @test length(selectid([:pregnant] => x->x==true, pop))==0

   # once more
   mgroup3 = copy_group(mgroup1)
   mgroup4 = copy_group(mgroup2)
   fgroup.dams .= sort(shuffle(17:24)[1:8])
   for i=1:4
      mating!(mgroup3, fgroup, "fgroup", n=1, method="dairy_standard_ai", plan="once_per_female", calving=false)
   end
   mating!(mgroup4, fgroup, "fgroup", n="all", method="dairy_standard_ai", plan="once_per_female", calving=false)

   @test all(fgroup.generation[17:24] .== 2)
   prog_id = selectid([:siregroup] => x->x==mgroup3.groupid, fgroup)
   @test all( (25 .<= prog_id) .& (prog_id .<= 28) )
   prog_id = selectid([:siregroup] => x->x==mgroup4.groupid, fgroup)
   @test all( (29 .<= prog_id) .& (prog_id .<= 32) )

   @test length(selectid([:pregnant] => x->x==true, pop))==8
   calving!(fgroup)
   @test length(selectid([:pregnant] => x->x==true, pop))==0
end

@testset "assign_phenotype!" begin
   par = PTParameters(50, 100, 0.5, 0, 0, 1.0)
   hp = generate_population(par,nm=10,nf=10)
   pop = generate_population(par)
   bulls = migrate_from_hp!(hp,pop,random_sampling(hp,4,male=true))
   cows = migrate_from_hp!(hp,pop,random_sampling(hp,4,female=true))
   mgroup = generate_group(pop, sires=bulls)
   fgroup = generate_group(pop, dams=cows)
   mating!(mgroup, fgroup, "fgroup", n="all", method="dairy_standard_ai", plan="once_per_female", calving=false, pmp=0.0)
   calving!(fgroup)
   assign_phenotype!(fgroup, gen="all", repeated=false)
   @test sum( map(x->length(x.cg), pop.animal) ) == sum(.!pop.df[fgroup.id,:male] )

   assign_phenotype!(fgroup, gen="all", repeated=false)
   assign_phenotype!(fgroup, gen="all", repeated=false)
   @test maximum(pop.df[:,:nrec]) == 1

   bulls2 = migrate_from_hp!(hp,pop,random_sampling(hp,4,male=true))
   mgroup2 = generate_group(pop, sires=bulls)
   fgroup.dams = selectid([:male]=>x->x!=true, fgroup, idlist=fgroup.id[fgroup.generation.==1])
   mating!(mgroup2, fgroup, "fgroup", n="all", method="dairy_standard_ai", plan="once_per_female", calving=false)
   calving!(fgroup)
   assign_phenotype!(fgroup, gen=2, repeated=false)
   @test length(selectid([:siregroup,:male]=>(x,y)->x==mgroup2.groupid && y==false,fgroup)) == sum(.!pop.df[fgroup.id[fgroup.generation .== 2],:male] )
end
