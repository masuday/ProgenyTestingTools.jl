using ProgenyTestingTools
using Test

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
end
