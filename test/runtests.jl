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

