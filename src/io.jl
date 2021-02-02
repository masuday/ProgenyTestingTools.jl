# functions for io

# not done yet
function show(io::IO, ::MIME"text/plain", group::PTGroup)
   print(io,"group_id = $(group.groupid)\n")
   print(io,"n = $(group.n)\n")
   print(io,"sires = $(group.sires)\n")
   print(io,"dams = $(group.dams)\n")
   if isempty(group.generation)
      print(io,"max generation = empty list")
   else
      print(io,"max generation = $(maximum(group.generation))")
   end
end

# write a data file
function write(io::IO, pop::PTPopulation; header::Bool=true, missing::String="NA", repeated::Bool=true, skipmissing::Bool=false)
   if header
      print(io,@sprintf("#%8s%8s%8s%7s%5s%3s%11s%11s%11s%3s%3s%8s%11s%11s%11s\n","animal","sire","dam","inb","yr","sx", "TBV","PedBV","QTLBV","n","mu","cg","pheno","pe","e"))
   end
   for i in 1:pop.maxAnimal
      # animal, sire, dam
      s = pop.df[i,:sire]
      d = pop.df[i,:dam]
      inb = pop.df[i,:inb]
      sex = ifelse(pop.df[i,:male],"1","2")
      year = pop.df[i,:year]
      pedigree = @sprintf("%8d%8d%8d%7.4f%5d%3s",i,s,d,inb,year,sex)
      bv = @sprintf("%11.3g%11.3g%11.3g",pop.df[i,:tbv],pop.df[i,:pbv],pop.df[i,:qbv])

      # phenotypes
      nrec = length(pop.animal[i].y)
      if skipmissing && nrec==0
         continue
      end
      if nrec==0
         y = @sprintf("%11s",missing)
         pe = 0.0
         e = 0.0
         cg = 0
         record = @sprintf("%3d%3d%8d%11s%11.3g%11.3g",nrec,1,cg,y,pe,e)
         print(io," " * pedigree * bv * record * "\n")
      else
         if !repeated
            nrec = 1
         end
         for rec in 1:nrec
            y = pop.animal[i].y[rec]
            pe = pop.animal[i].pe[rec]
            e = pop.animal[i].e[rec]
            cg = pop.animal[i].cg[rec]
            record = @sprintf("%3d%3d%8d%11.3g%11.3g%11.3g",rec,1,cg,y,pe,e)
            print(io," " * pedigree * bv * record * "\n")
         end
      end
   end
end

function write_pedigree(io::IO, pop::PTPopulation; header::Bool=true)
   if header
      print(io,@sprintf("#%8s%8s%8s%7s%5s%3s%7s%4s%11s%11s%11s\n","animal","sire","dam","inb","yr","sx","recdau","rec","tbv","tbvpoly","tbvqtl"))
   end
   for i in 1:pop.maxAnimal
      # animal, sire, dam
      s = pop.df[i,:sire]
      d = pop.df[i,:dam]
      inb = pop.df[i,:inb]
      sex = ifelse(pop.df[i,:male],"1","2")
      year = pop.df[i,:year]
      code = get_inbupg_code(s,d,pop.df.inb)
      pedigree = @sprintf("%8d%8d%8d%5d%7.4f%5d%3s%7d%4d%11.3g%11.3g%11.3g",i,s,d,code,inb,year,sex,pop.df[i,:nrecprog],pop.df[i,:nrec],pop.df[i,:tbv],pop.df[i,:pbv],pop.df[i,:qbv])
      print(io," " * pedigree * "\n")
   end
end

function get_inbupg_code(s::Int,d::Int,inb::Vector{Float64})
   if s>0
      ms = 0.0
      fs = inb[s]
   else
      ms = 1.0
      fs = 0.0
   end
   if d>0
      md = 0.0
      fd = inb[d]
   else
      md = 1.0
      fd = 0.0
   end
   return round( 4000/((1+ms)*(1-fs) + (1+md)*(1-fd)) )
end
