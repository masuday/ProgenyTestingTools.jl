# functions for pedigree

"""
    update_inbreeding!(pop)

Calculate inbreeding coeffcients for individuals in the population.
"""
function update_inbreeding!(sires::Vector{Int}, dams::Vector{Int}, inb::Vector{Float64}; first::Int=1)
   if first == 1
      f = get_inbreeding(sires, dams)
      if length(f) == length(inb)
         inb .= f
      else
         error("size mismatch: f and inb")
      end
   else
      get_inbreeding!(sires, dams, inb, first=first)
   end
   return
end

function update_inbreeding!(pop::PTPopulation; first::Int=1)
   update_inbreeding!(pop.df.sire, pop.df.dam, pop.df.inb, first=first)
end

"""
   f = get_inbreeding(sires::Vector{Int},dams::Vector{Int})
   get_inbreeding(sires::Vector{Int},dams::Vector{Int},inb::Vector{Float64})

Calculate inbreeding coeffcients out of pedigree arrays.
"""
function get_inbreeding(sires::Vector{Int}, dams::Vector{Int})
   pedlist = get_pedigree_list(sires,dams,ml=true)
   f = kernel_meuwissen_and_luo!(Float64, pedlist)
   return f
end

function get_inbreeding!(sires::Vector{Int}, dams::Vector{Int}, inb::Vector{Float64}; first::Int=1)
   pedlist = get_pedigree_list(sires,dams,ml=true)
   kernel_meuwissen_and_luo2!(pedlist, inb, first=first)
   return nothing
end

"""
    pedlist = get_pedigree_list(sires::Vector{Int},dams::Vector{Int}; ml=false)

Convert two pedigree vectors to a pedigree matrix.
If `ml=true`, sire ID and dam ID would be exchanged so that sire ID > dam ID for the Meuwissen and Luo algorithm.
"""
function get_pedigree_list(sires::Vector{Int}, dams::Vector{Int}; ml::Bool=false)
   if length(sires) != length(dams)
      throw(ArgumentError("unequal size of sire and dam vectors"))
   end
   n = length(sires)
   pedlist = zeros(Int,2,n)
   for i=1:n
      s = sires[i]
      d = dams[i]
      if ml==true
         pedlist[1,i] = max(s,d)
         pedlist[2,i] = min(s,d)
      else
         pedlist[1,i] = s
         pedlist[2,i] = d
      end
   end
   return pedlist
end

# kernel of Meuwissen and Luo (1992)
# the argument "ped" being rewritten
function kernel_meuwissen_and_luo!(Tv::DataType, ped::Matrix{Ti}) where Ti<:Integer
   n::Ti = size(ped,2)
   f = OffsetArray{Tv}(undef,0:n)
   i::Ti = 0
   j::Ti = 0
   k::Ti = 0
   s0::Ti = 0
   d0::Ti = 0
   ks::Ti = 0
   kd::Ti = 0
   fi::Tv = 0.0
   r::Tv = 0.0
   point = zeros(Ti,n)
   T = zeros(Tv,n)
   B = zeros(Tv,n)

   # start
   f[0] = -1
   for i=1:n
      s0 = ped[1,i]
      d0 = ped[2,i]
      ped[1,i] = max(s0,d0)
      ped[2,i] = min(s0,d0)
      B[i] = 0.5 - 0.25*(f[s0]+f[d0])

      if s0==0 || d0==0
         f[i] = 0.0
      else
         fi = -1.0
         T[i] = 1.0
         j = i
         while(j!=0)
            k = j
            r = 0.5 * T[k]
            ks = ped[1,k]
            kd = ped[2,k]
            if ks != 0
               while(point[k]>ks)
                  k = point[k]
               end
               T[ks] = T[ks] + r
               if ks != point[k]
                  point[ks] = point[k]
                  point[k] = ks
               end
               if kd != 0
                  while(point[k]>kd)
                     k = point[k]
                  end
                  T[kd] = T[kd] + r
                  if kd != point[k]
                     point[kd] = point[k]
                     point[k] = kd
                  end
               end
            end
            fi = fi + T[j]*T[j]*B[j]
            T[j] = 0.0
            k = j
            j = point[j]
            point[k] = 0
         end
         f[i] = fi
      end
   end
   return f[1:n]
end

# kernel of Meuwissen and Luo (1992)
# specialized for new animals in the same generation
# the argument "ped" not being rewritten
function kernel_meuwissen_and_luo2!(ped::Matrix{Ti}, inb::Vector{Tv}; first=1) where {Tv<:Real, Ti<:Integer}
   n::Ti = size(ped,2)
   if n<length(inb)
      throw(DimensionMismatch("inb shorter than pedigree size"))
   end
   f = OffsetArray{Tv}(undef,0:n)
   f[1:n] .= inb[1:n]
   point = zeros(Ti,n,Threads.nthreads())
   T = zeros(Tv,n,Threads.nthreads())
   B = zeros(Tv,n)

   # start
   f[0] = -1
   for i=1:first-1
      s0 = ped[1,i]
      d0 = ped[2,i]
      B[i] = 0.5 - 0.25*(f[s0]+f[d0])
   end
   Base.Threads.@threads for i=first:n
      local s0 = ped[1,i]
      local d0 = ped[2,i]
      B[i] = 0.5 - 0.25*(f[s0]+f[d0])

      if s0==0 || d0==0
         f[i] = 0.0
      else
         local fi = -1.0
         T[i,Threads.threadid()] = 1.0
         local j = i
         while(j!=0)
            local k = j
            local r = 0.5 * T[k,Threads.threadid()]
            local ks = ped[1,k]
            local kd = ped[2,k]
            if ks != 0
               while(point[k,Threads.threadid()]>ks)
                  k = point[k,Threads.threadid()]
               end
               T[ks,Threads.threadid()] = T[ks,Threads.threadid()] + r
               if ks != point[k,Threads.threadid()]
                  point[ks,Threads.threadid()] = point[k,Threads.threadid()]
                  point[k,Threads.threadid()] = ks
               end
               if kd != 0
                  while(point[k,Threads.threadid()]>kd)
                     k = point[k,Threads.threadid()]
                  end
                  T[kd,Threads.threadid()] = T[kd,Threads.threadid()] + r
                  if kd != point[k,Threads.threadid()]
                     point[kd,Threads.threadid()] = point[k,Threads.threadid()]
                     point[k,Threads.threadid()] = kd
                  end
               end
            end
            fi = fi + T[j,Threads.threadid()]*T[j,Threads.threadid()]*B[j]
            T[j,Threads.threadid()] = 0.0
            k = j
            j = point[j,Threads.threadid()]
            point[k,Threads.threadid()] = 0
         end
         f[i] = fi
      end
   end
   inb[first:n] .= f[first:n]
   return nothing
end
