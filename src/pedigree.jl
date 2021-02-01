# functions for pedigree

"""
    update_inbreeding!(sires, dams, f)
    update_inbreeding!(pop)

Interenal routine.
"""
function update_inbreeding!(sires::Vector{Int}, dams::Vector{Int}, inb::Vector{Float64})
   f = get_inbreeding(sires, dams)
   if length(f) == length(inb)
      inb .= f
   else
      error("size mismatch: f and inb")
   end
   return
end

function update_inbreeding!(pop::PTPopulation)
   update_inbreeding!(pop.df.sire, pop.df.dam, pop.df.inb)
end

"""
   f = get_inbreeding(sires::Vector{Int},dams::Vector{Int})

Calculate inbreeding coeffcients out of pedigree arrays.
"""
function get_inbreeding(sires::Vector{Int}, dams::Vector{Int})
   pedlist = get_pedigree_list(sires,dams)
   f = kernel_meuwissen_and_luo!(Float64, pedlist)
   return f
end

"""
    pedlist = get_pedigree_list(sires::Vector{Int},dams::Vector{Int})

Convert two pedigree vectors to a pedigree matrix.
"""
function get_pedigree_list(sires::Vector{Int}, dams::Vector{Int})
   if length(sires) != length(dams)
      throw(ArgumentError("unequal size of sire and dam vectors"))
   end
   n = length(sires)
   pedlist = zeros(Int,2,n)
   for i=1:n
      s = sires[i]
      d = dams[i]
      pedlist[1,i] = s
      pedlist[2,i] = d
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
