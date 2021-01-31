# function for parameters

"""
    check_parameters(par::PTParameters)

Check errors in `par`.
Return `nothing` if there is no error.
"""
function check_parameters(par::PTParameters, warn::Bool=false)
   if par.var_p<0
      error("negative phenotypic variance")
   end
   if par.var_g_sim<0
      error("negative QMSim genetic variance")
   end
   if par.h2_poly<0 || par.h2_poly>1 || par.h2_qtl<0 || par.h2_qtl>1 || (par.h2_poly+par.h2_qtl)<0 || (par.h2_poly+par.h2_qtl)>1
      error("heritability out of range")
   end
   if par.rep>0 && (par.rep < par.h2_poly || par.rep < par.h2_qtl || par.rep < (par.h2_poly + par.h2_qtl))
      error("repeatability>0 but smaller than heritability")
   end
   if par.rep < 0 || par.rep > 1
      error("repeatability out of range")
   end

   if warn && get_var_pe(par)<=0.0
      @warn "PE variance is set to 0 (rep==0 or rep==h2). Be aware that the repeatability model would not be appropriate if you have repeated records."
   end
end

function get_var_poly(par::PTParameters)
   return par.h2_poly * par.var_p
end

function get_var_qtl(par::PTParameters)
   return par.h2_qtl * par.var_p
end

function get_var_g_sim(par::PTParameters)
   return par.var_g_sim
end

function get_var_pe(par::PTParameters)
   var_u = get_var_poly(par) + get_var_qtl(par)
   return max((par.rep * par.var_p)-var_u, 0.0)
end

function get_var_error(par::PTParameters)
   return par.var_p - (get_var_poly(par) + get_var_qtl(par) + get_var_pe(par))
end
