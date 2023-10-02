import matlabparser as mpars
# --- Convert a matlab file
mend = '.m'
pend = '.py'
fdir = '/home/pete/PycharmProjects/pyroms_MI/Oforc_GLORYS_python/'
# Done
# # p2 = mpars.matlab2python(fdir + 'pw_calc_surfaceSWI' + mend, output=fdir + 'pw_calc_surfaceSWI' + pend)
# # p5 = mpars.matlab2python(fdir + 'pw_fnperiodic_distm' + mend, output=fdir + 'pw_fnperiodic_distm' + pend)
# # p8 = mpars.matlab2python(fdir + 'pw_Xplin' + mend, output=fdir + 'pw_Xplin' + pend)
# # p4 = mpars.matlab2python(fdir + 'pw_fn_biascorrect_projections_PW_Sep23' + mend,
# #                          output=fdir + 'pw_fn_biascorrect_projections_PW_Sep23' + pend)
# # p7 = mpars.matlab2python(fdir + 'pw_make_roho_biascorrected_forcings' + mend,
# #                          output=fdir + 'pw_make_roho_biascorrected_forcings' + pend)
# # p6 = mpars.matlab2python(fdir + 'pw_ks_regress' + mend, output=fdir + 'pw_ks_regress' + pend)

# Not needed vvvvv GO with more recent version of file sent by Phil in Iceland
# # p3 = mpars.matlab2python(fdir + 'pw_fn_biascorrect_projections' + mend,
# #                          output=fdir + 'pw_fn_biascorrect_projections' + pend)
#
# p8 = mpars.matlab2python(fdir + 'make_testdata_biascorrected_projections' + mend,
#                          output=fdir + 'make_testdata_biascorrected_projections' + pend)

p9 = mpars.matlab2python(fdir + 'nccopy_vars' + mend,
                         output=fdir + 'nccopy_vars' + pend)
p10 = mpars.matlab2python(fdir + 'squantile' + mend,
                          output=fdir + 'squantile' + pend)
p11 = mpars.matlab2python(fdir + 'optimize_theta' + mend,
                          output=fdir + 'optimize_theta' + pend)
p12 = mpars.matlab2python(fdir + 'fndatenum_to_decyr' + mend,
                          output=fdir + 'fndatenum_to_decyr' + pend)
p13 = mpars.matlab2python(fdir + 'fn_spherical_distm' + mend,
                          output=fdir + 'fn_spherical_distm' + pend)
p14 = mpars.matlab2python(fdir + 'J_fit_partially_linear_model' + mend,
                          output=fdir + 'J_fit_partially_linear_model' + pend)
