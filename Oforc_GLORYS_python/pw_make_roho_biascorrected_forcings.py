
import numpy as np
import matplotlib.pyplot as plt

model = 'roho800'
oceanmod = 1
small_screen = 1
model_proj = 'noresm'
scenstr1 = 'ssp585'
use_daily = 1
interp_to_gregorian = 1
if use_daily == 0:
    interp_to_gregorian = 0

test_method = 0
seasonal = 1
use_month = 0
if use_daily == 1:
    use_month = 0

correct_iavar = 0
yearmin0 = 2000
yearmax0 = 2019

tclimc = np.arange(0.5,365.5+1)
method = 2

yearminref = 2010
yearmaxref = 2019

legacy = 0
use_subdelta = 1
remove_deltascale_var = 1
correct_substd = 0
correct_subquantiles = 1
if method == 2:
    correct_subquantiles = 0

match_subdelta_deltamean = 0
match_subdelta_hourly = 0
use_calc_surfaceSWI = 1
if oceanmod == 1:
    # dir1 = 'D:\USERS\Phil\ROHO800\ROHO800_input\projections\'
    dir1 = ''
else:
    # dir0 = 'C:\Data\ROHO800\ROHO800_input\projections\'
    dir0 = ''
    dir1 = dir0

##Grid point subselection indices (used only for atmospheric forcings)
#latsel = 11; lonsel = 30; #(58.5N,8.25E) -- actually not in ROHO800 domain
#latsel = 15; lonsel = 13; #(59.5N,4E)
#latsel = 19; lonsel = 25; #(60.5N,7E) -- innermost fjord area simulated in ROHO800 domain bottom right, Eidfjord
#latsel = 23; lonsel = 2; #(61.5N,1.25E) -- windy place!
latsel = np.arange(1,25+1)
lonsel = np.arange(1,37+1)

#latsel = 25; lonsel = 37;
#latsel = 16:20; lonsel = 19:27; #5x9 block containing target point (60.5N,7E) -> 54s with match_subdelta_deltamean=1
#latsel = 19:21; lonsel = 25:30; #3x6 block containing target point (60.5N,7E) -> 44s with match_subdelta_deltamean=1 --- needed for decent rain projections

nlatsel = len(latsel)
nlonsel = len(lonsel)
ns = nlatsel * nlonsel
##Atmospheric forcing variables
#varstr = 'Tair';
#varstr = 'Pair';
#varstr = 'Qair';
#varstr = 'cloud';
#varstr = 'Uwind';
#varstr = 'Vwind';
#varstr = 'swrad';
#varstr = 'lwrad_down';
#varstr = 'rain';

##Boundary condition variables
#varstr = 'temp';
#varstr = 'salt';
#varstr = 'zeta';
#varstr = 'ubar';
#varstr = 'vbar';
#varstr = 'u';
#varstr = 'v';

#varstr = 'O3_c';
#varstr = 'O3_TA';
#varstr = 'O2_o';
#varstr = 'N1_p';
#varstr = 'N3_n';
varstr = 'N5_s'
#bcstrv = '_west';
#bcstrv = '_north';
#bcstrv = '_east';
#bcstrv = '_south';
bcstrv = np.array(['_west','_north','_east','_south'])
#bcstrv = {'_north','_east','_south'};
if not iscell(bcstrv) :
    bcstrv = np.array([bcstrv])

do_ocean_bcs = 0
# dir0 = 'D:\USERS\Phil\ROHO800\ROHO800_input\ERA5\concatenated\'
dir0 = ''
source0str = 'ERA5 hindcast'
if np.any(str(varstr) == str(np.array(['temp','salt','zeta','ubar','vbar','u','v','O3_c','O3_TA','O2_o','N1_p','N3_n','N5_s']))):
    do_ocean_bcs = 1
    # dir0 = 'D:\USERS\Phil\ROHO800\ROHO800_input\old_vertgrid_pre_v2e\'
    dir0 = ''
    if np.any(str(varstr) == str(np.array(['temp','salt','zeta','ubar','vbar','u','v']))):
        source0str = 'GLORYS12 reanalysis'
    else:
        source0str = 'NORESM reanalysis v1'

if (do_ocean_bcs == 1 and test_method == 1):
    yearminref = 2014
    yearmaxref = 2019

if do_ocean_bcs == 1:
    nbcstr = len(bcstrv)
else:
    nbcstr = 1

#NOTE: RCM output from EURO-CORDEX (12x12km resolution, CMIP5) has already been bias-corrected
#      to a 1x1km observational dataset, using empirical quantile mapping with trends subtracted
#      and added back on (Wong et al., 2016; Hanssen-Bauer et al., 2015; 2017 (English)).
#      These data are available for temperature and precipitation from https://nedlasting.nve.no/klimadata/kss.
#      However, it is acknowledged that there may be issues due to lack of multivariate treatment
#      (e.g. warm spells should correlate with dry spells) and lack of spatial treatment (below the 12km
#      resolution of the RCM). Also, we need to have it bias-corrected to ERA5 for comparability with hindcast.

#NOTE: As pointed out by Maraun (2016), most quantile mapping approaches have been direct (seeking to
#      to correct the ESM/RCM output), but there are applications that treat it as a delta-change
#      (e.g. Willems and Vrac, 2011). This latter is equivalent to the 'repeated hindcast correction'
#      method in fn_biascorrect_projections.m (method=2), and may be the most conservative approach
#      to bias correction (because only the quantile-trend information from the ESM/RCM is used).
#      A possible drawback here is that potential improvements in trends after bias correction may
#      be lost (see Comparison_trends_NORESM_vs_biascorrected_QM_v1.pptx) --- although Maraun (2016) states
#      that "Current bias correction methods cannot plausibly correct climate change trends".

for ibcstr in np.arange(1,nbcstr+1).reshape(-1):
    if do_ocean_bcs == 1:
        bcstr = bcstrv[ibcstr]
    XpcL = []
    XpcH = []
    XpcfL = []
    XpcfH = []
    minX = 0
    if (str(varstr) == str('Tair') == 1):
        varstrc = 'Tair'
        varstrt = 'Tair_time'
        unitstr = ' [\circC]'
        varstr1 = 'TREFHT'
        #tol_dtclim = 5; fractional = 0; fractional_subdelta = 0;
        tol_dtclim = 50
        fractional = 0
        fractional_subdelta = 0
        if method == 2:
            tol_dtclim = 50
            fractional = 0
            fractional_subdelta = 0
            #tol_dtclim=5 with method=2 appears to be unrobust (for nq=100) and produces dubious extreme cold spells at (60.5N,7E)
#v2: tol_dtclim=20 days gives best monthly quantile maxMAE (=0.988 degC) vs. validation data (2000-2009), with match_spatially=0
#v2: tol_dtclim=20 days gives best monthly quantile maxMAE (=1.094 degC) vs. validation data (2000-2009), with match_spatially=1
#v3: tol_dtclim=50 days gives best monthly quantile meanRMSE,meanMAE (=0.440,0.339 degC) vs. validation data (2000-2009), with match_spatially=0
        # #     opt_subq = struct('bymonth',0,'model','tanh','nparb',3,'nonnegb',0,'pL',0,'pH',1,'nr',1,...
# #         'resmodel','ksmooth','bdw',5,'ndclosest',100,'kfunc',1,'xp',(-40:5:40)','doplot',1); #Old model for remove_deltascale_var=0
# #The quadratic-ksmooth approach seems to get a similar result (looking at 60.5N,7E) to the ksmooth-ksmooth model,
# #but the latter is probably more robust in extrapolation in general.
#     opt_subq = struct('bymonth',1,'model','polynomial','nparb',3,'nonnegb',0,...
#         'resmodel','ksmooth','bdwres',5,'ndclosestres',200,'kfuncres',1,'lr_orderres',0,'xpres',(-40:5:40)','doplot',1);
        opt_subq = struct('bymonth',1,'model','ksmooth','bdw',5,'ndclosest',200,'kfunc',1,'lr_order',1,'xp',np.transpose((np.arange(- 40,40+2.5,2.5))),'resmodel','ksmooth','bdwres',5,'ndclosestres',200,'kfuncres',1,'lr_orderres',0,'xpres',np.transpose((np.arange(- 40,40+2.5,2.5))),'doplot',1)
    else:
        if (str(varstr) == str('Pair') == 1):
            varstrc = 'Pair'
            varstrt = 'pair_time'
            unitstr = ' [Pa]'
            varstr1 = 'PS'
            tol_dtclim = 50
            fractional = 0
            fractional_subdelta = 0
            if (method == 2):
                tol_dtclim = 50
                fractional = 0
                fractional_subdelta = 0
                #v2: tol_dtclim=20 days gives best monthly quantile maxMAE (=335 Pa) vs. validation data (2000-2009), with match_spatially=0
#v2: tol_dtclim=20 days gives best monthly quantile maxMAE (=352 Pa) vs. validation data (2000-2009), with match_spatially=1
#v3: tol_dtclim=50 days gives best monthly quantile meanRMSE,meanMAE (=207,170 Pa) vs. validation data (2000-2009), with match_spatially=0
            opt_subq = struct('bymonth',1,'model','polynomial','logX',0,'logY',0,'nparb',2,'resmodel',[],'doplot',1)
        else:
            if (str(varstr) == str('Qair') == 1):
                varstrc = 'Qair'
                varstrt = 'qair_time'
                unitstr = ' [kg/kg]'
                varstr1 = 'QREFHT'
                tol_dtclim = 5
                fractional = 0
                fractional_subdelta = 0
                if (method == 2):
                    tol_dtclim = 50
                    fractional = 0
                    fractional_subdelta = 0
                    #v2: tol_dtclim=100 days gives best monthly quantile maxMAE (=0.299e-3 kg/kg) vs. validation data (2000-2009), with match_spatially=0
#v2: tol_dtclim=100 days gives best monthly quantile maxMAE (=0.302e-3 kg/kg) vs. validation data (2000-2009), with match_spatially=1
#v3: tol_dtclim=50 days gives best monthly quantile meanRMSE,meanMAE (=0.153e-3,0.131e-3 kg/kg) vs. validation data (2000-2009), with match_spatially=0
                opt_subq = struct('bymonth',1,'model','polynomial','logX',0,'logY',0,'nparb',2,'resmodel',[],'doplot',1)
                XpcL = 0
                XpcfL = 0
            else:
                if (str(varstr) == str('cloud') == 1):
                    varstrc = 'cloud'
                    varstrt = 'cloud_time'
                    unitstr = ''
                    varstr1 = 'CLDTOT'
                    tol_dtclim = 5
                    fractional = 0
                    fractional_subdelta = 0
                    XpcL = 0
                    XpcfL = 0
                    XpcH = 1
                    XpcfH = 1
                else:
                    if (str(varstr) == str('Uwind') == 1):
                        varstrc = 'Uwind'
                        varstrt = 'wind_time'
                        unitstr = ' [m/s]'
                        varstr1 = 'U'
                        tol_dtclim = 5
                        fractional = 0
                        fractional_subdelta = 0
                        if (method == 2):
                            seasonal = 0
                            #v2: tol_dtclim=100 days gives best monthly quantile maxMAE (=1.50 m/s) vs. validation data (2000-2009), with match_spatially=0
#v3: non-seasonal gives best monthly quantile meanRMSE,meanMAE (=0.68,0.62 m/s) vs. validation data (2000-2009), with match_spatially=0
                        # #     opt_subq = struct('bymonth',1,'model','tanh','logX',0,'logY',0,'nparb',3,'nonnegb',0,'pL',0,'pH',1,'nr',1,...
# #         'resmodel','ksmooth','bdw',2,'ndclosest',50,'kfunc',1,'xp',(-30:30)','doplot',1); #Old model for remove_deltascale_var=0
                        #     match_subdelta_deltamean = 1; match_subdelta_hourly = 1; correct_subquantiles = 0;
#     #This achieves marginally-best hourly quantile skill, but poorer daily-average skill than the ksmooth-ksmooth approach.
#     #Also, we still end up with a mean abs diff that is ~3x higher between 2300 and 0000 (cf. ~4x higher with match_subdelta_hourly=0).
                        #     opt_subq = struct('bymonth',1,'model','ksmooth','bdw',3,'ndclosest',200,'kfunc',1,'lr_order',1,'xp',(-30:1:30)',...
#         'resmodel','ksmooth','bdwres',3,'ndclosestres',200,'kfuncres',1,'lr_orderres',0,'xpres',(-30:1:30)','doplot',1); #Current BEST
#    opt_subq.lr_order = 0; #This doesn't improve skill
#     opt_subq.lr_orderres = 1; #This didn't seem to help
#     opt_subq.bymonth = 0; #This doesn't improve skill
                        #     correct_substd = 1; correct_subquantiles = 0;
#     opt_subq = struct('bymonth',1,'model','ksmooth','bdw',10,'ndclosest',50,'kfunc',1,'lr_order',0,'xp',(-30:1:30)','doplot',1); #This gives a significantly poorer fit
                        #     opt_subq = struct('bymonth',1,'model','ksmooth','bdw',10,'ndclosest',50,'kfunc',1,'lr_order',0,'xp',(-30:1:30)',...
#         'resmodel','ksmooth','bdwres',10,'ndclosestres',50,'kfuncres',1,'lr_orderres',0,'xpres',(-30:1:30)','doplot',1); #This gives a poorer fit
#     opt_subq.bymonth = 0; #This doesn't help
                        #     opt_subq = struct('bymonth',1,'model','polynomial','nparb',3,'resmodel','polynomial','logYres',1,'nparbres',3,'doplot',1); #This is not robust -- leads to too-strong winds
#     opt_subq = struct('bymonth',1,'model','polynomial','nparb',3,'resmodel','ksmooth','bdwres',3,'ndclosestres',200,'kfuncres',1,'lr_orderres',0,'xpres',(-30:1:30)','doplot',1); #This is almost as good as the ksmooth-ksmooth solution
                        opt_subq = struct('bymonth',1,'model','polynomial','nparb',3,'resmodel','gaussian','presL',0,'presH',1,'doplot',1)
                    else:
                        if (str(varstr) == str('Vwind') == 1):
                            varstrc = 'Vwind'
                            varstrt = 'wind_time'
                            unitstr = ' [m/s]'
                            varstr1 = 'V'
                            tol_dtclim = 5
                            fractional = 0
                            fractional_subdelta = 0
                            if (method == 2):
                                tol_dtclim = 100
                                seasonal = 0
                                #v2: For Uwind, tol_dtclim=100 days gives best monthly quantile maxMAE (=1.50 ms-1) vs. validation data (2000-2009), with match_spatially=0
#v3: non-seasonal gives monthly quantile meanRMSE,meanMAE (=0.67,0.62 m/s) vs. validation data (2000-2009), with match_spatially=0
                            #     opt_subq = struct('bymonth',1,'model','ksmooth','bdw',3,'ndclosest',200,'kfunc',1,'lr_order',1,'xp',(-30:1:30)',...
#         'resmodel','ksmooth','bdwres',3,'ndclosestres',200,'kfuncres',1,'lr_orderres',0,'xpres',(-30:1:30)','doplot',1);
                            #    opt_subq.lr_order = 0; #This doesn't improve skill
#     opt_subq.lr_orderres = 1; #This gave unstable results
                            #     opt_subq = struct('bymonth',1,'model','polynomial','nparb',3,'resmodel','polynomial','logYres',1,'nparbres',3,'doplot',1); #This is not robust -- leads to too-strong winds
#     opt_subq = struct('bymonth',1,'model','polynomial','nparb',3,'resmodel','ksmooth','bdwres',3,'ndclosestres',200,'kfuncres',1,'lr_orderres',0,'xpres',(-30:1:30)','doplot',1); #This is almost as good as the ksmooth-ksmooth solution
#     opt_subq = struct('bymonth',1,'model','polynomial','nparb',3,'resmodel','gaussian','presL',0,'presH',1,'doplot',1); #This actually achieves lower MAEs than the ksmooth solution, but maximum excursions are more severe.
                            opt_subq = struct('bymonth',1,'model','ksmooth','bdw',3,'ndclosest',200,'kfunc',1,'lr_order',1,'xp',np.transpose((np.arange(- 30,30+1,1))),'resmodel','gaussian','presL',0,'presH',1,'doplot',0)
                        else:
                            if (str(varstr) == str('swrad') == 1):
                                varstrc = 'swrad'
                                varstrt = 'swrad_time'
                                unitstr = ' [W/m2]'
                                varstr1 = 'FSDS'
                                tol_dtclim = 5
                                fractional = 1
                                fractional_subdelta = 1
                                remove_deltascale_var = 0
                                match_subdelta_deltamean = 1
                                correct_subquantiles = 0
                                #For swrad it seems better not to remove the daily (deltascale) variation when
#defining the (fractional) subdelta variability, since swrad drops to zero at the
#end of each day.
                                match_subdelta_hourly = 0
                                if (method == 2):
                                    tol_dtclim = 100
                                    fractional = 0
                                    fractional_subdelta = 1
                                    #fractional=0 gave better validation (2000-2009) skill over all months, at (60.5N,7E) with tol_dtclim = 5
#v2: tol_dtclim=100 days gives best monthly quantile maxMAE (=6.67 W/m2) vs. validation data (2000-2009), with match_spatially=0
#v2: tol_dtclim=100 days gives best monthly quantile maxMAE (=7.09 W/m2) vs. validation data (2000-2009), with match_spatially=1
#v3: tol_dtclim=100 days gives ~best monthly quantile meanRMSE,meanMAE (=4.81,3.18 W/m2) vs. validation data (2000-2009), with match_spatially=0
                                    match_subdelta_deltamean = 0
                                #opt_subq = struct('bymonth',1,'model','powerlaw','nparb',2,'pL',0,'pH',2,'nr',1,'resmodel',[],'doplot',0); #Model for remove_deltascale_var=0
#     opt_subq = struct('bymonth',1,'model','ksmooth','bdw',50,'ndclosest',200,'kfunc',1,'lr_order',1,'xp',(0:20:400)',...
#         'resmodel','ksmooth','bdwres',50,'ndclosestres',200,'kfuncres',1,'lr_orderres',0,'xpres',(0:20:400)','doplot',1);
                                XpcL = 0
                                XpcfL = 0
                                minXsub = 0
                                frcrit_capSWI = 0.0
                            else:
                                if str(varstr) == str('lwrad_down') == 1:
                                    varstrc = 'lwrad_down'
                                    varstrt = 'lwrad_time'
                                    unitstr = ' [W/m2]'
                                    varstr1 = 'FLDS'
                                    tol_dtclim = 5
                                    fractional = 1
                                    fractional_subdelta = 0
                                    if (method == 2):
                                        tol_dtclim = 100
                                        fractional = 0
                                        fractional_subdelta = 0
                                        #v2: tol_dtclim=20 days gives best monthly quantile maxMAE (=5.72 W/m2) vs. validation data (2000-2009), with match_spatially=0.
#v2: tol_dtclim=20 days gives best monthly quantile maxMAE (=6.04 W/m2) vs. validation data (2000-2009), with match_spatially=1.
#v3: tol_dtclim=50 days gives best monthly quantile meanRMSE,meanMAE (=2.96,2.69 W/m2) vs. validation data (2000-2009), with match_spatially=0.
                                    opt_subq = struct('bymonth',1,'model','powerlaw','nparb',2,'pL',0,'pH',2,'nr',1,'resmodel',[],'doplot',1)
                                    XpcL = 0
                                    XpcfL = 0
                                else:
                                    if str(varstr) == str('rain') == 1:
                                        varstrc = 'rain'
                                        varstrt = 'rain_time'
                                        unitstr = ' [kg/m2/s]'
                                        varstr1 = 'PRECT'
                                        #tol_dtclim = 5; fractional = 1; fractional_subdelta = 0;
#tol_dtclim = 50; fractional = 0; fractional_subdelta = 1; #This does not work at all (is unstable) with method=1
                                        tol_dtclim = 50
                                        fractional = 1
                                        fractional_subdelta = 0
                                        if method == 2:
                                            tol_dtclim = 100
                                            fractional = 0
                                            fractional_subdelta = 1
                                            #v2: tol_dtclim=20 days gives best monthly quantile maxMAE (=0.0126 g/m2/s) vs. validation data (2000-2009), with match_spatially=0.
#v3: tol_dtclim=50 days gives best monthly quantile meanRMSE,meanMAE (=0.0107,0.0049 g/m2/s) vs. validation data (2000-2009), with match_spatially=0.
                                        # # OLD models for remove_deltascale_var=0
#     #opt_subq = struct('bymonth',1,'model','polynomial','nparb',2,'logX',1,'logY',1,'resmodel',[],'doplot',1);  #This gives poor results
#     #opt_subq = struct('bymonth',1,'model','tanh','nparb',3,'nonnegb',0,'pL',0,'pH',1e4,'nr',1,... #This is ok, but fit not optimal for lower quantiles
                                        #     opt_subq = struct('bymonth',1,'model','ksmooth','bdw',5e-4,'ndclosest',200,'kfunc',1,'lr_order',1,'xp',1e-3*(0:0.1:2)',...
#         'resmodel','ksmooth','bdwres',5e-4,'ndclosestres',200,'kfuncres',1,'lr_orderres',0,'xpres',1e-3*(0:0.1:2)','doplot',0);
#     #This gave generally too high rainfall over the ROHO800 region, and warnings about matrices close to singular
                                        opt_subq = struct('bymonth',1,'model','powerlaw','nparb',2,'nonnegb',1,'pL',0,'pH',2,'nr',1,'resmodel','ksmooth','bdwres',0.0005,'ndclosestres',200,'kfuncres',1,'lr_orderres',0,'xpres',0.001 * np.transpose((np.arange(0,2+0.1,0.1))),'doplot',0)
                                        XpcL = 0
                                        XpcfL = 0
                                        minXsub = - 0.0001
                                    else:
                                        if str(varstr) == str('temp') == 1:
                                            varstrc0 = 'temp'
                                            varstrt = 'ocean_time'
                                            unitstr = ' [\circC]'
                                            varstr1 = 'thetao'
                                            varstrc = np.array([varstrc0, bcstr])
                                            use_daily = 0
                                            use_month = 1
                                            use_subdelta = 0
                                            tol_dtclim = 1
                                            fractional = 0
                                            fractional_subdelta = 0
                                            if method == 2:
                                                tol_dtclim = 1
                                        else:
                                            if str(varstr) == str('salt') == 1:
                                                varstrc0 = 'salt'
                                                varstrt = 'ocean_time'
                                                unitstr = ' [psu]'
                                                varstr1 = 'so'
                                                varstrc = np.array([varstrc0,bcstr])
                                                use_daily = 0
                                                use_month = 1
                                                use_subdelta = 0
                                                tol_dtclim = 1
                                                fractional = 0
                                                fractional_subdelta = 0
                                                if (method == 2):
                                                    tol_dtclim = 1
                                            else:
                                                if np.any(str(varstr) == str(np.array(['zeta', 'ubar', 'vbar', 'u', 'v']))):
                                                    varstrc0 = varstr
                                                    varstrt = 'ocean_time'
                                                    varstr1 = []
                                                    if str(varstr) == str('zeta'):
                                                        unitstr = ' [m]'
                                                    else:
                                                        unitstr = ' [m/s]'
                                                    varstrc = np.array([varstrc0,bcstr])
                                                    use_daily = 0
                                                    use_month = 1
                                                    use_subdelta = 0
                                                    tol_dtclim = 1
                                                    fractional = 0
                                                    fractional_subdelta = 0
                                                    if method == 2:
                                                        tol_dtclim = 1
                                                else:
                                                    if str(varstr) == str('O3_c') == 1:
                                                        varstrc0 = 'O3_c'
                                                        varstrt = 'ocean_time'
                                                        unitstr = ' [mmolC/m3]'
                                                        varstr1 = 'dissic'
                                                        varstrc = np.array([varstrc0,bcstr])
                                                        use_daily = 0
                                                        use_month = 1
                                                        use_subdelta = 0
                                                        tol_dtclim = 1
                                                        fractional = 0
                                                        fractional_subdelta = 0
                                                        if (method == 2):
                                                            tol_dtclim = 1
                                                    else:
                                                        if (str(varstr) == str('O3_TA') == 1):
                                                            varstrc0 = 'O3_TA'
                                                            varstrt = 'ocean_time'
                                                            unitstr = ' [mmol/m3]'
                                                            varstr1 = 'talk'
                                                            varstrc = np.array([varstrc0,bcstr])
                                                            use_daily = 0
                                                            use_month = 1
                                                            use_subdelta = 0
                                                            tol_dtclim = 1
                                                            fractional = 0
                                                            fractional_subdelta = 0
                                                            if (method == 2):
                                                                tol_dtclim = 1
                                                        else:
                                                            if (str(varstr) == str('O2_o') == 1):
                                                                varstrc0 = 'O2_o'
                                                                varstrt = 'ocean_time'
                                                                unitstr = ' [mmolO2/m3]'
                                                                varstr1 = 'o2'
                                                                varstrc = np.array([varstrc0,bcstr])
                                                                use_daily = 0
                                                                use_month = 1
                                                                use_subdelta = 0
                                                                tol_dtclim = 1
                                                                fractional = 0
                                                                fractional_subdelta = 0
                                                                if (method == 2):
                                                                    tol_dtclim = 1
                                                            else:
                                                                if (str(varstr) == str('N1_p') == 1):
                                                                    varstrc0 = 'N1_p'
                                                                    varstrt = 'ocean_time'
                                                                    unitstr = ' [mmolP/m3]'
                                                                    varstr1 = 'po4'
                                                                    varstrc = np.array([varstrc0,bcstr])
                                                                    use_daily = 0
                                                                    use_month = 1
                                                                    use_subdelta = 0
                                                                    tol_dtclim = 1
                                                                    fractional = 1
                                                                    fractional_subdelta = 0
                                                                    if (method == 2):
                                                                        tol_dtclim = 1
                                                                else:
                                                                    if (str(varstr) == str('N3_n') == 1):
                                                                        varstrc0 = 'N3_n'
                                                                        varstrt = 'ocean_time'
                                                                        unitstr = ' [mmolN/m3]'
                                                                        varstr1 = 'no3'
                                                                        varstrc = np.array([varstrc0,bcstr])
                                                                        use_daily = 0
                                                                        use_month = 1
                                                                        use_subdelta = 0
                                                                        tol_dtclim = 1
                                                                        fractional = 1
                                                                        fractional_subdelta = 0
                                                                        if (method == 2):
                                                                            tol_dtclim = 1
                                                                    else:
                                                                        if (str(varstr) == str('N5_s') == 1):
                                                                            varstrc0 = 'N5_s'
                                                                            varstrt = 'ocean_time'
                                                                            unitstr = ' [mmolSi/m3]'
                                                                            varstr1 = 'si'
                                                                            varstrc = np.array([varstrc0,bcstr])
                                                                            use_daily = 0
                                                                            use_month = 1
                                                                            use_subdelta = 0
                                                                            tol_dtclim = 1
                                                                            fractional = 1
                                                                            fractional_subdelta = 0
                                                                            if (method == 2):
                                                                                tol_dtclim = 1
    tstrv = np.array(['January','February','March','April','May','June','July','August','September','October','November','December'])
    if (use_daily == 1):
        scstr2 = 'daily'
        if seasonal == 1:
            tolstr = np.array([', tol=',str(tol_dtclim),' day(s)'])
        else:
            tolstr = ', non-seasonal'
    else:
        scstr2 = 'monthly'
        if seasonal == 1:
            tolstr = np.array([', tol=',str(tol_dtclim),' month(s)'])
        else:
            tolstr = ', non-seasonal'
    if (use_daily == 1):
        fname11 = 'NorESM2-MM.atm.day_HIST_SSP585_2000-2100_ROHO800.nc'
        fname1 = np.array([dir1,fname11])
        if np.any(str(varstr1) == str(np.array(['U','V']))):
            X1 = np.squeeze(ncread(fname1,varstr1,np.array([lonsel(1),latsel(1),1,1]),np.array([nlonsel,nlatsel,1,Inf])))
        else:
            X1 = np.squeeze(ncread(fname1,varstr1,np.array([lonsel(1),latsel(1),1]),np.array([nlonsel,nlatsel,Inf])))
        if str(varstr) == str('Tair'):
            X1 = X1 - 273.15
        if str(varstr) == str('rain'):
            X1 = X1 * 1000
        time1 = ncread(fname1,'time')
        lat11 = np.squeeze(ncread(fname1,'lat',latsel(1),nlatsel))
        lon11 = np.squeeze(ncread(fname1,'lon',lonsel(1),nlonsel))
        td11 = np.transpose((np.arange(datenum(2001,1,1),datenum(2001,12,31)+1)))
        __,month11,day11 = datevec(td11)
        year1 = kron(np.transpose((np.arange(2000,2099+1))),np.ones((365,1)))
        month1 = kron(np.ones((100,1)),month11)
        day1 = kron(np.ones((100,1)),day11)
        year1 = year1(np.arange(2,end()+1))
        month1 = month1(np.arange(2,end()+1))
        day1 = day1(np.arange(2,end()+1))
        td1 = datenum(year1,month1,day1)
        #Note that they are skipping instances of February 29th, hence max(diff(td1))=2.
        print(np.array(['Loaded ',varstr1,' from: ',fname1]))
        fname01 = np.array(['roho800_',varstr,'_',str(yearmin0),'_',str(yearmax0),'_daymean.nc'])
        fname0 = np.array([dir0,'daymean',fname01])
        X0 = np.squeeze(ncread(fname0,varstrc,np.array([lonsel(1),latsel(1),1]),np.array([nlonsel,nlatsel,Inf])))
        time0 = ncread(fname0,varstrt)
        td0 = datenum(1948,1,1) + time0
        year0,month0,day0 = datevec(td0)
        yrday0 = td0 - datenum(year0,1,1)
        tyr0 = fndatenum_to_decyr(td0)
        print(np.array(['Loaded ',varstrc,' from: ',fname0]))
        dtdmean = unique(td0 - int(np.floor(td0)))
        #For cdo, timestamps of daily means use the mean of the timestamps for all hourly data used,
#which is the mean of 0000,0100,...,2300 in the case of non-flux data and 0030,0130,...,2330
#in the case of flux data. We correct the projection time stamps for consistency with this:
        td1 = td1 + dtdmean
        yrday1 = td1 - datenum(year1,1,1)
    else:
        if (do_ocean_bcs == 1):
            fname11 = 'NorESM2-MM.ocn.mon_HIST_SSP585_1990-2100_ROHO800.nc'
            fname1 = np.array([dir1,fname11])
            if np.any(str(varstr) == str(np.array(['temp','salt','zeta','ubar','vbar','u','v']))):
                fname01 = 'roho800_bry_GLORYS_20070115_to_20200115_N25_fix.nc'
            else:
                fname01 = 'roho800_v2ax_bry_NORESMOCv1p2reanal_20070115_to_20211215_closest.nc'
            fname0 = np.array([dir0,fname01])
            if np.any(str(varstr) == str(np.array(['zeta','ubar','vbar','u','v']))):
                if method != 2:
                    stop
                #For the circulation variables we assume no climatic change, but in order to preserve the
#eddy variability we do not force with climatology. Instead we force with repeated (and
#uncorrected) blocks of hindcast variability. This is achieved by applying the RHC method
#(method=2) with constant (unity-valued) projection data.
                X1 = 1
            else:
                X1 = ncread(fname1,varstr1)
            if np.any(str(varstr) == str(np.array(['O3_c','O3_TA','O2_o','N1_p','N3_n','N5_s']))):
                X1 = X1 * 1000
            lat1 = ncread(fname1,'latitude')
            lon1 = ncread(fname1,'longitude')
            depth_bnds = ncread(fname1,'lev_bnds')
            z11 = np.transpose(mean(depth_bnds))
            nz11 = len(z11)
            X0 = ncread(fname0,varstrc)
            fnameg = np.array([dir1,'../Grid/ROHO800_grid_fix5_unsmooth1.nc'])
            Lats0 = ncread(fnameg,'lat_rho')
            Lons0 = ncread(fnameg,'lon_rho')
            H0 = ncread(fnameg,'h')
            M0 = ncread(fnameg,'mask_rho')
            year11 = np.transpose((np.arange(1990,2100+1)))
            td00 = datenum(1948,1,1)
            timeunit = 1 / 86400
        else:
            fname1 = np.array([dir1,'NorESM2-MM.atm.day_HIST_SSP585_2000-2100_ROHO800_monmean.nc'])
            fname0 = np.array([dir0,'roho800_',varstr,'_',str(yearmin0),'_',str(yearmax0),'_monmean.nc'])
            if np.any(str(varstr1) == str(np.array(['U','V']))):
                X1 = np.squeeze(ncread(fname1,varstr1,np.array([lonsel(1),latsel(1),1,1]),np.array([nlonsel,nlatsel,1,Inf])))
            else:
                X1 = np.squeeze(ncread(fname1,varstr1,np.array([lonsel(1),latsel(1),1]),np.array([nlonsel,nlatsel,Inf])))
            if str(varstr) == str('Tair'):
                X1 = X1 - 273.15
            X0 = np.squeeze(ncread(fname0,varstrc,np.array([lonsel(1),latsel(1),1]),np.array([nlonsel,nlatsel,Inf])))
            lat11 = np.squeeze(ncread(fname1,'lat',latsel(1),nlatsel))
            lon11 = np.squeeze(ncread(fname1,'lon',lonsel(1),nlonsel))
            year11 = np.transpose((np.arange(2000,2099+1)))
            td00 = datenum(1948,1,1)
            timeunit = 1
        time1 = ncread(fname1,'time')
        year1 = kron(year11,np.ones((12,1)))
        month1 = kron(np.ones((len(year11),1)),np.transpose((np.arange(1,12+1))))
        day1 = 15 * np.ones((len(year11) * 12,1))
        td1 = datenum(year1,month1,day1)
        print(np.array(['Loaded ',varstr1,' from: ',fname1]))
        time0 = ncread(fname0,varstrt)
        td0 = td00 + time0 * timeunit
        year0,month0,day0 = datevec(td0)
        tyr0 = fndatenum_to_decyr(td0)
        print(np.array(['Loaded ',varstrc,' from: ',fname0]))
    if (use_subdelta == 1):
        fname0f1 = np.array(['roho800_',varstr,'_',str(yearmin0),'_',str(yearmax0),'.nc'])
        fname0f = np.array([dir0,fname0f1])
        X0f = np.squeeze(ncread(fname0f,varstrc,np.array([lonsel(1),latsel(1),1]),np.array([nlonsel,nlatsel,Inf])))
        time0f = ncread(fname0f,varstrt)
        td0f = datenum(1948,1,1) + time0f
        year0f,month0f,day0f,hour0f,minutes0f = datevec(td0f)
        hour0fc = hour0f + minutes0f / 60
        yrday0f = td0f - datenum(year0f,1,1)
        tyr0f = fndatenum_to_decyr(td0f)
        print(np.array(['Loaded ',varstrc,' from: ',fname0f]))
    if (do_ocean_bcs == 1):
        ##
#Match each boundary condition location in X0 to a nearest grid point in X1
#This should be consistent with the approach in make_model_inputs.m
        if str(bcstr) == str('_south'):
            lat01 = Lats0(:,1)
            lon01 = Lons0(:,1)
            h0 = H0(:,1)
            m0 = M0(:,1)
        if str(bcstr) == str('_north'):
            lat01 = Lats0(:,end())
            lon01 = Lons0(:,end())
            h0 = H0(:,end())
            m0 = M0(:,end())
        if str(bcstr) == str('_west'):
            lat01 = np.transpose(Lats0(1,:))
            lon01 = np.transpose(Lons0(1,:))
            h0 = np.transpose(H0(1,:))
            m0 = np.transpose(M0(1,:))
        if str(bcstr) == str('_east'):
            lat01 = np.transpose(Lats0(end(),:))
            lon01 = np.transpose(Lons0(end(),:))
            h0 = np.transpose(H0(end(),:))
            m0 = np.transpose(M0(end(),:))
        Tcline,theta_s,theta_b,Vtransform,Vstretching,N,hc = setROMSvgrid('roho800','_v2bg')
        if np.any(str(varstr) == str(np.array(['zeta','ubar','vbar']))):
            N = 1
        sz0 = X0.shape
        n0 = sz0(1)
        nt0 = len(td0)
        nt1 = len(td1)
        X11 = np.squeeze(X1(:,:,1,1))
        isnan11 = np.isnan(X11)
        nx1,ny1 = X11.shape
        ns = N * n0
        #Loop over n0 target profiles: interpolate closest NORESM grid point data over depth, store in X10
        interpolate_ESM = 1
        use_lin = 1
        if np.any(str(varstr) == str(np.array(['zeta','ubar','vbar','u','v']))):
            interpolate_ESM = 0
        print('Interpolating NORESM projections to the ROHO800 boundary grid points')
        extrap_constz = 1
        mindist10 = NaN * np.ones((n0,1))
        X10 = NaN * np.ones((n0,N,nt1))
        z10 = NaN * np.ones((n0,N))
        if (interpolate_ESM == 1):
            Fv = cell(1,nz11)
            for k in np.arange(1,nz11+1).reshape(-1):
                #Horizontal linear interpolation, set triangulation for this NORESM depth level
                v = np.squeeze(X1(:,:,k,1))
                sel1k = find(not np.isnan(v) )
                if (len(sel1k) > 2):
                    x = lon1(sel1k)
                    y = lat1(sel1k)
                    v = v(sel1k)
                    Fv[k] = scatteredInterpolant(x,y,v,'linear','none')
        for j in np.arange(1,n0+1).reshape(-1):
            if (m0(j) == 1):
                #Calculate depth levels
                z01 = roms_zdepth(h0(j),hc,theta_s,theta_b,N,Vtransform,Vstretching)
                z01 = - 1 * z01
                z10[j,:] = np.transpose(z01)
                if (interpolate_ESM == 1):
                    dist10 = fn_spherical_distm(lat1,lon1,lat01(j),lon01(j))
                    dist10m = reshape(dist10,nx1,ny1)
                    mindist10[j] = np.amin(dist10(isnan11 == 0))
                    selxy1 = find(dist10 == np.logical_and(mindist10(j),isnan11) == 0)
                    selxy1 = selxy1(1)
                    selx1,sely1 = find(dist10m == np.logical_and(mindist10(j),isnan11) == 0)
                    selx1 = selx1(1)
                    sely1 = sely1(1)
                    #Linearly interpolate to ROMS depths levels (z01) for point j
                    X11_nn = np.squeeze(X1(selx1,sely1,:,:))
                    if (use_lin == 1):
                        #Get linear scattered interpolants, looping over depth (shallow-to-deep) and time
                        X11_lin = NaN * X11_nn
                        xi = lon01(j)
                        yi = lat01(j)
                        for k in np.arange(1,nz11+1).reshape(-1):
                            #Horizontal linear interpolation, set triangulation for this NORESM depth level
                            v = np.squeeze(X1(:,:,k,1))
                            sel1k = find(not np.isnan(v) )
                            if (len(sel1k) > 2):
                                F = Fv[k]
                                for l in np.arange(1,nt1+1).reshape(-1):
                                    X111 = np.squeeze(X1(:,:,k,l))
                                    F.Values = X111(sel1k)
                                    X11_lin[k,l] = F(xi,yi)
                    X11 = X11_nn
                    if use_lin == 1:
                        sel_lin = find(not np.isnan(X11_lin) )
                        X11[sel_lin] = X11_lin(sel_lin)
                    selz = find(not np.isnan(X11(:,1)) )
                    X10[j,:,:] = Xplin(z01,z11(selz),extrap_constz) * X11(selz,:)
                else:
                    X10[j,:,:] = np.ones((N,nt1))
            if np.mod(j,10) == 0:
                print(np.array(['Done ',str(j),' of ',str(n0),' boundary grid points']))
        print('Done interpolation')
        #Rearrange (X1,X10) into 2D arrays (nt x nseries)
        X0c = X0
        X0 = NaN * np.ones((nt0,ns))
        X1 = NaN * np.ones((nt1,ns))
        lat1 = NaN * np.ones((1,ns))
        lon1 = lat1
        z1 = lat1
        for i in np.arange(1,N+1).reshape(-1):
            for j in np.arange(1,n0+1).reshape(-1):
                if N > 1:
                    X0[:,[i - 1] * n0 + j] = np.squeeze(X0c(j,i,:))
                X1[:,[i - 1] * n0 + j] = np.squeeze(X10(j,i,:))
                lat1[[i - 1] * n0 + j] = lat01(j)
                lon1[[i - 1] * n0 + j] = lon01(j)
                z1[[i - 1] * n0 + j] = z10(j,i)
        if N == 1:
            X0 = np.transpose(X0c)
        clear('X0c','X10')
    else:
        if (ns > 1):
            lat1 = NaN * np.ones((1,ns))
            lon1 = lat1
            X1c = NaN * np.ones((X1.shape[3-1],ns))
            X0c = NaN * np.ones((X0.shape[3-1],ns))
            if use_subdelta == 1:
                X0fc = NaN * np.ones((X0f.shape[3-1],ns))
            for i in np.arange(1,nlonsel+1).reshape(-1):
                for j in np.arange(1,nlatsel+1).reshape(-1):
                    lat1[[i - 1] * nlatsel + j] = lat11(j)
                    lon1[[i - 1] * nlatsel + j] = lon11(i)
                    X1c[:,[i - 1] * nlatsel + j] = np.squeeze(X1(i,j,:))
                    X0c[:,[i - 1] * nlatsel + j] = np.squeeze(X0(i,j,:))
                    if use_subdelta == 1:
                        X0fc[:,[i - 1] * nlatsel + j] = np.squeeze(X0f(i,j,:))
            X1 = X1c
            X0 = X0c
            clear('X1c','X0c')
            if use_subdelta == 1:
                X0f = X0fc
                clear('X0fc')
        else:
            lat1 = lat11
            lon1 = lon11
    if (use_daily == 1 and interp_to_gregorian == 1):
        td1o = td1
        X1o = X1
        td1 = np.transpose((np.arange(datenum(2000,1,2) + dtdmean,datenum(2099,12,31) + dtdmean+1)))
        year1,month1,day1 = datevec(td1)
        yrday1 = td1 - datenum(year1,1,1)
        tyr1 = fndatenum_to_decyr(td1)
        X1 = interp1(td1o,X1o,td1)
    #Correct bogus negative values
    if (str(varstr) == str('rain') == 1):
        Xmin = 1e-09
        X0 = np.amax(Xmin,X0)
        if use_subdelta == 1:
            X0f = np.amax(Xmin,X0f)
    ## Check formula for maximum (clear-sky) irradiance, if needed
    test_maxswrad_formulae = 0
    if (str(varstr) == str('swrad') == 1 and test_maxswrad_formulae == 1):
        sel0f = np.arange(1,len(year0f)+1)
        #sel0f = find(year0f>=2007);
        sels = np.arange(1,ns+1)
        #sels = find(lat1==60.5 & lon1==7);
        #optSWI = struct('ROMS_model',1); #Overshoot ~ 32.6 W/m2
        optSWI = struct('f_model',1,'decl_model',1,'use_eqtime',1,'Q_model',0,'cloud_model',0,'ndays_year',365.2422,'use_tday_lag',1,'year',year0f(sel0f))
        #Note: switching off the Equation-of-time correction (use_eqtime=0) exacerbates overshoot (~14.3 W/m2 at 60.5N,7E; ~32.8 W/m2 over full ROHO800 domain)
#Note: switching off leap year correction (use_tday_lag=0) makes no significant difference (~6.5 W/m2 at 60.5N,7E; ~20.0 W/m2 over full ROHO800 domain)
        maxswrad0f = calc_surfaceSWI(lat1(sels),lon1(sels),yrday0f(sel0f),optSWI)
        np.amax(np.amax(X0f(sel0f,sels) - maxswrad0f))
        #figure;plot(yrday0f,maxswrad0f,'r.',yrday0f,X0f,'k.')
#figure;plot(hour0f,maxswrad0f,'r.',hour0f,X0f,'k.')
        stop
    ## Apply delta change and quantile correction using function fn_biascorrect_projections.m
    tdp = td1
    Xp = X1
    yearp,monthp,dayp = datevec(tdp)
    yrdayp = tdp - datenum(yearp,1,1)
    tyrp = fndatenum_to_decyr(tdp)
    np = len(tyrp)
    opt = struct('seasonal',seasonal,'use_month',use_month,'use_subdelta',use_subdelta,'fractional',fractional,'correct_iavar',correct_iavar,'dyear_iavar',10,'XpcL',XpcL,'XpcH',XpcH,'legacy',legacy)
    if fractional == 1:
        opt.minX = minX
    if use_daily == 1:
        opt.tclimc = tclimc
        opt.tol_dtclim = tol_dtclim
        opt.tol_dtclimp = 1 / 24
    if (str(varstr) == str('swrad') == 1):
        opt.use_XpcH_SWI = use_calc_surfaceSWI
        opt.latp = lat1
        opt.lonp = lon1
        opt.optSWI = struct('f_model',1,'decl_model',1,'use_eqtime',1,'Q_model',0,'cloud_model',0,'ndays_year',365.2422,'use_tday_lag',1)
        opt.use_SWI_hav = 1
    #opt.yearpminv = 2000:20:2080; opt.yearpmaxv = 2019:20:2099;
    opt.yearpminv = np.arange(2000,2090+10,10)
    opt.yearpmaxv = np.arange(2009,2099+10,10)
    if (method == 2):
        opt.yearpminv = np.arange(2000,2090+10,10)
        opt.yearpmaxv = np.arange(2009,2099+10,10)
        #if test_method==1; opt.yearpminv = 2008:6:2026; opt.yearpmaxv = 2013:6:2031; end
    if (use_subdelta == 1):
        #opt.tdpfmax = datenum(2039,12,31);
        opt.Xhf = X0f
        opt.tdhf = td0f
        opt.remove_deltascale_var = remove_deltascale_var
        opt.correct_substd = correct_substd
        opt.correct_subquantiles = correct_subquantiles
        if (correct_substd == 1 or correct_subquantiles == 1):
            opt.opt_subq = opt_subq
        opt.match_subdelta_deltamean = match_subdelta_deltamean
        if match_subdelta_deltamean == 0:
            opt.yearminsub = yearminref
            opt.yearmaxsub = yearmaxref
        opt.match_subdelta_hourly = match_subdelta_hourly
        opt.fractional_subdelta = fractional_subdelta
        if fractional_subdelta == 1:
            opt.minXsub = minXsub
        opt.XpcfL = XpcfL
        opt.XpcfH = XpcfH
        if (str(varstr) == str('swrad') == 1):
            opt.use_XpcfH_SWI = 2 * use_calc_surfaceSWI
            opt.frcrit_capSWI = frcrit_capSWI
        opt.recalculate_Xpc = 1
        #opt.recalculate_Xpc = 0;
    doXpc = 0
    if doXpc == 1:
        Xpc,out = fn_biascorrect_projections(Xp,tdp,X0,td0,yearminref,yearmaxref,0,opt)
    doplot = 0
    if (doplot == 1):
        if use_daily == 1:
            tclim0 = yrday0
            tclimp = yrdayp
            tclimc = np.transpose(opt.tclimc)
        else:
            tclim0 = month0
            tclimp = monthp
            tclimc = np.transpose((np.arange(1,12+1)))
        ntclimc = len(tclimc)
        X0clim = NaN * np.ones((ntclimc,1))
        Xpcclim = X0clim
        Xhclimc = X0clim
        Xpclimc = X0clim
        for i in np.arange(1,ntclimc+1).reshape(-1):
            #These are the climatologies without smoothing, after correction
# REDO        X0clim(i) = mean(mean(X0(year0>=yearminref & year0<=yearmaxref & tclim0==tclimc(i),:)));
# REDO        Xpcclim(i) = mean(mean(Xpc(yearp>=yearminref & yearp<=yearmaxref & tclimp==tclimc(i),:)));
            pass
        plt.figure(20)
        clf
        plt.plot(tclimc,X0clim,'k-',tclimc,Xpcclim,'m-')
    opt2 = opt
    opt2.tol_dtclim = tol_dtclim
    #opt2.qsel = 0.005:0.01:0.995; #This doesn't really help for method=1
    if (method == 2):
        if use_daily == 1:
            opt2.qsel = np.arange(0.005,0.995+0.01,0.01)
        else:
            opt2.qsel = np.arange(0.05,0.95+0.1,0.1)
    Xpc2,out2 = fn_biascorrect_projections(Xp,tdp,X0,td0,yearminref,yearmaxref,method,opt2)
    if (use_subdelta == 1):
        tdpf = out2.tdpf
        yearpf,monthpf,daypf,hourpf,minutespf = datevec(tdpf)
        hourpfc = hourpf + minutespf / 60
        yrdaypf = tdpf - datenum(yearpf,1,1)
        tyrpf = fndatenum_to_decyr(tdpf)
        if doXpc == 1:
            Xpcf = out.Xpcf
        Xpcf2 = out2.Xpcf
        npf = len(tdpf)
        #Check: For e.g. swrad, the character of the diel variability should not be altered by
#       any capping using clear-sky irradiances, see:
#sel0f = find(year0f==2007 & month0f==1 & day0f==1); selpf = find(yearpf==2007 & monthpf==1 & daypf==1);
#sels = 1; [X0f(sel0f,sels) out2.Xpcfo(selpf,sels) Xpcf2(selpf,sels) out2.XpcfH_SWI(selpf,sels)]
        #Also check:
#sum(out2.Xpcfo(:)>Xpcf2(:))/length(Xpcf2(:)) #This should be a small fraction (e.g. <1#)
    ## Check for statistical continuity:
    doplot = 0
    if (doplot == 1):
        ##
        tfontsize = 11 - 2 * small_screen
        lfontsize = 11
        lfontweight = 'bold'
        showtitle = 1
        varstrp = varstr
        if str(varstr) == str('O3_c'):
            varstrp = 'DIC'
        if str(varstr) == str('O3_TA'):
            varstrp = 'TA'
        if str(varstr) == str('O2_o'):
            varstrp = 'DO'
        if str(varstr) == str('N1_p'):
            varstrp = 'PO4'
        if str(varstr) == str('N3_n'):
            varstrp = 'NO3+NO2'
        if str(varstr) == str('N5_s'):
            varstrp = 'Si'
        #sels = 4;
        sels = find(lat1 == np.logical_and(60.5,lon1) == 7)
        if do_ocean_bcs == 1:
            sels = 1 + (25 - 1) * n0
        if ns == 1:
            sels = 1
        if (use_subdelta == 1):
            # REDO        sel0f = find(year0f>=2000);
#selpf = find(yearpf>=2020 & yearpf<=2033);
#selpf = find(yearpf>=2020);
# REDO        selpf = find(yearpf>=2000);
            plt.figure(30 + correct_substd)
            clf
            show_Xpcf = 0
            show_Xpcf2 = 1
            if match_subdelta_deltamean == 1:
                matchstr = ', match_subdelta_deltamean=1'
            else:
                matchstr = []
            plt.plot(tyr0f(sel0f),X0f(sel0f,sels),'k-')
            tstr1 = np.array([varstrp,' hourly at (',num2str(lat1(sels),'%3.1f'),'N,',num2str(lon1(sels),'%3.1f'),'E) from: ERA5 hindcast (black)'])
            if show_Xpcf == 1:
                hold('on')
                plt.plot(tyrpf(selpf),Xpcf(selpf,sels),'m-')
                tstr1 = np.array([tstr1,'; NORESM2-SSP585 with daily delta-change (magenta)'])
            if show_Xpcf2 == 1:
                hold('on')
                plt.plot(tyrpf(selpf),Xpcf2(selpf,sels),'c-')
                tstr1 = np.array([tstr1,'; NORESM2-SSP585 with daily quantile correction to ERA5 (cyan)'])
            tstr1 = np.array([tstr1,tolstr,matchstr])
            plt.axis('tight')
            #if strcmp(varstr,'Tair'); axis([floor(min(tyr0f(sel0f))) ceil(max(tyrpf(selpf))) -40 40]); end
            if showtitle == 1:
                plt.title(tstr1,'FontSize',tfontsize,'FontWeight','bold','Interpreter','none')
            plt.xlabel('Year','FontSize',lfontsize,'FontWeight',lfontweight,'Interpreter','none')
            plt.ylabel(np.array([varstrp,unitstr]),'FontSize',lfontsize,'FontWeight',lfontweight,'Interpreter','none')
            box('on')
        #NOTE: The continuity plot can be deceptive (e.g. for Vwind at (60.5N,7E)), because the variability level in the projection
#      model may actually be lower during the 2000-2020 period compared to the rest of the century, in which case the
#      bias-corrected projection should show an increase in variability between hindcast and projection.
#      This can be confirmed by changing selpf to include 2000-2020 period.
        doplot_deltascale = 1
        if (doplot_deltascale == 1):
            plt.figure(40)
            clf
            show_Xp = 1
            show_Xpc = 0
            show_Xpc2 = 1
            # REDO        sel0 = find(year0>=2000); sel1 = find(yearp>=2000);
#selp = find(yearp>=2020);
# REDO        selp = find(yearp>=2000);
            if do_ocean_bcs == 1:
                zstr = np.array([str(np.round(z1(sels))),'m'])
            else:
                zstr = []
            tstr1 = np.array([varstrp,' ',scstr2,' averages at (',num2str(lat1(sels),'%3.1f'),'N,',num2str(lon1(sels),'%3.1f'),'E,',zstr,') from: ',source0str,' (black)'])
            if show_Xp == 1:
                plt.plot(tyrp(sel1),Xp(sel1,sels),'r-')
                tstr1 = np.array([tstr1,'; NORESM2-SSP585 (red)'])
            hold('on')
            plt.plot(tyr0(sel0),X0(sel0,sels),'k-')
            if show_Xpc == 1:
                hold('on')
                plt.plot(tyrp(selp),Xpc(selp,sels),'m-')
                tstr1 = np.array([tstr1,'; NORESM2-SSP585 with ',scstr2,' delta-change (magenta)'])
            if show_Xpc2 == 1:
                hold('on')
                plt.plot(tyrp(selp),Xpc2(selp,sels),'c-')
                tstr1 = np.array([tstr1,'; NORESM2-SSP585 with ',scstr2,' quantile correction to ',source0str,' (cyan)'])
            #hold on; plot(tyrp(selp),Xpc2A(selp,sels),'b-');
            tstr1 = np.array([tstr1,tolstr])
            plt.axis('tight')
            #set(gca,'XLim',[2007 2030])
            if showtitle == 1:
                plt.title(tstr1,'FontSize',tfontsize,'FontWeight','bold','Interpreter','none')
            plt.xlabel('Year','FontSize',lfontsize,'FontWeight',lfontweight,'Interpreter','none')
            plt.ylabel(np.array([varstrp,unitstr]),'FontSize',lfontsize,'FontWeight',lfontweight,'Interpreter','none')
            box('on')
        #Check: For method=2, we should have exact agreement (within rounding) error during the reference period
#       except for at the boundary days, where the interpolation between daily averages is affected by data
#       outside the reference period (if remove_deltascale_var=1):
#max(abs(X0f(td0f>datenum(2010,1,2) & td0f<datenum(2019,12,30),1)-out2.Xpcf(tdpf>datenum(2010,1,2) & tdpf<datenum(2019,12,30),1))) #Should be within rounding error
    ## Make monthly averages if req'd
    domav = 0
    if (domav == 1):
        if (use_daily == 1):
            #year0min = min(year0); year0max = max(year0);
            year0min = 2007
            year0max = 2019
            nyrs0 = year0max - year0min + 1
            n0m = nyrs0 * 12
            td0m = NaN * np.ones((n0m,1))
            X0m = NaN * np.ones((n0m,ns))
            X0mstd = X0m
            dX0sub = NaN * X0
            dX0subf = NaN * X0
            for i in np.arange(1,nyrs0+1).reshape(-1):
                for j in np.arange(1,12+1).reshape(-1):
                    sel01 = find(year0 == np.logical_and((year0min + (i - 1)),month0) == j)
                    X0m[[i - 1] * 12 + j,:] = meanan(X0(sel01,:))
                    X0mstd[[i - 1] * 12 + j,:] = stdan(X0(sel01,:))
                    dX0sub[sel01,:] = X0(sel01,:) - (np.ones((len(sel01),1)) * X0m((i - 1) * 12 + j,:))
                    dX0subf[sel01,:] = X0(sel01,:) / (np.ones((len(sel01),1)) * X0m((i - 1) * 12 + j,:))
                    td0m[[i - 1] * 12 + j] = datenum(year0min + (i - 1),j,15)
            year0m,month0m,day0m = datevec(td0m)
            tyr0m = fndatenum_to_decyr(td0m)
            yearpmin = np.amin(yearp)
            yearpmax = np.amax(yearp)
            nyrsp = yearpmax - yearpmin + 1
            npm = nyrsp * 12
            tdpm = NaN * np.ones((npm,1))
            Xpcm = NaN * np.ones((npm,ns))
            Xpcmstd = Xpcm
            Xpc2m = Xpcm
            Xpc2mstd = Xpcm
            for i in np.arange(1,nyrsp+1).reshape(-1):
                for j in np.arange(1,12+1).reshape(-1):
                    selp1 = find(yearp == np.logical_and((yearpmin + (i - 1)),monthp) == j)
                    Xpcm[[i - 1] * 12 + j,:] = meanan(Xpc(selp1,:))
                    Xpcmstd[[i - 1] * 12 + j,:] = stdan(Xpc(selp1,:))
                    Xpc2m[[i - 1] * 12 + j,:] = meanan(Xpc2(selp1,:))
                    Xpc2mstd[[i - 1] * 12 + j,:] = stdan(Xpc2(selp1,:))
                    tdpm[[i - 1] * 12 + j] = datenum(yearpmin + (i - 1),j,15)
            yearpm,monthpm,daypm = datevec(tdpm)
            tyrpm = fndatenum_to_decyr(tdpm)
    ## Make decadal monthly climatological statistics if req'd
    do_clim = 0
    if (do_clim == 1):
        print('Computing climatologies')
        tic
        X0_clim = NaN * np.ones((12,ns))
        X0_climstd = X0_clim
        X0_climL = X0_clim
        X0_climH = X0_clim
        X1_clim = X0_clim
        X1_climstd = X0_clim
        X1_climL = X0_clim
        X1_climH = X0_clim
        if doXpc == 1:
            Xpc_clim = X0_clim
            Xpc_climstd = X0_clim
            Xpc_climL = X0_clim
            Xpc_climH = X0_clim
        Xpc2_clim = X0_clim
        Xpc2_climstd = X0_clim
        Xpc2_climL = X0_clim
        Xpc2_climH = X0_clim
        if (use_subdelta == 1):
            X0f_clim = X0_clim
            X0f_climstd = X0_clim
            X0f_climL = X0_clim
            X0f_climH = X0_clim
            if doXpc == 1:
                Xpcf_clim = X0_clim
                Xpcf_climstd = X0_clim
                Xpcf_climL = X0_clim
                Xpcf_climH = X0_clim
            Xpcf2_clim = X0_clim
            Xpcf2_climstd = X0_clim
            Xpcf2_climL = X0_clim
            Xpcf2_climH = X0_clim
        yearminrefc = 2010
        yearmaxrefc = 2019
        qsel = np.array([0.025,0.975])
        # REDO    X0_meanm = mean(X0(year0>=yearminrefc & year0<=yearmaxrefc,:));
# REDO    X0_stdm = std(X0(year0>=yearminrefc & year0<=yearmaxrefc,:));
        if (use_subdelta == 1):
            # REDO        X0f_meanm = mean(X0f(year0f>=yearminrefc & year0f<=yearmaxrefc,:));
# REDO        X0f_stdm = std(X0f(year0f>=yearminrefc & year0f<=yearmaxrefc,:));
            pass
        for i in np.arange(1,12+1).reshape(-1):
            # REDO        selt0 = find(year0>=yearminrefc & year0<=yearmaxrefc & month0==i);
            X0_clim[i,:] = mean(X0(selt0,:))
            X0_climstd[i,:] = std(X0(selt0,:))
            X0_climL[i,:] = quantile(X0(selt0,:),qsel(1))
            X0_climH[i,:] = quantile(X0(selt0,:),qsel(2))
            # REDO        selt1 = find(year1>=yearminrefc & year1<=yearmaxrefc & month1==i);
            X1_clim[i,:] = mean(X1(selt1,:))
            X1_climstd[i,:] = std(X1(selt1,:))
            X1_climL[i,:] = quantile(X1(selt1,:),qsel(1))
            X1_climH[i,:] = quantile(X1(selt1,:),qsel(2))
            # REDO        seltp = find(yearp>=yearminrefc & yearp<=yearmaxrefc & monthp==i);
            if (doXpc == 1):
                Xpc_clim[i,:] = mean(Xpc(seltp,:))
                Xpc_climstd[i,:] = std(Xpc(seltp,:))
                Xpc_climL[i,:] = quantile(Xpc(seltp,:),qsel(1))
                Xpc_climH[i,:] = quantile(Xpc(seltp,:),qsel(2))
            Xpc2_clim[i,:] = mean(Xpc2(seltp,:))
            Xpc2_climstd[i,:] = std(Xpc2(seltp,:))
            Xpc2_climL[i,:] = quantile(Xpc2(seltp,:),qsel(1))
            Xpc2_climH[i,:] = quantile(Xpc2(seltp,:),qsel(2))
            if (use_subdelta == 1):
                # REDO            selt0f = find(year0f>=yearminrefc & year0f<=yearmaxrefc & month0f==i);
                X0f_clim[i,:] = mean(X0f(selt0f,:))
                X0f_climstd[i,:] = std(X0f(selt0f,:))
                X0f_climL[i,:] = quantile(X0f(selt0f,:),qsel(1))
                X0f_climH[i,:] = quantile(X0f(selt0f,:),qsel(2))
                # REDO            seltpf = find(yearpf>=yearminrefc & yearpf<=yearmaxrefc & monthpf==i);
                if (doXpc == 1):
                    Xpcf_clim[i,:] = mean(Xpcf(seltpf,:))
                    Xpcf_climstd[i,:] = std(Xpcf(seltpf,:))
                    Xpcf_climL[i,:] = quantile(Xpcf(seltpf,:),qsel(1))
                    Xpcf_climH[i,:] = quantile(Xpcf(seltpf,:),qsel(2))
                Xpcf2_clim[i,:] = mean(Xpcf2(seltpf,:))
                Xpcf2_climstd[i,:] = std(Xpcf2(seltpf,:))
                Xpcf2_climL[i,:] = quantile(Xpcf2(seltpf,:),qsel(1))
                Xpcf2_climH[i,:] = quantile(Xpcf2(seltpf,:),qsel(2))
        yearminv = np.arange(2010,2090+10,10)
        yearmaxv = np.arange(2019,2099+10,10)
        nperiods = len(yearminv)
        tyr_clim = NaN * np.ones((1,nperiods))
        X1_meanm = NaN * np.ones((nperiods,ns))
        X1_stdm = X1_meanm
        if doXpc == 1:
            Xpc_meanm = X1_meanm
            Xpc_stdm = X1_meanm
        Xpc2_meanm = X1_meanm
        Xpc2_stdm = X1_meanm
        if doXpc == 1:
            Xpcf_meanm = X1_meanm
            Xpcf_stdm = X1_meanm
        Xpcf2_meanm = X1_meanm
        Xpcf2_stdm = X1_meanm
        X1_climm = NaN * np.ones((12,nperiods,ns))
        X1_climstdm = X1_climm
        X1_climLm = X1_climm
        X1_climHm = X1_climm
        if doXpc == 1:
            Xpc_climm = X1_climm
            Xpc_climstdm = X1_climm
            Xpc_climLm = X1_climm
            Xpc_climHm = X1_climm
        Xpc2_climm = X1_climm
        Xpc2_climstdm = X1_climm
        Xpc2_climLm = X1_climm
        Xpc2_climHm = X1_climm
        if (use_subdelta == 1):
            if doXpc == 1:
                Xpcf_climm = X1_climm
                Xpcf_climstdm = X1_climm
                Xpcf_climLm = X1_climm
                Xpcf_climHm = X1_climm
            Xpcf2_climm = X1_climm
            Xpcf2_climstdm = X1_climm
            Xpcf2_climLm = X1_climm
            Xpcf2_climHm = X1_climm
        for j in np.arange(1,nperiods+1).reshape(-1):
            tyr_clim[j] = mean(np.arange(yearminv(j) + 0.5,yearmaxv(j) + 0.5+1))
            # REDO        selt1 = find(year1>=yearminv(j) & year1<=yearmaxv(j));
            X1_meanm[j,:] = mean(X1(selt1,:))
            X1_stdm[j,:] = std(X1(selt1,:))
            # REDO        seltp = find(yearp>=yearminv(j) & yearp<=yearmaxv(j));
            if doXpc == 1:
                Xpc_meanm[j,:] = mean(Xpc(seltp,:))
                Xpc_stdm[j,:] = std(Xpc(seltp,:))
            Xpc2_meanm[j,:] = mean(Xpc2(seltp,:))
            Xpc2_stdm[j,:] = std(Xpc2(seltp,:))
            if (use_subdelta == 1):
                # REDO            seltpf = find(yearpf>=yearminv(j) & yearpf<=yearmaxv(j));
                if doXpc == 1:
                    Xpcf_meanm[j,:] = mean(Xpcf(seltpf,:))
                    Xpcf_stdm[j,:] = std(Xpcf(seltpf,:))
                Xpcf2_meanm[j,:] = mean(Xpcf2(seltpf,:))
                Xpcf2_stdm[j,:] = std(Xpcf2(seltpf,:))
            for i in np.arange(1,12+1).reshape(-1):
                # REDO            selt1 = find(year1>=yearminv(j) & year1<=yearmaxv(j) & month1==i);
                X1_climm[i,j,:] = mean(X1(selt1,:))
                X1_climstdm[i,j,:] = std(X1(selt1,:))
                X1_climLm[i,j,:] = quantile(X1(selt1,:),qsel(1))
                X1_climHm[i,j,:] = quantile(X1(selt1,:),qsel(2))
                # REDO            seltp = find(yearp>=yearminv(j) & yearp<=yearmaxv(j) & monthp==i);
                if doXpc == 1:
                    Xpc_climm[i,j,:] = mean(Xpc(seltp,:))
                    Xpc_climstdm[i,j,:] = std(Xpc(seltp,:))
                    Xpc_climLm[i,j,:] = quantile(Xpc(seltp,:),qsel(1))
                    Xpc_climHm[i,j,:] = quantile(Xpc(seltp,:),qsel(2))
                Xpc2_climm[i,j,:] = mean(Xpc2(seltp,:))
                Xpc2_climstdm[i,j,:] = std(Xpc2(seltp,:))
                Xpc2_climLm[i,j,:] = quantile(Xpc2(seltp,:),qsel(1))
                Xpc2_climHm[i,j,:] = quantile(Xpc2(seltp,:),qsel(2))
                if (use_subdelta == 1):
                    # REDO                seltpf = find(yearpf>=yearminv(j) & yearpf<=yearmaxv(j) & monthpf==i);
                    if doXpc == 1:
                        Xpcf_climm[i,j,:] = mean(Xpcf(seltpf,:))
                        Xpcf_climstdm[i,j,:] = std(Xpcf(seltpf,:))
                        Xpcf_climLm[i,j,:] = quantile(Xpcf(seltpf,:),qsel(1))
                        Xpcf_climHm[i,j,:] = quantile(Xpcf(seltpf,:),qsel(2))
                    Xpcf2_climm[i,j,:] = mean(Xpcf2(seltpf,:))
                    Xpcf2_climstdm[i,j,:] = std(Xpcf2(seltpf,:))
                    Xpcf2_climLm[i,j,:] = quantile(Xpcf2(seltpf,:),qsel(1))
                    Xpcf2_climHm[i,j,:] = quantile(Xpcf2(seltpf,:),qsel(2))
        print('Done climatologies')
        toc
        #Plot hindcast and projection decadal climatologies if req'd
        doplotclims = 0
        if (doplotclims == 1):
            ##
            plt.figure(3)
            clf
            show_proj = 1
            show_qs = 1
            #fac = -273.15;
            fac = 0
            xlim = np.array([0.5,12.5])
            xlv = np.array([3.6,5])
            xlv2 = np.array([7,8])
            nperiodsc = 5
            #sels = 1;
#sels = find(lat1==60.5 & lon1==7);
            sels = np.arange(1,ns+1)
            if ns == 1:
                sels = 1
            if str(varstr) == str('Tair'):
                ylim = np.array([- 25,30])
                ytop = - 12
                dy = 2.5
                ytop2 = 23
            if (str(varstr) == str('Tair') and ns > 1):
                ylim = np.array([- 13,25])
                ytop = - 3
                dy = 2
                ytop2 = 23
            if str(varstr) == str('Pair'):
                ylim = np.array([0.88,1.06]) * 100000.0
                ytop = 93000.0
                dy = 1000.0
                ytop2 = 105000.0
            if str(varstr) == str('Qair'):
                xlv = np.array([4.6,6])
                ylim = np.array([0,12]) * 0.001
                ytop = 0.0025
                dy = 0.0005
                ytop2 = 0.011
            if str(varstr) == str('cloud'):
                ylim = np.array([- 0.1,1.5])
                ytop = 1.4
                dy = 0.07
                ytop2 = 1.05
            if str(varstr) == str('Uwind') or str(varstr) == str('Vwind'):
                ylim = np.array([- 15,15])
                ytop = - 8
                dy = 1.3
                ytop2 = 12
            if str(varstr) == str('swrad'):
                ylim = np.array([- 50,500])
                ytop = 450
                dy = 20
                ytop2 = 480
            if str(varstr) == str('lwrad'):
                xlv = np.array([4.6,6])
                ylim = np.array([140,400])
                ytop = 190
                dy = 10
                ytop2 = 390
            if str(varstr) == str('rain'):
                ylim = np.array([- 0.1,5]) * 0.0001
                ytop = 0.00045
                dy = 2e-05
                ytop2 = 0.00047
            if (show_qs == 1):
                X0_errL = X0_clim(:,sels) - X0_climL(:,sels)
                X0_errH = X0_climH(:,sels) - X0_clim(:,sels)
                X1_errL = X1_clim(:,sels) - X1_climL(:,sels)
                X1_errH = X1_climH(:,sels) - X1_clim(:,sels)
                if doXpc == 1:
                    Xpc_errL = Xpc_clim(:,sels) - Xpc_climL(:,sels)
                    Xpc_errH = Xpc_climH(:,sels) - Xpc_clim(:,sels)
                Xpc2_errL = Xpc2_clim(:,sels) - Xpc2_climL(:,sels)
                Xpc2_errH = Xpc2_climH(:,sels) - Xpc2_clim(:,sels)
                X1_errLm = mean(X1_climm(:,:,sels) - X1_climLm(:,:,sels),3)
                X1_errHm = mean(X1_climHm(:,:,sels) - X1_climm(:,:,sels),3)
                if doXpc == 1:
                    Xpc_errLm = mean(Xpc_climm(:,:,sels) - Xpc_climLm(:,:,sels),3)
                    Xpc_errHm = mean(Xpc_climHm(:,:,sels) - Xpc_climm(:,:,sels),3)
                Xpc2_errLm = mean(Xpc2_climm(:,:,sels) - Xpc2_climLm(:,:,sels),3)
                Xpc2_errHm = mean(Xpc2_climHm(:,:,sels) - Xpc2_climm(:,:,sels),3)
                errstr = np.array(['means\pm',str(100 * np.diff(qsel)),'%CIs from ',scstr2,' averages'])
            else:
                X0_errL = X0_climstd(:,sels)
                X0_errH = X0_climstd(:,sels)
                X1_errL = X1_climstd(:,sels)
                X1_errH = X1_climstd(:,sels)
                if doXpc == 1:
                    Xpc_errL = Xpc_climstd(:,sels)
                    Xpc_errH = Xpc_climstd(:,sels)
                Xpc2_errL = Xpc2_climstd(:,sels)
                Xpc2_errH = Xpc2_climstd(:,sels)
                X1_errLm = mean(X1_climstdm(:,:,sels),3)
                X1_errHm = mean(X1_climstdm(:,:,sels),3)
                if doXpc == 1:
                    Xpc_errLm = mean(Xpc_climstdm(:,:,sels),3)
                    Xpc_errHm = mean(Xpc_climstdm(:,:,sels),3)
                Xpc2_errLm = mean(Xpc2_climstdm(:,:,sels),3)
                Xpc2_errHm = mean(Xpc2_climstdm(:,:,sels),3)
                errstr = np.array(['means\pm1stdevs from ',scstr2,' averages'])
            X1_clim1m = mean(X1_climm(:,:,sels),3)
            if doXpc == 1:
                Xpc_clim1m = mean(Xpc_climm(:,:,sels),3)
            Xpc2_clim1m = mean(Xpc2_climm(:,:,sels),3)
            colm = jet(nperiodsc)
            x1 = np.transpose((np.arange(1,12+1)))
            optl = struct('lfontsize',9)
            xfac = 1.2
            yfac = 1.2
            xoffv = np.array([- 0.02,0.02])
            xstr = 'Month of year'
            ystr = np.array([varstr,unitstr])
            lfontsize = 13
            lfontweight = 'bold'
            fontsize = 12
            subplot(2,2,1)
            errorbar(x1 - 0.3,mean(X0_clim(:,sels),2) + fac,mean(X0_errL,2),mean(X0_errH,2),'ko','Color','k','LineWidth',2)
            hold('on')
            errorbar(x1 - 0.1,mean(X1_clim(:,sels),2) + fac,mean(X1_errL,2),mean(X1_errH,2),'^','Color','r','LineWidth',2)
            if doXpc == 1:
                hold('on')
                errorbar(x1 + 0.1,mean(Xpc_clim(:,sels),2) + fac,mean(Xpc_errL,2),mean(Xpc_errH,2),'^','Color','m','LineWidth',2)
            hold('on')
            errorbar(x1 + 0.3,mean(Xpc2_clim(:,sels),2) + fac,mean(Xpc2_errL,2),mean(Xpc2_errH,2),'^','Color','c','LineWidth',2)
            legend2(xlv2,ytop2,dy,np.array([errstr]),[],np.array(['k']),optl)
            #legend2([3.6 5],-5,2.2,{'ERA5 hindcast 2010-2019','NorESM2 2010-2019',...
            if (doXpc == 1):
                legstrv = np.array(['ERA5 hindcast 2010-2019','NorESM2 2010-2019',np.array(['NorESM2-ERA5 2010-2019, DC',tolstr]),np.array(['NorESM2-ERA5 2010-2019, QC',tolstr])])
                colv = np.array(['k','r','m','c'])
            else:
                legstrv = np.array(['ERA5 hindcast 2010-2019','NorESM2 2010-2019',np.array(['NorESM2-ERA5 2010-2019, QC',tolstr])])
                colv = np.array(['k','r','c'])
            legend2(xlv,ytop,dy,legstrv,[],colv,optl)
            plt.axis(np.array([xlim,ylim]))
            resizesubplot(gca,xfac,yfac,xoffv(1),0.03)
            plt.xlabel(xstr,'FontSize',lfontsize,'FontWeight',lfontweight)
            plt.ylabel(ystr,'FontSize',lfontsize,'FontWeight',lfontweight)
            box('on')
            if (show_proj == 1):
                dxv = (np.arange(- 0.4,- 0.4 + 0.1 * (nperiodsc - 1)+0.1,0.1))
                subplot(2,2,2)
                for j in np.arange(1,nperiodsc+1).reshape(-1):
                    hold('on')
                    errorbar(x1 + dxv(j),X1_clim1m(:,j) + fac,X1_errLm(:,j),X1_errHm(:,j),'^-','Color',colm(j,:),'LineWidth',1.5)
                strm = cell(nperiodsc,1)
                for j in np.arange(1,nperiodsc+1).reshape(-1):
                    strm[j] = np.array(['NorESM2 SSP585 ',str(yearminv(j)),'-',str(yearmaxv(j))])
                legend2(xlv,ytop,dy,strm,[],colm,optl)
                plt.axis(np.array([xlim,ylim]))
                resizesubplot(gca,xfac,yfac,xoffv(2),0.03)
                plt.xlabel(xstr,'FontSize',lfontsize,'FontWeight',lfontweight)
                plt.ylabel(ystr,'FontSize',lfontsize,'FontWeight',lfontweight)
                box('on')
                if (doXpc == 1):
                    subplot(2,2,3)
                    for j in np.arange(1,nperiodsc+1).reshape(-1):
                        hold('on')
                        errorbar(x1 + dxv(j),Xpc_clim1m(:,j) + fac,Xpc_errLm(:,j),Xpc_errHm(:,j),'^-','Color',colm(j,:),'LineWidth',1.5)
                    strm = cell(nperiodsc,1)
                    for j in np.arange(1,nperiodsc+1).reshape(-1):
                        strm[j] = np.array(['NorESM2-ERA5 SSP585 ',str(yearminv(j)),'-',str(yearmaxv(j)),', DC',tolstr])
                    legend2(xlv,ytop,dy,strm,[],colm,optl)
                    plt.axis(np.array([xlim,ylim]))
                    resizesubplot(gca,xfac,yfac,xoffv(1),0)
                    plt.xlabel(xstr,'FontSize',lfontsize,'FontWeight',lfontweight)
                    plt.ylabel(ystr,'FontSize',lfontsize,'FontWeight',lfontweight)
                    box('on')
                subplot(2,2,4)
                for j in np.arange(1,nperiodsc+1).reshape(-1):
                    hold('on')
                    errorbar(x1 + dxv(j),Xpc2_clim1m(:,j) + fac,Xpc2_errLm(:,j),Xpc2_errHm(:,j),'^-','Color',colm(j,:),'LineWidth',1.5)
                strm = cell(nperiodsc,1)
                for j in np.arange(1,nperiodsc+1).reshape(-1):
                    strm[j] = np.array(['NorESM2-ERA5 SSP585 ',str(yearminv(j)),'-',str(yearmaxv(j)),', QC',tolstr])
                legend2(xlv,ytop,dy,strm,[],colm,optl)
                plt.axis(np.array([xlim,ylim]))
                resizesubplot(gca,xfac,yfac,xoffv(2),0)
                plt.xlabel(xstr,'FontSize',lfontsize,'FontWeight',lfontweight)
                plt.ylabel(ystr,'FontSize',lfontsize,'FontWeight',lfontweight)
                box('on')
    ## Plot data by month if req'd
    doplot_bymonth = 0
    if (doplot_bymonth == 1):
        ##
        nrows = 3
        ncols = 4
        xstr = 'Year'
        lfontsize = 9
        #if strcmp(varstr,'Tair'); ystr = 'T [\circC]'; end
        ystr = np.array([varstr,unitstr])
        tfontsize = 8
        tfontweight = 'bold'
        fontsize = 9
        sgfontsize = 9
        col0 = 'k'
        colp = 'm'
        colp2 = 'c'
        fix_ylim = 2
        show_Xpc = 1
        if doXpc == 0:
            show_Xpc = 0
        show_Xpc2 = 1
        show_hourly = 1
        show_daily = 1
        if use_daily == 0:
            show_hourly = 0
            show_daily = 0
        show_skill = 1
        fmtstr = '%3.3f'
        if str(varstr) == str('rain'):
            fmtstr = '%3.4f'
        if str(varstr) == str('Qair'):
            fmtstr = '%3.6f'
        #sels = 8;
#sels = find(lat1==60.5 & lon1==7);
        sels = np.arange(1,ns+1)
        if ns == 1:
            sels = 1
        if use_daily == 0:
            show_daily = 0
        if str(varstr) == str('rain'):
            fac1 = 1000.0
            unitstrc = ' [g/m2/s]'
        else:
            fac1 = 1
            unitstrc = unitstr
        #if contains(varstr,'_'); interpstr = 'none'; else interpstr = 'tex'; end
        interpstr = 'tex'
        varstrc = varstr(varstr != '_')
        if str(varstr) == str('lwrad_down'):
            varstrc = 'lwrad'
        plot_type = 1
        #1 for quantile-quantile plots (hourly if show_hourly=1)
#2 for decadal means from daily/monthly averages
#3 for decadal standard deviations from daily/monthly averages
#4 for decadal means from hourly data
#5 for decadal standard deviations from hourly data
#6 for diel variability from hourly data (means +/- CIs)
        aggregation_type = 0
        if plot_type == 1:
            aggregation_type = 2
        match_spatially = 0
        #yearmin1 = yearminref; yearmax1 = yearmaxref;
        yearmin1 = 2000
        yearmax1 = 2009
        #if (do_ocean_bcs==1 && test_method==1); yearmin1 = 2008; yearmax1 = 2014; end
        xlim = np.array([2007,2040])
        # REDO    if (plot_type>=2 && plot_type<=5); xlim = [2010 2100]; end
        if plot_type == 4 or plot_type == 5:
            show_hourly = 1
        if plot_type == 2 or plot_type == 3:
            show_hourly = 0
        if plot_type == 6:
            show_hourly = 1
            aggregation_type = 2
        if (use_daily == 1):
            if (show_hourly == 1):
                tyr01 = tyr0f
                X01 = X0f(:,sels)
                year01 = year0f
                yrday01 = yrday0f
                month01 = month0f
                tyrp1 = tyrpf
                yearp1 = yearpf
                yrdayp1 = yrdaypf
                monthp1 = monthpf
                yrdaysel = np.array([15,45,75,105,135,165,195,225,255,285,315,345])
                if show_Xpc == 1:
                    Xpc1 = Xpcf(:,sels)
                if show_Xpc2 == 1:
                    Xpc21 = Xpcf2(:,sels)
                lwidth0 = 1
            else:
                if (show_daily == 1):
                    tyr01 = tyr0
                    X01 = X0(:,sels)
                    year01 = year0
                    yrday01 = yrday0
                    month01 = month0
                    tyrp1 = tyrp
                    yearp1 = yearp
                    yrdayp1 = yrdayp
                    monthp1 = monthp
                    yrdaysel = np.array([15,45,75,105,135,165,195,225,255,285,315,345])
                    if show_Xpc == 1:
                        Xpc1 = Xpc(:,sels)
                    if show_Xpc2 == 1:
                        Xpc21 = Xpc2(:,sels)
                    lwidth0 = 1
                else:
                    tyr01 = tyr0m
                    X01 = X0m(:,sels)
                    year01 = year0m
                    month01 = month0m
                    tyrp1 = tyrpm
                    yearp1 = yearpm
                    monthp1 = monthpm
                    if show_Xpc == 1:
                        Xpc1 = Xpcm(:,sels)
                    if show_Xpc2 == 1:
                        Xpc21 = Xpc2m(:,sels)
                    lwidth0 = 1.5
        else:
            tyr01 = tyr0
            X01 = X0(:,sels)
            year01 = year0
            month01 = month0
            tyrp1 = tyrp
            yearp1 = yearp
            monthp1 = monthp
            if show_Xpc == 1:
                Xpc21 = Xpc2(:,sels)
            if show_Xpc2 == 1:
                Xpc21 = Xpc2(:,sels)
            lwidth0 = 1.5
        plt.figure(1)
        clf
        tstr1 = []
        if fix_ylim == 1:
            ylim = np.array([np.amin(np.array([[X01],[Xpc1]])),np.amax(np.array([[X01],[Xpc1]]))])
        MAEv = NaN * np.ones((1,12))
        RMSEv = MAEv
        for j in np.arange(1,12+1).reshape(-1):
            subplot(nrows,ncols,j)
            if (show_daily == 1):
                if aggregation_type == 0:
                    sel0 = find(yrday01 == yrdaysel(j))
                    selp = find(yrdayp1 == yrdaysel(j))
                else:
                    if aggregation_type == 1:
                        dist1 = fnperiodic_distm(yrday01,yrdaysel(j),365)
                        # REDO                sel0 = find(dist1<=tol_dtclim);
                        dist1p = fnperiodic_distm(yrdayp1,yrdaysel(j),365)
                        # REDO                selp = find(dist1p<=tol_dtclim);
                    else:
                        if aggregation_type == 2:
                            sel0 = find(month01 == j)
                            selp = find(monthp1 == j)
                        else:
                            stop
            else:
                sel0 = find(month01 == j)
                selp = find(monthp1 == j)
            #if (fix_ylim==2 && show_Xpc==1); ylim = [min([X01(sel0);Xpc1(selp)]) max([X01(sel0);Xpc1(selp)])]; end
            if (plot_type == 0):
                plt.plot(tyr01(sel0),X01(sel0),'Color',col0,'LineStyle','-','LineWidth',lwidth0)
                hold('on')
                plt.plot(tyrp1(selp),Xpc1(selp),'Color',colp,'LineStyle','-')
                if show_Xpc2 == 1:
                    hold('on')
                    plt.plot(tyrp1(selp),Xpc21(selp),'Color',colp2,'LineStyle','-')
                set(gca,'XLim',xlim,'FontSize',fontsize)
                if fix_ylim == 1:
                    set(gca,'Ylim',ylim)
                if j > 8:
                    plt.xlabel(xstr,'FontSize',lfontsize)
                #if mod(j,4)==1; ylabel(ystr,'FontSize',lfontsize); end
                plt.ylabel(ystr,'FontSize',lfontsize)
            else:
                if (plot_type == 1):
                    # REDO            sel01 = sel0(year01(sel0)>=yearmin1 & year01(sel0)<=yearmax1 & ~isnan(X01(sel0,1)));
# REDO            selp1 = selp(yearp1(selp)>=yearmin1 & yearp1(selp)<=yearmax1 & ~isnan(Xpc21(selp,1)));
                    nsel1 = np.amin(len(sel01),len(selp1))
                    sel01 = sel01(np.arange(1,nsel1+1))
                    selp1 = selp1(np.arange(1,nsel1+1))
                    #if length(selp1)<length(sel01); selp1 = [selp1;selp1(1:(length(sel01)-length(selp1)))]; end ##ok<AGROW>
                    xlim = []
                    #xlim = [0 100];
                    X011 = fac1 * X01(sel01,:)
                    if match_spatially == 1:
                        X011 = __builtint__.sorted(X011)
                        X011 = X011
                    else:
                        X011 = __builtint__.sorted(X011)
                    Xf = X011
                    if (show_Xpc == 1):
                        Xpc11 = fac1 * Xpc1(selp1,:)
                        if match_spatially == 1:
                            Xpc11 = __builtint__.sorted(Xpc11)
                            Xpc11 = Xpc11
                        else:
                            Xpc11 = __builtint__.sorted(Xpc11)
                        Xf = np.array([[Xf],[Xpc11]])
                    if (show_Xpc2 == 1):
                        Xpc211 = fac1 * Xpc21(selp1,:)
                        if match_spatially == 1:
                            Xpc211 = __builtint__.sorted(Xpc211)
                            Xpc211 = Xpc211
                        else:
                            Xpc211 = __builtint__.sorted(Xpc211)
                        Xf = np.array([[Xf],[Xpc211]])
                    if len(xlim)==0:
                        xlim = np.array([np.amin(Xf),np.amax(Xf)])
                    plt.plot(xlim,xlim,'-','Color',np.array([0.7,0.7,0.7]))
                    if show_Xpc == 1:
                        hold('on')
                        plt.plot(Xpc11,X011,'.','Color',colp)
                    if show_Xpc2 == 1:
                        hold('on')
                        plt.plot(Xpc211,X011,'.','Color',colp2)
                    plt.axis(np.array([xlim,xlim]))
                    set(gca,'FontSize',fontsize)
                    plt.ylabel(np.array(['q',varstrc,' (',source0str,')',unitstrc]),'FontSize',lfontsize,'Interpreter',interpstr)
                    plt.xlabel(np.array(['q',varstrc,' (NORESM-',source0str,')',unitstrc]),'FontSize',lfontsize,'Interpreter',interpstr)
                    MAEv[j] = mean(np.abs(Xpc211 - X011))
                    RMSEv[j] = np.sqrt(mean((Xpc211 - X011) ** 2))
                    #if show_skill==1; tstr1 = [', MAE=',num2str(MAEv(j),fmtstr),unitstrc,', n=',str(length(X011))]; end
                    if show_skill == 1:
                        tstr1 = np.array([', RMSE=',num2str(RMSEv(j),fmtstr),unitstrc,', n=',str(len(X011))])
                else:
                    if (plot_type == 2):
                        plt.plot(tyr_clim(1),mean(X0_clim(j,sels),2),'+','Color',col0,'LineWidth',lwidth0)
                        if show_Xpc == 1:
                            hold('on')
                            plt.plot(tyr_clim,mean(np.squeeze(Xpc_climm(j,:,sels)),2),'Color',colp,'LineStyle','-','LineWidth',lwidth0)
                        if show_Xpc2 == 1:
                            hold('on')
                            plt.plot(tyr_clim,mean(np.squeeze(Xpc2_climm(j,:,sels)),2),'Color',colp2,'LineStyle','-','LineWidth',lwidth0)
                        set(gca,'XLim',xlim,'FontSize',fontsize)
                        if fix_ylim == 1:
                            set(gca,'Ylim',ylim)
                        if j > 8:
                            plt.xlabel(xstr,'FontSize',lfontsize)
                        plt.ylabel(ystr,'FontSize',lfontsize)
                    else:
                        if (plot_type == 3):
                            plt.plot(tyr_clim(1),mean(X0_climstd(j,sels),2),'+','Color',col0,'LineWidth',lwidth0)
                            if show_Xpc == 1:
                                hold('on')
                                plt.plot(tyr_clim,mean(np.squeeze(Xpc_climstdm(j,:,sels)),2),'Color',colp,'LineStyle','-','LineWidth',lwidth0)
                            if show_Xpc2 == 1:
                                hold('on')
                                plt.plot(tyr_clim,mean(np.squeeze(Xpc2_climstdm(j,:,sels)),2),'Color',colp2,'LineStyle','-','LineWidth',lwidth0)
                            set(gca,'XLim',xlim,'FontSize',fontsize)
                            if fix_ylim == 1:
                                set(gca,'Ylim',ylim)
                            if j > 8:
                                plt.xlabel(xstr,'FontSize',lfontsize)
                            plt.ylabel(ystr,'FontSize',lfontsize)
                        else:
                            if (plot_type == 4):
                                plt.plot(tyr_clim(1),mean(X0f_clim(j,sels),2),'+','Color',col0,'LineWidth',lwidth0)
                                if show_Xpc == 1:
                                    hold('on')
                                    plt.plot(tyr_clim,mean(np.squeeze(Xpcf_climm(j,:,sels)),2),'Color',colp,'LineStyle','-','LineWidth',lwidth0)
                                if show_Xpc2 == 1:
                                    hold('on')
                                    plt.plot(tyr_clim,mean(np.squeeze(Xpcf2_climm(j,:,sels)),2),'Color',colp2,'LineStyle','-','LineWidth',lwidth0)
                                set(gca,'XLim',xlim,'FontSize',fontsize)
                                if fix_ylim == 1:
                                    set(gca,'Ylim',ylim)
                                if j > 8:
                                    plt.xlabel(xstr,'FontSize',lfontsize)
                                plt.ylabel(ystr,'FontSize',lfontsize)
                            else:
                                if (plot_type == 5):
                                    plt.plot(tyr_clim(1),mean(X0f_climstd(j,sels),2),'+','Color',col0,'LineWidth',lwidth0)
                                    if show_Xpc == 1:
                                        hold('on')
                                        plt.plot(tyr_clim,mean(np.squeeze(Xpcf_climstdm(j,:,sels)),2),'Color',colp,'LineStyle','-','LineWidth',lwidth0)
                                    if show_Xpc2 == 1:
                                        hold('on')
                                        plt.plot(tyr_clim,mean(np.squeeze(Xpcf2_climstdm(j,:,sels)),2),'Color',colp2,'LineStyle','-','LineWidth',lwidth0)
                                    set(gca,'XLim',xlim,'FontSize',fontsize)
                                    if fix_ylim == 1:
                                        set(gca,'Ylim',ylim)
                                    if j > 8:
                                        plt.xlabel(xstr,'FontSize',lfontsize)
                                    plt.ylabel(ystr,'FontSize',lfontsize)
                                else:
                                    if (plot_type == 6):
                                        # REDO            sel01 = sel0(year01(sel0)>=yearmin1 & year01(sel0)<=yearmax1 & ~isnan(X01(sel0,1)));
# REDO            selp1 = selp(yearp1(selp)>=yearmin1 & yearp1(selp)<=yearmax1 & ~isnan(Xpc21(selp,1)));
                                        qselhourly = np.array([0.025,0.975])
                                        hourv = unique(hour0fc)
                                        X0f_hourly1 = NaN * np.ones((24,1))
                                        X0f_hourlyL1 = X0f_hourly1
                                        X0f_hourlyH1 = X0f_hourly1
                                        if show_Xpc2 == 1:
                                            Xpcf2_hourly1 = X0f_hourly1
                                            Xpcf2_hourlyL1 = X0f_hourly1
                                            Xpcf2_hourlyH1 = X0f_hourly1
                                        for i in np.arange(1,len(hourv)+1).reshape(-1):
                                            sel011 = sel01(hour0fc(sel01) == hourv(i))
                                            X011 = X01(sel011,:)
                                            X011 = X011
                                            X0f_hourly1[i] = mean(X011)
                                            X0f_hourlyL1[i] = quantile(X011,qselhourly(1))
                                            X0f_hourlyH1[i] = quantile(X011,qselhourly(2))
                                            if (show_Xpc2 == 1):
                                                selp11 = selp1(hourpfc(selp1) == hourv(i))
                                                Xpc211 = Xpc21(selp11,:)
                                                Xpc211 = Xpc211
                                                Xpcf2_hourly1[i] = mean(Xpc211)
                                                Xpcf2_hourlyL1[i] = quantile(Xpc211,qselhourly(1))
                                                Xpcf2_hourlyH1[i] = quantile(Xpc211,qselhourly(2))
                                        lwidth1 = 1
                                        ylim = []
                                        #ylim = [0 100]; #This to focus on low values where smoothing effect of ERA5 averages may be influential.
                                        errorbar(hourv,X0f_hourly1,X0f_hourly1 - X0f_hourlyL1,X0f_hourlyH1 - X0f_hourly1,'kx','LineWidth',lwidth1)
                                        if show_Xpc2 == 1:
                                            hold('on')
                                            errorbar(hourv,Xpcf2_hourly1,Xpcf2_hourly1 - Xpcf2_hourlyL1,Xpcf2_hourlyH1 - Xpcf2_hourly1,'co','LineWidth',lwidth1)
                                        if not len(ylim)==0 :
                                            set(gca,'YLim',ylim)
            if (show_daily == 1 and aggregation_type != 2 and plot_type < 2):
                plt.title(np.array(['yrday = ',str(yrdaysel(j)),tstr1]),'FontSize',tfontsize,'FontWeight',tfontweight)
            else:
                plt.title(np.array([tstrv[j],tstr1]),'FontSize',tfontsize,'FontWeight',tfontweight)
        if show_hourly == 1:
            avstr1 = 'hourly data'
        else:
            if show_daily == 1:
                avstr1 = 'daily averages'
            else:
                avstr1 = 'monthly averages'
        statstr1 = []
        symstr1 = []
        if plot_type == 2 or plot_type == 4:
            statstr1 = ' decadal means'
            symstr1 = ', cross'
        if plot_type == 3 or plot_type == 5:
            statstr1 = ' decadal stds'
            symstr1 = ', cross'
        if (show_hourly == 1 and match_subdelta_deltamean == 1):
            matchstr = ', match_subdelta_deltamean=1'
        else:
            matchstr = []
        if len(sels) == 1:
            locstr = np.array(['at (',num2str(lat1(sels),'%3.2f'),'N,',num2str(lon1(sels),'%3.2f'),'E)'])
        else:
            locstr = 'over ROHO800 region'
        if seasonal == 1:
            tolstr1 = tolstr(np.arange(3,end()+1))
        else:
            tolstr1 = 'non-seasonal'
        if (plot_type == 1):
            if (show_Xpc == 1):
                sgtitle(np.array([varstrc,' from ',avstr1,', ',source0str,' quantiles vs. NORESM2-SSP585 quantiles ',locstr,' using: ',scstr2,' delta-change (magenta); ',scstr2,' QC (cyan) (',tolstr1,matchstr,')']),'FontSize',sgfontsize,'FontWeight','bold','Interpreter',interpstr)
            else:
                sgtitle(np.array([varstrc,' from ',avstr1,', ',source0str,' quantiles vs. NORESM2-SSP585 quantiles ',locstr,'\newlineusing ',scstr2,' quantile correction (',tolstr1,matchstr,'), mean,maxRMSE=',num2str(mean(RMSEv),fmtstr),',',num2str(np.amax(RMSEv),fmtstr),', mean,maxMAE=',num2str(mean(MAEv),fmtstr),',',num2str(np.amax(MAEv),fmtstr),unitstrc]),'FontSize',sgfontsize,'FontWeight','bold','Interpreter',interpstr)
                #' using ',scstr2,' quantile correction (',tolstr(3:end),matchstr,'), mean,maxMAE=',num2str(mean(MAEv),fmtstr),',',num2str(max(MAEv),fmtstr),unitstrc],'FontSize',sgfontsize,'FontWeight','bold','Interpreter',interpstr)
        else:
            if (show_Xpc == 1):
                sgtitle(np.array([varstrc,statstr1,' from ',avstr1,' ',locstr,': ',source0str,' (black',symstr1,'); NORESM2-SSP585 ',scstr2,' delta-change (magenta); NORESM2-SSP585 ',scstr2,' QC (cyan) (',tolstr1,matchstr,')']),'FontSize',sgfontsize,'FontWeight','bold','Interpreter',interpstr)
            else:
                sgtitle(np.array([varstrc,statstr1,' from ',avstr1,' ',locstr,': ',source0str,' (black',symstr1,'); NORESM2-SSP585 ',scstr2,' QC (cyan lines) (',tolstr1,matchstr,')']),'FontSize',sgfontsize,'FontWeight','bold','Interpreter',interpstr)
        clear('X01','Xpc1','Xpc21','Xf','X011','Xpc11','Xpc211')
    ## Plot statistics derived from all months if req'd
    doplot_1 = 0
    if (doplot_1 == 1):
        ##
        nrows = 2
        ncols = 2
        xstr = 'Year'
        lfontsize = 13
        #if strcmp(varstr,'Tair'); ystr = 'T [\circC]'; end
        ystr = np.array([varstr,unitstr])
        tfontsize = 13
        tfontweight = 'bold'
        fontsize = 11
        legfontsize = 10
        col0 = 'k'
        colp = 'm'
        colp2 = 'c'
        msize = 12
        lwidth = 2
        fix_ylim = 2
        show_Xpc = 1
        if doXpc == 0:
            show_Xpc = 0
        show_Xpc2 = 1
        show_Xpc2o2 = 0
        show_hourly = 1
        show_daily = 1
        if use_daily == 0:
            show_hourly = 0
            show_daily = 0
        #sels = 8;
#sels = find(lat1==60.5 & lon1==7);
        sels = np.arange(1,ns+1)
        if ns == 1:
            sels = 1
        if use_daily == 0:
            show_daily = 0
        xfac = 1.2
        yfac = 1.05
        plot_type = 1
        #1 for quantile-quantile plots (daily and hourly)
        aggregation_type = 0
        if plot_type == 1:
            aggregation_type = 2
        #yearmin1 = yearminref; yearmax1 = yearmaxref;
        yearmin1 = 2000
        yearmax1 = 2009
        #if (do_ocean_bcs==1 && test_method==1); yearmin1 = 2008; yearmax1 = 2014; end
        xlim = np.array([2010,2100])
        if len(sels) == 1:
            locstr = np.array(['at (',num2str(lat1(sels),'%3.2f'),'N,',num2str(lon1(sels),'%3.2f'),'E)'])
        else:
            locstr = 'over ROHO800 region'
        plt.figure(10)
        clf
        if (plot_type == 0):
            subplot(nrows,ncols,1)
            plt.plot(tyr_clim(1),mean(X0_meanm(1,sels),2),'+','MarkerSize',msize,'Color',col0,'LineWidth',lwidth)
            if show_Xpc == 1:
                hold('on')
                plt.plot(tyr_clim,mean(Xpc_meanm(:,sels),2),'Color',colp,'LineStyle','-','LineWidth',lwidth)
            if show_Xpc2 == 1:
                hold('on')
                plt.plot(tyr_clim,mean(Xpc2_meanm(:,sels),2),'Color',colp2,'LineStyle','-','LineWidth',lwidth)
            set(gca,'XLim',xlim,'FontSize',fontsize)
            resizesubplot(gca,xfac,yfac,- 0.02,0)
            plt.xlabel(xstr,'FontSize',lfontsize)
            plt.ylabel(ystr,'FontSize',lfontsize)
            plt.title(np.array([varstr,' decadal means from daily averages ',locstr]),'FontSize',tfontsize,'FontWeight',tfontweight,'Interpreter','none')
            if show_Xpc == 1:
                legstrv = np.array(['ERA5 hindcast',np.array(['NORESM2-SSP585 ',scstr2,' DC (',tolstr(np.arange(3,end()+1)),')']),np.array(['NORESM2-SSP585 ',scstr2,' QC (',tolstr(np.arange(3,end()+1)),')'])])
            else:
                legstrv = np.array(['ERA5 hindcast',np.array(['NORESM2-SSP585 ',scstr2,' QC (',tolstr(np.arange(3,end()+1)),')'])])
            plt.legend(legstrv,'FontSize',legfontsize,'Location','SouthEast')
            box('on')
            subplot(nrows,ncols,2)
            plt.plot(tyr_clim(1),mean(X0f_meanm(1,sels),2),'+','MarkerSize',msize,'Color',col0,'LineWidth',lwidth)
            if show_Xpc == 1:
                hold('on')
                plt.plot(tyr_clim,mean(Xpcf_meanm(:,sels),2),'Color',colp,'LineStyle','-','LineWidth',lwidth)
            if show_Xpc2 == 1:
                hold('on')
                plt.plot(tyr_clim,mean(Xpcf2_meanm(:,sels),2),'Color',colp2,'LineStyle','-','LineWidth',lwidth)
            set(gca,'XLim',xlim,'FontSize',fontsize)
            resizesubplot(gca,xfac,yfac,0.02,0)
            plt.xlabel(xstr,'FontSize',lfontsize)
            plt.ylabel(ystr,'FontSize',lfontsize)
            plt.title(np.array([varstr,' decadal means from hourly data ',locstr]),'FontSize',tfontsize,'FontWeight',tfontweight,'Interpreter','none')
            if show_Xpc == 1:
                legstrv = np.array(['ERA5 hindcast',np.array(['NORESM2-SSP585 ',scstr2,' DC (',tolstr(np.arange(3,end()+1)),matchstr,')']),np.array(['NORESM2-SSP585 ',scstr2,' QC (',tolstr(np.arange(3,end()+1)),matchstr,')'])])
            else:
                legstrv = np.array(['ERA5 hindcast',np.array(['NORESM2-SSP585 ',scstr2,' QC (',tolstr(np.arange(3,end()+1)),matchstr,')'])])
            plt.legend(legstrv,'FontSize',legfontsize,'Location','SouthEast','Interpreter','none')
            box('on')
            subplot(nrows,ncols,3)
            plt.plot(tyr_clim(1),mean(X0_stdm(1,sels),2),'+','MarkerSize',msize,'Color',col0,'LineWidth',lwidth)
            if show_Xpc == 1:
                hold('on')
                plt.plot(tyr_clim,mean(Xpc_stdm(:,sels),2),'Color',colp,'LineStyle','-','LineWidth',lwidth)
            if show_Xpc2 == 1:
                hold('on')
                plt.plot(tyr_clim,mean(Xpc2_stdm(:,sels),2),'Color',colp2,'LineStyle','-','LineWidth',lwidth)
            set(gca,'XLim',xlim,'FontSize',fontsize)
            resizesubplot(gca,xfac,yfac,- 0.02,0)
            plt.xlabel(xstr,'FontSize',lfontsize)
            plt.ylabel(ystr,'FontSize',lfontsize)
            plt.title(np.array([varstr,' decadal stds from daily averages ',locstr]),'FontSize',tfontsize,'FontWeight',tfontweight,'Interpreter','none')
            #legend({'ERA5 hindcast',['NORESM2-SSP585 ',scstr2,' DC (',tolstr(3:end),')'],['NORESM2-SSP585 ',scstr2,' QC (',tolstr(3:end),')']},'Location','SouthWest')
            box('on')
            subplot(nrows,ncols,4)
            plt.plot(tyr_clim(1),mean(X0f_stdm(1,sels),2),'+','MarkerSize',msize,'Color',col0,'LineWidth',lwidth)
            if show_Xpc == 1:
                hold('on')
                plt.plot(tyr_clim,mean(Xpcf_stdm(:,sels),2),'Color',colp,'LineStyle','-','LineWidth',lwidth)
            if show_Xpc2 == 1:
                hold('on')
                plt.plot(tyr_clim,mean(Xpcf2_stdm(:,sels),2),'Color',colp2,'LineStyle','-','LineWidth',lwidth)
            set(gca,'XLim',xlim,'FontSize',fontsize)
            resizesubplot(gca,xfac,yfac,0.02,0)
            plt.xlabel(xstr,'FontSize',lfontsize)
            plt.ylabel(ystr,'FontSize',lfontsize)
            plt.title(np.array([varstr,' decadal stds from hourly data ',locstr]),'FontSize',tfontsize,'FontWeight',tfontweight,'Interpreter','none')
            #legend({'ERA5 hindcast',['NORESM2-SSP585 ',scstr2,' DC (',tolstr(3:end),')'],['NORESM2-SSP585 ',scstr2,' QC (',tolstr(3:end),')']},'Location','SouthWest')
            box('on')
        else:
            if (plot_type == 1):
                subplot(1,2,1)
                # REDO        sel01 = find(year0>=yearmin1 & year0<=yearmax1 & ~isnan(X0(:,1)));
#sel01 = sel01(2:end); #This is to avoid missing 01/01/2000 from tdp
#sel01 = sel01(~isnan(Xpc2(sel01-1,1))); #To avoid potential NaNs in Xpc2
#selp1 = sel01-1;
# REDO        selp1 = find(yearp>=yearmin1 & yearp<=yearmax1 & ~isnan(Xpc2(:,1)));
                nsel = np.amin(len(sel01),len(selp1))
                sel01 = sel01(np.arange(1,nsel+1))
                selp1 = selp1(np.arange(1,nsel+1))
                #if length(selp1)<length(sel01); selp1 = [selp1;selp1(1:(length(sel01)-length(selp1)))]; end
                X011 = X0(sel01,sels)
                X011 = __builtint__.sorted(X011)
                Xf = X011
                if show_Xpc == 1:
                    Xpc11 = Xpc(selp1,sels)
                    Xpc11 = __builtint__.sorted(Xpc11)
                    Xf = np.array([[X011],[Xpc11]])
                if show_Xpc2 == 1:
                    Xpc211 = Xpc2(selp1,sels)
                    Xpc211 = __builtint__.sorted(Xpc211)
                    Xf = np.array([[Xf],[Xpc211]])
                if show_Xpc2o2 == 1:
                    Xpc211 = out2.Xpco2(selp1,sels)
                    Xpc211 = __builtint__.sorted(Xpc211)
                    Xf = np.array([[Xf],[Xpc211]])
                xlim = np.array([np.amin(Xf),np.amax(Xf)])
                plt.plot(xlim,xlim,'-','Color',np.array([0.7,0.7,0.7]))
                if show_Xpc == 1:
                    hold('on')
                    plt.plot(Xpc11,X011,'.','Color',colp)
                if show_Xpc2 == 1 or show_Xpc2o2 == 1:
                    hold('on')
                    plt.plot(Xpc211,X011,'.','Color',colp2)
                plt.axis(np.array([xlim,xlim]))
                resizesubplot(gca,xfac,yfac,- 0.02,0)
                plt.ylabel(np.array(['q',varstr,' (ERA5, daily averages)',unitstr]),'FontSize',lfontsize,'Interpreter','none')
                plt.xlabel(np.array(['q',varstr,' (NORESM-ERA5, daily averages)',unitstr]),'FontSize',lfontsize,'Interpreter','none')
                fmtstr = '%3.3f'
                if str(varstr) == str('rain'):
                    fmtstr = '%3.8f'
                if str(varstr) == str('Qair'):
                    fmtstr = '%3.7f'
                MAE = mean(np.abs(Xpc211 - X011))
                RMSE = np.sqrt(mean((Xpc211 - X011) ** 2))
                #title(['MAE = ',num2str(MAE,fmtstr),unitstr,', n = ',str(length(X011))],'FontSize',tfontsize,'FontWeight',tfontweight)
                plt.title(np.array(['RMSE = ',num2str(RMSE,fmtstr),', MAE = ',num2str(MAE,fmtstr),unitstr,', n = ',str(len(X011))]),'FontSize',tfontsize,'FontWeight',tfontweight)
                if (show_hourly == 1):
                    subplot(1,2,2)
                    # REDO            sel01 = find(year0f>=yearmin1 & year0f<=yearmax1 & ~isnan(X0f(:,1)));
#sel01 = sel01(25:end); #This is to avoid missing 01/01/2000 from tdpf.
#td0f(sel01) is now aligned with tdpf (e.g. datestr(td0f(sel01(1:10))) agrees with datestr(tdpf((1:10))))
#sel01 = sel01(~isnan(Xpcf2(sel01-24,1))); #To avoid potential NaNs in Xpcf2
#selp1 = sel01-24;
# REDO            selp1 = find(yearpf>=yearmin1 & yearpf<=yearmax1 & ~isnan(Xpcf2(:,1)));
                    nself = np.amin(len(sel01),len(selp1))
                    sel01 = sel01(np.arange(1,nself+1))
                    selp1 = selp1(np.arange(1,nself+1))
                    #if length(selp1)>length(sel01); selp1 = selp1(1:length(sel01)); end #selp1 can end up 1 larger due to rounding error, datestr(max(tdpf(selp1))) = '31-Dec-2019 23:59:59'
#if length(selp1)<length(sel01); selp1 = [selp1;selp1(1:(length(sel01)-length(selp1)))]; end
                    X011 = X0f(sel01,sels)
                    X011 = __builtint__.sorted(X011)
                    Xf = X011
                    if show_Xpc == 1:
                        Xpc11 = Xpcf(selp1,sels)
                        Xpc11 = __builtint__.sorted(Xpc11)
                        Xf = np.array([[X011],[Xpc11]])
                    if show_Xpc2 == 1:
                        Xpc211 = Xpcf2(selp1,sels)
                        Xpc211 = __builtint__.sorted(Xpc211)
                        Xf = np.array([[Xf],[Xpc211]])
                    xlim = np.array([np.amin(Xf),np.amax(Xf)])
                    plt.plot(xlim,xlim,'-','Color',np.array([0.7,0.7,0.7]))
                    if show_Xpc == 1:
                        hold('on')
                        plt.plot(Xpc11,X011,'.','Color',colp)
                    if show_Xpc2 == 1:
                        hold('on')
                        plt.plot(Xpc211,X011,'.','Color',colp2)
                    plt.axis(np.array([xlim,xlim]))
                    resizesubplot(gca,xfac,yfac,0.02,0)
                    plt.ylabel(np.array(['q',varstr,' (ERA5, hourly)',unitstr]),'FontSize',lfontsize,'Interpreter','none')
                    plt.xlabel(np.array(['q',varstr,' (NORESM-ERA5, hourly)',unitstr]),'FontSize',lfontsize,'Interpreter','none')
                    MAEf = mean(np.abs(Xpc211 - X011))
                    RMSEf = np.sqrt(mean((Xpc211 - X011) ** 2))
                    #title(['MAE = ',num2str(MAEf,fmtstr),unitstr,', n = ',str(length(X011))],'FontSize',tfontsize,'FontWeight',tfontweight)
                    plt.title(np.array(['RMSE = ',num2str(RMSEf,fmtstr),', MAE = ',num2str(MAEf,fmtstr),unitstr,', n = ',str(len(X011))]),'FontSize',tfontsize,'FontWeight',tfontweight)
                clear('Xf','X011','Xpc11','Xpc211')
    ## Write to NetCDF if req'd
    write_nc = 1
    if (write_nc == 1):
        ##
        verstr = '_v3'
        if (do_ocean_bcs == 1):
            # REDO        sel1 = find(yearp>=2000 & yearp<=2099);
            if np.any(str(varstr) == str(np.array(['temp','salt','zeta','ubar','vbar','u','v']))):
                modelstr1 = 'NorESM2-GLORYS12reanal'
            else:
                modelstr1 = 'NorESM2-NorESMreanal1'
            fnamep1 = np.array(['roho800_bry_',modelstr1,'_',scenstr1,'_2000-2099',verstr,'.nc'])
            fname0f = fname0
            fname0str = fname01
            varstr1c = varstrc
            use_ncwriteschema = 0
            time_1day = 86400
            tdpf = tdp
            yearpf = yearp
        else:
            sel1 = np.arange(1,npf+1)
            fnamep1 = np.array(['roho800_',varstr,'_NorESM2-ERA5_',scenstr1,'_',str(np.amin(yearpf)),'-',str(np.amax(yearpf)),verstr,'.nc'])
            fname0str = np.array([fname0f1,', ',fname01])
            varstr1c = varstr
            use_ncwriteschema = 1
            time_1day = 1
        np1 = len(sel1)
        fnamep = np.array([dir1,fnamep1])
        print(np.array(['Writing bias-corrected projections to: ',fnamep1]))
        tic
        if (use_ncwriteschema == 1):
            finfo = ncinfo(fname0c)
            finfo.Filename = fnamep
            finfo.Dimensions(str(np.array([finfo.Dimensions.Name])) == str(varstrt)).Length = np1
            finfo.Variables(str(np.array([finfo.Variables.Name])) == str(varstrt)).Dimensions.Length = np1
            dimnames1 = np.array([finfo.Variables(str(np.array([finfo.Variables.Name])) == str(varstr1c)).Dimensions.Name])
            finfo.Variables(str(np.array([finfo.Variables.Name])) == str(varstr1c)).Dimensions(str(dimnames1) == str(varstrt)).Length = np1
            history1 = finfo.Attributes(str(np.array([finfo.Attributes.Name])) == str('history')).Value
            history1 = np.array([date,': Bias-corrected projections calculated using make_roho_biascorrected_forcings.m with method=',str(method),': ',newline,'             ',source0str,' file(s): ',fname0str,';',newline,'             ESM file: ',fname11,'.',newline,history1])
            history1 = np.array([history1,newline,'options used for fn_biascorrect_projections.m:',newline,evalc('disp(opt2)')])
            if isfield(opt2,'opt_subq'):
                history1 = np.array([history1,newline,'where opt_subq = ',newline,evalc('disp(opt_subq)'),newline,'(see opt2c_',varstr,'_ROHO800',verstr,'.mat)'])
            finfo.Attributes(str(np.array([finfo.Attributes.Name])) == str('history')).Value = history1
            ncwriteschema(fnamep,finfo)
        timepf = time_1day * (tdpf(sel1) - datenum(1948,1,1))
        ncwrite(fnamep,varstrt,timepf)
        if (do_ocean_bcs == 1 and str(varstr) == str('zeta') == 1):
            nccopy_vars(np.array(['theta_b','hc','theta_s','Tcline','s_rho','lat_psi','lon_psi','angle','h','lat_rho','lon_rho']),fname0,fnamep,struct('write_atts',0,'newfile',0))
        else:
            if (do_ocean_bcs == 0):
                ncwrite(fnamep,'lat',lat11)
                ncwrite(fnamep,'lon',lon11)
        yearpminv = np.arange(2000,2090+10,10)
        yearpmaxv = np.arange(2009,2099+10,10)
        nchunks = len(yearpminv)
        for i in np.arange(1,nchunks+1).reshape(-1):
            # REDO        self1 = find(yearpf(sel1)>=yearpminv(i) & yearpf(sel1)<=yearpmaxv(i));
            npf1 = len(self1)
            if (do_ocean_bcs == 1):
                X1fc = NaN * np.ones((n0,N,npf1))
                for j in np.arange(1,npf1+1).reshape(-1):
                    X1fc[np.arange[1,n0+1],np.arange[1,N+1],j] = reshape(Xpc2(sel1(self1(j)),:),n0,N)
            else:
                X1fc = NaN * np.ones((nlonsel,nlatsel,npf1))
                for j in np.arange(1,npf1+1).reshape(-1):
                    X1fc[np.arange[1,nlonsel+1],np.arange[1,nlatsel+1],j] = np.transpose(reshape(Xpcf2(sel1(self1(j)),:),nlatsel,nlonsel))
            if (do_ocean_bcs == 1 and N == 1):
                ncwrite(fnamep,varstr1c,np.squeeze(X1fc),np.array([1,self1(1)]))
            else:
                ncwrite(fnamep,varstr1c,X1fc,np.array([1,1,self1(1)]))
            print(np.array(['Done writing ',varstr,' bias-corrected projections for years ',str(yearpminv(i)),' to ',str(yearpmaxv(i))]))
        print(np.array(['Done writing bias-corrected projections to ',fnamep1]))
        toc
        clear('X1fc')

## Save figures
savefigure = 0

if (savefigure == 1):
    ##
#ifig = 3;
    #Remember to enlarge figure to full size BEFORE you run this section
    if oceanmod == 1:
        dir00 = 'E:/Users/Phil/'
    else:
        dir00 = 'C:/Data/'
    figdir1 = np.array([dir00,'figs_temp/'])
    figdir2 = np.array([dir00,'figs_temp/fig/'])
    #appstr = ['_tol_dtclim',str(tol_dtclim),'_bymonth1_ksmootndclosest200_lrorder10'];
#appstr = ['_tol_dtclim',str(tol_dtclim),'_bymonth0_ksmoothndclosest50'];
#appstr = ['_tol_dtclim',str(tol_dtclim),'_bymonth1_ksmoothndclosest200_lrorder00'];
#appstr = ['_tol_dtclim',str(tol_dtclim),'_bymonth1_nparb3_ksmooth200lrorderres0'];
#appstr = ['_tol_dtclim',str(tol_dtclim),'_bymonth1_nparb3_gaussianres'];
#appstr = ['_tol_dtclim',str(tol_dtclim),'_match_subdelta_deltamean_frcrit_capSWI0_use_SWI_hav1'];
#appstr = ['_tol_dtclim',str(tol_dtclim),'_match_subdelta_deltamean_hourly'];
#appstr = ['_tol_dtclim',str(tol_dtclim),'_bymonth1_powerlaw_ksmooth200lrorderres0'];
#appstr = ['_tol_dtclim',str(tol_dtclim),'_bymonth1_ksmooth200lrorder1_gaussian'];
    methodstr = np.array(['method',str(method)])
    if legacy == 1:
        methodstr = np.array([methodstr,'legacy'])
    if seasonal == 1:
        tolstr1 = np.array(['_tol_dtclim',str(tol_dtclim)])
    else:
        tolstr1 = '_nonseasonal'
    appstr = np.array(['_',methodstr,tolstr1])
    if not ('match_spatially' is not None) :
        match_spatially = 0
    if match_spatially == 1:
        mstr = '_matched_spatially'
    else:
        mstr = []
    #ifig = 1; figstr = ['q',varstr,'_daily_allROHO800_era5_NORESM2_ssp585_daily_QC',appstr];
#ifig = 1; figstr = ['q',varstr,'_hourly_allROHO800_era5_NORESM2_ssp585_daily_QC',mstr,appstr];
#ifig = 1; figstr = [varstr,'_daily_decmean_allROHO800_era5_NORESM2_ssp585_daily_QC',appstr];
#ifig = 1; figstr = [varstr,'_hourly_decmean_allROHO800_era5_NORESM2_ssp585_daily_QC',appstr];
#ifig = 1; figstr = [varstr,'_daily_decstd_allROHO800_era5_NORESM2_ssp585_daily_QC',appstr];
#ifig = 1; figstr = [varstr,'_hourly_decstd_allROHO800_era5_NORESM2_ssp585_daily_QC',appstr];
#ifig = 1; figstr = [varstr,'_diel_variability_allROHO800_era5_NORESM2_ssp585_daily_QC',appstr];
#ifig = 3; figstr = [varstr,'_daily_clim_allROHO800_era5_NORESM2_ssp585_daily_QC',appstr];
#ifig = 10; figstr = [varstr,'_daily_hourly_decmeanstd_allmonths_allROHO800_era5_NORESM2_ssp585_daily_QC',appstr];
#ifig = 10; figstr = ['q',varstr,'_daily_hourly_allmonths_allROHO800_era5_NORESM2_ssp585_daily_QC',appstr];
#ifig = 30; figstr = [varstr,'_hourly_60p5N_7E_era5_NORESM2_ssp585_daily_QC',appstr];
#ifig = 40; figstr = [varstr,'_',scstr2,'_averages_60p5N_7E_era5_NORESM2_ssp585_daily_QC',appstr];
#ifig = 40; figstr = [varstr,'_',scstr2,'_averages_58p4N_5p7E_3m_GLORYS12reanal_NORESM2_ssp585_monthly_QC',appstr];
#ifig = 40; figstr = [varstr,'_',scstr2,'_averages_58p4N_5p7E_100m_GLORYS12reanal_NORESM2_ssp585_monthly_QC',appstr];
#ifig = 40; figstr = [varstr,'_',scstr2,'_averages_58p4N_5p7E_100m_NORESMreanal1_NORESM2_ssp585_monthly_QC',appstr];
    ifig = 40
    figstr = np.array([varstr,'_',scstr2,'_averages_58p4N_5p7E_3m_NORESMreanal1_NORESM2_ssp585_monthly_QC',appstr])
    optfig = struct('figdir1',figdir1,'figdir2',figdir2,'savepng',1,'saveeps',0,'savefig',0)
    fn_savefig(ifig,figstr,optfig)

##Can save selected time series as a structure variable, e.g.:
#u10 = struct('td0',td0,'X0',X0(:,sels),'td0f',td0f,'X0f',X0f(:,sels),'tdp',tdp,'Xp',Xp(:,sels),'Xpc',Xpc(:,sels),'tdpf',tdpf,'Xpcf',Xpcf(:,sels),'Xpc2',Xpc2(:,sels),'Xpcf2',Xpcf2(:,sels));
#save('u10.mat','u10','-v7.3')

savemat = 0
if (savemat == 1):
    ##
    opt2c = rmfield(opt2,'Xhf')
    opt2c = rmfield(opt2c,'tdhf')
    str1 = np.array(['../ROHO800_mat/ROHO800_projections_mat/opt2c_',varstr,'_ROHO800',verstr,'.mat'])
    save(str1,'opt2c')
    print(np.array(['Saved options structure opt2c as: ',str1]))

##This highlights a potential problem with larger changes going from 2300->0000 each day, due to use of match_subdelta_deltamean=1
#diff1 = mean(abs(diff(Xpcf2(:,sels))),2);figure;plot(hourpf(1:end-1),diff1,'k.');
#for i=1:24; meandiff1(i) = meanan(diff1(hourpf(1:end-1)==(i-1)));end #MADs by hour
#NOTE: Even with match_subdelta_deltamean=0, the hour crossing to new day tends to have a larger
#      absolute difference (e.g. ~2x larger for rain).