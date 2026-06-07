



* REDH/ tmax from CPC
use /media/tack/Climate-7/CPC_daily/final_seasonXregion/final_CPC_RED_regridded_monthXregion_base1989_round.001_clean.dta, clear

sort id year month

* mean SAT
merge 1:1 id year month using /media/tack/Climate-7/GHCN_CAMS/final_seasonXregion/final_GHCN_tas_regridded_monthXregion_base1989_round.001.dta

keep if _merge == 3
drop _merge

* soil moisture
merge 1:1 id year month using /media/tack/Climate-7/VSW1/final_seasonXregion/final_ERA5_swvl1_regridded_monthXregion_base1989_round.001.dta

keep if _merge == 3
drop _merge


* water vapor
merge 1:1 id year month using /media/tack/Climate-7/TCWV/final_seasonXregion/final_ERA5_tcwv_regridded_monthXregion_base1989_round.001.dta

keep if _merge == 3
drop _merge

sort id year month

gen model  = "OBS"
ren swvl1 mrsos
replace mrsos = mrsos*(20/13)
ren redh_day redh_obs
ren redl_day redl_obs

sort id year month
order id year month

save /media/tack/Climate-7/CPC_daily/final_seasonXregion_obs_all/final_obs_monthXregion.dta, replace






** step 2.  SVAR regression **
local sp   = "CPC"
local md = "OBS"
local nloc = 69
local red = "redh_obs"
local yy1 = 1989
local yy2 = 2023

cd  /media/tack/Climate-7/CPC_daily/final_seasonXregion_obs_`red'_`yy1'_`yy2'/

forvalues rr = 1/`nloc' {

  di " -----------   region  `rr'  ---------- "
  qui {
  use   /media/tack/Climate-7/CPC_daily/final_seasonXregion_obs_all/final_obs_monthXregion.dta, clear
  
  * since RED in 1989 are always 1.
  drop if year == 1989
  keep if year >= `yy1' & year <= `yy2'
 
  keep if model ==  "`md'" & id == `rr'
  sort year month
  gen time = _n
  save tmp_`sp'_`md'_`rr'_`red'_`yy1'_`yy2'.dta, replace
 
  set seed 12345
  tsset time
  *varsoc `red' mrsos tas tcwv, maxlag(14)
  *matrix Z = r(stats)
  *svmat Z, name(col)
  *egen minh = min(AIC)
  *gen optimal_lag = lag if minh == AIC
  *sort optimal_lag
  *local optlag = optimal_lag[1]
  *sort time

  local optlag = 10
  var `red' mrsos tas tcwv, lags(1/`optlag') dfk small
 
  regsave using  /media/tack/Climate-7/CPC_daily/final_seasonXregion_obs_`red'_`yy1'_`yy2'/regsave/regsave_`sp'_`md'_region`rr'_`red'_`yy1'_`yy2'.dta, tstat ci detail(all)  replace

  irf create var`rr', step(`optlag') set(tmp_irf_`sp'_`md'_region`rr'_order1_`red'_`yy1'_`yy2') order(mrsos tas tcwv `red') replace
  irf create var`rr', step(`optlag') set(tmp_irf_`sp'_`md'_region`rr'_order2_`red'_`yy1'_`yy2') order(mrsos tcwv tas `red') replace
  irf create var`rr', step(`optlag') set(tmp_irf_`sp'_`md'_region`rr'_order3_`red'_`yy1'_`yy2') order(tas mrsos tcwv `red') replace
  irf create var`rr', step(`optlag') set(tmp_irf_`sp'_`md'_region`rr'_order4_`red'_`yy1'_`yy2') order(tas tcwv mrsos `red') replace
  irf create var`rr', step(`optlag') set(tmp_irf_`sp'_`md'_region`rr'_order5_`red'_`yy1'_`yy2') order(tcwv tas mrsos `red') replace
  irf create var`rr', step(`optlag') set(tmp_irf_`sp'_`md'_region`rr'_order6_`red'_`yy1'_`yy2') order(tcwv mrsos tas `red') replace
 
  irf create var`rr', step(`optlag') set(tmp_irf_`sp'_`md'_region`rr'_order7_`red'_`yy1'_`yy2')   order(tas tcwv `red' mrsos) replace
  irf create var`rr', step(`optlag') set(tmp_irf_`sp'_`md'_region`rr'_order8_`red'_`yy1'_`yy2')   order(tas `red' tcwv mrsos) replace
  irf create var`rr', step(`optlag') set(tmp_irf_`sp'_`md'_region`rr'_order9_`red'_`yy1'_`yy2')   order(tcwv tas `red' mrsos) replace
  irf create var`rr', step(`optlag') set(tmp_irf_`sp'_`md'_region`rr'_order10_`red'_`yy1'_`yy2') order(tcwv `red' tas mrsos) replace
  irf create var`rr', step(`optlag') set(tmp_irf_`sp'_`md'_region`rr'_order11_`red'_`yy1'_`yy2') order(`red' tcwv tas mrsos) replace
  irf create var`rr', step(`optlag') set(tmp_irf_`sp'_`md'_region`rr'_order12_`red'_`yy1'_`yy2') order(`red' tas tcwv mrsos) replace
 
  irf create var`rr', step(`optlag') set(tmp_irf_`sp'_`md'_region`rr'_order13_`red'_`yy1'_`yy2') order(mrsos tcwv `red' tas) replace
  irf create var`rr', step(`optlag') set(tmp_irf_`sp'_`md'_region`rr'_order14_`red'_`yy1'_`yy2') order(mrsos `red' tcwv tas) replace
  irf create var`rr', step(`optlag') set(tmp_irf_`sp'_`md'_region`rr'_order15_`red'_`yy1'_`yy2') order(tcwv mrsos `red' tas) replace
  irf create var`rr', step(`optlag') set(tmp_irf_`sp'_`md'_region`rr'_order16_`red'_`yy1'_`yy2') order(tcwv `red' mrsos tas) replace
  irf create var`rr', step(`optlag') set(tmp_irf_`sp'_`md'_region`rr'_order17_`red'_`yy1'_`yy2') order(`red' tcwv mrsos tas) replace
  irf create var`rr', step(`optlag') set(tmp_irf_`sp'_`md'_region`rr'_order18_`red'_`yy1'_`yy2') order(`red' mrsos tcwv tas) replace
 
  irf create var`rr', step(`optlag') set(tmp_irf_`sp'_`md'_region`rr'_order19_`red'_`yy1'_`yy2') order(mrsos tas `red' tcwv) replace
  irf create var`rr', step(`optlag') set(tmp_irf_`sp'_`md'_region`rr'_order20_`red'_`yy1'_`yy2') order(mrsos `red' tas tcwv) replace
  irf create var`rr', step(`optlag') set(tmp_irf_`sp'_`md'_region`rr'_order21_`red'_`yy1'_`yy2') order(tas mrsos `red' tcwv) replace
  irf create var`rr', step(`optlag') set(tmp_irf_`sp'_`md'_region`rr'_order22_`red'_`yy1'_`yy2') order(tas `red' mrsos tcwv) replace
  irf create var`rr', step(`optlag') set(tmp_irf_`sp'_`md'_region`rr'_order23_`red'_`yy1'_`yy2') order(`red' tas mrsos tcwv) replace
  irf create var`rr', step(`optlag') set(tmp_irf_`sp'_`md'_region`rr'_order24_`red'_`yy1'_`yy2') order(`red' mrsos tas tcwv) replace
  }
 
  di "      1. mrsos -> `red'"
  local list1 = "2 3 4 5 6 13 14 15 19 20 21"
  use tmp_irf_`sp'_`md'_region`rr'_order1_`red'_`yy1'_`yy2'.irf, clear
  gen order = 1
  foreach k in `list1' {
    append using tmp_irf_`sp'_`md'_region`rr'_order`k'_`red'_`yy1'_`yy2'.irf
    replace order = `k' if order == .
  }
  keep if impulse == "mrsos" & response == "`red'"
  gen model = "`md'"
  gen ssp =   "`sp'"
  gen region = `rr'
  save irf_`sp'_`md'_region`rr'_mrsos_to_`red'_`red'_`yy1'_`yy2'.dta, replace
 
  di "      2. tas -> `red'"
  local list1 = "2 3 4 5 6 7 8 9 19 21 22"
  use tmp_irf_`sp'_`md'_region`rr'_order1_`red'_`yy1'_`yy2'.irf, clear
  gen order = 1
  foreach k in `list1' {
    append using tmp_irf_`sp'_`md'_region`rr'_order`k'_`red'_`yy1'_`yy2'.irf
    replace order = `k' if order == .
  }
  keep if impulse == "tas" & response ==  "`red'"
  gen model = "`md'"
  gen ssp =   "`sp'"
  gen region = `rr'
  save irf_`sp'_`md'_region`rr'_tas_to_`red'_`red'_`yy1'_`yy2'.dta, replace
 
  di "      3. tcwv -> `red'"
  local list1 = "2 3 4 5 6 7 9 10 13 15 16"
  use tmp_irf_`sp'_`md'_region`rr'_order1_`red'_`yy1'_`yy2'.irf, clear
  gen order = 1
  foreach k in `list1' {
    append using tmp_irf_`sp'_`md'_region`rr'_order`k'_`red'_`yy1'_`yy2'.irf
    replace order = `k' if order == .
  }
  keep if impulse == "tcwv" & response == "`red'"
  gen model = "`md'"
  gen ssp =   "`sp'"
  gen region = `rr'
  save irf_`sp'_`md'_region`rr'_tcwv_to_`red'_`red'_`yy1'_`yy2'.dta, replace
 
  di "      4. mrsos -> tas"
  local list1 = "2 6 13 14 15 16 17 18 19 20 24"
  use tmp_irf_`sp'_`md'_region`rr'_order1_`red'_`yy1'_`yy2'.irf, clear
  gen order = 1
  foreach k in `list1' {
    append using tmp_irf_`sp'_`md'_region`rr'_order`k'_`red'_`yy1'_`yy2'.irf
    replace order = `k' if order == .
  }
  keep if impulse == "mrsos" & response == "tas"
  gen model = "`md'"
  gen ssp =   "`sp'"
  gen region = `rr'
  save irf_`sp'_`md'_region`rr'_mrsos_to_tas_`red'_`yy1'_`yy2'.dta, replace
 
  di "      5. tcwv -> tas"
  local list1 = "5 6 9 10 11 13 14 15 16 17 18"
  use tmp_irf_`sp'_`md'_region`rr'_order2_`red'_`yy1'_`yy2'.irf, clear
  gen order = 2
  foreach k in `list1' {
    append using tmp_irf_`sp'_`md'_region`rr'_order`k'_`red'_`yy1'_`yy2'.irf
    replace order = `k' if order == .
  }
  keep if impulse == "tcwv" & response == "tas"
  gen model = "`md'"
  gen ssp =   "`sp'"
  gen region = `rr'
  save irf_`sp'_`md'_region`rr'_tcwv_to_tas_`red'_`yy1'_`yy2'.dta, replace
 
  di "      6. `red' -> tas"
  local list1 = "11 12 13 14 15 16 17 18 20 23 24"
  use tmp_irf_`sp'_`md'_region`rr'_order10_`red'_`yy1'_`yy2'.irf, clear
  gen order = 10
  foreach k in `list1' {
    append using tmp_irf_`sp'_`md'_region`rr'_order`k'_`red'_`yy1'_`yy2'.irf
    replace order = `k' if order == .
  }
  keep if impulse == "`red'" & response == "tas"
  gen model = "`md'"
  gen ssp =   "`sp'"
  gen region = `rr'
  save irf_`sp'_`md'_region`rr'_`red'_to_tas_`red'_`yy1'_`yy2'.dta, replace
 
  di "      7. `red' -> mrsos"
  local list1 = "8 9 10 11 12 16 17 18 22 23 24"
  use tmp_irf_`sp'_`md'_region`rr'_order7_`red'_`yy1'_`yy2'.irf, clear
  gen order = 7
  foreach k in `list1' {
    append using tmp_irf_`sp'_`md'_region`rr'_order`k'_`red'_`yy1'_`yy2'.irf
    replace order = `k' if order == .
  }
  keep if impulse == "`red'" & response == "mrsos"
  gen model = "`md'"
  gen ssp =   "`sp'"
  gen region = `rr'
  save irf_`sp'_`md'_region`rr'_`red'_to_mrsos_`red'_`yy1'_`yy2'.dta, replace
 
  di "      8. tas -> mrsos"
  local list1 = "4 5 7 8 9 10 11 12  21 22 23"
  use tmp_irf_`sp'_`md'_region`rr'_order3_`red'_`yy1'_`yy2'.irf, clear
  gen order = 3
  foreach k in `list1' {
    append using tmp_irf_`sp'_`md'_region`rr'_order`k'_`red'_`yy1'_`yy2'.irf
    replace order = `k' if order == .
  }
  keep if impulse == "tas" & response == "mrsos"
  gen model = "`md'"
  gen ssp =   "`sp'"
  gen region = `rr'
  save irf_`sp'_`md'_region`rr'_tas_to_mrsos_`red'_`yy1'_`yy2'.dta, replace
 
  di "      9. tcwv ->  mrsos"
  local list1 = "5 6 7 8 9 10 11 12 15 16 17"
  use tmp_irf_`sp'_`md'_region`rr'_order4_`red'_`yy1'_`yy2'.irf, clear
  gen order = 4
  foreach k in `list1' {
    append using tmp_irf_`sp'_`md'_region`rr'_order`k'_`red'_`yy1'_`yy2'.irf
    replace order = `k' if order == .
  }
  keep if impulse == "tcwv" & response == "mrsos"
  gen model = "`md'"
  gen ssp =   "`sp'"
  gen region = `rr'
  save irf_`sp'_`md'_region`rr'_tcwv_to_mrsos_`red'_`yy1'_`yy2'.dta, replace
 
  di "      10. mrsos  -> tcwv"
  local list1 = "2 3 13 14 18 19 20 21 22 23 24"
  use tmp_irf_`sp'_`md'_region`rr'_order1_`red'_`yy1'_`yy2'.irf, clear
  gen order = 1
  foreach k in `list1' {
    append using tmp_irf_`sp'_`md'_region`rr'_order`k'_`red'_`yy1'_`yy2'.irf
    replace order = `k' if order == .
  }
  keep if impulse == "mrsos" & response == "tcwv"
  gen model = "`md'"
  gen ssp =   "`sp'"
  gen region = `rr'
  save irf_`sp'_`md'_region`rr'_mrsos_to_tcwv_`red'_`yy1'_`yy2'.dta, replace
 
  di "      11. tas -> tcwv"
  local list1 = "3 4 7 8 12 19 20 21 22 23 24"
  use tmp_irf_`sp'_`md'_region`rr'_order1_`red'_`yy1'_`yy2'.irf, clear
  gen order = 1
  foreach k in `list1' {
    append using tmp_irf_`sp'_`md'_region`rr'_order`k'_`red'_`yy1'_`yy2'.irf
    replace order = `k' if order == .
  }
  keep if impulse == "tas" & response == "tcwv"
  gen model = "`md'"
  gen ssp =   "`sp'"
  gen region = `rr'
  save irf_`sp'_`md'_region`rr'_tas_to_tcwv_`red'_`yy1'_`yy2'.dta, replace
 
  di "      12. `red' -> tcwv"
  local list1 = "11 12 14 17 18 19 20 21 22 23 24"
  use tmp_irf_`sp'_`md'_region`rr'_order8_`red'_`yy1'_`yy2'.irf, clear
  gen order = 8
  foreach k in `list1' {
    append using tmp_irf_`sp'_`md'_region`rr'_order`k'_`red'_`yy1'_`yy2'.irf
    replace order = `k' if order == .
  }
  keep if impulse == "`red'" & response == "tcwv"
  gen model = "`md'"
  gen ssp =   "`sp'"
  gen region = `rr'
  save irf_`sp'_`md'_region`rr'_`red'_to_tcwv_`red'_`yy1'_`yy2'.dta, replace
}



*** combine all over region X RED category

local sp   = "CPC"
local md = "OBS"
local nloc = 69
local rdlist "redh_obs redl_obs"
local yylist "1989_2023"

cd  /media/tack/Climate-7/CPC_daily/final_seasonXregion_obs_all/

clear
gen red = "" 
gen years = ""

save final_irf_obs.dta, replace


forvalues rr = 1/`nloc' {
   di  " ---------  region `rr' ---------- "
   qui {
   foreach rd in `rdlist' {
    local vvlist "tas_to_`rd'_`rd' mrsos_to_`rd'_`rd' tcwv_to_`rd'_`rd'  tas_to_mrsos_`rd' `rd'_to_mrsos_`rd' tcwv_to_mrsos_`rd' mrsos_to_tas_`rd' `rd'_to_tas_`rd' tcwv_to_tas_`rd' mrsos_to_tcwv_`rd' `rd'_to_tcwv_`rd' tas_to_tcwv_`rd'"
     foreach vv in `vvlist' {
       di "`vv'"
      foreach yy in `yylist' {
       di "`yy'"  
     use /media/tack/Climate-7/CPC_daily/final_seasonXregion_obs_`rd'_`yy'/irf_`sp'_`md'_region`rr'_`vv'_`yy'.dta, clear
     keep region order  impulse response step irf oirf fevd sfevd stdirf stdoirf stdfevd stdsfevd mse 
     gen red  = "`rd'"    
     gen years = "`yy'"  
     
     append using final_irf_obs.dta
     save                final_irf_obs.dta, replace
}
}
}
}
}



*** add information about latitude to compute weighted average ***

cd  /media/tack/Climate-7/CPC_daily/final_seasonXregion_obs_all/

use final_irf_obs.dta, clear

replace red =  "redh_day" if red == "redh_obs" 
replace red =  "redl_day" if red == "redl_obs" 
replace impulse =  "redh_day" if impulse == "redh_obs" 
replace impulse =  "redl_day"  if impulse == "redl_obs" 
replace response =  "redh_day" if response == "redh_obs" 
replace response =  "redl_day"  if response == "redl_obs" 
gen ssp = "OBS"
gen model = "OBS"

save final_irf_obs.dta, replace

sort region

merge n:1 region using  /media/tack/8GBlack/CMIP6_fromERA5/flexible_coarse_grid_mapping_merged_short60plus.dta

keep if _merge == 3
drop _merge

save final_irf_obs_latlon.dta, replace





* 1. compute global weighted-average of IRF for IR pairs.

cd  /media/tack/Climate-7/CPC_daily/final_seasonXregion_obs_all/

use final_irf_obs_latlon.dta, clear

gen wght = cos(coarse_center_lat*_pi/180)

sort red years impulse response step order region
bysort red years impulse response step order (region): egen t_wght = total(wght)
bysort red years impulse response step order (region): gen f_wght = (wght/t_wght)
*bysort  red years impulse response step order (region): egen ss_wght = total(f_wght)
*sum ss_wght
bysort red years impulse response step order (region):  gen nn = _N
tab nn

bysort red years impulse response step order (region):  gen f_oirf        = f_wght*oirf
bysort red years impulse response step order (region):  gen f_stdoirf   = f_wght*stdoirf
bysort red years impulse response step order (region):  gen f_fevd      = f_wght*fevd
bysort red years impulse response step order (region):  gen f_stdfevd = f_wght*stdfevd

bysort red years impulse response step order (region): egen s_oirf        = total(f_oirf)
bysort red years impulse response step order (region): egen s_stdoirf   = total(f_stdoirf)
bysort red years impulse response step order (region): egen s_fevd      = total(f_fevd)
bysort red years impulse response step order (region): egen s_stdfevd = total(f_stdfevd)
bysort red years impulse response step order (region): keep if _n == 1

* mean over model X order
*bysort red years ssp impulse response step: gen nn = _N
*tab nn

bysort red years impulse response step (order): egen m_oirf    = mean(s_oirf)
bysort red years impulse response step (order): egen sd_oirf   = mean(s_stdoirf)
bysort red years impulse response step (order): egen m_fevd  = mean(s_fevd)
bysort red years impulse response step (order): egen sd_fevd = mean(s_stdfevd)

bysort red years impulse response step (order): keep if _n == 1

keep  red years impulse response step m_oirf sd_oirf m_fevd sd_fevd 
order red years impulse response step m_oirf sd_oirf m_fevd sd_fevd 

save  final_irf_obs_global.dta, replace


* 2. regional mean

cd  /media/tack/Climate-7/CPC_daily/final_seasonXregion_obs_all/

use final_irf_obs_latlon.dta, clear

gen wght = cos(coarse_center_lat*_pi/180)

* breaking by latitude
gen group = 1
replace group = 2 if coarse_center_lat > 30 & coarse_center_lat <=60
replace group = 3 if coarse_center_lat > 0   & coarse_center_lat <=30
replace group = 4 if coarse_center_lat <=0

order ssp red years model group region order impulse response step
sort ssp red years model group region order impulse response step

sort red years ssp impulse response step group model order region
bysort red years ssp impulse response step group model order (region): gen nn = _N
tab nn
drop nn

bysort red years ssp impulse response step group model order (region): egen t_wght = total(wght)
bysort red years ssp impulse response step group model order (region): gen f_wght = (wght/t_wght)

bysort red years ssp impulse response step group model order (region): gen f_oirf        = f_wght*oirf
bysort red years ssp impulse response step group model order (region): gen f_stdoirf   = f_wght*stdoirf
bysort red years ssp impulse response step group model order (region): gen f_fevd      = f_wght*fevd
bysort red years ssp impulse response step group model order (region): gen f_stdfevd = f_wght*stdfevd

bysort red years ssp impulse response step group model order (region): egen s_oirf        = total(f_oirf)
bysort red years ssp impulse response step group model order (region): egen s_stdoirf   = total(f_stdoirf)
bysort red years ssp impulse response step group model order (region): egen s_fevd      = total(f_fevd)
bysort red years ssp impulse response step group model order (region): egen s_stdfevd = total(f_stdfevd)
bysort red years ssp impulse response step group model order (region): keep if _n == 1

bysort red years ssp impulse response step group: gen nn = _N
tab nn
drop nn

bysort red years ssp impulse response step group: egen m_oirf    = mean(s_oirf)
bysort red years ssp impulse response step group: egen sd_oirf   = mean(s_stdoirf)
bysort red years ssp impulse response step group: egen m_fevd  = mean(s_fevd)
bysort red years ssp impulse response step group: egen sd_fevd = mean(s_stdfevd)
bysort red years ssp impulse response step group: keep if _n == 1

bysort red years ssp impulse response step group: gen upper = m_oirf + 1.9*sd_oirf
bysort red years ssp impulse response step group: gen lower  = m_oirf - 1.9*sd_oirf
gen sig = (lower > 0 |  upper < 0)

keep  red years ssp group impulse response step m_oirf sd_oirf m_fevd sd_fevd upper lower 
order red years ssp group impulse response step m_oirf sd_oirf m_fevd sd_fevd upper lower 

save  final_irf_obs_group.dta, replace







** SVAR regression on self effect  **

local sp   = "CPC"
local md = "OBS"
local nloc = 69
local redlist = "redh_obs  redl_obs"
local yy1 = 1989
local yy2 = 2023

cd /media/tack/Climate-7/CPC_daily/final_seasonXregion_obs_all/

clear
gen red = ""
save  final_selfirf_obs.dta, replace

foreach red  in `redlist' {
 forvalues rr = 1/`nloc' {
  di " ------  red   `red'     |  region   `rr'    ------- "
  qui {
  forvalues or = 1/24 {
   use /media/tack/Climate-7/CPC_daily/final_seasonXregion_obs_`red'_`yy1'_`yy2'/tmp_irf_`sp'_`md'_region`rr'_order1_`red'_`yy1'_`yy2'.irf, clear

   order impulse response step oirf stdoirf fevd stdfevd mse
   keep  impulse response step oirf stdoirf fevd stdfevd mse

   keep if impulse == response 

   sort impulse response step

   gen red = "`red'"
   gen region = `rr' 
   gen order  = `or'

   append using final_selfirf_obs.dta
   save                final_selfirf_obs.dta, replace
 }
}
}
}


cd /media/tack/Climate-7/CPC_daily/final_seasonXregion_obs_all/

use  final_selfirf_obs.dta, clear

gen years = "1989_2023"
gen model = "OBS"
gen ssp = "OBS"

sort region

merge n:1 region using  /media/tack/8GBlack/CMIP6_fromERA5/flexible_coarse_grid_mapping_merged_short60plus.dta

keep if _merge == 3
drop _merge

save final_selfirf_obs.dta, replace


*1. make  global mean *

cd /media/tack/Climate-7/CPC_daily/final_seasonXregion_obs_all/

use final_selfirf_obs.dta, clear

gen wght = cos(coarse_center_lat*_pi/180)

duplicates report red years impulse response step order region

sort red years impulse response step order region
bysort red years impulse response step order (region): egen t_wght = total(wght)
bysort red years impulse response step order (region): gen f_wght = (wght/t_wght)
*bysort  red years impulse response step order (region): egen ss_wght = total(f_wght)
*sum ss_wght
bysort red years impulse response step order (region):  gen nn = _N
tab nn
drop nn

bysort red years impulse response step order (region):  gen f_oirf        = f_wght*oirf
bysort red years impulse response step order (region):  gen f_stdoirf   = f_wght*stdoirf
bysort red years impulse response step order (region):  gen f_fevd      = f_wght*fevd
bysort red years impulse response step order (region):  gen f_stdfevd = f_wght*stdfevd

bysort red years impulse response step order (region): egen s_oirf        = total(f_oirf)
bysort red years impulse response step order (region): egen s_stdoirf   = total(f_stdoirf)
bysort red years impulse response step order (region): egen s_fevd      = total(f_fevd)
bysort red years impulse response step order (region): egen s_stdfevd = total(f_stdfevd)
bysort red years impulse response step order (region): keep if _n == 1

* mean over order
bysort red years ssp impulse response step: gen nn = _N
tab nn
drop nn

bysort red years impulse response step (order): egen m_oirf    = mean(s_oirf)
bysort red years impulse response step (order): egen sd_oirf   = mean(s_stdoirf)
bysort red years impulse response step (order): egen m_fevd  = mean(s_fevd)
bysort red years impulse response step (order): egen sd_fevd = mean(s_stdfevd)

bysort red years impulse response step (order): keep if _n == 1

keep  red years ssp impulse response step m_oirf sd_oirf m_fevd sd_fevd
order red years ssp impulse response step m_oirf sd_oirf m_fevd sd_fevd

save  final_selfirf_obs_global.dta, replace


* 2. regional mean: just average over model X order 

cd  /media/tack/Climate-7/CPC_daily/final_seasonXregion_obs_all/

use final_selfirf_obs.dta, clear

duplicates report  red years impulse response step region order 

sort red years impulse response step region order 
bysort red years impulse response step region (order): egen m_oirf    = mean(oirf)
bysort red years impulse response step region (order): egen sd_oirf   = mean(stdoirf)
bysort red years impulse response step region (order): egen m_fevd  = mean(fevd)
bysort red years impulse response step region (order): egen sd_fevd = mean(stdfevd)
bysort red years impulse response step region (order): keep if _n == 1

keep  red years region impulse response step m_oirf sd_oirf m_fevd sd_fevd
order red years region impulse response step m_oirf sd_oirf m_fevd sd_fevd

save  final_selfirf_obs_regional.dta, replace


