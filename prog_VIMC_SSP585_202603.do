** step 1. combine all variables by model

local mdlist  "ACCESS-CM2 CMCC-ESM2 CNRM-CM6-1 GFDL-ESM4 GISS-E2-1-G MIROC6 MPI-ESM1-2-LR"
local vlist1  "mrsos tas tcwv RED"
local vlist2  "mrsos tas tcwv"
local splist  "ssp585"

local nloc  = 69
local rd = 0.001
local by = 1989

foreach sp in `splist' {

cd /media/tack/8GBlack/CMIP6_fromERA5/final_seasonXregion_`sp'/

foreach vv in `vlist1' {
 clear
 gen model = ""
 save tmp1_`vv'.dta, replace

 foreach md in `mdlist' {
  di " ---------  variable   `vv'   |   model   `md' --------------  "
  append using  final_`sp'_`md'_`vv'_regridded_monthXregion_base`by'_round`rd'.dta
  replace model = "`md'" if model == ""
 }
 sort model id year month
 save tmp1_`vv'.dta, replace
}

* merge all variables
use tmp1_RED.dta, clear

foreach vv in `vlist2' {
 di " ---------  merging  `vv'   --------------  "
 merge 1:1 model id year month using tmp1_`vv'.dta
 keep if _merge == 3
 drop _merge
}

sort model id year month

save  final_vimc_data_`sp'.dta, replace

}

tabstat redh_day, s(min p1 p5 p10 p50 p90 p95 p99 max)



** step 2.  SVAR regression for SSP585 **

local sp = "ssp585"
local mdlist  "ACCESS-CM2 CMCC-ESM2 CNRM-CM6-1 GFDL-ESM4 GISS-E2-1-G MIROC6 MPI-ESM1-2-LR"
local nloc = 69
local red = "redl_day"
local yy1 = 2024
local yy2 = 2099

cd  /media/tack/8GBlack/CMIP6_fromERA5/final_seasonXregion_`sp'_`red'_`yy1'_`yy2'/

foreach md in `mdlist' {
 forvalues rr = 1/`nloc' {

  di " ----------- ssp `sp'  &  model  `md'   :   region  `rr'  ---------- "
  qui {
  use  /media/tack/8GBlack/CMIP6_fromERA5/final_seasonXregion_`sp'/final_vimc_data_`sp'.dta, clear
  
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
 
  regsave using  /media/tack/8GBlack/CMIP6_fromERA5/final_seasonXregion_`sp'_`red'_`yy1'_`yy2'/regsave/regsave_`sp'_`md'_region`rr'_`red'_`yy1'_`yy2'.dta, tstat ci detail(all)  replace

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
}



*** collecting all IRF results ***

set more off
clear all

local splist "ssp585"
local vvlist "redl_day redh_day"
local yylist "1989_2023 1989_2099 2024_2099 2040_2069 2070_2099"

foreach sp in `splist' {
 foreach vv in `vvlist' {
  foreach yrange in `yylist' {

    di "--------------- ssp `sp' | variable `vv' | range `yrange' -----------------"

    * Set working directory
    cd "/media/tack/8GBlack/CMIP6_fromERA5/final_seasonXregion_`sp'_`vv'_`yrange'/"

    * Get file list
    local files: dir . files "irf_*.dta"
    local total : word count `files'

    di "Total files: `total'"

    * Skip if no files
    if `total' == 0 {
        di "No files found. Skipping..."
        continue
    }

    * Batch size (tune if needed)
    local batch = 200
    local start = 1
    local batch_id = 1

    * -------------------------------
    * Step 1: Create batch files
    * -------------------------------
    while `start' <= `total' {

        local end = `start' + `batch' - 1
        if `end' > `total' local end = `total'

        di "Processing batch `batch_id': files `start' to `end'"

        * Load first file
        local first : word `start' of `files'
        quietly use "`first'", clear

        * Append rest
        forvalues i = `=`start'+1'/`end' {
            local f : word `i' of `files'
            quietly append using "`f'"
        }

        * Save batch
        save "temp_batch_`batch_id'.dta", replace

        local start = `end' + 1
        local batch_id = `batch_id' + 1
    }

    * -------------------------------
    * Step 2: Combine batch files
    * -------------------------------
    local batch_files: dir . files "temp_batch_*.dta"
    local nb : word count `batch_files'

    di "Combining `nb' batch files..."

    local first : word 1 of `batch_files'
    use "`first'", clear

    forvalues i = 2/`nb' {
        local f : word `i' of `batch_files'
        quietly append using "`f'"
    }

    * -------------------------------
    * Step 3: Clean + organize
    * -------------------------------
    duplicates report ssp model region order impulse response step

    keep  ssp model region order impulse response step irf oirf fevd stdirf stdoirf stdfevd mse

    order ssp model region order impulse response step irf oirf fevd stdirf stdoirf stdfevd mse

    * -------------------------------
    * Step 4: Save output
    * -------------------------------
    *cap mkdir "/media/tack/8GBlack/CMIP6_fromERA5/final_seasonXregion_`sp'_all"

    save "/media/tack/8GBlack/CMIP6_fromERA5/final_seasonXregion_`sp'_all/final_irf_`sp'_`vv'_`yrange'.dta", replace

    * -------------------------------
    * Step 5: Clean temp files
    * -------------------------------
    foreach f of local batch_files {
        erase "`f'"
    }

    di "Finished: `sp' `vv' `yrange'"
  }
 }
}


local sp "ssp585"
local vvlist "redh_day redl_day"
local yylist "1989_2023 1989_2099 2024_2099 2040_2069 2070_2099"

cd /media/tack/8GBlack/CMIP6_fromERA5/final_seasonXregion_`sp'_all/

clear
gen years =  ""
gen red = ""
save final_irf_`sp'_all.dta, replace

foreach vv in `vvlist' {
 foreach yrange in `yylist' {
  
  append using final_irf_`sp'_`vv'_`yrange'.dta
  replace years = "`yrange'"  if years == ""
  replace red = "`vv'" if red == ""
  save final_irf_`sp'_all.dta, replace

 }
}



local sp ="ssp585"

cd /media/tack/8GBlack/CMIP6_fromERA5/final_seasonXregion_`sp'_all/

use  final_irf_`sp'_all.dta, clear

duplicates report ssp red years model region order impulse response step

sort region

merge n:1 region using  /media/tack/8GBlack/CMIP6_fromERA5/flexible_coarse_grid_mapping_merged_short60plus.dta

keep if _merge == 3
drop _merge

save final_irf_`sp'_all.dta, replace



*** 1. make global data ***

local sp ="ssp585"

cd /media/tack/8GBlack/CMIP6_fromERA5/final_seasonXregion_`sp'_all/

use  final_irf_`sp'_all.dta, clear

order ssp red years model region order impulse response step
sort ssp red years model region order impulse response step

gen wght = cos(coarse_center_lat*_pi/180)

sort red years ssp impulse response step model order region
bysort red years ssp impulse response step model order (region): egen t_wght = total(wght)
bysort red years ssp impulse response step model order (region): gen f_wght = (wght/t_wght)
*bysort red years ssp impulse response step model order (region): egen ss_wght = total(f_wght)
*sum ss_wght

bysort red years ssp impulse response step model order (region): gen f_oirf        = f_wght*oirf
bysort red years ssp impulse response step model order (region): gen f_stdoirf   = f_wght*stdoirf
bysort red years ssp impulse response step model order (region): gen f_fevd      = f_wght*fevd
bysort red years ssp impulse response step model order (region): gen f_stdfevd = f_wght*stdfevd

bysort red years ssp impulse response step model order (region): egen s_oirf        = total(f_oirf)
bysort red years ssp impulse response step model order (region): egen s_stdoirf   = total(f_stdoirf)
bysort red years ssp impulse response step model order (region): egen s_fevd      = total(f_fevd)
bysort red years ssp impulse response step model order (region): egen s_stdfevd = total(f_stdfevd)
bysort red years ssp impulse response step model order (region): keep if _n == 1

* mean over model X order
*bysort red years ssp impulse response step: gen nn = _N
*tab nn
bysort red years ssp impulse response step: egen m_oirf    = mean(s_oirf)
bysort red years ssp impulse response step: egen sd_oirf   = mean(s_stdoirf)
bysort red years ssp impulse response step: egen m_fevd  = mean(s_fevd)
bysort red years ssp impulse response step: egen sd_fevd = mean(s_stdfevd)

bysort red years ssp impulse response step: keep if _n == 1

bysort red years ssp impulse response step: gen upper = m_oirf + 1.9*sd_oirf
bysort red years ssp impulse response step: gen lower = m_oirf - 1.9*sd_oirf
gen sig = (lower > 0 |  upper < 0)


keep  red years ssp impulse response step m_oirf sd_oirf m_fevd sd_fevd upper lower 
order red years ssp impulse response step m_oirf sd_oirf m_fevd sd_fevd upper lower 

save final_irf_`sp'_global.dta, replace



*** 2. make group data ***

local sp ="ssp585"

cd /media/tack/8GBlack/CMIP6_fromERA5/final_seasonXregion_`sp'_all/

use  final_irf_`sp'_all.dta, clear

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

save  final_irf_`sp'_group.dta, replace



*** combine all SSPs ***

* global *

cd /media/tack/8GBlack/CMIP6_fromERA5/final_seasonXregion_all/

use                final_irf_ssp126_global.dta, clear
append using final_irf_ssp245_global.dta
append using final_irf_ssp585_global.dta

sort red years ssp impulse response step 

save final_irf_global.dta, replace

* regional *

cd /media/tack/8GBlack/CMIP6_fromERA5/final_seasonXregion_all/

use                 final_irf_ssp126_regional.dta, clear
append using final_irf_ssp245_regional.dta
append using final_irf_ssp585_regional.dta

sort red years ssp region impulse response step 

save final_irf_regional.dta, replace






** SVAR regression on self effect  **

local sp = "ssp585"
local mdlist  "ACCESS-CM2 CMCC-ESM2 CNRM-CM6-1 GFDL-ESM4 GISS-E2-1-G MIROC6 MPI-ESM1-2-LR"
local nloc = 69
local redlist "redl_day redh_day"
local yylist "1989_2023 1989_2099 2024_2099 2040_2069 2070_2099"

cd  /media/tack/8GBlack/CMIP6_fromERA5/final_seasonXregion_`sp'/

clear
gen red = ""
save  final_selfirf_`sp'.dta, replace

foreach red  in `redlist' {
 foreach yy in `yylist' {
 forvalues rr = 1/`nloc' {
  di " ------  red   `red'     |  years   `yy'   |  region   `rr'    ------- "
  qui {
  foreach md in `mdlist' {
   forvalues or = 1/24 {
    use  /media/tack/8GBlack/CMIP6_fromERA5/final_seasonXregion_`sp'_`red'_`yy'/tmp_irf_`sp'_`md'_region`rr'_order`or'_`red'_`yy'.irf, clear

    order impulse response step oirf stdoirf fevd stdfevd mse
    keep  impulse response step oirf stdoirf fevd stdfevd mse

    keep if impulse == response 

    sort impulse response step

    gen ssp = "`sp'"
    gen red = "`red'"
    gen model = "`md'"
    gen years = "`yy'"
    gen region = `rr' 
    gen order  = `or'

    append using final_selfirf_`sp'.dta
    save               final_selfirf_`sp'.dta, replace
 }
}
}
}
}
}

local sp  "ssp585"

cd  /media/tack/8GBlack/CMIP6_fromERA5/final_seasonXregion_`sp'/

use  final_selfirf_`sp'.dta, clear

sort region

merge n:1 region using  /media/tack/8GBlack/CMIP6_fromERA5/flexible_coarse_grid_mapping_merged_short60plus.dta

keep if _merge == 3
drop _merge

save final_selfirf_`sp'.dta, replace


*1. make  global mean *

cd  /media/tack/8GBlack/CMIP6_fromERA5/final_seasonXregion_`sp'/

use final_selfirf_`sp'.dta, clear

gen wght = cos(coarse_center_lat*_pi/180)

sort red years ssp impulse response step model order region
bysort red years ssp impulse response step model order (region): egen t_wght = total(wght)
bysort red years ssp impulse response step model order (region): gen f_wght = (wght/t_wght)
*bysort red years ssp impulse response step model order (region): egen ss_wght = total(f_wght)
*sum ss_wght

bysort red years ssp impulse response step model order (region): gen f_oirf        = f_wght*oirf
bysort red years ssp impulse response step model order (region): gen f_stdoirf   = f_wght*stdoirf
bysort red years ssp impulse response step model order (region): gen f_fevd      = f_wght*fevd
bysort red years ssp impulse response step model order (region): gen f_stdfevd = f_wght*stdfevd

bysort red years ssp impulse response step model order (region): egen s_oirf        = total(f_oirf)
bysort red years ssp impulse response step model order (region): egen s_stdoirf   = total(f_stdoirf)
bysort red years ssp impulse response step model order (region): egen s_fevd      = total(f_fevd)
bysort red years ssp impulse response step model order (region): egen s_stdfevd = total(f_stdfevd)
bysort red years ssp impulse response step model order (region): keep if _n == 1

* mean over model X order (7  x 24 =  168)
bysort red years ssp impulse response step: gen nn = _N
tab nn
drop nn

bysort red years ssp impulse response step: egen m_oirf    = mean(s_oirf)
bysort red years ssp impulse response step: egen sd_oirf   = mean(s_stdoirf)
bysort red years ssp impulse response step: egen m_fevd  = mean(s_fevd)
bysort red years ssp impulse response step: egen sd_fevd = mean(s_stdfevd)

bysort red years ssp impulse response step: keep if _n == 1

bysort red years ssp impulse response step: gen upper = m_oirf + 1.9*sd_oirf
bysort red years ssp impulse response step: gen lower  = m_oirf - 1.9*sd_oirf
gen sig = (lower > 0 |  upper < 0)


keep  red years ssp impulse response step m_oirf sd_oirf m_fevd sd_fevd upper lower 
order red years ssp impulse response step m_oirf sd_oirf m_fevd sd_fevd upper lower 

save  final_selfirf_`sp'_global.dta, replace


*** 2. make group data ***

local sp ="ssp585"

cd  /media/tack/8GBlack/CMIP6_fromERA5/final_seasonXregion_`sp'/

use final_selfirf_`sp'.dta, clear


* breaking by latitude
gen group = 1
replace group = 2 if coarse_center_lat > 30 & coarse_center_lat <=60
replace group = 3 if coarse_center_lat > 0   & coarse_center_lat <=30
replace group = 4 if coarse_center_lat <=0

order ssp red years model group region order impulse response step
sort ssp red years model group region order impulse response step

gen wght = cos(coarse_center_lat*_pi/180)

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

save  final_selfirf_`sp'_group.dta, replace





* nn is model x orders (7 x 12 =  84)
sort red years ssp impulse response step region order 
bysort red years ssp impulse response step region (order): gen nn = _N
tab nn
drop nn

bysort red years ssp impulse response step region (order): egen m_oirf    = mean(oirf)
bysort red years ssp impulse response step region (order): egen sd_oirf   = mean(stdoirf)
bysort red years ssp impulse response step region (order): egen m_fevd  = mean(fevd)
bysort red years ssp impulse response step region (order): egen sd_fevd = mean(stdfevd)
bysort red years ssp impulse response step region (order): keep if _n == 1

bysort red years ssp impulse response step region (order): gen upper = m_oirf + 1.9*sd_oirf
bysort red years ssp impulse response step region (order): gen lower = m_oirf - 1.9*sd_oirf
gen sig = (lower > 0 |  upper < 0)

keep  red years ssp region impulse response step m_oirf sd_oirf m_fevd sd_fevd upper lower 
order red years ssp region impulse response step m_oirf sd_oirf m_fevd sd_fevd upper lower 

save  final_selfirf_`sp'_regional.dta, replace





