
**********************  import data *******************************

******  SSP585 ******

local md = "CMCC-ESM2"
local dt = "ssp585"
local vv = "tmax"

cd /media/tack/8GBlack/CMIP6_fromERA5/original/`md'_`vv'/

foreach y in 2015 2040 2065 2090 {
 
di " -------------- year `y' ------------------ "

insheet using  tasmax_day_`md'_`dt'_`y'_1.txt, clear

duplicates report lat lon year month day

*tw scatter lat lon if year == 2015 & month == 1 & day == 1, msize(vtiny)

ren mean_tasmax `vv'

drop coarse_grid_id

ren group_center_lat coarse_center_lat
ren group_center_lon coarse_center_lon

sort lat lon

* merge coarse grid id *
merge n:1 lat lon using /media/tack/8GBlack/CMIP6_fromERA5/flexible_coarse_grid_mapping_merged.dta
 
keep if _merge == 3
drop _merge

drop if day  > 28
drop if year > 2099 | year < 1979

drop if coarse_center_lat < -60

*tw scatter coarse_center_lat coarse_center_lon if year == 2015 & month == 1 & day == 1, msize(tiny)

drop coarse_grid_id
egen coarse_grid_id = group(coarse_center_lat coarse_center_lon)

sort coarse_grid_id

local nloc = coarse_grid_id[_N]

di `nloc'

save tasmax_day_`md'_`dt'_`y'_1_regridded.dta, replace

qui {
forvalue i = 1/12 {
 forvalues r = 1/`nloc' {
  di " --- month `i' & region `r' ----- "
  use  tasmax_day_`md'_`dt'_`y'_1_regridded.dta, clear
  keep if month == `i' & coarse_grid_id ==`r'
  save tasmax_day_`md'_`dt'_`y'_regridded_month`i'_region`r'.dta, replace
  }
}
}
}




* historical *

local md = "CMCC-ESM2"
local vv = "tmax"
local dt = "historical"

cd /media/tack/8GBlack/CMIP6_fromERA5/original/`md'_`vv'/

foreach y in 1975 2000 {
 
di " -------------- year `y' ------------------ "

insheet using  tasmax_day_`md'_`dt'_`y'_1.txt, clear

duplicates report lat lon year month day

*tw scatter lat lon if year == 2015 & month == 1 & day == 1, msize(vtiny)

ren mean_tasmax `vv'

drop coarse_grid_id

ren group_center_lat coarse_center_lat
ren group_center_lon coarse_center_lon

sort lat lon

* merge coarse grid id *
merge n:1 lat lon using /media/tack/8GBlack/CMIP6_fromERA5/flexible_coarse_grid_mapping_merged.dta
 
keep if _merge == 3
drop _merge

drop if day  > 28
drop if year > 2099 | year < 1979

drop if coarse_center_lat < -60

*tw scatter coarse_center_lat coarse_center_lon if year == 2015 & month == 1 & day == 1, msize(tiny)

drop coarse_grid_id
egen coarse_grid_id = group(coarse_center_lat coarse_center_lon)

sort coarse_grid_id

local nloc = coarse_grid_id[_N]

di `nloc'

save tasmax_day_`md'_`dt'_`y'_1_regridded.dta, replace

qui {
forvalue i = 1/12 {
 forvalues r = 1/`nloc' {
  di " --- month `i' & region `r' ----- "
  use  tasmax_day_`md'_`dt'_`y'_1_regridded.dta, clear
  keep if month == `i' & coarse_grid_id ==`r'
  save tasmax_day_`md'_`dt'_`y'_regridded_month`i'_region`r'.dta, replace
  }
}
}
}








** adding ssp585 + historical to have the whole year **

local md = "CMCC-ESM2"
local vv = "tmax"

cd /media/tack/8GBlack/CMIP6_fromERA5/original/`md'_`vv'/

forvalues m = 1/12 {
 forvalues r = 1/`nloc' {

 di "  --------------  month `m'  & region `r' ----------------- "

use                   tasmax_day_`md'_historical_1975_regridded_month`m'_region`r'.dta, clear
append using  tasmax_day_`md'_historical_2000_regridded_month`m'_region`r'.dta
append using  tasmax_day_`md'_ssp585_2015_regridded_month`m'_region`r'.dta
append using  tasmax_day_`md'_ssp585_2040_regridded_month`m'_region`r'.dta
append using  tasmax_day_`md'_ssp585_2065_regridded_month`m'_region`r'.dta
append using  tasmax_day_`md'_ssp585_2090_regridded_month`m'_region`r'.dta

save          tasmax_day_`md'_all_month`m'_region`r'_regridded.dta, replace

}
}




************  remove outliers :  Tmax   ******************

local md = "CMCC-ESM2"
local vv = "tmax"
local nloc = 69

cd /media/tack/8GBlack/CMIP6_fromERA5/original/`md'_`vv'/

forvalues m = 1/12 {
 forvalues r = 1/`nloc' {

 di " ------------- month `m' & region `r' ---------------- "

 use tasmax_day_`md'_all_month`m'_region`r'_regridded.dta, clear
 
 replace `vv' = `vv' - 273.15

 sort lat lon year month day
 bysort lat lon year: gen c1 = _n
 bysort lat lon year: gen c2 = _N
 bysort lat lon year: gen l1_`vv' = `vv'[_n-1]
 bysort lat lon year: gen l2_`vv' = `vv'[_n-2]
 bysort lat lon year: gen l3_`vv' = `vv'[_n-3]
 bysort lat lon year: gen f1_`vv' = `vv'[_n+1]
 bysort lat lon year: gen f2_`vv' = `vv'[_n+2]
 bysort lat lon year: gen f3_`vv' = `vv'[_n+3]
 bysort lat lon year: gen mneigh     = (l1_`vv' + l2_`vv' + l3_`vv' + f1_`vv' + f2_`vv' + f3_`vv')/6

 bysort lat lon year: replace mneigh = (l1_`vv' + l2_`vv' + l3_`vv')/3 			   if _n == c2
 bysort lat lon year: replace mneigh = (l1_`vv' + l2_`vv' + l3_`vv' + f1_`vv')/4           if _n == c2-1
 bysort lat lon year: replace mneigh = (l1_`vv' + l2_`vv' + l3_`vv' + f1_`vv' + f2_`vv')/5 if _n == c2-2
 bysort lat lon year: replace mneigh = (f1_`vv' + f2_`vv' + f3_`vv')/3 			   if c1 == 1
 bysort lat lon year: replace mneigh = (f1_`vv' + f2_`vv' + f3_`vv' + l1_`vv')/4           if c1 == 2
 bysort lat lon year: replace mneigh = (f1_`vv' + f2_`vv' + f3_`vv' + l1_`vv' +l2_`vv')/5  if c1 == 3

 count if mneigh == .

 replace `vv' = . if (`vv' > mneigh + 25) | (`vv' < mneigh - 25)

 ** number of outliers **
 count if `vv' == .

 keep lat lon year month day tmax coarse_center_lat coarse_center_lon coarse_grid_id

 save tasmax_day_`md'_all_month`m'_region`r'_regridded_clean.dta, replace
}
}



** historical only for ssp 245 & 126 **

local md = "CMCC-ESM2"
local vv  = "tmax"
local nloc = 69

cd /media/tack/8GBlack/CMIP6_fromERA5/original/`md'_`vv'/

forvalues m = 1/12 {

 di "  --------------  month `m'  ----------------- "

 clear
 gen year =.
 save                   /media/tack/8GBlack/CMIP6_fromERA5/original/`md'_`vv'/regridded_new_hist/tasmax_day_`md'_historical_month`m'_regridded_new.dta, replace
 
 qui { 
 forvalues r = 1/`nloc' {

 use                   tasmax_day_`md'_all_month`m'_region`r'_regridded_clean.dta, clear
 keep if year <= 2014

 append using  /media/tack/8GBlack/CMIP6_fromERA5/original/`md'_`vv'/regridded_new_hist/tasmax_day_`md'_historical_month`m'_regridded_new.dta
 save                 /media/tack/8GBlack/CMIP6_fromERA5/original/`md'_`vv'/regridded_new_hist/tasmax_day_`md'_historical_month`m'_regridded_new.dta, replace
 }
}
}






** add shuffling **

local nsf =  100
local nloc =  69
local vv = "tmax"
local md = "CMCC-ESM2"

set seed 1234

cd /media/tack/12TBHelio/CMIP6/`md'/

forvalues m = 7/12 {
 forvalues r = 1/`nloc' {
 
 clear
 gen month = .
 save final_`vv'_`md'_regridded_month`m'_region`r'_addshuffle.dta, replace
 

 di " ------------- month = `m',  region = `r' ----------------  "

 use   tasmax_day_`md'_all_month`m'_region`r'_regridded_clean.dta, clear

 gen shuffle = 0

 append using  final_`vv'_`md'_regridded_month`m'_region`r'_addshuffle.dta
 save          final_`vv'_`md'_regridded_month`m'_region`r'_addshuffle.dta, replace

qui {
* Loop `nsf' times to generate the shuffled datasets
 forvalues s = 1/`nsf' {

    * Reload base data
    use  tasmax_day_`md'_all_month`m'_region`r'_regridded_clean.dta, clear
     
    * Create a random variable to shuffle day within each group
    gen double rand = runiform()
   
    * Sort by group and random variable
    sort lat lon year month rand
    bysort lat lon year month (rand): gen temp_day = _n

    replace day = temp_day

    * Drop temporary variables
    drop rand temp_day
   
    gen shuffle = `s'

    * Save shuffled dataset
    append using  final_`vv'_`md'_regridded_month`m'_region`r'_addshuffle.dta
    save          final_`vv'_`md'_regridded_month`m'_region`r'_addshuffle.dta, replace
}
}
}
}





****************** add record and create 10 year interval ***********************


local vv =  "tmax"
local md = "CMCC-ESM2"
local nloc = 69
local rd = 0.001
local by = 1989

cd /media/tack/12TBHelio/CMIP6/`md'/


* creating base year *

forvalues m = 1/12 {
 forvalues r = 1/`nloc' {
 
 di " ----------  month `m', region `r'  ----------- "
 
 use   tasmax_day_`md'_all_month`m'_region`r'_regridded_clean.dta, clear
 
 gen shuffle = 0

 drop if day > 28

 replace year = `by' if year < `by'

 * creating base year *
 sort lat lon month shuffle day year
 bysort lat lon month shuffle day year: egen m1 = mean(`vv')
 bysort lat lon month shuffle day year: keep if _n == 1
 drop `vv'
 ren m1 `vv'

 save /media/tack/12TBHelio/CMIP6/`md'_`vv'_seasonXregion/final_`vv'_`md'_regridded_month`m'_region`r'_base`by'.dta, replace

}
}



* create 10 year intervaled data *

local vv =  "tmax"
local md = "CMCC-ESM2"
local nloc = 69
local rd = 0.001
local by = 1989

cd /media/tack/8GBlack/CMIP6_fromERA5/original/`md'_`vv'_seasonXregion/

forvalues p = 1/11 {
 forvalues m = 1/12 {
  forvalues r = 1/`nloc' {

di " ----------  part `p',  month `m', region `r'  ----------- "
 
qui {

use  final_`vv'_`md'_regridded_month`m'_region`r'_base`by'.dta, clear

keep if year >= 1989 + (`p'-1)*10 & year <= 1989 + (`p')*10

* adding record indicator *
*gen `vv'_r = round(`vv',`rd')
*replace `vv' = `vv'_r
*drop  `vv'_r

**** further correction ****
gen `vv'_r1 = round(`vv',`rd')
gen `vv'_r2 = round(`vv',0.1)
gen `vv'_r3 = round(`vv',0.2)
gen `vv'_r4 = round(`vv',0.3)

replace `vv' = `vv'_r1 if year <= 2030
replace `vv' = `vv'_r2 if year > 2030 & year <= 2050
replace `vv' = `vv'_r3 if year > 2050 & year <= 2080
replace `vv' = `vv'_r4 if year > 2080

drop `vv'_r*

egen loc = group(lat lon)
egen loc_time = group(lat lon month shuffle day)

sort loc month shuffle day year

sort loc month shuffle day year
bysort loc month shuffle day (year): gen bestt = `vv' if _n == 1
bysort loc month shuffle day (year): replace bestt = max(`vv', bestt[_n-1]) if loc_time == loc_time[_n-1]
bysort loc month shuffle day (year): gen recordmax = 0
bysort loc month shuffle day (year): replace recordmax = 1 if loc_time != loc_time[_n-1] | (bestt > bestt[_n-1])

drop bestt

sort loc month shuffle day year
bysort loc month shuffle day (year): gen bestt = `vv' if _n == 1
bysort loc month shuffle day (year): replace bestt = min(`vv', bestt[_n-1]) if loc_time == loc_time[_n-1]
bysort loc month shuffle day (year): gen recordmin = 0
bysort loc month shuffle day (year): replace recordmin = 1 if loc_time != loc_time[_n-1] | (bestt < bestt[_n-1])

drop bestt loc_time
 
save final_`vv'_`md'_regridded_month`m'_region`r'_base`by'_round`rd'_record_part`p'.dta, replace
}
}
}
}



**************  all combined : create allmonth, season, allmonth X allregion, season x allregion  *****************


local vv =  "tmax"
local md = "CMCC-ESM2"
local nloc = 69
local rd = 0.001
local by = 1989

cd /media/tack/8GBlack/CMIP6_fromERA5/original/`md'_`vv'_seasonXregion/

forvalues p = 1/11 {
  di  " ---- part `p' ------"
  qui {
  * create allmonth data
  forvalues r = 1/`nloc' {
 
    use                  final_`vv'_`md'_regridded_month1_region`r'_base`by'_round`rd'_record_part`p'.dta, clear
    append using final_`vv'_`md'_regridded_month2_region`r'_base`by'_round`rd'_record_part`p'.dta
    append using final_`vv'_`md'_regridded_month12_region`r'_base`by'_round`rd'_record_part`p'.dta
    save                final_`vv'_`md'_regridded_season1_region`r'_base`by'_round`rd'_record_part`p'.dta, replace

   use                  final_`vv'_`md'_regridded_month3_region`r'_base`by'_round`rd'_record_part`p'.dta, clear
   append using final_`vv'_`md'_regridded_month4_region`r'_base`by'_round`rd'_record_part`p'.dta
   append using final_`vv'_`md'_regridded_month5_region`r'_base`by'_round`rd'_record_part`p'.dta
   save                final_`vv'_`md'_regridded_season2_region`r'_base`by'_round`rd'_record_part`p'.dta, replace

   use                  final_`vv'_`md'_regridded_month6_region`r'_base`by'_round`rd'_record_part`p'.dta, clear
   append using final_`vv'_`md'_regridded_month7_region`r'_base`by'_round`rd'_record_part`p'.dta
   append using final_`vv'_`md'_regridded_month8_region`r'_base`by'_round`rd'_record_part`p'.dta
   save                final_`vv'_`md'_regridded_season3_region`r'_base`by'_round`rd'_record_part`p'.dta, replace
   
   use                  final_`vv'_`md'_regridded_month9_region`r'_base`by'_round`rd'_record_part`p'.dta, clear
   append using final_`vv'_`md'_regridded_month10_region`r'_base`by'_round`rd'_record_part`p'.dta
   append using final_`vv'_`md'_regridded_month11_region`r'_base`by'_round`rd'_record_part`p'.dta
   save                final_`vv'_`md'_regridded_season4_region`r'_base`by'_round`rd'_record_part`p'.dta, replace
 
  append using final_`vv'_`md'_regridded_season1_region`r'_base`by'_round`rd'_record_part`p'.dta
  append using final_`vv'_`md'_regridded_season2_region`r'_base`by'_round`rd'_record_part`p'.dta
  append using final_`vv'_`md'_regridded_season3_region`r'_base`by'_round`rd'_record_part`p'.dta
  save                final_`vv'_`md'_regridded_allmonth_region`r'_base`by'_round`rd'_record_part`p'.dta, replace

  }
 
  * create allmonth X allregion data
  use                     final_`vv'_`md'_regridded_allmonth_region1_base`by'_round`rd'_record_part`p'.dta, clear
  forvalues rr = 2/`nloc' {
   append using  final_`vv'_`md'_regridded_allmonth_region`rr'_base`by'_round`rd'_record_part`p'.dta
  }
  save                  final_`vv'_`md'_regridded_allmonth_allregion_base`by'_round`rd'_record_part`p'.dta, replace
 
  * create season X allregion data
  forvalues s = 1/4 {
  use                     final_`vv'_`md'_regridded_season`s'_region1_base`by'_round`rd'_record_part`p'.dta, clear
  forvalues rr = 2/`nloc' {
   append using  final_`vv'_`md'_regridded_season`s'_region`rr'_base`by'_round`rd'_record_part`p'.dta
  }
  save                   final_`vv'_`md'_regridded_season`s'_allregion_base`by'_round`rd'_record_part`p'.dta, replace
  }
 }
}





*************** 1. compute RED : global (allmonth X allregion) ******************

local sp = "ssp585"
local vv =  "tmax"
local md = "CMCC-ESM2"
local nloc = 69
local rd   = 0.001
local by   = 1989

cd /media/tack/8GBlack/CMIP6_fromERA5/original/`md'_`vv'_seasonXregion/

* number of years in each part *
local ny  = 11

clear
gen year = .
save   /media/tack/8GBlack/CMIP6_fromERA5/final_seasonXregion_`sp'/final_`md'_RED_regridded_allmonthXallregion_base`by'_round`rd'.dta, replace


forvalues p = 1/11 {

use final_`vv'_`md'_regridded_allmonth_allregion_base`by'_round`rd'_record_part`p'.dta, clear

gen wght = cos(lat*_pi/180)

sort year
bysort year: egen t_wght = total(wght)
bysort year: gen  f_wght = (wght/t_wght)
bysort year: egen s_wght = total(f_wght)
*sum s_wght

bysort year: egen sf_wght = total(f_wght)

 
di " ----- part `p'  ------"

qui {

bysort year: gen b1 = (f_wght/sf_wght)*recordmax
bysort year: gen b2 = (f_wght/sf_wght)*recordmin
bysort year: gen b3 = (f_wght/sf_wght)*`vv'

bysort year: egen mh    = total(b1)
bysort year: egen ml    = total(b2)
bysort year: egen mean1 = total(b3)
bysort year: egen mean2 = mean(wght*`vv')

bysort year: keep if _n == 1

replace mh = 1 if _n == 1
replace ml = 1 if _n == 1

drop b1 b2 b3



*  s = c, t = lamda, x = RED
forvalues d = 1/`ny' {
 gen s`d' = 0
}
replace s1 = 1

forvalues d = 1/`ny' {
 gen x`d' = 0
}
replace x1 = 1 if _n == 1

forvalues d = 1/`ny' {
 gen t`d' = 0
}
replace t1 = s1


sort year

forvalues d = 2/`ny' {

 local d1 = `d' - 1
 replace x`d' =  mh*t`d1'/(1-mh) if _n == `d'
 replace s`d' = x`d'[`d'] if _n >= `d'

 forvalues i = 1/`d' {
  replace t`d' = t`d' + s`i'
 }
}

egen redh  = rowtotal(x1-x`ny')

drop s1-s`ny' x1-x`ny' t1-t`ny'


forvalues d = 1/`ny' {
 gen s`d' = 0
}
replace s1 = 1

forvalues d = 1/`ny' {
 gen x`d' = 0
}
replace x1 = 1 if _n == 1

forvalues d = 1/`ny' {
 gen t`d' = 0
}
replace t1 = s1


sort year

forvalues d = 2/`ny' {

 local d1 = `d' - 1
 replace x`d' =  ml*t`d1'/(1-ml) if _n == `d'
 replace s`d' = x`d'[`d'] if _n >= `d'

 forvalues i = 1/`d' {
  replace t`d' = t`d' + s`i'
 }
}

egen redl  = rowtotal(x1-x`ny')
drop s1-s`ny' x1-x`ny' t1-t`ny'

keep year mh ml mean1 mean2 redh redl

gen part = `p'

sort year

append using /media/tack/8GBlack/CMIP6_fromERA5/final_seasonXregion_`sp'/final_`md'_RED_regridded_allmonthXallregion_base`by'_round`rd'.dta
save                /media/tack/8GBlack/CMIP6_fromERA5/final_seasonXregion_`sp'/final_`md'_RED_regridded_allmonthXallregion_base`by'_round`rd'.dta, replace

}
}


cd /media/tack/8GBlack/CMIP6_fromERA5/final_seasonXregion_`sp'/

use final_`sp'_`md'_RED_regridded_allmonthXallregion_base`by'_round`rd'.dta, clear


sort year part
bysort year (part): keep if _n == 1

sort year
forvalues y = 1999(10)2089 {
 gen p`y'     = redh       if year >= `y'
 replace p`y' = p`y'[_n-1] if p`y'[_n-1] ! =.
 replace p`y' = 1          if p`y' ==. | year == `y'
}

forvalues y = 1999(10)2089 {
 gen q`y'     = redl       if year >= `y'
 replace q`y' = q`y'[_n-1] if q`y'[_n-1] ! =.
 replace q`y' = 1          if q`y' ==. | year == `y'
}

gen rh = redh*p1999*p2009*p2019*p2029*p2039*p2049*p2059*p2069*p2079*p2089
gen rl   = redl*q1999*q2009*q2019*q2029*q2039*q2049*q2059*q2069*q2079*q2089

replace redh = rh
replace redl = rl

duplicates report year part

ren mean1 tmax
ren redh redh_day
ren redl redl_day
drop mean2

keep  year tmax redh_day redl_day
order year tmax redh_day redl_day

save  final_`sp'_`md'_RED_regridded_allmonthXallregion_base`by'_round`rd'.dta, replace

tabstat redh_day, s(mean) by(year)







*************** 2. compute RED : by region (allmonth x region) ******************

local nloc = 69
local rd   = 0.001
local by   = 1989
local vv = "tmax"
local md = "CMCC-ESM2"

cd /media/tack/12TBHelio/CMIP6/`md'_`vv'_seasonXregion/

* number of years in each part *
local ny  = 11


clear
gen year = .
save  /media/tack/12TBHelio/CMIP6/final_seasonXregion/final_`md'_RED_regridded_allmonthXregion_base`by'_round`rd'.dta, replace


forvalues p = 1/11 {
 forvalues r = 1/`nloc' {

use final_`vv'_`md'_regridded_allmonth_region`r'_base`by'_round`rd'_record_part`p'.dta, clear

gen wght = cos(lat*_pi/180)

sort year
bysort year: egen t_wght = total(wght)
bysort year: gen  f_wght = (wght/t_wght)
bysort year: egen s_wght = total(f_wght)
*sum s_wght

bysort year: egen sf_wght = total(f_wght)

 
di " ----- part `p' , region `r'  ------"

qui {

bysort year: gen b1 = (f_wght/sf_wght)*recordmax
bysort year: gen b2 = (f_wght/sf_wght)*recordmin
bysort year: gen b3 = (f_wght/sf_wght)*`vv'

bysort year: egen mh    = total(b1)
bysort year: egen ml    = total(b2)
bysort year: egen mean1 = total(b3)
bysort year: egen mean2 = mean(wght*`vv')

bysort year: keep if _n == 1

replace mh = 1 if _n == 1
replace ml = 1 if _n == 1

drop b1 b2 b3



*  s = c, t = lamda, x = RED
forvalues d = 1/`ny' {
 gen s`d' = 0
}
replace s1 = 1

forvalues d = 1/`ny' {
 gen x`d' = 0
}
replace x1 = 1 if _n == 1

forvalues d = 1/`ny' {
 gen t`d' = 0
}
replace t1 = s1


sort year

forvalues d = 2/`ny' {

 local d1 = `d' - 1
 replace x`d' =  mh*t`d1'/(1-mh) if _n == `d'
 replace s`d' = x`d'[`d'] if _n >= `d'

 forvalues i = 1/`d' {
  replace t`d' = t`d' + s`i'
 }
}

egen redh  = rowtotal(x1-x`ny')

drop s1-s`ny' x1-x`ny' t1-t`ny'


forvalues d = 1/`ny' {
 gen s`d' = 0
}
replace s1 = 1

forvalues d = 1/`ny' {
 gen x`d' = 0
}
replace x1 = 1 if _n == 1

forvalues d = 1/`ny' {
 gen t`d' = 0
}
replace t1 = s1


sort year

forvalues d = 2/`ny' {

 local d1 = `d' - 1
 replace x`d' =  ml*t`d1'/(1-ml) if _n == `d'
 replace s`d' = x`d'[`d'] if _n >= `d'

 forvalues i = 1/`d' {
  replace t`d' = t`d' + s`i'
 }
}

egen redl  = rowtotal(x1-x`ny')
drop s1-s`ny' x1-x`ny' t1-t`ny'

keep year mh ml mean1 mean2 redh redl

gen part    = `p'
gen id      = `r'

sort year

append using /media/tack/12TBHelio/CMIP6/final_seasonXregion/final_`md'_RED_regridded_allmonthXregion_base`by'_round`rd'.dta
save                /media/tack/12TBHelio/CMIP6/final_seasonXregion/final_`md'_RED_regridded_allmonthXregion_base`by'_round`rd'.dta, replace

}
}
}

cd  /media/tack/12TBHelio/CMIP6/final_seasonXregion/

use final_`md'_RED_regridded_allmonthXregion_base`by'_round`rd'.dta, clear

sort id year part
bysort id year (part): keep if _n == 1

forvalues y = 1999(10)2089 {
 bysort id (year): gen p`y'     = redh       if year >= `y'
 bysort id (year): replace p`y' = p`y'[_n-1] if p`y'[_n-1] ! =.
 bysort id (year): replace p`y' = 1          if p`y' ==. | year == `y'
}

forvalues y = 1999(10)2089 {
 bysort id (year): gen q`y'     = redl       if year >= `y'
 bysort id (year): replace q`y' = q`y'[_n-1] if q`y'[_n-1] ! =.
 bysort id (year): replace q`y' = 1          if q`y' ==. | year == `y'
}

gen rh = redh*p1999*p2009*p2019*p2029*p2039*p2049*p2059*p2069*p2079*p2089
gen rl  =  redl*q1999*q2009*q2019*q2029*q2039*q2049*q2059*q2069*q2079*q2089

replace redh = rh
replace redl = rl

duplicates report id year part

ren mean1 tmax
ren redh redh_day
ren redl redl_day
drop mean2

keep  id year tmax redh_day redl_day
order id year tmax redh_day redl_day

save final_`md'_RED_regridded_allmonthXregion_base`by'_round`rd'_clean.dta, replace

tabstat redh_day, s(mean) by(year)





*************** 3. compute RED : season X allregion ******************

local nloc = 69
local rd   = 0.001
local by   = 1989
local vv = "tmax"
local md = "CMCC-ESM2"

cd /media/tack/12TBHelio/CMIP6/`md'_`vv'_seasonXregion/

* number of years in each part *
local ny  = 11

clear
gen year = .
save /media/tack/12TBHelio/CMIP6/final_seasonXregion/final_`md'_RED_regridded_seasonXallregion_base`by'_round`rd'.dta, replace


forvalues pp = 1/11 {
 forvalues ss = 1/4 {

use final_`vv'_`md'_regridded_season`ss'_allregion_base`by'_round`rd'_record_part`pp'.dta, clear

gen wght = cos(lat*_pi/180)

sort year
bysort year: egen t_wght = total(wght)
bysort year: gen  f_wght = (wght/t_wght)
bysort year: egen s_wght = total(f_wght)
*sum s_wght

bysort year: egen sf_wght = total(f_wght)

 
di " ----- part `pp', season `ss' ------"

qui {

bysort year: gen b1 = (f_wght/sf_wght)*recordmax
bysort year: gen b2 = (f_wght/sf_wght)*recordmin
bysort year: gen b3 = (f_wght/sf_wght)*`vv'

bysort year: egen mh    = total(b1)
bysort year: egen ml    = total(b2)
bysort year: egen mean1 = total(b3)
bysort year: egen mean2 = mean(wght*`vv')

bysort year: keep if _n == 1

replace mh = 1 if _n == 1
replace ml = 1 if _n == 1

drop b1 b2 b3



*  s = c, t = lamda, x = RED
forvalues d = 1/`ny' {
 gen s`d' = 0
}
replace s1 = 1

forvalues d = 1/`ny' {
 gen x`d' = 0
}
replace x1 = 1 if _n == 1

forvalues d = 1/`ny' {
 gen t`d' = 0
}
replace t1 = s1


sort year

forvalues d = 2/`ny' {

 local d1 = `d' - 1
 replace x`d' =  mh*t`d1'/(1-mh) if _n == `d'
 replace s`d' = x`d'[`d'] if _n >= `d'

 forvalues i = 1/`d' {
  replace t`d' = t`d' + s`i'
 }
}

egen redh  = rowtotal(x1-x`ny')

drop s1-s`ny' x1-x`ny' t1-t`ny'


forvalues d = 1/`ny' {
 gen s`d' = 0
}
replace s1 = 1

forvalues d = 1/`ny' {
 gen x`d' = 0
}
replace x1 = 1 if _n == 1

forvalues d = 1/`ny' {
 gen t`d' = 0
}
replace t1 = s1


sort year

forvalues d = 2/`ny' {

 local d1 = `d' - 1
 replace x`d' =  ml*t`d1'/(1-ml) if _n == `d'
 replace s`d' = x`d'[`d'] if _n >= `d'

 forvalues i = 1/`d' {
  replace t`d' = t`d' + s`i'
 }
}

egen redl  = rowtotal(x1-x`ny')
drop s1-s`ny' x1-x`ny' t1-t`ny'

keep year mh ml mean1 mean2 redh redl

gen part    = `pp'
gen season  = `ss'

sort year

append using /media/tack/12TBHelio/CMIP6/final_seasonXregion/final_`md'_RED_regridded_seasonXallregion_base`by'_round`rd'.dta
save         /media/tack/12TBHelio/CMIP6/final_seasonXregion/final_`md'_RED_regridded_seasonXallregion_base`by'_round`rd'.dta, replace

}
}
}



cd   /media/tack/12TBHelio/CMIP6/final_seasonXregion/

use final_`md'_RED_regridded_seasonXallregion_base`by'_round`rd'.dta, clear

sort season year part
bysort season year (part): keep if _n == 1

forvalues y = 1999(10)2089 {
 bysort season (year): gen p`y'     = redh       if year >= `y'
 bysort season (year): replace p`y' = p`y'[_n-1] if p`y'[_n-1] ! =.
 bysort season (year): replace p`y' = 1          if p`y' ==. | year == `y'
}

forvalues y = 1999(10)2089 {
 bysort season (year): gen q`y'     = redl       if year >= `y'
 bysort season (year): replace q`y' = q`y'[_n-1] if q`y'[_n-1] ! =.
 bysort season (year): replace q`y' = 1          if q`y' ==. | year == `y'
}

gen rh = redh*p1999*p2009*p2019*p2029*p2039*p2049*p2059*p2069*p2079*p2089
gen rl   = redl*q1999*q2009*q2019*q2029*q2039*q2049*q2059*q2069*q2079*q2089

replace redh = rh
replace redl = rl

duplicates report season year part

ren mean1 tmax
ren redh redh_day
ren redl redl_day
drop mean2

keep  season year tmax redh_day redl_day
order  season year tmax redh_day redl_day

save final_`md'_RED_regridded_seasonXallregion_base`by'_round`rd'_clean.dta, replace

tabstat redh_day, s(mean) by(year)







*************** 4. compute RED : season  x  region ******************

local nloc = 69
local rd   = 0.001
local by   = 1989
local vv = "tmax"
local md = "CMCC-ESM2"

cd /media/tack/12TBHelio/CMIP6/`md'_`vv'_seasonXregion/

* number of years in each part *
local ny  = 11

clear
gen year = .
save /media/tack/12TBHelio/CMIP6/final_seasonXregion/final_`md'_RED_regridded_seasonXregion_base`by'_round`rd'.dta, replace


forvalues pp = 1/11 {
 forvalues ss = 1/4 {
  forvalues rr = 1/`nloc' {

use final_`vv'_`md'_regridded_season`ss'_region`rr'_base`by'_round`rd'_record_part`pp'.dta, clear

gen wght = cos(lat*_pi/180)

sort year
bysort year: egen t_wght = total(wght)
bysort year: gen  f_wght = (wght/t_wght)
bysort year: egen s_wght = total(f_wght)
*sum s_wght

bysort year: egen sf_wght = total(f_wght)

di " ----- part `pp',  season `ss', region `rr'  ------"

qui {

bysort year: gen b1 = (f_wght/sf_wght)*recordmax
bysort year: gen b2 = (f_wght/sf_wght)*recordmin
bysort year: gen b3 = (f_wght/sf_wght)*`vv'

bysort year: egen mh    = total(b1)
bysort year: egen ml    = total(b2)
bysort year: egen mean1 = total(b3)
bysort year: egen mean2 = mean(wght*`vv')

bysort year: keep if _n == 1

replace mh = 1 if _n == 1
replace ml = 1 if _n == 1

drop b1 b2 b3



*  s = c, t = lamda, x = RED
forvalues d = 1/`ny' {
 gen s`d' = 0
}
replace s1 = 1

forvalues d = 1/`ny' {
 gen x`d' = 0
}
replace x1 = 1 if _n == 1

forvalues d = 1/`ny' {
 gen t`d' = 0
}
replace t1 = s1


sort year

forvalues d = 2/`ny' {

 local d1 = `d' - 1
 replace x`d' =  mh*t`d1'/(1-mh) if _n == `d'
 replace s`d' = x`d'[`d'] if _n >= `d'

 forvalues i = 1/`d' {
  replace t`d' = t`d' + s`i'
 }
}

egen redh  = rowtotal(x1-x`ny')

drop s1-s`ny' x1-x`ny' t1-t`ny'


forvalues d = 1/`ny' {
 gen s`d' = 0
}
replace s1 = 1

forvalues d = 1/`ny' {
 gen x`d' = 0
}
replace x1 = 1 if _n == 1

forvalues d = 1/`ny' {
 gen t`d' = 0
}
replace t1 = s1


sort year

forvalues d = 2/`ny' {

 local d1 = `d' - 1
 replace x`d' =  ml*t`d1'/(1-ml) if _n == `d'
 replace s`d' = x`d'[`d'] if _n >= `d'

 forvalues i = 1/`d' {
  replace t`d' = t`d' + s`i'
 }
}

egen redl  = rowtotal(x1-x`ny')
drop s1-s`ny' x1-x`ny' t1-t`ny'

keep year mh ml mean1 mean2 redh redl

gen part    = `pp'
gen season  = `ss'
gen id      = `rr'

sort year

append using /media/tack/12TBHelio/CMIP6/final_seasonXregion/final_`md'_RED_regridded_seasonXregion_base`by'_round`rd'.dta
save         /media/tack/12TBHelio/CMIP6/final_seasonXregion/final_`md'_RED_regridded_seasonXregion_base`by'_round`rd'.dta, replace

}
}
}
}


cd /media/tack/12TBHelio/CMIP6/final_seasonXregion/

use final_`md'_RED_regridded_seasonXregion_base`by'_round`rd'.dta, clear

sort id season year part
bysort id season year (part): keep if _n == 1

forvalues y = 1999(10)2089 {
 bysort id season (year): gen p`y'     = redh       if year >= `y'
 bysort id season (year): replace p`y' = p`y'[_n-1] if p`y'[_n-1] ! =.
 bysort id season (year): replace p`y' = 1          if p`y' ==. | year == `y'
}

forvalues y = 1999(10)2089 {
 bysort id season (year): gen q`y'     = redl       if year >= `y'
 bysort id season (year): replace q`y' = q`y'[_n-1] if q`y'[_n-1] ! =.
 bysort id season (year): replace q`y' = 1          if q`y' ==. | year == `y'
}

gen rh = redh*p1999*p2009*p2019*p2029*p2039*p2049*p2059*p2069*p2079*p2089
gen rl  = redl*q1999*q2009*q2019*q2029*q2039*q2049*q2059*q2069*q2079*q2089

replace redh = rh
replace redl = rl

duplicates report id season year part

ren mean1 tmax
ren redh redh_day
ren redl redl_day
drop mean2

keep  id season year tmax redh_day redl_day
order  id season year tmax redh_day redl_day

save final_`md'_RED_regridded_seasonXregion_base`by'_round`rd'_clean.dta, replace

tabstat redh_day, s(mean) by(year)








*************** 5. compute RED : month  x  region ******************

local vv = "tmax"
local md = "CMCC-ESM2"
local sp = "ssp585"
local nloc  = 69
local rd = 0.001
local by = 1989
*local mlist = "01 02 03 04 05 06 07 08 09 10 11 12"

cd /media/tack/8GBlack/CMIP6_fromERA5/original/`md'_`vv'_seasonXregion/

* number of years in each part *
local ny  = 11

clear
gen year = .
save  /media/tack/8GBlack/CMIP6_fromERA5/final_seasonXregion_`sp'/final_`sp'_`md'_RED_regridded_monthXregion_base`by'_round`rd'.dta, replace


forvalues pp = 1/11 {
 forvalues mm = 1/12 {
  forvalues rr = 1/`nloc' {

use final_`vv'_`md'_regridded_month`mm'_region`rr'_base`by'_round`rd'_record_part`pp'.dta, clear

gen wght = cos(lat*_pi/180)

sort year
bysort year: egen t_wght = total(wght)
bysort year: gen  f_wght = (wght/t_wght)
bysort year: egen s_wght = total(f_wght)
*sum s_wght

bysort year: egen sf_wght = total(f_wght)

di " ----- part `pp',  month `mm', region `rr'  ------"

qui {

bysort year: gen b1 = (f_wght/sf_wght)*recordmax
bysort year: gen b2 = (f_wght/sf_wght)*recordmin
bysort year: gen b3 = (f_wght/sf_wght)*`vv'

bysort year: egen mh    = total(b1)
bysort year: egen ml    = total(b2)
bysort year: egen mean1 = total(b3)
bysort year: egen mean2 = mean(wght*`vv')

bysort year: keep if _n == 1

replace mh = 1 if _n == 1
replace ml = 1 if _n == 1

drop b1 b2 b3



*  s = c, t = lamda, x = RED
forvalues d = 1/`ny' {
 gen s`d' = 0
}
replace s1 = 1

forvalues d = 1/`ny' {
 gen x`d' = 0
}
replace x1 = 1 if _n == 1

forvalues d = 1/`ny' {
 gen t`d' = 0
}
replace t1 = s1


sort year

forvalues d = 2/`ny' {

 local d1 = `d' - 1
 replace x`d' =  mh*t`d1'/(1-mh) if _n == `d'
 replace s`d' = x`d'[`d'] if _n >= `d'

 forvalues i = 1/`d' {
  replace t`d' = t`d' + s`i'
 }
}

egen redh  = rowtotal(x1-x`ny')

drop s1-s`ny' x1-x`ny' t1-t`ny'


forvalues d = 1/`ny' {
 gen s`d' = 0
}
replace s1 = 1

forvalues d = 1/`ny' {
 gen x`d' = 0
}
replace x1 = 1 if _n == 1

forvalues d = 1/`ny' {
 gen t`d' = 0
}
replace t1 = s1


sort year

forvalues d = 2/`ny' {

 local d1 = `d' - 1
 replace x`d' =  ml*t`d1'/(1-ml) if _n == `d'
 replace s`d' = x`d'[`d'] if _n >= `d'

 forvalues i = 1/`d' {
  replace t`d' = t`d' + s`i'
 }
}

egen redl  = rowtotal(x1-x`ny')
drop s1-s`ny' x1-x`ny' t1-t`ny'

keep year mh ml mean1 mean2 redh redl

gen part    = `pp'
gen month  = `mm'
gen id      = `rr'

sort year

append using /media/tack/8GBlack/CMIP6_fromERA5/final_seasonXregion_`sp'/final_`sp'_`md'_RED_regridded_monthXregion_base`by'_round`rd'.dta
save                /media/tack/8GBlack/CMIP6_fromERA5/final_seasonXregion_`sp'/final_`sp'_`md'_RED_regridded_monthXregion_base`by'_round`rd'.dta, replace
  
}
}
}
}


cd /media/tack/8GBlack/CMIP6_fromERA5/final_seasonXregion_`sp'/

use final_`sp'_`md'_RED_regridded_monthXregion_base`by'_round`rd'.dta, clear

sort id month year part
bysort id month year (part): keep if _n == 1

forvalues y = 1999(10)2089 {
 bysort id month (year): gen p`y'     = redh       if year >= `y'
 bysort id month (year): replace p`y' = p`y'[_n-1] if p`y'[_n-1] ! =.
 bysort id month (year): replace p`y' = 1          if p`y' ==. | year == `y'
}

forvalues y = 1999(10)2089 {
 bysort id month (year): gen q`y'     = redl       if year >= `y'
 bysort id month (year): replace q`y' = q`y'[_n-1] if q`y'[_n-1] ! =.
 bysort id month (year): replace q`y' = 1          if q`y' ==. | year == `y'
}

gen rh = redh*p1999*p2009*p2019*p2029*p2039*p2049*p2059*p2069*p2079*p2089
gen rl  = redl*q1999*q2009*q2019*q2029*q2039*q2049*q2059*q2069*q2079*q2089

replace redh = rh
replace redl = rl

duplicates report id month year part

ren mean1 tmax
ren redh redh_day
ren redl redl_day
drop mean2

keep  id month year tmax redh_day redl_day
order id month year tmax redh_day redl_day

save final_`sp'_`md'_RED_regridded_monthXregion_base`by'_round`rd'.dta, replace

tabstat redh_day, s(mean) by(year) 

tabstat redh_day, s(min p5 p10 p25 p50 p75 p90 p99 max)





*************** 6. compute RED : regional (allmonth x group) ******************

local sp = "ssp585"
local vv = "tmax"
local md = "CMCC-ESM2"
local nloc  = 69
local rd = 0.001
local by = 1989

cd /media/tack/8GBlack/CMIP6_fromERA5/original/`md'_`vv'_seasonXregion/

* number of years in each part *
local ny  = 11

clear
gen year = .
save  /media/tack/8GBlack/CMIP6_fromERA5/final_seasonXregion_`sp'/final_`sp'_`md'_RED_regridded_allmonthXgroup_base`by'_round`rd'.dta, replace


forvalues p = 1/11 {
 forvalues g = 1/4 {
 
qui { 
use final_`vv'_`md'_regridded_allmonth_allregion_base`by'_round`rd'_record_part`p'.dta, clear

* breaking by latitude
gen group = 1
replace group = 2 if coarse_center_lat > 30 & coarse_center_lat <=60
replace group = 3 if coarse_center_lat > 0   & coarse_center_lat <=30
replace group = 4 if coarse_center_lat <=0

keep if group == `g'

gen wght = cos(lat*_pi/180)

sort year
bysort year: egen t_wght = total(wght)
bysort year: gen  f_wght = (wght/t_wght)
bysort year: egen s_wght = total(f_wght)
*sum s_wght

bysort year: egen sf_wght = total(f_wght)
}
 
di " ----- part `p' , group`g'  ------"

qui {

bysort year: gen b1 = (f_wght/sf_wght)*recordmax
bysort year: gen b2 = (f_wght/sf_wght)*recordmin
bysort year: gen b3 = (f_wght/sf_wght)*`vv'

bysort year: egen mh   = total(b1)
bysort year: egen ml    = total(b2)
bysort year: egen mean1 = total(b3)
bysort year: egen mean2 = mean(wght*`vv')

bysort year: keep if _n == 1

replace mh = 1 if _n == 1
replace ml = 1 if _n == 1

drop b1 b2 b3



*  s = c, t = lamda, x = RED
forvalues d = 1/`ny' {
 gen s`d' = 0
}
replace s1 = 1

forvalues d = 1/`ny' {
 gen x`d' = 0
}
replace x1 = 1 if _n == 1

forvalues d = 1/`ny' {
 gen t`d' = 0
}
replace t1 = s1


sort year

forvalues d = 2/`ny' {

 local d1 = `d' - 1
 replace x`d' =  mh*t`d1'/(1-mh) if _n == `d'
 replace s`d' = x`d'[`d'] if _n >= `d'

 forvalues i = 1/`d' {
  replace t`d' = t`d' + s`i'
 }
}

egen redh  = rowtotal(x1-x`ny')

drop s1-s`ny' x1-x`ny' t1-t`ny'


forvalues d = 1/`ny' {
 gen s`d' = 0
}
replace s1 = 1

forvalues d = 1/`ny' {
 gen x`d' = 0
}
replace x1 = 1 if _n == 1

forvalues d = 1/`ny' {
 gen t`d' = 0
}
replace t1 = s1


sort year

forvalues d = 2/`ny' {

 local d1 = `d' - 1
 replace x`d' =  ml*t`d1'/(1-ml) if _n == `d'
 replace s`d' = x`d'[`d'] if _n >= `d'

 forvalues i = 1/`d' {
  replace t`d' = t`d' + s`i'
 }
}

egen redl  = rowtotal(x1-x`ny')
drop s1-s`ny' x1-x`ny' t1-t`ny'

keep year mh ml mean1 mean2 redh redl

gen part     = `p'
gen group  = `g'

sort year

append using  /media/tack/8GBlack/CMIP6_fromERA5/final_seasonXregion_`sp'/final_`sp'_`md'_RED_regridded_allmonthXgroup_base`by'_round`rd'.dta
save                 /media/tack/8GBlack/CMIP6_fromERA5/final_seasonXregion_`sp'/final_`sp'_`md'_RED_regridded_allmonthXgroup_base`by'_round`rd'.dta, replace

}
}
}


cd   /media/tack/8GBlack/CMIP6_fromERA5/final_seasonXregion_`sp'/

use final_`sp'_`md'_RED_regridded_allmonthXgroup_base`by'_round`rd'.dta, clear

sort group year part
bysort group year (part): keep if _n == 1

forvalues y = 1999(10)2089 {
 bysort group (year): gen p`y'     = redh            if year >= `y'
 bysort group (year): replace p`y' = p`y'[_n-1] if p`y'[_n-1] ! =.
 bysort group (year): replace p`y' = 1               if p`y' ==. | year == `y'
}

forvalues y = 1999(10)2089 {
 bysort group (year): gen q`y'     = redl             if year >= `y'
 bysort group (year): replace q`y' = q`y'[_n-1] if q`y'[_n-1] ! =.
 bysort group (year): replace q`y' = 1               if q`y' ==. | year == `y'
}

gen rh =  redh*p1999*p2009*p2019*p2029*p2039*p2049*p2059*p2069*p2079*p2089
gen rl  =  redl*q1999*q2009*q2019*q2029*q2039*q2049*q2059*q2069*q2079*q2089

replace redh = rh
replace redl = rl

duplicates report group year part

ren mean1 tmax
ren redh redh_day
ren redl redl_day
drop mean2

keep  group year tmax redh_day redl_day
order group year tmax redh_day redl_day

save final_`sp'_`md'_RED_regridded_allmonthXgroup_base`by'_round`rd'.dta, replace

tabstat redh_day, s(mean max) by(year)



******* single time series data and record  **********

local sp = "ssp585"
local vv =  "tmax"
local md = "CMCC-ESM2"
local nloc = 69
local rd = 0.001
local by = 1989

cd /media/tack/8GBlack/CMIP6_fromERA5/original/`md'_`vv'_seasonXregion/


forvalues m = 1/12 {
  forvalues r = 1/`nloc' {

di " ----------  month `m', region `r'  ----------- "
 
qui {

use  final_`vv'_`md'_regridded_month`m'_region`r'_base`by'.dta, clear

keep if year >= 1989  & year <= 2099

* adding record indicator *
gen `vv'_r = round(`vv',`rd')
replace `vv' = `vv'_r

egen loc = group(lat lon)
egen loc_time = group(lat lon month shuffle day)

sort loc month shuffle day year

sort loc month shuffle day year
bysort loc month shuffle day (year): gen bestt = `vv' if _n == 1
bysort loc month shuffle day (year): replace bestt = max(`vv', bestt[_n-1]) if loc_time == loc_time[_n-1]
bysort loc month shuffle day (year): gen recordmax = 0
bysort loc month shuffle day (year): replace recordmax = 1 if loc_time != loc_time[_n-1] | (bestt > bestt[_n-1])
bysort loc month shuffle day (year): gen recordmax_gap =  .
bysort loc month shuffle day (year): replace recordmax_gap =  bestt - bestt[_n-1] if loc_time != loc_time[_n-1] | (bestt > bestt[_n-1])

drop bestt

sort loc month shuffle day year
bysort loc month shuffle day (year): gen bestt = `vv' if _n == 1
bysort loc month shuffle day (year): replace bestt = min(`vv', bestt[_n-1]) if loc_time == loc_time[_n-1]
bysort loc month shuffle day (year): gen recordmin = 0
bysort loc month shuffle day (year): replace recordmin = 1 if loc_time != loc_time[_n-1] | (bestt < bestt[_n-1])
bysort loc month shuffle day (year): gen recordmin_gap =  .
bysort loc month shuffle day (year): replace recordmin_gap =  bestt - bestt[_n-1] if loc_time != loc_time[_n-1] | (bestt < bestt[_n-1])

drop bestt loc_time
 
save final_`vv'_`md'_regridded_month`m'_region`r'_base`by'_round`rd'_record_single.dta, replace
}
}
}




**************  all combined : create allmonth, season, allmonth X allregion, season x allregion  *****************

local sp = "ssp585"
local vv =  "tmax"
local md = "CMCC-ESM2"
local nloc = 69
local rd = 0.001
local by = 1989

cd /media/tack/8GBlack/CMIP6_fromERA5/original/`md'_`vv'_seasonXregion/


* create allmonth data
forvalues r = 1/`nloc' {
 
    use                  final_`vv'_`md'_regridded_month1_region`r'_base`by'_round`rd'_record_single.dta, clear
    append using final_`vv'_`md'_regridded_month2_region`r'_base`by'_round`rd'_record_single.dta
    append using final_`vv'_`md'_regridded_month12_region`r'_base`by'_round`rd'_record_single.dta
    save                final_`vv'_`md'_regridded_season1_region`r'_base`by'_round`rd'_record_single.dta, replace

   use                  final_`vv'_`md'_regridded_month3_region`r'_base`by'_round`rd'_record_single.dta, clear
   append using final_`vv'_`md'_regridded_month4_region`r'_base`by'_round`rd'_record_single.dta
   append using final_`vv'_`md'_regridded_month5_region`r'_base`by'_round`rd'_record_single.dta
   save                final_`vv'_`md'_regridded_season2_region`r'_base`by'_round`rd'_record_single.dta, replace

   use                  final_`vv'_`md'_regridded_month6_region`r'_base`by'_round`rd'_record_single.dta, clear
   append using final_`vv'_`md'_regridded_month7_region`r'_base`by'_round`rd'_record_single.dta
   append using final_`vv'_`md'_regridded_month8_region`r'_base`by'_round`rd'_record_single.dta
   save                final_`vv'_`md'_regridded_season3_region`r'_base`by'_round`rd'_record_single.dta, replace
   
   use                  final_`vv'_`md'_regridded_month9_region`r'_base`by'_round`rd'_record_single.dta, clear
   append using final_`vv'_`md'_regridded_month10_region`r'_base`by'_round`rd'_record_single.dta
   append using final_`vv'_`md'_regridded_month11_region`r'_base`by'_round`rd'_record_single.dta
   save                final_`vv'_`md'_regridded_season4_region`r'_base`by'_round`rd'_record_single.dta, replace
 
  append using final_`vv'_`md'_regridded_season1_region`r'_base`by'_round`rd'_record_single.dta
  append using final_`vv'_`md'_regridded_season2_region`r'_base`by'_round`rd'_record_single.dta
  append using final_`vv'_`md'_regridded_season3_region`r'_base`by'_round`rd'_record_single.dta
  save                final_`vv'_`md'_regridded_allmonth_region`r'_base`by'_round`rd'_record_single.dta, replace

}
 
  * create allmonth X allregion data
use                     final_`vv'_`md'_regridded_allmonth_region1_base`by'_round`rd'_record_single.dta, clear
forvalues rr = 2/`nloc' {
  append using  final_`vv'_`md'_regridded_allmonth_region`rr'_base`by'_round`rd'_record_single.dta
}
save                  final_`vv'_`md'_regridded_allmonth_allregion_base`by'_round`rd'_record_single.dta, replace
 
  * create season X allregion data
forvalues s = 1/4 {
  use                     final_`vv'_`md'_regridded_season`s'_region1_base`by'_round`rd'_record_single.dta, clear
  forvalues rr = 2/`nloc' {
   append using  final_`vv'_`md'_regridded_season`s'_region`rr'_base`by'_round`rd'_record_single.dta
  }
  save                   final_`vv'_`md'_regridded_season`s'_allregion_base`by'_round`rd'_record_single.dta, replace
}



****  compute the mean recordmax and min gaps ****

local sp = "ssp585"
local vv = "tmax"
local md = "CMCC-ESM2"
local nloc  = 69
local rd = 0.001
local by = 1989

cd /media/tack/8GBlack/CMIP6_fromERA5/original/`md'_`vv'_seasonXregion/

use  final_`vv'_`md'_regridded_allmonth_allregion_base`by'_round`rd'_record_single.dta, clear

tabstat recordmax_gap, s(min p5 p10 p50 p75 p90 p95 p99 max)

* breaking by latitude
gen group = 1
replace group = 2 if coarse_center_lat > 30 & coarse_center_lat <=60
replace group = 3 if coarse_center_lat > 0   & coarse_center_lat <=30
replace group = 4 if coarse_center_lat <=0

*tabstat recordmax_gap if year >= 2070, s(mean sd p5 p50 p75 p90 p95) by(group)
*tabstat recordmax_gap if year < 2024, s(mean sd p5 p50 p75 p90 p95) by(group)
*tabstat recordmax_gap if group == 1, s(mean) by(year)

*egen id2 =  group(lat lon)
*sort group year  month day id2
*bysort group year  month day (id2): gen dd = _N
*tabstat dd, s(mean) by(group)


* 1. average recordmax_gap and recordmin_gap over locations, month, and day 
*sort group year 
*bysort group year: egen m_rmax_gap = mean(recordmax_gap)
*bysort group year: egen m_rmin_gap  = mean(recordmin_gap)
*bysort group year: keep if _n == 1

*tabstat m_rmax_gap, s(min p5 p10 p50 p75 p90 p95 p99 max)

* 2. weighted average
gen wght = cos(lat*_pi/180)

sort group year month day lat lon
bysort group year month day: egen t_wght = total(wght)
bysort group year month day: gen   f_wght = (wght/t_wght)
bysort group year month day: egen s_wght = total(f_wght)
*sum s_wght
bysort group year month day: egen sf_wght = total(f_wght)
sum sf_wght

bysort group year month day: gen b1 = (f_wght/sf_wght)*recordmax_gap
bysort group year month day: gen b2 = (f_wght/sf_wght)*recordmin_gap
bysort group year month day: egen s1 = total(b1)
bysort group year month day: egen s2 = total(b2)
bysort group year: egen m_rmax_gap = mean(s1)
bysort group year: egen m_rmin_gap  = mean(s2)
bysort group year: keep if _n == 1

keep group year m_rmax_gap m_rmin_gap

gen ssp = "`sp'"
gen model = "`md'"

save  /media/tack/8GBlack/CMIP6_fromERA5/final_seasonXregion_`sp'_gaps/final_`sp'_`vv'_`md'_regridded_base`by'_round`rd'_gaps.dta, replace



