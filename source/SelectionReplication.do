/* SelectionReplication.do        KTS/DCC/NLB               yyyy-mm-dd:2024-11-10
----|----1----|----2----|----3----|----4----|----5----|----6----|----7----|----8

This do file replicates the data needed for the selection analysis
*/

*-------------------------------------------------------------------------------
*--- Preamble
*-------------------------------------------------------------------------------

clear all
set more off
*capture log close
*set scheme s1mono
graph set window fontface "Times New Roman"
set matsize 11000
set seed 30081985

*-------------------------------------------------------------------------------
*--- Directories and Log 
*-------------------------------------------------------------------------------

global DIR "F:\intergenerational\Replication"
global DAT "$DIR\data"
global SEL "$DIR\data\Selection"
cap mkdir $SEL
global OUT "$DIR\results\Selection"
cap mkdir $OUT
*global LOG "$DIR\log\Selection"
*cap mkdir $LOG

* Log file:
*log using "$LOG\SelectionReplication.log", replace

*-------------------------------------------------------------------------------
*--- Variables
*-------------------------------------------------------------------------------

global g1blnctrls edmom EDAD_MADRE married doc_aten bregion_?? byear_????
global g2blnctrls edgmom EDAD_ABUELA marriedm doc_atenm bregion_??m byear_????m

*-------------------------------------------------------------------------------
*--- Mortality by 15
*-------------------------------------------------------------------------------

use "$DAT/workingdata_age_work.dta", clear

gen g1smpl = SEXO == 2 ///
	& EDAD_MADRE >= 15 & EDAD_MADRE <= 45 ///
	& mrg_mbdata2NAC == 1 ///
	& g1ok_blnctrls == 1
	
keep if g1smpl == 1
	
sum ANO_NAC
scalar max_ano_nac = r(max)
scalar min_ano_nac = r(min)

sum ANO_DEF
scalar max_ano_def = r(max)
scalar min_ano_def = r(min)

scalar max_age = max_ano_def - min_ano_nac - 1

 if 15 >= 1 & 15 <= `=max_age' {
  egen byte dead_by_a15 = rowtotal(dead_at_a00-dead_at_a14), m
  replace dead_by_a15 = . if ANO_NAC + 14 > max_ano_def

  #delimit ;
  rdrobust dead_by_a15 PESO if sem32 == 1, c(1500) scalepar(-1) all
           covs(${g1blnctrls} dbw1??0) vce(cluster PESO);
  #delimit cr
	 	  
  estimates save "$OUT/g1f_mortbya15_o32_hps.ster", replace
 } 
estimates clear

*-------------------------------------------------------------------------------
*--- Number of Births at Age
*-------------------------------------------------------------------------------

use "$DAT/workingdata_age_work.dta", clear

gen g1smpl = SEXO == 2 ///
	& EDAD_MADRE >= 15 & EDAD_MADRE <= 45 ///
	& mrg_mbdata2NAC == 1 ///
	& g1ok_blnctrls == 1
	
keep if g1smpl == 1

sum ANO_NAC
scalar max_ano_nac = r(max)
scalar min_ano_nac = r(min)

scalar max_age = max_ano_nac - min_ano_nac

foreach y of numlist 15/26 {

  if 15 >= 1 & 15 <= `=max_age' {
 
  #delimit ;
  rdrobust nbirths_at_`y' PESO if sem32 == 1, c(1500) scalepar(-1) all
  covs(${g1blnctrls} dbw1??0) vce(cluster PESO);
  #delimit cr

  estimates save "$OUT/g1f_nbirthsat`y'_o32_hps.ster", replace
  }
 }

*-------------------------------------------------------------------------------
*--- Counterfactual Mothers
*-------------------------------------------------------------------------------

* Load birth data:
use "$DAT\NAC_1992_2018_NOGLOSAS_NODUPS_NONAS.dta", clear

* Keep only girls born before 2003 (so that they could have been 15 by 2018):
keep if SEXO == 2 & ANO_NAC <= 2003

* Create 50 gram bins according to birth weight:
gen int bin50g = 50*floor(PESO / 50)
label var bin50g "Birth weight 50 gram bins"

********************************************************************************
* CALCULATE MORTALITY RATES:
* Merge in mortality variables (for girls born before 2003):
merge 1:1 ID_RECIEN_NACIDO using "$DAT\MORTALITY_VARIABLES.dta", ///
	nogen keepusing(dead_at_a??) keep(master match)

* Construct variable of death by age 15:
egen dead_by_a15 = rowtotal(dead_at_a00-dead_at_a14), m
label var dead_by_a15 "Dead before 15th birthday."

* Calculate mortality rates before age 15 by birth year and weight bin:
bys ANO_NAC bin50g: egen double mort_by_cell = mean(dead_by_a15)
label var mort_by_cell "Mortality rates before age 15 by birth year and weight bin."

********************************************************************************
* CALCULATE FERTILITY RATES:
* Merge in fertility variables:
merge 1:1 ID_RECIEN_NACIDO using "$DAT\FERTILITY_VARIABLES.dta", ///
	nogen keep(master match) keepusing(nbirths_at_??)

* Calculate birth rates for survivors by birth year, and weight bin:
foreach var of varlist nbirths_at_?? {
	local aa = subinstr("`var'", "nbirths_at_", "", .)
	local a = real("`aa'")
	bys ANO_NAC bin50g: egen double fert_by_cell_at_a`aa' = mean(cond(dead_by_a15 == 0, `var', .))
	label var fert_by_cell_at_a`aa' "Fertility (mean births) by birth year and weight bin of survivors at age `a'"
}

********************************************************************************
* SELECT SAMPLE AND VARIABLES TO CREATE COUNTERFACTUAL MOTHERS DATABASE:
* Keep only deceased girls:
keep if dead_by_a15 == 1

* Drop individual mortality and fertility variables:
drop dead_at_* dead_by_a15 nbirths_at_??

* Rename some variables:
rename bin50g bin50gm
label var bin50gm "Mother: Birth weight 50 gram bins"

********************************************************************************
* RENAME NAC VARIABLES:
* Rename father variables to grandfather:
foreach var of varlist *_PADRE {
	local varlbl : variable label `var'
	local nvar = subinstr("`var'", "_PADRE", "_ABUELO", .)
	rename `var' `nvar'
	label var `nvar' "Abuelo: `varlbl'"
}

* Rename mother variables to grandmother:
foreach var of varlist *_MADRE {
	local varlbl : variable label `var'
	local nvar = subinstr("`var'", "_MADRE", "_ABUELA", .)
	rename `var' `nvar'
	label var `nvar' "Abuela: `varlbl'"
}

* Rename the rest of the variables as mother variables:
ds *_ABUELO *_ABUELA ID_RECIEN_NACIDO fert_by_cell_at_a?? mort_by_cell bin50gm, not
foreach var of varlist `r(varlist)' {
	local varlbl : variable label `var'
	rename `var' `var'_MADRE
	label var `var'_MADRE "Madre: `varlbl'"
}
********************************************************************************
* EDUCATION MERGE AND RENAME:
* Merge in education variables:
merge 1:1 ID_RECIEN_NACIDO using "$DAT\EDUCATION_VARIABLES.dta", ///
	keep(match) nogen keepusing(ed?lvls_* highest_* ed*o?)

replace edgmom = edmom
replace edmom = .

replace edgpop = edpop
replace edpop = .

replace ed3lvls_aa = ed3lvls_m
replace ed3lvls_m = .
replace ed3lvls_ao = ed3lvls_p
replace ed3lvls_p = .

replace highest_ed3lvls_aaoo = highest_ed3lvls_mp
replace highest_ed3lvls_mp = .

replace ed2lvls_aa = ed2lvls_m
replace ed2lvls_m = .
replace ed2lvls_ao = ed2lvls_p
replace ed2lvls_p = .

replace highest_ed2lvls_aaoo = highest_ed2lvls_mp
replace highest_ed2lvls_mp = .

********************************************************************************
* HEAPING MERGE AND RENAME:
* Merge in heaping variables:
merge 1:1 ID_RECIEN_NACIDO using "$DAT\HEAPING_VARIABLES.dta", ///
	keep(match) nogen keepusing(round* dbw*)
	
foreach var of varlist dbw???? {
    local mvar = "`var'm"
    capture confirm var `mvar'
    if _rc == 0 {
        replace `mvar' = `var'
    }
    else {
        local var_lbl : variable label `var'
        local mvar_lbl  = subinstr("`var_lbl'", "Birth weight", "Mother's Birth weight", .)
        gen `mvar' = `var'
        label var `mvar' "`mvar_lbl'"		
    }
    drop `var'
}

replace round050m = round050
drop round050
replace round100m = round100
drop round100

********************************************************************************
* LBW INDICATORS MERGE AND RENAME:
* Merge in low birth weight indicators:
merge 1:1 ID_RECIEN_NACIDO using "$DAT\LBW_INDICATORS.dta", ///
keep(match) nogen keepusing(bw_below_* vlbw* bw2k* lbw*)

foreach var of varlist bw_below_???? {
    local mvar = "`var'm"
    capture confirm var `mvar'
    if _rc == 0 {
        replace `mvar' = `var'
    }
    else {
        local var_lbl : variable label `var'
        local mvar_lbl  = subinstr("`var_lbl'", "Child's", "Mother's", .)
        gen `mvar' = `var'
        label var `mvar' "`mvar_lbl'"		
    }
    drop `var'
}

replace vlbwm = vlbw
drop vlbw
replace bw2km = bw2k
drop bw2k
replace lbwm = lbw
drop lbw

********************************************************************************
* BIRTH TIMING MERGE AND RENAME:
* Merge in birth timing variables:
merge 1:1 ID_RECIEN_NACIDO using "$DAT\BIRTH_TIMING_VARIABLES.dta", ///
	keep(match) nogen keepusing(sem32* gw_below_* dayofweek* wkndbirth*)
	
replace sem32m = sem32
drop sem32
replace dayofweekm = dayofweek
drop dayofweek
replace wkndbirthm = wkndbirth
drop wkndbirth
replace wkndbirth2m = wkndbirth2
drop wkndbirth2

foreach var of varlist gw_below_?? {
    local mvar = "`var'm"
    capture confirm var `mvar'
    if _rc == 0 {
        replace `mvar' = `var'
    }
    else {
        local var_lbl : variable label `var'
        local mvar_lbl  = subinstr("`var_lbl'", "Child's", "Mother's", .)
        gen `mvar' = `var'
        label var `mvar' "`mvar_lbl'"		
    }
    drop `var'
}

********************************************************************************
* BLN CONTROLS MERGE AND RENAME
* Merge in BLN (2013) control variables:
merge 1:1 ID_RECIEN_NACIDO using "$DAT\BLN2013_CONTROLS.dta", ///
	keep(match) nogen keepusing(married* doc_aten* bregion_* byear_* g?ok_blnctrls)
	
replace marriedm = married
drop married

replace doc_atenm = doc_aten
drop doc_aten

foreach var of varlist bregion_?? {
	local mvar = "`var'm"
	replace `mvar' = `var'
	drop `var'
}

foreach var of varlist byear_???? {
	local mvar = "`var'm"
	capture confirm var `mvar'
	if _rc == 0 {
		replace `mvar' = `var'
	}
	drop `var'
}

* Rename ID variable:
rename ID_RECIEN_NACIDO ID_MADRE
label var ID_MADRE "Identificador único y anónimo de la madre del recién nacido vivo"

********************************************************************************
* SIMULATE SURVIVAL:
* Obtain point estimates of the effect on mortality for gen 1 females by age 15:
estimates use "$OUT/g1f_mortbya15_o32_hps.ster" 
scalar tau_bc_mortbya15_o32_hps = e(tau_bc)

* Create survival probability:
gen surv_prob_o32_hps = max(min(1 - (mort_by_cell + tau_bc_mortbya15_o32_hps), .9999999), 0.000001)

* Create survival dummy:
gen surv_o32_hps = rbinomial(1, surv_prob_o32_hps) if PESO >= 1500 & SEMANAS >= 32
replace surv_o32_hps = 0 if PESO < 1500 & SEMANAS >= 32

* Compress, label, sign, and save:
compress
label data "Counterfactual mothers: girls that died by age 15 with mother variables"
notes drop _all
note: This dataset contains observations of girls that died by age 15, as if they had become mothers.
save "$SEL\COUNTERFACTUAL_MOTHERS.dta", replace

*-------------------------------------------------------------------------------
*--- Expected Births Dead Girls
*-------------------------------------------------------------------------------

use "$SEL\COUNTERFACTUAL_MOTHERS.dta", clear

********************************************************************************
* SIMULATE FERTILITY:
* Obtain point estimates of the effect on fertility for gen 1 females at age:
local o32_list : dir "$OUT" files "g1f_nbirthsat??_o32_hps.ster", respect
foreach est of local o32_list {
    local aa = subinstr(subinstr("`est'", "g1f_nbirthsat", "", .), "_o32_hps.ster", "", .)
    estimates use "$OUT/`est'"
    scalar tau_bc_nbirthsat`aa'_o32_hps = e(tau_bc)
}

* Create variable for number of counterfactual births:
foreach var of varlist fert_by_cell_at_a?? {
    local aa = subinstr("`var'", "fert_by_cell_at_a", "", .)
    local a = real("`aa'")
    local b = `a' + 1
    local bb = string(`b', "%02.0f")
    local z = `a' - 1
    local zz = string(`z', "%02.0f")
    
    * Calculate probability of having a child at a given age:
	gen fert_prob_at_`aa'_o32_hps = max(min(`var' + tau_bc_nbirthsat`aa'_o32_hps, .9999999), 0.000001)
	
	* Determine expected births at a given age:
	gen ebirths_at_`aa'_o32_hps = rbinomial(1, fert_prob_at_`aa'_o32_hps) if surv_o32_hps == 1 & PESO >= 1500 & SEMANAS >= 32 & ANO_NAC + `a' <= 2018
	
	* Correct for recent birth:
	if `a' >= 16 {
		replace ebirths_at_`aa'_o32_hps = 0 if ebirths_at_`zz'_o32_hps == 1
	}
	
	* Calculate running sum of how many expected children by a given age:
	egen ebirths_by_`bb'_o32_hps = rowtotal(ebirths_at_??_o32_hps) if surv_o32_hps == 1 & PESO >= 1500 & SEMANAS >= 32 & ANO_NAC + `a' <= 2018, m
}

* Calculate total number of expected births:
egen total_ebirths_o32_hps = rowtotal(ebirths_at_??_o32_hps), m

* Keep only relevant variables:
keep ID_MADRE total_ebirths_???_hps ebirths_??_??_???_hps

* Compress, label, sign, and save:
compress
label data "Expected births for simulated survivals."
notes drop _all
note: This dataset contains observations of girls that died by age 15, with the number of expected births if they had survived.
save "$SEL\expected_births_dead_girls.dta", replace

*-------------------------------------------------------------------------------
*--- Cell Stats - Bin 50 grams
*-------------------------------------------------------------------------------

* Load data:
use "$DAT\workingdata_age.dta", clear

* Base sample for second generation:
/* 
-> children for whom we have mother's birth data: 
	mrg_mbdata2main == 3
-> children of mothers who were between 15 and 26 years old at time of birth: 
	EDAD_MADRE >= 15 & EDAD_MADRE <= 26
-> children born from 2007 onwards (only 9 births between 2001-2006 given above
 constraints): 
	ANO_NAC >= 2007
*/
gen g2smpl = mrg_mbdata2NAC == 3 ///
	& EDAD_MADRE >= 15 & EDAD_MADRE <= 26 ///
	& ANO_NAC >= 2007

* Tabulate sample by birth year:
tab ANO_NAC g2smpl

* Drop unused observations:
keep if g2smpl == 1

* Construct some pending variables:
gen fgrate = PESO/SEMANAS
label var fgrate "Fetal growth rate (PESO/SEMANAS)"

gen sex_male = SEXO == 1 if SEXO == 1 | SEXO == 2
label var sex_male "Sex is male"

* Create 50 gram bins according to mother's birth weight:
gen int bin50gm = 50*floor(PESO_MADRE / 50)
label var bin50gm "Mother's birth weight 50 gram bins"

* Construct variable of death by age 15:
egen dead_by_a15 = rowtotal(dead_at_a00-dead_at_a14), m
label var dead_by_a15 "Dead before 15th birthday."

* Construct local for collapse:
local dvar_list PESO SEMANAS TALLA dead_at_a00 gw_below_36 vlbw fgrate sex_male sgabw sgasz

* Calculate percentiles 10 to 90 of dependent variables by mother's 50g bins:
gen one = 1
macro drop _filelist
forval p = 5(5)95 {
	preserve
	collapse (p`p') `dvar_list' (sum) cellcount = one, by(bin50gm)
	drop if bin50gm == .
	gen p = `p'
	tempfile p`p'
	compress
	save `p`p'', replace
	local filelist `filelist' `p`p''
	restore
}

clear

append using `filelist'

macro drop _dvar_list_p
foreach var of varlist `dvar_list' {
	rename `var' `var'_p
	local dvar_list_p `dvar_list_p' `var'_p
}

reshape wide `dvar_list_p', i(bin50gm cellcount) j(p)

save "$SEL\cell_stats_bin50gm_5pc.dta", replace

*-------------------------------------------------------------------------------
*--- Counterfactual Babies
*-------------------------------------------------------------------------------

* Load data:
use "$SEL\expected_births_dead_girls.dta", clear

foreach ggg in o32 {
	preserve
	
	* Use only o32 estimates:
	keep *_`ggg'_hps ID_MADRE

	* Keep only girls that are expected to have at least one child:
	keep if total_ebirths_`ggg'_hps > 0 & total_ebirths_`ggg'_hps != .

	foreach var of varlist ebirths_* {
		local t = word(subinstr("`var'", "_", " ", .), 2)
		local aa = word(subinstr("`var'", "_", " ", .), 3)
		rename `var' ebirths_`ggg'_hps_`t'_`aa'
	}
	drop ebirths_`ggg'_hps_by_??
	reshape long ebirths_`ggg'_hps_at_, i(ID_MADRE) j(EDAD_MADRE)
	keep if ebirths_`ggg'_hps_at_ == 1

	* Drop unnecessary variables:
	drop ebirths_`ggg'_hps_at_ total_ebirths_`ggg'_hps

	* Create birth order variable:
	sort ID_MADRE EDAD_MADRE
	by ID_MADRE: egen birth_order_by_mother_t = rank(EDAD_MADRE), track
	
	* Merge in counterfactual mothers' data:
	merge m:1 ID_MADRE using "$SEL\COUNTERFACTUAL_MOTHERS.dta", nogen keep(master match) // should be a perfect match (I checked!)
	
	* Create year of birth of the expected child:
	gen ANO_NAC = ANO_NAC_MADRE + EDAD_MADRE
	
	* Merge in cell stats (non matches are because mother has missing weight):
	merge m:1 bin50gm using "$SEL\cell_stats_bin50gm_5pc.dta", nogen keep(match)

	* Compress and save:
	compress
	save "$SEL\counterfactual_babies_`ggg'_hps.dta", replace
	
	restore
}

*-------------------------------------------------------------------------------
*--- Counterfactual Babies - Fixed
*-------------------------------------------------------------------------------

* Load surviving mothers based on pooled mortality estimates (18,930):
use "$SEL\COUNTERFACTUAL_MOTHERS.dta", clear

foreach ggg in o32 {
    * Loop over inputed number of children
    foreach nbs of numlist 1 2 3 4 5 6 7 8 9 10 11 12 {
        preserve
        * Keep only surviving mothers:
	keep if surv_`ggg'_hps == 1
		
	* Expand dataset by fixed number of children:
        if `nbs'==1|`nbs'==2 {
            expand `nbs'
	}	
        else if `nbs'< 12 {
            local val = (`nbs'-2)/10
            gen keeper = runiform()
            keep if keeper < `val'
            drop keeper
        }
        if `nbs'==12 {
            local val = 0
            gen keeper = runiform()
            keep if keeper < `val'
            drop keeper
        }

        * Merge in cell stats:
 	merge m:1 bin50gm using "$SEL\cell_stats_bin50gm_5pc.dta", nogen keep(match)
	
	* Compress and save:
	compress
	save "$SEL\counterfactual_babies_`ggg'_fixed`nbs'.dta", replace
	restore
    }
}


*log close
