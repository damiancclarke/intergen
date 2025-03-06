/* AppendixReplication.do        KTS/DCC/NLB               yyyy-mm-dd:2025-01-18
----|----1----|----2----|----3----|----4----|----5----|----6----|----7----|----8

This do file replicates the results provided in the paper "Estimating 
Inter­generational Returns to Medical Care: New Evidence from At­Risk Newborns" 
written by Damian Clarke, Nicolas Lillo Bustos and Kathya Tapia-Schythe.  
In certain cases these results will require the user-written, estout, 
reghdfe, rdbwselect, rdplot, rdrobust, binsreg, blindschemes, swindex 
schemepack, rdpow commands.  
*/

*-------------------------------------------------------------------------------
*--- Preamble
*-------------------------------------------------------------------------------

clear all
set more off
capture log close
*set scheme s1mono
graph set window fontface "Times New Roman"
set matsize 11000
set seed 30081985

*-------------------------------------------------------------------------------
*--- Directories and Log 
*-------------------------------------------------------------------------------

global DIR "F:\intergenerational\Replication"
global DO  "$DIR\source"
global DAT "$DIR\data"
global EXP "$DIR\data\Appendix"
cap mkdir $EXP
global LOG "$DIR\log\Appendix"
cap mkdir $LOG
global OUT "$DIR\results\Appendix"
cap mkdir $OUT
global FIG "$DIR\figures\Appendix"
cap mkdir $FIG
global TAB "$DIR\tables\Appendix"
cap mkdir $TAB

* Log file:
log using "$LOG\AppendixReplication.log", replace

*-------------------------------------------------------------------------------
*--- Variables
*-------------------------------------------------------------------------------

global g1blnctrls edmom EDAD_MADRE married doc_aten bregion_?? byear_????
global g2blnctrls edgmom EDAD_ABUELA marriedm doc_atenm bregion_??m byear_????m

*-------------------------------------------------------------------------------
*--- Tables
*-------------------------------------------------------------------------------

*** Table A2: Temporal Links between Mother-Child Matched Birth Years ***
use "$DAT/workingdata_age_work.dta", clear
tab ANO_NAC_MADRE ANO_NAC

*11 Corrections between 2001-2006 
drop if ANO_NAC_MADRE==1999 & ANO_NAC==2001 & EDAD_MADRE==26 
drop if ANO_NAC_MADRE==1999 & ANO_NAC==2002 & EDAD_MADRE==23 
drop if ANO_NAC_MADRE==1999 & ANO_NAC==2002 & EDAD_MADRE==26
drop if ANO_NAC_MADRE==1999 & ANO_NAC==2002 & EDAD_MADRE==27
drop if ANO_NAC_MADRE==1994 & ANO_NAC==2003 & EDAD_MADRE==19 
drop if ANO_NAC_MADRE==1997 & ANO_NAC==2003 & EDAD_MADRE==16 
drop if ANO_NAC_MADRE==1994 & ANO_NAC==2004 & EDAD_MADRE==17 
drop if ANO_NAC_MADRE==1999 & ANO_NAC==2004 & EDAD_MADRE==30 
drop if ANO_NAC_MADRE==1996 & ANO_NAC==2005 & EDAD_MADRE==16 
drop if ANO_NAC_MADRE==1995 & ANO_NAC==2006 & EDAD_MADRE==20 

tab ANO_NAC_MADRE ANO_NAC 


*** Table A3: Summary Statistics - Births local to the 1500 gram threshold ***
use "$DAT/workingdata_age_work.dta", clear

*Labels
lab var SEMANAS     "Gestation Weeks"
lab var sem32       "$\geq$ 32 Gestation Weeks"
lab var PESO        "Birth Weight in Grams"
lab var vlbw        "Birth Weight < 1,500"
lab var TALLA       "Birth Length in cm"
lab var dead_at_a00 "Death Within 1st Year of Birth"   /*check observations*/
lab var days_y00    "Days Spent in Hospital by Year 1" 
lab var nadmssn_y00 "Number of Admissions to Hospital by Year 1"
lab var EDAD_MADRE  "Mother's Age"
lab var edmom       "Mother's Education Years"

* Panel A: Bandwidth 134.4 IMR Gen1 all sample
preserve 
*Muestra madres
keep if EDAD_MADRE>=15 & EDAD_MADRE<=45

keep if PESO<=1500+134.4
keep if PESO>=1500-134.4

#delimit ;
estpost sum SEMANAS sem32 PESO vlbw TALLA dead_at_a00 days_y00 nadmssn_y00 EDAD_MADRE edmom;
estout using "$TAB/SummaryB_G1.tex", replace label style(tex)
cells("count(fmt(%15.0gc)) mean(fmt(2)) sd(fmt(2)) min(fmt(2)) max(fmt(2))")
collabels(, none) mlabels(, none);
#delimit cr
restore

* Panel B Bandwidth 245.5 Anderson index Gen2
preserve
gen	g2smpl = (mrg_mbdata2main==3 & EDAD_MADRE>=15 & EDAD_MADRE<=45 & ANO_NAC>=2007)
keep if g2smpl==1

keep if PESO_MADRE<=1500+245.5
keep if PESO_MADRE>=1500-245.5

#delimit ;
estpost sum SEMANAS sem32 PESO vlbw TALLA dead_at_a00 EDAD_MADRE edmom;
estout using "$TAB/SummaryB_G2.tex", replace label style(tex)
cells("count(fmt(%15.0gc)) mean(fmt(2)) sd(fmt(2)) min(fmt(2)) max(fmt(2))")
collabels(, none) mlabels(, none);
#delimit cr
restore
estimates clear


*** Table B1: Initial Policy Impacts on Infant Mortality Rates ***
use "$DAT/workingdata_age_work.dta", clear

keep if SEMANAS>=32& SEMANAS!=.

foreach year of numlist 2007 2018 2001 {
 preserve
 keep if ANO_NAC >= 1992 & ANO_NAC<=`year' 
         *& g2ok_blnctrls == 1 & g2smpl = mrg_mbdata2NAC == 3 ///
         *& EDAD_MADRE >= 15 & EDAD_MADRE <= 45
 
 * Odd Columns
 eststo: rdrobust dead_at_a00 PESO, c(1500) scalepar(-1) all covs(${g1blnctrls} dbw1??0) vce(cluster PESO) h(100)
 estadd scalar Hopt = e(h_l)
 estadd scalar Nl   = e(N_h_l)
 estadd scalar Nr   = e(N_h_r)
 local ef = e(N_h_l) + e(N_h_r)
 estadd scalar effopt = `ef'
 sum dead_at_a00 if PESO>=1400 & PESO<=1600
 estadd scalar dvmean = r(mean)
 
 * Even Columns
 rdbwselect dead_at_a00 PESO, c(1500) covs(${g1blnctrls} dbw1??0) vce(cluster PESO)
 local RBW = e(h_mserd)
 local LBW = e(h_mserd)
 local maxwt = 1500+`RBW'
 local minwt = 1500-`LBW'
	
 eststo: rdrobust dead_at_a00 PESO, c(1500) scalepar(-1) all covs(${g1blnctrls} dbw1??0) vce(cluster PESO)
 estadd scalar Hopt = e(h_l)
 estadd scalar Nl   = e(N_h_l)
 estadd scalar Nr   = e(N_h_r)
 local ef = e(N_h_l) + e(N_h_r)
 estadd scalar effopt = `ef'
 sum dead_at_a00 if PESO>=`minwt' & PESO<=`maxwt'
 estadd scalar dvmean = r(mean)
    
 restore
}

gen Robust =.
lab var Robust "Robust"
gen Conventional =.
lab var Conventional "Conventional"
#delimit ;
esttab est1 est2 est3 est4 est5 est6 using "$TAB\tableIMR.tex",
stats(dvmean N h_r effopt Nl Nr, fmt(%8.3f %15.0gc %5.1f %8.0gc %8.0gc %8.0gc)
labels("\\ Mean of Dep. Var." "Observations" "Optimal Bandwidth"
       "Effective Observations" "Observations (left)"
       "Observations (right)")) keep(Conventional Robust) label
b(%-9.4f) se(%-9.4f) noobs nonotes nogaps mlabels(, none) nonumbers style(tex)
fragment replace noline starlevel("*" 0.10 "**" 0.05 "***" 0.01);
#delimit cr
estimates clear


*** Table B2: Balance Tests - Generation 1 Family and Birth Characteristics ***
use "$DAT/workingdata_age_work.dta", clear

keep if EDAD_MADRE >= 15 & EDAD_MADRE <= 45 
drop if PESO_MADRE != . // Drop Gen2

gen marrieda = EST_CIV_MADRE==2 

gen ageMum  = EDAD_MADRE
gen educMum = CURSO_MADRE if NIVEL_MADRE==4 | NIVEL_MADRE==5 
replace educMum = CURSO_MADRE+8 if NIVEL_MADRE==2 | NIVEL_MADRE==3
replace educMum = CURSO_MADRE+12 if NIVEL_MADRE==1
gen employedMum = activ_m==1

gen ageDad  = EDAD_PADRE
gen educDad = CURSO_PADRE if NIVEL_PADRE==4 | NIVEL_PADRE==5 
replace educDad = CURSO_PADRE+8 if NIVEL_PADRE==2 | NIVEL_PADRE==3
replace educDad = CURSO_PADRE+12 if NIVEL_PADRE==1
gen employedDad = activ_p==1

gen urbano= URBANO_RURAL==1

gen difedad = abs(ageDad-ageMum)

gen obs_dad = (EDAD_PADRE!=.)

foreach var of varlist ageMum educMum employedMum ageDad educDad employedDad marrieda difedad obs_dad urbano* {
 keep if SEMANAS>=32& SEMANAS!=.
	
 eststo: rdrobust `var' PESO, c(1500) scalepar(-1) vce(cluster PESO) all 
 local ef = e(N_h_l)+e(N_h_r)
 estadd scalar effopt = `ef'
 local BW = e(h_l)
 local maxwt = 1500+`BW'
 local minwt = 1500-`BW'
 sum `var' if PESO >=`minwt' & PESO<=`maxwt'
 estadd scalar dvmean = r(mean)
	
 estimates save "$OUT/g1_`var'_o32.ster", replace	
}

clear 
gen Robust =.
lab var Robust "Birth weight $<$ 1,500"
gen po32 = .
set obs 10
local i = 1
foreach outcome in ageMum educMum employedMum ageDad educDad employedDad marrieda difedad obs_dad urbano { 
 estimates use "$OUT/g1_`outcome'_o32.ster"
 eststo est`i'
 replace po32=e(pv_rb) in `i'
 local ++i
}
rename po32 pval
qsharpenedp pval

local l = 1
foreach outcome in ageMum educMum employedMum ageDad educDad employedDad marrieda difedad obs_dad urbano {
    sum bky06 in `l'
    estadd scalar qpv = r(mean): est`l'
    local ++l
}

local ests est1 est2 est3 est4 est5 est6 est7 est8 est9 est10
#delimit ;
esttab `ests' using "$TAB/tablefigureA3.tex",
replace booktabs cells(b(fmt(%-9.3f) star) se(fmt(%-9.3f) par(( )))) label
stats(dvmean N h_r effopt N_h_l N_h_r qpv, 
      fmt(%05.3f %12.0gc %5.1f %9.0gc %9.0gc %9.0gc %05.3f) 
      label("\\ Mean of Dep. Var." "Observations" "Optimal Bandwidth" 
            "Effective Observations" "Observations (left)" "Observations (right)"
			"q-sharpened p-value"))
nonotes nogaps mlabels(, none) nonumbers style(tex) fragment noline keep(Robust) varlabels(Robust "Birth weight $<$ 1,500")
collabels(none) starlevel("*" 0.1 "**" 0.05 "***" 0.01);
#delimit cr
estimates clear


*** Table B3: Testing for birth heaping at the cut-off ***
use "$DAT/workingdata_age_work.dta", clear

keep if PESO>=1300&PESO<=1700
sum PESO

gen A = 1
collapse TALLA SEMANAS (count) nt=TALLA ns=SEMANAS (sum) A, by(PESO)

local sigstar starlevel("*" 0.1 "**" 0.05 "***" 0.01)
local nstar star
local ests est1 est2 est3 est4 est5 est6 est7 est8 est9 est10
lab var ns "Number of Births"
gen below = PESO<1500 if PESO!=1500
lab var below "$<$ 1,500 gram threshold"

* Panel A:
gen absDif = abs(PESO-1500)
foreach num of numlist 10 20 30 40 50 75 100 125 150 200 {
    eststo: reg A below if absDif<=`num'
    sum A if e(sample)==1
    local TOT = r(N)*r(mean)
    estadd scalar total = `TOT'    
}

#delimit ;
esttab `ests' using "$TAB/AlmondTest.tex",
replace booktabs cells(b(fmt(%-9.3f) `nstar') se(fmt(%-9.3f) par(( )) )) label
stats(total N r2, fmt(%9.0gc %9.0gc %05.3f)
      label("\\ Total Births" "Gram-Specific cells" R-Squared))
nonotes nogaps mlabels(, none) nonumbers style(tex) fragment noline keep(below)
collabels(none) `sigstar';
#delimit cr
estimates clear

* Panel B:
gen absDifBelow = absDif*below
lab var absDif "Absolute Difference (grams)"
lab var absDifBelow "Absolute Difference (grams) $\times <$ 1,500 grams"
foreach num of numlist 10 20 30 40 50 75 100 125 150 200 {
    eststo: reg A below absDif absDifBelow if absDif<=`num'
    sum A if e(sample)==1
    local TOT = r(N)*r(mean)
    estadd scalar total = `TOT'  
	estadd scalar am    = `num'
}
#delimit ;
esttab `ests' using "$TAB/AlmondTest_trend.tex",
replace booktabs cells(b(fmt(%-9.3f) `nstar') se(fmt(%-9.3f) par(( )))) label
stats(total N r2 am, fmt(%9.0gc %9.0gc %05.3f %9.0gc )
      label("\\ Total Births" "Gram-Specific cells" "R-Squared" "Absolute margin (grams)"))
nonotes nogaps mlabels(, none) nonumbers style(tex) fragment noline keep(below absDif absDifBelow)
collabels(none) `sigstar';
#delimit cr
estimates clear


*** Table D1 and D2: Alternative Specifications ***
use "$DAT/workingdata_age_work.dta", clear
gen g2smpl = mrg_mbdata2main == 3 ///
	& EDAD_MADRE >= 15 & EDAD_MADRE <= 45 & ANO_NAC >= 2007 & g2ok_blnctrls == 1

keep if g2smpl == 1

gen premature = SEMANAS<37
gen fgrate = PESO/SEMANAS

gen not_vlbw = 1 - vlbw
gen not_dead_at_a00 = 1 - dead_at_a00
gen not_gw_below_36 = 1 - gw_below_36
gen control_group = PESO_MADRE >= 1500 if PESO_MADRE != .
swindex PESO SEMANAS TALLA not_dead_at_a00 not_vlbw not_gw_below_36 fgrate, normby(control_group) gen(aindex)

foreach outcome in SEMANAS PESO TALLA dead_at_a00 premature vlbw fgrate aindex {

if ("`outcome'"=="SEMANAS")     local lab "gwks"
if ("`outcome'"=="PESO")        local lab "peso"
if ("`outcome'"=="TALLA")       local lab "talla"
if ("`outcome'"=="dead_at_a00") local lab "deadata00"
if ("`outcome'"=="premature")   local lab "prem36"
if ("`outcome'"=="vlbw")        local lab "vlbw"
if ("`outcome'"=="fgrate")      local lab "fgrate"
if ("`outcome'"=="aindex")      local lab "aindex"

local j = 1
foreach p of numlist 1 2 {

eststo: rdrobust `outcome' PESO_MADRE if sem32m == 1, c(1500) scalepar(-1) all ///
          covs(${g2blnctrls} dbw1??0m) vce(cluster PESO_MADRE) p(`p') 
estadd scalar esti = _b[Robust]
local se = string(_se[Robust], "%05.3f")
estadd local se "(`se')": est`j'
local ++j

eststo: rdrobust `outcome' PESO_MADRE if sem32m == 1, c(1500) scalepar(-1) all ///
          covs(${g2blnctrls} dbw1??0m) vce(cluster PESO_MADRE) bwselect(cerrd) p(`p') 
estadd scalar esti = _b[Robust]
local se = string(_se[Robust], "%05.3f")
estadd local se "(`se')": est`j'
local ++j

eststo: rdrobust `outcome' PESO_MADRE if sem32m == 1, c(1500) scalepar(-1) all ///
          covs(${g2blnctrls} dbw1??0m) vce(cluster PESO_MADRE) p(`p') 
estadd scalar esti = _b[Conventional]
local se = string(_se[Robust], "%05.3f")
estadd local se "(`se')": est`j'
local ++j

eststo: rdrobust `outcome' PESO_MADRE if sem32m == 1, c(1500) scalepar(-1) all ///
          covs(${g2blnctrls} dbw1??0m) vce(cluster PESO_MADRE) p(`p')
estadd scalar esti = _b[Bias-corrected]
local se = string(_se[Robust], "%05.3f")
estadd local se "(`se')": est`j'
local ++j
}

* Panel A-H:		  
#delimit ;
esttab est1 est2 est3 est4 est5 est6 est7 est8 using "$TAB\_`lab'.tex",
drop(Conventional Bias-corrected Robust) b(%-9.3f) se(%-9.3f)
stats(esti se N h_r, fmt(%-9.3f %05.3f %9.0gc %5.1f) 
labels("Birth weight $<$1500" " " "Observations" "Optimal Bandwidth") 
star(esti)) noobs nonotes nogaps mlabels(, none) nonumbers style(tex) 
fragment replace noline starlevel("*" 0.10 "**" 0.05 "***" 0.01);
estimates clear;
#delimit cr
}
estimates clear


*** Table D3: Identificacion Considerations – Donut RD and Placebo Treatments ***

* Panel A:
use "$DAT/workingdata_age_work.dta", clear
gen g2smpl = mrg_mbdata2main == 3 ///
	& EDAD_MADRE >= 15 & EDAD_MADRE <= 45 & ANO_NAC >= 2007 & g2ok_blnctrls == 1

keep if g2smpl == 1

gen premature = SEMANAS<37
gen fgrate = PESO/SEMANAS

foreach var in SEMANAS PESO TALLA dead_at_a00 premature vlbw fgrate {

 if "`var'"=="SEMANAS" local lab "gwks"
 if "`var'"=="PESO" local lab "peso"
 if "`var'"=="TALLA" local lab "talla"
 if "`var'"=="dead_at_a00" local lab "deadata00"
 if "`var'"=="premature" local lab "prem36"
 if "`var'"=="vlbw" local lab "vlbw"
 if "`var'"=="fgrate" local lab "fgrate"

 rdrobust `var' PESO_MADRE if sem32m == 1, c(1500) scalepar(-1) all ///
          covs(${g2blnctrls} dbw1??0m) vce(cluster PESO_MADRE)

 estimates store `lab'_o32
}

#delimit ;
esttab gwks_o32 peso_o32 talla_o32 deadata00_o32 prem36_o32 vlbw_o32  
       fgrate_o32 using "$TAB/T6A_o32_hps.tex",
replace keep(Robust) cells(b(fmt(%-9.3f) star) se(fmt(%-9.3f) par(( )) ))
starlevel("*" 0.1 "**" 0.05 "***" 0.01) tex nomtitles nonumbers fragment nodepvars
coeflabels(Robust "Birth weight < 1,500") collabels(none) nolines noobs;
#delimit cr

* Panel B: 
use "$DAT/workingdata_age_work.dta", clear
gen g2smpl = mrg_mbdata2main == 3 ///
	& EDAD_MADRE >= 15 & EDAD_MADRE <= 45 & ANO_NAC >= 2007 & g2ok_blnctrls == 1

keep if g2smpl == 1

gen premature = SEMANAS<37
gen fgrate = PESO/SEMANAS

foreach d in 0 5 10 20 {
 if `d'==0  local dlab "00"
 if `d'==5  local dlab "05"
 if `d'==10 local dlab "10"
 if `d'==20 local dlab "20"
 
 keep if abs(PESO_MADRE - 1500) >= `d'
 
 foreach var in SEMANAS PESO TALLA dead_at_a00 premature vlbw fgrate {
  if "`var'"=="SEMANAS" local lab "gwks"
  if "`var'"=="PESO" local lab "peso"
  if "`var'"=="TALLA" local lab "talla"
  if "`var'"=="dead_at_a00" local lab "deadata00"
  if "`var'"=="premature" local lab "prem36"
  if "`var'"=="vlbw" local lab "vlbw"
  if "`var'"=="fgrate" local lab "fgrate"
  
  rdrobust `var' PESO_MADRE if sem32m == 1, c(1500) scalepar(-1) all ///
            covs(${g2blnctrls} dbw1??0m) vce(cluster PESO_MADRE)

  estimates store `lab'`d'_o32
 }

 #delimit ;
 esttab gwks`d'_o32 peso`d'_o32 talla`d'_o32 deadata00`d'_o32 prem36`d'_o32 
        vlbw`d'_o32 fgrate`d'_o32  using "$TAB/T6B`dlab'_o32_hps.tex",
 replace keep(Robust) cells(b(fmt(%-9.3f) star) se(fmt(%-9.3f) par(( )) ))
 starlevel("*" 0.1 "**" 0.05 "***" 0.01) tex nomtitles nonumbers fragment nodepvars
 coeflabels(Robust "Donut = `d'") collabels(none) nolines noobs;
 #delimit cr
}

* Panel C:
use "$DAT/workingdata_age_work.dta", clear
gen g2smpl = mrg_mbdata2main == 3 ///
	& EDAD_MADRE >= 15 & EDAD_MADRE <= 45 & ANO_NAC >= 2007 & g2ok_blnctrls == 1

keep if g2smpl == 1

gen premature = SEMANAS<37
gen fgrate = PESO/SEMANAS

foreach c in 1250 1750 2000 2250 2500 {
 if `c'==1250 local clab "1,250"
 if `c'==1750 local clab "1,750"
 if `c'==2000 local clab "2,000"
 if `c'==2250 local clab "2,250"
 if `c'==2500 local clab "2,500"
  
 foreach var in SEMANAS PESO TALLA dead_at_a00 premature vlbw fgrate {
  if "`var'"=="SEMANAS" local lab "gwks"
  if "`var'"=="PESO" local lab "peso"
  if "`var'"=="TALLA" local lab "talla"
  if "`var'"=="dead_at_a00" local lab "deadata00"
  if "`var'"=="premature" local lab "prem36"
  if "`var'"=="vlbw" local lab "vlbw"
  if "`var'"=="fgrate" local lab "fgrate"
  
  rdrobust `var' PESO_MADRE if sem32m == 1, c(`c') scalepar(-1) all ///
            covs(${g2blnctrls} dbw1??0m) vce(cluster PESO_MADRE)

  estimates store `lab'`c'_o32
 }

 #delimit ;
 esttab gwks`c'_o32 peso`c'_o32 talla`c'_o32 deadata00`c'_o32 prem36`c'_o32 
        vlbw`c'_o32 fgrate`c'_o32 using "$TAB/T6C`c'_o32_hps.tex",
 replace keep(Robust) cells(b(fmt(%-9.3f) star) se(fmt(%-9.3f) par(( )) ))
 starlevel("*" 0.1 "**" 0.05 "***" 0.01) tex nomtitles nonumbers fragment nodepvars
 coeflabels(Robust "Birth weight < `clab'") collabels(none) nolines noobs;
 #delimit cr
}
estimates clear


*** Table D4: Intensive Health Investments and Birth Outcomes of the Second Generation (All births) ***
use "$DAT/workingdata_age_work.dta", clear
gen g2smpl = mrg_mbdata2main == 3 ///
	& EDAD_MADRE >= 15 & EDAD_MADRE <= 45 & ANO_NAC >= 2007 & g2ok_blnctrls == 1

keep if g2smpl == 1

gen premature = SEMANAS<37
gen fgrate = PESO/SEMANAS

gen not_vlbw = 1 - vlbw
gen not_dead_at_a00 = 1 - dead_at_a00
gen not_gw_below_36 = 1 - gw_below_36
gen control_group = PESO_MADRE >= 1500 if PESO_MADRE != .
swindex PESO SEMANAS TALLA not_dead_at_a00 not_vlbw not_gw_below_36 fgrate, normby(control_group) gen(aindex)

foreach outcome in SEMANAS PESO TALLA dead_at_a00 premature vlbw fgrate aindex {

if "`outcome'"=="SEMANAS"     local lab "gwks"
if "`outcome'"=="PESO"        local lab "peso"
if "`outcome'"=="TALLA"       local lab "talla"
if "`outcome'"=="dead_at_a00" local lab "imr_ata00"
if "`outcome'"=="premature"   local lab "gw36"
if "`outcome'"=="vlbw"        local lab "bw1500"
if "`outcome'"=="fgrate"      local lab "fgrate"
if "`outcome'"=="aindex"      local lab "aindex"

rdbwselect `outcome' PESO_MADRE, c(1500) covs(${g2blnctrls} dbw1??0m) vce(cluster PESO_MADRE)
scalar RBW = e(h_mserd)
scalar LBW = e(h_mserd)
scalar maxwt = 1500 + RBW
scalar minwt = 1500 - LBW
gen g2_`outcome'_obw = PESO_MADRE >= minwt & PESO_MADRE <= maxwt

rdrobust `outcome' PESO_MADRE, c(1500) scalepar(-1) all ///
          covs(${g2blnctrls} dbw1??0m) vce(cluster PESO_MADRE)
		  
eststo `outcome'_all
 if ("`lab'"=="imr_ata00") local pv = normal(_b[Robust]/_se[Robust])
 if ("`lab'"=="gw36")      local pv = normal(_b[Robust]/_se[Robust])
 if ("`lab'"=="bw1500")    local pv = normal(_b[Robust]/_se[Robust])
 if ("`lab'"=="gwks")      local pv = 1-normal(_b[Robust]/_se[Robust])
 if ("`lab'"=="peso")      local pv = 1-normal(_b[Robust]/_se[Robust])
 if ("`lab'"=="talla")     local pv = 1-normal(_b[Robust]/_se[Robust])
 if ("`lab'"=="fgrate")    local pv = 1-normal(_b[Robust]/_se[Robust])
 
 local p_value = string(`pv', "%05.3f")
 estadd local p_value "[`p_value']": `outcome'_all
  
estadd scalar Hopt = e(h_l)
estadd scalar Nl   = e(N_h_l)
estadd scalar Nr   = e(N_h_r)
local ef = e(N_h_l) + e(N_h_r)
estadd scalar effopt = `ef'
sum `outcome' if g2_`outcome'_obw == 1 & `outcome' !=. & PESO_MADRE !=.
estadd scalar dvmean = r(mean)

estimates restore `outcome'_all
estimates save "$OUT/g2_`lab'_blnctrls_all_hps.ster", replace
}

clear 
local end _blnctrls
gen Robust =.
lab var Robust "Birth weight $<$ 1,500"
gen pall = .
set obs 7
local i = 1
foreach outcome in gwks peso talla imr_ata00 gw36 bw1500 fgrate { 
 estimates use "$OUT/g2_`outcome'`end'_all_hps.ster"
 eststo `outcome'_all
 replace pall=e(pv_rb) in `i'
 local ++i
}
rename pall pval
qsharpenedp pval

local l = 1
foreach outcome in gwks peso talla imr_ata00 gw36 bw1500 fgrate {
    sum bky06 in `l'
    estadd scalar qpv = r(mean): `outcome'_all
    local ++l
}

#delimit ;
esttab gwks_all peso_all talla_all imr_ata00_all using "$TAB\T2A_all_hps.tex",
stats(p_value dvmean N h_r effopt Nl Nr qpv,
fmt(%05.3f %8.3f %9.0gc %5.1f %8.0gc %8.0gc %8.0gc %05.3f)
labels(" " "\\ Mean of Dep.\ Var." "Observations" "Optimal Bandwidth"
       "Effective Observations" "Observations (left)"
       "Observations (right)" "q-sharpened p-value")) keep("Robust") label
b(%-9.3f) se(%-9.3f) noobs nonotes nogaps mlabels(, none) nonumbers style(tex)
fragment replace noline starlevel("*" 0.10 "**" 0.05 "***" 0.01);
    
esttab gw36_all bw1500_all fgrate_all aindex_all using "$TAB\T2B_all_hps.tex",
stats(p_value dvmean N h_r effopt Nl Nr qpv,
fmt(%05.3f %8.3f %9.0gc %5.1f %8.0gc %8.0gc %8.0gc %05.3f)
labels(" " "\\ Mean of Dep.\ Var." "Observations" "Optimal Bandwidth"
       "Effective Observations" "Observations (left)"
       "Observations (right)" "q-sharpened p-value")) keep("Robust") label
b(%-9.3f) se(%-9.3f) noobs nonotes nogaps mlabels(, none) nonumbers style(tex)
fragment replace noline starlevel("*" 0.10 "**" 0.05 "***" 0.01);
#delimit cr
estimates clear


*** Table D5: Effects Conditioning on Education and Partnership Choices ***

* Panel A
use "$DAT/workingdata_age_work.dta", clear
gen g2smpl = mrg_mbdata2main == 3 ///
	& EDAD_MADRE >= 15 & EDAD_MADRE <= 45 & ANO_NAC >= 2007 & g2ok_blnctrls == 1

keep if g2smpl == 1

gen obs_dad = (EDAD_PADRE!=.)

gen premature = SEMANAS<37
gen fgrate = PESO/SEMANAS

gen not_vlbw = 1 - vlbw
gen not_dead_at_a00 = 1 - dead_at_a00
gen not_gw_below_36 = 1 - gw_below_36
gen control_group = PESO_MADRE >= 1500 if PESO_MADRE != .
swindex PESO SEMANAS TALLA not_dead_at_a00 not_vlbw not_gw_below_36 fgrate, normby(control_group) gen(aindex)

foreach outcome in SEMANAS PESO TALLA dead_at_a00 premature vlbw fgrate aindex {

if "`outcome'"=="SEMANAS"     local lab "gwks"
if "`outcome'"=="PESO"        local lab "peso"
if "`outcome'"=="TALLA"       local lab "talla"
if "`outcome'"=="dead_at_a00" local lab "imr_ata00"
if "`outcome'"=="premature"   local lab "gw36"
if "`outcome'"=="vlbw"        local lab "bw1500"
if "`outcome'"=="fgrate"      local lab "fgrate"
if "`outcome'"=="aindex"      local lab "aindex"

rdbwselect `outcome' PESO_MADRE if sem32m == 1, c(1500) ///
            covs(${g2blnctrls} dbw1??0m edmom) vce(cluster PESO_MADRE)
scalar RBW_o32 = e(h_mserd)
scalar LBW_o32 = e(h_mserd)
scalar maxwt_o32 = 1500 + RBW_o32
scalar minwt_o32 = 1500 - LBW_o32
gen g2_`outcome'_obw_o32 = PESO_MADRE >= minwt_o32 & PESO_MADRE <= maxwt_o32

rdrobust `outcome' PESO_MADRE if sem32m == 1, c(1500) scalepar(-1) all ///
          covs(${g2blnctrls} dbw1??0m edmom) vce(cluster PESO_MADRE)
		  
eststo `outcome'o32
  
estadd scalar Hopt = e(h_l)
estadd scalar Nl   = e(N_h_l)
estadd scalar Nr   = e(N_h_r)
local ef = e(N_h_l) + e(N_h_r)
estadd scalar effopt = `ef'
sum `outcome' if sem32m == 1 & g2_`outcome'_obw_o32 == 1 & `outcome' !=. & PESO_MADRE !=.
estadd scalar dvmean = r(mean)

estimates restore `outcome'o32
estimates save "$OUT/g2_`lab'_blnctrls_o32_hps_controls.ster", replace
}

clear 
local end _blnctrls
gen Robust =.
lab var Robust "Birth weight $<$ 1,500"
gen po32 = .
set obs 7
local i = 1
foreach outcome in gwks peso talla imr_ata00 gw36 bw1500 fgrate { 
    estimates use "$OUT/g2_`outcome'`end'_o32_hps_controls.ster"
    eststo `outcome'o32
    local ++i
}

local ests 
#delimit ;
esttab gwkso32 pesoo32 tallao32 imr_ata00o32 
        gw36o32 bw1500o32 fgrateo32 aindexo32
using "$TAB/T2A_o32_controls.tex", 
stats(dvmean N h_r, fmt(%8.3f %9.0gc %5.1f)
labels("\\ Mean of Dep.\ Var." "Observations" "Optimal Bandwidth")) keep("Robust") label
b(%-9.3f) se(%-9.3f) noobs nonotes nogaps mlabels(, none) nonumbers style(tex)
fragment replace noline starlevel("*" 0.10 "**" 0.05 "***" 0.01);
#delimit cr
estimates clear

* Panel B
use "$DAT/workingdata_age_work.dta", clear
gen g2smpl = mrg_mbdata2main == 3 ///
	& EDAD_MADRE >= 15 & EDAD_MADRE <= 45 & ANO_NAC >= 2007 & g2ok_blnctrls == 1

keep if g2smpl == 1

sum edgmom, d
keep if edgmom>=r(p50)

gen premature = SEMANAS<37
gen fgrate = PESO/SEMANAS

gen not_vlbw = 1 - vlbw
gen not_dead_at_a00 = 1 - dead_at_a00
gen not_gw_below_36 = 1 - gw_below_36
gen control_group = PESO_MADRE >= 1500 if PESO_MADRE != .
swindex PESO SEMANAS TALLA not_dead_at_a00 not_vlbw not_gw_below_36 fgrate, normby(control_group) gen(aindex)

foreach outcome in SEMANAS PESO TALLA dead_at_a00 premature vlbw fgrate aindex {

if "`outcome'"=="SEMANAS"     local lab "gwks"
if "`outcome'"=="PESO"        local lab "peso"
if "`outcome'"=="TALLA"       local lab "talla"
if "`outcome'"=="dead_at_a00" local lab "imr_ata00"
if "`outcome'"=="premature"   local lab "gw36"
if "`outcome'"=="vlbw"        local lab "bw1500"
if "`outcome'"=="fgrate"      local lab "fgrate"
if "`outcome'"=="aindex"      local lab "aindex"

rdbwselect `outcome' PESO_MADRE if sem32m == 1, c(1500) ///
            covs(${g2blnctrls} dbw1??0m edmom) vce(cluster PESO_MADRE)
scalar RBW_o32 = e(h_mserd)
scalar LBW_o32 = e(h_mserd)
scalar maxwt_o32 = 1500 + RBW_o32
scalar minwt_o32 = 1500 - LBW_o32
gen g2_`outcome'_obw_o32 = PESO_MADRE >= minwt_o32 & PESO_MADRE <= maxwt_o32

rdrobust `outcome' PESO_MADRE if sem32m == 1, c(1500) scalepar(-1) all ///
          covs(${g2blnctrls} dbw1??0m edmom) vce(cluster PESO_MADRE)
		  
eststo `outcome'o32
  
estadd scalar Hopt = e(h_l)
estadd scalar Nl   = e(N_h_l)
estadd scalar Nr   = e(N_h_r)
local ef = e(N_h_l) + e(N_h_r)
estadd scalar effopt = `ef'
sum `outcome' if sem32m == 1 & g2_`outcome'_obw_o32 == 1 & `outcome' !=. & PESO_MADRE !=.
estadd scalar dvmean = r(mean)

estimates restore `outcome'o32
estimates save "$OUT/g2_`lab'_blnctrls_o32_hps_highEd.ster", replace
}

clear 
local end _blnctrls
gen Robust =.
lab var Robust "Birth weight $<$ 1,500"
gen po32 = .
set obs 7
local i = 1
foreach outcome in gwks peso talla imr_ata00 gw36 bw1500 fgrate { 
    estimates use "$OUT/g2_`outcome'`end'_o32_hps_highEd.ster"
    eststo `outcome'o32
    local ++i
}

local ests 
#delimit ;
esttab gwkso32 pesoo32 tallao32 imr_ata00o32 
        gw36o32 bw1500o32 fgrateo32 aindexo32
using "$TAB/T2A_o32_highEd.tex", 
stats(dvmean N h_r, fmt(%8.3f %9.0gc %5.1f)
labels("\\ Mean of Dep.\ Var." "Observations" "Optimal Bandwidth")) keep("Robust") label
b(%-9.3f) se(%-9.3f) noobs nonotes nogaps mlabels(, none) nonumbers style(tex)
fragment replace noline starlevel("*" 0.10 "**" 0.05 "***" 0.01);
#delimit cr
estimates clear

* Panel C
use "$DAT/workingdata_age_work.dta", clear
gen g2smpl = mrg_mbdata2main == 3 ///
	& EDAD_MADRE >= 15 & EDAD_MADRE <= 45 & ANO_NAC >= 2007 & g2ok_blnctrls == 1

keep if g2smpl == 1

gen obs_dad = (EDAD_PADRE!=.)

gen premature = SEMANAS<37
gen fgrate = PESO/SEMANAS

gen not_vlbw = 1 - vlbw
gen not_dead_at_a00 = 1 - dead_at_a00
gen not_gw_below_36 = 1 - gw_below_36
gen control_group = PESO_MADRE >= 1500 if PESO_MADRE != .
swindex PESO SEMANAS TALLA not_dead_at_a00 not_vlbw not_gw_below_36 fgrate, normby(control_group) gen(aindex)

foreach outcome in SEMANAS PESO TALLA dead_at_a00 premature vlbw fgrate aindex {

if "`outcome'"=="SEMANAS"     local lab "gwks"
if "`outcome'"=="PESO"        local lab "peso"
if "`outcome'"=="TALLA"       local lab "talla"
if "`outcome'"=="dead_at_a00" local lab "imr_ata00"
if "`outcome'"=="premature"   local lab "gw36"
if "`outcome'"=="vlbw"        local lab "bw1500"
if "`outcome'"=="fgrate"      local lab "fgrate"
if "`outcome'"=="aindex"      local lab "aindex"

rdbwselect `outcome' PESO_MADRE if sem32m == 1, c(1500) ///
            covs(${g2blnctrls} dbw1??0m obs_dad) vce(cluster PESO_MADRE)
scalar RBW_o32 = e(h_mserd)
scalar LBW_o32 = e(h_mserd)
scalar maxwt_o32 = 1500 + RBW_o32
scalar minwt_o32 = 1500 - LBW_o32
gen g2_`outcome'_obw_o32 = PESO_MADRE >= minwt_o32 & PESO_MADRE <= maxwt_o32

rdrobust `outcome' PESO_MADRE if sem32m == 1, c(1500) scalepar(-1) all ///
          covs(${g2blnctrls} dbw1??0m obs_dad) vce(cluster PESO_MADRE)
		  
eststo `outcome'o32
  
estadd scalar Hopt = e(h_l)
estadd scalar Nl   = e(N_h_l)
estadd scalar Nr   = e(N_h_r)
local ef = e(N_h_l) + e(N_h_r)
estadd scalar effopt = `ef'
sum `outcome' if sem32m == 1 & g2_`outcome'_obw_o32 == 1 & `outcome' !=. & PESO_MADRE !=.
estadd scalar dvmean = r(mean)

estimates restore `outcome'o32
estimates save "$OUT/g2_`lab'_blnctrls_o32_hps_control_relationship.ster", replace
}

clear 
local end _blnctrls
gen Robust =.
lab var Robust "Birth weight $<$ 1,500"
gen po32 = .
set obs 7
local i = 1
foreach outcome in gwks peso talla imr_ata00 gw36 bw1500 fgrate { 
    estimates use "$OUT/g2_`outcome'`end'_o32_hps_control_relationship.ster"
    eststo `outcome'o32
    local ++i
}

local ests 
#delimit ;
esttab gwkso32 pesoo32 tallao32 imr_ata00o32 
       gw36o32 bw1500o32 fgrateo32 aindexo32
using "$TAB/T2A_o32_control_relationship.tex", 
stats(dvmean N h_r, fmt(%8.3f %9.0gc %5.1f)
labels("\\ Mean of Dep.\ Var." "Observations" "Optimal Bandwidth")) keep("Robust") label
b(%-9.3f) se(%-9.3f) noobs nonotes nogaps mlabels(, none) nonumbers style(tex)
fragment replace noline starlevel("*" 0.10 "**" 0.05 "***" 0.01);
#delimit cr
estimates clear

* Panel D
use "$DAT/workingdata_age_work.dta", clear
gen g2smpl = mrg_mbdata2main == 3 ///
	& EDAD_MADRE >= 15 & EDAD_MADRE <= 45 & ANO_NAC >= 2007 & g2ok_blnctrls == 1

keep if g2smpl == 1

gen obs_dad = (EDAD_PADRE!=.)
keep if obs_dad==1

gen premature = SEMANAS<37
gen fgrate = PESO/SEMANAS

gen not_vlbw = 1 - vlbw
gen not_dead_at_a00 = 1 - dead_at_a00
gen not_gw_below_36 = 1 - gw_below_36
gen control_group = PESO_MADRE >= 1500 if PESO_MADRE != .
swindex PESO SEMANAS TALLA not_dead_at_a00 not_vlbw not_gw_below_36 fgrate, normby(control_group) gen(aindex)

foreach outcome in SEMANAS PESO TALLA dead_at_a00 premature vlbw fgrate aindex {

if "`outcome'"=="SEMANAS"     local lab "gwks"
if "`outcome'"=="PESO"        local lab "peso"
if "`outcome'"=="TALLA"       local lab "talla"
if "`outcome'"=="dead_at_a00" local lab "imr_ata00"
if "`outcome'"=="premature"   local lab "gw36"
if "`outcome'"=="vlbw"        local lab "bw1500"
if "`outcome'"=="fgrate"      local lab "fgrate"
if "`outcome'"=="aindex"      local lab "aindex"

rdbwselect `outcome' PESO_MADRE if sem32m == 1, c(1500) ///
            covs(${g2blnctrls} dbw1??0m) vce(cluster PESO_MADRE)
scalar RBW_o32 = e(h_mserd)
scalar LBW_o32 = e(h_mserd)
scalar maxwt_o32 = 1500 + RBW_o32
scalar minwt_o32 = 1500 - LBW_o32
gen g2_`outcome'_obw_o32 = PESO_MADRE >= minwt_o32 & PESO_MADRE <= maxwt_o32

rdrobust `outcome' PESO_MADRE if sem32m == 1, c(1500) scalepar(-1) all ///
          covs(${g2blnctrls} dbw1??0m) vce(cluster PESO_MADRE)
		  
eststo `outcome'o32
  
estadd scalar Hopt = e(h_l)
estadd scalar Nl   = e(N_h_l)
estadd scalar Nr   = e(N_h_r)
local ef = e(N_h_l) + e(N_h_r)
estadd scalar effopt = `ef'
sum `outcome' if sem32m == 1 & g2_`outcome'_obw_o32 == 1 & `outcome' !=. & PESO_MADRE !=.
estadd scalar dvmean = r(mean)

estimates restore `outcome'o32
estimates save "$OUT/g2_`lab'_blnctrls_o32_hps_relationship.ster", replace
}

clear 
local end _blnctrls
gen Robust =.
lab var Robust "Birth weight $<$ 1,500"
gen po32 = .
set obs 7
local i = 1
foreach outcome in gwks peso talla imr_ata00 gw36 bw1500 fgrate { 
    estimates use "$OUT/g2_`outcome'`end'_o32_hps_relationship.ster"
    eststo `outcome'o32
    local ++i
}

local ests 
#delimit ;
esttab gwkso32 pesoo32 tallao32 imr_ata00o32 
       gw36o32 bw1500o32 fgrateo32 aindexo32
using "$TAB/T2A_o32_relationship.tex", 
stats(dvmean N h_r, fmt(%8.3f %9.0gc %5.1f)
labels("\\ Mean of Dep.\ Var." "Observations" "Optimal Bandwidth")) keep("Robust") label
b(%-9.3f) se(%-9.3f) noobs nonotes nogaps mlabels(, none) nonumbers style(tex)
fragment replace noline starlevel("*" 0.10 "**" 0.05 "***" 0.01);
#delimit cr
estimates clear


*** Table E1: Selective survival and second generation outcomes - counterfactual analysis ***
do "$DO/SelectionReplication.do" // Create data

* Panel A-G:
use "$DAT/workingdata_age_work.dta", clear
gen g2smpl = mrg_mbdata2main == 3 ///
	& EDAD_MADRE >= 15 & EDAD_MADRE <= 45 & ANO_NAC >= 2007 & g2ok_blnctrls == 1

keep if g2smpl == 1

gen premature = SEMANAS<37
gen fgrate = PESO/SEMANAS

foreach var in SEMANAS PESO TALLA dead_at_a00 premature vlbw fgrate {

if "`var'"=="SEMANAS" local lab "gwks"
if "`var'"=="PESO" local lab "peso"
if "`var'"=="TALLA" local lab "talla"
if "`var'"=="dead_at_a00" local lab "deadata00"
if "`var'"=="premature" local lab "prem36"
if "`var'"=="vlbw" local lab "vlbw"
if "`var'"=="fgrate" local lab "fgrate"

rdrobust `var' PESO_MADRE if sem32m == 1, c(1500) scalepar(-1) all ///
          covs(${g2blnctrls} dbw1??0m) vce(cluster PESO_MADRE)
eststo `lab'o32

local ef = e(N_h_l) + e(N_h_r)
estadd scalar effopt = `ef'

estadd scalar Nr_i = 0

estimates restore `lab'o32
estimates save "$OUT/g2_ip05`lab'_o32_hps.ster", replace
}

use "$DAT/workingdata_age_work.dta", clear

gen g2smpl = mrg_mbdata2NAC == 3 ///
	& EDAD_MADRE >= 15 & EDAD_MADRE <= 45 & ANO_NAC >= 2007 & g2ok_blnctrls == 1
	
keep if g2smpl == 1

gen fgrate = PESO/SEMANAS

append using "$DAT\Selection\counterfactual_babies_o32_hps.dta", gen(app_cbs)

foreach var in SEMANAS PESO TALLA dead_at_a00 gw_below_36 vlbw fgrate {
 rename `var'_p5 `var'_p05
}

foreach var in SEMANAS PESO TALLA dead_at_a00 gw_below_36 vlbw fgrate {
 if "`var'"=="SEMANAS" local lab "gwks"
 if "`var'"=="PESO" local lab "peso"
 if "`var'"=="TALLA" local lab "talla"
 if "`var'"=="dead_at_a00" local lab "deadata00"
 if "`var'"=="gw_below_36" local lab "prem36"
 if "`var'"=="vlbw" local lab "vlbw"
 if "`var'"=="fgrate" local lab "fgrate"
 
 foreach pp of numlist 10 30 50 70 95 {
 gen `var'_i`pp' = `var' if app_cbs == 0
 replace `var'_i`pp' = `var'_p`pp' if app_cbs == 1 & SEMANAS_MADRE >= 32 & PESO_MADRE >= 1500
    
 rdrobust `var'_i`pp' PESO_MADRE if sem32m == 1, c(1500) scalepar(-1) all ///
           covs(${g2blnctrls} dbw1??0m) vce(cluster PESO_MADRE)
 estimates store select`pp'_`lab'o32
 
 local ef = e(N_h_l) + e(N_h_r)
 estadd scalar effopt = `ef'
 
 tempvar obw_all
 scalar maxwt_all = 1500 + e(h_r)
 scalar minwt_all = 1500 - e(h_l)
 qui gen `obw_all' = PESO_MADRE >= minwt_all & PESO_MADRE <= maxwt_all

 
 count if app_cbs == 1 & `obw_all' == 1 & `var'_i`pp' != . & PESO_MADRE != . & PESO_MADRE < 1500
 estadd scalar Nl_i = r(N)
 count if app_cbs == 1 & `obw_all' == 1 & `var'_i`pp' != . & PESO_MADRE != . & PESO_MADRE >= 1500
 estadd scalar Nr_i = r(N)
		   
 estimates restore select`pp'_`lab'o32
 estimates save "$OUT/g2_ip`pp'`lab'_o32_hps.ster", replace
 }
}

foreach var in gwks peso talla deadata00 prem36 vlbw fgrate {
 if ("`var'"=="peso") local clabel "fmt(%-9.2f)"
 else if ("`var'"=="deadata00") local clabel "fmt(%-9.3f)"
 else local clabel "fmt(%-9.4f)"
 
 clear
 set obs 19
 gen n = _n*5
 gen `var' = .

 cap estimates use "$OUT/g2_ip05`var'_o32_hps.ster"
 if _rc==0 estimates store `var'base 
 foreach num of numlist 10 30 50 70 95 {
  estimates use "$OUT/g2_ip`num'`var'_o32_hps.ster"
  replace `var' = _b[Robust] if n==`num'
  estimates store `var'`num'
 }
 
 #delimit ;
 esttab `var'base `var'10 `var'30 `var'50 `var'70 `var'95 using "$TAB/select_`var'.tex",
 replace keep(Robust) cells(b(`clabel' star) se(`clabel' par(( )) ))
 stats(effopt Nr_i, fmt(%9.0gc %9.0gc) label("\\ Effective Obs." "Imputed Obs."))
 starlevel("*" 0.1 "**" 0.05 "***" 0.01) tex nomtitles nonumbers fragment nodepvars
 coeflabels(Robust "Birth weight < 1,500") collabels(none) nolines;
 #delimit cr
}
estimates clear


*** Table E2: Impacts of Early Life Health Interventions on Number of Children ***
use "$DAT/workingdata_age_work.dta", clear

gen g1smpl = SEXO == 2 ///
  & EDAD_MADRE >= 15 & EDAD_MADRE <= 45 ///
  & mrg_mbdata2NAC == 1 ///
  & g1ok_blnctrls == 1
	
keep if g1smpl == 1

local list nbirths_at_15 nbirths_at_16 nbirths_at_17 nbirths_at_18 nbirths_at_19 /// 
           nbirths_at_20 nbirths_at_21 nbirths_at_22 nbirths_at_23 nbirths_at_24 ///
		   nbirths_at_25 nbirths_at_26 
								  
egen nbirths_all=rowtotal(`list'), missing					  
 
keep if sem32 == 1
 
rdbwselect nbirths_all PESO, c(1500) ///
           covs(${g1blnctrls} dbw1??0) vce(cluster PESO)
		   
scalar RBW_o32 = e(h_mserd)
scalar LBW_o32 = e(h_mserd)
scalar maxwt_o32 = 1500 + RBW_o32
scalar minwt_o32 = 1500 - LBW_o32
gen g1_nbirths_all = PESO >= minwt_o32 & PESO <= maxwt_o32
 
rdrobust nbirths_all PESO, c(1500) ///
scalepar(-1) all covs(${g1blnctrls} dbw1??0) vce(cluster PESO)
eststo all_all
     
estadd scalar Hopt = e(h_l)
estadd scalar Nl   = e(N_h_l)
estadd scalar Nr   = e(N_h_r)
local ef = e(N_h_l) + e(N_h_r)
estadd scalar effopt = `ef'
estadd scalar Ntot = e(N)

sum nbirths_all if sem32 == 1 & g1_nbirths_all == 1 & nbirths_all !=. & PESO !=.
estadd scalar dvmean = r(mean)   // All
sum nbirths_all if sem32 == 1 & PESO >= minwt_o32 & PESO < 1500 & nbirths_all !=. & PESO !=.
estadd scalar dvmean_l = r(mean) // Left
sum nbirths_all if sem32 == 1 & PESO >= 1500 & PESO <= maxwt_o32 & nbirths_all !=. & PESO !=.
estadd scalar dvmean_r = r(mean) // Right
 
gen Robust = .
lab var Robust "Birth weight $<$ 1,500"
 
#delimit ;
esttab all_all using "$TAB\TA18_o32.tex",
stats(Ntot dvmean h_r effopt dvmean_l Nl dvmean_r Nr,
fmt( %15.0gc %8.3f %5.1f %8.0gc %8.3f %8.0gc %8.3f %8.0gc)
labels("\midrule Observations" "Mean of Dep.\ Var." "Optimal Bandwidth"
       "Effective Observations" "\midrule Mean of Dep.\ Var. (left)" 
	   "Observations (left)" "\midrule Mean of Dep.\ Var. (right)"
       "Observations (right)")) keep(Robust) label
b(%-9.3f) se(%-9.3f) noobs nonotes nogaps mlabels(, none) nonumbers style(tex)
fragment replace noline starlevel("*" 0.10 "**" 0.05 "***" 0.01);
#delimit cr
estimates clear


*** Table E3: Is there Policy-Driven Selection intro Childbirth? RD Estimates ***
use "$DAT/workingdata_age_work.dta", clear

keep if SEMANAS_MADRE>=32& SEMANAS_MADRE!=.

gen marrieda = EST_CIV_MADRE==2 

gen ageMum  = EDAD_MADRE
gen educMum = CURSO_MADRE if NIVEL_MADRE==4 | NIVEL_MADRE==5 
replace educMum = CURSO_MADRE+8 if NIVEL_MADRE==2 | NIVEL_MADRE==3
replace educMum = CURSO_MADRE+12 if NIVEL_MADRE==1
gen employedMum = activ_m==1

gen ageDad  = EDAD_PADRE
gen educDad = CURSO_PADRE if NIVEL_PADRE==4 | NIVEL_PADRE==5 
replace educDad = CURSO_PADRE+8 if NIVEL_PADRE==2 | NIVEL_PADRE==3
replace educDad = CURSO_PADRE+12 if NIVEL_PADRE==1
gen employedDad = activ_p==1

gen urbano= URBANO_RURAL==1

gen difedad = abs(ageDad-ageMum)

gen obs_dad = (EDAD_PADRE!=.)

local i=1
foreach var of varlist ageMum educMum employedMum ageDad educDad employedDad marrieda difedad obs_dad urbano {
 keep if SEMANAS_MADRE>=32& SEMANAS_MADRE!=.
	
 eststo: rdrobust `var' PESO_MADRE, c(1500) scalepar(-1) covs(${g2blnctrls} dbw1??0m) vce(cluster PESO_MADRE) all 
 
 if ("`var'"=="ageMum") 	 local pv = normal(_b[Robust]/_se[Robust])
 if ("`var'"=="educMum")     local pv = 1-normal(_b[Robust]/_se[Robust])
 if ("`var'"=="employedMum") local pv = 1-normal(_b[Robust]/_se[Robust])
 if ("`var'"=="ageDad")      local pv = normal(_b[Robust]/_se[Robust])
 if ("`var'"=="educDad")     local pv = 1-normal(_b[Robust]/_se[Robust])
 if ("`var'"=="employedDad") local pv = 1-normal(_b[Robust]/_se[Robust])
 if ("`var'"=="marrieda")    local pv = 1-normal(_b[Robust]/_se[Robust])
 if ("`var'"=="difedad")     local pv = normal(_b[Robust]/_se[Robust])
 if ("`var'"=="obs_dad")     local pv = 1-normal(_b[Robust]/_se[Robust])
 if ("`var'"=="urbano")      local pv = 1-normal(_b[Robust]/_se[Robust])
 
 local p_value = string(`pv', "%05.3f")
 estadd local p_value "[`p_value']": est`i'
 
 estadd scalar Nl   = e(N_h_l)
 estadd scalar Nr   = e(N_h_r)
 local ef = e(N_h_l)+e(N_h_r)
 estadd scalar effopt = `ef'
 local BW = e(h_l)
 local maxwt = 1500+`BW'
 local minwt = 1500-`BW'
 sum `var' if PESO_MADRE >=`minwt' & PESO_MADRE<=`maxwt'
 estadd scalar dvmean = r(mean)
	
 estimates save "$OUT/g2_`var'_blnctrls_o32_hps.ster", replace
 local ++i
}

clear 
local end _blnctrls
gen Robust =.
lab var Robust "Birth weight $<$ 1,500"
gen po32 = .
set obs 10
local i = 1
foreach outcome in ageMum educMum employedMum ageDad educDad employedDad marrieda difedad obs_dad urbano { 
 estimates use "$OUT/g2_`outcome'`end'_o32_hps.ster"
 eststo est`i'
 replace po32=e(pv_rb) in `i'
 local ++i
}
rename po32 pval
qsharpenedp pval

local l = 1
foreach outcome in ageMum educMum employedMum ageDad educDad employedDad marrieda difedad obs_dad urbano {
    sum bky06 in `l'
    estadd scalar qpv = r(mean): est`l'
    local ++l
}

local ests est1 est2 est3 est4 est5 est6 est7 est8 est9 est10
#delimit ;
esttab `ests' using "$TAB/fig7_o32_hps.tex",
replace booktabs cells(b(fmt(%-9.3f) star) se(fmt(%-9.3f) par(( )))) label
stats(p_value dvmean N h_r effopt Nl Nr qpv, 
      fmt(%05.3f %05.3f %12.0gc %5.1f %9.0gc %9.0gc %9.0gc %05.3f) 
      label(" " "\\ Mean of Dep. Var." "Observations" "Optimal Bandwidth" 
            "Effective Observations" "Observations (left)" "Observations (right)" 
			"q-sharpened p-value"))
nonotes nogaps mlabels(, none) nonumbers style(tex) fragment noline keep(Robust) varlabels(Robust "Birth weight $<$ 1,500")
collabels(none) starlevel("*" 0.1 "**" 0.05 "***" 0.01);
#delimit cr
estimates clear


*** Table E4: Mother’s Age Baseline Test by Parents’ Characteristics – Generation 2 ***
use "$DAT/workingdata_age_work.dta", clear

* G2: Bandwidth 245.5 Anderson index
gen	g2smpl = (mrg_mbdata2main==3 & EDAD_MADRE>=15 & EDAD_MADRE<=45 & ANO_NAC>=2007)
keep if g2smpl==1

keep if sem32m==1

keep if PESO_MADRE<=1500+245.5
keep if PESO_MADRE>=1500-245.5

* Groups
* Mother's Education: Above vs. Below Median (12 years)
sum edmom, d
local median=r(p50)
gen aoa_mother=.
replace aoa_mother=0 if  edmom <  `median'
replace aoa_mother=1 if  edmom >= `median'
* Grandmother's Education: Above vs. Below Median (9 years)
sum edgmom, d
local median=r(p50)
gen aoa_grandmother=.
replace aoa_grandmother=0 if  edgmom <  `median'
replace aoa_grandmother=1 if  edgmom >= `median'
* Parent's Age Difference: 0-5 Years vs. +6 Years
gen ageMum  = EDAD_MADRE
gen ageDad  = EDAD_PADRE
gen difedad = abs(ageDad-ageMum) // absolute gap
gen difedad_plus6 = .
replace difedad_plus6 = 0 if difedad<=5 // 0-5 
replace difedad_plus6 = 1 if difedad>5  // 6+
* Observed Father
gen obs_dad = (EDAD_PADRE!=.)


keep EDAD_MADRE vlbwm aoa_mother aoa_grandmother difedad_plus6 obs_dad

foreach y of numlist 15/26 {

 gen age`y' = (EDAD_MADRE==`y')
 eststo: reg age`y' vlbwm, level(99)

 foreach var of varlist aoa_mother aoa_grandmother difedad_plus6 obs_dad {
 
  foreach n of numlist 0 1 {
  
   eststo: reg age`y' vlbwm if `var'==`n', level(99)
   
  }
 }
}

* Estimates
foreach j of numlist 1/12 {
 foreach num of numlist 1/9 {
  local n=`num'+`j'*9-9
  local ests`j' `ests`j'' est`n'
 }
 local ++j
}

foreach n of numlist 1/11 {
 local age=14+`n'
 #delimit ;
 esttab `ests`n''
 using "$TAB\TA20_y`age'_o32.tex", 
 keep(vlbwm) coeflabels(vlbwm "Mother's age: `age'") 
 b(%-9.3f) se(%-9.3f) noobs nonotes nogaps mlabels(, none) nonumbers style(tex)
 fragment replace noline starlevel("*" 0.10 "**" 0.05 "***" 0.01);
 #delimit cr
}

#delimit ;
esttab `ests12'
using "$TAB\TA20_y26_o32.tex", 
keep(vlbwm) coeflabels(vlbwm "Mother's age: `age'") 
stats(N,fmt(%9.0gc) labels("\\ Observations")) 
b(%-9.3f) se(%-9.3f) noobs nonotes nogaps mlabels(, none) nonumbers style(tex)
fragment replace noline starlevel("*" 0.10 "**" 0.05 "***" 0.01);
#delimit cr
estimates clear

*-------------------------------------------------------------------------------
*--- Figures
*-------------------------------------------------------------------------------

*** Figure A2: Age of Mothers at Chilbirth ***

* Panel A and B:
use "$DAT/workingdata_age_work.dta", clear

graph set window fontface "Times New Roman"
set scheme plotplainblind

foreach y of numlist 1992 2018 {

 preserve
 keep if ANO_NAC==`y'

 #delimit ;
 histogram EDAD_MADRE if EDAD_MADRE<=26, freq discrete color(red%100)
           ylabel(,format("%9.0f")) 
		   ytitle("Number of births") xtitle("Age of mother at birth") 
	 	   addplot(histogram EDAD_MADRE if EDAD_MADRE>26, freq discrete 
		   color(green%100)) legend(off) xlabel(10(10)60);
 #delimit cr
 graph export "$FIG/nacimientos_edadmadre`y'.eps", replace
 restore
 
}

* Panel C 
#delimit ;
histogram EDAD_MADRE if EDAD_MADRE<=26, freq discrete color(red%100)
          ylabel(,format("%9.0f"))  
		  ytitle("Number of births") 
		  xtitle("Age of mother at birth") 
		  addplot(histogram EDAD_MADRE if EDAD_MADRE>26, freq discrete 
		  color(green%100)) legend(off) xlabel(10(10)60);
#delimit cr
graph export "$FIG/nacimientos_edadmadre.eps", replace

* Panel D 
preserve
keep if PESO_MADRE!=.
#delimit ;
histogram EDAD_MADRE if EDAD_MADRE<=26, freq discrete color(red%100)
          ylabel(,format("%9.0f"))  
		  ytitle("Number of births") 
		  xtitle("Age of mother at birth") 
		  addplot(histogram EDAD_MADRE if EDAD_MADRE>26, freq discrete 
		  color(green%100)) legend(off) xlabel(10(10)60);
#delimit cr
graph export "$FIG/nacimientos_edadmadreG1.eps", replace
restore
estimates clear


*** Figure A3: Observed Birth Spacing Between Children and their Future Siblings ***
use "$DAT/MERGED_DEF2NAC.dta", clear

drop if mrg_DEF2NAC==2
drop if NAC_ID_MADRE=="NA"

#delimit ;
preserve;
gen nbirth = 1;
collapse NAC_FECHA_NACIMIENTO_SIF (sum) nbirth,
               by(NAC_ID_MADRE NAC_ANO_NAC NAC_MES_NAC);
bys NAC_ID_MADRE: gen N=_N;
bys NAC_ID_MADRE (NAC_ANO_NAC NAC_MES_NAC): gen n=_n;
bys NAC_ID_MADRE (NAC_ANO_NAC NAC_MES_NAC):
      gen birthSpacing = NAC_FECHA_NACIMIENTO_SIF[_n+1]-NAC_FECHA_NACIMIENTO_SIF;
#delimit cr
tempfile birthorders
save `birthorders'
restore

merge m:1 NAC_ID_MADRE NAC_ANO_NAC NAC_MES_NAC using `birthorders'
gen followingKids = N - n

* Panel A:
replace birthSpacing = birthSpacing/365
local hopts xlabel(1(1)20) scheme(lean2) xtitle("Spacing to Following Birth (Years)")
hist birthSpacing, `hopts' ylabel(, format(%04.2f))
graph export "$FIG/birthSpacingAll.eps", replace

*Panel B:
hist birthSpacing if NAC_PESO>=1300&NAC_PESO<=1700, `hopts' ylabel(, format(%04.2f))
graph export "$FIG/birthSpacingClose.eps", replace
estimates clear


*** Figure B1: Power Analysis - RD Models of Intergenerational Impacts ***
use "$DAT/workingdata_age_work.dta", clear
gen g2smpl = mrg_mbdata2main == 3 ///
	& EDAD_MADRE >= 15 & EDAD_MADRE <= 45 & ANO_NAC >= 2007 & g2ok_blnctrls == 1

keep if g2smpl == 1

*Graph range 
local range_PESO    "-300 300"
local range_TALLA   "-1.5 1.5"
local range_SEMANAS "-1.5 1.5"
*Graph step 
local step_PESO    100
local step_TALLA   0.5
local step_SEMANAS 0.5

foreach var in SEMANAS PESO TALLA {
 
 #delimit ;
 rdpow `var' PESO_MADRE if sem32m==1, c(1500) scalepar(-1) all plot
	   graph_range(`range_`var'') graph_step(`step_`var'') covs(${g2blnctrls})
	   graph_options(scheme(lean2) ytitle("Power") xtitle("Tau")
	   legend(rows(1) position(6)) xline(0, lcolor(red) lpattern(shortdash))); 
 graph export "$FIG/NAC_`var'_o32_Power.eps", replace;
		
 rdpow `var' PESO_MADRE if sem32m==1, c(1500) scalepar(-1) all plot
	   graph_range(`range_`var'') graph_step(`step_`var'') covs(${g2blnctrls} dbw1??0m)
	   graph_options(scheme(lean2) ytitle("Power") xtitle("Tau")
	   legend(rows(1) position(6)) xline(0, lcolor(red) lpattern(shortdash)));
 graph export "$FIG/NAC_`var'_o32_heaping_Power.eps", replace;
 #delimit cr
		
}
estimates clear


*** Figure B2: Birth weight Assignment Thresholds and Infant Mortality ***

* Panel A: 
use "$DAT/workingdata_age_work.dta", clear

keep if ANO_NAC >= 1992 & ANO_NAC<=2007 
       *& g2ok_blnctrls == 1 & g2smpl = mrg_mbdata2NAC == 3 ///
	   *& EDAD_MADRE >= 15 & EDAD_MADRE <= 45
keep if SEMANAS>=32& SEMANAS!=.
keep if PESO>=1385 & PESO<=1615

preserve
    drop if mrg_DEF2NAC==2
    replace dead_at_a00=. if PESO==1500
	replace dead_at_a00=. if PESO==1400
	replace dead_at_a00=. if PESO==1600
    keep if PESO>=1500 & PESO<=1600 & SEMANAS>=32
    gen wt=abs(1600-PESO)/100
    qui reg dead_at_a00 PESO [pw=wt]
    predict dead_at_a00_hat
    tempfile D_o32
    keep dead_at_a00_hat PESO
    drop if dead_at_a00_hat==. & PESO==.
    gen side="D"
    save "`D_o32'"
restore

preserve
    drop if mrg_DEF2NAC==2
    replace dead_at_a00=. if PESO==1500
	replace dead_at_a00=. if PESO==1400
	replace dead_at_a00=. if PESO==1600
    keep if PESO<1500 & PESO>=1400 & SEMANAS>=32
    gen wt=abs(1400-PESO)/100
    qui reg dead_at_a00 PESO [pw=wt]
    predict dead_at_a00_hat
    tempfile I_o32
    keep dead_at_a00_hat PESO
    drop if dead_at_a00_hat==. & PESO==.
    gen side="I"
    save "`I_o32'"
    append using "`D_o32'"
    tempfile B_o32
    save "`B_o32'"
restore

    preserve
    drop if mrg_DEF2NAC==2
    drop if PESO==1500
    drop if PESO==1400
    drop if PESO==1600
    keep if SEMANAS>=32
	foreach n of numlist 1400(10)1480 1520(10)1600{
        gen bin_`n'=.
        qui replace bin_`n'=1 if abs(PESO-`n')<=15 
    }
    gen bin_1490=.
    replace bin_1490=1 if abs(PESO-1490)<=15 & PESO<1500
    gen bin_1510=.
    replace bin_1510=1 if abs(PESO-1510)<=15 & PESO>=1500
    foreach n of numlist 1400(10)1490 1510(10)1600{
        qui sum dead_at_a00 if bin_`n'==1
        local mean_`n'=r(mean)
    }
    gen bin=.
    gen mean_imrt=.
    local i=1
    foreach n of numlist 1400(10)1490 1510(10)1600{
        qui replace bin=`n'              in `i'
        qui replace mean_imrt=`mean_`n'' in `i'
        local ++i
    }
    keep bin mean_imrt
    keep if bin!=. & mean_imrt!=.
    append using "`B_o32'"
	save "$DAT/blnfinal.dta", replace

local farea fcolor(gs12%50) clcolor(orange) clstyle(solid)
twoway line dead_at_a00_hat PESO if side=="I", `farea' clwidth(medthick) ///
|| line dead_at_a00_hat PESO if side=="D", `farea' clwidth(medthick) ///
|| scatter mean_imrt bin, ms(Oh) mcolor(blue) msize(medlarge) ///
ylabel(, format(%03.2f)) ///
xline(1500, lcolor(red)) scheme(plotplain) legend(off) ///
ytitle("Infant Mortality") xtitle("Mother's Birth Weight")
graph export "$FIG/imrt_o32_BLN_19922007_despues.eps", replace
restore

* Panel B and C:
use "$DAT/workingdata_age_work.dta", clear

keep if SEMANAS>=32& SEMANAS!=.

foreach year of numlist 2018 2001 {
 preserve
 keep if ANO_NAC >= 1992 & ANO_NAC<=`year' 
         *& g2ok_blnctrls == 1 & g2smpl = mrg_mbdata2NAC == 3 ///
         *& EDAD_MADRE >= 15 & EDAD_MADRE <= 45
	
 rdrobust dead_at_a00 PESO, c(1500) scalepar(-1) all covs(${g1blnctrls} dbw1??0) vce(cluster PESO)
 local opt = e(h_l)
    
 gen cond = PESO-1500
 keep if abs(cond)<=`opt'
 gen gram10 = floor(PESO/20)*20

 sum PESO
 local min = floor(r(min)/20)*20
 if mod(`opt', 20)<10 replace gram10 = gram10+20 if gram10==`min'
    
 bys gram10: egen meanWt = mean(dead_at_a00)
 bys gram10: gen N=_n
 bys gram10: gen wt=_N
    
 local farea fcolor(gs12%50) clcolor(orange) clstyle(solid)
 twoway qfitci dead_at_a00 PESO if PESO<1500, `farea' ///
     || qfitci dead_at_a00 PESO if PESO>1500, `farea' ///
     || scatter meanWt gram10 [aw=wt] if N==1, ms(Oh) mcolor(blue) ///
     xline(1500, lcolor(red)) scheme(plotplain) legend(off) ///
     ylabel(, format(%03.2f)) ///
     ytitle("Infant Mortality") xtitle("Weight in grams")
	
 graph export "$FIG/imrt_o32_optimal_1992`year'.eps", replace
 restore
}
estimates clear


*** Figure B3: Birth Weight Assignment Thresholds and Infant Mortality (Early Years Only) ***
use "$DAT/workingdata_age_work.dta", clear

keep if SEMANAS>=32& SEMANAS!=.

foreach year of numlist 2001 2004 2005 2006 {
preserve
 keep if ANO_NAC >= 1992 & ANO_NAC<=`year' 
		 *& g2ok_blnctrls == 1 & g2smpl = mrg_mbdata2NAC == 3 ///
		 *& EDAD_MADRE >= 15 & EDAD_MADRE <= 45
	
 rdrobust dead_at_a00 PESO, c(1500) scalepar(-1) all covs(${g1blnctrls} dbw1??0) vce(cluster PESO)
 local opt = e(h_l)
    
 gen cond = PESO-1500
 keep if abs(cond)<=`opt'
 gen gram10 = floor(PESO/20)*20
 ***AQUI ANTES TUVIMOS sum PESO_MADRE, y lo cambie por sum PESO. Ahora todo funciona
 sum PESO
 local min = floor(r(min)/20)*20
 if mod(`opt', 20)<10 replace gram10 = gram10+20 if gram10==`min'
    
 bys gram10: egen meanWt = mean(dead_at_a00)
 bys gram10: gen N=_n
 bys gram10: gen wt=_N
    
 local farea fcolor(gs12%50) clcolor(orange) clstyle(solid)
 twoway qfitci dead_at_a00 PESO if PESO<1500, `farea' ///
     || qfitci dead_at_a00 PESO if PESO>1500, `farea' ///
     || scatter meanWt gram10 [aw=wt] if N==1, ms(Oh) mcolor(blue) ///
     xline(1500, lcolor(red)) scheme(plotplain) legend(off) ///
     ylabel(, format(%03.2f)) ///
     ytitle("Infant Mortality") xtitle("Weight in grams")
	
 graph export "$FIG/imrt_o32_optimal_1992`year'.eps", replace
 restore
}
estimates clear


*** Figure B4: Discontinuities in Hospitalization in Babies Born at <32 Gestational Weeks ***
use "$DAT/workingdata_age_work.dta", clear

gen g1smpl = EDAD_MADRE >= 15 & EDAD_MADRE <= 45 ///
	& mrg_mbdata2NAC == 1 ///
	& g1ok_blnctrls == 1

keep if g1smpl == 1

sum ANO_NAC
scalar max_ano_nac = r(max)
scalar min_ano_nac = r(min)
scalar max_age = max_ano_nac - min_ano_nac

drop days_y25 days_y26 days_y27

foreach var of varlist days_y?? {

 dis "`var'"
 
 local aa = substr("`var'", -2, 2)
 local a = real("`aa'")
	
 gen days0_y`aa' = `var'
 replace days0_y`aa' = 0 if days_y`aa' == . & nadmssn_y`aa' != . & ndischrg_y`aa' != .
 label var days0_y`aa' "Number of days in hospital at age `a' (inputing 0 for no hospitalization)"
	
  if `a' >= 0 & `a' <= `=max_age'  {
   #delimit ;
   cap rdrobust days0_y`aa' PESO if sem32 == 0, c(1500) scalepar(-1) all
           covs(${g1blnctrls} dbw1??0) vce(cluster PESO);
   #delimit cr
		
	
   cap estimates save "$OUT/g1_hdays0at`aa'_u32_hps.ster", replace	
	}
}

clear 
set more off
set obs 25
gen age  = _n-1
gen EST  = .
gen LB   = .
gen UB   = .
gen LB90 = .
gen UB90 = .
estimates clear
foreach age of numlist 0(1)9 {
 estimates use "$OUT/g1_hdays0at0`age'_u32_hps.ster"
 eststo hosp
   qui: replace EST = _b[Robust] if age==`age'
   qui: replace LB  = _b[Robust] + invnormal(0.025)*_se[Robust] if  age==`age'
   qui: replace UB  = _b[Robust] + invnormal(0.975)*_se[Robust] if  age==`age'
   qui: replace LB90  = _b[Robust] + invnormal(0.05)*_se[Robust] if age==`age'
   qui: replace UB90  = _b[Robust] + invnormal(0.95)*_se[Robust] if age==`age'
}
foreach age of numlist 10(1)24 {
 estimates use "$OUT/g1_hdays0at`age'_u32_hps.ster"
 eststo hosp
   qui: replace EST = _b[Robust] if age==`age'
   qui: replace LB  = _b[Robust] + invnormal(0.025)*_se[Robust] if  age==`age'
   qui: replace UB  = _b[Robust] + invnormal(0.975)*_se[Robust] if  age==`age'
   qui: replace LB90  = _b[Robust] + invnormal(0.05)*_se[Robust] if age==`age'
   qui: replace UB90  = _b[Robust] + invnormal(0.95)*_se[Robust] if age==`age'
}
graph set window fontface "Times New Roman"
set scheme plotplainblind
#delimit ;
twoway rcap LB90 UB90 age, ylabel(,format("%3.0f")) lwidth(thick)
    || rcap LB UB age, ylabel(,format("%3.0f")) lcolor(black)
    || scatter EST age, mc(black) ms(S) ysize(6) xsize(6.5)
yline(0, lpattern(dash)) xtitle("Exposed Child's Age") legend(off)
xlabel(0(2)24) ytitle("Days of hospitalization");
#delimit cr
graph export "$FIG/hospitalization_u32.eps", replace
estimates clear


*** Figure B5: Balance Tests - Generation 1 family and Birth Characteristics ***
use "$DAT/workingdata_age_work.dta", clear

keep if EDAD_MADRE >= 15 & EDAD_MADRE <= 45 
drop if PESO_MADRE != . // Drop Gen2

gen marrieda = EST_CIV_MADRE==2 

gen ageMum  = EDAD_MADRE
gen educMum = CURSO_MADRE if NIVEL_MADRE==4 | NIVEL_MADRE==5 
replace educMum = CURSO_MADRE+8 if NIVEL_MADRE==2 | NIVEL_MADRE==3
replace educMum = CURSO_MADRE+12 if NIVEL_MADRE==1
gen employedMum = activ_m==1 

gen ageDad  = EDAD_PADRE
gen educDad = CURSO_PADRE if NIVEL_PADRE==4 | NIVEL_PADRE==5 
replace educDad = CURSO_PADRE+8 if NIVEL_PADRE==2 | NIVEL_PADRE==3
replace educDad = CURSO_PADRE+12 if NIVEL_PADRE==1
gen employedDad = activ_p==1

gen urbano= URBANO_RURAL==1

gen difedad = abs(ageDad-ageMum)

foreach var of varlist ageMum educMum employedMum ageDad educDad employedDad marrieda difedad urbano {
    if `"`var'"'=="marrieda" {
		local xx "Married"
		local ff  %04.2f
	}
	if `"`var'"'=="difedad" {
		local xx "Parent's Difference in Age"
		local ff  %04.2f
	}
	if `"`var'"'=="ageMum" {
		local xx "Mother´s Age"
		local ff  %03.1f
	}
	if `"`var'"'=="educMum" {
		local xx "Mother´s Education"
		local ff  %03.1f
	}
	if `"`var'"'=="employedMum" {
		local xx "Mother´s Employment"
		local ff  %04.2f
	}
	if `"`var'"'=="ageDad" {
		local xx "Father´s Age"
		local ff  %03.1f
	}
	if `"`var'"'=="educDad" {
		local xx "Father´s Education"
		local ff  %03.1f
	}
	if `"`var'"'=="employedDad" {
		local xx "Father´s Employment"
		local ff  %04.2f
	}
	if `"`var'"'=="urbano" {
		local xx "Urban Status"
		local ff  %04.2f
	}
	
	keep if SEMANAS>=32& SEMANAS!=.
	
    preserve

    #delimit ;
    rdrobust `var' PESO, c(1500) scalepar(-1) all
    covs(${g1blnctrls} dbw1??0) vce(cluster PESO);
    #delimit cr
    local opt = e(h_l)
    display `opt'
    gen cond = PESO-1500
    keep if abs(cond)<=`opt'
    gen gram10 = floor(PESO/20)*20
	sum PESO
	local min = floor(r(min)/20)*20
	if mod(`opt', 20)<10 replace gram10 = gram10+20 if gram10==`min'
    *drop if gram10==1500
    bys gram10: egen meanWt = mean(`var')
    bys gram10: gen N=_n
    bys gram10: gen wt=_N
    
    local farea fcolor(gs12%50) clcolor(orange) clstyle(solid)
    #delimit ;
    twoway qfitci `var' PESO if PESO<1500, `farea'
    || qfitci `var' PESO if PESO>1500, `farea'
    || scatter meanWt gram10 [aw=wt] if N==1, ms(Oh) mcolor(blue)
	ylabel(, format(`ff'))
    xline(1500, lcolor(red)) scheme(plotplain) legend(off)
    ytitle("`xx'") xtitle("Weight in grams");
    #delimit cr
    graph export "$FIG/figureA3_RDplot_`var'.eps", replace
    restore
}
estimates clear


*** Figure B6: Observable Birth Outcomes by Birth Weight
use "$DAT/workingdata_age_work.dta", clear

keep if PESO>=1300&PESO<=1700

gen A = 1
collapse TALLA SEMANAS (count) nt=TALLA ns=SEMANAS (sum) A, by(PESO)

* Panel A:
#delimit ;
twoway scatter TALLA PESO [aw=nt], xline(1500) scheme(plotplainblind) ms(Oh) mc(black)
 ||    qfit  TALLA PESO if PESO<1500, lpattern(dash) lcolor(orange) lwidth(thick)   
 ||    qfit  TALLA PESO if PESO>1500, lpattern(dash) lcolor(orange) lwidth(thick)   
ytitle("Size at birth (cms)") xtitle("Birth weight (grams)")
legend(order(2 "Split quadratic fit" 1 "Scatter") position(6) rows(1));
graph export "$FIG/sizeByWeight.eps", replace;
#delimit cr

* Panel B:
#delimit ;
twoway scatter SEMANAS PESO [aw=nt], xline(1500) scheme(plotplainblind) ms(Oh) mc(black)
 ||    qfit  SEMANAS PESO if PESO<1500, lpattern(dash) lcolor(orange) lwidth(thick)
 ||    qfit  SEMANAS PESO if PESO>1500, lpattern(dash) lcolor(orange) lwidth(thick)
ytitle("Gestational length (weeks)") xtitle("Birth weight (grams)")
legend(order(2 "Split quadratic fit" 1 "Scatter") position(6) rows(1));
graph export "$FIG/gestByWeight.eps", replace;
#delimit cr
estimates clear


*** Figure B7: Birth Weight Frequency in Adminstrative Records ***
use "$DAT/workingdata_age_work.dta", clear
keep if  ANO_NAC < 2007
graph set window fontface "Times New Roman"
set scheme plotplainblind

* Panel A:
histogram PESO, xtitle("Birth weights (grams)") fcolor(gray) ylabel(0(.0002).001)
graph export "$FIG/histFull2007.eps", replace

* Panel B:
preserve
drop if PESO<1300
drop if PESO>1700
histogram PESO, width(10) xtitle("Birth weights (grams)") fcolor(gray)
graph export "$FIG/histShort2007.eps", replace
restore
estimates clear


*** Figure D1: RD Estimates of Intergenerational Impacts Varying Estimation Bandwidths ***
use "$DAT/workingdata_age_work.dta", clear
gen g2smpl = mrg_mbdata2main == 3 ///
	& EDAD_MADRE >= 15 & EDAD_MADRE <= 45 & ANO_NAC >= 2007 & g2ok_blnctrls == 1

keep if g2smpl == 1

gen premature = SEMANAS<37
gen fgrate = PESO/SEMANAS

gen not_vlbw = 1 - vlbw
gen not_dead_at_a00 = 1 - dead_at_a00
gen not_gw_below_36 = 1 - gw_below_36
gen control_group = PESO_MADRE >= 1500 if PESO_MADRE != .
swindex PESO SEMANAS TALLA not_dead_at_a00 not_vlbw not_gw_below_36 fgrate, normby(control_group) gen(aindex)

foreach outcome in SEMANAS PESO TALLA dead_at_a00 premature vlbw fgrate aindex {
if "`outcome'"=="SEMANAS"     local lab "gweeks"
if "`outcome'"=="PESO"        local lab "peso"
if "`outcome'"=="TALLA"       local lab "talla"
if "`outcome'"=="dead_at_a00" local lab "deadata00"
if "`outcome'"=="premature"   local lab "prem36"
if "`outcome'"=="vlbw"        local lab "vlbw"
if "`outcome'"=="fgrate"      local lab "fgrate"
if "`outcome'"=="aindex"      local lab "aindex"

 foreach h of numlist 90(10)300 {
  rdrobust `outcome' PESO_MADRE if sem32m == 1, c(1500) scalepar(-1) all ///
           covs(${g2blnctrls} dbw1??0m) vce(cluster PESO_MADRE) h(`h')
  eststo `outcome'o32
 
  estimates save "$OUT/g2_`lab'h`h'_blnctrls_o32_hps.ster", replace
 }
}

clear
set obs 22
gen bw   = .
local i=1
foreach h of numlist 90(10)300 {
 replace bw = `h' in `i'
 local ++i
}
gen EST  = .
gen LB   = .
gen UB   = .
gen LB90 = .
gen UB90 = .
estimates clear

foreach outcome in SEMANAS PESO TALLA dead_at_a00 premature vlbw fgrate aindex {
 if "`outcome'"=="SEMANAS"     local lab "gweeks"
 if "`outcome'"=="PESO"        local lab "peso"
 if "`outcome'"=="TALLA"       local lab "talla"
 if "`outcome'"=="dead_at_a00" local lab "deadata00"
 if "`outcome'"=="premature"   local lab "prem36"
 if "`outcome'"=="vlbw"        local lab "vlbw"
 if "`outcome'"=="fgrate"      local lab "fgrate"
 if "`outcome'"=="aindex"      local lab "aindex"
 
 preserve
 foreach h of numlist 90(10)300 {

  estimates use "$OUT/g2_`lab'h`h'_blnctrls_o32_hps.ster"
  eststo bandwidth
  replace EST   = _b[Robust] if bw==`h'
  replace LB    = _b[Robust] + invnormal(0.025)*_se[Robust] if bw==`h'
  replace UB    = _b[Robust] + invnormal(0.975)*_se[Robust] if bw==`h'
  replace LB90  = _b[Robust] + invnormal(0.05)*_se[Robust]  if bw==`h'
  replace UB90  = _b[Robust] + invnormal(0.95)*_se[Robust]  if bw==`h'
 }
 graph set window fontface "Times New Roman"
 set scheme plotplainblind
 #delimit ;
 twoway connected EST bw, lcolor(black) lpattern(solid) msymbol(Oh) mcolor(color)
     || rcap LB UB bw, lwidth(medthick) lcolor(gray)
     || rcap LB90 UB90 bw, lwidth(thick) lcolor(midblue)
 yline(0, lpattern(dash)) xtitle("Bandwidth") legend(off)
 xlabel(90(20)300) ytitle("RDD Bias-Corrected Estimate (Robust CIs)");
 #delimit cr
 graph export "$FIG/g2_`lab'_bandwith_o32.eps", replace
 restore
}
estimates clear


*** Figure D2: RD Estimates of Intergenerational Impacts Varying Estimation Bandwidths with Constant Relative Bias***
use "$DAT/workingdata_age_work.dta", clear
gen g2smpl = mrg_mbdata2main == 3 ///
	& EDAD_MADRE >= 15 & EDAD_MADRE <= 45 & ANO_NAC >= 2007 & g2ok_blnctrls == 1

keep if g2smpl == 1

gen premature = SEMANAS<37
gen fgrate = PESO/SEMANAS

gen not_vlbw = 1 - vlbw
gen not_dead_at_a00 = 1 - dead_at_a00
gen not_gw_below_36 = 1 - gw_below_36
gen control_group = PESO_MADRE >= 1500 if PESO_MADRE != .
swindex PESO SEMANAS TALLA not_dead_at_a00 not_vlbw not_gw_below_36 fgrate, normby(control_group) gen(aindex)

foreach outcome in SEMANAS PESO TALLA dead_at_a00 premature vlbw fgrate aindex {
if "`outcome'"=="SEMANAS"     local lab "gweeks"
if "`outcome'"=="PESO"        local lab "peso"
if "`outcome'"=="TALLA"       local lab "talla"
if "`outcome'"=="dead_at_a00" local lab "deadata00"
if "`outcome'"=="premature"   local lab "prem36"
if "`outcome'"=="vlbw"        local lab "vlbw"
if "`outcome'"=="fgrate"      local lab "fgrate"
if "`outcome'"=="aindex"      local lab "aindex"

 rdbwselect `outcome' PESO_MADRE if sem32m == 1, c(1500) ///
            covs(${g2blnctrls} dbw1??0m) vce(cluster PESO_MADRE)
  local rho = e(h_mserd)/e(b_mserd)      

 foreach h of numlist 90(10)300 {
 			
  rdrobust `outcome' PESO_MADRE if sem32m == 1, c(1500) scalepar(-1) all ///
           covs(${g2blnctrls} dbw1??0m) vce(cluster PESO_MADRE) h(`h') rho(`rho')
  eststo `outcome'o32
 
  estimates save "$OUT/g2_`lab'h`h'b_blnctrls_o32_hps.ster", replace
 }
}

clear
set obs 22
gen bw   = .
local i=1
foreach h of numlist 90(10)300 {
 replace bw = `h' in `i'
 local ++i
}
gen EST  = .
gen LB   = .
gen UB   = .
gen LB90 = .
gen UB90 = .
estimates clear

foreach outcome in SEMANAS PESO TALLA dead_at_a00 premature vlbw fgrate aindex {
 if "`outcome'"=="SEMANAS"     local lab "gweeks"
 if "`outcome'"=="PESO"        local lab "peso"
 if "`outcome'"=="TALLA"       local lab "talla"
 if "`outcome'"=="dead_at_a00" local lab "deadata00"
 if "`outcome'"=="premature"   local lab "prem36"
 if "`outcome'"=="vlbw"        local lab "vlbw"
 if "`outcome'"=="fgrate"      local lab "fgrate"
 if "`outcome'"=="aindex"      local lab "aindex"
 
 preserve
 foreach h of numlist 90(10)300 {

  estimates use "$OUT/g2_`lab'h`h'b_blnctrls_o32_hps.ster"
  eststo bandwidth
  replace EST   = _b[Robust] if bw==`h'
  replace LB    = _b[Robust] + invnormal(0.025)*_se[Robust] if bw==`h'
  replace UB    = _b[Robust] + invnormal(0.975)*_se[Robust] if bw==`h'
  replace LB90  = _b[Robust] + invnormal(0.05)*_se[Robust]  if bw==`h'
  replace UB90  = _b[Robust] + invnormal(0.95)*_se[Robust]  if bw==`h'
 }
 graph set window fontface "Times New Roman"
 set scheme plotplainblind
 #delimit ;
 twoway connected EST bw, lcolor(black) lpattern(solid) msymbol(Oh) mcolor(color)
     || rcap LB UB bw, lwidth(medthick) lcolor(gray)
     || rcap LB90 UB90 bw, lwidth(thick) lcolor(midblue)
 yline(0, lpattern(dash)) xtitle("Bandwidth") legend(off)
 xlabel(90(20)300) ytitle("RDD Bias-Corrected Estimate (Robust CIs)");
 #delimit cr
 graph export "$FIG/g2_`lab'_biasbandwith_o32.eps", replace
 restore
}
estimates clear


*** Figure D3: Descriptive Plots of Parental Policy Receipt and Child Health Measures (All births) ***
use "$DAT/workingdata_age_work.dta", clear

gen g2smpl = mrg_mbdata2main == 3 ///
	& EDAD_MADRE >= 15 & EDAD_MADRE <= 45 & ANO_NAC >= 2007 & g2ok_blnctrls == 1

keep if g2smpl == 1

gen premature = SEMANAS<37
gen fgrate = PESO/SEMANAS

foreach var of varlist PESO SEMANAS TALLA vlbw premature fgrate  {
 if `"`var'"'=="PESO" local xx "Child's birth weight"
 if `"`var'"'=="TALLA" local xx "Child's gestational length"
 if `"`var'"'=="SEMANAS" local xx "Child's size at birth"
 if `"`var'"'=="premature" local xx "Prematurity (child)"
 if `"`var'"'=="vlbw" local xx "Very low birth weight (child)"
 if `"`var'"'=="fgrate" local xx "Fetal growth rate"
     
 preserve
 
 keep if SEMANAS_MADRE!=.
    
 #delimit ;
 rdrobust `var' PESO_MADRE, c(1500) scalepar(-1) all
           covs(${g2blnctrls} dbw1??0m) vce(cluster PESO_MADRE);
 #delimit cr
 local opt = e(h_l)
    
 gen cond = PESO_MADRE-1500
 keep if abs(cond)<=`opt'
 gen gram10 = floor(PESO_MADRE/20)*20
 gen kernel = (PESO_MADRE-1500)/`opt'
 replace kernel = (1-abs(kernel))
    
 bys gram10: egen meanWt = mean(`var')
 bys gram10: gen N=_n
 bys gram10: gen wt=_N
    
    
 local farea fcolor(gs12%50) clcolor(orange) clstyle(solid)
 #delimit ;
 twoway qfitci `var' PESO_MADRE if PESO_MADRE<1500, `farea'
 || qfitci `var' PESO_MADRE if PESO_MADRE>1500, `farea'
 || scatter meanWt gram10 [aw=wt] if N==1, ms(Oh) mcolor(blue)
 xline(1500, lcolor(red)) scheme(plotplain) legend(off)
 ytitle("`xx'") xtitle("Mother's Birth Weight");
 #delimit cr
 graph export "$FIG/Figure2ALL_RDplot_`var'.eps", replace
 restore
}
estimates clear


*** Figure D4: Distributional Impacts of Early Life Health Interventions on Second Generation Health Stocks (All births)-0.1 ***
use "$DAT/workingdata_age_work.dta", clear
gen g2smpl = mrg_mbdata2main == 3 ///
	& EDAD_MADRE >= 15 & EDAD_MADRE <= 45 & ANO_NAC >= 2007 ///
	& g2ok_blnctrls == 1

keep if g2smpl == 1

foreach outcome in PESO SEMANAS {

if "`outcome'"=="SEMANAS" {
 local lab "gw"
 local sample "27(1)39"
}

else if "`outcome'"=="PESO" {
 local lab "bw"
 local sample "1000(250)4000"
 }
 
 foreach wt of numlist `sample' { 
 
  rdrobust `lab'_below_`wt' PESO_MADRE, c(1500) scalepar(-1) all ///
            covs(${g2blnctrls} dbw1??0m) vce(cluster PESO_MADRE)
		  
  eststo `outcome'o32_`wt'
  estimates save "$OUT/g2_`lab'`wt'_all_hps.ster", replace
 }
}

* Panel A:
clear
set obs 13
gen J    =_n
gen wt   = 1000+(_n-1)*250
gen EST  = .
gen LB   = .
gen UB   = .
gen LB90 = .
gen UB90 = .
estimates clear
 foreach wt of numlist 1000(250)4000 {
  estimates use "$OUT/g2_bw`wt'_all_hps.ster"
  eststo bw1000
  replace EST = _b[Robust] if wt==`wt'
  replace LB  = _b[Robust] + invnormal(0.025)*_se[Robust] if wt==`wt'
  replace UB  = _b[Robust] + invnormal(0.975)*_se[Robust] if wt==`wt'
  replace LB90  = _b[Robust] + invnormal(0.05)*_se[Robust] if wt==`wt'
  replace UB90  = _b[Robust] + invnormal(0.95)*_se[Robust] if wt==`wt'
 }
graph set window fontface "Times New Roman"
set scheme plotplainblind
#delimit ;
twoway rarea LB UB wt , ylabel(,format("%03.1f")) color(gs13%20)
    || rarea LB90 UB90 wt , ylabel(,format("%03.1f")) color(gs10%20)
    || scatter EST wt, mc(black) ms(S) ysize(6) xsize(6.5)
    yline(0, lpattern(dash)) xtitle("Child's birth weight threshold") legend(off)
    xlabel(1000(500)4000) ytitle("Effect of Maternal Policy Receipt");
#delimit cr
graph export "$FIG/DistBW_all.pdf", replace

*Panel B:
clear
set obs 13
gen J    =_n
gen wt   = _n+26
gen EST  = .
gen LB   = .
gen UB   = .
gen LB90 = .
gen UB90 = .
estimates clear
 foreach wt of numlist 27/39 {
  estimates use "$OUT/g2_gw`wt'_all_hps.ster"
  eststo bw1000
  replace EST = _b[Robust] if wt==`wt'
  replace LB  = _b[Robust] + invnormal(0.025)*_se[Robust] if wt==`wt'
  replace UB  = _b[Robust] + invnormal(0.975)*_se[Robust] if wt==`wt'
  replace LB90  = _b[Robust] + invnormal(0.05)*_se[Robust] if wt==`wt'
  replace UB90  = _b[Robust] + invnormal(0.95)*_se[Robust] if wt==`wt'
 }
graph set window fontface "Times New Roman"
set scheme plotplainblind
#delimit ;
twoway rarea LB UB wt , ylabel(,format("%03.1f")) color(gs13%20)
    || rarea LB90 UB90 wt , ylabel(,format("%03.1f")) color(gs10%20)
    || scatter EST wt, mc(black) ms(S) ysize(6) xsize(6.5)
    yline(0, lpattern(dash)) xtitle("Child's gestational length threshold") legend(off)
    xlabel(27(2)39) ytitle("Effect of Maternal Policy Receipt");
#delimit cr
graph export "$FIG/DistWeeks_all.pdf", replace
estimates clear


*** Figure D5: Impacts of Early Life Health Interventions on Fertility and Spontaneous Abortions (All births) ***

* Panel A:
use "$DAT/workingdata_age_work.dta", clear

gen g1smpl = SEXO == 2 ///
	& EDAD_MADRE >= 15 & EDAD_MADRE <= 45 ///
	& mrg_mbdata2NAC == 1 ///
	& g1ok_blnctrls == 1
	
keep if g1smpl == 1
keep if SEMANAS!=.

foreach y of numlist 15/26 {
 
 #delimit ;
  rdrobust nbirths_by_`y' PESO, c(1500) scalepar(-1) all
  covs(${g1blnctrls} dbw1??0) vce(cluster PESO);
 #delimit cr

 estimates save "$OUT/g1f_nbirths_by`y'_blnctrls_all_hps.ster", replace
 }

clear
set obs 12
gen age  = _n+14
gen EST  = .
gen LB   = .
gen UB   = .
gen LB90 = .
gen UB90 = .
estimates clear
foreach age of numlist 15/26 {
 estimates use "$OUT/g1f_nbirths_by`age'_blnctrls_all_hps.ster"
 eststo fert
 replace EST = _b[Robust] if age==`age'
 replace LB  = _b[Robust] + invnormal(0.025)*_se[Robust] if  age==`age'
 replace UB  = _b[Robust] + invnormal(0.975)*_se[Robust] if  age==`age'
 replace LB90  = _b[Robust] + invnormal(0.05)*_se[Robust] if age==`age'
 replace UB90  = _b[Robust] + invnormal(0.95)*_se[Robust] if age==`age'
}
graph set window fontface "Times New Roman"
set scheme plotplainblind
#delimit ;
twoway rcap LB90 UB90 age, ylabel(,format("%03.2f")) lwidth(thick)
    || rcap LB UB age, ylabel(,format("%03.2f")) lcolor(black)
    || scatter EST age, mc(black) ms(S) ysize(6) xsize(6.5)
yline(0, lpattern(dash)) xtitle("Exposed Girl's Age") legend(off)
xlabel(15(1)26) ytitle("Number of Additional Births");
#delimit cr
graph export "$FIG/fertility_all.eps", replace

* Panel B:
use "$DAT/workingdata_age_work.dta", clear

gen g1smpl = SEXO == 2 ///
	& EDAD_MADRE >= 15 & EDAD_MADRE <= 45 ///
	& mrg_mbdata2NAC == 1 ///
	& g1ok_blnctrls == 1
	
keep if g1smpl == 1
keep if SEMANAS!=.

foreach y of numlist 15/26 {

 * Any birth by age:
 gen had_birth_by_`y' = (nbirths_by_`y' > 0) if nbirths_by_`y'!= .
 label var had_birth_by_`y' "Any birth by age `num'"
  
 #delimit ;
  rdrobust had_birth_by_`y' PESO, c(1500) scalepar(-1) all
  covs(${g1blnctrls} dbw1??0) vce(cluster PESO);
 #delimit cr

 estimates save "$OUT/g1f_nbirth_by`y'_blnctrls_o32_hps.ster", replace
 }

clear
set obs 12
gen age  = _n+14
gen EST  = .
gen LB   = .
gen UB   = .
gen LB90 = .
gen UB90 = .
estimates clear
foreach age of numlist 15/26 {
 estimates use "$OUT/g1f_nbirth_by`age'_blnctrls_o32_hps.ster"
 eststo fert
 replace EST = _b[Robust] if age==`age'
 replace LB  = _b[Robust] + invnormal(0.025)*_se[Robust] if  age==`age'
 replace UB  = _b[Robust] + invnormal(0.975)*_se[Robust] if  age==`age'
 replace LB90  = _b[Robust] + invnormal(0.05)*_se[Robust] if age==`age'
 replace UB90  = _b[Robust] + invnormal(0.95)*_se[Robust] if age==`age'
}
graph set window fontface "Times New Roman"
set scheme plotplainblind
#delimit ;
twoway rcap LB90 UB90 age, ylabel(,format("%03.2f")) lwidth(thick)
    || rcap LB UB age, ylabel(,format("%03.2f")) lcolor(black)
    || scatter EST age, mc(black) ms(S) ysize(6) xsize(6.5)
yline(0, lpattern(dash)) xtitle("Exposed Girl's Age") legend(off)
xlabel(15(1)26) ytitle("Any Births");
#delimit cr
graph export "$FIG/fertility2_all.eps", replace

* Panel C:
use "$DAT/workingdata_age_work.dta", clear

gen g1smpl = SEXO == 2 ///
	& EDAD_MADRE >= 15 & EDAD_MADRE <= 45 ///
	& mrg_mbdata2NAC == 1 ///
	& g1ok_blnctrls == 1

keep if g1smpl == 1

* Maximum horizon:
sum ANO_NAC
scalar max_ano_nac = r(max)
scalar min_ano_nac = r(min)
scalar max_age = max_ano_nac - min_ano_nac

* EEHH data coverage:
scalar min_ano_eehh = 2001
scalar max_ano_eehh = 2019

keep if SEMANAS!=.

foreach y of numlist 15/26 {

 * Number of abortions by age:
	egen numabrtns_by_y`y' = rowtotal(nadmssn_y00_abrtn-nadmssn_y`y'_abrtn), m
	replace numabrtns_by_y`y' = . if ANO_NAC + `y' < min_ano_eehh
	replace numabrtns_by_y`y' = . if ANO_NAC + `y' > max_ano_eehh
	replace numabrtns_by_y`y' = . if ANO_NAC + `y' > ANO_DEF & mrg_DEF2NAC == 3
	label var numabrtns_by_y`y' "Number of abortions by age `y'"

  * Any abortion by age:
	gen anyabrtn_by_y`y' = numabrtns_by_y`y' > 0 if numabrtns_by_y`y'!= .
	label var anyabrtn_by_y`y' "Any abortion by age `y'"
 
 #delimit ;
  rdrobust anyabrtn_by_y`y' PESO, c(1500) scalepar(-1) all
  covs(${g1blnctrls} dbw1??0) vce(cluster PESO);
 #delimit cr

 estimates save "$OUT/g1f_anyabrtn_by_y`y'_blnctrls_all_hps.ster", replace
 }

clear
set obs 12
gen age  = _n+14
gen EST  = .
gen LB   = .
gen UB   = .
gen LB90 = .
gen UB90 = .
estimates clear
foreach age of numlist 15/26 {
 estimates use "$OUT/g1f_anyabrtn_by_y`age'_blnctrls_all_hps.ster"
 eststo fert
 replace EST = _b[Robust] if age==`age'
 replace LB  = _b[Robust] + invnormal(0.025)*_se[Robust] if  age==`age'
 replace UB  = _b[Robust] + invnormal(0.975)*_se[Robust] if  age==`age'
 replace LB90  = _b[Robust] + invnormal(0.05)*_se[Robust] if age==`age'
 replace UB90  = _b[Robust] + invnormal(0.95)*_se[Robust] if age==`age'
}
graph set window fontface "Times New Roman"
set scheme plotplainblind
#delimit ;
twoway rcap LB90 UB90 age, ylabel(,format("%03.2f")) lwidth(thick)
    || rcap LB UB age, ylabel(,format("%03.2f")) lcolor(black)
    || scatter EST age, mc(black) ms(S) ysize(6) xsize(6.5)
yline(0, lpattern(dash)) xtitle("Exposed Girl's Age") legend(off)
xlabel(15(1)26) ytitle("Number of Abortions");
#delimit cr
graph export "$FIG/abortions_all.eps", replace
estimates clear

* Panel D:
use "$DAT/workingdata_age_work.dta", clear

keep if SEMANAS!=.
gen sex_male = SEXO == 1 if SEXO == 1 | SEXO == 2

gen age  = _n+16 in 1/10
gen EST  = .
gen LB   = .
gen UB   = .
gen LB90 = .
gen UB90 = .

foreach age of numlist 18(1)26 {
    #delimit ;
    rdrobust sex_male PESO_MADRE if EDAD_MADRE<=`age', c(1500) scalepar(-1) all
    covs(${g2blnctrls} dbw1??0m) vce(cluster PESO_MADRE);
    #delimit cr
    
    replace EST = _b[Robust] if age==`age'
    replace LB  = _b[Robust] + invnormal(0.025)*_se[Robust] if  age==`age'
    replace UB  = _b[Robust] + invnormal(0.975)*_se[Robust] if  age==`age'
    replace LB90  = _b[Robust] + invnormal(0.05)*_se[Robust] if age==`age'
    replace UB90  = _b[Robust] + invnormal(0.95)*_se[Robust] if age==`age'
}    
graph set window fontface "Times New Roman"
set scheme plotplainblind
#delimit ;
twoway rcap LB90 UB90 age, ylabel(,format("%03.2f")) lwidth(thick)
||     rcap LB UB age, ylabel(,format("%03.2f")) lcolor(black)
||    scatter EST age, mc(black) ms(S) ysize(6) xsize(6.5)
yline(0, lpattern(dash)) xtitle("Exposed Girl's Age") legend(off)
xlabel(18(1)26) ytitle("Child is Male");
#delimit cr
graph export "$FIG/sexratio_all.eps", replace
estimates clear


*** Figure D6: Estimated RDD Effects by Second Generation Mother’s Age ***
use "$DAT/workingdata_age_work.dta", clear

gen g2smpl = mrg_mbdata2NAC == 3 ///
			& EDAD_MADRE >= 15 & EDAD_MADRE <= 45 ///
			& ANO_NAC >= 2007 ///
			& g2ok_blnctrls == 1

keep if g2smpl == 1
keep if sem32m == 1

* Panel A and B:
foreach yvar of varlist PESO SEMANAS {
 if `"`yvar'"'=="PESO"    local yl "Child's birth weight"
 if `"`yvar'"'=="SEMANAS" local yl "Child's gestational period"

 gen group=_n in 1/6
 gen EST= .
 gen LB = .
 gen UB = .
 
 local opts c(1500) scalepar(-1) all covs(${g2blnctrls} dbw1??0m) vce(cluster PESO_MADRE)
    
 rdrobust `yvar' PESO_MADRE if EDAD_MADRE<=16, `opts'
 replace EST = _b[Robust] if group==1
 replace LB  = _b[Robust] + invnormal(0.025)*_se[Robust] if group==1
 replace UB  = _b[Robust] + invnormal(0.975)*_se[Robust] if group==1
    
 rdrobust `yvar' PESO_MADRE if EDAD_MADRE>=16&EDAD_MADRE<=18, `opts'
 replace EST = _b[Robust] if group==2
 replace LB  = _b[Robust] + invnormal(0.025)*_se[Robust] if group==2
 replace UB  = _b[Robust] + invnormal(0.975)*_se[Robust] if group==2
    
 rdrobust `yvar' PESO_MADRE if EDAD_MADRE>=18&EDAD_MADRE<=20, `opts'
 replace EST = _b[Robust] if group==3
 replace LB  = _b[Robust] + invnormal(0.025)*_se[Robust] if group==3
 replace UB  = _b[Robust] + invnormal(0.975)*_se[Robust] if group==3
    
 rdrobust `yvar' PESO_MADRE if EDAD_MADRE>=20&EDAD_MADRE<=22, `opts'
 replace EST = _b[Robust] if group==4
 replace LB  = _b[Robust] + invnormal(0.025)*_se[Robust] if group==4
 replace UB  = _b[Robust] + invnormal(0.975)*_se[Robust] if group==4
    
 rdrobust `yvar' PESO_MADRE if EDAD_MADRE>=22&EDAD_MADRE<=24, `opts'
 replace EST = _b[Robust] if group==5
 replace LB  = _b[Robust] + invnormal(0.025)*_se[Robust] if group==5
 replace UB  = _b[Robust] + invnormal(0.975)*_se[Robust] if group==5
    
 rdrobust `yvar' PESO_MADRE if EDAD_MADRE>=24&EDAD_MADRE<=26, `opts'
 replace EST = _b[Robust] if group==6
 replace LB  = _b[Robust] + invnormal(0.025)*_se[Robust] if group==6
 replace UB  = _b[Robust] + invnormal(0.975)*_se[Robust] if group==6
    
 graph set window fontface "Times New Roman"
 set scheme plotplainblind
 #delimit ;
 twoway rcap LB UB group, ylabel(,format("%8.0gc")) lcolor(black)
     || scatter EST group, mc(black) ms(S) ysize(6) xsize(6.5)
 yline(0, lpattern(dash)) xtitle("Age at Birth") 
 xlabel(1 "{&le} 16" 2 "16-18" 3 "18-20" 4 "20-22" 5 "22-24" 6 "24-26")
 ytitle("`yl'") legend(order(2 "Effect Size" 1 "95% CI") pos(6) rows(1));
 #delimit cr
 graph export "$FIG/effectsByAge_`yvar'.eps", replace
 drop group EST LB UB
}
estimates clear


*** Figure E1: Intergenerational Transmission Under Alternative Health and Fertility Counterfactuals ***
do "$DO/SelectionReplication.do" // Create data

* Panel A-D:
use "$DAT/workingdata_age_work.dta", clear

gen g2smpl = mrg_mbdata2NAC == 3 ///
	& EDAD_MADRE >= 15 & EDAD_MADRE <= 45 & ANO_NAC >= 2007 & g2ok_blnctrls == 1
	
keep if g2smpl == 1

gen fgrate = PESO/SEMANAS

append using "$DAT\Selection\counterfactual_babies_o32_hps.dta", gen(app_cbs)

cap mkdir "$OUT\selection"

foreach var in PESO SEMANAS TALLA fgrate {
 rename `var'_p5 `var'_p05
}

foreach var in PESO SEMANAS TALLA fgrate {
 if "`var'"=="PESO" local lab "peso"
 if "`var'"=="SEMANAS" local lab "gwks"
 if "`var'"=="TALLA" local lab "talla"
 if "`var'"=="fgrate" local lab "fgrate"
 
 foreach pp of numlist 5(5)95  {
  if `pp'==5 local ft 05
  else       local ft=`pp'
  
 gen `var'_i`ft' = `var' if app_cbs == 0
 replace `var'_i`ft' = `var'_p`ft' if app_cbs == 1 & SEMANAS_MADRE >= 32 & PESO_MADRE >= 1500
    
 rdrobust `var'_i`ft' PESO_MADRE if sem32m == 1, c(1500) scalepar(-1) all ///
           covs(${g2blnctrls} dbw1??0m) vce(cluster PESO_MADRE)
 estimates store select`ft'_`lab'o32

 estimates save "$OUT/Selection/g2_ip`ft'`lab'_o32_hps.ster", replace
 }
 
}

use "$DAT/workingdata_age_work.dta", clear

gen g2smpl = mrg_mbdata2NAC == 3 ///
	& EDAD_MADRE >= 15 & EDAD_MADRE <= 45 ///
	& ANO_NAC >= 2007 ///
	& g2ok_blnctrls == 1

keep if g2smpl == 1

gen fgrate = PESO/SEMANAS

foreach s of numlist 1 2 3 5 12 4 6 7 8 9 10 {

 cap mkdir "$OUT\selection`s'"

 preserve
 cap append using "$DAT\Selection\counterfactual_babies_o32_fixed`s'.dta", gen(app_cbs)	
 
 foreach var in PESO SEMANAS TALLA fgrate {
 rename `var'_p5 `var'_p05
 }
 
 foreach var in  PESO SEMANAS TALLA fgrate {
 if "`var'"=="PESO" local lab "peso"
 if "`var'"=="SEMANAS" local lab "gwks"
 if "`var'"=="TALLA" local lab "talla"
 if "`var'"=="fgrate" local lab "fgrate"
 
 foreach pp of numlist 5(5)95  {
  if `pp'==5 local ft 05
  else       local ft=`pp'
  
  gen `var'_i`ft' = `var' if app_cbs == 0
  replace `var'_i`ft' = `var'_p`ft' if app_cbs == 1 & SEMANAS_MADRE >= 32 & PESO_MADRE >= 1500
    
  rdrobust `var'_i`ft' PESO_MADRE if sem32m == 1, c(1500) scalepar(-1) all ///
            covs(${g2blnctrls} dbw1??0m) vce(cluster PESO_MADRE)
  estimates store select`ft'_`lab'o32

  estimates save "$OUT/selection`s'/g2_ip`ft'`lab'_o32_hps.ster", replace
 }
}
 restore
}

foreach var in peso gwks talla fgrate {

 local lring 0
 local ys
 
 local clabel fmt(%-9.4f)
 if `"`var'"'=="peso" {
  local top = 155
  local textloc 77.5
  local lpos 6
  local clabel fmt(%-9.2f)
 }
 if `"`var'"'=="gwks" {
  local top = 0.63
  local textloc 0.315
  local lpos 7
  local ys ylabel(-1(0.5)0.5, format(%04.2f))
 }
 if `"`var'"'=="talla" {
  local top = 0.74
  local textloc 0.37
  local lpos 7
 }
 if `"`var'"'=="fgrate" {
  local top = 4
  local textloc 2
  local lpos 7
 }
 
 clear   
 set obs 19
 gen n = _n*5
 gen `var' = .
 
 foreach num of numlist 5(5)95 {
   if `num'==5 local ft 05
   else        local ft=`num'
   cap estimates use "$OUT/Selection/g2_ip`ft'`var'_o32_hps.ster"
   if _rc==0 replace `var' = _b[Robust] if n==`num'
 }
 
 gen high = `top'
 gen low  = 0
 gen n2 = (_n-1)*105.5/19

 foreach s of numlist 1 2 3 5 12 4 6 7 8 9 10 {
  gen `var'`s' = .
  foreach num of numlist 5(5)95 {
   if `num'==5 local ft 05
   else        local ft=`num'
   cap estimates use "$OUT/selection`s'/g2_ip`ft'`var'_o32_hps.ster"
   if _rc==0 replace `var'`s' = _b[Robust] if n==`num'
  }
 }
 set scheme white_jet
 #delimit ;
 twoway rarea high low n2, color(gs14)
 ||     connected `var'   n, ms(sh)
 ||     connected `var'12 n, ms(oh)
 ||     connected `var'4  n, ms(dh)
 ||     connected `var'6  n, ms(th)
 ||     connected `var'8  n, ms(x)
 ||     connected `var'10 n, ms(d)
 ||     connected `var'1  n, ms(s)
 ||     connected `var'2  n, ms(o)
 ytitle("Estimated Intergenerational Returns")
 xtitle("Imputed Centile: Health")
 text(`textloc' 50 "Positive Intergenerational Transmission") yline(0)
 legend(order(2 "E[fertility | BW]" 3 "0" 4 "0.2"
              5 "0.4" 6 "0.6" 7 "0.8" 8 "1" 9 "2")
 ring(`lring') pos(`lpos') rows(2)) `ys';
 graph export "$FIG/`var'100.eps", replace;
 #delimit cr
}


*** Figure E2: Weight at Birth by Mothers's Age ***
use "$DAT/MERGED_DEF2NAC.dta", clear

keep if NAC_EDAD_MADRE >= 15 & NAC_EDAD_MADRE <= 45

gen vlbw = NAC_PESO<1500
gen premature = NAC_SEMANAS<37

set scheme plotplainblind

* Panel A: 
preserve
tempfile FA20_a
#delimit ;
binsreg NAC_PESO NAC_EDAD_MADRE, nbins(31) dots(0,0) line(2,2) ci(2,2) cb(2,2)
ytitle("Child's Weight at Birth") xtitle("Mother's Age") savedata(`FA20_a')
dotsplotopt(ms(Oh) mcolor(black)) ciplotopt(lcolor(gs3)) lineplotopt(lcolor(lavender))
cbplotopt(color(blue%40)) level(99) xlabel(15(5)45) ylabel(3100(50)3350);
graph export "$FIG/birthWeightAge.pdf", replace;
#delimit cr
use `FA20_a', clear
keep dots_x dots_fit CI_l CI_r
save "$EXP/FA20_a.dta", replace
restore

* Panel B: Manual
preserve
drop if NAC_ID_MADRE=="NA"

reghdfe NAC_PESO i.NAC_EDAD_MADRE, absorb(NAC_ID_MADRE) level(99) 
scalar df = e(df_r)

margins i.NAC_EDAD_MADRE, post level(99)
mat v = e(V)

gen age=_n+14 in 1/31 
gen b=.
gen se=.
gen lb=.
gen ub=.

local i=1
foreach n of numlist 15/45 {
 replace b = _b[`n'.NAC_EDAD_MADRE] in `i'
 replace se = sqrt(v[`i',`i'])      in `i'
 replace lb = _b[`n'.NAC_EDAD_MADRE] - invttail(df,0.005)*sqrt(v[`i',`i']) in `i'
 replace ub = _b[`n'.NAC_EDAD_MADRE] + invttail(df,0.005)*sqrt(v[`i',`i']) in `i'
 local ++i
}

keep age b lb ub
keep if b!=.
save "$EXP/FA20_b.dta", replace

#delimit ;
twoway scatter b age, ms(Oh) mcolor(black) || 
       rcap ub lb age, lcolor(black)
	   xlabel(15(5)45) ytitle("Child's Weight at Birth") xtitle("Mother's Age")
	   legend(off);
graph export "$FIG/birthWeightAge_FE.pdf", replace;
#delimit cr 
restore

* Panel C: 
preserve
tempfile FA20_c
#delimit ;
binsreg vlbw NAC_EDAD_MADRE, nbins(31) dots(0,0) line(2,2) ci(2,2) cb(2,2)
ytitle("Child Very Low Birth Weight") xtitle("Mother's Age") savedata(`FA20_c')
dotsplotopt(ms(Oh) mcolor(black)) ciplotopt(lcolor(gs3)) lineplotopt(lcolor(lavender))
cbplotopt(color(blue%40)) ylabel(.005(.005).03, format(%04.3f)) level(99) xlabel(15(5)45);
graph export "$FIG/vlbwAge.pdf", replace;
#delimit cr
use `FA20_c', clear
keep dots_x dots_fit CI_l CI_r
save "$EXP/FA20_c.dta", replace
restore

* Panel D: 
preserve
drop if NAC_ID_MADRE=="NA"

reghdfe vlbw i.NAC_EDAD_MADRE, absorb(NAC_ID_MADRE) level(99) 
scalar df = e(df_r)

margins i.NAC_EDAD_MADRE, post level(99)
mat v = e(V)

gen age=_n+14 in 1/31 
gen b=.
gen se=.
gen lb=.
gen ub=.

local i=1
foreach n of numlist 15/45 {
 replace b = _b[`n'.NAC_EDAD_MADRE] in `i'
 replace se = sqrt(v[`i',`i'])      in `i'
 replace lb = _b[`n'.NAC_EDAD_MADRE] - invttail(df,0.005)*sqrt(v[`i',`i']) in `i'
 replace ub = _b[`n'.NAC_EDAD_MADRE] + invttail(df,0.005)*sqrt(v[`i',`i']) in `i'
 local ++i
}

keep age b lb ub
keep if b!=.
save "$EXP/FA20_d.dta", replace

#delimit ;
twoway scatter b age, ms(Oh) mcolor(black) || 
       rcap ub lb age, lcolor(black)
	   xlabel(15(5)45) ytitle("Child Very Low Birth Weight") xtitle("Mother's Age")
	   legend(off) ylabel(, format(%04.3f));
graph export "$FIG/vlbwAge_FE.pdf", replace;
#delimit cr 
restore

* Panel E: 
preserve
tempfile FA20_e
#delimit ;
binsreg premature NAC_EDAD_MADRE, nbins(31) dots(0,0) line(2,2) ci(2,2) cb(2,2)
ytitle("Child Prematurity") xtitle("Mother's Age") savedata(`FA20_e')
dotsplotopt(ms(Oh) mcolor(black)) ciplotopt(lcolor(gs3)) lineplotopt(lcolor(lavender))
cbplotopt(color(blue%40)) ylabel(.05(.05).2, format(%04.2f)) level(99) xlabel(15(5)45);
graph export "$FIG/prematureAge.pdf", replace;
#delimit cr
use `FA20_e', clear
keep dots_x dots_fit CI_l CI_r
save "$EXP/FA20_e.dta", replace
restore

* Panel F: 
preserve
drop if NAC_ID_MADRE=="NA"

reghdfe premature i.NAC_EDAD_MADRE, absorb(NAC_ID_MADRE) level(99) 
scalar df = e(df_r)

margins i.NAC_EDAD_MADRE, post level(99)
mat v = e(V)

gen age=_n+14 in 1/31 
gen b=.
gen se=.
gen lb=.
gen ub=.

local i=1
foreach n of numlist 15/45 {
 replace b = _b[`n'.NAC_EDAD_MADRE] in `i'
 replace se = sqrt(v[`i',`i'])      in `i'
 replace lb = _b[`n'.NAC_EDAD_MADRE] - invttail(df,0.005)*sqrt(v[`i',`i']) in `i'
 replace ub = _b[`n'.NAC_EDAD_MADRE] + invttail(df,0.005)*sqrt(v[`i',`i']) in `i'
 local ++i
}

keep age b lb ub
keep if b!=.
save "$EXP/FA20_f.dta", replace

#delimit ;
twoway scatter b age, ms(Oh) mcolor(black) || 
       rcap ub lb age, lcolor(black)
	   xlabel(15(5)45) ytitle("Child Prematurity") xtitle("Mother's Age")
	   legend(off) ylabel(, format(%04.2f));
graph export "$FIG/prematureAge_FE.pdf", replace;
#delimit cr 
restore
estimates clear


*** Figure E3: Childbirth Around the 1500 Grams Threshold ***

* Panel A: Average childbirths by age
use "$DAT/workingdata_age_work.dta", clear

keep if PESO<=1500+134.4 // same as in Table A3
keep if PESO>=1500-134.4

gen g1smpl = SEXO == 2 ///
	& EDAD_MADRE >= 15 & EDAD_MADRE <= 45 ///
	& mrg_mbdata2NAC == 1
	*& g1ok_blnctrls == 1
	
keep if g1smpl == 1

collapse (mean)nbirths_by_15 (mean)nbirths_by_16 (mean)nbirths_by_17 ///
         (mean)nbirths_by_18 (mean)nbirths_by_19 (mean)nbirths_by_20 ///
		 (mean)nbirths_by_21 (mean)nbirths_by_22 (mean)nbirths_by_23 ///
		 (mean)nbirths_by_24 (mean)nbirths_by_25 (mean)nbirths_by_26, by(vlbw)

reshape long nbirths_by_, i(vlbw) j(age)	
	 
graph set window fontface "Times New Roman"
set scheme plotplainblind

#delimit ;
twoway connected nbirths_by_ age if vlbw==1 ||
       connected nbirths_by_ age if vlbw==0, 
       ylabel(,format("%9.1f")) ytitle("Average Births") 
       xtitle("Age of mother") xlabel(15(1)26)
	   legend(order(1 "Below Cutoff" 2 "Above Cutoff")
	   position(6) rows(1));
#delimit cr
graph export "$FIG/avgbirthsby_vlbw.eps", replace
estimates clear

* Panel B: Probability of childbirth at age
use "$DAT/workingdata_age_work.dta", clear

keep if PESO<=1500+134.4 // same as in Table A3
keep if PESO>=1500-134.4

gen g1smpl = SEXO == 2 ///
	& EDAD_MADRE >= 15 & EDAD_MADRE <= 45 ///
	& mrg_mbdata2NAC == 1
	*& g1ok_blnctrls == 1
	
keep if g1smpl == 1

foreach y of numlist 15/26 {
 gen had_birth_at_`y' = (nbirths_at_`y' > 0) if nbirths_at_`y'!= .
 label var had_birth_at_`y' "Any birth at age `num'"
}

collapse (mean)had_birth_at_15 (mean)had_birth_at_16 (mean)had_birth_at_17 ///
         (mean)had_birth_at_18 (mean)had_birth_at_19 (mean)had_birth_at_20 ///
		 (mean)had_birth_at_21 (mean)had_birth_at_22 (mean)had_birth_at_23 ///
		 (mean)had_birth_at_24 (mean)had_birth_at_25 (mean)had_birth_at_26, by(vlbw)

reshape long had_birth_at_, i(vlbw) j(age)	
	 
graph set window fontface "Times New Roman"
set scheme plotplainblind

drop if age==26
#delimit ;
twoway connected had_birth_at_ age if vlbw==1 ||
       connected had_birth_at_ age if vlbw==0, 
       ylabel(,format("%9.2f")) ytitle("Birth probability") 
       xtitle("Age of mother") xlabel(15(1)25)
	   legend(order(1 "Below Cutoff" 2 "Above Cutoff")
	   position(6) rows(1));
#delimit cr
graph export "$FIG/anybirthsat_vlbw.eps", replace
estimates clear

* Panel C: Probability of childbirth by age
use "$DAT/workingdata_age_work.dta", clear

keep if PESO<=1500+134.4 // same as in Table A3
keep if PESO>=1500-134.4

gen g1smpl = SEXO == 2 ///
	& EDAD_MADRE >= 15 & EDAD_MADRE <= 45 ///
	& mrg_mbdata2NAC == 1
	*& g1ok_blnctrls == 1
	
keep if g1smpl == 1

foreach y of numlist 15/26 {
 gen had_birth_by_`y' = (nbirths_by_`y' > 0) if nbirths_by_`y'!= .
 label var had_birth_by_`y' "Any birth by age `num'"
}

collapse (mean)had_birth_by_15 (mean)had_birth_by_16 (mean)had_birth_by_17 ///
         (mean)had_birth_by_18 (mean)had_birth_by_19 (mean)had_birth_by_20 ///
		 (mean)had_birth_by_21 (mean)had_birth_by_22 (mean)had_birth_by_23 ///
		 (mean)had_birth_by_24 (mean)had_birth_by_25 (mean)had_birth_by_26, by(vlbw)

reshape long had_birth_by_, i(vlbw) j(age)	
	 
graph set window fontface "Times New Roman"
set scheme plotplainblind

#delimit ;
twoway connected had_birth_by_ age if vlbw==1 ||
       connected had_birth_by_ age if vlbw==0, 
       ylabel(,format("%9.1f")) ytitle("Birth probability") 
       xtitle("Age of mother") xlabel(15(1)26)
	   legend(order(1 "Below Cutoff" 2 "Above Cutoff")
	   position(6) rows(1));
#delimit cr
graph export "$FIG/anybirthsby_vlbw.eps", replace
estimates clear

* Panel D: Probability of first childbirth at age
use "$DAT/workingdata_age_work.dta", clear

keep if PESO<=1500+134.4 // same as in Table A3
keep if PESO>=1500-134.4

gen g1smpl = SEXO == 2 ///
	& EDAD_MADRE >= 15 & EDAD_MADRE <= 45 ///
	& mrg_mbdata2NAC == 1
	*& g1ok_blnctrls == 1
	
keep if g1smpl == 1

foreach y of numlist 15/26 {
 gen first_birth_at_`y' = (nbirths_at_`y' == 1 & nbirths_by_`y' == 0) if nbirths_at_`y'!=. & nbirths_by_`y'!= .
 label var first_birth_at_`y' "First birth at age `num'"
}


collapse (mean)first_birth_at_15 (mean)first_birth_at_16 (mean)first_birth_at_17 ///
         (mean)first_birth_at_18 (mean)first_birth_at_19 (mean)first_birth_at_20 ///
		 (mean)first_birth_at_21 (mean)first_birth_at_22 (mean)first_birth_at_23 ///
		 (mean)first_birth_at_24 (mean)first_birth_at_25 (mean)first_birth_at_26, by(vlbw)

reshape long first_birth_at_, i(vlbw) j(age)	
	 
graph set window fontface "Times New Roman"
set scheme plotplainblind

drop if age==26
#delimit ;
twoway connected first_birth_at_ age if vlbw==1 ||
       connected first_birth_at_ age if vlbw==0, 
       ylabel(,format("%9.2f")) ytitle("First birth probability") 
       xtitle("Age of mother") xlabel(15(1)25)
	   legend(order(1 "Below Cutoff" 2 "Above Cutoff")
	   position(6) rows(1));
#delimit cr
graph export "$FIG/firstbirthat_vlbw.eps", replace
estimates clear


*** Figure E4: Probability of Childbirth Relative to Median ***

use "$DAT/workingdata_age_work.dta", clear

gen g1smpl = SEXO == 2 ///
	& EDAD_MADRE >= 15 & EDAD_MADRE <= 45 ///
	& mrg_mbdata2NAC == 1
	*& g1ok_blnctrls == 1
	
keep if g1smpl == 1

sum PESO, d
gen median=(PESO<r(p50)) // 3,300 grams


* Panel A: Average childbirths by age
preserve
collapse (mean)nbirths_by_15 (mean)nbirths_by_16 (mean)nbirths_by_17 ///
         (mean)nbirths_by_18 (mean)nbirths_by_19 (mean)nbirths_by_20 ///
		 (mean)nbirths_by_21 (mean)nbirths_by_22 (mean)nbirths_by_23 ///
		 (mean)nbirths_by_24 (mean)nbirths_by_25 (mean)nbirths_by_26, by(median)

reshape long nbirths_by_, i(median) j(age)	
	 
graph set window fontface "Times New Roman"
set scheme plotplainblind

#delimit ;
twoway connected nbirths_by_ age if median==1 ||
       connected nbirths_by_ age if median==0, 
       ylabel(,format("%9.1f")) ytitle("Average Births") 
       xtitle("Age of mother") xlabel(15(1)26)
	   legend(order(1 "Below Median" 2 "Above Median")
	   position(6) rows(1));
#delimit cr
graph export "$FIG/avgbirthsby_median.eps", replace
restore

* Panel B: Probability of childbirth at age
preserve
foreach y of numlist 15/26 {
 gen had_birth_at_`y' = (nbirths_at_`y' > 0) if nbirths_at_`y'!= .
 label var had_birth_at_`y' "Any birth at age `num'"
}

collapse (mean)had_birth_at_15 (mean)had_birth_at_16 (mean)had_birth_at_17 ///
         (mean)had_birth_at_18 (mean)had_birth_at_19 (mean)had_birth_at_20 ///
		 (mean)had_birth_at_21 (mean)had_birth_at_22 (mean)had_birth_at_23 ///
		 (mean)had_birth_at_24 (mean)had_birth_at_25 (mean)had_birth_at_26, by(median)

reshape long had_birth_at_, i(median) j(age)	
	 
graph set window fontface "Times New Roman"
set scheme plotplainblind

drop if age==26
#delimit ;
twoway connected had_birth_at_ age if median==1 ||
       connected had_birth_at_ age if median==0, 
       ylabel(,format("%9.2f")) ytitle("Birth probability") 
       xtitle("Age of mother") xlabel(15(1)25)
	   legend(order(1 "Below Median" 2 "Above Median")
	   position(6) rows(1));
#delimit cr
graph export "$FIG/anybirthsat_median.eps", replace
restore

* Panel C: Probability of childbirth by age
preserve
foreach y of numlist 15/26 {
 gen had_birth_by_`y' = (nbirths_by_`y' > 0) if nbirths_by_`y'!= .
 label var had_birth_by_`y' "Any birth by age `num'"
}

collapse (mean)had_birth_by_15 (mean)had_birth_by_16 (mean)had_birth_by_17 ///
         (mean)had_birth_by_18 (mean)had_birth_by_19 (mean)had_birth_by_20 ///
		 (mean)had_birth_by_21 (mean)had_birth_by_22 (mean)had_birth_by_23 ///
		 (mean)had_birth_by_24 (mean)had_birth_by_25 (mean)had_birth_by_26, by(median)

reshape long had_birth_by_, i(median) j(age)	
	 
graph set window fontface "Times New Roman"
set scheme plotplainblind

#delimit ;
twoway connected had_birth_by_ age if median==1 ||
       connected had_birth_by_ age if median==0, 
       ylabel(,format("%9.1f")) ytitle("Birth probability") 
       xtitle("Age of mother at birth") xlabel(15(1)26)
	   legend(order(1 "Below Median" 2 "Above Median")
	   position(6) rows(1));
#delimit cr
graph export "$FIG/anybirthsby_median.eps", replace
restore


* Panel D: Probability of first childbirth at age
preserve
foreach y of numlist 15/26 {
 gen first_birth_at_`y' = (nbirths_at_`y' == 1 & nbirths_by_`y' == 0) if nbirths_at_`y'!=. & nbirths_by_`y'!= .
 label var first_birth_at_`y' "First birth at age `num'"
}


collapse (mean)first_birth_at_15 (mean)first_birth_at_16 (mean)first_birth_at_17 ///
         (mean)first_birth_at_18 (mean)first_birth_at_19 (mean)first_birth_at_20 ///
		 (mean)first_birth_at_21 (mean)first_birth_at_22 (mean)first_birth_at_23 ///
		 (mean)first_birth_at_24 (mean)first_birth_at_25 (mean)first_birth_at_26, by(median)

reshape long first_birth_at_, i(median) j(age)	
	 
graph set window fontface "Times New Roman"
set scheme plotplainblind

drop if age==26
#delimit ;
twoway connected first_birth_at_ age if median==1 ||
       connected first_birth_at_ age if median==0, 
       ylabel(,format("%9.2f")) ytitle("First birth probability") 
       xtitle("Age of mother") xlabel(15(1)25)
	   legend(order(1 "Below Median" 2 "Above Median")
	   position(6) rows(1));
#delimit cr
graph export "$FIG/firstbirthat_median.eps", replace
restore
estimates clear


*** Figure E5: Weight at Birth by Mothers's Age and Parent Characteristics ***

* Panel A:
use "$DAT/workingdata_age_work.dta", clear

keep if EDAD_MADRE >= 15 & EDAD_MADRE <= 45

gen heduc_parents=.
replace heduc_parents=2 if  highest_ed3lvls_mp <= 2
replace heduc_parents=3 if  highest_ed3lvls_mp == 3

set scheme plotplainblind

foreach n of numlist 2 3 {
 preserve
 tempfile heduc`n'
 drop if ID_MADRE=="NA"

 reghdfe PESO i.EDAD_MADRE if heduc_parents==`n', absorb(ID_MADRE) level(99) 
 scalar df = e(df_r)

 margins i.EDAD_MADRE, post level(99)
 mat v = e(V)

 gen age`n'=_n+14 in 1/31 
 gen b`n'=.
 gen se`n'=.
 gen lb`n'=.
 gen ub`n'=.

 local i=1
 foreach j of numlist 15/45 {
  replace b`n'  = _b[`j'.EDAD_MADRE] in `i'
  replace se`n' = sqrt(v[`i',`i'])   in `i'
  replace lb`n' = _b[`j'.EDAD_MADRE] - invttail(df,0.005)*sqrt(v[`i',`i']) in `i'
  replace ub`n' = _b[`j'.EDAD_MADRE] + invttail(df,0.005)*sqrt(v[`i',`i']) in `i'
  local ++i
 }

 keep age`n' b`n' lb`n' ub`n'
 keep if b`n'!=.
 if `n'==3 replace age`n'=age`n'+.5
 save `heduc`n'', replace
 restore
}

use `heduc2', clear
merge 1:1 _n using `heduc3'
save "$EXP/FA30_b.dta", replace

#delimit ;
twoway scatter b2    age2, ms(Oh) mcolor(black) ||
       rcap lb2 ub2  age2, lcolor(black)        ||
	   scatter b3    age3, ms(Oh) mcolor(gs10)  ||
       rcap lb3 ub3  age3, lcolor(gs10)        
legend(order(2 "Secondary or Below" 4 "University") position(6) rows(1))
ytitle("Child's Weight at Birth") xtitle("Mother's Age")
xlabel(15(5)45);
graph export "$FIG/birthWeightAge_educ_FE.pdf", replace;
#delimit cr

* Panel B: 
use "$DAT/workingdata_age_work.dta", clear

keep if EDAD_MADRE >= 15 & EDAD_MADRE <= 45

gen heduc_parents=.
replace heduc_parents=2 if  highest_ed3lvls_mp <= 2
replace heduc_parents=3 if  highest_ed3lvls_mp == 3

set scheme plotplainblind

foreach n of numlist 2 3 {
 preserve
 tempfile heduc`n'
 drop if ID_MADRE=="NA"

 reghdfe vlbw i.EDAD_MADRE if heduc_parents==`n', absorb(ID_MADRE) level(99) 
 scalar df = e(df_r)

 margins i.EDAD_MADRE, post level(99)
 mat v = e(V)

 gen age`n'=_n+14 in 1/31 
 gen b`n'=.
 gen se`n'=.
 gen lb`n'=.
 gen ub`n'=.

 local i=1
 foreach j of numlist 15/45 {
  replace b`n'  = _b[`j'.EDAD_MADRE] in `i'
  replace se`n' = sqrt(v[`i',`i'])   in `i'
  replace lb`n' = _b[`j'.EDAD_MADRE] - invttail(df,0.005)*sqrt(v[`i',`i']) in `i'
  replace ub`n' = _b[`j'.EDAD_MADRE] + invttail(df,0.005)*sqrt(v[`i',`i']) in `i'
  local ++i
 }

 keep age`n' b`n' lb`n' ub`n'
 keep if b`n'!=.
 if `n'==3 replace age`n'=age`n'+.5
 save `heduc`n'', replace
 restore
}

use `heduc2', clear
merge 1:1 _n using `heduc3'
save "$EXP/FA30_d.dta", replace

#delimit ;
twoway scatter b2    age2, ms(Oh) mcolor(black) ||
       rcap lb2 ub2  age2, lcolor(black)        ||
	   scatter b3    age3, ms(Oh) mcolor(gs10)  ||
       rcap lb3 ub3  age3, lcolor(gs10)        
legend(order(2 "Secondary or Below" 4 "University") position(6) rows(1))
ytitle("Child Very Low Birth Weight") xtitle("Mother's Age")
xlabel(15(5)45) ylabel(, format(%04.2f));
graph export "$FIG/vlbwAge_educ_FE.pdf", replace;
#delimit cr

* Panel C: 
use "$DAT/workingdata_age_work.dta", clear

keep if EDAD_MADRE >= 15 & EDAD_MADRE <= 45

gen premature = SEMANAS<37

gen heduc_parents=.
replace heduc_parents=2 if  highest_ed3lvls_mp <= 2
replace heduc_parents=3 if  highest_ed3lvls_mp == 3

set scheme plotplainblind

foreach n of numlist 2 3 {
 preserve
 tempfile heduc`n'
 drop if ID_MADRE=="NA"

 reghdfe premature i.EDAD_MADRE if heduc_parents==`n', absorb(ID_MADRE) level(99) 
 scalar df = e(df_r)

 margins i.EDAD_MADRE, post level(99)
 mat v = e(V)

 gen age`n'=_n+14 in 1/31 
 gen b`n'=.
 gen se`n'=.
 gen lb`n'=.
 gen ub`n'=.

 local i=1
 foreach j of numlist 15/45 {
  replace b`n'  = _b[`j'.EDAD_MADRE] in `i'
  replace se`n' = sqrt(v[`i',`i'])   in `i'
  replace lb`n' = _b[`j'.EDAD_MADRE] - invttail(df,0.005)*sqrt(v[`i',`i']) in `i'
  replace ub`n' = _b[`j'.EDAD_MADRE] + invttail(df,0.005)*sqrt(v[`i',`i']) in `i'
  local ++i
 }

 keep age`n' b`n' lb`n' ub`n'
 keep if b`n'!=.
 if `n'==3 replace age`n'=age`n'+.5
 save `heduc`n'', replace
 restore
}

use `heduc2', clear
merge 1:1 _n using `heduc3'
save "$EXP/FA30_f.dta", replace

#delimit ;
twoway scatter b2    age2, ms(Oh) mcolor(black) ||
       rcap lb2 ub2  age2, lcolor(black)        ||
	   scatter b3    age3, ms(Oh) mcolor(gs10)  ||
       rcap lb3 ub3  age3, lcolor(gs10)        
legend(order(2 "Secondary or Below" 4 "University") position(6) rows(1))
ytitle("Child Premature") xtitle("Mother's Age")
xlabel(15(5)45) ylabel(, format(%04.2f));
graph export "$FIG/prematureAge_educ_FE.pdf", replace;
#delimit cr
estimates clear

* Panel D: 
use "$DAT/workingdata_age_work.dta", clear

keep if EDAD_MADRE >= 15 & EDAD_MADRE <= 45

gen ageMum  = EDAD_MADRE
gen ageDad  = EDAD_PADRE
gen difedad = abs(ageDad-ageMum) // absolute gap

gen lev_difedad = .
replace lev_difedad = 1 if difedad<=5 // 0-5 
replace lev_difedad = 2 if difedad>5  // 6+

set scheme plotplainblind

foreach n of numlist 1 2 {
 preserve
 tempfile difedad`n'
 drop if ID_MADRE=="NA"

 reghdfe PESO i.EDAD_MADRE if lev_difedad==`n', absorb(ID_MADRE) level(99) 
 scalar df = e(df_r)

 margins i.EDAD_MADRE, post level(99)
 mat v = e(V)

 gen age`n'=_n+14 in 1/31 
 gen b`n'=.
 gen se`n'=.
 gen lb`n'=.
 gen ub`n'=.

 local i=1
 foreach j of numlist 15/45 {
  replace b`n'  = _b[`j'.EDAD_MADRE] in `i'
  replace se`n' = sqrt(v[`i',`i'])   in `i'
  replace lb`n' = _b[`j'.EDAD_MADRE] - invttail(df,0.005)*sqrt(v[`i',`i']) in `i'
  replace ub`n' = _b[`j'.EDAD_MADRE] + invttail(df,0.005)*sqrt(v[`i',`i']) in `i'
  local ++i
 }

 keep age`n' b`n' lb`n' ub`n'
 keep if b`n'!=.
 if `n'==2 replace age`n'=age`n'+.5
 save `difedad`n'', replace
 restore
}

use `difedad1', clear
merge 1:1 _n using `difedad2'
save "$EXP/FA31_b.dta", replace

#delimit ;
twoway scatter b1    age1, ms(Oh) mcolor(black) ||
       rcap lb1 ub1  age1, lcolor(black)        ||
	   scatter b2    age2, ms(Oh) mcolor(gs10)  ||
       rcap lb2 ub2  age2, lcolor(gs10)        
legend(order(2 "0-5  Years" 4 ">5 Years") position(6) rows(1))
ytitle("Child's Weight at Birth") xtitle("Mother's Age")
xlabel(15(5)45);
graph export "$FIG/birthWeightAge_agediff_FE.pdf", replace;
#delimit cr

* Panel E: 
use "$DAT/workingdata_age_work.dta", clear

keep if EDAD_MADRE >= 15 & EDAD_MADRE <= 45

gen ageMum  = EDAD_MADRE
gen ageDad  = EDAD_PADRE
gen difedad = abs(ageDad-ageMum) // absolute gap

gen lev_difedad = .
replace lev_difedad = 1 if difedad<=5 // 0-5 
replace lev_difedad = 2 if difedad>5  // 6+

set scheme plotplainblind

foreach n of numlist 1 2 {
 preserve
 tempfile difedad`n'
 drop if ID_MADRE=="NA"

 reghdfe vlbw i.EDAD_MADRE if lev_difedad==`n', absorb(ID_MADRE) level(99) 
 scalar df = e(df_r)

 margins i.EDAD_MADRE, post level(99)
 mat v = e(V)

 gen age`n'=_n+14 in 1/31 
 gen b`n'=.
 gen se`n'=.
 gen lb`n'=.
 gen ub`n'=.

 local i=1
 foreach j of numlist 15/45 {
  replace b`n'  = _b[`j'.EDAD_MADRE] in `i'
  replace se`n' = sqrt(v[`i',`i'])   in `i'
  replace lb`n' = _b[`j'.EDAD_MADRE] - invttail(df,0.005)*sqrt(v[`i',`i']) in `i'
  replace ub`n' = _b[`j'.EDAD_MADRE] + invttail(df,0.005)*sqrt(v[`i',`i']) in `i'
  local ++i
 }

 keep age`n' b`n' lb`n' ub`n'
 keep if b`n'!=.
 if `n'==2 replace age`n'=age`n'+.5
 save `difedad`n'', replace
 restore
}

use `difedad1', clear
merge 1:1 _n using `difedad2'
save "$EXP/FA31_d.dta", replace

#delimit ;
twoway scatter b1    age1, ms(Oh) mcolor(black) ||
       rcap lb1 ub1  age1, lcolor(black)        ||
	   scatter b2    age2, ms(Oh) mcolor(gs10)  ||
       rcap lb2 ub2  age2, lcolor(gs10)        
legend(order(2 "0-5  Years" 4 ">5 Years") position(6) rows(1))
ytitle("Child Very Low Birth Weight") xtitle("Mother's Age")
xlabel(15(5)45) ylabel(, format(%04.2f));
graph export "$FIG/vlbwAge_agediff_FE.pdf", replace;
#delimit cr

* Panel F: 
use "$DAT/workingdata_age_work.dta", clear

keep if EDAD_MADRE >= 15 & EDAD_MADRE <= 45

gen premature = SEMANAS<37

gen ageMum  = EDAD_MADRE
gen ageDad  = EDAD_PADRE
gen difedad = abs(ageDad-ageMum) // absolute gap

gen lev_difedad = .
replace lev_difedad = 1 if difedad<=5 // 0-5 
replace lev_difedad = 2 if difedad>5  // 6+

set scheme plotplainblind

foreach n of numlist 1 2 {
 preserve
 tempfile difedad`n'
 drop if ID_MADRE=="NA"

 reghdfe premature i.EDAD_MADRE if lev_difedad==`n', absorb(ID_MADRE) level(99) 
 scalar df = e(df_r)

 margins i.EDAD_MADRE, post level(99)
 mat v = e(V)

 gen age`n'=_n+14 in 1/31 
 gen b`n'=.
 gen se`n'=.
 gen lb`n'=.
 gen ub`n'=.

 local i=1
 foreach j of numlist 15/45 {
  replace b`n'  = _b[`j'.EDAD_MADRE] in `i'
  replace se`n' = sqrt(v[`i',`i'])   in `i'
  replace lb`n' = _b[`j'.EDAD_MADRE] - invttail(df,0.005)*sqrt(v[`i',`i']) in `i'
  replace ub`n' = _b[`j'.EDAD_MADRE] + invttail(df,0.005)*sqrt(v[`i',`i']) in `i'
  local ++i
 }

 keep age`n' b`n' lb`n' ub`n'
 keep if b`n'!=.
 if `n'==2 replace age`n'=age`n'+.5
 save `difedad`n'', replace
 restore
}

use `difedad1', clear
merge 1:1 _n using `difedad2'
save "$EXP/FA31_f.dta", replace

#delimit ;
twoway scatter b1    age1, ms(Oh) mcolor(black) ||
       rcap lb1 ub1  age1, lcolor(black)        ||
	   scatter b2    age2, ms(Oh) mcolor(gs10)  ||
       rcap lb2 ub2  age2, lcolor(gs10)        
legend(order(2 "0-5  Years" 4 ">5 Years") position(6) rows(1))
ytitle("Child Premature") xtitle("Mother's Age")
xlabel(15(5)45) ylabel(, format(%04.2f));
graph export "$FIG/prematureAge_agediff_FE.pdf", replace;
#delimit cr
estimates clear

*-------------------------------------------------------------------------------
*--- Log Close
*-------------------------------------------------------------------------------

log close
