/* IntergenReplication.do        KTS/DCC/NLB               yyyy-mm-dd:2025-01-18
----|----1----|----2----|----3----|----4----|----5----|----6----|----7----|----8

This do file replicates the results provided in the paper "Estimating 
Inter­generational Returns to Medical Care: New Evidence from At­Risk Newborns" 
written by Damian Clarke, Nicolas Lillo Bustos and Kathya Tapia-Schythe.  
In certain cases these results will require the user-written, estout, 
reghdfe, rdbwselect, rdplot, rdrobust, binsreg, blindschemes, swindex 
personage, autofmt commands.  
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
global LOG "$DIR\log"
cap mkdir $LOG
global OUT "$DIR\results"
cap mkdir $OUT
global FIG "$DIR\figures"
cap mkdir $FIG
global TAB "$DIR\tables"
cap mkdir $TAB

* Log file:
log using "$LOG\IntergenReplication.log", replace

*-------------------------------------------------------------------------------
*--- Variables
*-------------------------------------------------------------------------------

global g1blnctrls edmom EDAD_MADRE married doc_aten bregion_?? byear_????
global g2blnctrls edgmom EDAD_ABUELA marriedm doc_atenm bregion_??m byear_????m

*-------------------------------------------------------------------------------
*--- Data Corrections 
*-------------------------------------------------------------------------------

use "$DAT/workingdata.dta", clear

*a.- Age
gen	g2smpl = (mrg_mbdata2main==3 & EDAD_MADRE>=15 & EDAD_MADRE<=45 & ANO_NAC>=2007)

personage FECHA_NACIMIENTO_SIF_MADRE FECHA_NACIMIENTO_SIF if g2smpl==1, gen(mother_age) 
replace EDAD_MADRE=mother_age if EDAD_MADRE>26 & g2smpl==1

drop g2smpl
drop if EDAD_MADRE<10 // 1 corrected observation

*save "$DAT/workingdata_age.dta", replace

*b.- Empleo
tab ACTIV_PADRE activ_p, m
replace activ_p=0 if ACTIV_PADRE=="0"
replace activ_p=1 if ACTIV_PADRE=="1"
replace activ_p=1 if ACTIV_PADRE=="2"
replace activ_p=. if ACTIV_PADRE=="3"
replace activ_p=. if ACTIV_PADRE==""
tab ACTIV_PADRE activ_p, m

tab ACTIV_MADRE activ_m, m
replace activ_m=0 if ACTIV_MADRE=="0"
replace activ_m=1 if ACTIV_MADRE=="1"
replace activ_m=1 if ACTIV_MADRE=="2"
replace activ_m=. if ACTIV_MADRE=="9"
replace activ_m=. if ACTIV_MADRE==""
tab ACTIV_MADRE activ_m, m

save "$DAT/workingdata_age_work.dta", replace

*-------------------------------------------------------------------------------
*--- Main Tables
*-------------------------------------------------------------------------------

*** Table 1: Summary Statistics - All Births ***
use "$DAT/workingdata_age_work.dta", clear

*Labels
lab var SEMANAS     "Gestation Weeks"
lab var sem32       "$\geq$ 32 Gestation Weeks"
lab var PESO        "Birth Weight in Grams"
lab var vlbw        "Birth Weight < 1,500"
lab var TALLA       "Birth Length in cm"
lab var dead_at_a00 "Death Within 1st Year of Birth"   
lab var days_y00    "Days Spent in Hospital by Year 1" 
lab var nadmssn_y00 "Number of Admissions to Hospital by Year 1"
lab var EDAD_MADRE  "Mother's Age"
lab var edmom       "Mother's Education Years"

* Panel A:
preserve
keep if  EDAD_MADRE>=15 & EDAD_MADRE<=45
#delimit ;
estpost sum SEMANAS sem32 PESO vlbw TALLA dead_at_a00 days_y00 nadmssn_y00 EDAD_MADRE edmom;
estout using "$TAB/Summary_G1.tex", replace label style(tex)
cells("count(fmt(%15.0gc)) mean(fmt(2)) sd(fmt(2)) min(fmt(2)) max(fmt(2))")
collabels(, none) mlabels(, none);
#delimit cr
restore

* Panel B:
preserve
gen	g2smpl = (mrg_mbdata2main==3 & EDAD_MADRE>=15 & EDAD_MADRE<=45 & ANO_NAC>=2007)
keep if g2smpl==1
#delimit ;
estpost sum SEMANAS sem32 PESO vlbw TALLA dead_at_a00 EDAD_MADRE edmom;
estout using "$TAB/Summary_G2.tex", replace label style(tex)
cells("count(fmt(%15.0gc)) mean(fmt(2)) sd(fmt(2)) min(fmt(2)) max(fmt(2))")
collabels(, none) mlabels(, none);
#delimit cr
restore
estimates clear


*** Table 2: Intensive Health Investments and Birth Outcomes of the Second Generation ***
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
 if ("`lab'"=="imr_ata00") local pv = normal(_b[Robust]/_se[Robust])
 if ("`lab'"=="gw36")      local pv = normal(_b[Robust]/_se[Robust])
 if ("`lab'"=="bw1500")    local pv = normal(_b[Robust]/_se[Robust])
 if ("`lab'"=="gwks")      local pv = 1-normal(_b[Robust]/_se[Robust])
 if ("`lab'"=="peso")      local pv = 1-normal(_b[Robust]/_se[Robust])
 if ("`lab'"=="talla")     local pv = 1-normal(_b[Robust]/_se[Robust])
 if ("`lab'"=="fgrate")    local pv = 1-normal(_b[Robust]/_se[Robust])
 if ("`lab'"=="aindex")    local pv = 1-normal(_b[Robust]/_se[Robust])
 
 local p_value = string(`pv', "%05.3f")
 estadd local p_value "[`p_value']": `outcome'o32
  
estadd scalar Hopt = e(h_l)
estadd scalar Nl   = e(N_h_l)
estadd scalar Nr   = e(N_h_r)
local ef = e(N_h_l) + e(N_h_r)
estadd scalar effopt = `ef'
sum `outcome' if sem32m == 1 & g2_`outcome'_obw_o32 == 1 & `outcome' !=. & PESO_MADRE !=.
estadd scalar dvmean = r(mean)

estimates restore `outcome'o32
estimates save "$OUT/g2_`lab'_blnctrls_o32_hps.ster", replace
}

clear 
local end _blnctrls
gen Robust =.
lab var Robust "Birth weight $<$ 1,500"
gen po32 = .
set obs 7
local i = 1
foreach outcome in gwks peso talla imr_ata00 gw36 bw1500 fgrate { 
 estimates use "$OUT/g2_`outcome'`end'_o32_hps.ster"
 eststo `outcome'o32
 replace po32=e(pv_rb) in `i'
 local ++i
}
rename po32 pval
qsharpenedp pval

local l = 1
foreach outcome in gwks peso talla imr_ata00 gw36 bw1500 fgrate {
    sum bky06 in `l'
    estadd scalar qpv = r(mean): `outcome'o32
    local ++l
}

#delimit ;
esttab gwkso32 pesoo32 tallao32 imr_ata00o32 using "$TAB\T2A_o32.tex",
stats(p_value dvmean N h_r effopt Nl Nr qpv,
fmt(%05.3f %8.3f %9.0gc %5.1f %8.0gc %8.0gc %8.0gc %05.3f)
labels(" " "\\ Mean of Dep.\ Var." "Observations" "Optimal Bandwidth"
       "Effective Observations" "Observations (left)"
       "Observations (right)" "q-sharpened p-value")) keep("Robust") label
b(%-9.3f) se(%-9.3f) noobs nonotes nogaps mlabels(, none) nonumbers style(tex)
fragment replace noline starlevel("*" 0.10 "**" 0.05 "***" 0.01);
    
esttab gw36o32 bw1500o32 fgrateo32 aindexo32 using "$TAB\T2B_o32.tex",
stats(p_value dvmean N h_r effopt Nl Nr qpv,
fmt(%05.3f %8.3f %9.0gc %5.1f %8.0gc %8.0gc %8.0gc %05.3f)
labels(" " "\\ Mean of Dep.\ Var." "Observations" "Optimal Bandwidth"
       "Effective Observations" "Observations (left)"
       "Observations (right)" "q-sharpened p-value")) keep("Robust") label
b(%-9.3f) se(%-9.3f) noobs nonotes nogaps mlabels(, none) nonumbers style(tex)
fragment replace noline starlevel("*" 0.10 "**" 0.05 "***" 0.01);
#delimit cr
estimates clear


*** Table 3: Parental Responses to Treatment Receipt ***

* Panel A:
use "$DAT/workingdata_age_work.dta", clear

gen g1smpl = EDAD_MADRE >= 15 & EDAD_MADRE <= 45 ///
	& ANO_NAC >= 2001 ///
	& mrg_mbdata2NAC == 1 ///
	& nbirths_by_mother > birth_order_by_mother_t ///
	& g1ok_blnctrls == 1
	
keep if g1smpl == 1

foreach var of varlist leavlf_m joinlf_m next_hi_MADRE next_ryai_MADRE ///
                       leavlf_p joinlf_p next_hi_PADRE next_ryai_PADRE {

 rdbwselect `var' PESO if sem32 == 1, c(1499.9) ///
             covs(${g1blnctrls} dbw1??0) vce(cluster PESO)
 scalar RBW_o32 = e(h_mserd)
 scalar LBW_o32 = e(h_mserd)
 scalar maxwt_o32 = 1500 + RBW_o32
 scalar minwt_o32 = 1500 - LBW_o32
 gen g1_`var'_obw_o32 = PESO >= minwt_o32 & PESO <= maxwt_o32


 rdrobust `var' PESO if sem32 == 1, c(1500) scalepar(-1) all ///
           covs(${g1blnctrls} dbw1??0) vce(cluster PESO)
 eststo `var'o32

 local pv = normal(_b[Robust]/_se[Robust])
 local p_value = string(`pv', "%05.3f")
 estadd local p_value "[`p_value']": `var'o32 

 estadd scalar Hopt = e(h_l)
 estadd scalar Nl   = e(N_h_l)
 estadd scalar Nr   = e(N_h_r)
 local ef = e(N_h_l) + e(N_h_r)
 estadd scalar effopt = `ef'
 sum `var' if sem32 == 1 & g1_`var'_obw_o32 == 1 & `var'!= . & PESO != .
 estadd scalar dvmean = r(mean)

 estimates restore `var'o32
 estimates save "$OUT/g1_`var'_blnctrls_o32_hps.ster", replace
}

clear 
local end _blnctrls
gen Robust =.
lab var Robust "Birth weight $<$ 1,500"
gen po32 = .
set obs 8
local i = 1
foreach outcome in leavlf_m joinlf_m next_hi_MADRE next_ryai_MADRE ///
                   leavlf_p joinlf_p next_hi_PADRE next_ryai_PADRE {
 estimates use "$OUT/g1_`outcome'`end'_o32_hps.ster"
 eststo `outcome'o32
 replace po32=e(pv_rb) in `i'
 local ++i
}
rename po32 pval
qsharpenedp pval

local l = 1
foreach outcome in leavlf_m joinlf_m next_hi_MADRE next_ryai_MADRE ///
                   leavlf_p joinlf_p next_hi_PADRE next_ryai_PADRE {
    sum bky06 in `l'
    estadd scalar qpv = r(mean): `outcome'o32
    local ++l
}

#delimit ;
esttab leavlf_mo32 joinlf_mo32 next_hi_MADREo32 next_ryai_MADREo32 
       leavlf_po32 joinlf_po32 next_hi_PADREo32 next_ryai_PADREo32
       using "$TAB\T3A_o32_hps_edit.tex",
stats(p_value dvmean N h_r effopt Nl Nr qpv,
fmt(%05.3f %15.3f %15.0gc %5.1f %8.0gc %8.0gc %8.0gc %05.3f)
labels(" " "\\ Mean of Dep.\ Var." "Observations" "Optimal Bandwidth"
       "Effective Observations" "Observations (left)"
       "Observations (right)" "q-sharpened p-value")) keep("Robust") label
b(%-9.3f) se(%-9.3f) noobs nonotes nogaps mlabels(, none) nonumbers style(tex)
fragment replace noline starlevel("*" 0.10 "**" 0.05 "***" 0.01);
#delimit cr

* Panel B
use "$DAT/workingdata_age_work.dta", clear

foreach var of varlist days_y02_ispr days_y04_ispr days_y06_ispr ///
					   days_y02_priv days_y04_priv days_y06_priv {
					   
 preserve					
 gen g1smpl = ANO_NAC >= 2001 ///
   & EDAD_MADRE >= 15 & EDAD_MADRE <= 45 ///
   & mrg_mbdata2NAC == 1 ///
   & g1ok_blnctrls == 1
	
 keep if g1smpl == 1
 
 * Maximum horizon:
 sum ANO_NAC
 scalar max_ano_nac = r(max)
 scalar min_ano_nac = r(min)
 scalar max_age = max_ano_nac - min_ano_nac
 
 local aa = substr(word(subinstr("`var'", "_", " ", .), 2), -2, 2)
 local a = real("`aa'")
 
 gen share_`var' = `var' / days_y`aa' if `var' != . & days_y`aa' > 0 & days_y`aa' != .
 replace share_`var' = 0 if `var' == . & days_y`aa' > 0 & days_y`aa' != .
	
 if `a' >= 0 & `a' <= `=max_age' {
 
 rdbwselect share_`var' PESO if sem32 == 1, c(1500) ///
            covs(${g1blnctrls} dbw1??0) vce(cluster PESO)
 scalar RBW_o32 = e(h_mserd)
 scalar LBW_o32 = e(h_mserd)
 scalar maxwt_o32 = 1500 + RBW_o32
 scalar minwt_o32 = 1500 - LBW_o32
 gen g1_`var'_obw_o32 = PESO >= minwt_o32 & PESO <= maxwt_o32

 rdrobust share_`var' PESO if sem32 == 1, c(1500) scalepar(-1) all ///
          covs(${g1blnctrls} dbw1??0) vce(cluster PESO)
 eststo `var'o32
 
 local pv = normal(_b[Robust]/_se[Robust])
 local p_value = string(`pv', "%05.3f")
 estadd local p_value "[`p_value']": `var'o32

 estadd scalar Hopt = e(h_l)
 estadd scalar Nl   = e(N_h_l)
 estadd scalar Nr   = e(N_h_r)
 local ef = e(N_h_l) + e(N_h_r)
 estadd scalar effopt = `ef'
 sum share_`var' if sem32 == 1 & g1_`var'_obw_o32 == 1 & share_`var'!= . & PESO != .
 estadd scalar dvmean = r(mean)

 estimates restore `var'o32
 estimates save "$OUT/g1_`var'_blnctrls_o32_hps.ster", replace
 
 }
 restore
}

foreach var of varlist another_birth bspcngf {

 preserve					
 gen g1smpl = EDAD_MADRE >= 15 & EDAD_MADRE <= 45 ///
	& ANO_NAC >= 2001     ///
	& ANO_NAC <= 2017     ///
	& mrg_mbdata2NAC == 1 ///
	& g1ok_blnctrls == 1
  
 keep if g1smpl == 1

 rdbwselect `var' PESO if sem32 == 1, c(1500) ///
             covs(${g1blnctrls} dbw1??0) vce(cluster PESO)
 scalar RBW_o32 = e(h_mserd)
 scalar LBW_o32 = e(h_mserd)
 scalar maxwt_o32 = 1500 + RBW_o32
 scalar minwt_o32 = 1500 - LBW_o32
 gen g1_`var'_obw_o32 = PESO >= minwt_o32 & PESO <= maxwt_o32


 rdrobust `var' PESO if sem32 == 1, c(1500) scalepar(-1) all ///
           covs(${g1blnctrls} dbw1??0) vce(cluster PESO)
 eststo `var'o32
 
 local pv = normal(_b[Robust]/_se[Robust])
 local p_value = string(`pv', "%05.3f")
 estadd local p_value "[`p_value']": `var'o32

 estadd scalar Hopt = e(h_l)
 estadd scalar Nl   = e(N_h_l)
 estadd scalar Nr   = e(N_h_r)
 local ef = e(N_h_l) + e(N_h_r)
 estadd scalar effopt = `ef'
 sum `var' if sem32 == 1 & g1_`var'_obw_o32 == 1 & `var'!= . & PESO != .
 estadd scalar dvmean = r(mean)

 estimates restore `var'o32
 estimates save "$OUT/g1_`var'_blnctrls_o32_hps.ster", replace
 restore
}

clear 
local end _blnctrls
gen Robust =.
lab var Robust "Birth weight $<$ 1,500"
gen po32 = .
set obs 8
local i = 1
foreach outcome in days_y02_ispr days_y04_ispr days_y06_ispr ///
				   days_y02_priv days_y04_priv days_y06_priv ///
				   another_birth bspcngf {
 estimates use "$OUT/g1_`outcome'`end'_o32_hps.ster"
 eststo `outcome'o32
 replace po32=e(pv_rb) in `i'
 local ++i
}
rename po32 pval
qsharpenedp pval

local l = 1
foreach outcome in days_y02_ispr days_y04_ispr days_y06_ispr ///
				   days_y02_priv days_y04_priv days_y06_priv ///
				   another_birth bspcngf {
    sum bky06 in `l'
    estadd scalar qpv = r(mean): `outcome'o32
    local ++l
}


#delimit ;
esttab days_y02_ispro32 days_y04_ispro32 days_y06_ispro32 
	   days_y02_privo32 days_y04_privo32 days_y06_privo32 
	   another_birtho32  bspcngfo32 
       using "$TAB\T3B_o32_hps_edit.tex",
stats(p_value dvmean N h_r effopt Nl Nr qpv,
fmt(%05.3f %15.3f %15.0gc %5.1f %8.0gc %8.0gc %8.0gc %05.3f)
labels(" " "\\ Mean of Dep.\ Var." "Observations" "Optimal Bandwidth"
       "Effective Observations" "Observations (left)"
       "Observations (right)" "q-sharpened p-value")) keep("Robust") label
b(%-9.3f) se(%-9.3f) noobs nonotes nogaps mlabels(, none) nonumbers style(tex)
fragment replace noline starlevel("*" 0.10 "**" 0.05 "***" 0.01);
#delimit cr
estimates clear
kñkñk

*** Table 4: The Effect of Compositional Change on Health at Birth ***
use "$DAT/workingdata_age_work.dta", clear
gen premature = SEMANAS<37
gen fgrate = PESO/SEMANAS

bys ID_MADRE: gen N=_N
drop if N>10
sum SEMANAS PESO TALLA dead_at_a00 premature vlbw fgrate
foreach y of varlist SEMANAS PESO TALLA dead_at_a00 premature vlbw fgrate {
    reghdfe `y' i.EDAD_MADRE, absorb(ID_MADRE)
    margins i.EDAD_MADRE, post level(99)
    foreach age of numlist 15(1)26 {
        local `y'`age' = _b[`age'.EDAD_MADRE] 
    }
}

use "$DAT/workingdata_age_work.dta", clear
gen premature = SEMANAS<37
gen fgrate = PESO/SEMANAS

gen g2smpl = mrg_mbdata2main == 3 ///
	& EDAD_MADRE >= 15 & EDAD_MADRE <= 45 & ANO_NAC >= 2007 & g2ok_blnctrls == 1

keep if g2smpl == 1
gen motherTreated = PESO_MADRE < 1500
foreach age of numlist 15(1)26 {
    gen age`age' = EDAD_MADRE == `age'
    reg age`age' motherTreated  if PESO_MADRE >=1300& PESO_MADRE <= 1700
}
local Deltas
local proportions
foreach y of varlist SEMANAS PESO TALLA dead_at_a00 premature vlbw fgrate {
    rdrobust `y' PESO_MADRE if sem32m == 1, c(1500) scalepar(-1) all ///
       covs(${g2blnctrls} dbw1??0m) vce(cluster PESO_MADRE)
    local Hopt = e(h_l)
    local down = 1500 - `Hopt' 
    local up   = 1500 + `Hopt' 
    local b    = _b[Robust]

    foreach age of numlist 15(1)26 {
        reg age`age' motherTreated if PESO_MADRE >=`down'& PESO_MADRE <= `up'
    }
    preserve
    collapse age15-age26 if PESO_MADRE >=1300& PESO_MADRE<=1700, by(motherTreated)
    reshape long age, i(motherTreated) j(a)
    rename (age a) (prop age)
    gen y = .
    foreach age of numlist 15(1)26 {
        replace y = ``y'`age'' if age == `age'
    }
    gen composition = y * prop 
    collapse (sum) composition , by(motherTreated)
    list
    local ceffect = composition[2] - composition[1] 
    dis "ceffect is `ceffect'"
    autofmt, input(`ceffect') dec(3)
    local D`y' = r(output1)
    local pcnt = `D`y''/`b'*100
    dis "peffect is `pcnt'"
    autofmt, input(`pcnt') dec(3)
    local peffect = r(output1)

    local Deltas      `Deltas' & `D`y''
    local proportions `proportions' & `peffect'
    dis "`Deltas'"
    dis "`proportions'"


    restore
}

// Panel A
dis "`Deltas'"
dis "`proportions'"

use "$DAT/workingdata_age_work.dta", clear
gen premature = SEMANAS<37
gen fgrate = PESO/SEMANAS

bys ID_MADRE: gen N=_N
drop if N>10

gen difedad = abs(EDAD_PADRE-EDAD_MADRE)
gen highEd  = highest_ed3lvls_mp==2|highest_ed3lvls_mp==3
gen g1 = difedad<=5
gen g2 = difedad>0
gen g3 = highEd==0
gen g4 = highEd==1

foreach y of varlist SEMANAS PESO TALLA dead_at_a00 premature vlbw fgrate {
    foreach g of varlist g1 g2 g3 g4 {
        reghdfe `y' i.EDAD_MADRE if `g'==1, absorb(ID_MADRE)
        margins i.EDAD_MADRE if `g'==1, post level(99)
        foreach age of numlist 15(1)26 {
            local `y'`age'_`g' = _b[`age'.EDAD_MADRE] 
        }
    }
}

use "$DAT/workingdata_age_work.dta", clear
gen premature = SEMANAS<37
gen fgrate = PESO/SEMANAS

gen g2smpl = mrg_mbdata2main == 3 ///
	& EDAD_MADRE >= 15 & EDAD_MADRE <= 45 & ANO_NAC >= 2007 & g2ok_blnctrls == 1

keep if g2smpl == 1
gen motherTreated = PESO_MADRE < 1500

gen difedad = abs(EDAD_PADRE-EDAD_MADRE)
gen relQual = difedad<=5
gen highEd  = highest_ed3lvls_aaoo==2|highest_ed3lvls_aaoo==3
gen g1 = difedad<=5
gen g2 = difedad>0
gen g3 = highEd==0
gen g4 = highEd==1


foreach age of numlist 15(1)26 {
    gen age`age' = EDAD_MADRE == `age'
    reg age`age' motherTreated if PESO_MADRE >=1300& PESO_MADRE <= 1700
    reg age`age' motherTreated if highEd==1 & PESO_MADRE >=1300& PESO_MADRE <= 1700
    reg age`age' motherTreated if highEd==0 & PESO_MADRE >=1300& PESO_MADRE <= 1700
}
local DeltasEd
local proportionsEd
local DeltasPart
local proportionsPart
local DeltasBoth
local proportionsBoth
foreach y of varlist SEMANAS PESO TALLA dead_at_a00 premature vlbw fgrate {
    rdrobust `y' PESO_MADRE if sem32m == 1, c(1500) scalepar(-1) all ///
       covs(${g2blnctrls} dbw1??0m) vce(cluster PESO_MADRE)
    local Hopt = e(h_l)
    local down = 1500 - `Hopt' 
    local up   = 1500 + `Hopt' 
    local b    = _b[Robust]

    // EDUC
    preserve
    collapse age15-age26 if PESO_MADRE >=1300& PESO_MADRE<=1700, by(motherTreated highEd)
    reshape long age, i(motherTreated highEd) j(a)
    rename (age a) (prop age)
    gen g3 = highEd==0
    gen g4 = highEd==1
    gen y = .
    foreach age of numlist 15(1)26 {
        replace y = ``y'`age'_g3' if age == `age'&g3==1
        replace y = ``y'`age'_g4' if age == `age'&g4==1
    }
    gen composition = y * prop 
    collapse (sum) composition , by(motherTreated)
    list
    local ceffect = composition[2] - composition[1] 
    dis "ceffect is `ceffect'"
    autofmt, input(`ceffect') dec(3)
    local D`y' = r(output1)
    local pcnt = `D`y''/`b'*100
    dis "peffect is `pcnt'"
    autofmt, input(`pcnt') dec(3)
    local peffect = r(output1)

    local DeltasEd      `DeltasEd' & `D`y''
    local proportionsEd `proportionsEd' & `peffect'
    restore

    // PARTNERSHIP
    preserve
    collapse age15-age26 if PESO_MADRE >=1300& PESO_MADRE<=1700, by(motherTreated relQual)
    reshape long age, i(motherTreated relQual) j(a)
    rename (age a) (prop age)
    gen g1 = relQual==1
    gen g2 = relQual==0
    gen y = .
    foreach age of numlist 15(1)26 {
        replace y = ``y'`age'_g1' if age == `age'&g1==1
        replace y = ``y'`age'_g2' if age == `age'&g2==1
    }
    gen composition = y * prop 
    collapse (sum) composition , by(motherTreated)
    list
    local ceffect = composition[2] - composition[1] 
    dis "ceffect is `ceffect'"
    autofmt, input(`ceffect') dec(3)
    local D`y' = r(output1)
    local pcnt = `D`y''/`b'*100
    dis "peffect is `pcnt'"
    autofmt, input(`pcnt') dec(3)
    local peffect = r(output1)

    local DeltasPart      `DeltasPart' & `D`y''
    local proportionsPart `proportionsPart' & `peffect'
    restore

    // BOTH
    preserve
    collapse age15-age26 if PESO_MADRE >=1300& PESO_MADRE<=1700, by(motherTreated highEd relQual)
    reshape long age, i(motherTreated highEd relQual) j(a)
    rename (age a) (prop age)
    gen g1 = relQual==1
    gen g2 = relQual==0
    gen g3 = highEd==0
    gen g4 = highEd==1
    gen y = .
    foreach age of numlist 15(1)26 {
        replace y = ``y'`age'_g1' if age == `age'&g1==1
        replace y = ``y'`age'_g2' if age == `age'&g2==1
        replace y = ``y'`age'_g3' if age == `age'&g3==1
        replace y = ``y'`age'_g4' if age == `age'&g4==1
    }
    gen composition = y * prop 
    collapse (sum) composition , by(motherTreated)
    list
    local ceffect = composition[2] - composition[1] 
    dis "ceffect is `ceffect'"
    autofmt, input(`ceffect') dec(3)
    local D`y' = r(output1)
    local pcnt = `D`y''/`b'*100
    dis "peffect is `pcnt'"
    autofmt, input(`pcnt') dec(3)
    local peffect = r(output1)

    local DeltasBoth      `DeltasBoth' & `D`y''
    local proportionsBoth `proportionsBoth' & `peffect'
    dis "`DeltasBoth'"
    dis "`proportionsBoth'"
    restore
}

// Panel B-D
dis "`DeltasEd'"
dis "`proportionsEd'"

dis "`DeltasPart'"
dis "`proportionsPart'"

dis "`DeltasBoth'"
dis "`proportionsBoth'"

*-------------------------------------------------------------------------------
*--- Main Figures
*-------------------------------------------------------------------------------

*** Figure 1: Intergenerational Transmission of Early Life Health Outcomes ***
use "$DAT/MERGED_DEF2NAC.dta", clear

gen LBW       = NAC_PESO<2500

keep if NAC_PESO_MADRE>=1000&NAC_PESO_MADRE<=5000

* Panel A: 
#delimit ;
binsreg NAC_PESO NAC_PESO_MADRE, nbins(30) dots(0,0) line(3,3) ci(3,3) cb(3,3)
ytitle("Child's Weight at Birth") xtitle("Mother's Weight at Birth")
dotsplotopt(ms(Oh) mcolor(black)) ciplotopt(lcolor(gs3)) lineplotopt(lcolor(lavender))
cbplotopt(color(blue%40)) level(99) scheme(plotplain);
graph export "$FIG/birthWeightIntergen.pdf", replace ;
#delimit cr

* Panel B: 
format LBW %05.2f
#delimit ;
binsreg LBW NAC_PESO_MADRE, nbins(30) dots(0,0) line(3,3) ci(3,3) cb(3,3)
ytitle("Child Low Birth Weight") xtitle("Mother's Weight at Birth")
dotsplotopt(ms(Oh) mcolor(black)) ciplotopt(lcolor(gs3)) lineplotopt(lcolor(lavender))
cbplotopt(color(blue%40))  ylabel(, format(%04.2f)) level(99) scheme(plotplain);
graph export "$FIG/lbwIntergen.pdf", replace;
#delimit cr
estimates clear


*** Figure 2: Descriptive Plots of Parental Policy Receipt and Child Health Measures ***
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
 
 keep if SEMANAS_MADRE>=32& SEMANAS_MADRE!=.
    
 #delimit ;
 rdrobust `var' PESO_MADRE if sem32m == 1, c(1500) scalepar(-1) all
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
 graph export "$FIG/RDplot_`var'.eps", replace
 restore
}
estimates clear


*** Figure 3: Distributional Inmpacts of Early Life Health Interventions on Second Generation Health ***
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
 
  rdrobust `lab'_below_`wt' PESO_MADRE if sem32m == 1, c(1500) scalepar(-1) all ///
            covs(${g2blnctrls} dbw1??0m) vce(cluster PESO_MADRE)
		  
  eststo `outcome'o32_`wt'
  estimates save "$OUT/g2_`lab'`wt'_o32_hps.ster", replace
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
  estimates use "$OUT/g2_bw`wt'_o32_hps.ster"
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
graph export "$FIG/DistBW_o32.pdf", replace

* Panel B:
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
  estimates use "$OUT/g2_gw`wt'_o32_hps.ster"
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
graph export "$FIG/DistWeeks_o32.pdf", replace
estimates clear


*** Figure 4: Long-Term Health Stocks and Early Life Interventions ***

* Panel A:
use "$DAT/workingdata_age_work.dta", clear

gen g1smpl = EDAD_MADRE >= 15 & EDAD_MADRE <= 45 ///
	& mrg_mbdata2NAC == 1 ///
	& g1ok_blnctrls == 1

keep if g1smpl == 1

* Maximum horizon:
sum ANO_NAC
scalar max_ano_nac = r(max)
scalar min_ano_nac = r(min)
scalar max_age = max_ano_nac - min_ano_nac

drop days_y25 days_y26 days_y27

foreach var of varlist days_y?? {
 local aa = substr("`var'", -2, 2)
 local a = real("`aa'")
	
 gen days0_y`aa' = `var'
 replace days0_y`aa' = 0 if days_y`aa' == . & nadmssn_y`aa' != . & ndischrg_y`aa' != .
 label var days0_y`aa' "Number of days in hospital at age `a' (inputing 0 for no hospitalization)"
	
 if `a' >= 0 & `a' <= `=max_age'  {
  #delimit ;
  rdrobust days0_y`aa' PESO if sem32 == 1, c(1500) scalepar(-1) all
           covs(${g1blnctrls} dbw1??0) vce(cluster PESO);
  #delimit cr
		
 estimates save "$OUT/g1_hdays0at`aa'_o32_hps.ster", replace	
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
 estimates use "$OUT/g1_hdays0at0`age'_o32_hps.ster"
 eststo hosp
   qui: replace EST = _b[Robust] if age==`age'
   qui: replace LB  = _b[Robust] + invnormal(0.025)*_se[Robust] if  age==`age'
   qui: replace UB  = _b[Robust] + invnormal(0.975)*_se[Robust] if  age==`age'
   qui: replace LB90  = _b[Robust] + invnormal(0.05)*_se[Robust] if age==`age'
   qui: replace UB90  = _b[Robust] + invnormal(0.95)*_se[Robust] if age==`age'
}
foreach age of numlist 10(1)24 {
 estimates use "$OUT/g1_hdays0at`age'_o32_hps.ster"
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
graph export "$FIG/hospitalization_o32.eps", replace

* Panel B:
use "$DAT/workingdata_age_work.dta", clear

gen g1smpl = EDAD_MADRE >= 15 & EDAD_MADRE <= 45 ///
	& mrg_mbdata2NAC == 1 ///
	& g1ok_blnctrls == 1

keep if g1smpl == 1

* Maximum horizon:
sum ANO_NAC
scalar max_ano_nac = r(max)
scalar min_ano_nac = r(min)
scalar max_age = max_ano_nac - min_ano_nac

local bb0 = substr(string(min_ano_nac), -2, 2)
local bb1 = substr(string(max_ano_nac), -2, 2)

drop days_y25_notext days_y26_notext days_y27_notext

foreach var of varlist days_y??_notext {
 local aa = substr(word(subinstr("`var'", "_", " ", .), 2), -2, 2)
 local a = real("`aa'")
 
 gen days0_y`aa'_notext = `var'
 replace days0_y`aa'_notext = 0 if days_y`aa' == . & nadmssn_y`aa' != . & ndischrg_y`aa' != .
 label var days0_y`aa'_notext "Number of days in hospital at age `a' (inputing 0 for no hospitalization)"
	
 if `a' >= 0 & `a' <= `=max_age'  {
  #delimit ;
  rdrobust days0_y`aa'_notext PESO if sem32 == 1, c(1500) scalepar(-1) all
           covs(${g1blnctrls} dbw1??0) vce(cluster PESO);
  #delimit cr
		
  estimates save "$OUT/g1_hdaysnex0at`aa'_o32_hps.ster", replace	
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
 estimates use "$OUT/g1_hdaysnex0at0`age'_o32_hps.ster"
 eststo hosp
   qui: replace EST = _b[Robust] if age==`age'
   qui: replace LB  = _b[Robust] + invnormal(0.025)*_se[Robust] if  age==`age'
   qui: replace UB  = _b[Robust] + invnormal(0.975)*_se[Robust] if  age==`age'
   qui: replace LB90  = _b[Robust] + invnormal(0.05)*_se[Robust] if age==`age'
   qui: replace UB90  = _b[Robust] + invnormal(0.95)*_se[Robust] if age==`age'
}
foreach age of numlist 10(1)24 {
 estimates use "$OUT/g1_hdaysnex0at`age'_o32_hps.ster"
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
graph export "$FIG/hospitalizationNoExt_o32.eps", replace
estimates clear


*** Figure 5: Fertility and Stocks of Health at Birth ***
use "$DAT/workingdata_age_work.dta", clear

*Panel A-D:
gen weightBin=ceil((PESO-1500)/50)
foreach num of numlist 16 19 22 25 {
 local ytit "Birth by age `num'"
 local svnm "childBy"
 
 * Any birth by age:
 gen had_birth_by_`num' = (nbirths_by_`num' > 0) if nbirths_by_`num'!= .
 label var had_birth_by_`num' "Any birth by age `num'"

 preserve
 keep if had_birth_by_`num'!=.
 gen NUM=1
 keep if weightBin>-18 & weightBin<66
 collapse had_birth_by_`num' (sum) NUM, by(weightBin)
 tab NUM
 keep if NUM>3
 replace weightBin = weightBin*50+1500
 #delimit ;
 twoway scatter had_birth_by_`num' weightBin [fw=NUM],
 scheme(plotplainblind) xline(1525) ylabel(, format(%04.2f))
 xtitle("Birth weight bin") ytitle("`ytit'");
 #delimit cr
 graph export "$FIG/`svnm'`num'.eps", replace
 restore
}
estimates clear


*** Figure 6: Impacts of Early Life Health Interventions on Fertility and Spontaneous Abortions ***

* Panel A:
use "$DAT/workingdata_age_work.dta", clear

gen g1smpl = SEXO == 2 ///
	& EDAD_MADRE >= 15 & EDAD_MADRE <= 45 ///
	& mrg_mbdata2NAC == 1 ///
	& g1ok_blnctrls == 1
	
keep if g1smpl == 1
keep if SEMANAS>=32& SEMANAS!=.

foreach y of numlist 15/26 {
 
 #delimit ;
  rdrobust nbirths_by_`y' PESO, c(1500) scalepar(-1) all
  covs(${g1blnctrls} dbw1??0) vce(cluster PESO);
 #delimit cr

 estimates save "$OUT/g1f_nbirths_by`y'_blnctrls_o32_hps.ster", replace
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
 estimates use "$OUT/g1f_nbirths_by`age'_blnctrls_o32_hps.ster"
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
graph export "$FIG/fertility_o32.eps", replace

* Panel B:
use "$DAT/workingdata_age_work.dta", clear

gen g1smpl = SEXO == 2 ///
	& EDAD_MADRE >= 15 & EDAD_MADRE <= 45 ///
	& mrg_mbdata2NAC == 1 ///
	& g1ok_blnctrls == 1
	
keep if g1smpl == 1
keep if SEMANAS>=32& SEMANAS!=.

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
graph export "$FIG/fertility2_o32.eps", replace

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

keep if SEMANAS>=32& SEMANAS!=.

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

 estimates save "$OUT/g1f_anyabrtn_by_y`y'_blnctrls_o32_hps.ster", replace
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
 estimates use "$OUT/g1f_anyabrtn_by_y`age'_blnctrls_o32_hps.ster"
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
xlabel(15(1)26) ytitle("Any Abortion");
#delimit cr
graph export "$FIG/abortions_o32.eps", replace

* Panel D:
use "$DAT/workingdata_age_work.dta", clear

keep if SEMANAS_MADRE>=32& SEMANAS_MADRE!=.
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
graph export "$FIG/sexratio_o32.eps", replace
estimates clear


*** Figure 7: Is There Policy-Driven Selection Into Childbirth ***
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
	
    preserve

    #delimit ;
    rdrobust `var' PESO_MADRE, c(1500) scalepar(-1) all
    covs(${g2blnctrls} dbw1??0m) vce(cluster PESO_MADRE);
    #delimit cr
    local opt = e(h_l)
    display `opt'
    gen cond = PESO_MADRE-1500
    keep if abs(cond)<=`opt'
    gen gram10 = floor(PESO_MADRE/20)*20
	sum PESO_MADRE
	local min = floor(r(min)/20)*20
	if mod(`opt', 20)<10 replace gram10 = gram10+20 if gram10==`min'
    bys gram10: egen meanWt = mean(`var')
    bys gram10: gen N=_n
    bys gram10: gen wt=_N
    
    local farea fcolor(gs12%50) clcolor(orange) clstyle(solid)
    #delimit ;
    twoway qfitci `var' PESO_MADRE if PESO_MADRE<1500, `farea'
    || qfitci `var' PESO_MADRE if PESO_MADRE>1500, `farea'
    || scatter meanWt gram10 [aw=wt] if N==1, ms(Oh) mcolor(blue)
	ylabel(, format(`ff'))
    xline(1500, lcolor(red)) scheme(plotplain) legend(off)
    ytitle("`xx'") xtitle("Weight in grams");
    #delimit cr
    graph export "$FIG/figure7_RDplot_`var'.eps", replace
    restore
}
estimates clear

*-------------------------------------------------------------------------------
*--- Log Close
*-------------------------------------------------------------------------------

log close
