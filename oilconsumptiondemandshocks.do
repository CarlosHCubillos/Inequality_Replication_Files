clear all
set more off
set matsize 800

local path0="C:\Users\Carlos H Cubillos\Documents\Javeriana PhD\A. Seminario Investigación\Inequality_CGKS_replication_folder\Workfiles_figure3"
capture cd "`path0'"

local path1="C:\Users\Carlos H Cubillos\Documents\Javeriana PhD\A. Seminario Investigación\Inequality_CGKS_replication_folder\Workfiles_figure3"
capture cd "`path1'"

local path_graphs="C:\Users\Carlos H Cubillos\Documents\Javeriana PhD\A. Seminario Investigación\Inequality_CGKS_replication_folder\Workfiles_figure3"

local ARn=2
local MAn=12

local ymin=-0.005
local ymax=0.01

local counter=1

quiet foreach vvx in  D_SD_LNYBTIMP2_SA D_SD_LNSALARYIMP_SA D_SD_LNTOTALEXP3_SA D_SD_LNCONS_SA ///
				D_GINI_YBTIMP2_SA D_GINI_SALARYIMP_SA D_GINI_TOTALEXP3_SA D_GINI_CONS_SA ///
				P9010_LNYBTIMP2_SA P9010_LNSALARYIMP_SA P9010_LNTOTALEXP3_SA P9010_LNCONS_SA {


noisily di "`vvx'"				

use "C:\Users\Carlos H Cubillos\Documents\Javeriana PhD\A. Seminario Investigación\Inequality_CGKS_replication_folder\replication_folder\workfiles\oilconsumptiondemandshocks.dta", clear

tsset t

local DP="`vvx'"

*** create lags of the shock series
forvalues i=0(1)`MAn' {
	gen sh_rr_L`i'=l`i'.sh_rr
}
*** create lags of the depedent variable
forvalues i=1(1)`ARn' {
	gen `DP'_L`i'=l`i'.`DP'
}

*** create leads of the depedent variable
forvalues i=0(1)20 {
	gen `DP'_F`i'=f`i'.`DP'
}

keep year quarter t time sh_rr* `DP'_*

*** stack data by horizon
reshape long `DP'_F, i(t) j(hor)

tsset hor t

*** interact all lags with horizon dummy
*** lags of the dep variable
forvalues i=0(1)20 {
	forvalues j=1(1)`ARn' {

		gen _`DP'_L`j'_H`i'=`DP'_L`j'*(hor==`i')
	}
}
*** lags of the shock series
forvalues i=0(1)20 {
	forvalues j=0(1)`MAn' {

		gen _sh_rr_L`j'_H`i'=sh_rr_L`j'*(hor==`i')
	}
	
}

*** run the pooled regression
xtscc `DP'_F _sh* _`DP'* , fe


egen tx=fill(1(1)10)
capture drop t0
gen t0=tx-1

*====================================================================
*			construct cumulative impulse response
*====================================================================
capture drop b_
capture drop se_
gen b_=.
gen se_=.

quiet forvalues i=0(1)20 {

	local str0="_sh_rr_L0_H0"
	forvalues j=1(1)`i' {
		local str0="`str0'" + "+_sh_rr_L0_H`j'"
	}
	lincom "`str0'"
	replace b_=r(estimate) if t0==`i'
	replace se_=r(se) if t0==`i'
	
	
}

*====================================================================
*			hypothesis testing for path being euqal to zero
*====================================================================		
quiet forvalues i=0(1)20 {
	local str0="_b[_sh_rr_L0_H0]"
	forvalues j=1(1)`i' {
		local j1=`j'-1
		local str`j'="`str`j1''" + "+_b[_sh_rr_L0_H`j']"
	}	
	 di "`str`i''"
	 di " "
}

local strM=" "
quiet forvalues i=0(1)20 {
	local strM="`strM'"+"`str`i''="
}
local strM="`strM'"+"0"
test "`strM'"
	
quiet {	
	noisily di "F-stat =  " r(F)
	noisily di "p-value = " r(p)
}
local rp=string(r(p),"%4.3f")

capture drop b_hi
capture drop b_lo
gen b_lo=b_+se_
gen b_hi=b_-se_

capture drop b_hi90
capture drop b_lo90
gen b_lo90=b_+1.65*se_
gen b_hi90=b_-1.65*se_

if `counter'==1 {
	local title0="Income"
	local yx="st.dev."
	}
if `counter'==2 {
	local title0="Earnings"
	local yx="st.dev."
	}	
if `counter'==3 {
	local title0="Expenditure"
	local yx="st.dev."
	}	
if `counter'==4 {
	local title0="Consumption"
	local yx="st.dev."
	}	
if `counter'==5 {
	local title0="Income"
	local yx="Gini"
	}
if `counter'==6 {
	local title0="Earnings"
	local yx="Gini"
	}	
if `counter'==7 {
	local title0="Expenditure"
	local yx="Gini"
	}	
if `counter'==8 {
	local title0="Consumption"
	local yx="Gini"
	}	
if `counter'==9 {
	local title0="Income"
	local yx="90-10"
	}
if `counter'==10 {
	local title0="Earnings"
	local yx="90-10"
	}	
if `counter'==11 {
	local title0="Expenditure"
	local yx="90-10"
	}	
if `counter'==12 {
	local title0="Consumption"
	local yx="90-10"
	}	
	

twoway 	(rarea b_lo90 b_hi90 t0 if t0<=20 & hor==0, color(gs14)) ///
		(rarea b_lo b_hi t0 if t0<=20 & hor==0, color(gs6)) ///
		(line b_ t0 if t0<=20 & hor==0, lcolor(black) lwidth(thick)) /// 
		, legend(off) yline(0) ///
		xtitle("") ytitle("`yx'") ///
		subtitle("`title0' (p-val = `rp')") ///
		saving("`path_graphs'\`DP'.gph", replace) ///
		name("`DP'", replace)
	
		
di "`path_graphs'\`DP'.gph"		
local counter=`counter'+1
}
				
graph combine 	D_SD_LNYBTIMP2_SA D_SD_LNSALARYIMP_SA D_SD_LNTOTALEXP3_SA D_SD_LNCONS_SA ///
				D_GINI_YBTIMP2_SA D_GINI_SALARYIMP_SA D_GINI_TOTALEXP3_SA D_GINI_CONS_SA ///
				P9010_LNYBTIMP2_SA P9010_LNSALARYIMP_SA P9010_LNTOTALEXP3_SA P9010_LNCONS_SA ///
	, imargin(tiny) rows(3) saving(Figure3, replace)

