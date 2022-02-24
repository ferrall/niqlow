loc rf RustEmet1987_type1_col3_ALL
loc rp ruspy_all
cap log close
log using Rust1987infinitehorizon, replace
use `rf', clear
format %5.0f _all
format %8.0f mi
sort id s
rename t tniqlow
gen n = _n
sort id n
qui by id: gen t = _n-1
sort id t
merge 1:1 id s using `rp' 
tab d d_ruspy
qui by id : gen dx90_ruspy = x90_ruspy[_n+1] - (1-d_ruspy)*x90_ruspy
tab dx90 if group==4
tab dx90_ruspy if group==4
tab dx90 dx90_ruspy if group==4



