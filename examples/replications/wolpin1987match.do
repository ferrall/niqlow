infile t      accept        done    hasoffer using  wolpin1997.txt, clear
replace t = t-61
gen s = 1
replace s = s[_n-1]*(1-accept[_n-1]) if _n>1
gen mtch = 1-s if t==0
loc lo 62
loc hi 74
foreach k of num 1/14 {
    replace mtch = (s[`lo']-s[`hi'])/s[`lo'] if _n==`hi'
    loc lo `hi'
    loc hi = `hi'+13
    }
keep if mtch!=.
list t mtch, noobs
save wolpin1987match, replace
