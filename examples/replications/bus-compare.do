loc dno 0
forvalues n = 1(1)4329 {
	 if (dx90[`n']!=dx90EL[`n']) {
	 	loc dno = `dno' + 1
		di "Case `dno'"
		loc l = `n'-2
		loc u = `n'+2
		list id busid s mi dmi x90 d dEL dx90 dx90EL  in `l'/`u'
		di "--------------------------------"
	 }
}