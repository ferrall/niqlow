#include "oxstd.oxh"

main() {
	decl x;
	x=loadmat("aguirregabiria_mira/bus1234.dht");
	x ~= range(1,rows(x))';
	savemat("bus1234x.dta",x,{"busid","model","yr","mth","d","mi","n"});
	println(x[:10][]);
}