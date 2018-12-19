#include "oxstd.oxh"

main() {
	decl x;
	savemat("bus1234.dta",x=loadmat("aguirregabiria_mira/bus1234.dht"));
	println(x[:10][]);
}