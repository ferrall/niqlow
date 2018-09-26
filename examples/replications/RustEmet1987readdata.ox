#include "oxstd.h"

/** Labels for the head lines for each bus.**/
enum{inid,r1m,r1y,r2m,r2y,r2mi,r3m,r3y,r3mi,fm,fy,Nheader}

/** Code for columns of ouptut data. **/
enum{id,mm,i,x,Nrows}

main () {
decl f = shape(loadmat("a530875.mat",1),128,37),i,
	k,bus,myid,repm,prev,Ntot=0,smp;
decl fo = fopen("RustEmet1987.data","w");

smp = <>;
for(i=0;i<columns(f);++i) {
	myid = f[inid][i];
	prev = Nheader+1;	 //why extra observation??
	println("** ",i);
	if (f[r2m][i]) { // one rebuild occurred
		repm = prev+12*(f[r2y][i]-f[r1y][i]) + (f[r2m][i]-f[r1m][i])-1;
		bus = myid~f[prev:repm][i]~(zeros(repm-prev,1)|1);
		if (f[r3m][i]) { // a second occurred
			prev = repm;
			repm += 12*(f[r3y][i]-f[r2y][i]) + (f[r3m][i]-f[r2m][i])-1;
			bus |= myid~(f[prev+1:repm][i]-f[r2mi][i])~(zeros(repm-prev-1,1)|1);
			bus |= myid~(f[repm+1:][i]-f[r3mi][i])~0;	//
			}
		else   // only one rebuid
			bus |= myid~(f[repm+1:][i]-f[r2mi][i])~0;
		}
	else {	 // engine was not rebuilt
		bus = f[inid][i]~f[prev:][i]~0;
		}
//	println("%cf",{"%5.0f","%12.0f","%3.0f"},bus);
	Ntot += rows(bus);
	bus ~= bus[][mm];
	for(k=0;k<89;++k)  //convert odometer reading into 5K mileage category.
		bus[][x] = bus[][x].>=5000*k .&& bus[][x].<5000*(k+1) .? k .: bus[][x];
	println(i,bus);
	fprintln(fo,"%v",bus);
	smp |= bus;	  // append bus to the data set
	}
fclose(fo);
println("Total observations= ",Ntot," ",Ntot-37);
savemat("RustEmet1987.dta",smp,{"id","mi","d","x"});
}
