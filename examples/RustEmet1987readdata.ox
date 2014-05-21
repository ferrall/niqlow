/** Read bus data from Rust (1987) to be read into niqlow.

This program reads the original ascii data for Group 4 buses in the paper.

The file is a vectorized version of a 128x37 matrix, each column is a bus.

Each column has a header of 11 rows followed by odometer readings for each month.

The header contains up to two engine replacement months and years (after the bus
was brought into service, which is coded as a first replacement.

This data is then used to convert the raw odometer readings into a combination of
the model <var>(i,x)</var> pair, consisting of the rebuild decision and odometer bin category.

The conversion is not fully verified to be accurate.  A slight discrepancy with
the reported log-likelihood could be due to slight differences with the original data.

<dd>The first four rows of the data set look like this:<pre>
         id         mileage        i             x
       5297.0       6299.0      0.00000       1.0000
       5297.0       10479.      0.00000       2.0000
       5297.0       15201.      0.00000       3.0000
</pre>
Note that this is from the output of the program.  The data are saved directly
to the <code>.dta</code> file. The first reading of Bus 5297 had an odometer of
6299, which is the second category (<code>x=1</code>).
Bus 5297 engine was replaced once, and farther down in the data this looks like:
<pre>
       5297.0  1.5256e+005      0.00000       30.000
       5297.0  1.5510e+005       1.0000       31.000
       5297.0       4770.0      0.00000      0.00000 </pre>
The bus is replaced when <code>x=31</code> and the next odometer reading, after
subtracting the reading at replacement (which is in the header at row <code>r2mi</code>,
the sixth row), is 47770 miles more, which means <code>x = 0</code>.
	   </dd>

**/
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
