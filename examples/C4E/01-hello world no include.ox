acf(const ma, const ilag);
acos(const ma);
aggregatec(const ma, const istep);
aggregater(const ma, const istep);
any(const ma);
arglist();
asin(const ma);
atan(const ma);
atan2(const my, const mx);

bessel(const ma, const type, const n01);
betafunc(const da, const db, const dx);
binand(const ia, const ib, ...);
bincomp(const ia);
binor(const ia, const ib, ...);
binpop(const ia, ...);
binomial(const n, const k);
binvec(const ia);
binxor(const ia, const ib, ...);

cabs(const ma);
cdiv(const ma, const mb);
ceil(const ma);
cerf(const ma);
cexp(const ma);
chdir(const s);
choleski(const ma);
classname(const obj);
clog(const ma);
clone(const obj, ...);
cmul(const ma, const mb);
columns(const ma);
constant(const dval, const r, ...);
correlation(const ma);
cos(const ma);
cosh(const ma);
countc(const ma, const va);
countr(const ma, const va);
csqrt(const ma);
cumprod(const mfac, ...);
cumsum(const mx, const vp, ...);
cumulate(const ma, ...);

date();
dawson(const ma);
dayofcalendar(...);
dayofeaster(const year);
dayofmonth(const year, const month, const dayofweek, const nth);
dayofweek(const year, ...);
decldl(const ma, const aml, const amd);
decldlband(const ma, const aml, const amd);
declu(const ma, const aml, const amu, const amp);
decqr(const ma, const amht, const amr, const amp);
decqrmul(const mht, ...);
decqrupdate(const amq, const amr, const i1, ...);
decschur(const ma, const amval, const ams, ...);
decschurgen(const ma, const mb, const amalpha, const ambeta, const ams, const amt, ...);
decsvd(const ma, ...);
deletec(const ma, ...);
deleteifc(const ma, const mif);
deleteifr(const ma, const mif);
deleter(const ma, ...);
denschi(const ma, const df);
densf(const ma, const df1, const df2);
densn(const ma);
denst(const ma, const df);
determinant(const ma);
dfft(const mx, ...);
diag(const ma);
diagcat(const ma, const mb);
diagonal(const ma, ...);
diagonalize(const ma);
diff(const ma, ...);
diff0(const ma, ...);
discretize(const vy, const dmin, const dmax, const icount, const ioption);
dropc(const ma, const midx, ...);
dropr(const ma, const midx, ...);

eigen(const ma, const amval, ...);
eigensym(const ma, const amval, ...);
eigensymgen(const ma, const mb, const amval, const amvec);
eprint(const a, ...);
erf(const ma);
exclusion(const ma, const mb, ...);
exit(const iexit);
exp(const ma);
expint(const x);

fabs(const ma);
factorial(const n);
fclose(const file);
feof(const file);
fflush(const file);
fft(const mx, ...);
fft1d(const mx, ...);
find(const where, const what, ...);
findsample(const mdata, const vvarsel, const vlagsel,
    const it1, const it2, const imode, const ait1, const ait2);
floor(const ma);
fmod(const da, const db);
fopen(const sfile, ...);
format(const sfmt);
fprint(const file, const a, ...);
fprintln(const file, const a, ...);
fread(const file, const type, ...);
fremove(const sfile);
fscan(const file, const a, ...);
fseek(const file, ...);
fsize(const file);
ftime(const file);
fuzziness(const deps);
fwrite(const file, const a);

gammafact(const ma);
gammafunc(const dx, const dp);
getcwd();
getenv(const s);
getfiles(const sfilemask);
getfolders(const sfilemask);

hyper_2F1(const a, const b, const c, const z);

idiv(const ia, const ib);
imod(const ia, const ib);
insertc(const ma, const c, const cadd, ...);
insertr(const ma, const r, const radd, ...);
intersection(const ma, const mb, ...);
invert(const ma, ...);
inverteps(const deps);
invertgen(const ma, ...);
invertsym(const ma, ...);
isaddress(const a);
isarray(const a);
isclass(const a, ...);
isdotfeq(const a1, const a2);
isdotinf(const ma);
isdotmissing(const ma);
isdotnan(const ma);
isdouble(const a);
iseq(const a1, const a2);
isfeq(const a1, const a2);
isfile(const a);
isfunction(const a);
isint(const a);
ismatrix(const a);
ismember(const obj, const smember);
ismissing(const ma);
isnan(const ma);
isstatic(const a);
isstring(const a);

lag(const ma, ...);
lag0(const ma, ...);
limits(const mx, ...);
loadmat(const sname, ...);
loadsheet(const sname, ...);
log(const ma);
log10(const ma);
logdet(const ma, const asign);
loggamma(const ma);
lower(const ma);

max(const a, ...);
maxc(const a);
maxcindex(const a);
maxr(const a);
meanc(const ma);
meanr(const ma);
min(const a, ...);
minc(const a);
mincindex(const a);
minr(const a);
moments(const ma, ...);

nans(const r, ...);
norm(const ma, ...);
nullspace(const ma);

ols2c(const my, const mx, const amb, ...);
ols2r(const my, const mx, const amb, ...);
olsc(const my, const mx, const amb, ...);
olsr(const my, const mx, const amb, ...);
ones(const r, ...);
outer(const mx, const ms, ...);
oxfilename(const itype);
oxprintlevel(const ilevel);
oxrunerror(const smsg, ...);
oxversion();
oxversionispro();
oxwarning(const flset);

periodogram(const ma, ...);
polydiv(const ma, const mb, const cc);
polyeval(const ma, const mx);
polygamma(const ma, const n);
polymake(const roots);
polymul(const ma, const mb);
polyroots(const ma, const amroots);
print(const a, ...);
println(const a, ...);
probchi(const ma, const df, ...);
probf(const ma, const df1, const df2, ...);
probn(const ma);
probt(const ma, const df, ...);
prodc(const ma);
prodr(const ma);
pow(const ma, const p);

quanchi(const mp, const df);
quanf(const mp, const df1, const df2);
quann(const mp);
quant(const mp, const df);
quantilec(const ma, ...);
quantiler(const ma, ...);

range(const imin, const imax, ...);
rank(const ma, ...);
ranloopseed(const iloop, const istage);
rann(const r, const c);
ranseed(const iseed);
ranu(const r, const c);
reflect(const ma);
replace(inval, oldval, newval, ...);
reshape(const ma, const r, const c);
reversec(const ma);
reverser(const ma);
round(const ma);
rows(const ma);

savemat(const sname, const ma, ...);
scan(const a, ...);
selectc(const ma, ...);
selectifc(const ma, const mif);
selectifr(const ma, const mif);
selectr(const ma, ...);
selectrc(const ma, const mr, const mc);
setbounds(const ma, const dlo, const dhi);
setdiagonal(const ma, const mdiag);
setlower(const ma, const ml, ...);
setupper(const ma, const mu, ...);
shape(const ma, const r, const c);
sin(const ma);
sinh(const ma);
sizec(const ma);
sizeof(const ma);
sizer(const ma);
sizerc(const ma);
solveldl(const ml, const md, const mb);
solveldlband(const ml, const md, const mb);
solvelu(const ml, const mu, const mp, const mb);
solvetoeplitz(const mr, const cm, const mb, ...);
sortbyc(const ma, const icol);
sortbyr(const ma, const irow);
sortc(const ma);
sortcindex(const ma);
sortr(const ma);
spline(const my, const mx, const alpha, ...);
sprint(const a, ...);
sprintbuffer(const len);
sqr(const ma);
sqrt(const ma);
sscan(const str, const a, ...);
standardize(const ma);
strfind(const as, const s, ...);
strfindr(const as, const s, ...);
strifind(const as, const s, ...);
strifindr(const as, const s, ...);
strlwr(const s);
strtrim(const s);
strupr(const s);
submat(const ma, const r1, const r2, const c1, const c2, ...);
sumc(const ma);
sumr(const ma);
sumsqrc(const ma);
sumsqrr(const ma);
systemcall(const s);

tailchi(const ma, const df);
tailf(const ma, const df1, const df2);
tailn(const ma);
tailt(const ma, const df);
tan(const ma);
tanh(const ma);
thinc(const ma, const c);
thinr(const ma, const r);
time();
timeofday(...);
timer();
timespan(const time, ...);
timestr(const time);
timing(const mdates, ...);
today();
toeplitz(const ma, ...);
trace(const ma);
trunc(const ma);
truncf(const ma);

union(const ma, const mb);
unique(const ma);
unit(const rc, ...);
unvech(const ma);
upper(const ma);

va_arglist();
varc(const ma);
variance(const ma);
varr(const ma);
vec(const ma);
vech(const ma);
vecindex(const ma, ...);
vecr(const ma);
vecrindex(const ma, ...);

zeros(const r, ...);

#pragma _ox_stdlib(0)

enum
{   FALSE, TRUE
};
enum
{   WFL_DECFAILED = 1, WFL_ITMAX = 2, WFL_CONCAT = 4, WR_VECIDXMAT = 8,
	WFL_DETERMINANT = 16, WFL_USER = 64
};
enum
{	SAM_ALLVALID,	SAM_ENDSVALID, 	SAM_ANYVALID
};
enum
{	SCAN_EOF		=   -1,
	SCAN_IDENTIFIER	=	0,
	SCAN_LITERAL	=	1,
	SCAN_SYMBOL		=	2,
	SCAN_SPACE		=	3
};	

main() {
	println("hello world");
	}