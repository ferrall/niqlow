/** Multinomial Probit using <span class="n">FiveO</span>.

This is an example of how to implement an econometric objective.  FiveO is not designed to handle sample selection and other standard tasks for the user.  However,
built in features of Ox are readily available for those purposes.


**/
#include "MNP.h"

/** Create a new MNP model.
@param fn string, a file to load the data from using Ox <code>Database.Load()</code>
@param Yname string, name or label of the column in the file that contains <var>Y</var><br>
Yname can contain any integers.  MNP will translate the unique sorted values into 0...Jvals-1
@param Xnames a string of the form <q>var1 var2 ... varN</q>
@comments Observations with any missing data are deleted.<br> A constant column is appended at the end of the X matrix as in Stata.<br> Summary statistics are reported.
**/
MNP::MNP(fn,const Yname,const Xnames)	{
	decl j,data, sample;
	BlackBox::BlackBox("Multinomial Probit");
	namearray = varlist(Xnames);
	data = new Database();
	data.Load(fn);
	sample = deleter(data.GetVar({Yname}|namearray));
	delete data;
	Y =sample[][0];
	X = sample[][1:]~1;
	namearray |= "Cons";
	Jvals = unique(Y);
	J = columns(Jvals);
	NvfuncTerms = rows(Y);
	NN = range(0,NvfuncTerms-1);
	nX = columns(X);
	indY = maxcindex( (Y.==Jvals)' )' ;
	betas = array(zeros(nX,1));
	for (j=1;j<J;++j){
		indY ~= (j-1 .<indY[][0]) .? j-1 .: j;
		betas |= new Coefficients("Y="+sprint(Jvals[j]),nX,namearray);
		Parameters(betas[j]);
		}
	MyMoments(Y~X,{Yname}|namearray);
   }

/** Estimate the model using Simplex and BHHH, print results.
**/
MNP::Estimate() {
	decl j,se,bhhh =new BHHH(this);
	bhhh->Iterate(0);
	println("Default value of Y = ","%2.0f",Jvals[0]);
	for (j=1;j<J;++j)  {
		se = vcur.SE[(j-1)*nX:j*nX-1]';
		println("Y = ","%2.0f",Jvals[j],"%c",{"Estimate","Std.Err","t-ratio"},"%r",namearray,
			betas[j].v~se~(betas[j].v./se));
		}
	}

GQMNP::GQMNP(Npts,const fn,const Yname,const Xnames)	 {
	MNP(fn,Yname,Xnames);	
	GQH::Initialize(Npts);
	this.Npts = Npts;
	Encode(0);
	}
	
/**  Compute and return the vector of log-likelihoods at the current parameters.
**/
GQMNP::vfunc() {
	decl j, D, lk, myD;
	for (j=0,D=<>;j<J;++j)	D ~= X*CV(betas[j]);
	myD = selectrc(D,NN,indY[][0]);
	for (j=1,lk=ones(Npts,NvfuncTerms);j<J;++j) lk .*=  probn(GQH::nodes+myD-selectrc(D,NN,indY[][j]) );
	return log(GQH::wght * lk / M_SQRT2PI)' ;
	}

GHKMNP::GHKMNP(R,const iSigma,const fn,const Yname,const Xnames) {
	MNP(fn,Yname,Xnames);
	ghk = new GHK(R,J,0);
	if (isint(iSigma)) {
		sigfree=identity;
	   }
	else if (sizerc(iSigma)==J) {
		sigfree = onlydiag;
		SigLT = new Coefficients("Sigma",iSigma,0);
		Block(SigLT);
		}
	else if (sizerc(iSigma)==J*J) {
		sigfree = lowertriangle;
		SigLT = new Coefficients("Sigma",vech(iSigma),0);
		Block(SigLT);
		}
	else oxrunerror("Initial Variance matrix must be JxJ");
	Encode(0);
	}
	
/**  Compute and return the vector of log-likelihoods at the current parameters.
**/
GHKMNP::vfunc() {
	decl i, j, D, lk,Sigma;
	ranseed(-1);
	for (j=0,D=<>;j<J;++j)	D ~= X*CV(betas[j]);
	Sigma = sigfree==identity
						? unit(J)
						: sigfree==onlydiag
							? diag(SigLT.v)
							: unvech(SigLT.v);	
	for (i=0,lk=zeros(NvfuncTerms,1);i<NvfuncTerms;++i)	lk[i] = ghk->SimProb(indY[i][0],D[i][],Sigma);
	return log(lk) ;
	}
	
