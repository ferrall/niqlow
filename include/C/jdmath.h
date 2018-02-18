/*--------------------------------------------------------------------------
 * jdmath.h - definitions/declarations for mathematics library
 *
 *       Jurgen Doornik 1990-1995
 *
 *--------------------------------------------------------------------------*/

#ifndef INC_JDMATH
#define INC_JDMATH

/*=========================== typedefs =====================================*/
#ifndef INC_JDTYPES
#include "jdtypes.h"
#endif
#ifndef INC_JDMATRIX
#include "jdmatrix.h"			      /* includes default matrix allocators */
#ifdef JDDEBUG_
#include "jddebug.h"		/* replaces matrix allocators by debug versions */
#endif
#endif

#ifdef __cplusplus
extern "C" {
#endif

/*=========================== constants ====================================*/
#ifndef DBL_RADIX
#define DBL_RADIX 2
#endif
/* Threading is designed to be activated when it roughly matches the serial */
/* speed; this factor makes it more conservative */
#define OMP_FACTOR 			2
#define OMP_FACTOR_LINEAR 	8

#define JDRAN_MWC_R  	256			/* George Marsaglia MWC state dimension */
#define JDRAN_ZIGNOR_C 	128						        /* number of blocks */

/*============================ typedefs ====================================*/
struct ranState;									 /* forward declaration */

typedef double 			(JDCALL * DRANFUN)(struct ranState *pRan);
typedef int    			(JDCALL * IRANFUN)(struct ranState *pRan);
typedef unsigned int	(JDCALL * UIRANFUN)(struct ranState *pRan);
typedef MATRIX 			(JDCALL * MATRANFUN)(struct ranState *pRan, MATRIX, int, int);
typedef void   			(JDCALL * RANSETSEEDFUN)(struct ranState *pRan, int *, int);
typedef int    			(JDCALL * RANGETSEEDFUN)(struct ranState *pRan, int *, int);
typedef void   			(JDCALL * RANCLRSEEDFUN)(struct ranState *pRan, struct ranState *pRanInit, int);

struct ranState
{
	int 			iSeedPM;
	unsigned int 	uiSeed1GM, uiSeed2GM;
	unsigned int 	uiSeed1LE, uiSeed2LE, uiSeed3LE, uiSeed4LE;
	unsigned int 	uiStateMWC, uiCarryMWC, auiStateMWC[JDRAN_MWC_R + 1];
	double 			adZigX[JDRAN_ZIGNOR_C + 1], adZigR[JDRAN_ZIGNOR_C];
	BOOL 			fNormalInStore;
	UIRANFUN 		fnUIRanu;
	DRANFUN 		fnDRanu, fnDRann;
	MATRANFUN  		fnMatRanu, fnMatRann;
	RANSETSEEDFUN 	fnRanSetSeed;
	RANGETSEEDFUN 	fnRanGetSeed;
	RANCLRSEEDFUN 	fnRanClrSeed;
};

/*================== function prototypes & macro's =========================*/
#ifndef max
#define max(a,b)  (((a) > (b)) ? (a) : (b))
#endif
#ifndef min
#define min(a,b)  (((a) < (b)) ? (a) : (b))
#endif

#ifndef SIGNTRANSFER                       /* Fortran style sign() function */
#define SIGNTRANSFER(a, b)   ((b) >= 0 ? fabs(a) : -fabs(a))
#endif


/* mtarma.c */
BOOL   JDCALL FGetPartAcf(VECTOR vAcf, int cAcf, VECTOR vPartAcf);
MATRIX JDCALL MatPartAcf(MATRIX mPartAcf, MATRIX mAcf, int cAcf, MATRIX mY, int cY, double *pdLogDet, BOOL bFilter);
BOOL   JDCALL FArmaVar(VECTOR vP, int p, int q, double dVar, VECTOR vSigma, int cT);

/* mtmatrix.c */
void   JDCALL SetInvertEps(double dEps);
double JDCALL DGetInvertEps(void);
double JDCALL DDiagXSXt(int iT, MATRIX mX, MATRIX mS, int cX);
double JDCALL DDiagXtSXtt(int iT, MATRIX mXt, MATRIX mS, int cX);
MATRIX JDCALL MatSumOuter(MATRIX mQ, MATRIX mX, int cT, int cX);
double JDCALL DTrace(MATRIX, int);
double JDCALL DTraceAB(MATRIX, MATRIX, int, int);
double JDCALL DVecsum(VECTOR, int);
double JDCALL DVecXtY(VECTOR vX, VECTOR vY, int cC);
MATRIX JDCALL MatAB(MATRIX, int, int, MATRIX, int, MATRIX);
MATRIX JDCALL MatABt(MATRIX, int, int, MATRIX, int, MATRIX);
MATRIX JDCALL MatAtB(MATRIX, int, int, MATRIX, int, MATRIX);
MATRIX JDCALL MatAdd(MATRIX mA, int cM, int cN, MATRIX mB, double dFac, MATRIX mAplusB);
MATRIX JDCALL MatCopyTranspose(MATRIX mB, MATRIX mA, int cM, int cN);
MATRIX JDCALL MatReflect(MATRIX, int);
MATRIX JDCALL MatTranspose(MATRIX, int);
VECTOR JDCALL VecTranspose(VECTOR vA, int cM, int cN);
MATRIX JDCALL MatS1S2(MATRIX, MATRIX, int);
MATRIX JDCALL MatXtXtt(MATRIX mXt, int cX, int iT1, int iT2, MATRIX mXtX);
MATRIX JDCALL MatXXt(MATRIX mX, int iT1, int iT2, int cX, MATRIX mXtX);
MATRIX JDCALL MatXtZtt(MATRIX mXt, int cX, int iT1, int iT2, MATRIX mZt, int cZ, MATRIX mXtZ);
MATRIX JDCALL MatXXtVec(MATRIX mX, int iT1, int iT2, int cX, VECTOR vY, MATRIX mXtX);
MATRIX JDCALL MatBSBt(MATRIX, int, MATRIX, int, MATRIX);
MATRIX JDCALL MatBtSB(MATRIX, int, MATRIX, int, MATRIX);
MATRIX JDCALL MatBtB(MATRIX mB, int rB, int cB, MATRIX mBtB);
MATRIX JDCALL MatBtBVec(MATRIX mB, int cT, int cB, VECTOR vY, MATRIX mBtB);
MATRIX JDCALL MatBBt(MATRIX mB, int rB, int cB, MATRIX mBBt);
void   JDCALL QSortVecMat(VECTOR, VECTOR, MATRIX, int, int, int);

/* mtband.h */
int    JDCALL ToeplitzSolve(VECTOR, int, int, MATRIX, int, VECTOR v_1);
int    JDCALL ILDLbandDec(MATRIX, VECTOR, int, int);
void   JDCALL LDLbandSolve(MATRIX, VECTOR, VECTOR, VECTOR, int, int);

/* mtbessel.h */
double JDCALL DBessel01(double x, int type, int n);
/* mtbesselfrac.h */
double JDCALL DExpInt(double x);
double JDCALL DExpInte(double x);
double JDCALL DExpInt1(double x);
double JDCALL DDawson(double x);
double JDCALL DBesselNu(double x, int type, double nu);
double JDCALL DDensGIG(double dX, double dNu, double dDelta, double dGamma);
double JDCALL DDensGH(double dX, double dNu, double dDelta, double dGamma, double dBeta);


/* mtcomplx.c */
double JDCALL c_abs(double xr, double xi);
BOOL   JDCALL c_div(double xr, double xi, double yr, double yi, double *zr, double *zi);
void   JDCALL c_mul(double xr, double xi, double yr, double yi, double *zr, double *zi);
void   JDCALL c_sqrt(double xr, double xi, double *yr, double *yi);
void   JDCALL c_exp(double xr, double xi, double *yr, double *yi);
void   JDCALL c_log(double xr, double xi, double *yr, double *yi);
void   JDCALL c_pow(double xr, double xi, double a, double *yr, double *yi);
BOOL   JDCALL FftDiscrete(VECTOR vXr, VECTOR vXi, int cN, int iForward);
void   JDCALL c_erf(double x, double y, double *erfx, double *erfy);
void   JDCALL c_sin(double xr, double xi, double *yr, double *yi);
void   JDCALL c_cos(double xr, double xi, double *yr, double *yi);
void   JDCALL c_tan(double xr, double xi, double *yr, double *yi);
void   JDCALL c_cot(double xr, double xi, double *yr, double *yi);
void   JDCALL c_sinh(double xr, double xi, double *yr, double *yi);
void   JDCALL c_cosh(double xr, double xi, double *yr, double *yi);
void   JDCALL c_tanh(double xr, double xi, double *yr, double *yi);
void   JDCALL c_asin(double xr, double xi, double *yr, double *yi);
void   JDCALL c_acos(double xr, double xi, double *yr, double *yi);
void   JDCALL c_atan(double xr, double xi, double *yr, double *yi);

/* mtchol.h */
int    JDCALL ILDLdec(MATRIX, VECTOR, int);
BOOL   JDCALL FPPtDec(MATRIX, int);
void   JDCALL LDLsolve(MATRIX, VECTOR, VECTOR, VECTOR, int );
void   JDCALL LDLsolveInv(MATRIX mLDLt, MATRIX mAinv, int cA);
int    JDCALL ISymInv(MATRIX, int);
int    JDCALL ISymInvDet(MATRIX mA, int cA, double *pdLogDet);
double JDCALL DGetInvertEpsNorm(MATRIX mA, int cA);
int    JDCALL ILUPlogdet(MATRIX mU, int cA, int *piPiv, double dNormEps, double *pdLogDet);
int    JDCALL ILUPdec(MATRIX mA, int cA, int *piPiv, double *pdLogDet,
              int *piSignDet, MATRIX mUt);
int    JDCALL ILUPdecNoPiv(MATRIX mA, int cA);
void   JDCALL LUPsolve(MATRIX mL, MATRIX mU, int *piPiv, VECTOR vB, int cA);
void   JDCALL LUPsolveEx(MATRIX mL, MATRIX mU, int *piPiv, MATRIX mB, int cA, int cB, BOOL fUseDiagL, BOOL fBisI);
void   JDCALL LUPsolveInv(MATRIX mL, MATRIX mU, int *piPiv, MATRIX mAinv, int cA);
int    JDCALL IInvert(MATRIX, int);
int    JDCALL IInvDet(MATRIX mA, int cA, double *pdLogDet, int *piSignDet);

/* mtdenest.c */
int    JDCALL IDenest(double *, int, double, double, double, double, double, double *,
      double *, int, int);

/* mteigsym.c */
int    JDCALL IEigenSym(MATRIX mA, int cA, VECTOR vEval, int fDoVectors);
int    JDCALL IEigValSym(MATRIX, VECTOR, int);
int    JDCALL IEigVecSym(MATRIX, VECTOR, int);
int    JDCALL IGenEigVecSym(MATRIX, MATRIX, VECTOR, VECTOR, int);
void   JDCALL HHsym(MATRIX, VECTOR, VECTOR, int, double, BOOL);
int    JDCALL IQLsym(MATRIX, VECTOR, VECTOR, int, int);

/* mteigval.c */
void   JDCALL Balance(MATRIX, VECTOR, int, int *, int *);
void   JDCALL DeBalance(MATRIX, VECTOR, int, int, int, int);
int    JDCALL IEigValPoly(VECTOR vPoly, VECTOR evr, VECTOR evi, int dim);
int    JDCALL IEigValReal(MATRIX, VECTOR, VECTOR, int);
void   JDCALL HHreal(MATRIX, VECTOR, int, int lo, int hi, double);
void   JDCALL QLreal(MATRIX, MATRIX, VECTOR, VECTOR, int, int, int, int, int *);

/* mteigvec.c */
int    JDCALL IEigen(MATRIX mA, int cA, VECTOR vRe, VECTOR vIm, MATRIX mEvec);
void   JDCALL EigVecDiv(MATRIX, VECTOR, VECTOR, int);
int    JDCALL IEigVecReal(MATRIX, VECTOR, VECTOR, int);
void   JDCALL HHrealBak(MATRIX, MATRIX, VECTOR, int);
void   JDCALL HHrealTrans(MATRIX, MATRIX, VECTOR, int, int, int);
void   JDCALL QLrealBak(MATRIX, MATRIX, VECTOR, VECTOR, int, int, int);

/* mtfft.c */
void   JDCALL FftComplex(VECTOR vXr, VECTOR vXi, int iPower, int iForward);
void   JDCALL FftReal(VECTOR vXr, VECTOR vXi, int iPower, int iForward);
int    JDCALL FFT1d(MATRIX mDest, MATRIX mSrc, int n, int iReverse, int isComplex);

/* mtfilter.c */
double JDCALL DEpKernelBandWidth(VECTOR vX, int iT1, int iT2);
BOOL   JDCALL FEpKernel(VECTOR vY, VECTOR vT, int cT, double *pdAlpha, VECTOR vG,
	VECTOR vX, double *pdCV, double *pdPar, BOOL fAuto, int iDesiredPar);
BOOL   JDCALL FCubicSpline(VECTOR vY, VECTOR vT, int cT, double *pdAlpha, VECTOR vG,
    VECTOR vX, double *pdCV, double *pdEDF, BOOL fAuto, int iDesiredDf);
BOOL   JDCALL FCubicSplineTime(VECTOR vY, int cT, double dAlpha, VECTOR vG, BOOL fHP);
BOOL   JDCALL FEWMA(VECTOR vY, int cT, double dLambda0, double dLambda1, VECTOR vG);
BOOL   JDCALL FEWMC(VECTOR vY, VECTOR vX, int cT, double dLambda, BOOL fMean0, VECTOR vG);

/* mtgamma.c */
double JDCALL DBetaFunc(double x, double a, double b);
double JDCALL DBetaFunc1(double x, double a, double b);
double JDCALL DGammaFunc(double x, double p);
double JDCALL DGammaFuncUpper(double x, double p);
double JDCALL DLogGamma(double z);
double JDCALL DPolyGamma(double z, int n);
double JDCALL DGamma(double z);
double JDCALL DFactorial(double z);
double JDCALL DBinomial(double n, double k);
void   JDCALL CLogGamma(double zr, double zi, double *yr, double *yi);
void   JDCALL CGamma(double zr, double zi, double *yr, double *yi);
void   JDCALL CPolyGamma(double zr, double zi, int n, double *yr, double *yi);
double JDCALL DPochhammer(double a, int n);
double JDCALL DPochhammerRatio_22(double a1, double a2, double b1, double b2, int n);
double JDCALL DGammaRatio_22(double a1, double a2, double b1, double b2);
double JDCALL DGammaRatio_23(double a1, double a2, double b1, double b2, double b3);

/* mthyper.c */
double JDCALL Hyper_2F1(double a, double b, double c, double z);

/* mtorth.c */
int    JDCALL IOrthMGS(MATRIX, int, int, int, MATRIX);
void   JDCALL DGivens(double a, double b, double *pdC, double *pdS);
MATRIX JDCALL GivensApplyRow(MATRIX mA, int i1, int i2, int c1, int c2, double c, double s);
MATRIX JDCALL GivensApplyCol(MATRIX mA, int i1, int i2, int r1, int r2, double c, double s);
void   JDCALL GivensSweepLZDiag(MATRIX mL, int rL, int cL, MATRIX mZ, int cZ, int jFrom, int jTo);
void   JDCALL GivensSweepLZRow(MATRIX mL, int rL, int cL, MATRIX mZ, int cZ, int iFrom, int iTo);
void   JDCALL GivensSweepQRDiag(MATRIX mR, int rR, int cR, MATRIX mQ, int rQ, int iFrom, int iTo);
void   JDCALL GivensSweepQRCol(MATRIX mR, int rR, int cR, MATRIX mQ, int rQ, int jFrom, int jTo);

/* mtprob.c */
double JDCALL DErf(double x);

double JDCALL DDensBeta(double x, double a, double b);
double JDCALL DDensBinomial(double dX, double dN, double p);
double JDCALL DDensCauchy(double x);
double JDCALL DDensChi(double chi, double dDf);
double JDCALL DDensExp(double x, double dLambda);
double JDCALL DDensExtremeValue(double x, double dAlpha, double dBeta);
double JDCALL DDensF(double f, double n1, double n2);
double JDCALL DDensGamma(double, double, double);
double JDCALL DDensGeometric(double dX, double dP);
double JDCALL DDensHypergeometric(double dX, double dN, double dK, double dM);
double JDCALL DDensInvGaussian(double x, double dMu, double dLambda);
double JDCALL DDensKernel(double x, int iType);
double JDCALL DDensLogNormal(double x);
double JDCALL DDensLogarithmic(double dX, double dA);
double JDCALL DDensLogistic(double x, double dAlpha, double dBeta);
double JDCALL DDensMises(double x, double dMu, double dKappa);
double JDCALL DDensNegBin(double dX, double dK, double dP);
double JDCALL DDensNormal(double x);
double JDCALL DDensPareto(double dX, double dK, double dA);
double JDCALL DDensPoisson(double dMu, int iK);
double JDCALL DDensPoissonEx(double dX, double dMu);
double JDCALL DDensT(double t, double n1);
double JDCALL DDensWeibull(double dX, double dA, double dB);

double JDCALL DProbBeta(double x, double a, double b);
double JDCALL DProbBinomial(double dX, double dN, double dP);
double JDCALL DProbCauchy(double x);
double JDCALL DProbChi(double, double);
double JDCALL DProbExp(double x, double dLambda);
double JDCALL DProbExtremeValue(double x, double dAlpha, double dBeta);
double JDCALL DProbF(double, double, double);
double JDCALL DProbGamma(double, double, double);
double JDCALL DProbGeometric(double dX, double dP);
double JDCALL DProbHypergeometric(double dX, double dN, double dK, double dM);
double JDCALL DProbInvGaussian(double x, double dMu, double dLambda);
double JDCALL DProbLogNormal(double x);
double JDCALL DProbLogarithmic(double dX, double dA);
double JDCALL DProbLogistic(double x, double dAlpha, double dBeta);
double JDCALL DProbMises(double x, double dMu, double dKappa);
double JDCALL DProbNegBin(double dX, double dK, double dP);
double JDCALL DProbNormal(double);
double JDCALL DProbPareto(double dX, double dK, double dA);
double JDCALL DProbPoisson(double dMu, int iK);
double JDCALL DProbPoissonEx(double dX, double dMu);
double JDCALL DProbT(double, int);
double JDCALL DProbWeibull(double dX, double dA, double dB);

double JDCALL DTailProbChi(double, double);
double JDCALL DTailProbF(double, double, double);
double JDCALL DTailProbNormal(double);
double JDCALL DTailProbT(double, int);

double JDCALL DProbChiNc(double, double, double dnc);
double JDCALL DProbFNc(double, double, double, double dnc);
double JDCALL DProbTNc(double, double, double dnc);

double JDCALL DProbBVN(double dLo1, double dLo2, double dRho);
double JDCALL DProbMVN(int n, VECTOR vX, MATRIX mSigma);

/* mtquan.c */
double JDCALL DFval(double, double, double, int *);
double JDCALL DQuanBeta(double x, double a, double b);
double JDCALL DQuanBinomial(double p, double dN, double dP);
double JDCALL DQuanCauchy(double p);
double JDCALL DQuanChi(double p, double v);
double JDCALL DQuanExp(double x, double dLambda);
double JDCALL DQuanExtremeValue(double p, double dAlpha, double dBeta);
double JDCALL DQuanF(double, double, double);
double JDCALL DQuanGamma(double p, double r, double a);
double JDCALL DQuanGeometric(double p, double dP);
double JDCALL DQuanHypergeometric(double p, double dN, double dK, double dM);
double JDCALL DQuanInvGaussian(double x, double dMu, double dLambda);
double JDCALL DQuanLogNormal(double x);
double JDCALL DQuanLogarithmic(double p, double dA);
double JDCALL DQuanLogistic(double p, double dAlpha, double dBeta);
double JDCALL DQuanMises(double p, double dMu, double dKappa);
double JDCALL DQuanNegBin(double p, double dK, double dP);
double JDCALL DQuanNormal(double dProb);
double JDCALL DQuanPareto(double p, double dK, double dA);
double JDCALL DQuanPoisson(double p, double dMu);
double JDCALL DQuanT(double p, int v);
double JDCALL DQuanTD(double p, double v);
double JDCALL DQuanWeibull(double p, double dA, double dB);

/* mtrandom.c */
void   JDCALL SetRanLocker(struct ranState * (JDCALL * pfLocker)(int lockit));
void   JDCALL RanInit(struct ranState *pRan);
void   JDCALL RanCopy(struct ranState *pRan);
unsigned int JDCALL IRanU(void);
double JDCALL DRanU(void);
void   JDCALL RanSetSeed(int *piSeed, int cSeed);
int    JDCALL RanGetSeed(int *piSeed, int cSeed);
void   JDCALL RanSetLoopSeed(int iLoop, int iStage);
void   JDCALL RanSetRan(const char *sRan);
const char *JDCALL RanGetRan(void);
void   JDCALL RanNewRNG(UIRANFUN fnUIRanu, RANSETSEEDFUN fnRanSetSeed,
	RANGETSEEDFUN fnRanGetSeed, RANCLRSEEDFUN fnRanClrSeed);
MATRIX JDCALL MatRan(MATRIX, int, int);
MATRIX JDCALL MatRann(MATRIX mA, int cR, int cC);
double JDCALL DRanBeta(double a, double b);
int    JDCALL IRanBinomial(int, double);
void   JDCALL RanBrownianMotion(VECTOR vY, int cY, VECTOR vEvalTimes);
double JDCALL DRanCauchy(void);
double JDCALL DRanChi(double);
void   JDCALL RanDirichlet(VECTOR vX, VECTOR vAlpha, int cAlpha);
double JDCALL DRanExp(double);
double JDCALL DRanExtremeValue(double dAlpha, double dBeta);
double JDCALL DRanF(double, double);
double JDCALL DRanGIG(double dNu, double dDelta, double dGamma);
double JDCALL DRanGamma(double, double);
int    JDCALL IRanGeometric(double dP);
double JDCALL DRanHypergeometric(double dN, double dK, double dM);
double JDCALL DRanInvGaussian(double dMu, double dLambda);
double JDCALL DRanLogNormal(void);
int    JDCALL IRanLogarithmic(double dA);
double JDCALL DRanLogistic(void);
double JDCALL DRanMises(double dKappa);
int    JDCALL IRanNegBin(int iN, double dP);
double JDCALL DRanNormal(void);
double JDCALL DRanNormalPM(void);
double JDCALL DRanPareto(double dK, double dA);
int    JDCALL IRanPoisson(double);
void   JDCALL RanPoissonProcess(VECTOR vY, int cY, VECTOR vEvalTimes, double dMu);
double JDCALL DRanStable(double dA, double dB);
double JDCALL DRanStudentT(double dDf);
void   JDCALL RanSubSample(VECTOR vU, int cU, int cN);
void   JDCALL RanSubSampleVec(VECTOR vU, int cU, VECTOR vX, int cX);
double JDCALL DRanT(int);
void   JDCALL RanUorder(VECTOR vU, int cU);
double JDCALL DRanWeibull(double dA, double dB);
void   JDCALL RanWishart(MATRIX mX, int cX, int cT);

/* mtstats.c */
BOOL   JDCALL FGetAcf(VECTOR vX, int cT, int cLag, VECTOR vAcf);
int    JDCALL IGetAcf(VECTOR, int, int, VECTOR, BOOL bCov);
int    JDCALL IGetAcfCross(VECTOR vX, VECTOR vY, int cT, int cLag, VECTOR vAcf, double *pdAcf0);
BOOL   JDCALL FGetAcfRun(VECTOR, int, int, VECTOR);
int    JDCALL IGetAcfRun(VECTOR vX, int cT, int cLag, VECTOR vAcf, BOOL bCov);
MATRIX JDCALL MatAcf(MATRIX mAcf, MATRIX mX, int cT, int cX, int mxLag);
MATRIX JDCALL MatVariance(MATRIX mXtX, MATRIX mX, int cT, int cX, BOOL fCorr);
MATRIX JDCALL MatStandardize(MATRIX mXdest, MATRIX mX, int cT, int cX);
BOOL   JDCALL FPeriodogram(VECTOR vX, int cT, int iTrunc, int cS, VECTOR vS, int iMode);
BOOL   JDCALL FPeriodogramAcf(VECTOR vAcf, int cT, int iTrunc, int cS, VECTOR vS, int iMode, int cTwin);
void   JDCALL OlsQRacc(MATRIX mXt, int cX, int cT, int *piPiv, int cR, VECTOR vTau,
			MATRIX mYt, int cY,	MATRIX mB, MATRIX mXtXinv, MATRIX mXtX);
int    JDCALL IOlsQR(MATRIX mXt, int cX, int cT, MATRIX mYt, int cY,
            MATRIX mB, MATRIX mXtXinv, MATRIX mXtX, VECTOR vW);
int    JDCALL IOlsNorm(MATRIX mXt, int cX, int cT, MATRIX mYt, int cY, 
			MATRIX mB, MATRIX mXtXinv, MATRIX mXtX, BOOL fInRows);
int	   JDCALL IDecQRtRank(MATRIX mXt, int cX, int cT, int *pcR);
int    JDCALL IDecQRtEx(MATRIX mXt, int cX, int cT, int *piPiv, VECTOR vTau);
int    JDCALL IDecQRt(MATRIX mXt, int cX, int cT, int *piPiv, int *pcR);
void   JDCALL DecQRtMult(MATRIX mQt, int cX, int cT, MATRIX mYt, int cY, int cR);
void   JDCALL DecQRtMul(MATRIX mQt, int cX, int cT, MATRIX mY, int cY, int cR);
VECTOR JDCALL VecDiscretize(VECTOR vY, int cY, double dMin, double dMax,
			VECTOR vDisc, int cM, VECTOR vT, int iOption);
void   JDCALL Moments(MATRIX mX, int cT, int cX, BOOL fCheckNaN, BOOL fRatios,
			int cMom, MATRIX mMom);

/* mtsvd.c */
int    JDCALL IDecSVD(MATRIX mA, int cM, int cN, VECTOR vW, int fDoU, MATRIX mU,
    		int fDoV, MATRIX mV, int fSort);
int    JDCALL ISVDdec(MATRIX, int, int, VECTOR, BOOL, MATRIX, BOOL, MATRIX, VECTOR, BOOL);
int    JDCALL INullSpace(MATRIX, int, int, BOOL);
MATRIX JDCALL MatGenInvert(MATRIX mA, int cM, int cN, MATRIX mRes, VECTOR vSval);
MATRIX JDCALL MatGenInvertSym(MATRIX mA, int cA, MATRIX mRes, VECTOR vSval);
int    JDCALL ISVDdecVals(MATRIX mA, int cM, int cN, VECTOR vW);
int	   JDCALL IMatRank(MATRIX mA, int cM, int cN, double dEps, BOOL bAbsolute);
double JDCALL DMatNorm(MATRIX mA, int cM, int cN, int iType);


#ifdef __cplusplus
}
#endif

#endif  /* INC_JDMATH */

