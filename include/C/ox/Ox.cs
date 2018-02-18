using System;
using System.Runtime.InteropServices;

static class Ox
{
// basic mathematical and statistical functions and allocation/deallocation
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern double c_abs(double xr, double xi);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern Int32 c_div(double xr, double xi, double yr, double yi, ref double zr, ref double zi);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern void c_erf(double x, double y, ref double erfx, ref double erfy);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern void c_exp(double xr, double xi, ref double yr, ref double yi);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern void c_log(double xr, double xi, ref double yr, ref double yi);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern void c_mul(double xr, double xi, double yr, double yi, ref double zr, ref double zi);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern void c_sqrt(double xr, double xi, ref double yr, ref double yi);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern double DBessel01(double x, Int32 iType, Int32 n);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern double DBesselNu(double x, Int32 iType, double n);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern double DBetaFunc(double dX, double dA, double dB);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern double DDawson(double dX);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern double DDensBeta(double x, double a, double b);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern double DDensChi(double x, double dDf);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern double DDensF(double x, double dDf1, double dDf2);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern double DDensGamma(double g, double r, double a);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern double DDensGH(double dX, double dNu, double dDelta, double dGamma, double dBeta);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern double DDensGIG(double dX, double dNu, double dDelta, double dGamma);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern double DDensMises(double x, double dMu, double dKappa);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern double DDensNormal(double x);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern double DDensPoisson(double dMu, Int32 k);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern double DDensT(double x, double dDf);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern double DDiagXSXt(Int32 iT, IntPtr mX, IntPtr mS, Int32 cS);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern double DDiagXtSXtt(Int32 cX, IntPtr mXt, Int32 mS, Int32 cS);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern void DecQRtMul(IntPtr mQt, Int32 cX, Int32 cT, IntPtr mY, Int32 cY, Int32 cR);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern void DecQRtMult(IntPtr mQt, Int32 cX, Int32 cT, IntPtr mYt, Int32 cY, Int32 cR);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern double DErf(double z);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern double DExpInt(double z);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern double DExpInt1(double z);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern double DExpInte(double z);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern double DGamma(double z);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern double DGammaFunc(double dX, double dR);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern double DGetInvertEps();
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern double DGetInvertEpsNorm(IntPtr mA, double cA);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern double DLogGamma(double dA);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern double DPolyGamma(double dA, Int32 n);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern double DProbBeta(double x, double a, double b);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern double DProbBVN(double dLo1, double dLo2, double dRho);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern double DProbChi(double x, double dDf);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern double DProbChiNc(double x, double df, double dNc);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern double DProbF(double x, double dDf1, double dDf2);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern double DProbFNc(double x, double dDf1, double dDf2, double dNc);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern double DProbGamma(double x, double dR, double dA);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern double DProbMises(double x, double dMu, double dKappa);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern double DProbNormal(double x);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern double DProbPoisson(double dMu, Int32 k);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern double DProbT(double x, Int32 iDf);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern double DProbTNc(double x, double dDf, double dNc);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern double DQuanBeta(double x, double a, double b);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern double DQuanChi(double p, double dDf);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern double DQuanF(double p, double dDf1, double dDf2);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern double DQuanGamma(double p, double r, double a);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern double DQuanMises(double p, double dMu, double dKappa);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern double DQuanNormal(double p);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern double DQuanT(double p, Int32 iDf);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern double DQuanTD(double p, double dDf);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern double DRanBeta(double a, double b);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern double DRanChi(double dDf);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern double DRanExp(double dLambda);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern double DRanF(double dDf1, double dDf2);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern double DRanGamma(double dR, double dA);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern double DRanGIG(double dNu, double dDelta, double dGamma);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern double DRanInvGaussian(double dMu, double dLambda);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern double DRanLogistic();
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern double DRanLogNormal();
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern double DRanMises(double dKappa);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern double DRanNormalPM();
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern double DRanStable(double dA, double dB);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern double DRanStudentT(double dDf);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern double DRanT(Int32 iDf);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern double DRanU();
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern double DTailProbChi(double x, double dDf);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern double DTailProbF(double x, double dDf1, double dDf2);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern double DTailProbGamma(double x, double dR, double dA);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern double DTailProbNormal(double x);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern double DTailProbT(double x, Int32 iDf);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern double DTrace(IntPtr mat, Int32 cA);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern double DTraceAB(IntPtr mA, IntPtr mB, Int32 cM, Int32 cN);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern double DVecsum(double[] vA, Int32 cA);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern void EigVecDiv(ref IntPtr mmE, double[] vEr, double[] vEi, Int32 cA);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern Int32 FCubicSpline(double[] vY, double[] vT, Int32 cT, ref double pdAlpha, double[] vG, double[] vX, ref double pdCV, ref double pdPar, Int32 fAuto, Int32 iDesiredPar);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern Int32 FCubicSplineTime(double[] vY, double cT, double dAlpha, ref Int32 vG, Int32 fHP);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern Int32 FFT1d(IntPtr mDest, IntPtr mSrc, Int32 iDir, Int32 isComplex);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern void FftComplex(double[] vXr, double[] vXi, Int32 iPower, Int32 iDir);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern void FftReal(double[] vXr, double[] vXi, Int32 iPower, Int32 iDir);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern Int32 FftDiscrete(double[] vXr, double[] vXi, Int32 cN, Int32 iDir);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern Int32 FGetAcf(double[] vX, Int32 cT, Int32 cLag, double[] vAcf);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern Int32 FGetAcfRun(double[] vX, Int32 cT, Int32 cLag, double[] vAcf);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern Int32 FIsInf(double d);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern Int32 FIsNaN(double d);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern Int32 FPeriodogram(double[] vX, Int32 cT, Int32 iTrunc, Int32 cS, double[] vS, Int32 iMode);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern Int32 FPeriodogramAcf(double[] vAcf, Int32 cT, Int32 iTrunc, Int32 cS1, double[] vS, Int32 iMode, Int32 cTwin);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern object FPPtDec(IntPtr mA, Int32 cA);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern Int32 IDecQRt(IntPtr mXt, Int32 cX, Int32 cT, ref Int32 piPiv, ref Int32 pcR);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern Int32 IDecQRtEx(IntPtr mXt, Int32 cX, Int32 cT, ref Int32 piPiv, double[] vTau);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern Int32 IDecQRtRank(IntPtr mXt, Int32 cX, Int32 cT, ref Int32 pcR);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern Int32 IDecSVD(IntPtr mA, Int32 cM, Int32 cN, double[] vW, Int32 fDoU, IntPtr mU, Int32 fDoV, IntPtr mV, Int32 fSort);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern Int32 IEigen(IntPtr mA, Int32 cA, double[] vEr, double[] vEi, IntPtr mmE);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern Int32 IEigenSym(IntPtr mA, Int32 cA, double[] vEval, Int32 fDoVectors);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern Int32 IEigValPoly(double[] vPoly, double[] vEr, double[] vEi, Int32 cA);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern Int32 IEigValReal(IntPtr mA, double[] vEr, double[] vEi, Int32 cA);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern Int32 IEigValSym(IntPtr mA, double[] vEv, Int32 cA);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern Int32 IEigVecReal(IntPtr mA, double[] vEr, double[] vEi, Int32 cA);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern Int32 IEigVecSym(IntPtr mA, double[] vEv, Int32 cA);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern Int32 IGenEigVecSym(IntPtr mA, IntPtr mB, double[] vEval, double[] vSubd, Int32 cA);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern Int32 IGetAcf(double[] vX, Int32 cT, Int32 cLag, double[] vAcf, Int32 bCov);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern Int32 IInvDet(IntPtr mA, Int32 cA, ref double pdLogDet, ref Int32 piSignDet);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern Int32 IInvert(IntPtr mA, Int32 cA);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern Int32 ILDLbandDec(IntPtr mA, double[] vD, Int32 cB, Int32 cA);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern Int32 ILDLdec(IntPtr mA, double[] vD, Int32 cA);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern Int32 ILUPdec(IntPtr mA, Int32 cA, ref Int32 piPiv, ref double pdLogDet, ref Int32 piSignDet, IntPtr mUt);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern Int32 ILUPlogdet(IntPtr mU, Int32 cA, ref Int32 piPiv, double dNormEps, ref double pdLogDet);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern Int32 IMatRank(IntPtr mA, Int32 cM, Int32 cN, double dEps, Int32 bAbsolute);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern IntPtr IntMatAllocBlock(Int32 cR, Int32 cC);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern void IntMatFreeBlock(IntPtr m);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern Int32 INullSpace(IntPtr mA, Int32 cM, Int32 cN, Int32 fAppend);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern Int32 IOlsNorm(IntPtr mXt, Int32 cX, Int32 cT, IntPtr mYt, Int32 cY, IntPtr mB, IntPtr mXtXinv, IntPtr mXtX, Int32 fInRows);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern Int32 IOlsQR(IntPtr mXt, Int32 cX, Int32 cT, IntPtr mYt, Int32 cY, IntPtr mB, IntPtr mXtXinv, IntPtr mXtX, ref Int32 vW);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern Int32 IRanBinomial(Int32 n, double p);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern Int32 IRanLogarithmic(double dA);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern Int32 IRanNegBin(Int32 iiN, double dP);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern Int32 IRanPoisson(double dMu);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern Int32 ISVDdec(IntPtr mA, Int32 cM, Int32 cN, double[] vW, Int32 fDoU, IntPtr mU, Int32 fDoV, IntPtr mV, double[] v_1, Int32 fSort);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern Int32 ISymInv(IntPtr mA, Int32 cA);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern Int32 ISymInvDet(IntPtr mA, Int32 cA, ref double pdLogDet);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern void LDLbandSolve(IntPtr mL, double[] vD, double[] vX, double[] vB, Int32 cB, Int32 cA);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern void LDLsolve(IntPtr mL, double[] vD, double[] vX, double[] vB, Int32 cA);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern void LDLsolveInv(IntPtr mLDLt, IntPtr mAinv, Int32 cA);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern void LUPsolve(IntPtr mL, IntPtr mU, ref Int32 piPiv, double[] vB, Int32 cA);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern void LUPsolveInv(IntPtr mL, IntPtr mU, ref Int32 piPiv, IntPtr mAinv, Int32 cA);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern Int32 MatAB(IntPtr mA, Int32 cA, Int32 cC, IntPtr mB, Int32 cB, IntPtr mAB);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern Int32 MatABt(IntPtr mA, Int32 cA, Int32 cC, IntPtr mB, Int32 cB, IntPtr mABt);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern Int32 MatAcf(IntPtr mAcf, IntPtr mX, Int32 cT, Int32 cX, Int32 mxLag);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern Int32 MatAdd(IntPtr mA, Int32 cM, Int32 cN, IntPtr mB, double dFac, IntPtr mAplusB);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern IntPtr MatAllocBlock(Int32 cR, Int32 cC);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern Int32 MatAtB(IntPtr mA, Int32 cA, Int32 cC, IntPtr mB, Int32 cB, IntPtr mAtB);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern Int32 MatBBt(IntPtr mB, Int32 cB, Int32 cS, IntPtr mBBt);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern Int32 MatBSBt(IntPtr mB, Int32 cB, IntPtr mS, Int32 cS, IntPtr mBSBt);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern Int32 MatBtB(IntPtr mB, Int32 cB, Int32 cS, Int32 mBtB);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern Int32 MatBtBVec(IntPtr mB, Int32 cB, Int32 cS, double[] vY, IntPtr mBtB);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern Int32 MatBtSB(IntPtr mB, Int32 cB, IntPtr mS, Int32 cS, IntPtr mBtSB);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern Int32 MatCopy(IntPtr mDest, IntPtr mSrc, Int32 cM, Int32 cN);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern Int32 MatCopyTranspose(IntPtr mDestT, IntPtr mSrc, Int32 cM, Int32 cN);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern void MatCopyVecc(IntPtr mDest, double[] vSrc_c, Int32 cM, Int32 cN);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern void MatCopyVecr(IntPtr mDest, double[] vSrc_r, Int32 cM, Int32 cN);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern Int32 MatDup(IntPtr mSrc, Int32 cM, Int32 cN);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern void MatFreeBlock(IntPtr m);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern Int32 MatGenInvert(IntPtr mA, Int32 cM, Int32 cN, IntPtr mRes, double[] vSval);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern double MatGetAt(IntPtr mSrc, Int32 i, Int32 j);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern Int32 MatI(IntPtr mDest, Int32 cM);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern Int32 MatNaN(IntPtr mDest, Int32 cM, Int32 cN);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern Int32 MatPartAcf(IntPtr mPartAcf, IntPtr mAcf, Int32 cAcf, IntPtr mY, Int32 cY, ref double pdLogDet, Int32 bFilter);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern Int32 MatRan(IntPtr mA, Int32 cR, Int32 cC);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern Int32 MatRann(IntPtr mA, Int32 cR, Int32 cC);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern Int32 MatReflect(IntPtr mA, Int32 cA);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern void MatSetAt(IntPtr mDest, double d, Int32 i, Int32 j);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern Int32 MatStandardize(IntPtr mXdest, IntPtr mX, Int32 cT, Int32 cX);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern Int32 MatTranspose(IntPtr mA, Int32 cA);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern Int32 MatVariance(IntPtr mXtX, IntPtr mX, Int32 cT, Int32 cX, Int32 fCorr);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern Int32 MatZero(IntPtr mDest, Int32 cM, Int32 cN);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern void OlsQRacc(IntPtr mXt, Int32 cX, Int32 cT, Int32[] piPiv, Int32 cR, double[] vTau, IntPtr mYt, Int32 cY, IntPtr mB, IntPtr mXtXinv,
	IntPtr mXtX);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern void RanDirichlet(double[] vX, double[] vAlpha, Int32 cALpha);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern Int32 RanGetSeed(Int32[] piSeed, Int32 cSeed);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern void RanInit();
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern void RanNewRan(IntPtr fnDRanu, IntPtr fnRanSetSeed, IntPtr fnRanGetSeed);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern void RanSetRan(string sRan);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern void RanSetSeed(Int32[] piSeed, Int32 cSeed);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern void RanSubSample(double[] vU, Int32 cU, Int32 cN);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern void RanUorder(double[] vU, Int32 cU);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern void RanWishart(IntPtr mX, Int32 cX);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern void SetFastMath(Int32 fYes);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern void SetInf(ref double pd);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern void SetInvertEps(double dEps);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern void SetNaN(ref double pd);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern void ToeplitzSolve(double[] vR, Int32 cR, Int32 cM, IntPtr mB, Int32 cB, double[] v_1);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern void VeccCopyMat(double[] vDest_c, IntPtr mSrc, Int32 cM, Int32 cN);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern Int32 VecDiscretize(double[] vY, Int32 cY, double dMin, double dMax, double[] vDisc, Int32 cM, double[] vT, Int32 iOption);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern Int32 VecDup(double[] vSrc, Int32 cM);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern void VecrCopyMat(double[] vDest_r, IntPtr mSrc, Int32 cM, Int32 cN);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern Int32 VecTranspose(double[] vA, Int32 cM, Int32 cN);

// Ox run-time functionality
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern Int32 FOxCallBack(IntPtr pvFunc, IntPtr rtn, IntPtr pv, Int32 cArg);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern Int32 FOxCallBackMember(IntPtr pvClass, string sMember, IntPtr rtn, IntPtr pv, Int32 cArg);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern Int32 FOxCreateObject(string sClass, IntPtr rtn, IntPtr pv, Int32 cArg);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern Int32 FOxGetDataMember(IntPtr pvClass, string sMember, IntPtr rtn);

	public delegate void NewOxFunc(IntPtr rtn, IntPtr pv, int cArg);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern Int32 FOxLibAddFunction(string sFunc, NewOxFunc pFunc, Int32 fVarArg);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern Int32 FOxLibAddFunctionEx(string sFunc, NewOxFunc pFunc, Int32 cArgs, Int32 flFlags);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern Int32 FOxRun(Int32 iMainIP, string sFunc);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern Int32 FOxSetDataMember(IntPtr pvClass, string sMember, IntPtr pv);

	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern Int32 IOxRunInit();
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern Int32 IOxVersion();
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern Int32 IOxVersionIsProfessional();
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern Int32 IOxVersionOxo();

	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern void OxCloneObject(IntPtr rtn, IntPtr pvObject, Int32 bDeep);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern void OxDeleteObject(IntPtr pvClass);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern void OxFnDouble(IntPtr rtn, IntPtr pv, IntPtr fn);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern void OxFnDouble2(IntPtr rtn, IntPtr pv, IntPtr fn);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern void OxFnDouble3(IntPtr rtn, IntPtr pv, IntPtr fn);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern void OxFnDouble4(IntPtr rtn, IntPtr pv, IntPtr fn);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern void OxFnDoubleInt(IntPtr rtn, IntPtr pv, IntPtr fn);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern void OxFreeByValue(IntPtr pv);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern void OxGetMainArgs(Int32[] pcArgc, ref string[] pasArgv);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern void OxGetOxArgs(Int32[] pcArgc, ref string[] pasArgv);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern Int32 OxGetPrintlevel();
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern Int32 OxGetUserExitCode();

	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern void OxLibArgError(Int32 iArg);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern void OxLibArgTypeError(Int32 iArg, Int32 iExpected, Int32 iFound);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern void OxLibCheckArrayMatrix(IntPtr pv, Int32 iFirst, Int32 iLast, IntPtr m);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern void OxLibCheckMatrixSize(IntPtr pv, Int32 iFirst, Int32 iLast, Int32 r, Int32 c);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern void OxLibCheckSquareMatrix(IntPtr pv, Int32 iFirst, Int32 iLast);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern void OxLibCheckType(Int32 iType, IntPtr pv, Int32 iFirst, Int32 iLast);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern void OxLibValArrayCalloc(IntPtr pv, Int32 c);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern void OxLibValMatDup(IntPtr pv, IntPtr mSrc, Int32 r, Int32 c);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern void OxLibValMatMalloc(IntPtr pv, Int32 r, Int32 c);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern void OxLibValStrMalloc(IntPtr pv, Int32 c);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern void OxLibValZero(IntPtr pv);

	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern Int32 OxMain(Int32 argc, string[] argv);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern Int32 OxMainCmd(string sCommand);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern void OxMainExit();
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern void OxMainInit();
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern void OxMakeByValue(IntPtr pv);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern void OxMessage(string s);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern void OxPuts(string s);

	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern void OxRunAbort(Int32 i);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern void OxRunError(Int32 iErno, string sToken);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern void OxRunErrorMessage(string s);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern void OxRunExit();
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern void OxRunMainExitCall(IntPtr fn);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern void OxRunMessage(string s);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern void OxRunWarningMessage(string sFunc, string sMsg);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern void OxSetMainArgs(Int32 argc, string[] argv);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern void OxSetOxArgs(Int32 cArgc, string[] asArgv);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern void OxSetPrintlevel(Int32 iSet);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern IntPtr OxStoreCreate(Int32 c);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern void OxStoreDelete(IntPtr pv, Int32 c);

	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern Int32 OxValColumns(IntPtr pv);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern Int32 OxValDuplicate2(IntPtr pvDest, IntPtr pvSrc);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern IntPtr OxValGetArray(IntPtr pv);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern Int32 OxValGetArrayLen(IntPtr pv);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern IntPtr OxValGetArrayVal(IntPtr pv, Int32 i);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern IntPtr OxValGetBlob(IntPtr pv, ref Int32 pI1, ref Int32 pI2);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern string OxValGetClassName(IntPtr pv);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern Int32 OxValGetDouble(IntPtr pv, ref double pdVal);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern Int32 OxValGetInt(IntPtr pv, ref Int32 piVal);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern IntPtr OxValGetMat(IntPtr pv);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern Int32 OxValGetMatc(IntPtr pv);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern Int32 OxValGetMatr(IntPtr pv);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern Int32 OxValGetMatrc(IntPtr pv);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern IntPtr OxValGetStaticObject(IntPtr pv);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern string OxValGetString(IntPtr pv);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern Int32 OxValGetStringCopy(IntPtr pv, string s, Int32 mxlen);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern Int32 OxValGetStringLen(IntPtr pv);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern IntPtr OxValGetVal(IntPtr pv, Int32 i);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern Int32 OxValGetVecc(IntPtr pv, double[] vX);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern Int32 OxValGetVecr(IntPtr pv, double[] vX);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern Int32 OxValHasFlag(IntPtr pv, Int32 iFlag);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern Int32 OxValHasType(IntPtr pv, Int32 iType);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern Int32 OxValRows(IntPtr pv);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern void OxValSetBlob(IntPtr pv, Int32 i1, Int32 i2, IntPtr p);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern void OxValSetDouble(IntPtr pv, double dVal);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern void OxValSetInt(IntPtr pv, Int32 iVal);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern void OxValSetMat(IntPtr pv, IntPtr mVal, int r, int c);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern void OxValSetMatZero(IntPtr pv, int r, int c);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern void OxValSetNull(IntPtr pv);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern void OxValSetString(IntPtr pv, string sVal);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern void OxValSetVecc(IntPtr pv, double[] vX, int r, int c);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern void OxValSetVecr(IntPtr pv, double[] vX, int r, int c);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern void OxValSetZero(IntPtr pv);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern Int32 OxValSizec(IntPtr pv);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern Int32 OxValSizer(IntPtr pv);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern Int32 OxValSizerc(IntPtr pv);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern Int32 OxValType(IntPtr pv);

	// these expect StdCall function pointers: Reverse P/Invoke
	public delegate void NewOxExit(int code);
	public delegate String NewOxGets(string str, int n);
	public delegate void NewOxMessage(string str);
	public delegate void NewOxPuts(string str);
	public delegate void NewOxRunMessage(string str);
	public delegate int  NewOxTextWindow(string sAction, string sArg);
	public delegate int  NewOxDraw(string sAction, int iArea, double[] vY, int cY, string sY, double[] vX, int cX, string sX, double[] vZ, int cZ, string sZ, int iSymbol, int iIndex, int[] piArgInt, int cArgInt, double[] pdArgDbl, int cArgDbl);
	public delegate int  NewOxDrawWindow(string, string);

	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern void SetOxExit(NewOxExit pfnNewOxExit);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern void SetOxGets(NewOxGets pfnNewOxGets);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern void SetOxMessage(NewOxMessage pfnNewOxMessage);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern void SetOxPipe(Int32 cPipe);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern void SetOxPuts(NewOxPuts pfnNewOxPuts);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern void SetOxRunMessage(NewOxRunMessage pfnNewOxRunMessage);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern void SetOxTextWindow(NewOxTextWindow pfnNewOxTextWindow);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern void SetOxDraw(NewOxDraw pfnNewOxDraw);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern void SetOxDrawWindow(NewOxDrawWindow pfnNewOxDrawWindow);

	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern string SOxGetTypeName(Int32 iType);
	[DllImport("OxWin", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern string SOxIntFunc();

	public enum OxTypes
	{
		OX_INT = 1,
		OX_DOUBLE = 2,
		OX_MATRIX = 3,
		OX_STRING = 4,
		OX_ARRAY = 5,
		OX_FUNCTION = 6,
		OX_CLASS = 7,
		OX_VECTOR = 8,
		OX_INTFUNC = 9,
		OX_RANGE = 10,
		OX_FILE = 11,
		OX_IMPORT = 12,
		OX_LAMBDA = 13,
		OX_BLOB = 14,
		OX_RETURN = 64,

		OX_NULL = 0x100,
		OX_VALUE = 0x200,
		OX_CONST = 0x400,
		OX_RESERVED = 0x800,
		OX_EXTERN = 0x1000,
		OX_GLOBAL = 0x2000,
		OX_STATDECL = 0x4000,
		OX_INLINE = 0x8000,
		OX_KEYWORD = 0x8000,
		OX_MEMBER = 0x10000,
		OX_STATIC = 0x20000,
		OX_VIRTUAL = 0x40000,
		OX_PUBLIC = 0x80000,
		OX_INDEX = 0x100000,
		OX_ADDRESS = 0x200000,
		OX_ARGUMENT = 0x400000,
		OX_VARARGS = 0x800000,
		OX_SERIAL = 0x1000000,
		OX_VECMAT = 0x2000000,
		OX_VECRANGE = 0x4000000,
		OX_IDXSCALAR = 0x8000000,
		OX_INTERNAL = 0x10000000
	}

	public enum OxErrors
	{
		ER_RUNTIME,
		ER_FATAL,
		ER_LIBFUNC,
		ER_STACKOVERFLOW,
		ER_OM,
		ER_NOFUNC,
		ER_DIVIDE,
		ER_EXPOBJECT,
		ER_UNDECLMEM,
		ER_SINGULAR,
		ER_BADOP,
		ER_INCTYPE,
		ER_NULL,
		ER_ARGS,
		ER_BADTYPE,
		ER_IDXBOUND,
		ER_IDXTYPE,
		ER_EXPSCALAR,
		ER_EXPSCALMAT,
		ER_EXPSCALSTR,
		ER_CONSTAS,
		ER_EXPCLASS,
		ER_UNDEFMEMFU,
		ER_NOMEMFUNC,
		ER_UNDEFFU,
		ER_INVALIDCLASS,
		ER_CONSTDEL,
		ER_ADDRESS,
		ER_NEGATIVE,
		ER_SVD,
		ER_CAST,
		ER_IDXRANGE,
		ER_IDXARRAY,
		ER_POWOVERFLOW,
		ER_ARGSAME,
		ER_NOCALLMEM,
		ER_EXPFUNC,
		ER_ARGUMENT,
		ER_IDXMISSING,
		ER_SIZENEGATIVE,
		ER_ARGERROR,
		ER_CALLBACK,
		ER_INVALID,
		ER_OVERFLOW,
		ER_FPEXCEPTION,
		ER_CASTINT,
		ER_CASTDBL,
		ER_CASTMAT,
		ER_CASTVEC,
		ER_USER1,
		ER_USER2,
		ER_USER3,
		ER_USER4,
		ER_USER5,
		ER_USER6,
		ER_EXPVECTOR,
		ER_JSTINVALID,
		ER_ARRAYMULTI,
		ER_NOPIPE,
		ER_USERABORT,
		ER_EXPPUBLIC,
		ER_EXSTRING,
		ER_CONSTUSE,
		ER_LAMBDA,
		ER_VARCHANGED,
		ER_STACK,
		ER_FOR,
		ER_FOR_INF,
		ER_LAST
	}

	public enum OxWarnings
	{
		WR_DECFAILED,
		WR_ITMAX,
		WR_CONCAT,
		WR_CASTINT,
		WR_VECIDXMAT,
		WR_DETERMINANT,
		WR_USER
	}
}
