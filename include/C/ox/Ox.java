package ox;

import com.sun.jna.Callback;
import com.sun.jna.Library;
import com.sun.jna.Native;
import com.sun.jna.Pointer;
import com.sun.jna.PointerType;
import com.sun.jna.ptr.DoubleByReference;
import com.sun.jna.ptr.IntByReference;
import com.sun.jna.ptr.PointerByReference;
import com.sun.jna.win32.StdCallLibrary;
import com.sun.jna.Platform;

public class Ox implements StdCallLibrary {
	static {
		Native.register(Platform.isWindows() ? "oxwin" : Platform.isMac() ? "libox.7.dylib" : "libox.so.7");
	}
	public static native double c_abs(double xr, double xi);
	public static native int c_div(double xr, double xi, double yr, double yi, DoubleByReference zr, DoubleByReference zi);
	public static native void c_erf(double x, double y, DoubleByReference erfx, DoubleByReference erfy);
	public static native void c_exp(double xr, double xi, DoubleByReference yr, DoubleByReference yi);
	public static native void c_log(double xr, double xi, DoubleByReference yr, DoubleByReference yi);
	public static native void c_mul(double xr, double xi, double yr, double yi, DoubleByReference zr, DoubleByReference zi);
	public static native void c_sqrt(double xr, double xi, DoubleByReference yr, DoubleByReference yi);
	public static native double DBessel01(double x, int type, int n);
	public static native double DBesselNu(double x, int type, double nu);
	public static native double DBetaFunc(double x, double a, double b);
	public static native double DDawson(double x);
	public static native double DDensBeta(double x, double a, double b);
	public static native double DDensChi(double chi, double dDf);
	public static native double DDensF(double f, double n1, double n2);
	public static native double DDensGamma(double double1, double double2, double double3);
	public static native double DDensGH(double dX, double dNu, double dDelta, double dGamma, double dBeta);
	public static native double DDensGIG(double dX, double dNu, double dDelta, double dGamma);
	public static native double DDensMises(double x, double dMu, double dKappa);
	public static native double DDensNormal(double x);
	public static native double DDensPoisson(double dMu, int iK);
	public static native double DDensT(double t, double n1);
	public static native double DDiagXSXt(int iT, Pointer mX, Pointer mS, int cX);
	public static native double DDiagXtSXtt(int iT, Pointer mXt, Pointer mS, int cX);
	public static native void DecQRtMul(Pointer mQt, int cX, int cT, Pointer mY, int cY, int cR);
	public static native void DecQRtMult(Pointer mQt, int cX, int cT, Pointer mYt, int cY, int cR);
	public static native double DEpKernelBandWidth(double[] vX, int iT1, int iT2);
	public static native double DErf(double x);
	public static native double DExpInt(double x);
	public static native double DExpInt1(double x);
	public static native double DExpInte(double x);
	public static native double DGamma(double z);
	public static native double DGammaFunc(double x, double p);
	public static native double DGetInvertEps();
	public static native double DGetInvertEpsNorm(Pointer mA, int cA);
	public static native double DLogGamma(double z);
	public static native double DPolyGamma(double z, int n);
	public static native double DProbBeta(double x, double a, double b);
	public static native double DProbBVN(double dLo1, double dLo2, double dRho);
	public static native double DProbChi(double double1, double double2);
	public static native double DProbChiNc(double double1, double double2, double dnc);
	public static native double DProbF(double double1, double double2, double double3);
	public static native double DProbFNc(double double1, double double2, double double3, double dnc);
	public static native double DProbGamma(double double1, double double2, double double3);
	public static native double DProbMises(double x, double dMu, double dKappa);
	public static native double DProbMVN(int n, double[] vX, Pointer mSigma);
	public static native double DProbNormal(double x);
	public static native double DProbPoisson(double dMu, int iK);
	public static native double DProbT(double x, int iNc);
	public static native double DProbTNc(double double1, double double2, double dnc);
	public static native double DQuanBeta(double x, double a, double b);
	public static native double DQuanChi(double p, double v);
	public static native double DQuanF(double double1, double double2, double double3);
	public static native double DQuanGamma(double p, double r, double a);
	public static native double DQuanMises(double p, double dMu, double dKappa);
	public static native double DQuanNormal(double dProb);
	public static native double DQuanT(double p, int v);
	public static native double DQuanTD(double p, double v);
	public static native double DRanBeta(double a, double b);
	public static native double DRanChi(double double1);
	public static native double DRanExp(double double1);
	public static native double DRanF(double double1, double double2);
	public static native double DRanGamma(double double1, double double2);
	public static native double DRanGIG(double dNu, double dDelta, double dGamma);
	public static native double DRanInvGaussian(double dMu, double dLambda);
	public static native double DRanLogistic();
	public static native double DRanLogNormal();
	public static native double DRanMises(double dKappa);
	public static native double DRanNormal();
	public static native double DRanNormalPM();
	public static native double DRanStable(double dA, double dB);
	public static native double DRanStudentT(double dDf);
	public static native double DRanT(int int1);
	public static native double DRanU();
	public static native double DTailProbChi(double double1, double double2);
	public static native double DTailProbF(double double1, double double2, double double3);
	public static native double DTailProbNormal(double double1);
	public static native double DTailProbT(double double1, int int1);
	public static native double DTrace(Pointer MATRIX1, int int1);
	public static native double DTraceAB(Pointer MATRIX1, Pointer MATRIX2, int int1, int int2);
	public static native double DVecsum(double[] vX, int cX);
	public static native void EigVecDiv(Pointer MATRIX1, double[] VECTOR1, double[] VECTOR2, int int1);
	public static native int FCubicSpline(double[] vY, double[] vT, int cT, double[] pdAlpha, double[] vG, double[] vX, DoubleByReference pdCV, DoubleByReference pdEDF, int fAuto, int iDesiredDf);
	public static native int FCubicSplineTime(double[] vY, int cT, double dAlpha, double[] vG, int fHP);
	public static native int FEpKernel(double[] vY, double[] vT, int cT, DoubleByReference pdAlpha, double[] vG, double[] vX, double[] pdCV, DoubleByReference pdPar, int fAuto, int iDesiredPar);
	public static native int FEWMA(double[] vY, int cT, double dLambda0, double dLambda1, double[] vG);
	public static native int FEWMC(double[] vY, double[] vX, int cT, double dLambda, int fMean0, double[] vG);
	public static native int FFT1d(Pointer mDest, Pointer mSrc, int n, int iReverse, int isComplex);
	public static native void FftComplex(double[] vXr, double[] vXi, int iPower, int iForward);
	public static native int FftDiscrete(double[] vXr, double[] vXi, int cN, int iForward);
	public static native void FftReal(double[] vXr, double[] vXi, int iPower, int iForward);
	public static native int FGetAcf(double[] vX, int cT, int cLag, double[] vAcf);
	public static native int FGetAcfRun(double[] VECTOR1, int int1, int int2, double[] VECTOR2);
	public static native int FGetPartAcf(double[] vAcf, int cAcf, double[] vPartAcf);
	public static native int FIsInf(double d);
	public static native int FIsInfOrNaN(double d);
	public static native int FIsNaN(double d);
	public static native int FOxCallBack(Pointer pvFunc, Pointer rtn, Pointer pv, int cArg);
	public static native int FOxCallBackMember(Pointer pvClass, String sMember, Pointer rtn, Pointer pv, int cArg);
	public static native int FOxCreateObject(String sClass, Pointer rtn, Pointer pv, int cArg);
	public static native int FOxGetDataMember(Pointer pvClass, String sMember, Pointer rtn);

	public interface CallbackOxFunc extends StdCallCallback {
	    void NewOxFunc(Pointer rtn, Pointer pv, int cArg);
	}
	public static native int FOxLibAddFunction(String sFunc, Ox.CallbackOxFunc pFunc, int fVarArg);
	public static native int FOxLibAddFunctionEx(String sFunc, Ox.CallbackOxFunc pFunc, int cArgs, int flFlags);

	public static native int FOxRun(int iMainIP, String sFunc);
	public static native int FOxSetDataMember(Pointer pvClass, String sMember, Pointer pv);
	public static native int FPeriodogram(double[] vX, int cT, int iTrunc, int cS, double[] vS, int iMode);
	public static native int FPeriodogramAcf(double[] vAcf, int cT, int iTrunc, int cS, double[] vS, int iMode, int cTwin);
	public static native int FPPtDec(Pointer MATRIX1, int int1);
	public static native double Hyper_2F1(double a, double b, double c, double z);
	public static native int IDecQRt(Pointer mXt, int cX, int cT, int[] piPiv, IntByReference pcR);
	public static native int IDecQRtEx(Pointer mXt, int cX, int cT, int[] piPiv, double[] vTau);
	public static native int IDecQRtRank(Pointer mXt, int cX, int cT, IntByReference pcR);
	public static native int IDecSVD(Pointer mA, int cM, int cN, double[] vW, int fDoU, Pointer mU, int fDoV, Pointer mV, int fSort);
	public static native int IDenest(double[] doublePtr1, int int1, double double1, double double2, double double3, double double4, double double5, double[] doublePtr2, double[] doublePtr3, int int2, int int3);
	public static native int IEigen(Pointer mA, int cA, double[] vRe, double[] vIm, Pointer mEvec);
	public static native int IEigenSym(Pointer mA, int cA, double[] vEval, int fDoVectors);
	public static native int IEigValPoly(double[] vPoly, double[] evr, double[] evi, int dim);
	public static native int IEigValReal(Pointer MATRIX1, double[] VECTOR1, double[] VECTOR2, int int1);
	public static native int IEigValSym(Pointer MATRIX1, double[] VECTOR1, int int1);
	public static native int IEigVecReal(Pointer MATRIX1, double[] VECTOR1, double[] VECTOR2, int int1);
	public static native int IEigVecSym(Pointer MATRIX1, double[] VECTOR1, int int1);
	public static native int IGenEigVecSym(Pointer MATRIX1, Pointer MATRIX2, double[] VECTOR1, double[] VECTOR2, int int1);
	public static native int IGetAcf(double[] VECTOR1, int int1, int int2, double[] VECTOR2, int bCov);
	public static native int IGetAcfCross(double[] vX, double[] vY, int cT, int cLag, double[] vAcf, double[] pdAcf0);
	public static native int IGetAcfRun(double[] vX, int cT, int cLag, double[] vAcf, int bCov);
	public static native int IInvDet(Pointer mA, int cA, DoubleByReference pdLogDet, IntByReference piSignDet);
	public static native int IInvert(Pointer MATRIX1, int int1);
	public static native int ILDLbandDec(Pointer MATRIX1, double[] VECTOR1, int int1, int int2);
	public static native int ILDLdec(Pointer MATRIX1, double[] VECTOR1, int int1);
	public static native int ILUPdec(Pointer mA, int cA, int[] piPiv, DoubleByReference pdLogDet, IntByReference piSignDet, Pointer mUt);
	public static native int ILUPlogdet(Pointer mU, int cA, int[] piPiv, double dNormEps, DoubleByReference pdLogDet);
	public static native int IMatRank(Pointer mA, int cM, int cN, double dEps, int bAbsolute);
	public static native Pointer IntMatAllocBlock(int cR, int cC);
	public static native void IntMatFreeBlock(Pointer mX);
	public static native int INullSpace(Pointer MATRIX1, int int1, int int2, int int3);
	public static native int IOlsNorm(Pointer mXt, int cX, int cT, Pointer mYt, int cY, Pointer mB, Pointer mXtXinv, Pointer mXtX, int fInRows);
	public static native int IOlsQR(Pointer mXt, int cX, int cT, Pointer mYt, int cY, Pointer mB, Pointer mXtXinv, Pointer mXtX, double[] vW);
	public static native int IOrthMGS(Pointer MATRIX1, int int1, int int2, int int3, Pointer MATRIX2);
	public static native int IOxRunInit();
	public static native int IOxVersion();
	public static native int IOxVersionIsProfessional();
	public static native int IOxVersionOxo();
	public static native int IRanBinomial(int int1, double double1);
	public static native int IRanLogarithmic(double dA);
	public static native int IRanNegBin(int iN, double dP);
	public static native int IRanPoisson(double double1);
	public static native int IRanU();
	public static native int IsDbMissing(double x);
	public static native int ISVDdec(Pointer MATRIX1, int int1, int int2, double[] VECTOR1, int int3, Pointer MATRIX2, int int4, Pointer MATRIX3, double[] VECTOR2, int int5);
	public static native int ISymInv(Pointer MATRIX1, int int1);
	public static native int ISymInvDet(Pointer mA, int cA, DoubleByReference pdLogDet);
	public static native void LDLbandSolve(Pointer MATRIX1, double[] VECTOR1, double[] VECTOR2, double[] VECTOR3, int int1, int int2);
	public static native void LDLsolve(Pointer MATRIX1, double[] VECTOR1, double[] VECTOR2, double[] VECTOR3, int int1);
	public static native void LDLsolveInv(Pointer mLDLt, Pointer mAinv, int cA);
	public static native void LUPsolve(Pointer mL, Pointer mU, int[] piPiv, double[] vB, int cA);
	public static native void LUPsolveInv(Pointer mL, Pointer mU, int[] piPiv, Pointer mAinv, int cA);
	public static native Pointer MatAB(Pointer MATRIX1, int int1, int int2, Pointer MATRIX2, int int3, Pointer MATRIX3);
	public static native Pointer MatABt(Pointer MATRIX1, int int1, int int2, Pointer MATRIX2, int int3, Pointer MATRIX3);
	public static native Pointer MatAcf(Pointer mAcf, Pointer mX, int cT, int cX, int mxLag);
	public static native Pointer MatAdd(Pointer mA, int cM, int cN, Pointer mB, double dFac, Pointer mAplusB);
	public static native Pointer MatAllocBlock(int cR, int cC);
	public static native Pointer MatAtB(Pointer MATRIX1, int int1, int int2, Pointer MATRIX2, int int3, Pointer MATRIX3);
	public static native Pointer MatBBt(Pointer mB, int rB, int cB, Pointer mBBt);
	public static native Pointer MatBSBt(Pointer MATRIX1, int int1, Pointer MATRIX2, int int2, Pointer MATRIX3);
	public static native Pointer MatBtB(Pointer mB, int rB, int cB, Pointer mBtB);
	public static native Pointer MatBtBVec(Pointer mB, int cT, int cB, double[] vY, Pointer mBtB);
	public static native Pointer MatBtSB(Pointer MATRIX1, int int1, Pointer MATRIX2, int int2, Pointer MATRIX3);
	public static native Pointer MatCopy(Pointer mDest, Pointer mSrc, int cR, int cC);
	public static native Pointer MatCopyTranspose(Pointer mB, Pointer mA, int cM, int cN);
	public static native void MatCopyVecc(Pointer mX, double[] vY, int cR, int cC);
	public static native void MatCopyVecr(Pointer mX, double[] vY, int cR, int cC);
	public static native Pointer MatCpy(Pointer mDest, Pointer mSrc, int cR, int cC);
	public static native Pointer MatDup(Pointer mSrc, int cR, int cC);
	public static native void MatFreeBlock(Pointer mX);
	public static native Pointer MatGenInvert(Pointer mA, int cM, int cN, Pointer mRes, double[] vSval);
	public static native Pointer MatGenInvertSym(Pointer mA, int cA, Pointer mRes, double[] vSval);
	public static native double MatGetAt(Pointer mX, int i, int j);
	public static native Pointer MatI(Pointer mX, int cX);
	public static native Pointer MatNaN(Pointer mX, int cR, int cC);
	public static native Pointer MatPartAcf(Pointer mPartAcf, Pointer mAcf, int cAcf, Pointer mY, int cY, DoubleByReference pdLogDet, int bFilter);
	public static native Pointer MatRan(Pointer MATRIX1, int int1, int int2);
	public static native Pointer MatRann(Pointer mA, int cR, int cC);
	public static native Pointer MatReflect(Pointer MATRIX1, int int1);
	public static native Pointer MatS1S2(Pointer MATRIX1, Pointer MATRIX2, int int1);
	public static native void MatSetAt(Pointer mX, double d, int i, int j);
	public static native Pointer MatStandardize(Pointer mXdest, Pointer mX, int cT, int cX);
	public static native Pointer MatSumOuter(Pointer mQ, Pointer mX, int cT, int cX);
	public static native Pointer MatTranspose(Pointer MATRIX1, int int1);
	public static native Pointer MatVariance(Pointer mXtX, Pointer mX, int cT, int cX, int fCorr);
	public static native Pointer MatXtXtt(Pointer mXt, int cX, int iT1, int iT2, Pointer mXtX);
	public static native Pointer MatXtZtt(Pointer mXt, int cX, int iT1, int iT2, Pointer mZt, int cZ, Pointer mXtZ);
	public static native Pointer MatXXt(Pointer mX, int iT1, int iT2, int cX, Pointer mXtX);
	public static native Pointer MatXXtVec(Pointer mX, int iT1, int iT2, int cX, double[] vY, Pointer mXtX);
	public static native Pointer MatZero(Pointer mX, int cR, int cC);
	public static native void OlsQRacc(Pointer mXt, int cX, int cT, int[] piPiv, int cR, double[] vTau, Pointer mYt, int cY, Pointer mB, Pointer mXtXinv, Pointer mXtX);
	public static native void OxCloneObject(Pointer rtn, Pointer pvObject, int bDeep);
	public static native void OxDeleteObject(Pointer pvClass);
	public static native void OxFreeByValue(Pointer pv);
	public static native void OxGetMainArgs(IntByReference pcArgc, PointerByReference pasArgv);
	public static native void OxGetOxArgs(IntByReference pcArgc, PointerByReference pasArgv);
	public static native int OxGetPrintlevel();
	public static native int OxGetUserExitCode();
	public static native void OxLibArgError(int iArg);
	public static native void OxLibArgTypeError(int iArg, int iExpected, int iFound);
	public static native void OxLibCheckArrayMatrix(Pointer pv, int iFirst, int iLast, Pointer m);
	public static native void OxLibCheckMatrixSize(Pointer pv, int iFirst, int iLast, int r, int c);
	public static native void OxLibCheckSquareMatrix(Pointer pv, int iFirst, int iLast);
	public static native void OxLibCheckType(int iType, Pointer pv, int iFirst, int iLast);
	public static native void OxLibValArrayCalloc(Pointer pv, int c);
	public static native void OxLibValMatDup(Pointer pv, Pointer mSrc, int r, int c);
	public static native void OxLibValMatMalloc(Pointer pv, int r, int c);
	public static native void OxLibValStrMalloc(Pointer pv, int c);
	public static native int OxMain(int argc, PointerByReference argv);
	public static native int OxMain_T(int argc, PointerByReference argv);
	public static native int OxMainCmd(String sCommand);
	public static native void OxMainExit();
	public static native void OxMainInit();
	public static native void OxMakeByValue(Pointer pv);
	public static native void OxMessage(String s);
	public static native void OxPuts(String s);
	public static native void OxRunAbort(int i);
	public static native void OxRunError(int iErno, String sToken);
	public static native void OxRunErrorMessage(String s);
	public static native void OxRunExit();
	public static native void OxRunMessage(String s);
	public static native void OxRunWarningMessage(String sFunc, String sMsg);
	public static native void OxSetMainArgs(int argc, PointerByReference argv);
	public static native void OxSetOxArgs(int cArgc, PointerByReference asArgv);
	public static native void OxSetPrintlevel(int iSet);
	public static native void OxSetUserExitCode(int iError);
	public static native Pointer OxStoreCreate(int c);
	public static native void OxStoreDelete(Pointer pv, int c);
	public static native int OxValColumns(Pointer pv);
	public static native void OxValDuplicate2(Pointer pvDest, Pointer pvSrc);
	public static native Pointer OxValGetArray(Pointer pv);
	public static native int OxValGetArrayLen(Pointer pv);
	public static native Pointer OxValGetArrayVal(Pointer pv, int i);
	public static native Pointer OxValGetBlob(Pointer pv, IntByReference pI1, IntByReference pI2);
	public static native Pointer OxValGetClassName(Pointer pv);
	public static native int OxValGetDouble(Pointer pv, DoubleByReference pdVal);
	public static native int OxValGetInt(Pointer pv, IntByReference piVal);
	public static native Pointer OxValGetMat(Pointer pv);
	public static native int OxValGetMatc(Pointer pv);
	public static native int OxValGetMatr(Pointer pv);
	public static native int OxValGetMatrc(Pointer pv);
	public static native Pointer OxValGetStaticObject(Pointer pv);
	public static native Pointer OxValGetString(Pointer pv);
	public static native int OxValGetStringCopy(Pointer pv, byte[] s, int mxLen);
	public static native int OxValGetStringLen(Pointer pv);
	public static native Pointer OxValGetVal(Pointer pv, int i);
	public static native int OxValGetVecc(Pointer pv, double[] vX);
	public static native int OxValGetVecr(Pointer pv, double[] vX);
	public static native int OxValHasFlag(Pointer pv, int iFlag);
	public static native int OxValHasType(Pointer pv, int iType);
	public static native int OxValRows(Pointer pv);
	public static native void OxValSetBlob(Pointer pv, int i1, int i2, Pointer p);
	public static native void OxValSetDouble(Pointer pv, double dVal);
	public static native void OxValSetInt(Pointer pv, int iVal);
	public static native void OxValSetMat(Pointer pv, Pointer mVal, int r, int c);
	public static native void OxValSetMatZero(Pointer pv, int r, int c);
	public static native void OxValSetNull(Pointer pv);
	public static native void OxValSetString(Pointer pv, String sVal);
	public static native void OxValSetZero(Pointer pv);
	public static native void OxValSetVecc(Pointer pv, double[] vX, int r, int c);
	public static native void OxValSetVecr(Pointer pv, double[] vX, int r, int c);
	public static native int OxValSizec(Pointer pv);
	public static native int OxValSizer(Pointer pv);
	public static native int OxValSizerc(Pointer pv);
	public static native void OxValTransfer(Pointer pv);
	public static native int OxValType(Pointer pv);
	public static native void QSortVecMat(double[] VECTOR1, double[] VECTOR2, Pointer MATRIX1, int int1, int int2, int int3);
	public static native void RanDirichlet(double[] vX, double[] vAlpha, int cAlpha);
	public static native int RanGetSeed(int[] piSeed, int cSeed);
	public static native void RanInit(Pointer pRan);
	public static native void RanSetRan(String sRan);
	public static native void RanSetSeed(int[] piSeed, int cSeed);
	public static native void RanSubSample(double[] vU, int cU, int cN);
	public static native void RanUorder(double[] vU, int cU);
	public static native void RanWishart(Pointer mX, int cX, int cT);
	public static native void SetInf(DoubleByReference pd);
	public static native void SetInvertEps(double dEps);
	public static native void SetNaN(DoubleByReference pd);
	public static native void SetOxPipe(int cPipe);
	public static native int SortMatCol(Pointer mX, int iCol, int cX);
	public static native int SortmXByCol(int iCol, Pointer mX, int cT, int cX);
	public static native int SortmXtByVec(int cT, double[] vBy, Pointer mXt, int cX);
	public static native int SortVec(double[] vX, int cX);
	public static native void SortVecByVec(int cT, double[] vBy, double[] vAc);
	public static native Pointer SOxGetTypeName(int iType);
	public static native Pointer SOxIntFunc();
	public static native int ToeplitzSolve(double[] VECTOR1, int int1, int int2, Pointer MATRIX1, int int3, double[] v_1);
	public static native DoubleByReference VecAllocBlock(int cX);
	public static native void VeccCopyMat(double[] vDest, Pointer mSrc, int cR, int cC);
	public static native DoubleByReference VecCopy(double[] vDest, double[] vSrc, int cX);
	public static native DoubleByReference VecDiscretize(double[] vY, int cY, double dMin, double dMax, double[] vDisc, int cM, double[] vT, int iOption);
	public static native DoubleByReference VecDupBlock(double[] vSrc, int cX);
	public static native void VecFreeBlock(double[] vX);
	public static native void VecrCopyMat(double[] vDest, Pointer mSrc, int cR, int cC);
	public static native DoubleByReference VecTranspose(double[] vA, int cM, int cN);

	public interface CallbackOxExit extends StdCallCallback {
	    void NewOxExit(int code);
	}
	public static native void SetOxExit(Ox.CallbackOxExit pfnNewOxExit);
	
	public static native void SetOxDraw(Pointer pfnNewOxDraw);
	public static native void SetOxDrawWindow(Pointer pfnNewOxDrawWindow);
	public static native void SetOxGets(Pointer pfnNewOxGets);
	public static native void SetOxMessage(Pointer pfnNewOxMessage);
	public static native void SetOxPuts(Pointer pfnNewOxPuts);
	public static native void SetOxRunMessage(Pointer pfnNewOxRunMessage);
	public static native void SetOxTextWindow(Pointer pfnNewOxTextWindow);
	
	public enum OxTypes
	{
		OX_INT(1),
		OX_DOUBLE(2),
		OX_MATRIX(3),
		OX_STRING(4),
		OX_ARRAY(5),
		OX_FUNCTION(6),
		OX_CLASS(7),
		OX_VECTOR(8),
		OX_INTFUNC(9),
		OX_RANGE(10),
		OX_FILE(11),
		OX_IMPORT(12),
		OX_LAMBDA(13),
		OX_BLOB(14),
		OX_RETURN(64),

		OX_NULL(0x100),
		OX_VALUE(0x200),
		OX_CONST(0x400),
		OX_RESERVED(0x800),
		OX_EXTERN(0x1000),
		OX_GLOBAL(0x2000),
		OX_STATDECL(0x4000),
		OX_INLINE(0x8000),
		OX_KEYWORD(0x8000),
		OX_MEMBER(0x10000),
		OX_STATIC(0x20000),
		OX_VIRTUAL(0x40000),
		OX_PUBLIC(0x80000),
		OX_INDEX(0x100000),
		OX_ADDRESS(0x200000),
		OX_ARGUMENT(0x400000),
		OX_VARARGS(0x800000),
		OX_SERIAL(0x1000000),
		OX_VECMAT(0x2000000),
		OX_VECRANGE(0x4000000),
		OX_IDXSCALAR(0x8000000),
		OX_INTERNAL(0x10000000);

		private int intValue;
		private static java.util.HashMap<Integer, OxTypes> mappings;
		private static java.util.HashMap<Integer, OxTypes> getMappings()
		{
			if (mappings == null)
			{
				synchronized (OxTypes.class)
				{
					if (mappings == null)
					{
						mappings = new java.util.HashMap<Integer, OxTypes>();
					}
				}
			}
			return mappings;
		}

		private OxTypes(int value)
		{
			intValue = value;
			OxTypes.getMappings().put(value, this);
		}

		public int getValue()
		{
			return intValue;
		}

		public static OxTypes forValue(int value)
		{
			return getMappings().get(value);
		}
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
		ER_LAST;

		public int getValue()
		{
			return this.ordinal();
		}

		public static OxErrors forValue(int value)
		{
			return values()[value];
		}
	}

	public enum OxWarnings
	{
		WR_DECFAILED,
		WR_ITMAX,
		WR_CONCAT,
		WR_CASTINT,
		WR_VECIDXMAT,
		WR_DETERMINANT,
		WR_USER;

		public int getValue()
		{
			return this.ordinal();
		}

		public static OxWarnings forValue(int value)
		{
			return values()[value];
		}
	}

}
 