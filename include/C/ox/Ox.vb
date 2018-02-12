'converted from Ox.cs by http://www.developerfusion.com/tools/convert/csharp-to-vb/
'
Imports System.Runtime.InteropServices

NotInheritable Class Ox
	Private Sub New()
	End Sub

	Const OXLIB As String = "oxwin"

	' basic mathematical and statistical functions and allocation/deallocation
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function c_abs(xr As Double, xi As Double) As Double
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function c_div(xr As Double, xi As Double, yr As Double, yi As Double, ByRef zr As Double, ByRef zi As Double) As Int32
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Sub c_erf(x As Double, y As Double, ByRef erfx As Double, ByRef erfy As Double)
	End Sub
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Sub c_exp(xr As Double, xi As Double, ByRef yr As Double, ByRef yi As Double)
	End Sub
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Sub c_log(xr As Double, xi As Double, ByRef yr As Double, ByRef yi As Double)
	End Sub
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Sub c_mul(xr As Double, xi As Double, yr As Double, yi As Double, ByRef zr As Double, ByRef zi As Double)
	End Sub
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Sub c_sqrt(xr As Double, xi As Double, ByRef yr As Double, ByRef yi As Double)
	End Sub
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function DBessel01(x As Double, iType As Int32, n As Int32) As Double
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function DBesselNu(x As Double, iType As Int32, n As Double) As Double
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function DBetaFunc(dX As Double, dA As Double, dB As Double) As Double
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function DDawson(dX As Double) As Double
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function DDensBeta(x As Double, a As Double, b As Double) As Double
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function DDensChi(x As Double, dDf As Double) As Double
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function DDensF(x As Double, dDf1 As Double, dDf2 As Double) As Double
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function DDensGamma(g As Double, r As Double, a As Double) As Double
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function DDensGH(dX As Double, dNu As Double, dDelta As Double, dGamma As Double, dBeta As Double) As Double
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function DDensGIG(dX As Double, dNu As Double, dDelta As Double, dGamma As Double) As Double
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function DDensMises(x As Double, dMu As Double, dKappa As Double) As Double
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function DDensNormal(x As Double) As Double
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function DDensPoisson(dMu As Double, k As Int32) As Double
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function DDensT(x As Double, dDf As Double) As Double
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function DDiagXSXt(iT As Int32, mX As IntPtr, mS As IntPtr, cS As Int32) As Double
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function DDiagXtSXtt(cX As Int32, mXt As IntPtr, mS As Int32, cS As Int32) As Double
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Sub DecQRtMul(mQt As IntPtr, cX As Int32, cT As Int32, mY As IntPtr, cY As Int32, cR As Int32)
	End Sub
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Sub DecQRtMult(mQt As IntPtr, cX As Int32, cT As Int32, mYt As IntPtr, cY As Int32, cR As Int32)
	End Sub
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function DErf(z As Double) As Double
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function DExpInt(z As Double) As Double
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function DExpInt1(z As Double) As Double
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function DExpInte(z As Double) As Double
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function DGamma(z As Double) As Double
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function DGammaFunc(dX As Double, dR As Double) As Double
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function DGetInvertEps() As Double
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function DGetInvertEpsNorm(mA As IntPtr, cA As Double) As Double
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function DLogGamma(dA As Double) As Double
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function DPolyGamma(dA As Double, n As Int32) As Double
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function DProbBeta(x As Double, a As Double, b As Double) As Double
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function DProbBVN(dLo1 As Double, dLo2 As Double, dRho As Double) As Double
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function DProbChi(x As Double, dDf As Double) As Double
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function DProbChiNc(x As Double, df As Double, dNc As Double) As Double
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function DProbF(x As Double, dDf1 As Double, dDf2 As Double) As Double
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function DProbFNc(x As Double, dDf1 As Double, dDf2 As Double, dNc As Double) As Double
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function DProbGamma(x As Double, dR As Double, dA As Double) As Double
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function DProbMises(x As Double, dMu As Double, dKappa As Double) As Double
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function DProbNormal(x As Double) As Double
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function DProbPoisson(dMu As Double, k As Int32) As Double
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function DProbT(x As Double, iDf As Int32) As Double
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function DProbTNc(x As Double, dDf As Double, dNc As Double) As Double
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function DQuanBeta(x As Double, a As Double, b As Double) As Double
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function DQuanChi(p As Double, dDf As Double) As Double
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function DQuanF(p As Double, dDf1 As Double, dDf2 As Double) As Double
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function DQuanGamma(p As Double, r As Double, a As Double) As Double
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function DQuanMises(p As Double, dMu As Double, dKappa As Double) As Double
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function DQuanNormal(p As Double) As Double
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function DQuanT(p As Double, iDf As Int32) As Double
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function DQuanTD(p As Double, dDf As Double) As Double
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function DRanBeta(a As Double, b As Double) As Double
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function DRanChi(dDf As Double) As Double
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function DRanExp(dLambda As Double) As Double
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function DRanF(dDf1 As Double, dDf2 As Double) As Double
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function DRanGamma(dR As Double, dA As Double) As Double
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function DRanGIG(dNu As Double, dDelta As Double, dGamma As Double) As Double
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function DRanInvGaussian(dMu As Double, dLambda As Double) As Double
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function DRanLogistic() As Double
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function DRanLogNormal() As Double
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function DRanMises(dKappa As Double) As Double
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function DRanNormalPM() As Double
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function DRanStable(dA As Double, dB As Double) As Double
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function DRanStudentT(dDf As Double) As Double
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function DRanT(iDf As Int32) As Double
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function DRanU() As Double
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function DTailProbChi(x As Double, dDf As Double) As Double
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function DTailProbF(x As Double, dDf1 As Double, dDf2 As Double) As Double
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function DTailProbGamma(x As Double, dR As Double, dA As Double) As Double
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function DTailProbNormal(x As Double) As Double
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function DTailProbT(x As Double, iDf As Int32) As Double
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function DTrace(mat As IntPtr, cA As Int32) As Double
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function DTraceAB(mA As IntPtr, mB As IntPtr, cM As Int32, cN As Int32) As Double
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function DVecsum(vA As Double(), cA As Int32) As Double
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Sub EigVecDiv(ByRef mmE As IntPtr, vEr As Double(), vEi As Double(), cA As Int32)
	End Sub
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function FCubicSpline(vY As Double(), vT As Double(), cT As Int32, ByRef pdAlpha As Double, vG As Double(), vX As Double(), _
		ByRef pdCV As Double, ByRef pdPar As Double, fAuto As Int32, iDesiredPar As Int32) As Int32
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function FCubicSplineTime(vY As Double(), cT As Double, dAlpha As Double, ByRef vG As Int32, fHP As Int32) As Int32
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function FFT1d(mDest As IntPtr, mSrc As IntPtr, iDir As Int32, isComplex As Int32) As Int32
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Sub FftComplex(vXr As Double(), vXi As Double(), iPower As Int32, iDir As Int32)
	End Sub
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Sub FftReal(vXr As Double(), vXi As Double(), iPower As Int32, iDir As Int32)
	End Sub
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function FftDiscrete(vXr As Double(), vXi As Double(), cN As Int32, iDir As Int32) As Int32
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function FGetAcf(vX As Double(), cT As Int32, cLag As Int32, vAcf As Double()) As Int32
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function FGetAcfRun(vX As Double(), cT As Int32, cLag As Int32, vAcf As Double()) As Int32
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function FIsInf(d As Double) As Int32
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function FIsNaN(d As Double) As Int32
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function FPeriodogram(vX As Double(), cT As Int32, iTrunc As Int32, cS As Int32, vS As Double(), iMode As Int32) As Int32
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function FPeriodogramAcf(vAcf As Double(), cT As Int32, iTrunc As Int32, cS1 As Int32, vS As Double(), iMode As Int32, _
		cTwin As Int32) As Int32
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function FPPtDec(mA As IntPtr, cA As Int32) As Object
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function IDecQRt(mXt As IntPtr, cX As Int32, cT As Int32, ByRef piPiv As Int32, ByRef pcR As Int32) As Int32
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function IDecQRtEx(mXt As IntPtr, cX As Int32, cT As Int32, ByRef piPiv As Int32, vTau As Double()) As Int32
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function IDecQRtRank(mXt As IntPtr, cX As Int32, cT As Int32, ByRef pcR As Int32) As Int32
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function IDecSVD(mA As IntPtr, cM As Int32, cN As Int32, vW As Double(), fDoU As Int32, mU As IntPtr, _
		fDoV As Int32, mV As IntPtr, fSort As Int32) As Int32
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function IEigen(mA As IntPtr, cA As Int32, vEr As Double(), vEi As Double(), mmE As IntPtr) As Int32
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function IEigenSym(mA As IntPtr, cA As Int32, vEval As Double(), fDoVectors As Int32) As Int32
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function IEigValPoly(vPoly As Double(), vEr As Double(), vEi As Double(), cA As Int32) As Int32
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function IEigValReal(mA As IntPtr, vEr As Double(), vEi As Double(), cA As Int32) As Int32
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function IEigValSym(mA As IntPtr, vEv As Double(), cA As Int32) As Int32
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function IEigVecReal(mA As IntPtr, vEr As Double(), vEi As Double(), cA As Int32) As Int32
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function IEigVecSym(mA As IntPtr, vEv As Double(), cA As Int32) As Int32
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function IGenEigVecSym(mA As IntPtr, mB As IntPtr, vEval As Double(), vSubd As Double(), cA As Int32) As Int32
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function IGetAcf(vX As Double(), cT As Int32, cLag As Int32, vAcf As Double(), bCov As Int32) As Int32
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function IInvDet(mA As IntPtr, cA As Int32, ByRef pdLogDet As Double, ByRef piSignDet As Int32) As Int32
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function IInvert(mA As IntPtr, cA As Int32) As Int32
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function ILDLbandDec(mA As IntPtr, vD As Double(), cB As Int32, cA As Int32) As Int32
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function ILDLdec(mA As IntPtr, vD As Double(), cA As Int32) As Int32
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function ILUPdec(mA As IntPtr, cA As Int32, ByRef piPiv As Int32, ByRef pdLogDet As Double, ByRef piSignDet As Int32, mUt As IntPtr) As Int32
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function ILUPlogdet(mU As IntPtr, cA As Int32, ByRef piPiv As Int32, dNormEps As Double, ByRef pdLogDet As Double) As Int32
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function IMatRank(mA As IntPtr, cM As Int32, cN As Int32, dEps As Double, bAbsolute As Int32) As Int32
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function IntMatAllocBlock(cR As Int32, cC As Int32) As IntPtr
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Sub IntMatFreeBlock(m As IntPtr)
	End Sub
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function INullSpace(mA As IntPtr, cM As Int32, cN As Int32, fAppend As Int32) As Int32
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function IOlsNorm(mXt As IntPtr, cX As Int32, cT As Int32, mYt As IntPtr, cY As Int32, mB As IntPtr, _
		mXtXinv As IntPtr, mXtX As IntPtr, fInRows As Int32) As Int32
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function IOlsQR(mXt As IntPtr, cX As Int32, cT As Int32, mYt As IntPtr, cY As Int32, mB As IntPtr, _
		mXtXinv As IntPtr, mXtX As IntPtr, ByRef vW As Int32) As Int32
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function IRanBinomial(n As Int32, p As Double) As Int32
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function IRanLogarithmic(dA As Double) As Int32
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function IRanNegBin(iiN As Int32, dP As Double) As Int32
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function IRanPoisson(dMu As Double) As Int32
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function ISVDdec(mA As IntPtr, cM As Int32, cN As Int32, vW As Double(), fDoU As Int32, mU As IntPtr, _
		fDoV As Int32, mV As IntPtr, v_1 As Double(), fSort As Int32) As Int32
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function ISymInv(mA As IntPtr, cA As Int32) As Int32
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function ISymInvDet(mA As IntPtr, cA As Int32, ByRef pdLogDet As Double) As Int32
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Sub LDLbandSolve(mL As IntPtr, vD As Double(), vX As Double(), vB As Double(), cB As Int32, cA As Int32)
	End Sub
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Sub LDLsolve(mL As IntPtr, vD As Double(), vX As Double(), vB As Double(), cA As Int32)
	End Sub
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Sub LDLsolveInv(mLDLt As IntPtr, mAinv As IntPtr, cA As Int32)
	End Sub
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Sub LUPsolve(mL As IntPtr, mU As IntPtr, ByRef piPiv As Int32, vB As Double(), cA As Int32)
	End Sub
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Sub LUPsolveInv(mL As IntPtr, mU As IntPtr, ByRef piPiv As Int32, mAinv As IntPtr, cA As Int32)
	End Sub
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function MatAB(mA As IntPtr, cA As Int32, cC As Int32, mB As IntPtr, cB As Int32, mAB As IntPtr) As Int32
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function MatABt(mA As IntPtr, cA As Int32, cC As Int32, mB As IntPtr, cB As Int32, mABt As IntPtr) As Int32
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function MatAcf(mAcf As IntPtr, mX As IntPtr, cT As Int32, cX As Int32, mxLag As Int32) As Int32
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function MatAdd(mA As IntPtr, cM As Int32, cN As Int32, mB As IntPtr, dFac As Double, mAplusB As IntPtr) As Int32
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function MatAllocBlock(cR As Int32, cC As Int32) As IntPtr
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function MatAtB(mA As IntPtr, cA As Int32, cC As Int32, mB As IntPtr, cB As Int32, mAtB As IntPtr) As Int32
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function MatBBt(mB As IntPtr, cB As Int32, cS As Int32, mBBt As IntPtr) As Int32
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function MatBSBt(mB As IntPtr, cB As Int32, mS As IntPtr, cS As Int32, mBSBt As IntPtr) As Int32
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function MatBtB(mB As IntPtr, cB As Int32, cS As Int32, mBtB As Int32) As Int32
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function MatBtBVec(mB As IntPtr, cB As Int32, cS As Int32, vY As Double(), mBtB As IntPtr) As Int32
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function MatBtSB(mB As IntPtr, cB As Int32, mS As IntPtr, cS As Int32, mBtSB As IntPtr) As Int32
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function MatCopy(mDest As IntPtr, mSrc As IntPtr, cM As Int32, cN As Int32) As Int32
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function MatCopyTranspose(mDestT As IntPtr, mSrc As IntPtr, cM As Int32, cN As Int32) As Int32
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Sub MatCopyVecc(mDest As IntPtr, vSrc_c As Double(), cM As Int32, cN As Int32)
	End Sub
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Sub MatCopyVecr(mDest As IntPtr, vSrc_r As Double(), cM As Int32, cN As Int32)
	End Sub
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function MatDup(mSrc As IntPtr, cM As Int32, cN As Int32) As Int32
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Sub MatFreeBlock(m As IntPtr)
	End Sub
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function MatGenInvert(mA As IntPtr, cM As Int32, cN As Int32, mRes As IntPtr, vSval As Double()) As Int32
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function MatGetAt(mSrc As IntPtr, i As Int32, j As Int32) As Double
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function MatI(mDest As IntPtr, cM As Int32) As Int32
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function MatNaN(mDest As IntPtr, cM As Int32, cN As Int32) As Int32
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function MatPartAcf(mPartAcf As IntPtr, mAcf As IntPtr, cAcf As Int32, mY As IntPtr, cY As Int32, ByRef pdLogDet As Double, _
		bFilter As Int32) As Int32
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function MatRan(mA As IntPtr, cR As Int32, cC As Int32) As Int32
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function MatRann(mA As IntPtr, cR As Int32, cC As Int32) As Int32
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function MatReflect(mA As IntPtr, cA As Int32) As Int32
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Sub MatSetAt(mDest As IntPtr, d As Double, i As Int32, j As Int32)
	End Sub
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function MatStandardize(mXdest As IntPtr, mX As IntPtr, cT As Int32, cX As Int32) As Int32
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function MatTranspose(mA As IntPtr, cA As Int32) As Int32
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function MatVariance(mXtX As IntPtr, mX As IntPtr, cT As Int32, cX As Int32, fCorr As Int32) As Int32
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function MatZero(mDest As IntPtr, cM As Int32, cN As Int32) As Int32
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Sub OlsQRacc(mXt As IntPtr, cX As Int32, cT As Int32, piPiv As Int32(), cR As Int32, vTau As Double(), _
		mYt As IntPtr, cY As Int32, mB As IntPtr, mXtXinv As IntPtr, mXtX As IntPtr)
	End Sub
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Sub RanDirichlet(vX As Double(), vAlpha As Double(), cALpha As Int32)
	End Sub
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function RanGetSeed(piSeed As Int32(), cSeed As Int32) As Int32
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Sub RanInit()
	End Sub
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Sub RanNewRan(fnDRanu As IntPtr, fnRanSetSeed As IntPtr, fnRanGetSeed As IntPtr)
	End Sub
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Sub RanSetRan(sRan As String)
	End Sub
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Sub RanSetSeed(piSeed As Int32(), cSeed As Int32)
	End Sub
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Sub RanSubSample(vU As Double(), cU As Int32, cN As Int32)
	End Sub
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Sub RanUorder(vU As Double(), cU As Int32)
	End Sub
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Sub RanWishart(mX As IntPtr, cX As Int32)
	End Sub
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Sub SetFastMath(fYes As Int32)
	End Sub
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Sub SetInf(ByRef pd As Double)
	End Sub
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Sub SetInvertEps(dEps As Double)
	End Sub
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Sub SetNaN(ByRef pd As Double)
	End Sub
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Sub ToeplitzSolve(vR As Double(), cR As Int32, cM As Int32, mB As IntPtr, cB As Int32, v_1 As Double())
	End Sub
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Sub VeccCopyMat(vDest_c As Double(), mSrc As IntPtr, cM As Int32, cN As Int32)
	End Sub
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function VecDiscretize(vY As Double(), cY As Int32, dMin As Double, dMax As Double, vDisc As Double(), cM As Int32, _
		vT As Double(), iOption As Int32) As Int32
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function VecDup(vSrc As Double(), cM As Int32) As Int32
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Sub VecrCopyMat(vDest_r As Double(), mSrc As IntPtr, cM As Int32, cN As Int32)
	End Sub
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function VecTranspose(vA As Double(), cM As Int32, cN As Int32) As Int32
	End Function

	' Ox run-time functionality
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function FOxCallBack(pvFunc As IntPtr, rtn As IntPtr, pv As IntPtr, cArg As Int32) As Int32
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function FOxCallBackMember(pvClass As IntPtr, sMember As String, rtn As IntPtr, pv As IntPtr, cArg As Int32) As Int32
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function FOxCreateObject(sClass As String, rtn As IntPtr, pv As IntPtr, cArg As Int32) As Int32
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function FOxGetDataMember(pvClass As IntPtr, sMember As String, rtn As IntPtr) As Int32
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function FOxLibAddFunction(sFunc As String, pFunc As IntPtr, fVarArg As Int32) As Int32
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function FOxLibAddFunctionEx(sFunc As String, pFunc As IntPtr, cArgs As Int32, flFlags As Int32) As Int32
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function FOxRun(iMainIP As Int32, sFunc As String) As Int32
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function FOxSetDataMember(pvClass As IntPtr, sMember As String, pv As IntPtr) As Int32
	End Function

	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function IOxRunInit() As Int32
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function IOxVersion() As Int32
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function IOxVersionIsProfessional() As Int32
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function IOxVersionOxo() As Int32
	End Function

	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Sub OxCloneObject(rtn As IntPtr, pvObject As IntPtr, bDeep As Int32)
	End Sub
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Sub OxDeleteObject(pvClass As IntPtr)
	End Sub
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Sub OxFnDouble(rtn As IntPtr, pv As IntPtr, fn As IntPtr)
	End Sub
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Sub OxFnDouble2(rtn As IntPtr, pv As IntPtr, fn As IntPtr)
	End Sub
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Sub OxFnDouble3(rtn As IntPtr, pv As IntPtr, fn As IntPtr)
	End Sub
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Sub OxFnDouble4(rtn As IntPtr, pv As IntPtr, fn As IntPtr)
	End Sub
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Sub OxFnDoubleInt(rtn As IntPtr, pv As IntPtr, fn As IntPtr)
	End Sub
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Sub OxFreeByValue(pv As IntPtr)
	End Sub
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Sub OxGetMainArgs(pcArgc As Int32(), ByRef pasArgv As String())
	End Sub
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Sub OxGetOxArgs(pcArgc As Int32(), ByRef pasArgv As String())
	End Sub
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function OxGetPrintlevel() As Int32
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function OxGetUserExitCode() As Int32
	End Function

	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Sub OxLibArgError(iArg As Int32)
	End Sub
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Sub OxLibArgTypeError(iArg As Int32, iExpected As Int32, iFound As Int32)
	End Sub
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Sub OxLibCheckArrayMatrix(pv As IntPtr, iFirst As Int32, iLast As Int32, m As IntPtr)
	End Sub
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Sub OxLibCheckMatrixSize(pv As IntPtr, iFirst As Int32, iLast As Int32, r As Int32, c As Int32)
	End Sub
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Sub OxLibCheckSquareMatrix(pv As IntPtr, iFirst As Int32, iLast As Int32)
	End Sub
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Sub OxLibCheckType(iType As Int32, pv As IntPtr, iFirst As Int32, iLast As Int32)
	End Sub
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Sub OxLibValArrayCalloc(pv As IntPtr, c As Int32)
	End Sub
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Sub OxLibValMatDup(pv As IntPtr, mSrc As IntPtr, r As Int32, c As Int32)
	End Sub
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Sub OxLibValMatMalloc(pv As IntPtr, r As Int32, c As Int32)
	End Sub
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Sub OxLibValStrMalloc(pv As IntPtr, c As Int32)
	End Sub
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Sub OxLibValZero(pv As IntPtr)
	End Sub

	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function OxMain(argc As Int32, argv As String()) As Int32
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function OxMainCmd(sCommand As String) As Int32
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Sub OxMainExit()
	End Sub
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Sub OxMainInit()
	End Sub
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Sub OxMakeByValue(pv As IntPtr)
	End Sub
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Sub OxMessage(s As String)
	End Sub
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Sub OxPuts(s As String)
	End Sub

	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Sub OxRunAbort(i As Int32)
	End Sub
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Sub OxRunError(iErno As Int32, sToken As String)
	End Sub
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Sub OxRunErrorMessage(s As String)
	End Sub
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Sub OxRunExit()
	End Sub
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Sub OxRunMainExitCall(fn As IntPtr)
	End Sub
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Sub OxRunMessage(s As String)
	End Sub
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Sub OxRunWarningMessage(sFunc As String, sMsg As String)
	End Sub
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Sub OxSetMainArgs(argc As Int32, argv As String())
	End Sub
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Sub OxSetOxArgs(cArgc As Int32, asArgv As String())
	End Sub
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Sub OxSetPrintlevel(iSet As Int32)
	End Sub
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function OxStoreCreate(c As Int32) As IntPtr
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Sub OxStoreDelete(pv As IntPtr, c As Int32)
	End Sub

	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function OxValColumns(pv As IntPtr) As Int32
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function OxValDuplicate2(pvDest As IntPtr, pvSrc As IntPtr) As Int32
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function OxValGetArray(pv As IntPtr) As IntPtr
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function OxValGetArrayLen(pv As IntPtr) As Int32
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function OxValGetArrayVal(pv As IntPtr, i As Int32) As IntPtr
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function OxValGetBlob(pv As IntPtr, ByRef pI1 As Int32, ByRef pI2 As Int32) As IntPtr
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function OxValGetClassName(pv As IntPtr) As String
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function OxValGetDouble(pv As IntPtr, ByRef pdVal As Double) As Int32
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function OxValGetInt(pv As IntPtr, ByRef piVal As Int32) As Int32
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function OxValGetMat(pv As IntPtr) As IntPtr
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function OxValGetMatc(pv As IntPtr) As Int32
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function OxValGetMatr(pv As IntPtr) As Int32
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function OxValGetMatrc(pv As IntPtr) As Int32
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function OxValGetStaticObject(pv As IntPtr) As IntPtr
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function OxValGetString(pv As IntPtr) As String
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function OxValGetStringCopy(pv As IntPtr, s As String, mxlen As Int32) As Int32
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function OxValGetStringLen(pv As IntPtr) As Int32
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function OxValGetVal(pv As IntPtr, i As Int32) As IntPtr
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function OxValGetVecc(pv As IntPtr, vX As Double()) As Int32
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function OxValGetVecr(pv As IntPtr, vX As Double()) As Int32
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function OxValHasFlag(pv As IntPtr, iFlag As Int32) As Int32
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function OxValHasType(pv As IntPtr, iType As Int32) As Int32
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function OxValRows(pv As IntPtr) As Int32
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Sub OxValSetBlob(pv As IntPtr, i1 As Int32, i2 As Int32, p As IntPtr)
	End Sub
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Sub OxValSetDouble(pv As IntPtr, dVal As Double)
	End Sub
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Sub OxValSetInt(pv As IntPtr, iVal As Int32)
	End Sub
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Sub OxValSetMat(pv As IntPtr, mVal As IntPtr, r As Integer, c As Integer)
	End Sub
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Sub OxValSetMatZero(pv As IntPtr, r As Integer, c As Integer)
	End Sub
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Sub OxValSetNull(pv As IntPtr)
	End Sub
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Sub OxValSetString(pv As IntPtr, sVal As String)
	End Sub
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Sub OxValSetVecc(pv As IntPtr, vX As Double(), r As Integer, c As Integer)
	End Sub
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Sub OxValSetVecr(pv As IntPtr, vX As Double(), r As Integer, c As Integer)
	End Sub
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Sub OxValSetZero(pv As IntPtr)
	End Sub
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function OxValSizec(pv As IntPtr) As Int32
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function OxValSizer(pv As IntPtr) As Int32
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function OxValSizerc(pv As IntPtr) As Int32
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function OxValType(pv As IntPtr) As Int32
	End Function


	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Sub SetOxExit(pfnNewOxExit As IntPtr)
	End Sub
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Sub SetOxGets(pfnNewOxGets As IntPtr)
	End Sub
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Sub SetOxMessage(pfnNewOxMessage As IntPtr)
	End Sub
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Sub SetOxPipe(cPipe As Int32)
	End Sub
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Sub SetOxPuts(pfnNewOxPuts As IntPtr)
	End Sub
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Sub SetOxRunMessage(pfnNewOxRunMessage As IntPtr)
	End Sub
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Sub SetOxTextWindow(pfnNewOxTextWindow As IntPtr)
	End Sub

	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function SOxGetTypeName(iType As Int32) As String
	End Function
	<DllImport(OXLIB, CharSet := CharSet.Ansi, ExactSpelling := True)> _
	Public Shared Function SOxIntFunc() As String
	End Function

	Public Enum OxTypes
		OX_INT = 1
		OX_DOUBLE = 2
		OX_MATRIX = 3
		OX_STRING = 4
		OX_ARRAY = 5
		OX_FUNCTION = 6
		OX_CLASS = 7
		OX_VECTOR = 8
		OX_INTFUNC = 9
		OX_RANGE = 10
		OX_FILE = 11
		OX_IMPORT = 12
		OX_LAMBDA = 13
		OX_BLOB = 14
		OX_RETURN = 64

		OX_NULL = &H100
		OX_VALUE = &H200
		OX_CONST = &H400
		OX_RESERVED = &H800
		OX_EXTERN = &H1000
		OX_GLOBAL = &H2000
		OX_STATDECL = &H4000
		OX_INLINE = &H8000
		OX_KEYWORD = &H8000
		OX_MEMBER = &H10000
		OX_STATIC = &H20000
		OX_VIRTUAL = &H40000
		OX_PUBLIC = &H80000
		OX_INDEX = &H100000
		OX_ADDRESS = &H200000
		OX_ARGUMENT = &H400000
		OX_VARARGS = &H800000
		OX_SERIAL = &H1000000
		OX_VECMAT = &H2000000
		OX_VECRANGE = &H4000000
		OX_IDXSCALAR = &H8000000
		OX_INTERNAL = &H10000000
	End Enum

	Public Enum OxErrors
		ER_RUNTIME
		ER_FATAL
		ER_LIBFUNC
		ER_STACKOVERFLOW
		ER_OM
		ER_NOFUNC
		ER_DIVIDE
		ER_EXPOBJECT
		ER_UNDECLMEM
		ER_SINGULAR
		ER_BADOP
		ER_INCTYPE
		ER_NULL
		ER_ARGS
		ER_BADTYPE
		ER_IDXBOUND
		ER_IDXTYPE
		ER_EXPSCALAR
		ER_EXPSCALMAT
		ER_EXPSCALSTR
		ER_CONSTAS
		ER_EXPCLASS
		ER_UNDEFMEMFU
		ER_NOMEMFUNC
		ER_UNDEFFU
		ER_INVALIDCLASS
		ER_CONSTDEL
		ER_ADDRESS
		ER_NEGATIVE
		ER_SVD
		ER_CAST
		ER_IDXRANGE
		ER_IDXARRAY
		ER_POWOVERFLOW
		ER_ARGSAME
		ER_NOCALLMEM
		ER_EXPFUNC
		ER_ARGUMENT
		ER_IDXMISSING
		ER_SIZENEGATIVE
		ER_ARGERROR
		ER_CALLBACK
		ER_INVALID
		ER_OVERFLOW
		ER_FPEXCEPTION
		ER_CASTINT
		ER_CASTDBL
		ER_CASTMAT
		ER_CASTVEC
		ER_USER1
		ER_USER2
		ER_USER3
		ER_USER4
		ER_USER5
		ER_USER6
		ER_EXPVECTOR
		ER_JSTINVALID
		ER_ARRAYMULTI
		ER_NOPIPE
		ER_USERABORT
		ER_EXPPUBLIC
		ER_EXSTRING
		ER_CONSTUSE
		ER_LAMBDA
		ER_VARCHANGED
		ER_STACK
		ER_FOR
		ER_FOR_INF
		ER_LAST
	End Enum

	Public Enum OxWarnings
		WR_DECFAILED
		WR_ITMAX
		WR_CONCAT
		WR_CASTINT
		WR_VECIDXMAT
		WR_DETERMINANT
		WR_USER
	End Enum
End Class
