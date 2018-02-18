/*--------------------------------------------------------------------------
 * oxmetrics.h - API for client interface with OxMetrics
 *
 *       (C) Jurgen Doornik 1995-2000
 *
 *--------------------------------------------------------------------------*/

#ifndef INC_OXMETRICS_H
#define INC_OXMETRICS_H

#include "jdtypes.h"

#ifdef __cplusplus
extern "C" {
#endif

enum DrawQQTypes
{   QQ_CHI, QQ_F, QQ_N, QQ_T, QQ_U, QQ_N_SE
};

/*================== function prototypes ===================================*/
/* old names */
BOOL JDCALL GiveWinStart(const TCHAR *sClient, int iBatchIndex);
BOOL JDCALL GiveWinStartA(const char *sClient, int iBatchIndex);
BOOL JDCALL GiveWinStartEx(const TCHAR *sClient, int iBatchIndex, int iGiveWinVersion);
BOOL JDCALL GiveWinStartExA(const char *sClient, int iBatchIndex, int iGiveWinVersion);
void JDCALL GiveWinInit();
void JDCALL GiveWinFinish();
int  JDCALL GiveWinDLLVersion(void);

/* shell provided */
BOOL JDCALL OxMetricsStart(const TCHAR *sClient, int iBatchIndex);
BOOL JDCALL OxMetricsStartA(const char *sClient, int iBatchIndex);
BOOL JDCALL OxMetricsStartEx(const TCHAR *sClient, int iBatchIndex, int iOxMetricsVersion);
BOOL JDCALL OxMetricsStartExA(const char *sClient, int iBatchIndex, int iOxMetricsVersion);
BOOL JDCALL OxMetricsStartAdvise(void);
BOOL JDCALL OxMetricsStopAdvise(void);
void JDCALL OxMetricsInit();
void JDCALL OxMetricsFinish();
int  JDCALL OxMetricsDLLVersion(void);
BOOL JDCALL ClientIsRegistered(const TCHAR *sClient, int iVersion100);
BOOL JDCALL ClientIsRegisteredA(const char *sClient, int iVersion100);
void JDCALL SetBatchClient(void);
void JDCALL SetOxMetricsBlock(int iBlock);
void JDCALL SetOxMetricsBlockOnShowDialog(int iShow);
void JDCALL SendOxMetricsMsg(const TCHAR *sMessage);
void JDCALL SendOxMetricsMsgA(const char *sMessage);
void JDCALL SetClassName(const TCHAR *sClassName);
void JDCALL SetClassNameA(const char *sClassName);
void JDCALL SetPackageName(const TCHAR *sPackageName);
void JDCALL SetPackageNameA(const char *sPackageName);

int  JDCALL FindFirstText(const TCHAR *sText, int iOption);
int  JDCALL FindFirstTextA(const char *sText, int iOption);
void JDCALL FocusTextWindow(void);
void JDCALL SetTextMarker(int iMode);
void JDCALL SetTextShowMode(int iMode);
void JDCALL SetTextWindow(const TCHAR *sTitle);
void JDCALL SetTextWindowA(const char *sTitle);
void JDCALL ShowTextWindow(void);
void JDCALL PrintText(const TCHAR *sResults);
void JDCALL PrintTextA(const char *sResults);
void JDCALL SpoolText(void);

void JDCALL BatchDone(const TCHAR *sResults);
void JDCALL BatchDoneA(const char *sResults);
void JDCALL SetBatchCommands(const TCHAR *sResults);
void JDCALL SetBatchCommandsA(const char *sResults);
void JDCALL SetBatchCode(const TCHAR *sResults);
void JDCALL SetBatchCodeA(const char *sResults);
void JDCALL AddBatchCode(const TCHAR *sResults);
void JDCALL AddBatchCodeA(const char *sResults);
void JDCALL SetOxCode(const TCHAR *sResults);
void JDCALL SetOxCodeA(const char *sResults);
void JDCALL AddOxCode(const TCHAR *sResults);
void JDCALL AddOxCodeA(const char *sResults);

void JDCALL ShowHelp(const TCHAR *sPackage, const TCHAR *sPath, const TCHAR *sCommand);
void JDCALL ShowHelpA(const char *sPackage, const char *sPath, const char *sCommand);
int  JDCALL RequestAdvise(TCHAR *sResults, int cResults);
int  JDCALL RequestAdviseA(char *sResults, int cResults);
void JDCALL OxMetricsFileOpen(const TCHAR *sFileName);
void JDCALL OxMetricsFileOpenA(const char *sFileName);

BOOL JDCALL CreateDb(const TCHAR *sDbName, int iFreq, int iYear1, int iPeriod1, int iYear2, int iPeriod2, int iOptions);
BOOL JDCALL CreateDbA(const char *sDbName, int iFreq, int iYear1, int iPeriod1, int iYear2, int iPeriod2, int iOptions);
int  JDCALL GetDbSample(SAMPLE *pSample);
BOOL JDCALL GetDbVarData(VECTOR vX, int cX, const TCHAR *sX, int *pcT);
BOOL JDCALL GetDbVarDataA(VECTOR vX, int cX, const char *sX, int *pcT);
BOOL JDCALL GetDbVar(VECTOR vX, int iX, int iT1, int iT2);
int  JDCALL GetDbVarRaw(const TCHAR *sX, TCHAR *sXraw, double *pdScale,
	int *piDiff, int *piSeasDiff, int *piTrans);
int  JDCALL GetDbVarRawA(const char *sX, char *sXraw, double *pdScale,
	int *piDiff, int *piSeasDiff, int *piTrans);
int  JDCALL GetDbVarIndex(const TCHAR *sX, int* piLagAlways0, TCHAR *sXCopy);
int  JDCALL GetDbVarIndexA(const char *sX, int* piLagAlways0, char *sXCopy);
int  JDCALL GetDbVarIndexDirect(const TCHAR *sX);
int  JDCALL GetDbVarIndexDirectA(const char *sX);
int  JDCALL GetDbVarIndexWithLag(const TCHAR *sX, int* piLag, TCHAR *sVarNameInDb);
int  JDCALL GetDbVarIndexWithLagA(const char *sX, int* piLag, char *sVarNameInDb);
int  JDCALL GetDbVarType(int iX);
void JDCALL GetDefaultDbName(TCHAR *sDbNameDest);
void JDCALL GetDefaultDbNameA(char *sDbNameDest);
BOOL JDCALL GetFirstDbName(TCHAR *sDb);
BOOL JDCALL GetFirstDbNameA(char *sDb);
BOOL JDCALL GetNextDbName(TCHAR *sDb);
BOOL JDCALL GetNextDbNameA(char *sDb);
BOOL JDCALL GetFirstDbVarName(TCHAR *sVar);
BOOL JDCALL GetFirstDbVarNameA(char *sVar);
BOOL JDCALL GetNextDbVarName(TCHAR *sVar);
BOOL JDCALL GetNextDbVarNameA(char *sVar);
BOOL JDCALL SelectDb(const TCHAR *sDbName, SAMPLE *pSample);
BOOL JDCALL SelectDbA(const char *sDbName, SAMPLE *pSample);
void JDCALL SetDbVar(VECTOR vX, int cX, int iT1, const TCHAR *sX, BOOL fQuery);
void JDCALL SetDbVarA(VECTOR vX, int cX, int iT1, const char *sX, BOOL fQuery);
void JDCALL SetDbVarEx(VECTOR vX, int cX, int iT1, const TCHAR *sX, BOOL bQuery, TCHAR *sVarNameInDb);
void JDCALL SetDbVarExA(VECTOR vX, int cX, int iT1, const char *sX, BOOL bQuery, char *sVarNameInDb);
void JDCALL GetDbFullPath(TCHAR *sDbFullPathDest);
void JDCALL GetDbFullPathA(char *sDbFullPathDest);

void JDCALL CloseDrawWindow(void);
void JDCALL FocusDrawWindow(void);
void JDCALL OpenDrawWindow(void);
void JDCALL SaveDrawWindow(const TCHAR *sFile);
void JDCALL SaveDrawWindowA(const char *sFile);
void JDCALL SetDrawShowMode(int iMode);
void JDCALL SetDrawWindow(const TCHAR *sTitle);
void JDCALL SetDrawWindowA(const char *sTitle);
void JDCALL ShowDrawWindowAdd(void);
void JDCALL ShowDrawWindow(void);

int  JDCALL Draw(const TCHAR *sAction, int iArea, VECTOR vY, int cY, const TCHAR *sY,
	VECTOR vX, int cX, const TCHAR *sX, VECTOR vZ, int cZ, const TCHAR *sZ,
	int iSymbol, int iIndex, int *piArgInt, int cArgInt, double *pdArgDbl, int cArgDbl);
int  JDCALL DrawA(const char *sAction, int iArea, VECTOR vY, int cY, const char *sY,
	VECTOR vX, int cX, const char *sX, VECTOR vZ, int cZ, const char *sZ,
	int iSymbol, int iIndex, int *piArgInt, int cArgInt, double *pdArgDbl, int cArgDbl);
int  JDCALL DrawA_NoRet(const char *sAction, int iArea, VECTOR vY, int cY, const char *sY,
	VECTOR vX, int cX, const char *sX, VECTOR vZ, int cZ, const char *sZ,
	int iSymbol, int iIndex, int *piArgInt, int cArgInt, double *pdArgDbl, int cArgDbl);

void JDCALL DrawAdjust(int iObjectNo, int iType, double d1, double d2, double d3, double d4);
void JDCALL SetDrawOptions(int iOption, int i1, int i2, int i3, int i4, int i5);

int  JDCALL IDrawAcf(int iArea, VECTOR vY, int cT, int iLength, const TCHAR *sY,
	BOOL bAcf, BOOL bPacf, BOOL bErrorBand, int iLineno, BOOL fBar);
int  JDCALL IDrawAcfA(int iArea, VECTOR vY, int cT, int iLength, const char *sY,
	BOOL bAcf, BOOL bPacf, BOOL bErrorBand, int iLineno, BOOL fBar);
void JDCALL DrawBoxPlot(int iArea, VECTOR vY, int cY, const TCHAR *sY, int iLineno);
void JDCALL DrawBoxPlotA(int iArea, VECTOR vY, int cY, const char *sY, int iLineno);
void JDCALL DrawCorrelogram(int iArea, VECTOR vY, int cY, const TCHAR *sY, int iOrder, int iLineno);
void JDCALL DrawCorrelogramA(int iArea, VECTOR vY, int cY, const char *sY, int iOrder, int iLineno);
void JDCALL DrawDensity(int iArea, VECTOR vY, int cY, const TCHAR *sY,
	BOOL fDens, BOOL fHist, BOOL fNormal, BOOL fCdf, BOOL fStand, int cBar, int iIndex);
void JDCALL DrawDensityA(int iArea, VECTOR vY, int cY, const char *sY,
	BOOL fDens, BOOL fHist, BOOL fNormal, BOOL fCdf, BOOL fStand, int cBar, int iIndex);
int  JDCALL IDrawHistogram(int iArea, VECTOR vBar, int cBar, 
    double dMin, double dStep, int iIndex, int iColorIn);	 
void JDCALL DrawQQ(int iArea, VECTOR vY, int cY, const TCHAR *sY,
	int iDist, double dDf1, double dDf2);
void JDCALL DrawQQA(int iArea, VECTOR vY, int cY, const char *sY,
	int iDist, double dDf1, double dDf2);
void JDCALL DrawSpectrum(int iArea, VECTOR vY, int cY, const TCHAR *sY, int iOrder, int iLineno);
void JDCALL DrawSpectrumA(int iArea, VECTOR vY, int cY, const char *sY, int iOrder, int iLineno);

int  JDCALL IDrawVector(int iArea, VECTOR vY, int cY, const TCHAR *sY,
        double dXfirst, double dXstep, int iSymbol, int iIndex);
int  JDCALL IDrawVectorA(int iArea, VECTOR vY, int cY, const char *sY,
        double dXfirst, double dXstep, int iSymbol, int iIndex);
int  JDCALL IDrawTimeVector(int iArea, VECTOR vY, int cY, const TCHAR *sY,
        int mnYear, int mnPeriod, int iFreq, int iSymbol, int iIndex);
int  JDCALL IDrawTimeVectorA(int iArea, VECTOR vY, int cY, const char *sY,
        int mnYear, int mnPeriod, int iFreq, int iSymbol, int iIndex);
int  JDCALL IDrawDateVector(int iArea, VECTOR vY, int cY, const TCHAR *sY,
        VECTOR vX, const TCHAR *sX, int iSymbol, int iIndex);
int  JDCALL IDrawDateVectorA(int iArea, VECTOR vY, int cY, const char *sY,
        VECTOR vX, const char *sX, int iSymbol, int iIndex);
int  JDCALL IDrawXVector(int iArea, VECTOR vY, int cY, const TCHAR *sY,
        VECTOR vX, const TCHAR *sX, int iSymbol, int iIndex);
int  JDCALL IDrawXVectorA(int iArea, VECTOR vY, int cY, const char *sY,
        VECTOR vX, const char *sX, int iSymbol, int iIndex);
int  JDCALL IDraw3D(int iArea, VECTOR vY, int cY, const TCHAR *sY,
	VECTOR vX, int cX, const TCHAR *sX, VECTOR vZ, int cZ,
	const TCHAR *sZ, int iType, int iPalette, int iIndex, int iMode);
int  JDCALL IDraw3DA(int iArea, VECTOR vY, int cY, const char *sY,
	VECTOR vX, int cX, const char *sX, VECTOR vZ, int cZ,
	const char *sZ, int iType, int iPalette, int iIndex, int iMode);
void JDCALL DrawZ(VECTOR vZ, int cZ, const TCHAR *sZ, int iMode, double dFac, int iIndex);
void JDCALL DrawZA(VECTOR vZ, int cZ, const char *sZ, int iMode, double dFac, int iIndex);
int  JDCALL IDrawLine(int iArea, double dX1, double dY1, double dX2, double dY2, int iIndex);
int  JDCALL IDrawSymbol(int iArea, double dX1, double dY1, double dX2, double dY2, int iSymType, int iIndex);
int  JDCALL IDrawSymbol3(int iArea, double dX1, double dY1, double dZ1, double dX2, double dY2, double dZ2, int iSymType, int iIndex);
int  JDCALL IDrawPLine(int iArea, int iX1, int iY1, int iX2, int iY2, int iIndex);
int  JDCALL IDrawPSymbol(int iArea, int iX1, int iY1, int iX2, int iY2, int iSymType, int iIndex);
int  JDCALL IDrawText(int iArea, const TCHAR *sText, double dX1, double dY1, int iFontNo, int iFontSize, int iTitle);
int  JDCALL IDrawTextA(int iArea, const char *sText, double dX1, double dY1, int iFontNo, int iFontSize, int iTitle);
int  JDCALL IDrawText3(int iArea, const TCHAR *sText, double dX1, double dY1, double dZ1, int iFontNo, int iFontSize, int iTitle, int iRotation);
int  JDCALL IDrawText3A(int iArea, const char *sText, double dX1, double dY1, double dZ1, int iFontNo, int iFontSize, int iTitle, int iRotation);
int  JDCALL IDrawPText(int iArea, const char *sText, int iX1, int iY1, int iFontNo, int iFontSize, int iTitle, int iRotation);
int  JDCALL IDrawPTextA(int iArea, const TCHAR *sText, int iX1, int iY1, int iFontNo, int iFontSize, int iTitle, int iRotation);
int  JDCALL IDrawTitle(int iArea, const TCHAR *sText);
int  JDCALL IDrawTitleA(int iArea, const char *sText);
int  JDCALL IDrawPLegend(int iArea, int iOffsetX, int iOffsetY, BOOL fHidden);
int  JDCALL IDrawAxis(int iArea, BOOL fIsXAxis, double dAnchor, double dAxmin, double dAxmax, double dFirstLarge, double dLargeStep, double dSmallStep, int iFreq);
int  JDCALL IDrawAxis3(int iArea, int iIsXAxis, double dAnchor, double dAnchor2, double dAxmin, double dAxmax, double dFirstLarge, double dLargeStep, double dSmallStep, int iFreq);
int  JDCALL IDrawAxisAuto(int iArea, BOOL fIsXAxis, BOOL fShow, int iAnchor, double dAnchor);
int  JDCALL IDrawAxisAuto3(int iArea, int iIsXAxis, BOOL fShow, int iAnchor, double dAnchor, double dAnchor2);

/* others */
void JDCALL BufferSWprintf(BOOL fStartBuffering);
void JDCALL FlushSWprintf(void);
TCHAR * JDCALLC SWprintf(const TCHAR *, ...);
void JDCALL SWputs(const TCHAR *sText);
void JDCALL SWputsA(const char *sText);

/* backward compatibility */
void JDCALL SWputsResults(const TCHAR *sResults);
void JDCALL SWputsResultsA(const char *sResults);
void JDCALL SWspoolResults(void);

#ifdef __cplusplus
}
#endif

#endif  /* INC_OXMETRICS_H */

