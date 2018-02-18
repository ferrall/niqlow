/*--------------------------------------------------------------------------
 * oxexport.h - definitions/declarations Ox exported functions
 *
 *       (C) Jurgen Doornik 1995
 *
 *--------------------------------------------------------------------------*/

#ifndef INC_OXEXPORT
#define INC_OXEXPORT

#include "jdtypes.h"
#include "jdmatrix.h"
#include "oxtypes.h"

#ifdef __cplusplus
extern "C" {
#endif

/*=========================== constants ====================================*/
/*---------------------- error codes and messages --------------------------*/
enum RunErrors                                       /* runtime error codes */
{	ER_RUNTIME,				 ER_FATAL,
	ER_LIBFUNC,
    ER_STACKOVERFLOW,        ER_OM,                   ER_NOFUNC,
    ER_DIVIDE,               ER_EXPOBJECT,            ER_UNDECLMEM,
    ER_SINGULAR,             ER_BADOP,                ER_INCTYPE,
    ER_NULL,                 ER_ARGS,                 ER_BADTYPE,
    ER_IDXBOUND,             ER_IDXTYPE,              ER_EXPSCALAR,
    ER_EXPSCALMAT,           ER_EXPSCALSTR,           ER_CONSTAS,
    ER_EXPCLASS,             ER_UNDEFMEMFU,           ER_NOMEMFUNC,
    ER_UNDEFFU,              ER_INVALIDCLASS,         ER_CONSTDEL,
    ER_ADDRESS,              ER_NEGATIVE,             ER_SVD,
    ER_CAST,                 ER_IDXRANGE,             ER_IDXARRAY,
    ER_POWOVERFLOW,          ER_ARGSAME,              ER_NOCALLMEM,
    ER_EXPFUNC,              ER_ARGUMENT,             ER_IDXMISSING,
	ER_SIZENEGATIVE,		 ER_ARGERROR,			  ER_CALLBACK,
	ER_INVALID,				 ER_OVERFLOW,			  ER_FPEXCEPTION,
	ER_CASTINT,				 ER_CASTDBL,			  ER_CASTMAT,
	ER_CASTVEC,
	ER_USER1,				 ER_USER2,				  ER_USER3,
	ER_USER4,				 ER_USER5,				  ER_USER6,
	ER_EXPVECTOR,			 ER_JSTINVALID,			  ER_ARRAYMULTI,
	ER_NOPIPE,				 ER_USERABORT,			  ER_EXPPUBLIC,
	ER_EXSTRING,			 ER_CONSTUSE,			  ER_LAMBDA,
	ER_VARCHANGED,			 ER_STACK,				  ER_FOR,
	ER_FOR_INF,
	ER_LAST
};
enum RunWarnings
{	WR_DECFAILED,            WR_ITMAX,                WR_CONCAT,
	WR_CASTINT,				 WR_VECIDXMAT,			  WR_DETERMINANT,
	WR_USER
};

/*========================== macro's =======================================*/

/*===================== typedef for handler functions ======================*/
typedef int	 OXCALL	FN_NewOxDraw(const char *sAction, int iArea, VECTOR vY, int cY, const char *sY, VECTOR vX, int cX, const char *sX, VECTOR vZ, int cZ, const char *sZ, int iSymbol, int iIndex, int *piArgInt, int cArgInt, double *pdArgDbl, int cArgDbl);
typedef int	 OXCALL	FN_NewOxDrawWindow(const char *, const char *);
typedef void OXCALL	FN_NewOxMessage(char *);
typedef void OXCALL	FN_NewOxPuts(char *);
typedef void OXCALL	FN_NewOxRunMessage(char *);
typedef int	 OXCALL	FN_NewOxTextWindow(const char *, const char *);

/*========================= functions ======================================*/
/* exported functions */
BOOL	OXCALL	FOxCallBack(OxVALUE *pvFunc, OxVALUE *rtn, OxVALUE *pv, int cArg);
BOOL	OXCALL	FOxCallBackMember(OxVALUE *pvClass, const char *sMember, OxVALUE *rtn, OxVALUE *pv, int cArg);
BOOL	OXCALL	FOxCreateObject(const char *sClass, OxVALUE *rtn, OxVALUE *pv, int cArg);
BOOL	OXCALL	FOxGetDataMember(OxVALUE *pvClass, const char *sMember, OxVALUE *rtn);
BOOL	OXCALL	FOxLibAddFunction(char *sFunc, OxFUNCP pFunc, BOOL fVarArg);
BOOL 	OXCALL  FOxLibAddFunctionEx(char *sFunc, OxFUNCP pFunc, int cArgs, int flFlags);
BOOL	OXCALL	FOxRun(int iMainIP, char *sFunc);
BOOL	OXCALL	FOxSetDataMember(OxVALUE *pvClass, const char *sMember, OxVALUE *pv);
int 	OXCALL	IOxRunInit(void);
int 	OXCALL	IOxVersion(void);
int  	OXCALL  IOxVersionIsProfessional(void);
int 	OXCALL	IOxVersionOxo(void);
void 	OXCALL  OxCloneObject(OxVALUE *rtn, OxVALUE *pvObject, BOOL bDeep);
void	OXCALL	OxDeleteObject(OxVALUE *pvClass);
void	OXCALL	OxFnDouble(OxVALUE *rtn, OxVALUE *pv, double (OXCALL * fn)(double) );
void	OXCALL	OxFnDouble2(OxVALUE *rtn, OxVALUE *pv, double (OXCALL * fn)(double,double) );
void	OXCALL	OxFnDouble3(OxVALUE *rtn, OxVALUE *pv, double (OXCALL * fn)(double,double,double) );
void	OXCALL	OxFnDouble4(OxVALUE *rtn, OxVALUE *pv, double (OXCALL * fn)(double,double,double,double) );
void	OXCALL	OxFnDoubleInt(OxVALUE *rtn, OxVALUE *pv, double (OXCALL * fn)(double,int) );
void 	OXCALL 	OxFnDoubleVec(OxVALUE *rtn, OxVALUE *pv, double (OXCALL * fn)(double), int iVecAt);
void 	OXCALL 	OxFnDouble2Vec(OxVALUE *rtn, OxVALUE *pv, double (OXCALL * fn)(double,double), int iVecAt);
void 	OXCALL 	OxFnDouble3Vec(OxVALUE *rtn, OxVALUE *pv, double (OXCALL * fn)(double,double,double), int iVecAt);
void	OXCALL	OxFreeByValue(OxVALUE *pv);
void	OXCALL	OxGetMainArgs(int *pcArgc, char ***pasArgv);
void	OXCALL	OxGetOxArgs(int *pcArgc, char ***pasArgv);
BOOL	OXCALL	OxGetOxEditMode(void);
int 	OXCALL	OxGetPrintlevel(void);
int 	OXCALL 	OxGetUserExitCode(void);
void	OXCALL	OxLibArgError(int iArg);
void	OXCALL	OxLibArgTypeError(int iArg, int iExpected, int iFound);
void	OXCALL	OxLibCheckArrayMatrix(OxVALUE *pv, int iFirst, int iLast, MATRIX m);
void	OXCALL	OxLibCheckMatrixSize(OxVALUE *pv, int iFirst, int iLast, int r, int c);
void	OXCALL	OxLibCheckSquareMatrix(OxVALUE *pv, int iFirst, int iLast);
void	OXCALL	OxLibCheckType(int iType, OxVALUE *pv, int iFirst, int iLast);
void	OXCALL	OxLibValArrayCalloc(OxVALUE *pv, int c);
void	OXCALL	OxLibValMatDup(OxVALUE *pv, MATRIX mSrc, int r, int c);
void	OXCALL	OxLibValMatMalloc(OxVALUE *pv, int r, int c);
void	OXCALL	OxLibValStrMalloc(OxVALUE *pv, int c);
int 	OXCALL	OxMain(int argc, char *argv[]);
int 	OXCALL	OxMain_T(int argc, TCHAR *argv[]);
int 	OXCALL	OxMainCmd(const char *sCommand);
void	OXCALL	OxMainExit(void);
void	OXCALL	OxMainInit(void);
void	OXCALL	OxMakeByValue(OxVALUE *pv);
void	OXCALL	OxMessage(char *s);
void	OXCALL	OxNothing(int);
void	OXCALL	OxPuts(char *s);
void	OXCALL	OxRunAbort(int i);
void	OXCALL	OxRunError(int iErno, char *sToken);
void	OXCALL	OxRunErrorMessage(char *s);
void	OXCALL	OxRunExit(void);
void	OXCALL	OxRunMainExitCall(void (OXCALL * fn)(void));
void	OXCALL	OxRunMessage(char *s);
void	OXCALL	OxRunWarningMessage(char *sFunc, char *sMsg);
void	OXCALL	OxSetMainArgs(int argc, char *argv[]);
void	OXCALL	OxSetOxArgs(int cArgc, char **asArgv);
void	OXCALL	OxSetPrintlevel(int iSet);
void 	OXCALL 	OxSetUserExitCode(int iError);
OxVALUE *OXCALL OxStoreCreate(int c);
void 	OXCALL  OxStoreDelete(OxVALUE *pv, int c);
int 	OXCALL	OxValColumns(OxVALUE *pv);
OxVALUE OXCALL  OxValDuplicate(OxVALUE *pv);
void 	OXCALL  OxValDuplicate2(OxVALUE *pvDest, OxVALUE *pvSrc);
OxVALUE *OXCALL	OxValGetArray(OxVALUE *pv);
int		OXCALL	OxValGetArrayLen(OxVALUE *pv);
OxVALUE *OXCALL	OxValGetArrayVal(OxVALUE *pv, int i);
void * 	OXCALL  OxValGetBlob(OxVALUE *pv, int *pI1, int *pI2);
const char * OXCALL OxValGetClassName(OxVALUE *pv);
BOOL	OXCALL	OxValGetDouble(OxVALUE *pv, double *pdVal);
BOOL	OXCALL	OxValGetInt(OxVALUE *pv, int *piVal);
MATRIX	OXCALL	OxValGetMat(OxVALUE *pv);
int		OXCALL	OxValGetMatc(OxVALUE *pv);
int		OXCALL	OxValGetMatr(OxVALUE *pv);
int		OXCALL	OxValGetMatrc(OxVALUE *pv);
char *	OXCALL	OxValGetString(OxVALUE *pv);
BOOL	OXCALL	OxValGetStringCopy(OxVALUE *pv, char *s, int mxLen);
int		OXCALL	OxValGetStringLen(OxVALUE *pv);
BOOL 	OXCALL  OxValGetVecc(OxVALUE *pv, VECTOR vX);
BOOL 	OXCALL  OxValGetVecr(OxVALUE *pv, VECTOR vX);
OxVALUE *OXCALL	OxValGetVal(OxVALUE *pv, int i);
OxVALUE *OXCALL OxValGetStaticObject(OxVALUE *pv);
BOOL	OXCALL	OxValHasFlag(OxVALUE *pv, int iFlag);
BOOL	OXCALL	OxValHasType(OxVALUE *pv, int iType);
int 	OXCALL	OxValRows(OxVALUE *pv);
void 	OXCALL  OxValSetArray(OxVALUE *pv, int c);
void 	OXCALL  OxValSetBlob(OxVALUE *pv, int i1, int i2, void *p);
void	OXCALL	OxValSetDouble(OxVALUE *pv, double dVal);
void	OXCALL	OxValSetInt(OxVALUE *pv, int iVal);
void 	OXCALL  OxValSetMat(OxVALUE *pv, MATRIX mVal, int r, int c);
void 	OXCALL  OxValSetMatZero(OxVALUE *pv, int r, int c);
void	OXCALL	OxValSetNull(OxVALUE *pv);
void	OXCALL	OxValSetString(OxVALUE *pv, const char *sVal);
void 	OXCALL  OxValSetVecc(OxVALUE *pv, VECTOR vX, int r, int c);
void 	OXCALL  OxValSetVecr(OxVALUE *pv, VECTOR vX, int r, int c);
void	OXCALL	OxValSetZero(OxVALUE *pv);
int 	OXCALL	OxValSizec(OxVALUE *pv);
int 	OXCALL	OxValSizer(OxVALUE *pv);
int 	OXCALL	OxValSizerc(OxVALUE *pv);
OxVALUE OXCALL  OxValTransfer(OxVALUE *pv);
int 	OXCALL	OxValType(OxVALUE *pv);
char *	OXCALL	SOxGetTypeName(int iType);
char *	OXCALL	SOxIntFunc(void);
void	OXCALL	SetOxDataWindow(int (OXCALL * pfnNewOxDataWindow)(const char *, const char *) );
void	OXCALL	SetOxDraw(FN_NewOxDraw *pfnNewOxDraw);
void	OXCALL	SetOxDrawWindow(FN_NewOxDrawWindow *pfnNewOxDrawWindow);
void	OXCALL	SetOxExit(void (OXCALL * pfnNewOxExit)(int) );
void	OXCALL	SetOxGets(char * (OXCALL * pfnNewOxGets)(char *, int) );
void	OXCALL	SetOxMessage(FN_NewOxMessage *pfnNewOxMessage);
void	OXCALL	SetOxPipe(int cPipe);
void	OXCALL	SetOxPuts(FN_NewOxPuts *pfnNewOxPuts);
void	OXCALL	SetOxRunMessage(FN_NewOxRunMessage *pfnNewOxRunMessage);
void	OXCALL	SetOxTextWindow(FN_NewOxTextWindow *pfnNewOxTextWindow);

/*===================== typedef for kernel functions =======================*/
typedef int	OXCALL	FN_OxKnlCholeski(MATRIX mA, int cA);
typedef int	OXCALL	FN_OxKnlInvertsym(MATRIX mA, int cA, double *pdLogDet);
typedef int	OXCALL	FN_OxKnlInvert(MATRIX mA, int cA, double *pdLogDet, int *piSignDet);
typedef int	OXCALL	FN_OxKnlEigensym(MATRIX mA, int cA, VECTOR vEval, int iDoVectors);
typedef int	OXCALL	FN_OxKnlEigen(MATRIX mA, int cA, VECTOR vRe, VECTOR vIm, MATRIX mEvec);
typedef int	OXCALL	FN_OxKnlDecqr(MATRIX mXt, int cX, int cT, int *piPiv, VECTOR vTau);
typedef int	OXCALL	FN_OxKnlDecsvd(MATRIX mA, int cM, int cN, VECTOR vW, int fDoU, MATRIX mU,
    int fDoV, MATRIX mV, int fSort);
typedef void OXCALL FN_OxKnlMatmul(MATRIX mC, MATRIX mA, int rA, int cA, MATRIX mB,
	int rB, int cB,	int fTransA, int fTransB);


#ifdef __cplusplus
}
#endif

#endif  /* INC_OXEXPORT */


