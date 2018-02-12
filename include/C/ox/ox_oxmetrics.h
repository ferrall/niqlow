/*--------------------------------------------------------------------------
 * ox_oxmetrics.h - definitions/declarations for linking Ox to OxMetrics
 *
 *       (C) Jurgen Doornik 1995-2004
 *
 *--------------------------------------------------------------------------*/

#ifndef INC_OX_OXMETRICS
#define INC_OX_OXMETRICS

#ifdef __cplusplus
extern "C" {
#endif

BOOL OXCALL FOxOxMetricsStart(const TCHAR *OxModuleName, const TCHAR *OxWindowName, BOOL bUseStdHandles);
BOOL OXCALL FOxOxMetricsStartA(const char *OxModuleName, const char *OxWindowName, BOOL bUseStdHandles);
BOOL OXCALL FOxOxMetricsStartBatch(const TCHAR *OxModuleName, const TCHAR *OxWindowName, BOOL bUseStdHandles, int iBatch);
BOOL OXCALL FOxOxMetricsStartBatchA(const char *OxModuleName, const char *OxWindowName, BOOL bUseStdHandles, int iBatch);
void OXCALL OxOxMetricsFinish(BOOL bFocusText);
void OXCALL OxOxMetricsRestart();
void OXCALL OxOxMetricsSetRunMessageLevel(int iLevel);
const TCHAR * OXCALL OxOxMetricsGetRunMessage();
BOOL OXCALL FOxOxMetricsStartChecked(const TCHAR *OxModuleName, const TCHAR *OxWindowName, BOOL bUseStdHandles, int iBatch,
	int iVersion100);
BOOL OXCALL FOxOxMetricsStartCheckedA(const char *OxModuleName, const char *OxWindowName, BOOL bUseStdHandles, int iBatch,
	int iVersion100);
void OXCALL OxOxMetricsStopOutput();
void * OXCALL OxMetricsGetOxHandlerA(const char *sName);
BOOL OXCALL OxOxMetricsUsingStdHandles();

#ifdef __cplusplus
}
#endif

#endif  /* INC_OX_OXMETRICS */
