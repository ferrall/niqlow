/*--------------------------------------------------------------------------
 * ox_oxmetrics_bridge.cpp - definitions/declarations for bridging Ox - OxMetrics
 *
 *       (C) Jurgen Doornik 2012
 *
 *--------------------------------------------------------------------------*/

#include "oxexport.h"
#include "ox_oxmetrics.h"
#include "ox_oxmetrics_bridge.h"


/*-------------------------- OxOxMetricsSetHandlers ------------------------*/
extern "C"
void OXCALL OxOxMetricsSetHandlers(BOOL bRestart)
{
	if (bRestart)
	{
		OxOxMetricsRestart();
		OxMainExit();
	}
	OxMainInit();			// call first, then replace io and drawing functions

	if (!OxOxMetricsUsingStdHandles())
	{
		SetOxPipe(0);
		SetOxMessage((FN_NewOxMessage *)OxMetricsGetOxHandlerA("OxMessage"));
		SetOxRunMessage((FN_NewOxRunMessage *)OxMetricsGetOxHandlerA("OxRunMessage"));
		SetOxPuts((FN_NewOxPuts *)OxMetricsGetOxHandlerA("OxPuts"));
		SetOxTextWindow((FN_NewOxTextWindow *)OxMetricsGetOxHandlerA("OxTextWindow"));
	}
	SetOxDrawWindow((FN_NewOxDrawWindow *)OxMetricsGetOxHandlerA("OxDrawWindow"));
	SetOxDraw((FN_NewOxDraw *)OxMetricsGetOxHandlerA("OxDraw"));
}
/*------------------------ END OxOxMetricsSetHandlers ----------------------*/
