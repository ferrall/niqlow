using System;
using System.Runtime.InteropServices;

static class OxMetrics
{
// Ox-OxMetrics run-time functionality
	[DllImport("OxMetrics7", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern int FOxOxMetricsStartA(string OxModuleName, string OxWindowName, int bUseStdHandles);
	[DllImport("OxMetrics7", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern int FOxOxMetricsStartBatchA(string OxModuleName, string OxWindowName, int bUseStdHandles, int iBatch);
	[DllImport("OxMetrics7", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern void OxOxMetricsFinish(int bFocusText);
	[DllImport("OxMetrics7", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern void OxOxMetricsRestart();
	[DllImport("OxMetrics7", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern IntPtr OxMetricsGetOxHandlerA(string sName);
	[DllImport("OxMetrics7", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern int OxOxMetricsUsingStdHandles();

// OxMetrics run-time functionality
	[DllImport("OxMetrics7", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern int OxMetricsStartA(string sClient, int iBatchIndex);
	[DllImport("OxMetrics7", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern int OxMetricsStartExA(string sClient, int iBatchIndex, int iOxMetricsVersion);
	[DllImport("OxMetrics7", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern int OxMetricsStartAdvise();
	[DllImport("OxMetrics7", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern int OxMetricsStopAdvise();
	[DllImport("OxMetrics7", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern void OxMetricsInit();
	[DllImport("OxMetrics7", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern void OxMetricsFinish();
	[DllImport("OxMetrics7", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern int OxMetricsDLLVersion();
	[DllImport("OxMetrics7", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern int ClientIsRegisteredA(string sClient, int iVersion100);
	[DllImport("OxMetrics7", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern void SetBatchClient();
	[DllImport("OxMetrics7", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern void SetOxMetricsBlock(int iBlock);
	[DllImport("OxMetrics7", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern void SetOxMetricsBlockOnShowDialog(int iShow);
	[DllImport("OxMetrics7", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern void SendOxMetricsMsgA(string sMessage);

// OxMetrics text functions
	[DllImport("OxMetrics7", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern int FindFirstTextA(string sText, int iOption);
	[DllImport("OxMetrics7", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern void FocusTextWindow();
	[DllImport("OxMetrics7", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern void SetTextMarker(int iMode);
	[DllImport("OxMetrics7", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern void SetTextShowMode(int iMode);
	[DllImport("OxMetrics7", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern void SetTextWindowA(string sTitle);
	[DllImport("OxMetrics7", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern void ShowTextWindow();
	[DllImport("OxMetrics7", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern void PrintTextA(string sResults);
	[DllImport("OxMetrics7", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern void SpoolText();
	[DllImport("OxMetrics7", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern void BatchDoneA(string sResults);
	[DllImport("OxMetrics7", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern void SetBatchCodeA(string sResults);
	[DllImport("OxMetrics7", CharSet = CharSet.Ansi, ExactSpelling = true)]
		public static extern void SetOxCodeA(string sResults);
}