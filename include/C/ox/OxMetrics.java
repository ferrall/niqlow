package oxmetrics;

import com.sun.jna.Library;
import com.sun.jna.Native;
import com.sun.jna.Pointer;
import com.sun.jna.PointerType;
import com.sun.jna.ptr.DoubleByReference;
import com.sun.jna.ptr.IntByReference;
import com.sun.jna.ptr.PointerByReference;
import com.sun.jna.win32.StdCallLibrary;
import com.sun.jna.Platform;


public class OxMetrics implements StdCallLibrary {
	static {
		Native.register(Platform.isWindows() ? "oxmetrics7" : Platform.isMac() ? "liboxmetrics.7.dylib" : "liboxmetrics.so.7");
	}
// Ox-OxMetrics run-time functionality
	public static native int FOxOxMetricsStartA(String OxModuleName, String OxWindowName, int bUseStdHandles);
	public static native int FOxOxMetricsStartBatchA(String OxModuleName, String OxWindowName, int bUseStdHandles, int iBatch);
	public static native void OxOxMetricsFinish(int bFocusText);
	public static native void OxOxMetricsRestart();
	public static native Pointer OxMetricsGetOxHandlerA(String sName);
	public static native int  OxOxMetricsUsingStdHandles();

// OxMetrics run-time functionality
	public static native int OxMetricsStartA(String sClient, int iBatchIndex);
	public static native int OxMetricsStartExA(String sClient, int iBatchIndex, int iOxMetricsVersion);
	public static native int OxMetricsStartAdvise();
	public static native int OxMetricsStopAdvise();
	public static native void OxMetricsInit();
	public static native void OxMetricsFinish();
	public static native int OxMetricsDLLVersion();
	public static native int ClientIsRegisteredA(String sClient, int iVersion100);
	public static native void SetBatchClient();
	public static native void SetOxMetricsBlock(int iBlock);
	public static native void SetOxMetricsBlockOnShowDialog(int iShow);
	public static native void SendOxMetricsMsgA(String sMessage);

// OxMetrics text functions
	public static native int FindFirstTextA(String sText, int iOption);
	public static native void FocusTextWindow();
	public static native void SetTextMarker(int iMode);
	public static native void SetTextShowMode(int iMode);
	public static native void SetTextWindowA(String sTitle);
	public static native void ShowTextWindow();
	public static native void PrintTextA(String sResults);
	public static native void SpoolText();
	public static native void BatchDoneA(String sResults);
	public static native void SetBatchCodeA(String sResults);
	public static native void SetOxCodeA(String sResults);
}