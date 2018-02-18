#include "net.h"
#ifdef OX_Windows
    const decl _platform = "Win";
#else
    const decl _platform = "Unix";
#endif

WGET::Base(burl) {
    if (!ready) {
        oxwarning("Can't set base URL before WGET is ready");
        return;
        }
    base = burl;
    }
WGET::Run(urls) {
    if (!ready) init();
    decl retval = systemcall(cmd+" -q -o _wglogXXX "+" "+urls);
    if (retval) oxwarning("If wget is not recognized, WGET::init(FINDCMD) or WGET::init(INSTALL)");
    return retval;
    }

WGET::init(findcmd) {
    if (ready) {
        if (ready>1) oxwarning("net commands already initialized");
        ready = 2;
        return;
        }
    ready = 1;
    base = "\"\"";
    iswin = _platform=="Win";
    cmd = "wget";
    decl notfnd;
    if (findcmd) {
        notfnd = (iswin)
                    ? systemcall("where /R \"%ProgramFiles(x86)%\" wget.exe > _wgcmdXXX")
                    :systemcall("which wget > _wgcmdXXX");
        if (notfnd) oxrunerror("Could find wget command");
        decl wf = fopen("_wgcmdXXX");
        fscan(wf,"%z",&cmd);
        cmd = "\""+replace(cmd,"\\","/","a")+"\"";
        cmd = cmd;
        fclose(wf);
        }
    }
