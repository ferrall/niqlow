#ifndef curlHEADED
#define curHEADED
#import "Shared"
	extern "CFcurl,fget"   			get(url,file);
#endif

enum{FINDCMD=1,INSTALL=2}
class net {
    static decl iswin;
    }

class WGET : net {
    static decl cmd, ready, base;
    static init(incmd=0);
    static Base(burl);
    static Run(urls);
    }
