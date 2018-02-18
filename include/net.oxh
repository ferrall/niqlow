#include <oxstd.h>
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

class CGI : net {
    static const decl keys ={
        "AUTH_TYPE",
        "CONTENT_LENGTH",
        "CONTENT_TYPE",
        "GATEWAY_INTERFACE",
        "PATH_INFO",
        "PATH_TRANSLATED",
        "QUERY_STRING",
        "REMOTE_HOST",
        "REMOTE_IDENT",
        "REMOTE_USER",
        "REQUEST_METHOD",
        "SCRIPT_NAME",
        "SERVER_NAME",
        "SERVER_PORT",
        "SERVER_PROTOCOL",
        "SERVER_SOFTWARE"};
    static const decl eq='=', amp='&', dnvsuffix = "-dnv", ivalsuffix="-val", cgiopt = "-DCGI",
            pfopt = "post=", htopt="html=";
    static decl iscgi, post, out, kvals, query;
    static Initialize(title="OX CGI");
    static ParseQ();
    static Parse();
    static GetVar(key);
    static Finalize();
    static VolumeCtrl(pref="",Volume=0);
    static CheckBox(nm,val,checked);
    static ReadForm(list);
    static CreateForm(list);
    }
