#include <oxstd.h>
class CGI {
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
    static decl http, kvals, query;
    static Initialize(title="OX CGI");
    static ParseQ();
    static GetVar(key);
    static Finalize();
    }
