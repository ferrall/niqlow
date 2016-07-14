#include "CGI.oxh"
CGI::Initialize(title) {
    kvals = {};
    decl key;																						
    foreach (key in keys) kvals |= getenv(key);
    decl aa = arglist();
    post = fopen(aa[2],"r");
    out = fopen(aa[1],"a");
    fprintln(out,"<!doctype html><html xml:lang=\"en\">",
        "<head><meta charset=\"utf-8\"><title>",title,
        "</title><meta name=\"description\" content=\"CGI for Ox\">",
        "<meta name=\"author\" content=\"Christopher Ferrall\">",
        "<script type=\"text/x-mathjax-config\">",
        "MathJax.Hub.Config({tex2jax: {inlineMath:",
        "[[\"$\",\"$\"],[\"\\(\",\"\\)\"]]}});</script>",
        "<script type=\"text/javascript\"",
        "src=\"http://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML\"></script>",
        "</head><body>");
    }

CGI::Finalize() {
    fprintln(out,"<h2>Ox Output</h2><code><pre>");
    fclose(out);
    fclose(post);
    }

CGI::ParseQ() {
    decl q = GetVar("QUERY_STRING");
    }

CGI::GetVar(key) {
    if (!isstring(key)) {
        oxwarning("key must be a string");
        return -1;
        }
    decl ind = strfind(keys,key);
    if (ind==-1) return -1;
    return kvals[ind];
    }

CGI::Parse() {
    if (!isfile(post)) {
        oxwarning("post data file not found.");
        return {};
        }
	decl instr, loc,nms,vals, val;
    nms = {}; vals = <>;
    fscan(post,"%s",&instr);
    do {
        loc=strfind(instr,eq);
        if (loc>0) {
            nms |= instr[:loc-1];
            instr = instr[loc:];
            loc = strfind(instr,amp);
            if ((loc>0)) {
                sscan(instr[1:loc-1],"%g",&val);
                instr = instr[loc+1:];
                }
            else {
                sscan(instr[1:],"%g",&val);
                instr = "";
                }
            vals |= val;
            }
        }  while (loc>0);
    return {nms,vals};
	}

CGI::VolumeCtrl(fp,pref,Volume) {
    fprint  (fp,"<input type=\"radio\" name=\"",pref,"Volume\""," value=\"",-1,"\" ",Volume==-1 ? "checked>" : ">","SILENT&nbsp;");
    fprint  (fp,"<input type=\"radio\" name=\"",pref,"Volume\""," value=\"", 0,"\" ",Volume== 0 ? "checked>" : ">","QUIET&nbsp;");
    fprint  (fp,"<input type=\"radio\" name=\"",pref,"Volume\""," value=\"", 1,"\" ",Volume== 1 ? "checked>" : ">","LOUD&nbsp;");
    fprintln(fp,"<input type=\"radio\" name=\"",pref,"Volume\""," value=\"", 2,"\" ",Volume== 2 ? "checked>" : ">","NOISY; &emsp;");
    }
