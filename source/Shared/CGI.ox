#include "CGI.oxh"
CGI::Initialize(title) {
    kvals = {};
    decl key;
    foreach (key in keys) kvals |= getenv(key);
    http = fopen(arglist()[0],"a");
    fprintln(http,"<!doctype html><html xml:lang=\"en\">",
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
    fprintln(http,"<h2>Ox Output</h2><code><pre>");
    fclose(http);
    println("</code></pre></body></html>");
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
