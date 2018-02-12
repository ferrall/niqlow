#include <oxstd.h>

extern "threes,FnThrees" Threes(const r, const c);
extern "threes,FnGet"    Get(url,file);

main()
{
    print(Threes(3,3));
//    Get("http://econ.queensu.ca/index.html","test.html");
}
