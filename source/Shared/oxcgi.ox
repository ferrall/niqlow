#import "CGI"

main() {
    CGI::Initialize();
    println(CGI::GetVar("Query_String"));
    CGI::Finalize();
    }
