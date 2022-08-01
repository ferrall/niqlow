#include "oxstd.h"

/** A menu object.
**/
struct Menu {
	enum{prompt,call,ItemSize}
	enum{INVALID=-4,EXIT,HELP,QUIT,DOALL,CHOICE}
	static const decl sep  = "\n---------------------------\n";
	static decl
        /** Folder for logs.**/                    logdir = "./";
	decl
    /** log calls from this menu.**/        keeplog,
	/** string: menu name.**/ 				name,
	/** array containing menu.**/			items;
    static SetLogDir(folder="./");
    static itemparse(token,items);
	Menu(name="Menu",keeplog=FALSE);
    ~Menu();
	virtual add(...);
    virtual CmdLine();
	}

/** A class for choosing from a menu of functions and submenus interactively or from the command line. **/
struct CallMenu : Menu {
    decl
	/** TRUE if top menu.**/				IamMain,
	/** string printed if HELP selected.**/ help_text;
	CallMenu(name="Menu",keeplog=FALSE,IamMain=TRUE);
	make_the_call(item);
	add(...);
    CmdLine(alist=1);  //program name is 0 element of arglist
	Run();
	}

/** A class for reading in parameters interactively.**/
struct ParamMenu : Menu {
    enum{RESET=-2,SEND,Value=1}
    decl pars,TargFunc;
    add(...);
    ParamMenu(name="Set Params",keeplog=FALSE,TargFunc=0);
    SetPars(TargFunc=0);
    CmdLine(TargFunc=0,args=1);
    }
