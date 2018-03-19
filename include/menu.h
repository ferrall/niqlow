#include "oxstd.h"

/** A menu object.
**/
struct Menu {
	enum{prompt,call,ItemSize}
	enum{EXIT = -3,HELP,QUIT,DOALL,CHOICE}
	static const decl sep  = "\n---------------------------\n";
	static decl
        /** Folder for logs.**/                    logdir = "./",
	   /** TRUE: create a logfile for this run.**/ logoutput = FALSE;
	decl
	/** string: menu name.**/ 				name,
	/** array containing menu.**/			items,
	/** TRUE if top menu.**/				IamMain,
	/** string printed if HELP selected.**/ help_text;
	Menu(name="Menu",IamMain=TRUE);
	~Menu();
	make_the_call(item);
	add(...);
	Run();
    CmdLine();
	}
