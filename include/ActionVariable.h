#import "Shared"

/** An element of the action vector &alpha;.
@example
<pre>MyModel : DPparent {
    static decl a;
    &hellip;
    }
&hellip;
MyModel::Initialize() {
    DPparent::Initialize();
    &hellip;
    a = new ActionVariable("choice",2);
    Actions(a);
    &hellip;
    CreateSpaces();
    }</pre>
</dd>
**/
struct ActionVariable : Discrete	{
	ActionVariable(L="a",N=1);
	}

/** Easy way to create a binary choice.
@example
<pre>MyModel : DPparent {
    static decl a;
    &hellip;
    }
&hellip;
MyModel::Initialize() {
    DPparent::Initialize();
    &hellip;
    a = new BinaryChoice();
    Actions(a);
    &hellip;
    CreateSpaces();
    }</pre>
</dd>
**/
struct BinaryChoice : ActionVariable {
    BinaryChoice();
    }
