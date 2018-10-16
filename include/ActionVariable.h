#import "DDPShared"

/** An element of the action vector &alpha;.
@example
<pre>MyModel : DPparent {
    static decl a;
    &vellip;
    }
&vellip;
MyModel::Initialize() {
    DPparent::Initialize(new MyModel());
    &vellip;
    a = new ActionVariable("choice",2);
    Actions(a);
    &vellip;
    CreateSpaces();
    }</pre>
</dd>
**/
struct ActionVariable : Discrete	{
    decl vL;
	ActionVariable(L="a",N=1);
    virtual myAV();
    virtual myCV();
	}

/** Easy way to create a binary choice.
@example
<pre>MyModel : DPparent {
    static decl a;
    &vellip;
    }
&hellip;
MyModel::Initialize() {
    DPparent::Initialize(new MyModel());
    &vellip;
    a = new BinaryChoice();
    Actions(a);
    &vellip;
    CreateSpaces();
    }</pre>
</dd>
**/
struct BinaryChoice : ActionVariable {
    BinaryChoice(L="a");
    }
