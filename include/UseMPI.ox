/** Include this file to use CFMPI and either fake or real message passing.

If you wish to use either real or fake message passing in the execution of your program then include
this file in your program.
<DD><pre>
#include "UseMPI.ox"
</pre></dd>

**/
#include "CFMPI.ox"
#ifdef CFMPIfakeDEFINED
#include "MPIinterface.ox"
#endif
