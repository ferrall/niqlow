
tool=&Ox
command=C:\Program Files\OxMetrics9\ox\oxl.exe
arguments=-iC:\Users\Chris\Documents\OFFICE\software\Microeconometrics\niqlow\include "$(FilePath)"
folder=$(FileDir)
defaultext=ox
capture=1
output=Ox Output
help=C:\Program Files\OxMetrics9\ox\doc\oxind*.html
clear=0
beep=0
warn=1
header=1
send=1
activate=1

tool=Ox&Run
command=C:\Program Files\OxMetrics9\apps\OxRun.exe
arguments="$(FilePath)"
folder=$(FileDir)
defaultext=ox
capture=0
output=
help=C:\Program Files\OxMetrics9\..\doc\oxind*.html
clear=0
beep=0
warn=1
header=1
send=0
activate=0

tool=Ox - &serial
command=C:\Program Files\OxMetrics9\ox\oxl.exe
arguments="-rp1 $(FilePath)"
folder=$(FileDir)
defaultext=ox
capture=1
output=Ox Output
help=C:\Program Files\OxMetrics9\ox\doc\oxind*.html
clear=0
beep=0
warn=1
header=1
send=0
activate=1

tool=Ox - &compile
command=C:\Program Files\OxMetrics9\ox\oxl.exe
arguments=-iC:\Users\Chris\Documents\OFFICE\software\Microeconometrics\niqlow\include -c "$(FilePath)"
folder=$(FileDir)
defaultext=ox
capture=1
output=Ox Output
help=C:\Program Files\OxMetrics9\ox\doc\oxind*.html
clear=0
beep=0
warn=1
header=1
send=0
activate=1

tool=Ox - &interactive
command=C:\Program Files\OxMetrics9\ox\oxl.exe
arguments=-q
folder=
defaultext=ox
capture=1
output=Session.ox
help=C:\Program Files\OxMetrics9\ox\doc\oxind*.html
clear=0
beep=0
warn=1
header=1
send=1
activate=1

tool=Ox - &debug
command=C:\Program Files\OxMetrics9\ox\oxl.exe
arguments=-d "$(FilePath)"
folder=$(FileDir)
defaultext=ox
capture=1
output=Debug.ox
help=C:\Program Files\OxMetrics9\ox\doc\oxind*.html
clear=0
beep=0
warn=1
header=1
send=1
activate=1
