include Makefile.inc

INC := ..\include
DOC := ../docs
OxMenu = ..\..\OxMenu\include

vpath %.ox .:Shared:DDP:FiveO:$(INC)
vpath %.h $(INC)
vpath %.oxo $(INC)
vpath %.ox.html $(DOC)

#Shared.oxo 
oxobjects = DDP.oxo FiveO.oxo CFMPI.ox

niqlow.oxo  : $(oxobjects)
#	-for d in $(DIRS); do ($(CD) $$d; $(MAKE) $(MKFLAGS); $(ECHO) "\n\n\n Build attempt for directory: $$d complete\n\n\n"); done

%.oxo : %.ox %.h
	$(OX) $(OXFLAGS) -i$(INC) -d -c $<
	$(COPY) $@ $(INC)
	$(ERASE) $@

.PHONY : document
document:
	$(ERASE) ..\docs\default.html
	$(OXDOC) -include $(INC) niqlow.ox
	${MAKE} -C About document
	${MAKE} -C DDP document
	PAUSE "DDP"
	${MAKE} -C FiveO document
	PAUSE "FiveO"
	${MAKE} -C CFMPI document
	${MAKE} -C $(XMP) document
	${MAKE} -C $(DOC) tweak
	$(COPY) ..\docs\index.tmp ..\docs\default.html
	$(COPY) ..\docs\index.tmp ..\docs\index.html
#	$(COPY) ..\docs\DDP\default.html ..\docs\DDP\list.html
#	$(COPY) ..\docs\DDP\default.tmp ..\docs\DDP\default.html
#	$(COPY) ..\docs\FIVEO\default.html ..\docs\FIVEO\list.html
#	$(COPY) ..\docs\FIVEO\default.tmp ..\docs\FIVEO\default.html

.PHONY : install
install:
	$(ERASE) $(INC)\*.oxo
	$(COPY) $(OxMenu)\menu.* $(INC)
	${MAKE} -C Shared
	${MAKE} -C DDP
	${MAKE} -C FiveO
	$(COPY) FiveO\ParallelObjective.ox $(INC)
	$(ERASE) $(INC)\ParallelObjective.oxo
#	${MAKE} -C $(INC)
	${MAKE} -C .
	${MAKE} -C $(XMP)
	$(ERASE) /S $(XMP)\*.oxo
