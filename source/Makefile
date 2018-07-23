include Makefile.inc

INC := ..\include
DOC := ../docs

vpath %.ox .:Shared:DDP:FiveO:Hybrids:$(INC)
vpath %.h $(INC)
vpath %.oxo $(INC)
vpath %.ox.html $(DOC)

oxobjects = Shared.oxo DDP.oxo FiveO.oxo CFMPI.ox Hybrids.oxo

niqlow.oxo  : $(oxobjects)
#	-for d in $(DIRS); do ($(CD) $$d; $(MAKE) $(MKFLAGS); $(ECHO) "\n\n\n Build attempt for directory: $$d complete\n\n\n"); done

%.oxo : %.ox %.h
	$(OX) $(OXFLAGS) -i$(INC) -d -c $<
	$(COPY) $@ $(INC)
	$(ERASE) $@

.PHONY : clean
clean:
	$(info $(CURDIR))
	$(ERASE) *.oxo
	$(ERASE) $(INC)\*.oxo

.PHONY : document
document:
#	$(ERASE) ..\doc\default.html
#	$(OXDOC) -uplevel -include $(INC) niqlow.ox "About niqlow"
	${MAKE} -C $(DOC) tweak
	${MAKE} -C About document
	${MAKE} -C DDP document
	${MAKE} -C FiveO document
	${MAKE} -C Hybrids document
	${MAKE} -C CFMPI document
	${MAKE} -C $(XMP) document
	
.PHONY : install
install:
	$(ERASE) $(INC)\*.oxo
	${MAKE} -C DDP
	${MAKE} -C FiveO
	$(COPY) FiveO\ParallelObjective.ox $(INC)
	$(ERASE) $(INC)\ParallelObjective.oxo
	${MAKE} -C Hybrids
	${MAKE} -C $(INC)
	${MAKE} -C .
	$(ERASE) $(XMP)\*.oxo
	${MAKE} -C $(XMP)
	$(ERASE) $(XMP)\PeerCommunicationTestA.oxo
	$(ERASE) $(XMP)\ClientServerTestA.oxo
#   $(ERASE) $(INC)\C\*.so
#	${MAKE} -C C
