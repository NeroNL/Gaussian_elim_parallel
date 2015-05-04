include $(PUB)/Arch/arch.gnu-c++11.generic
WARNINGS += -Wall -pedantic

WARNINGS += -Wall -pedantic


# This generates output about how the
# compiler vectorized the code
# We  suggest using level 2 (the integer after "verbose=")
# See the gcc manual for the other levels of output: levels 1-7
# http://gcc.gnu.org/onlinedocs/gcc-4.7.3/gcc/Debugging-Options.html#Debugging-Options
#OPTIMIZATION += -ftree-vectorizer-verbose=2
# OPTIMIZATION += -ftree-vectorize


# If you want to specify the logarithmic barrier
# specify log=1 on the "make" command line
ifeq ($(log), 1)
        C++FLAGS += -DLOG_BARRIER
endif   


app:	lu


OBJ	= lu.o elim.o verify.o init.o cmdLine.o Report.o Timer.o

lu:	$(OBJ)
	$(C++LINK) $(LDFLAGS) -o $@ $(OBJ)  $(LDLIBS)

.PHONY: clean
clean:	
	$(RM) *.o lu
	$(RM) core
