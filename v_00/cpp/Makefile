# code version:
VER = 00


CXXFLAGS := -I.$(pwd)/include

CXXFLAGS += -O2 -std=c++11

#CXXFLAGS += -Wall -Wclobbered -Wempty-body -Wignored-qualifiers -Wmissing-field-initializers -Wsign-compare -Wtype-limits -Wno-unused -Wno-uninitialized

#CXXFLAGS += -mtune=generic -msse2 -mfpmath=sse


# Archiving options:
AR = ar -crs

VPATH = . src build bin lib

code_src = main

code_lib =

code_exe = bg.$(VER).x


# targets:
all:            libs exes

libs:           $(code_lib)

exes:           $(code_exe)

$(code_exe):    $(addsuffix .o, $(code_src))
		$(CXX) -o $(@F) $(CXXFLAGS) $(LDFLAGS) $(addprefix build/, $(^F)) $(LIBS)
		mv -f $(@F) bin/

%.o:            %.cpp
%.o:            %.cpp Makefile
		$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) $<
		mv -f $(@F) build/

.PHONY:         clean

clean:
		rm -f *.o build/* lib/*

#add automatic generation of the source code dependencies:
build/%.d: %.cpp
	@set -e; rm -f $@; \
	$(CXX) $(CXXFLAGS) -MM $< > $@.$$$$; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$

include $(addprefix build/, $(addsuffix .d, $(code_src)))
