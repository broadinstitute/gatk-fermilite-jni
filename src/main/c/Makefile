ifeq ($(JAVA_HOME),)
$(error $$(JAVA_HOME) is not defined)
endif

JNI_INCLUDE_DIRS=$(addprefix -I,$(shell find $(JAVA_HOME)/include -type d))
CFLAGS=-O2 -Wall -fPIC $(JNI_INCLUDE_DIRS)
CC=gcc

#OS-dependent extension lookup
UNAME := $(shell uname)
ifeq ($(UNAME),Darwin)
LIB_EXT=Darwin.dylib
else
LIB_EXT=Linux.so
endif


all: libfml.$(LIB_EXT)

libfml.$(LIB_EXT): org_broadinstitute_hellbender_utils_fermi_FermiLiteAssembler.o fermi-lite/libfml.a
	$(CC) -dynamiclib -shared -o $@ $^ -lm -lz -lpthread

fermi-lite/libfml.a:
	git clone https://github.com/lh3/fermi-lite
	$(MAKE) CFLAGS="$(CFLAGS)" -C fermi-lite libfml.a

JNI_BASE_NAME=org_broadinstitute_hellbender_utils_fermi_FermiLiteAssembler
$(JNI_BASE_NAME).o: $(JNI_BASE_NAME).c fermi-lite/libfml.a

test: test.o fermi-lite/libfml.a
	$(CC) -o $@ $^ -lm -lz -lpthread

clean:
	rm -rf fermi-lite *.o *.$(LIB_EXT)

.PHONY: all clean
