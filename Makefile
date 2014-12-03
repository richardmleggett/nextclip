BIN = bin
MAXK = 63

ifdef MAC
define n


endef
$(warning $n$nNOTE On Mac OS X, you must compile using GCC. If you have not already done so, please change the line CC= in this Makefile to point to your GCC compiler. $nThis must be a full path, as Apple aliases gcc to point to it's own compiler.$n)

# Change the following to point to your GCC binary
CC=gcc

# On my MAC, I use CC=~/gcc/bin/gcc

# On older versions of XCode, it was necessary to include the following
#MACFLAG = -fnested-functions -L/opt/local/lib/ 
endif

ifeq ($(MAXK),31)
   BITFIELDS = 1
endif

ifeq ($(MAXK),63)
   BITFIELDS = 2
endif

ifeq ($(MAXK),95)
   BITFIELDS = 3
endif

ifeq ($(MAXK),127)
   BITFIELDS = 4
endif

OPT	= $(ARCH) $(MACFLAG) -Wall -O3 -DNUMBER_OF_BITFIELDS_IN_BINARY_KMER=$(BITFIELDS) -pthread -g

CFLAGS_NEXTCLIP = -Iinclude

NEXTCLIP_OBJ = obj/nextclip.o obj/hash_table.o obj/hash_value.o obj/logger.o obj/binary_kmer.o obj/element.o

all:remove_objects $(NEXTCLIP_OBJ)
	mkdir -p $(BIN); $(CC) $(OPT) -o $(BIN)/nextclip $(NEXTCLIP_OBJ) -lm -lz

clean:
	rm obj/*
	rm -rf $(BIN)/nextclip

remove_objects:
	rm -rf obj

obj/%.o : src/%.c
	mkdir -p obj; $(CC) $(CFLAGS_NEXTCLIP) $(OPT) -c $< -o $@

