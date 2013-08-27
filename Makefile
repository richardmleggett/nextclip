BIN = bin
MAXK = 63

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

ifdef MAC
MACFLAG = -fnested-functions -L/opt/local/lib/ 
endif

OPT	= $(ARCH) -Wall -O3 $(MACFLAG) -DNUMBER_OF_BITFIELDS_IN_BINARY_KMER=$(BITFIELDS) -pthread -g

CFLAGS_NEXTCLIP = -Iinclude

NEXTCLIP_OBJ = obj/nextclip.o obj/hash_table.o obj/hash_value.o obj/logger.o obj/binary_kmer.o obj/element.o

all:remove_objects $(NEXTCLIP_OBJ)
	mkdir -p $(BIN); $(CC) -lm $(OPT) -o $(BIN)/nextclip $(NEXTCLIP_OBJ)

clean:
	rm obj/*
	rm -rf $(BIN)/nextclip

remove_objects:
	rm -rf obj

obj/%.o : src/%.c
	mkdir -p obj; $(CC) $(CFLAGS_NEXTCLIP) $(OPT) -c $< -o $@

