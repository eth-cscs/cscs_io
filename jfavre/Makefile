CC     = cc
CFLAGS = -w
H := $(shell echo $(HOST) | tr [a-z] [A-Z] | tr -d [0-9] )

ifeq ($(PE_ENV),GNU)
        CFLAGS +=
endif

ifeq ($(PE_ENV),PGI)
        CFLAGS +=
endif

ifeq ($(PE_ENV),CRAY)
        CFLAGS +=
endif

ifeq ($(PE_ENV),INTEL)
        CFLAGS += 
endif

LDFLAGS = $(CFLAGS)
OBJ     = benchmark.o
LIB     = -lz
.SUFFIXES: .o .c

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $(PE_ENV)_$@

$(PE_ENV).$(H): $(OBJ)
	$(CC) $(PE_ENV)_$< $(LIB) $(LDFLAGS) -o $@
	@echo $(PE_ENV) / $(H) / c executable ready

clean:
	rm -f *.o o_* runme.slurm.* *.$(H) 
