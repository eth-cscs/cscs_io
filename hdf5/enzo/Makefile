CC       = CC
FLAGS    = -O3 # -g
DDTFLAGS = 
H := $(shell echo $(HOST) | tr [a-z] [A-Z] | tr -d [0-9] )

# Are you uts ?
ifneq ($(CRAY_HDF5_DIR),/opt/cray/hdf5-parallel/1.8.12)
  FLAGS += YOU MUST LOAD cray-hdf5-parallel/1.8.12
endif

# ifeq ($(PE_ENV),)
# 	echo "Setting PE_ENV=GNU"
#	PE_ENV=GNU
# endif

LDFLAGS        = $(FLAGS)
OBJ            = ReadInitialNew.o
OBJ2 	       = $(subst ReadInitialNew,$(PE_ENV)_ReadInitialNew,$(OBJ))
LIB            = 
.SUFFIXES: .o .C

%.o: %.C
	$(CC) $(FLAGS) -c $< -o $(PE_ENV)_$@ $(DDTFLAGS) 

$(PE_ENV).$(H): $(OBJ)
	$(CC) $(LDFLAGS) $(OBJ2) $(LIB) -o $@ $(DDTFLAGS)
	@echo $(PE_ENV) / $(H) / c executable ready

clean:
	rm -f *.o o_* runme.slurm.* *.$(H) *.mod cuda_profile* nvprof.output* *+pat
	rm -fr *.ptx *.opt *.cg *.ap2 *.rpt *.cub *.lst *.xf
	rm -fr scorep-* *.z *.thumb *.otf

