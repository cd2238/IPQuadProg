SRCPATH   = ../src/
LIBPATH   = ../lib/
LIB       = optim
STATUT=debug
OPTS      =  -Wall -fPIC 
ifeq ($(STATUS), debug)
OPTS     += -g -DDEBUG
else
OPTS     += -O2
endif
QPLP = $(SRCPATH)solvesysqp.o $(SRCPATH)quadprosimp.o $(SRCPATH)initqp.o \
       $(SRCPATH)modchol.o  $(SRCPATH)modchol2.o   $(SRCPATH)calcalpha.o
$(QPLP):$(FRC)

      

FRC:
	@FRC=$(FRC)	

all: EXE

$(LIB): $(QPLP)
	ar cr $(LIBPATH)$@.a $(QPLP)
	ranlib $(LIBPATH)$@.a
	
	
EXE: $(LIB)
	gfortran $(OPTS) -o ../bin/quadpro  $(SRCPATH)main.F90 ../lib/optim.a  -llapack
	gfortran $(OPTS) -o ../bin/modchol  $(SRCPATH)main_modchol.F90 ../lib/optim.a  -llapack

clean:
	rm -rf  ../src/*.o
	
cleanall: clean 
	rm -rf  $(LIB).a $(LIB).so ../bin/quadpro ../bin/modchol 
	
%.o: %.F90
	gfortran $(OPTS) -c $< -o $@		
