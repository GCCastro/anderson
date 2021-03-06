CC = g++
SPECTRA = -I ~/SOFT/spectra/include
CFLAGS  = -c -Wall -g
RTNORM  = -I ./src
LIB     = -lm -lgsl -lgslcblas

SRC1 = src/rtnorm.cpp 
HDRS = src/rtnorm.hpp

All:laplace 

laplace: laplace_solver.o rtnorm.o
	$(CC) -o a $^ $(SPECTRA) $(RTNORM) $(LIB)

esp: laplace_solver_esp.o rtnorm.o
	$(CC) -o esp $^ $(SPECTRA) $(RTNORM) $(LIB)

1d: loc_1d.o rtnorm.o
	$(CC) -o 1d $^ $(SPECTRA) $(RTNORM) $(LIB)

laplace_solver.o: laplace_solver.C
	$(CC) $(CFLAGS) $^ $(SPECTRA) $(RTNORM)

loc_1d.o: loc_1d.C
	$(CC) $(CFLAGS) $^ $(SPECTRA) $(RTNORM)

square: laplace_solver_square.o
	$(CC) -o square $^ $(SPECTRA)

laplace_solver_square.o: laplace_solver_square.C
	$(CC) $(CFLAGS) $^ $(SPECTRA)

laplace_solver_esp.o: laplace_solver_esp.C
	$(CC) $(CFLAGS) $^ $(SPECTRA) $(RTNORM)


rtnorm.o: $(SRC1)
	$(CC) $(CFLAGS) $^

gen: koch estranho quad

koch: koch_gen.o
	$(CC) -o koch_gen $^ 

koch2: koch2_gen.o
	$(CC) -o koch2_gen $^ 


estranho: estranho_gen.o
	$(CC) -o estranho_gen $^ 

quad: square_gen.o
	$(CC) -o square_gen $^ 

koch_gen.o: koch_gen.C
	$(CC) $(CFLAGS) $^ $(SPECTRA)

koch2_gen.o: koch2_gen.C
	$(CC) $(CFLAGS) $^ $(SPECTRA)

estranho_gen.o: estranho_gen.C
	$(CC) $(CFLAGS) $^ $(SPECTRA)

square_gen.o: square_gen.C
	$(CC) $(CFLAGS) $^ $(SPECTRA)



.C.o:
	$(CC) $(CFLAGS) $<

clean:
	rm -f *.o
