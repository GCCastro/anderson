CC = g++
SPECTRA = -I ~/SOFT/spectra/include
CFLAGS = -c -Wall -g

All:laplace 

laplace: laplace_solver.o
	$(CC) -o a $^ $(SPECTRA)

laplace_solver.o: laplace_solver.C
	$(CC) $(CFLAGS) $^ $(SPECTRA)

square: laplace_solver_square.o
	$(CC) -o square $^ $(SPECTRA)

laplace_solver_square.o: laplace_solver_square.C
	$(CC) $(CFLAGS) $^ $(SPECTRA)

gen: koch estranho

koch: koch_gen.o
	$(CC) -o koch_gen $^ 

estranho: estranho_gen.o
	$(CC) -o estranho_gen $^ 

koch_gen.o: koch_gen.C
	$(CC) $(CFLAGS) $^ $(SPECTRA)

estranho_gen.o: estranho_gen.C
	$(CC) $(CFLAGS) $^ $(SPECTRA)




.C.o:
	$(CC) $(CFLAGS) $<

clean:
	rm -f *.o
