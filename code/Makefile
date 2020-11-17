all:
	gcc -O3 -std=c11 msmpb.c forcefield?.c -lfftw3 -lm -o msmpb
	gcc -O3 -std=c11 pme.c pme?.c -lfftw3 -lm -o pme
demo:
	./msmpb ../data/NaClN8
