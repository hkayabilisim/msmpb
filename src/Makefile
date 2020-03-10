all:
	gcc -O2 -std=c99 msmpb.c forcefield?.c -lfftw3 -lm -o msmpb
	gcc -O2 -std=c99 pme.c pme?.c -lfftw3 -lm -o pme
demo:
	./msmpb ../data/NaClN8
