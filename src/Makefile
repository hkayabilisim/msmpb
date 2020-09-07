all:
	gcc -O2 -g -fno-omit-frame-pointer -fsanitize=address -std=c11 msmpb.c forcefield?.c -lfftw3 -lm -o msmpb
	gcc -O2 -g -fno-omit-frame-pointer -fsanitize=address -std=c11 pme.c pme?.c -lfftw3 -lm -o pme
	gcc -O2 -g -fno-omit-frame-pointer -fsanitize=address -std=c11 testpme.c pme?.c -lfftw3 -lm -o testpme
	gcc -O2 -g -fno-omit-frame-pointer -fsanitize=address -std=c11 testmsm.c forcefield?.c -lfftw3 -lm -o testmsm
	gcc -O2 -g -fno-omit-frame-pointer -fsanitize=address -std=c11 fudge.c forcefield?.c -lfftw3 -lm -o fudge
demo:
	./msmpb ../data/NaClN8
