shell:shell.o
	gcc -o a shell.o -fopenmp

shell.o:shell.c
	gcc -c -O3 shell.c -fopenmp
	
clean:
	rm -f a shell.o
