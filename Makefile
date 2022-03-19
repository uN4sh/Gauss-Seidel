CC=gcc
CFLAGS=-g -Wall

full_matrix: compil
	./pows full full_matrix.txt

sparse_matrix: compil
	./pows sparse sparse_matrix.txt


debug: compil
	valgrind ./pows sparse sparse_matrix.txt
	
compil: lecteur.o
	$(CC) -o pows lecteur.o 

clean:
	rm -f *.o
	rm -f pows
	ls -l
