CC=gcc
CFLAGS=-g -Wall

all_matrices: compil
	./pows all full_matrix.txt sparse_matrix.txt

full_matrix: compil
	./pows full full_matrix.txt

sparse_matrix: compil
	./pows sparse sparse_matrix.txt


debug: compil
	valgrind --leak-check=full ./pows all full_matrix.txt sparse_matrix.txt

compil: lecteur.o 
	$(CC) -o pows lecteur.o linked_list_manager.h

clean:
	rm -f *.o
	rm -f pows
	ls -l
