CC=gcc
CFLAGS=-g -Wall

### Exécution ###

gauss_seidel: compil
	./pows gauss_seidel resources/wb-cs-stanford.txt
# sparse_matrix.txt

all_matrices: compil
	./pows all full_matrix.txt sparse_matrix.txt

full_matrix: compil
	./pows full full_matrix.txt

sparse_matrix: compil
	./pows sparse sparse_matrix.txt

pagerank: compil
	./pows pagerank sparse_matrix.txt


### Édition de lien ###

debug: compil
	valgrind --leak-check=full ./pows all full_matrix.txt sparse_matrix.txt

compil: pagerank.o lecteur.o main.o linked_list_manager.o gauss_seidel.o
	$(CC) -o pows main.o pagerank.o gauss_seidel.o lecteur.o linked_list_manager.o

### Compilation des fichiers ###

main.o: main.c lecteur.h pagerank.h gauss_seidel.h
	$(CC) $(CFLAGS) -c main.c

pagerank.o: pagerank.c pagerank.h lecteur.h linked_list_manager.h
	$(CC) $(CFLAGS) -c pagerank.c

gauss_seidel.o: gauss_seidel.c gauss_seidel.h lecteur.h linked_list_manager.h
	$(CC) $(CFLAGS) -c gauss_seidel.c

lecteur.o: lecteur.c lecteur.h linked_list_manager.h
	$(CC) $(CFLAGS) -c lecteur.c

linked_list_manager.o: linked_list_manager.c linked_list_manager.h
	$(CC) $(CFLAGS) -c linked_list_manager.c


### Utilitaires ###

clean:
	rm -f *.o
	rm -f pows
	ls -l
