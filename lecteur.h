#ifndef __lecteur_h
#define __lecteur_h

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "linked_list_manager.h"

#define PRINT_ON_ERR(s)  { fprintf(stderr, "ERR -- %s\n", s); exit(1); }


typedef struct vector {
    int size;
    double *vect; 
} vector;

typedef struct matrix {
    int size;
    double **mat; 
} matrix;

typedef struct sparse_matrix {
    int size;
    int edges;
    int *f; // Vecteur f
    node **lists; 
} sparse_matrix;

typedef struct pagerank_struct {
    sparse_matrix P;
    vector f;
} pagerank_struct;

void print_vector(char *name, vector v);
void print_full_matrix(matrix P);
void print_sparse_matrix(sparse_matrix P);

matrix read_full_matrix(const char *path);
sparse_matrix read_sparse_matrix(const char *path);


#endif
