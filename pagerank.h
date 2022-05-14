#ifndef __pagerank_h
#define __pagerank_h

#include <time.h>
#include "lecteur.h"

double calculate_norm(vector v);
vector initialize_opi(int size);

void full_matrix_vector_product(vector *res, vector v, matrix P);
void sparse_matrix_vector_product(vector *res, vector v, sparse_matrix P);

vector full_pows_algorithm(const char* path, double epsilon);
vector sparse_pows_algorithm(const char* path, double epsilon);

double scalar_product(vector a, int *b);
vector pagerank_98(const char* path, double epsilon, double alpha);


#endif
