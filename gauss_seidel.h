#ifndef __gauss_seidel_h
#define __gauss_seidel_h

#include "pagerank.h"

double truncated_scalar_product(vector a, int *b, int start, int end);

/**
 * @brief Produit vecteur matrice en tronquant la matrice de a à b
 * 
 * @param res 
 * @param v 
 * @param P 
 * @param a
 * @param b  
 */
double truncated_vector_colmatrix_product(vector v, sparse_matrix P, int i, int a, int b);

double truncated_x_ColG_product(vector npi, sparse_matrix P, double alpha, int i, int a, int b);

vector gauss_seidel(const char* path, double epsilon, double alpha);

#endif
