#include "gauss_seidel.h"
#include <unistd.h>

double truncated_scalar_product(vector a, int *b, int start, int end)  {
    // if (a.size != b.size)   PRINT_ON_ERR("Vector sizes should be equal for a scalar product");
    if (start < 0 || end > a.size) PRINT_ON_ERR("Start and end indices should be coherent for a product");
    

    double res = 0;
    for (size_t i = start; i <= end; i++)  {
        res += a.vect[i] * b[i];
    }
    return res;
}

double truncated_vector_colmatrix_product(vector v, sparse_matrix P, int i, int a, int b) {
    double res = 0;
    node *ptr = NULL;
    
    ptr = P.lists[i];
    while (ptr != NULL)  {
        // On ne fait le produit que des lignes entre a et b (compris)
        if (ptr->row-1 < a) { 
            ptr = ptr->next;
            continue; 
        }
        if (ptr->row-1 > b) {
            ptr = ptr->next;
            continue;
        }
        res += v.vect[ptr->row-1] * ptr->val;
        ptr = ptr->next;
    }
    return res;
}


double truncated_x_ColG_product(vector npi, sparse_matrix P, double alpha, int i, int a, int b) {
    if (npi.size != P.size)    PRINT_ON_ERR("Matrix and vector sizes should be equal for a product");
    if (a < 0 || b > P.size) PRINT_ON_ERR("Start and end indices should be coherent for a product");
    
    if (a > b)      return 0;

    double res = 0;

    // vecteur xP (npi)
    res = truncated_vector_colmatrix_product(npi, P, i, a, b); // npi2 = npi*P

    // vecteur alpha*xP (alpha*npi)
    res *= alpha;
    
    // x*f = opi*f
    double xf;
    xf = truncated_scalar_product(npi, P.f, a, b);
    // Partie droite de la formule
    xf = (1-alpha)*(1.0/P.size) + alpha*(1.0/P.size)*xf;

    // Somme gauche et droite
    res += xf;
    
    return res;
}



vector gauss_seidel(const char* path, double epsilon, double alpha) {
    printf("\n\033[01;32mGauss-Seidel: Pows algorithm on a sparse matrix\033[m\n");

    // Lecture de la matrice creuse 
    sparse_matrix P;
    P = read_sparse_matrix(path);
    // print_sparse_matrix(P);
    
    vector npi = initialize_opi(P.size);
    
    // Itérations en do ... while
    vector npi2, sum;
    npi2.size = P.size;
    sum.size = P.size;
    npi2.vect = malloc(npi.size * sizeof(double));
    sum.vect = malloc(sum.size * sizeof(double));

    for (size_t i = 0; i < npi2.size; i++) {
        npi2.vect[i] = 0.0;
    }
    
    int it = 0;
    double tmp;
    do {
        it++;
        
        for (int i = 0; i < (int) npi2.size; i++) { // Pour chaque i de npi2
            npi2.vect[i] = truncated_x_ColG_product(npi2, P, alpha, i, 0, i-1)
                         + truncated_x_ColG_product(npi, P, alpha, i, i+1, P.size-1);
                        
            // ToDo : Diviser par 1-G[i,i]
        }
        
        
        // Normalisation du le vecteur
        tmp = 0;
        for (size_t i = 0; i < npi2.size; i++)  {
            tmp += npi2.vect[i];
        }
        // printf("sum avant normalisation: %lf\n", tmp);
        for (size_t i = 0; i < npi2.size; i++)  {
            npi2.vect[i] = npi2.vect[i] / tmp;
        }

        // Calcul de npi2 - npi
        for (size_t i = 0; i < sum.size; i++)  {
            sum.vect[i] = npi2.vect[i] - npi.vect[i];
        }

        // Copie de np2 dans npi
        for (size_t i = 0; i < P.size; i++)  {
            npi.vect[i] = npi2.vect[i];
        }
    } while(calculate_norm(sum) > epsilon);
    
    for (size_t i = 0; i < P.size; i++)
        free_list(P.lists[i]);
    free(P.lists);

    printf("Gauss-Seidel algorithm executed in %d iterations with ε = %f", it, epsilon);
    free(npi2.vect);
    free(sum.vect);
    return npi;
}

