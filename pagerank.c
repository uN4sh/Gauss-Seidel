#include "pagerank.h"

void full_matrix_vector_product(vector *res, vector v, matrix P)  {
    if (v.size != P.size)   PRINT_ON_ERR("Matrix and vector sizes should be equal for a product");

    for (size_t i = 0; i < res->size; i++)  {
        res->vect[i] = 0;
        for (size_t j = 0; j < res->size; j++)  {
            res->vect[i] += v.vect[j] * P.mat[j][i];
        }
    }
}

void sparse_matrix_vector_product(vector *res, vector v, sparse_matrix P)  {
    if (v.size != P.size)   PRINT_ON_ERR("Matrix and vector sizes should be equal for a product");

    node *ptr = NULL;
    for (size_t i = 0; i < P.size; i++)  {
        res->vect[i] = 0;
        ptr = P.lists[i];
        while (ptr != NULL)  {
            res->vect[i] += v.vect[ptr->row-1] * ptr->val;
            ptr = ptr->next;
        }
    }
}

// La norme 1 d’un vecteur x est la somme de ses valeurs (en valeur absolue)
double calculate_norm(vector v)  {
    double norm = 0;
    for (size_t i = 0; i < v.size; i++)  {
        if (v.vect[i] > 0)
            norm += v.vect[i];
        else
            norm -= v.vect[i];
    }
    return norm;
}

vector initialize_opi(int size)  {
    // Vecteur e plein de 1
    vector e;
    e.size = size;
    e.vect = malloc (e.size * sizeof(double));
    for (size_t i = 0; i < e.size; i++)  {
        e.vect[i] = 1;
    }
    
    // Initialisation (opi = e/N)
    vector opi;
    opi.size = size;
    opi.vect =  malloc (opi.size * sizeof(double));
    for (size_t i = 0; i < opi.size; i++)  {
        opi.vect[i] = (double) e.vect[i] / size;
    }
    free(e.vect);
    return opi;
}


vector full_pows_algorithm(const char* path, double epsilon) {
    printf("\n\033[01;32mPows algorithm on a full matrix\033[m\n");

    // Lecture de la matrice pleine 
    matrix P;
    P = read_full_matrix(path);
    print_full_matrix(P);

    vector opi = initialize_opi(P.size);
    
    // Itérations en do ... while
    vector npi, npi2, sum;
    npi2 = opi;
    npi.size = P.size;
    sum.size = P.size;
    npi.vect = malloc(npi.size * sizeof(double));
    sum.vect = malloc(sum.size * sizeof(double));
    int it = 0;
    do {
        it++;
        // Copie de pi-1 dans pi
        for (size_t i = 0; i < P.size; i++)  {
            npi.vect[i] = npi2.vect[i];
        }

        full_matrix_vector_product(&npi2, npi, P); // npi = opi*P

        // Calcul de npi2 -npi
        for (size_t i = 0; i < sum.size; i++)  {
            sum.vect[i] = npi2.vect[i] - npi.vect[i];
        }
    } while (calculate_norm(sum) > epsilon);

    for (size_t i = 0; i < P.size; i++)
        free(P.mat[i]);
    free(P.mat);
    
    printf("Pows algorithm executed in %d iterations with ε = %f", it, epsilon);
    print_vector("opi", opi);
    free(npi.vect);
    free(sum.vect);
    return npi2;
}

vector sparse_pows_algorithm(const char* path, double epsilon) {
    printf("\n\033[01;32mPows algorithm on a sparse matrix\033[m\n");

    // Lecture de la matrice creuse 
    sparse_matrix P;
    P = read_sparse_matrix(path);
    print_sparse_matrix(P);

    vector opi = initialize_opi(P.size);
    
    // Itérations en do ... while
    vector npi, npi2, sum;
    npi2 = opi;
    npi.size = P.size;
    sum.size = P.size;
    npi.vect = malloc(npi.size * sizeof(double));
    sum.vect = malloc(sum.size * sizeof(double));
    int it = 0;
    do {
        it++;
        // Copie de pi-1 dans pi
        for (size_t i = 0; i < P.size; i++)  {
            npi.vect[i] = npi2.vect[i];
        }

        sparse_matrix_vector_product(&npi2, npi, P); // npi = opi*P
        
        // Calcul de npi2 -npi
        for (size_t i = 0; i < sum.size; i++)  {
            sum.vect[i] = npi2.vect[i] - npi.vect[i];
        }
    } while (calculate_norm(sum) > epsilon);

    for (size_t i = 0; i < P.size; i++)
        free_list(P.lists[i]);
    free(P.lists);

    printf("Pows algorithm executed in %d iterations with ε = %f", it, epsilon);
    print_vector("opi", opi);
    free(npi.vect);
    free(sum.vect);
    return npi2;
}


double scalar_product(vector a, int *b)  {
    // if (a.size != b.size)   PRINT_ON_ERR("Vector sizes should be equal for a scalar product");

    double res = 0;
    for (size_t i = 0; i < a.size; i++)  {
        res += a.vect[i] * b[i];
    }
    return res;
}

vector pagerank_98(const char* path, double epsilon, double alpha) {
    printf("\n\033[01;32mPageRank 98: Pows algorithm on a sparse matrix\033[m\n");

    // Lecture de la matrice creuse 
    sparse_matrix P;
    P = read_sparse_matrix(path);
    // print_sparse_matrix(P);
    
    vector opi = initialize_opi(P.size);
    
    // Itérations en do ... while
    vector npi, npi2, sum;
    npi2 = opi;
    npi.size = P.size;
    sum.size = P.size;
    npi.vect = malloc(npi.size * sizeof(double));
    sum.vect = malloc(sum.size * sizeof(double));
    int it = 0;
    do {
        it++;
        // Copie de pi-1 dans pi
        for (size_t i = 0; i < P.size; i++)  {
            npi.vect[i] = npi2.vect[i];
        }

        // vecteur xP (npi)
        sparse_matrix_vector_product(&npi2, npi, P); // npi = opi*P

        // vecteur alpha*xP (alpha*npi)
        for (size_t i = 0; i < P.size; i++)  {
            npi2.vect[i] *= alpha;
        }

        // x*f = opi*f
        double xf;
        xf = scalar_product(npi, P.f);
        // Partie droite de la formule
        xf = (1-alpha)*(1.0/P.size) + alpha*(1.0/P.size)*xf;
        
        // Somme gauche et droite
        for (size_t i = 0; i < P.size; i++)  {
            npi2.vect[i] = npi2.vect[i] + xf;
        }
        
        // Calcul de npi2 -npi
        for (size_t i = 0; i < sum.size; i++)  {
            sum.vect[i] = npi2.vect[i] - npi.vect[i];
        }
    } while (calculate_norm(sum) > epsilon);

    for (size_t i = 0; i < P.size; i++)
        free_list(P.lists[i]);
    free(P.lists);

    printf("Pows algorithm executed in %d iterations with ε = %f", it, epsilon);
    free(npi.vect);
    free(sum.vect);
    return npi2;
}
