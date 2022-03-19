#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define epsilon 1e-6

typedef struct matrix {
    int size;
    double **mat; 
} matrix;

typedef struct vector {
    int size;
    double *vect; 
} vector;

void print_vector(char *name, vector v);
void print_matrix(matrix P);
matrix read_full_matrix(const char *path);
void matrix_vector_product(vector *res, vector v, matrix P);
double calculate_norm(vector v);
vector pows_algorithm(const char* path);


void print_vector(char *name, vector v)  {
    printf("\n");
    printf("%s = [ ", name);
    for (size_t i = 0; i < v.size; i++)  {
        printf("%.7lf ", v.vect[i]);
    }
    printf("]\n");
}

void print_matrix(matrix P)  {
    printf("Size: %dx%d\n", P.size, P.size);
    for (size_t i = 0; i < P.size; i++)  {
        for (size_t j = 0; j < P.size; j++)  {
            printf("%.2lf ", P.mat[i][j]);
        }
        printf("\n");
    }
}

matrix read_full_matrix(const char *path)  {
    FILE* ptr;
    ptr = fopen(path, "r");
 
    if (NULL == ptr)  {
        printf("file can't be opened \n");
        exit(1);
    }
 
    matrix P;
    fscanf(ptr, "%d", &P.size);

    // Malloc matrix vertex*vertex
    P.mat = malloc(P.size*sizeof(double*));
    for (size_t i = 0; i < P.size; i++)  {
        P.mat[i] = malloc(P.size*sizeof(double));
    }
    
    // Read matrix line by line
    for (size_t i = 0; i < P.size; i++)  {
        fscanf(ptr, "%lf %lf %lf %lf %lf %lf", &P.mat[i][0], &P.mat[i][1], &P.mat[i][2], &P.mat[i][3], &P.mat[i][4], &P.mat[i][5]);
    }

    fclose(ptr);
    return P;
}

matrix read_sparse_matrix(const char *path)  {

}

void matrix_vector_product(vector *res, vector v, matrix P)  {
    // ToDo: ajouter des checks de sizes ?

    for (size_t i = 0; i < res->size; i++)  {
        res->vect[i] = 0;
        // printf("[%ld] = ", i);
        for (size_t j = 0; j < res->size; j++)  {
            // printf("%lf * %lf + ", vect[j], P.mat[j][i]);
            res->vect[i] += v.vect[j] * P.mat[j][i];
        }
        // printf("\n");
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

vector pows_algorithm(const char* path) {
    // Lecture de la matrice pleine 
    matrix P;
    P = read_full_matrix(path);
    print_matrix(P);

    // Vecteur e plein de 1
    vector e;
    e.size = P.size;
    e.vect = malloc (e.size * sizeof(double));
    for (size_t i = 0; i < e.size; i++)  {
        e.vect[i] = 1;
    }
    
    // Initialisation (opi = e/N)
    vector opi;
    opi.size = P.size;
    opi.vect = malloc (opi.size * sizeof(double));
    for (size_t i = 0; i < opi.size; i++)  {
        opi.vect[i] = (double) e.vect[i] / P.size;
    }
    
    print_vector("opi", opi);
    
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
        npi.vect = malloc(npi.size * sizeof(double));
        for (size_t i = 0; i < P.size; i++)  {
            npi.vect[i] = npi2.vect[i];
        }

        matrix_vector_product(&npi2, npi, P); // npi = opi*P
        // print_vector("npi2", npi2);
        
        // Calcul de npi2 -npi
        for (size_t i = 0; i < sum.size; i++)  {
            sum.vect[i] = npi2.vect[i] - npi.vect[i];
        }

        // printf("norm: %.7lf\n", calculate_norm(sum));
        
    } while ( calculate_norm(sum) > epsilon);

    printf("Pows algorithm executed in %d iterations with ε = %f", it, epsilon);
    return npi2;
}



int main(int argc, char const *argv[]) {
    
    if (argc < 2) {
        printf("Syntax: ./lecteur nom_fichier\n");
        return -1;
    }

    vector npi = pows_algorithm(argv[1]);  
    print_vector("res", npi);

    // ToDo : modifier l'algorithme pour matrices creuses

    return 0;
}
