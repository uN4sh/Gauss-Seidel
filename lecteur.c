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

void full_matrix_vector_product(vector *res, vector v, matrix P);
void sparse_matrix_vector_product(vector *res, vector v, sparse_matrix P);

double calculate_norm(vector v);
vector initialize_opi(int size);

vector full_pows_algorithm(const char* path, double epsilon);
vector sparse_pows_algorithm(const char* path, double epsilon);


void print_vector(char *name, vector v)  {
    printf("\n");
    printf("%s = [ ", name);
    for (size_t i = 0; i < v.size; i++)  {
        printf("%lf ", v.vect[i]);
    }
    printf("]\n");
}

void print_full_matrix(matrix P)  {
    printf("%dx%d matrix\n", P.size, P.size);
    for (size_t i = 0; i < P.size; i++)  {
        printf(" [ ");
        for (size_t j = 0; j < P.size; j++)  {
            printf("%.2lf ", P.mat[i][j]);
        }
        printf("]\n");
    }
}

void print_sparse_matrix(sparse_matrix P)  {
    printf("%d vertices and %d edges\n", P.size, P.edges);
    for (size_t i = 0; i < P.size; i++)  {
        printf(" Colonne %ld: ", i+1);
        print_list(P.lists[i]);
    }
}

matrix read_full_matrix(const char *path)  {
    FILE* ptr;
    ptr = fopen(path, "r");
 
    if (NULL == ptr)    PRINT_ON_ERR("File can't be opened");

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

sparse_matrix read_sparse_matrix(const char *path)  {
    FILE* fp;
    fp = fopen(path, "r");
 
    if (NULL == fp)     PRINT_ON_ERR("File can't be opened");
    
    sparse_matrix P;
   
    fscanf(fp, "%d\n", &P.size);
    fscanf(fp, "%d\n", &P.edges);
    
    // Format: ind_ligne deg_ligne ind_col val_col ind_col val_col ...
    int row, col, deg;
    double val;
    char *line;
    size_t len = 0;
    size_t read;
    char *strToken;
    int flip; // If flip: read index, else: read value

    // Création des têtes de liste pour chaque colonne de la matrice
    P.lists = malloc(P.size * sizeof(node*));
    for (size_t i = 0; i < P.size; i++)  {
        P.lists[i] = NULL; 
    }

    // ToDo: passer le vecteur f en vecteur d'entiers
    P.f = malloc(P.size * sizeof(int));
    
    while ((read = getline(&line, &len, fp)) != -1)  {
        strToken = strtok (line, " ");
        row = atoi(strToken);
        strToken = strtok (NULL, " ");
        deg = atoi(strToken); 

        // Création du vecteur f : f[i] = 1 si deg = 0
        if (deg)
            P.f[row-1] = 0;
        else
            P.f[row-1] = 1;
        
        flip = 1;
        while ( (strToken = strtok ( NULL, " " )) != NULL ) {
            if (flip)
                col = atoi(strToken);
            else {
                val = strtod(strToken, NULL);
                insert_first(row, col, val, &P.lists[col-1]);
            }
            flip = 1 - flip;
        }
    }

    fclose(fp);
    if (line)
        free(line);
    
    return P;
}

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


int main(int argc, char const *argv[])  {
    if (argc < 3) {
        printf("Syntax: ./lecteur mode nom_fichier\n");
        return -1;
    }

    double epsilon = 1e-6;
    double alpha = 0.85;

    if (!strcmp(argv[1], "full")) {
        vector npi = full_pows_algorithm(argv[2], epsilon);  
        print_vector("res", npi);
    }
    else if (!strcmp(argv[1], "sparse")) {
        vector npi = sparse_pows_algorithm(argv[2], epsilon);  
        print_vector("res", npi);
    }
    else if (!strcmp(argv[1], "all")) {
        vector npi = full_pows_algorithm(argv[2], epsilon);  
        print_vector("res_full", npi);
        free(npi.vect);
        npi = sparse_pows_algorithm(argv[3], epsilon);  
        print_vector("res_sparse", npi);
        free(npi.vect);
    }
    else if (!strcmp(argv[1], "pagerank")) {
        vector npi = pagerank_98(argv[2], alpha, epsilon);
        print_vector("res", npi);
    }
    else {
        printf("Syntax: ./lecteur mode nom_fichier\n");
        return -1;
    }

    return 0;
}
