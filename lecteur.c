#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define epsilon 1e-6

typedef struct matrix {
    int size;
    double **mat; 
} matrix;

typedef struct node {
   int row;
   int col;
   double val;
   struct node *next;
} node;

typedef struct sparse_matrix {
    int size;
    int edges;
    node **lists; 
} sparse_matrix;

typedef struct vector {
    int size;
    double *vect; 
} vector;

void print_vector(char *name, vector v);
void print_matrix(matrix P);
matrix read_full_matrix(const char *path);
sparse_matrix read_sparse_matrix(const char *path);

void matrix_vector_product(vector *res, vector v, matrix P);
double calculate_norm(vector v);
vector pows_algorithm(const char* path);


void print_vector(char *name, vector v)  {
    printf("\n");
    printf("%s = [ ", name);
    for (size_t i = 0; i < v.size; i++)  {
        printf("%lf ", v.vect[i]);
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

void print_node(node n)  {
    printf("Node:\n Row: %d\n Col: %d\n val: %lf\n", n.row, n.col, n.val);
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


void print_list(node *head) {
   node *ptr = head;
   
    printf("[ ");
    while(ptr != NULL) {
        printf("(%d,%d,%lf) -> ",ptr->row,ptr->col,ptr->val);
        ptr = ptr->next;
    }
    printf(" ]\n");
}

void free_list(node *head)  {
    struct node* tmp;

    while (head != NULL)  {
        tmp = head;
        head = head->next;
        free(tmp);
    }

}

void insert_first(int row, int col, double val, node **head) {
    node *new = (node*) malloc(sizeof(node));
    new->row = row;
    new->col = col;
    new->val = val;
    
    new->next = *head;
    *head = new;
}

sparse_matrix read_sparse_matrix(const char *path)  {
    FILE* fp;
    fp = fopen(path, "r");
 
    if (NULL == fp)  {
        printf("file can't be opened \n");
        exit(1);
    }
 
    sparse_matrix P;
    fscanf(fp, "%d\n", &P.size);
    fscanf(fp, "%d\n", &P.edges);
    printf("\n%d vertex and %d edges\n\n", P.size, P.edges);

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

    while ((read = getline(&line, &len, fp)) != -1)  {
        printf("\n%s", line);
        strToken = strtok (line, " ");
        row = atoi(strToken);
        strToken = strtok (NULL, " ");
        deg = atoi(strToken);
        
        printf ("Ligne: %d, degré: %d\n", row, deg);
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


int main(int argc, char const *argv[])  {
    if (argc < 3) {
        printf("Syntax: ./lecteur mode nom_fichier\n");
        return -1;
    }

    if (!strcmp(argv[1], "full")) {
        vector npi = pows_algorithm(argv[2]);  
        print_vector("res", npi);
    }
    else if (!strcmp(argv[1], "sparse")) {
        sparse_matrix P;
        P = read_sparse_matrix(argv[2]);
        
        // Affichage des listes
        for (size_t i = 0; i < P.size; i++)  {
            printf("\n> Colonne %ld: ", i+1);
            print_list(P.lists[i]);
        }

        for (size_t i = 0; i < P.size; i++)
            free_list(P.lists[i]);
        free(P.lists);
    }
    else {
        printf("Syntax: ./lecteur mode nom_fichier\n");
        return -1;
    }

    // ToDo : modifier l'algorithme pour matrices creuses

    return 0;
}
