#include "lecteur.h"

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

