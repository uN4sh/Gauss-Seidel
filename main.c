#include "pagerank.h"
#include "gauss_seidel.h"

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
        vector npi = pagerank_98(argv[2], epsilon, alpha);
        print_vector("res", npi);
    }
    else if (!strcmp(argv[1], "gauss_seidel")) {
        pagerank_98(argv[2], epsilon, alpha);
        // print_vector("res", npi);

        gauss_seidel(argv[2], epsilon, alpha);
        // print_vector("res", npi);
    }
    else {
        printf("Syntax: ./lecteur mode nom_fichier\n");
        return -1;
    }

    return 0;
}
