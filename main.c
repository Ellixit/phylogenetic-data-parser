#include <stdio.h>
#include <stdlib.h>

#include "global.h"
#include "debug.h"

int main(int argc, char **argv)
{
    if(validargs(argc, argv))
        USAGE(*argv, EXIT_FAILURE);
    if(global_options == HELP_OPTION)
        USAGE(*argv, EXIT_SUCCESS);
    
    if (read_distance_data(stdin)) {
        exit(EXIT_FAILURE);
    }
    
    //Testing
    FILE *out = stdout;
    
    if (global_options == MATRIX_OPTION) {
        if (build_taxonomy(NULL)) {
            exit(EXIT_FAILURE);
        }
        if (emit_distance_matrix(out)) {
            exit(EXIT_FAILURE);
        } else {
            exit(EXIT_SUCCESS);
        }
    }
    
    if (global_options == NEWICK_OPTION) {
        if (build_taxonomy(NULL)) {
            exit(EXIT_FAILURE);
        }
        if (emit_newick_format(out)) {
            exit(EXIT_FAILURE);
        } else {
            exit(EXIT_SUCCESS);
        }
    }
    
    if (build_taxonomy(out)) {
        exit(EXIT_FAILURE);
    } else {
        exit(EXIT_SUCCESS);
    }
    
    return EXIT_FAILURE; 
}

/*
 * Just a reminder: All non-main functions should
 * be in another file not named main.c
 */
