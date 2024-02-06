#include <stdlib.h>

#include "global.h"
#include "debug.h"

/**
 * @brief Validates command line arguments passed to the program.
 * @details This function will validate all the arguments passed to the
 * program, returning 0 if validation succeeds and -1 if validation fails.
 * Upon successful return, the various options that were specified will be
 * encoded in the global variable 'global_options', where it will be
 * accessible elsewhere in the program.  For details of the required
 * encoding, see the assignment handout.
 *
 * @param argc The number of arguments passed to the program from the CLI.
 * @param argv The argument strings passed to the program from the CLI.
 * @return 0 if validation succeeds and -1 if validation fails.
 * @modifies global variable "global_options" to contain an encoded representation
 * of the selected program options.
 */

int validargs(int argc, char **argv) {

    int ctr = argc - 1;
    argv++;
    char *currentArg = *argv;
    char *temp = *argv;
    char flag = ' ';
   
    if (ctr == 0) {
        return 0;
    }
   
    if (*temp == '-') { //Checks for -h flag
        temp++;
        flag = *temp;
        temp++;
        if ((flag == 'h') && (*temp == '\0')) {
            global_options |= HELP_OPTION;
            return 0;
        }
    }
   
    while (ctr > 0) {
        if (*currentArg == '-') {
            currentArg++;
            flag = *currentArg;
            currentArg++;
            if ((flag == 'm') && (*currentArg == '\0') && (global_options == 0x0)) { //Verifies -m flag
                global_options |= MATRIX_OPTION;
            } else if ((flag == 'n') && (*currentArg == '\0') && (global_options == 0x0)) { //Verifies -n flag
                global_options |= NEWICK_OPTION;
            } else if ((flag == 'o') && (*currentArg == '\0') && (global_options == 0x2)) { //Verifies -o flag follows -n flag
                ctr--;
                if (ctr == 0) {
                    return -1;
                }
                argv++;
                outlier_name = *argv;
            } else { //Returns -1 if an invalid flag is passed
                global_options &= 0;
                return -1;
            }
        } else { //Returns -1 if an invalid argument is passed
            global_options &= 0;
            return -1;
        }
        argv++;
        ctr--;
        currentArg = *argv;
    }
    return 0;
    abort();
}
