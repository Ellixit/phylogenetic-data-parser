#include <stdlib.h>
#include <math.h>

#include "global.h"
#include "debug.h"

// Helper Function Prototypes
static void copyString(char *source, char *destination);
static int equalString(char *str1, char *str2);
static int duplicateName(char *str, char *nameArray, int numTaxa);
static double calculateNum(char *source);
static double qIJ(int indexI, int indexJ, int i, int j);
static double dUK(int indexF, int indexG, int f, int g);
static double accessIndices(int i, int j);
static void writeNodeName();
static int logn(int n);
static void printNewick(NODE *node, NODE *parent, FILE *out);
static int findIndexName(char *inputName);
static int findOutlier();
static int compareNode(NODE *node, NODE *parent);
static int isLeaf(NODE *node);

// Global Variables
int invalidNum = 0;
int firstLeafComma = 1;
int firstOpenParen = 0;

/**
 * @brief  Read genetic distance data and initialize data structures.
 * @details  This function reads genetic distance data from a specified
 * input stream, parses and validates it, and initializes internal data
 * structures.
 *
 * The input format is a simplified version of Comma Separated Values
 * (CSV).  Each line consists of text characters, terminated by a newline.
 * Lines that start with '#' are considered comments and are ignored.
 * Each non-comment line consists of a nonempty sequence of data fields;
 * each field iis terminated either by ',' or else newline for the last
 * field on a line.  The constant INPUT_MAX specifies the maximum number
 * of data characters that may be in an input field; fields with more than
 * that many characters are regarded as invalid input and cause an error
 * return.  The first field of the first data line is empty;
 * the subsequent fields on that line specify names of "taxa", which comprise
 * the leaf nodes of a phylogenetic tree.  The total number N of taxa is
 * equal to the number of fields on the first data line, minus one (for the
 * blank first field).  Following the first data line are N additional lines.
 * Each of these lines has N+1 fields.  The first field is a taxon name,
 * which must match the name in the corresponding column of the first line.
 * The subsequent fields are numeric fields that specify N "distances"
 * between this taxon and the others.  Any additional lines of input following
 * the last data line are ignored.  The distance data must form a symmetric
 * matrix (i.e. D[i][j] == D[j][i]) with zeroes on the main diagonal
 * (i.e. D[i][i] == 0).
 *
 * If 0 is returned, indicating data successfully read, then upon return
 * the following global variables and data structures have been set:
 *   num_taxa - set to the number N of taxa, determined from the first data line
 *   num_all_nodes - initialized to be equal to num_taxa
 *   num_active_nodes - initialized to be equal to num_taxa
 *   node_names - the first N entries contain the N taxa names, as C strings
 *   distances - initialized to an NxN matrix of distance values, where each
 *     row of the matrix contains the distance data from one of the data lines
 *   nodes - the "name" fields of the first N entries have been initialized
 *     with pointers to the corresponding taxa names stored in the node_names
 *     array.
 *   active_node_map - initialized to the identity mapping on [0..N);
 *     that is, active_node_map[i] == i for 0 <= i < N.
 *
 * @param in  The input stream from which to read the data.
 * @return 0 in case the data was successfully read, otherwise -1
 * if there was any error.  Premature termination of the input data,
 * failure of each line to have the same number of fields, and distance
 * fields that are not in numeric format should cause a one-line error
 * message to be printed to stderr and -1 to be returned.
 */

int read_distance_data(FILE *in) {
    
    // File parsing
    int currentChar = 0;
    int prevChar = 0;
    int inputLength = INPUT_MAX;
    
    // Logical flags
    int startLineFlag = 1;
    int skipLineFlag = 0;
    int emptyNameFieldFlag = 0;
    int dataInputFlag = 0;
    
    // Data transcription
    int numTaxa = 0;
    int rowLeft = -1;
    int colLeft = 0;
    int numTranscribed = 0;
    int periodCount = 0;
    int ctr = 0;
    int tracker = 0;
    
    // Pointers
    char *bufferPtr = input_buffer;
    double *distancePtr = *distances;
    char *namePtr = *node_names;
    NODE *nodePtr = nodes;
    
    while (1) {
        currentChar = fgetc(in);
        if (feof(in)) {
            break;
        }
        if (rowLeft == 0) {
            break;
        }
        
        // SKIP COMMENT LINES
        if ((startLineFlag == 1) && (currentChar == '#')) { // '#' at start of line
            skipLineFlag = 1;
            startLineFlag = 0;
            continue;
        } else if ((currentChar == '\n') && (skipLineFlag == 1)) { // '\n' at end of line
            skipLineFlag = 0;
            startLineFlag = 1;
            continue;
        } else if (skipLineFlag == 1) { // In comment line
            continue;
        }
        
        // SKIP FIRST FIELD
        if (emptyNameFieldFlag == 0) { // ',' at end of first field
            if (currentChar == ',') {
                emptyNameFieldFlag = 1;
                prevChar = currentChar;
                startLineFlag = 0;
                continue;
            } else {
                fprintf(stderr, "ERROR: NON-EMPTY FIRST FIELD");
                return -1;
            }
        }
        
        // TRANSCRIBING FIRST DATA LINE: TAXA NAMES
        if (dataInputFlag == 0) {
            if ((currentChar == ',') || (currentChar == '\n')) {
                goto TRANSCRIBE_TAXON_NAMES;
                COMPLETE_FIRST_LINE:
                    namePtr = *node_names;
                    rowLeft = numTaxa;
                    colLeft = numTaxa;
                    startLineFlag = 1;
                    dataInputFlag = 1;
                    continue;
            } else {
                *bufferPtr = currentChar;
                bufferPtr++;
                inputLength--;
            }
            goto NEXT_ITERATION;
        }
        
        // TRANSCRIBING DISTANCE DATA
        if (startLineFlag == 1) {
            if (currentChar == ',') {
                *bufferPtr = '\0';
                bufferPtr = input_buffer;
                if (!equalString(bufferPtr, namePtr)) {
                    fprintf(stderr, "ERROR: ROW AND COLUMN NAME MISMATCH");
                    return -1;
                }
                startLineFlag = 0;
                namePtr += (INPUT_MAX+1);
                inputLength = INPUT_MAX;
                prevChar = currentChar;
                continue;
            } else {
                *bufferPtr = currentChar;
                bufferPtr++;
                inputLength--;
                goto NEXT_ITERATION;
            }
        }
        if ((currentChar == ',') || (currentChar == '\n')) { 
            if (prevChar == ',') {
                fprintf(stderr, "ERROR: MISSING FIELD");
                return -1;
            }
            goto TRANSCRIBE_NUMERIC_FIELD;
            COMPLETE_DATA_LINE:
                if (colLeft != 0) {
                    fprintf(stderr, "ERROR: MISSING FIELD");
                    return -1;
                }
                distancePtr += (MAX_NODES - numTaxa);
                colLeft = numTaxa;
                rowLeft--;
                startLineFlag = 1;
                continue;
        } else if ((currentChar == '.') && (periodCount == 0)) {
            *bufferPtr = currentChar;
            bufferPtr++;
            periodCount++;
            inputLength--;
        } else if ((currentChar >= 48) && (currentChar <= 57)) {
            *bufferPtr = currentChar;
            bufferPtr++;
            inputLength--;
        } else {
            fprintf(stderr, "ERROR: UNACCEPTED NUMERIC CHARACTER");
            return -1;
        }
        goto NEXT_ITERATION;
        
        TRANSCRIBE_TAXON_NAMES:
            if (prevChar == ',') {
                fprintf(stderr, "ERROR: MISSING TAXON NAME");
                return -1;
            }
            *bufferPtr = '\0';
            bufferPtr = input_buffer;
            copyString(bufferPtr, namePtr);
            if (duplicateName(namePtr, *node_names, numTaxa)) {
                fprintf(stderr, "ERROR: DUPLICATE TAXON NAME");
                return -1;
            }
            nodePtr->name = namePtr;
            nodePtr++;
            numTaxa++;
            inputLength = INPUT_MAX;
            namePtr += (INPUT_MAX+1);
            prevChar = currentChar;
            if (currentChar == '\n') goto COMPLETE_FIRST_LINE;
            continue;
        
        TRANSCRIBE_NUMERIC_FIELD:
            *bufferPtr = '\0';
            bufferPtr = input_buffer;
            *distancePtr = calculateNum(bufferPtr);
            if (invalidNum == 1) { 
                fprintf(stderr, "ERROR: INVALID NUMERIC FORMAT");
                return -1;
            }
            distancePtr++;
            inputLength = INPUT_MAX;
            periodCount = 0;
            colLeft--;
            numTranscribed++;
            prevChar = currentChar;
            if (currentChar == '\n') goto COMPLETE_DATA_LINE;
            if (colLeft == 0) {
                fprintf(stderr, "ERROR: INCORRECT NUMBER OF FIELDS");
                return -1;
            }
            continue;
        
        NEXT_ITERATION:
            if (inputLength == -1) {
                fprintf(stderr, "ERROR: MAXIMUM INPUT SIZE EXCEEDED");
                return -1;
            }
            prevChar = currentChar;
            continue;   
    }
    if ((rowLeft != 0) || (colLeft != numTaxa)) {
        fprintf(stderr, "ERROR: PREMATURE TERMINATION OF INPUT DATA");
        return -1;
    }
    if (numTranscribed == 0) {
        fprintf(stderr, "ERROR: NO TAXONOMIC DATA");
        return -1;
    }
    
    double *ptr1 = *distances;
    double *ptr2 = *distances;
    tracker = numTaxa;
    ctr = (MAX_NODES - numTaxa);
    for (int i = 0; i < numTaxa; i++) { // Outer loop iterates through ROWs
        for (int j = 0; j < tracker; j++) { // Inner loop iterates through COLs
            if ((ptr1 == ptr2) && (fabs(*ptr1) > 0.00000000000000000000001)) { // Validates 0 diagonal
                fprintf(stderr, "ERROR: NON-ZERO DIAGONAL");
                return -1;
            }
            if (fabs(*ptr1 - *ptr2) > 0.00000000000000000000001) { // Validates symmetricity
                fprintf(stderr, "ERROR: ASYMMETRICAL MATRIX");
                return -1;
            }
            ptr1++;
            ptr2 += MAX_NODES;
        }
        tracker--;
        ctr++;
        ptr1 += ctr;
        ptr2 = ptr1;
    }
    num_taxa = numTaxa;
    num_all_nodes = num_taxa;
    num_active_nodes = num_taxa;
    return 0;
    abort();
}

/**
 * @brief  Emit a representation of the phylogenetic tree in Newick
 * format to a specified output stream.
 * @details  This function emits a representation in Newick format
 * of a synthesized phylogenetic tree to a specified output stream.
 * See (https://en.wikipedia.org/wiki/Newick_format) for a description
 * of Newick format.  The tree that is output will include for each
 * node the name of that node and the edge distance from that node
 * its parent.  Note that Newick format basically is only applicable
 * to rooted trees, whereas the trees constructed by the neighbor
 * joining method are unrooted.  In order to turn an unrooted tree
 * into a rooted one, a root will be identified according by the
 * following method: one of the original leaf nodes will be designated
 * as the "outlier" and the unique node adjacent to the outlier
 * will serve as the root of the tree.  Then for any other two nodes
 * adjacent in the tree, the node closer to the root will be regarded
 * as the "parent" and the node farther from the root as a "child".
 * The outlier node itself will not be included as part of the rooted
 * tree that is output.  The node to be used as the outlier will be
 * determined as follows:  If the global variable "outlier_name" is
 * non-NULL, then the leaf node having that name will be used as
 * the outlier.  If the value of "outlier_name" is NULL, then the
 * leaf node having the greatest total distance to the other leaves
 * will be used as the outlier.
 *
 * @param out  Stream to which to output a rooted tree represented in
 * Newick format.
 * @return 0 in case the output is successfully emitted, otherwise -1
 * if any error occurred.  If the global variable "outlier_name" is
 * non-NULL, then it is an error if no leaf node with that name exists
 * in the tree.
 */

int emit_newick_format(FILE *out) {
    
    // Calculations
    int outlierIndex = 0;
    int parentIndex = 0;
    
    // Pointers
    NODE *outlier = nodes;
    NODE *parent = nodes;
    
    // Setting outlier node
    if (outlier_name != NULL) {
        outlierIndex = findIndexName(outlier_name);
        if (outlierIndex == -1) {
            fprintf(stderr, "ERROR: NO NODE WITH NAME %s EXISTS", outlier_name);
            return -1;
        }
    } else {
        outlierIndex = findOutlier();
    }
    outlier += outlierIndex;
    
    
    parent = *(outlier->neighbors);
    printNewick(parent, outlier, out);
    fprintf(out, "\n");
    
    return 0;
    abort();
}

/**
 * @brief  Emit the synthesized distance matrix as CSV.
 * @details  This function emits to a specified output stream a representation
 * of the synthesized distance matrix resulting from the neighbor joining
 * algorithm.  The output is in the same CSV form as the program input.
 * The number of rows and columns of the matrix is equal to the value
 * of num_all_nodes at the end of execution of the algorithm.
 * The submatrix that consists of the first num_leaves rows and columns
 * is identical to the matrix given as input.  The remaining rows and columns
 * contain estimated distances to internal nodes that were synthesized during
 * the execution of the algorithm.
 *
 * @param out  Stream to which to output a CSV representation of the
 * synthesized distance matrix.
 * @return 0 in case the output is successfully emitted, otherwise -1
 * if any error occurred.
 */

int emit_distance_matrix(FILE *out) {
    
    // Pointers
    double *distancePtr = *distances;
    char *namePtr = *node_names;
    
    // Output node_names
    for (int i = 0; i < num_all_nodes; i++) { // Iterate through node_names
        fprintf(out, ",%s", namePtr);
        namePtr += (INPUT_MAX+1);
    }
    namePtr = *node_names;
    fprintf(out, "\n");
    fflush(out);
    
    // Output distances
    for (int i = 0; i < num_all_nodes; i++) { // Iterate through ROWs of distances
        fprintf(out, "%s", namePtr);
        namePtr += (INPUT_MAX+1);
        for (int j = 0; j < num_all_nodes; j++) { // Iterate through COLs of distances
            fprintf(out, ",%.2f", *distancePtr);
            distancePtr++;
        }
        distancePtr += MAX_NODES - num_all_nodes;
        fprintf(out, "\n");
        fflush(out);
    }
    
    return 0;
    abort();
}

/**
 * @brief  Build a phylogenetic tree using the distance data read by
 * a prior successful invocation of read_distance_data().
 * @details  This function assumes that global variables and data
 * structures have been initialized by a prior successful call to
 * read_distance_data(), in accordance with the specification for
 * that function.  The "neighbor joining" method is used to reconstruct
 * phylogenetic tree from the distance data.  The resulting tree is
 * an unrooted binary tree having the N taxa from the original input
 * as its leaf nodes, and if (N > 2) having in addition N-2 synthesized
 * internal nodes, each of which is adjacent to exactly three other
 * nodes (leaf or internal) in the tree.  As each internal node is
 * synthesized, information about the edges connecting it to other
 * nodes is output.  Each line of output describes one edge and
 * consists of three comma-separated fields.  The first two fields
 * give the names of the nodes that are connected by the edge.
 * The third field gives the distance that has been estimated for
 * this edge by the neighbor-joining method.  After N-2 internal
 * nodes have been synthesized and 2*(N-2) corresponding edges have
 * been output, one final edge is output that connects the two
 * internal nodes that still have only two neighbors at the end of
 * the algorithm.  In the degenerate case of N=1 leaf, the tree
 * consists of a single leaf node and no edges are output.  In the
 * case of N=2 leaves, then no internal nodes are synthesized and
 * just one edge is output that connects the two leaves.
 *
 * Besides emitting edge data (unless it has been suppressed),
 * as the tree is built a representation of it is constructed using
 * the NODE structures in the nodes array.  By the time this function
 * returns, the "neighbors" array for each node will have been
 * initialized with pointers to the NODE structure(s) for each of
 * its adjacent nodes.  Entries with indices less than N correspond
 * to leaf nodes and for these only the neighbors[0] entry will be
 * non-NULL.  Entries with indices greater than or equal to N
 * correspond to internal nodes and each of these will have non-NULL
 * pointers in all three entries of its neighbors array.
 * In addition, the "name" field each NODE structure will contain a
 * pointer to the name of that node (which is stored in the corresponding
 * entry of the node_names array).
 *
 * @param out  If non-NULL, an output stream to which to emit the edge data.
 * If NULL, then no edge data is output.
 * @return 0 in case the output is successfully emitted, otherwise -1
 * if any error occurred.
 */

int build_taxonomy(FILE *out) {
    
    // Calculations
    double rowTemp = 0;
    double minPair = 0;
    double qTemp = 0;
    int indexI = 0;
    int indexJ = 0;
    int indexF = 0; // Used to index node f in distances
    int indexG = 0; // Used to index node g in distances
    int indexK = 0;
    int indexU = 0;
    int activeNodeI = 0; // Used to index node f in active_node_map
    int activeNodeJ = 0; // Used to index node g in active_node_map
    int currentIndex = 0;
    double distanceVal = 0;
    
    // Pointers
    double *rowSumPtr = row_sums;
    double *uCol = *distances;
    double *uRow = *distances;
    int *nodeMapPtr = active_node_map;
    int *nodeMapI = active_node_map;
    int *nodeMapJ = active_node_map;
    NODE *nodePtrF = nodes;
    NODE *nodePtrG = nodes;
    NODE *nodePtrU = nodes;
    
    // Flags
    int initializeMinPairFlag = 1;
    int outputEdgeDataFlag = 1;
    
    // Sets flag for outputting edge data
    if ((out == NULL) || (num_taxa <= 1)) {
        outputEdgeDataFlag = 0;
    }
    
    // Initialize active_node_map 
    for (int i = 0; i < num_active_nodes; i++) { // Assigns indices (0, 1... num_taxa) to leaf nodes
        *nodeMapPtr = i;
        nodeMapPtr++;
    }
    
    // Neighbor Joining Algorithm
    while (num_active_nodes > 2) { // Runs until two nodes remain
    
        // Calculates row_sums
        for (int i = 0; i < num_active_nodes; i++) { // Iterates through ROWs of active nodes in distances
            indexI = *nodeMapI;
            for (int j = 0; j < num_active_nodes; j++) { // Iterates through COLs of active nodes in distances
                indexJ = *nodeMapJ;
                rowTemp += accessIndices(indexI, indexJ);
                nodeMapJ++;
            }
            *rowSumPtr = rowTemp;
            rowSumPtr++;
            nodeMapJ = active_node_map;
            nodeMapI++;
            rowTemp = 0;
        }
        nodeMapI = active_node_map;
        nodeMapJ = active_node_map;
        rowSumPtr = row_sums;
        
        // Calculates and determines node pair to be joined
        for (int i = 0; i < num_active_nodes; i++) { // Iterates through ROWs of active nodes in distances
            indexI = *nodeMapI;
            for (int j = 0; j < num_active_nodes; j++) { // Iterates through COLs of active nodes in distances
                indexJ = *nodeMapJ;
                if (indexI == indexJ) { // Skip main diagonal
                    nodeMapJ++;
                    continue;
                }
                qTemp = qIJ(indexI, indexJ, i, j);
                if (initializeMinPairFlag == 1) { // Initialize minPair, distances indices, active_node_map indices for current iteration
                    minPair = qTemp;
                    indexF = indexI;
                    activeNodeI = i;
                    indexG = indexJ;
                    activeNodeJ = j;
                    initializeMinPairFlag = 0;
                } else if (qTemp < minPair) { // If computed qIJ is less than minPair, reassign minPair, distances indices, active_node_map indices
                    minPair = qTemp;
                    indexF = indexI;
                    activeNodeI = i;
                    indexG = indexJ;
                    activeNodeJ = j;
                }
                nodeMapJ++;
            }
            nodeMapJ = active_node_map;
            nodeMapI++;
        }
        nodeMapI = active_node_map;
        nodeMapJ = active_node_map;
        initializeMinPairFlag = 1; // Reset initialization flag
        
        // Creates new node in nodes and creates name in node_names
        writeNodeName();
        
        // Calculate distances for new node
        indexU = num_all_nodes;
        nodeMapPtr = active_node_map;
        for (int i = 0; i < num_active_nodes; i++) { 
            indexK = *nodeMapPtr;
            uCol = *distances + num_all_nodes + (indexK * MAX_NODES);
            uRow = *distances + indexK+ (num_all_nodes * MAX_NODES);
            if (indexK == indexU) { // If current node is new node, assign distance of zero
                distanceVal = 0;
            } else if (indexK == indexF) { // If current node is child of new node, perform distance calculation
                distanceVal = dUK(indexF, indexG, activeNodeI, activeNodeJ);
            } else if (indexK == indexG) { // If current node is child of new node, perform distance calculation
                distanceVal = dUK(indexG, indexF, activeNodeJ, activeNodeI);
            } else { // If current node is not new node or child or new node, perform distance calculation
                distanceVal = ((accessIndices(indexF, indexK) + accessIndices(indexG, indexK) - accessIndices(indexF, indexG)) / 2);
            }
            *uCol = distanceVal;
            *uRow = distanceVal;
            nodeMapPtr++;
        }
        
        // Update node neighbor information
        nodePtrU = nodes + indexU;
        nodePtrF = nodes + indexF;
        nodePtrG = nodes + indexG;
        *(nodePtrU->neighbors + 1) = nodePtrF; // Assign F to the 1st index of U->neighbors (Left child)
        *(nodePtrU->neighbors + 2) = nodePtrG; // Assign G to the 2nd index of U->neighbors (Right child)
        *(nodePtrF->neighbors) = nodePtrU; // Assign U to the 0th index of F->neighbors (Parent)
        *(nodePtrG->neighbors) = nodePtrU; // Assign U to the 0th index of G->neighbors (Parent)
        
        // Output edge data
        if (outputEdgeDataFlag == 1) {
            fprintf(out, "%d,%d,%.2f\n", indexF, indexU, accessIndices(indexF, indexU)); 
            fprintf(out, "%d,%d,%.2f\n", indexG, indexU, accessIndices(indexG, indexU));
        }
        
        // Deactivate nodes
        *(active_node_map + activeNodeI) = indexU; 
        *(active_node_map + activeNodeJ) = *(active_node_map + num_active_nodes - 1);
        num_active_nodes--;
        num_all_nodes++;
    }
    
    // Neighbor join remaining nodes
    nodeMapPtr = active_node_map;
    indexF = *nodeMapPtr;
    indexG = *(nodeMapPtr + 1);
    nodePtrF = nodes + indexF;
    nodePtrG = nodes + indexG;
    *(nodePtrF->neighbors) = nodePtrG; // Assign G to 0th index of F->neighbors
    *(nodePtrG->neighbors) = nodePtrF; // Assign F to 0th index of G->neighbors
    
    // Output edge data for last nodes
    if (outputEdgeDataFlag == 1) {
        if (indexG < indexF) { // If the second node is smaller, print it's index first
            fprintf(out, "%d,%d,%.2f\n", indexG, indexF, accessIndices(indexF, indexG));
        } else {
            fprintf(out, "%d,%d,%.2f\n", indexF, indexG, accessIndices(indexF, indexG));
        }
    }
    
    return 0;
    abort();
}

static void copyString(char *source, char *destination) {
    while (*source != '\0') {
        *destination = *source;
        source++;
        destination++;
    }
    *destination = '\0';
}

static int equalString(char *str1, char *str2) {
    while (*str1 != '\0') {
        if (*str1 != *str2) {
            return 0;
        }
        str1++;
        str2++;
    }
    if (*str1 != *str2) {
        return 0;
    }
    return 1;
}

static int duplicateName(char *str, char *nameArray, int numTaxa) {
    for (int i = 0; i < numTaxa; i++) {
        if (equalString(str, nameArray)) {
            return -1;
        }
        nameArray += (INPUT_MAX+1);
    }
    return 0;
}

static double calculateNum(char *source) {
    int decimalFlag = 0, digitCounter = 0, decimalCounter = 1;
    double integerCalc = 0, decimalCalc = 0;
    while (*source != '\0') {
        if (*source == '.') {
            if ((integerCalc == 0) && (digitCounter != 1)) {
                invalidNum = 1;
                return 0;
            } 
            decimalFlag++;
            source++;
            continue;
        }
        if (decimalFlag == 0) {
            integerCalc *= 10;
            integerCalc += (*source - 48);
            digitCounter++;
        } else {
            decimalCalc *= 10;
            decimalCalc += (*source - 48);
            decimalCounter *= 10;
        }
        source++;
    }
    source--;
    if (*source == '.') {
        invalidNum = 1;
        return 0;
    }
    if (decimalFlag > 0) {
        decimalCalc /= decimalCounter;
        integerCalc += decimalCalc;
    }
    if (integerCalc == 0) {
        if ((decimalFlag == 1) && (digitCounter > 2)) {
            invalidNum = 1;
        } else if ((decimalFlag == 0) && (digitCounter > 1)) {
            invalidNum = 1;
        }
    }
    return integerCalc;
}

static double qIJ(int indexI, int indexJ, int i, int j) {
    double dIJ = accessIndices(indexI, indexJ);
    double *rowSumPtr = row_sums;
    double sI = *(rowSumPtr + i);
    rowSumPtr = row_sums;
    double sJ = *(rowSumPtr + j);
    return (num_active_nodes - 2) * dIJ - sI - sJ;
}

static double dUK(int indexF, int indexG, int f, int g) {
    double dFG = accessIndices(indexF, indexG);
    double *rowSumPtr = row_sums;
    double sF = *(rowSumPtr + f);
    rowSumPtr = row_sums;
    double sG = *(rowSumPtr + g);
    return (dFG + (sF - sG) / (num_active_nodes - 2)) / 2;
}

static double accessIndices(int i, int j) {
    double *distancePtr = *distances;
    distancePtr += i;
    distancePtr += j * MAX_NODES;
    return *distancePtr;
}

static void writeNodeName() {
    char *namePtr = *node_names + (num_all_nodes * (INPUT_MAX+1));
    NODE *nodePtr = nodes + (num_all_nodes);
    int temp = num_all_nodes;
    *namePtr = '#';
    namePtr++;
    for (int i = logn(num_all_nodes); i >= 0; i--) {
        for (int j = 0; j < i; j++) {
            temp /= 10;
        }
        temp %= 10;
        *namePtr = temp + 48;
        namePtr++;
        temp = num_all_nodes;
    }
    namePtr = *node_names + (num_all_nodes * (INPUT_MAX+1));
    nodePtr->name = namePtr;
}

static int logn(int n) {
    if (n > 10 - 1) {
        return 1 + logn(n / 10);
    } else {
        return 0;
    }
}

static void printNewick(NODE *node, NODE *parent, FILE *out) {
    int parentIndex = compareNode(node, parent);
    int leaf = isLeaf(node);
    if (!leaf) {
        if (firstOpenParen) {
            fprintf(out, ",");
            firstOpenParen--;
        }
        fprintf(out, "(");
        firstLeafComma = 1;
        firstOpenParen = 0;
    }
    if ((*(node->neighbors) != NULL) && (parentIndex != 0)) {
        printNewick(*(node->neighbors), node, out);
    }
    if ((*(node->neighbors + 1) != NULL) && (parentIndex != 1)) {
        printNewick(*(node->neighbors + 1), node, out);
    }
    if ((*(node->neighbors + 2) != NULL) && (parentIndex != 2)) {
        printNewick(*(node->neighbors + 2), node, out);
    }
    int indexI = findIndexName(node->name);
    int indexJ = findIndexName(parent->name);
    if (leaf) {
        if (firstLeafComma) {
            fprintf(out, "%s:%.2f", node->name, accessIndices(indexI, indexJ));
            firstLeafComma--;
        } else {
            fprintf(out, ",%s:%.2f", node->name, accessIndices(indexI, indexJ));
        }
    } else {
        fprintf(out, ")%s:%.2f", node->name, accessIndices(indexI, indexJ));
        firstOpenParen++;
    }
}

static int findIndexName(char *inputName) {
    NODE *nodePtr = nodes;
    char *namePtr = nodePtr->name;
    for (int i = 0; i < num_all_nodes; i++) {
        if (equalString(namePtr, inputName)) {
            return i;
        }
        nodePtr++;
        namePtr = nodePtr->name;
    }
    return -1;
}

static int findOutlier() {
    NODE *nodePtr = nodes;
    double *distancePtr = *distances;
    double rowSum = 0;
    double rowMax = 0;
    int rowIndex = 0;
    for (int i = 0; i < num_taxa; i++) {
        for (int j = 0; j < num_taxa; j++) {
            rowSum += *distancePtr;
            distancePtr++;
        }
        if (rowSum > rowMax) {
            rowMax = rowSum;
            rowIndex = i;
        }
        distancePtr += (MAX_NODES - num_taxa);
        rowSum = 0;
    }
    return rowIndex;
}

static int compareNode(NODE *node, NODE *parent) {
    char *parentName = parent->name;
    if (((*(node->neighbors) != NULL)) && (equalString(parentName, (*(node->neighbors))->name))) {
        return 0;
    } else if (((*(node->neighbors + 1) != NULL)) && (equalString(parentName, (*(node->neighbors + 1))->name))) {
        return 1;
    } else if (((*(node->neighbors + 2) != NULL)) && (equalString(parentName, (*(node->neighbors + 2))->name))) {
        return 2;
    }
    return -1;
}

static int isLeaf(NODE *node) {
    int numNeighbors = 0;
    if (*(node->neighbors) != NULL) {
        numNeighbors++;
    }
    if (*(node->neighbors + 1) != NULL) {
        numNeighbors++;
    }
    if (*(node->neighbors + 2) != NULL) {
        numNeighbors++;
    }
    if (numNeighbors == 1) {
        return 1;
    }
    return 0;
}
