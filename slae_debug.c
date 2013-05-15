#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define COMMW MPI_COMM_WORLD
#define HEIGHT 10
#define WIDTH 10
#define EPS 0.01
#define T 0.01
#define DBG 0
#define CHECK 0
#define FROM_FILE 0

double norma(double*);

int main(int argc, char *argv[]){

    FILE *file;
    int size = 0, rank = 0, i = 0, j = 0, iter = 0;
    int root = 0;
    int div_h = 0, mod_h = 0;
    double **partOfMatr = NULL, **matr = NULL;
    double *vector = NULL;
    double *loc_res = NULL, *result = NULL;
    int *loc_hghts = NULL, *starts_v = NULL;
    double *b_vector = NULL, *partOfB_vector = NULL;
    MPI_Status status;
    double prescission = 0.0;

    matr = (double **) malloc(sizeof(double*) * HEIGHT);
    vector = (double *) malloc(sizeof(double) * WIDTH);
    result = (double *) malloc(sizeof(double*) * HEIGHT);
    b_vector = (double *) malloc(sizeof(double*) * WIDTH);

    if((NULL == matr) || (NULL == vector) || (NULL == result) || (NULL == b_vector)){
        fprintf(stderr, "memory  wasn't allocated\n");
        exit(1);
    }

    for(i = 0; i < HEIGHT; i++){
      matr[i] = (double*) malloc(sizeof(double) * WIDTH);
      if(NULL == matr[i]){
        fprintf(stderr, "memory  wasn't allocated\n");
        exit(1);
      }
    }

    if(NULL == (file = fopen("matrix", "rb"))) {
       fprintf(stderr, "file with source matrix wasn't opened\n");
        exit(1); 
    }
                    //fill the source matrix and vector with numbers
    #if FROM_FILE
    for(i = 0; i < HEIGHT; i++){
	for(j = 0; j < WIDTH; j++){
    	    fscanf(file, "%lf", &matr[i][j]);
	}
    }

    for(j = 0; j < WIDTH; j++){
        fscanf(file, "%lf", &b_vector[j]);
	vector[j] = b_vector[j];
    }
    fclose(file);
    #else
    for(i = 0; i < HEIGHT; i++){
	for(j = 0; j < WIDTH; j++){
    	    matr[i][j] = (double) (i == j) ? 2.0 : 1.0;
	}
    	b_vector[i] = (double) (WIDTH + 1);
    	vector[i] = b_vector[i];
    }
    #endif

    prescission = norma(b_vector) * EPS;
    #if DBG
	printf("%f - prescission\n", prescission);
    #endif

    MPI_Init(&argc, &argv);
    MPI_Comm_size(COMMW, &size);
    MPI_Comm_rank(COMMW, &rank);

    loc_hghts = (int *) malloc(sizeof(int) * size); //heights of matrx's parts
    starts_v = (int *) malloc(sizeof(int) * (size + 1)); //the same navigation in vertical

    if((NULL == loc_hghts) || (NULL == starts_v)){
        fprintf(stderr, "memory  wasn't allocated\n");
        exit(1);
      }

    div_h = HEIGHT / size;
    mod_h = HEIGHT % size;
    #if DBG 
	if(CHECK == rank) printf("div_h mod_h\n"); 
	if(CHECK == rank) printf("%d     %d\n", div_h, mod_h); 
    #endif

    for(i = 0; i < size; i++){
    loc_hghts[i] = (i < mod_h) ? div_h + 1 : div_h; //first mod_h take +1 more row from the matrix
    starts_v[i] = (i <= mod_h) ? (div_h + 1) * i : (div_h + 1) * mod_h + div_h * (i - mod_h); // next process starts righter from previous
    }						//first mod_w shifts on steps (div_w + 1), the rest of other steps == (div_w + 1) and (div_w)
    starts_v[size] = HEIGHT; //bottom range

    partOfMatr = (double **) malloc(sizeof(double*) * loc_hghts[rank]);//allocating memory for working parts of matrix
    loc_res = (double *) malloc(sizeof(double) * (div_h + 1));	//and local result: size (div_h + 1) for any prcess (+1 - reserve)
								//so we take wich reserving (div_h + 1) for each process
    partOfB_vector = (double *) malloc(sizeof(double) * loc_hghts[rank]); //particular loc_hghts[rank] elements from source b_vector

    if((NULL == partOfMatr) || (NULL == loc_res) || (NULL == partOfB_vector)){
        fprintf(stderr, "memory  wasn't allocated\n");
        exit(1);
      }

    for(i = 0; i < loc_hghts[rank]; i++){
      partOfMatr[i] = (double*) malloc(sizeof(double) * WIDTH);
      if(NULL == partOfMatr[i]){
        fprintf(stderr, "memory  wasn't allocated\n");
        exit(1);
      }
    }

    while(1){

		    //fill numbers from the source matrix
    for(i = 0; i < loc_hghts[rank]; i++){
	for(j = 0; j < WIDTH; j++){
	    partOfMatr[i][j] = matr[i + starts_v[rank]][j];
	    #if DBG
	    if(CHECK == rank){ printf("%f ", partOfMatr[i][j]); }
	    #endif
	}
	partOfB_vector[i] = b_vector[i + starts_v[rank]]; //corect for square matrixs
	#if DBG
	if(CHECK == rank){ printf("%f ", partOfB_vector[i]); }
        printf("filling the partOfMatr\n");
	#endif
    }

		    //production part of matrix on whole vector
    for(i = 0; i < loc_hghts[rank]; i++){
    for(j = 0; j < WIDTH; j++){
        partOfMatr[i][j] *= vector[j];
	#if DBG 
	    if(CHECK == rank) printf("%f ", partOfMatr[i][j]); 
        #endif
    }
    #if DBG 
        if(CHECK == rank) printf(" part of producted matrix\n"); 
    #endif
    }

    memset(loc_res, 0, sizeof(double) * (div_h + 1)); //length for local result with reserving 1 element
    for(i = 0; i < loc_hghts[rank]; i++){	//but fill only sense elements
	for(j = 0; j < WIDTH; j ++){
    	    loc_res[i] += partOfMatr[i][j]; //summing up for local result
	}
	#if DBG
	printf("Ax[i] - b[i] = %f - %f\n", loc_res[i], partOfB_vector[i]);
	#endif
	loc_res[i] -= partOfB_vector[i];	// Axk - b

    }

    MPI_Allgatherv(loc_res, loc_hghts[rank], MPI_DOUBLE, //send particular loc_hghts[rank] - sometime less than div_h + 1
	    result, loc_hghts, starts_v, // loc_hghts contains lenghts of local results, starts contains their positions in result vector
	    MPI_DOUBLE, COMMW); // result will have exactly HEIGHT elements

    #if DBG 
    if(CHECK == rank){
	for(i = 0; i < WIDTH; i++){	 printf("%f Ax - b\n", result[i]);	}
    }
    #endif


    iter++;
    //chek norma of result <= EPS* || b || to finish algorithm or to continue

	if((norma(result) <= prescission) || (200 < iter)) {
	    
	    break;	
	}
	    else {
	    	for(i = 0; i < WIDTH; i++){
		    vector[i] = vector[i] - T * result[i];
		    #if DBG 
    			if(CHECK == rank) printf("%f for %d\n", vector[i], iter + 1); 
		    #endif
		}
	}
    }

    if(root == rank){
        for(i = 0; i < WIDTH; i++){
	    printf("%f \n", vector[i]);
	}
	printf("the result is above \n");
    }

    for(i = 0; i < loc_hghts[rank]; i++){
      free(partOfMatr[i]);
    }
    free(partOfMatr);

    free(loc_hghts);
    free(starts_v);
    free(loc_res);
    free(partOfB_vector);

    MPI_Finalize();

    for(i = 0; i < HEIGHT; i++){
      free(matr[i]);
    }
    free(matr);
    free(vector);
    free(result);
    free(b_vector);

return 0;
}

double norma(double *v){

    double max = 0.0, compare = 0.0;
    int i = 0;

    for(i = 0; i < HEIGHT; i++){
	compare = (0 < v[i]) ? v[i] : -v[i];
        if(compare > max){
	     max = compare;
	}
    }

return max;
}