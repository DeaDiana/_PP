#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define COMMW MPI_COMM_WORLD
#define DBG 0
#define HEIGHT 12
#define WIDTH 18
//#define HEIGHT 800
//#define WIDTH 1000
#define CHECK 1

int main(int argc, char* argv[]){

    int size, rank, i, j, iter;
    int root = 0;
    int div_h = 0, mod_h = 0;
    double **partOfMatr = NULL, **matr = NULL;
    double *vector = NULL;
    double *loc_res = NULL, *result = NULL;
    int *loc_hghts = NULL, *starts_v = NULL;

    MPI_Status status;
							//allocating memory for a source matrix and vector
    matr = (double **) malloc(sizeof(double*) * HEIGHT);
    vector = (double *) malloc(sizeof(double) * WIDTH);
    result = (double *) malloc(sizeof(double*) * HEIGHT);


    if((NULL == matr) || (NULL == vector)){
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
					//fill the source matrix and vector with numbers
    for(i = 0; i < HEIGHT; i++){
	for(j = 0; j < WIDTH; j++){
	    matr[i][j] = (i == j) ? 2.0 : 1.0;
	    vector[j] = 1.0;
	}
    }
//---------------------------------------parallel computation
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

    for(i = 0; i < size; i++){
	loc_hghts[i] = (i < mod_h) ? div_h + 1 : div_h; //first mod_h take +1 more row from the matrix
	starts_v[i] = (i <= mod_h) ? (div_h + 1) * i : (div_h + 1) * mod_h + div_h * (i - mod_h); // next process starts righter from previous
    }						//first mod_w shifts on steps (div_w + 1), the rest of other steps == (div_w + 1) and (div_w)
    starts_v[size] = HEIGHT; //bottom range

    partOfMatr = (double **) malloc(sizeof(double*) * loc_hghts[rank]);//allocating memory for working parts of matrix
    loc_res = (double *) malloc(sizeof(double) * (div_h + 1));	//and local result: size (div_h + 1) for any prcess (+1 - reserve)
								//so we take wich reserving (div_h + 1) for each process

    if((NULL == partOfMatr) || (NULL == loc_res)){
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
					//fill numbers from the source matrix
    for(i = 0; i < loc_hghts[rank]; i++){
	for(j = 0; j < WIDTH; j++){
	    partOfMatr[i][j] = matr[i + starts_v[rank]][j];
	}
    }
					//production part of matrix on whole vector
    for(i = 0; i < loc_hghts[rank]; i++){
	for(j = 0; j < WIDTH; j++){
	    partOfMatr[i][j] *= vector[j];
	}
    }

    memset(loc_res, 0, sizeof(double) * (div_h + 1)); //length for local result with reserving 1 element
    for(i = 0; i < loc_hghts[rank]; i++){	//but fill only sense elements
	for(j = 0; j < WIDTH; j ++){
	    loc_res[i] += partOfMatr[i][j]; //summing up for local result
	}
    }
    
    MPI_Allgatherv(loc_res, loc_hghts[rank], MPI_DOUBLE, //send particular loc_hghts[rank] - sometime less than div_h + 1
		    result, loc_hghts, starts_v, // loc_hghts contains lenghts of local results, starts contains their positions in result vector
		    MPI_DOUBLE, COMMW); // result will have exactly HEIGHT elements
    
    if(root == rank){
	for(i = 0; i < HEIGHT; i++){
	    printf("%.0f ", result[i]);
	}
	printf("\n");
    }
				//free all allocated memory
    for(i = 0; i < loc_hghts[rank]; i++){
      free(partOfMatr[i]);
    }
    free(partOfMatr);

    free(loc_hghts);
    free(starts_v);
    free(loc_res);

    MPI_Finalize();

    for(i = 0; i < HEIGHT; i++){
      free(matr[i]);
    }
    free(matr);
    free(vector);
    free(result);

return 0;
}