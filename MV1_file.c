#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define COMMW MPI_COMM_WORLD
#define DBG 1
#define HEIGHT 3
#define WIDTH 3
#define CHECK 0

int main(int argc, char* argv[]){

    FILE *file;
    int size, rank, i, j, iter, k = 0;
    int root = 0, start = 0, shift = 0;
    int div_h = 0, mod_h = 0, div_w = 0, mod_w = 0;
    double **partOfMatr = NULL, **matr = NULL;
    double *vector = NULL, *partOfVector = NULL;
    double *loc_res = NULL, *raw_res = NULL, *result = NULL;
    int *loc_hghts = NULL, *starts_v = NULL, *starts = NULL, *loc_wdths = NULL;

    MPI_Status status;
							//allocating memory for a source matrix and vector
    matr = (double **) malloc(sizeof(double*) * HEIGHT);
    vector = (double *) malloc(sizeof(double) * WIDTH);
    result = (double *) malloc(sizeof(double*) * HEIGHT);


    if((NULL == matr) || (NULL == vector) || (NULL == result)){
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
    for(i = 0; i < HEIGHT; i++){
	for(j = 0; j < WIDTH; j++){
		fscanf(file, "%lf", &matr[i][j]);		
	}
    }
	for(j = 0; j < WIDTH; j++){
		fscanf(file, "%lf", &vector[j]);
	}
    fclose(file);
//---------------------------------------parallel computation
    MPI_Init(&argc, &argv);
    MPI_Comm_size(COMMW, &size);
    MPI_Comm_rank(COMMW, &rank);

    loc_wdths = (int *) malloc(sizeof(int) * size); //widths of the vector's parts and vertical blocks in the matrix
    loc_hghts = (int *) malloc(sizeof(int) * size); //heights of matrx's parts
    starts = (int *) malloc(sizeof(int) * (size + 1)); //each process starts from its own position
    starts_v = (int *) malloc(sizeof(int) * (size + 1)); //the same navigation in vertical


    if((NULL == loc_wdths) || (NULL == loc_hghts) || (NULL == starts)){
        fprintf(stderr, "memory  wasn't allocated\n");
        exit(1);
      }

    div_h = HEIGHT / size;
    mod_h = HEIGHT % size;
    div_w = WIDTH / size;
    mod_w = WIDTH % size;
    #if DBG 
	if(CHECK == rank) printf("div_h mod_h div_w mod_w\n"); 
	if(CHECK == rank) printf("%d     %d     %d     %d\n", div_h, mod_h, div_w, mod_w); 
	if(CHECK == rank) printf("l_w l_h st st_v\n");
    #endif

    for(i = 0; i < size; i++){
	loc_wdths[i] = (i < mod_w) ? div_w + 1 : div_w;	//alloting WIDTH of the matrix to [size] of processes
	loc_hghts[i] = (i < mod_h) ? div_h + 1 : div_h; //first mod_h take +1 more row from the matrix
	starts[i] = (i <= mod_w) ? (div_w + 1) * i : (div_w + 1) * mod_w + div_w * (i - mod_w); // next process starts righter from previous
	starts_v[i] = (i <= mod_h) ? (div_h + 1) * i : (div_h + 1) * mod_h + div_h * (i - mod_h); // next process starts righter from previous
    #if DBG 
	if(CHECK == rank) printf("%d   %d   %d   %d\n", loc_wdths[i], loc_hghts[i], starts[i], starts_v[i]); 
    #endif
    }						//first mod_w shifts on steps (div_w + 1), the rest of other steps == (div_w + 1) and (div_w)
    starts[size] = WIDTH; //right range
    starts_v[size] = HEIGHT; //bottom range

    partOfMatr = (double **) malloc(sizeof(double*) * loc_hghts[rank]);
    partOfVector = (double *) malloc(sizeof(double) * (div_w + 1)); //allocating memory for working parts of matrix and vector
    loc_res = (double *) malloc(sizeof(double) * (div_h + 1));	//and local result: size (div_h + 1) for any prcess (+1 - reserve)
    raw_res = (double *) malloc(sizeof(double) * (div_h + 1) * size); //length more than HEIGHT as gathering demands common size for any buf
								//so we take wich reserving (div_h + 1) for each process

    if((NULL == partOfMatr) || (NULL == partOfVector) || (NULL == loc_res) || (NULL == raw_res)){
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
					//fill numbers from the source matrix and vector
    for(i = 0; i < loc_hghts[rank]; i++){
	for(j = 0; j < WIDTH; j++){
	    partOfMatr[i][j] = matr[i + starts_v[rank]][j];
    #if DBG 
	if(CHECK == rank) printf("%.0f ", partOfMatr[i][j]); 
    #endif
	}
    #if DBG 
	if(CHECK == rank) printf(" part of matrix\n"); 
    #endif
    }
    for(i = 0; i < loc_wdths[rank]; i++){
	    partOfVector[i] = vector[i + starts[rank]];
    }

    for(iter = 0; iter < size; iter++){	//number of iterations - [size]
        for(i = 0; i < loc_hghts[rank]; i++){
	    start = starts[(rank + iter) % size];//production in proper place
	    j = start;
	    while(j < start + loc_wdths[(rank + iter) % size]){ //j should be less than next blok starts
		partOfMatr[i][j] *= partOfVector[j - start];
    #if DBG 
	if(CHECK == rank) printf("%.0f ", partOfMatr[i][j]); 
    #endif
		j++;
	    }
    #if DBG 
	if(CHECK == rank) printf(" part of producted matrix\n"); 
    #endif
        }

	MPI_Sendrecv_replace(partOfVector, (div_w + 1), MPI_DOUBLE, //sending current part of vector, getting next part
			    (0 == rank) ? size - 1 : rank - 1, 1, (rank + 1) % size, 1, //all local vectors have (div_h + 1) length
			    COMMW, &status);
    #if DBG 
	if(CHECK == rank){
	    printf("part of vector after send recieve %d: ", iter);
	    for(i = 0; i < div_h + 1; i++){
		 printf("%.0f ", partOfVector[i]); 
	    }
	    printf("\n");
	}
    #endif
    }

    memset(loc_res, 0, sizeof(double) * (div_h + 1));
    for(i = 0; i < loc_hghts[rank]; i++){
	for(j = 0; j < WIDTH; j ++){
	    loc_res[i] += partOfMatr[i][j]; //summing up for local result
	}
    #if DBG 
	if(CHECK == rank) printf("%.0f local result\n", loc_res[i]); 
    #endif
    }
    
    MPI_Gather(loc_res, (div_h + 1), MPI_DOUBLE, 
		raw_res, (div_h + 1), MPI_DOUBLE, 
		root, COMMW); // raw result will have (div_h + 1) * size elements - more than HEIGHT
    
    if(root == rank){
	shift = 0;
	for(i = 0; i < size; i++){ //getting pure result from special places from raw result
	    memcpy((result + shift), (raw_res + i * (div_h + 1)), sizeof(double) * loc_hghts[i]);
	    shift += loc_hghts[i];
	}
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
    free(partOfVector);

    free(loc_wdths);
    free(loc_hghts);
    free(loc_res);
    free(raw_res);
    free(starts);
    free(starts_v);

    MPI_Finalize();

    for(i = 0; i < HEIGHT; i++){
      free(matr[i]);
    }
    free(matr);
    free(vector);
    free(result);

return 0;
}