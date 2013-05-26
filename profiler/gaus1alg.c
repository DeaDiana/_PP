#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define COMMW MPI_COMM_WORLD
#define DBG 0
#define HEIGHT 9
#define WIDTH 9
#define CHECK 1

int main(int argc, char* argv[]){

    FILE *file;
    int size, rank, i, j, iter, h = 0, k = 0, shift = 0, sub = 0, d = 0;
    int root = 0;
    int div_h = 0, mod_h = 0;
    double **partOfMatr = NULL, **matr = NULL;
    double *row = NULL;
    double *loc_res = NULL, *result = NULL;
    int *loc_hghts = NULL, *starts_v = NULL;
    double ijCoef = 0.0, iiCoef = 0.0l;
    MPI_Status status;
							//allocating memory for a source matrix and vector
    matr = (double **) malloc(sizeof(double*) * HEIGHT);
    result = (double *) malloc(sizeof(double*) * HEIGHT);
    row = (double *) malloc(sizeof(double) * (WIDTH + 1));

    if((NULL == matr) || (NULL == result) || (NULL == row)){
        fprintf(stderr, "memory  wasn't allocated\n");
        exit(1);
      }

    for(i = 0; i < HEIGHT; i++){
      matr[i] = (double*) malloc(sizeof(double) * (WIDTH + 1));
      if(NULL == matr[i]){
        fprintf(stderr, "memory  wasn't allocated\n");
        exit(1);
      }
    }
					
    if(NULL == (file = fopen("matrix2", "rb"))) {
       fprintf(stderr, "file with source matrix wasn't opened\n");
        exit(1); 
    }
                    //fill the source matrix and vector with numbers
    for(i = 0; i < HEIGHT; i++){
        for(j = 0; j <= WIDTH; j++){
	    fscanf(file, "%lf", &matr[i][j]);
	}
    }
    fclose(file);

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
      partOfMatr[i] = (double*) malloc(sizeof(double) * (WIDTH + 1));
      if(NULL == partOfMatr[i]){
        fprintf(stderr, "memory  wasn't allocated\n");
        exit(1);
      }
    }
					//fill numbers from the source matrix
    for(i = 0; i < loc_hghts[rank]; i++){
	for(j = 0; j <= WIDTH; j++){
	    partOfMatr[i][j] = matr[i + starts_v[rank]][j];
	}
    }

//gauss algorithm
    //1st pass																					
	//receive row
    
    for(iter = 0; iter < rank; iter++){
        for(sub = 0; sub < loc_hghts[iter]; sub++){
	    MPI_Recv(row, WIDTH + 1, MPI_DOUBLE, iter, 169, COMMW, &status);
	    #if DBG
	    if(CHECK == rank){
		for(d = 0; d <= WIDTH; d++){
		printf("%.1f ", row[d]);
		}
		printf("\n");
	    }
	    #endif
				// * and -
	    for(i = 0; i < loc_hghts[rank]; i++){
		ijCoef = partOfMatr[i][sub + starts_v[iter]];		//think better here!
		#if DBG
		//if(CHECK == rank) printf("%.1f ijCoef\n", ijCoef);
		#endif
		for(j = 0; j <= WIDTH; j++){
		    #if DBG
		    //if(CHECK == rank) printf("%.1f -= %.1f * %.1f:\n ", partOfMatr[i][j], ijCoef, row[j]);
		    #endif
		    partOfMatr[i][j] -= ijCoef * row[j];
		    #if DBG
		    //if(CHECK == rank) printf("%.1f ", partOfMatr[i][j]);
		    #endif
		}
		#if DBG
		//if(CHECK == rank) printf("calculating %d\n", i);
		#endif
	    }
	    shift++;
	}
    }

	shift = 0;
			//divide on diagonal element
    for(i = 0; i < loc_hghts[rank]; i++){
        iiCoef = partOfMatr[i][i + starts_v[rank]];
        for(j = i; j <= WIDTH; j++){
	    if(iiCoef != 0)
	    partOfMatr[i][j] /= (double) iiCoef;
	}
	    #if DBG
	    if(CHECK == rank){
	    for(d = 0; d <= WIDTH; d++){
	    printf("%.1f SEND %d\n", partOfMatr[i][d], i);
	    }
	    }
	    #endif
			//send for each process below
	for(j = (rank + 1); j < size; j++){
	    MPI_Send(partOfMatr[i], WIDTH + 1, MPI_DOUBLE, j, 169, COMMW);
	}
			// * and - to all row in this block
	for(h = i+1; h < loc_hghts[rank]; h++){
	    ijCoef = partOfMatr[h][i + starts_v[rank]];		//think better here
	    for(k = 0; k <= WIDTH; k++){
		partOfMatr[h][k] -= ijCoef * partOfMatr[i][k];
		#if DBG
		if(CHECK  == rank) printf("%.1f ", partOfMatr[h][k]);
		#endif
	    }
		#if DBG
		if(CHECK  == rank) printf("\n-------------\n");
		#endif
	}
	shift++;
    }

    //back pass
	for(iter = size - 1; iter > rank; iter--){ //more iterartions indeed
	    for(sub = (loc_hghts[iter] - 1); sub >= 0; sub--){
	        MPI_Recv(row, WIDTH + 1, MPI_DOUBLE, iter, 196, COMMW, &status);
	    #if DBG
	    if (0 == rank) printf("received\n");
	    #endif
		shift = 0;
				// * and - bottom row
		for(i = 0; i < loc_hghts[rank]; i++){
		    ijCoef = partOfMatr[i][sub + starts_v[iter]];
		    for(j = 0; j <= WIDTH; j++){	//could be optimized
			partOfMatr[i][j] -= ijCoef * row[j];
		    }
		}
		    shift++;
	    }
	}
		//works in current block
	for(i = (loc_hghts[rank] - 1); i >= 0; i--){
		//send row for all other above processes
	    for(j = rank - 1; j >= 0; j--){
		MPI_Send(partOfMatr[i], WIDTH + 1, MPI_DOUBLE, j, 196, COMMW);
	        #if DBG
		if (3 == rank) printf("sended\n");
		#endif
	    }
	    shift = 0;
	    for(h = 0; h < i; h++){
		ijCoef = partOfMatr[h][i + starts_v[rank]];
		for(k = 0; k <= WIDTH; k++){
		    partOfMatr[h][k] -= ijCoef * partOfMatr[i][k];
		}
	    }
	        #if DBG
		if (3 == rank) printf("calculated and readry for next send\n");
		#endif
		shift++;
	}

    memset(loc_res, 0, sizeof(double) * (div_h + 1)); //length for local result with reserving 1 element
    for(i = 0; i < loc_hghts[rank]; i++){	
	loc_res[i] = partOfMatr[i][WIDTH]; //result is in last column
    }
    
    MPI_Allgatherv(loc_res, loc_hghts[rank], MPI_DOUBLE, //send particular loc_hghts[rank] - sometime less than div_h + 1
		    result, loc_hghts, starts_v, // loc_hghts contains lenghts of local results, starts contains their positions in result vector
		    MPI_DOUBLE, COMMW); // result will have exactly HEIGHT elements
    
    if(root == rank){
	for(i = 0; i < HEIGHT; i++){
	    printf("%.0f ", result[i]);
	}
	printf("the result is above\n");
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
    free(result);
    free(row);

return 0;
}