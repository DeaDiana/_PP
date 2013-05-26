#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define COMMW MPI_COMM_WORLD
#define HEIGHT 9
#define WIDTH 21
#define EPS 0.00001
#define CHECK 1
#define DBG 0
#define PRINT_RESULT 1
#define PRINT_ORIGIN 1

int main(int argc, char *argv[])
{
	double ** matr = NULL;
	double ** current_partOfMatr = NULL;
	double ** last_partOfMatr = NULL;
	double ** swap = NULL;
	int * loc_wdths = NULL;
	int * starts_h = NULL;
	int rank = 0;
	int size = 0;
	FILE * file;
	MPI_Status status;
	double hx_step = 1.0;
	double hy_step = 1.0;
	double local_max = 0.0;
	double global_max = 0.0;
	double local_differ = 0.0;
	MPI_Request request_send_right;
	MPI_Request request_recv_right;
	MPI_Request request_send_left;
	MPI_Request request_recv_left;
	int root = 0;

	//allocating memory for source matrix
	matr = (double **) malloc(sizeof(double*) * HEIGHT);
	if(NULL == matr) {
        fprintf(stderr, "memory  wasn't allocated\n");
        exit(1);
    }
    for(int i = 0; i < HEIGHT; i++){
      matr[i] = (double *) malloc(sizeof(double) * WIDTH);
      if(NULL == matr[i]) {
      	fprintf(stderr, "memory  wasn't allocated\n");
        exit(1);
      }
    }

    if(NULL == (file = fopen("matrix.txt", "rb"))) {
       fprintf(stderr, "file with source matrix wasn't opened\n");
        exit(1); 
    }
	//fill the source matrix and vector with numbers
    for(int i = 0; i < HEIGHT; i++){
        for(int j = 0; j < WIDTH; j++){
	    	fscanf(file, "%lf", &matr[i][j]);
		}
    }
    fclose(file);

//---------------------------------------parallel computation
    MPI_Init(&argc, &argv);
    MPI_Comm_size(COMMW, &size);
    MPI_Comm_rank(COMMW, &rank);

    //widths of matrx's parts and navigation in horizontal
    loc_wdths = (int *) malloc(sizeof(int) * size); 
    starts_h = (int *) malloc(sizeof(int) * (size + 1));

    if((NULL == loc_wdths) || (NULL == starts_h)) {
        fprintf(stderr, "memory  wasn't allocated\n");
        exit(1);
    }

    int div_w = WIDTH / size;
    int mod_w = WIDTH % size;

    //fill navigation information
    for(int i = 0; i < size; i++){
		loc_wdths[i] = (i < mod_w) ? div_w + 1 : div_w; //first mod_w take +1 more row from the matrix
		starts_h[i] = (i <= mod_w) ? (div_w + 1) * i : (div_w + 1) * mod_w + div_w * (i - mod_w); // next process starts righter from previous
		if((0 == i) || ((size - 1) == i))
		{
			loc_wdths[i]++;
			if((size - 1) == i) 
			{
				starts_h[i]--;
			}
		} else {
			loc_wdths[i] += 2;
			starts_h[i]--;
		}
    }						//first mod_w shifts on steps (div_w + 1), the rest of other steps == (div_w + 1) and (div_w)
    starts_h[size] = WIDTH; //bottom range

    //allocating memory for working parts of matrix
    //algorthm demands 2 matrix to work - old and current
    current_partOfMatr = (double **) malloc(sizeof(double*) * HEIGHT);
    last_partOfMatr = (double **) malloc(sizeof(double*) * HEIGHT);

    if((NULL == current_partOfMatr) || (NULL == last_partOfMatr)) {
        fprintf(stderr, "memory  wasn't allocated\n");
        exit(1);
    }

    for(int i = 0; i < HEIGHT; i++) {
      current_partOfMatr[i] = (double *) malloc(sizeof(double) * loc_wdths[rank]);
      last_partOfMatr[i] = (double *) malloc(sizeof(double) * loc_wdths[rank]);
      if((NULL == current_partOfMatr[i]) || (NULL == last_partOfMatr[i])) {
        fprintf(stderr, "memory  wasn't allocated\n");
        exit(1);
      }
    }
	//fill numbers from the source matrix
    for(int i = 0; i < HEIGHT; i++){
		for(int j = 0; j < loc_wdths[rank]; j++){
		    current_partOfMatr[i][j] = matr[i][j + starts_h[rank]];
		    last_partOfMatr[i][j] = matr[i][j + starts_h[rank]];
		}
    }

   	#if PRINT_ORIGIN
		if(root == rank)
		{
			for(int i = 0; i < HEIGHT; i++)
			{
				for(int j = 0; j < WIDTH; j++)
				{
					printf("%.0f ", matr[i][j]);
				}
				printf("\n");
			}
		}
	#endif

    //grid-algorithm started
    double * edge_left_column_to_send = (double *) malloc(sizeof(double) * HEIGHT);
    double * ear_left_column_to_recv = (double *) malloc(sizeof(double) * HEIGHT);
    double * edge_right_column_to_send = (double *) malloc(sizeof(double) * HEIGHT);
    double * ear_right_column_to_recv = (double *) malloc(sizeof(double) * HEIGHT);
    memset(ear_left_column_to_recv, 0, sizeof(double) * HEIGHT);
    memset(ear_right_column_to_recv, 0, sizeof(double) * HEIGHT);
    //prepearings for 1st iteration
    for(int i = 0; i < HEIGHT; i++)
    {
    	edge_left_column_to_send[i] = current_partOfMatr[i][1];
    	edge_right_column_to_send[i] = current_partOfMatr[i][loc_wdths[rank] - 2];
    }
    if((size - 1) != rank)
    {
    	MPI_Isend(edge_right_column_to_send, (HEIGHT - 0), MPI_DOUBLE, rank + 1, 169, COMMW, &request_send_right);
		MPI_Irecv(ear_right_column_to_recv, (HEIGHT - 0), MPI_DOUBLE, rank + 1, 170, COMMW, &request_recv_right);
    }
    if(0 != rank)
    {
    	MPI_Isend(edge_left_column_to_send, (HEIGHT - 0), MPI_DOUBLE, rank - 1, 170, COMMW, &request_send_left);
		MPI_Irecv(ear_left_column_to_recv, (HEIGHT - 0), MPI_DOUBLE, rank - 1, 169, COMMW, &request_recv_left);
    }
    do
	{
	    //update right ear
		if((size - 1) != rank)
	    {
	    	MPI_Wait(&request_send_right, &status);
	    	MPI_Wait(&request_recv_right, &status);
	    	for(int i = 0; i < HEIGHT; i++)
	    	{
	    		current_partOfMatr[i][loc_wdths[rank] - 1] = ear_right_column_to_recv[i];
	    		last_partOfMatr[i][loc_wdths[rank] - 1] = ear_right_column_to_recv[i];
	    	}
	    }
	    //update left ear
	    if(0 != rank)
	    {
	    	MPI_Wait(&request_send_left, &status);
	    	MPI_Wait(&request_recv_left, &status);
	    	for(int i = 0; i < HEIGHT; i++)
	    	{
	    		current_partOfMatr[i][0] = ear_left_column_to_recv[i];
	    		last_partOfMatr[i][0] = ear_left_column_to_recv[i];
	    	}
	    }
	    #if DBG
		    if(CHECK == rank)
			{
		        for(int i = 0; i < (HEIGHT - 0); i++)
			    {
		    		printf("%.1f ", ear_right_column_to_recv[i]);
		    	}
		    	printf("\n");
		    }
		#endif
		local_max = 0.0;
		for(int i = 1; i < (HEIGHT - 1); i++)
		{
			// as column (WIDTH - 1) should be received from another process
			for(int j = 1; j < (loc_wdths[rank] - 1); j++) 
			{
				current_partOfMatr[i][j] = (1 / (2 / (hx_step * hx_step) + 2 / (hy_step * hy_step))) * ((current_partOfMatr[i][j - 1] + last_partOfMatr[i][j + 1]) / (hx_step * hx_step) + (current_partOfMatr[i - 1][j] + last_partOfMatr[i + 1][j]) / (hy_step * hy_step));
				//current_partOfMatr[i][j] = 0.25 * ((current_partOfMatr[i][j - 1] + last_partOfMatr[i][j + 1]) / (hx_step * hx_step) + (current_partOfMatr[i - 1][j] + last_partOfMatr[i + 1][j]) / (hy_step * hy_step));
				local_differ = abs(current_partOfMatr[i][j] - last_partOfMatr[i][j]);
				if(local_differ > local_max)
				{
					local_max = local_differ;
				}
				#if DBG
					if(CHECK == rank)
					{
						printf("%.1f ", current_partOfMatr[i][j]);
					}
				#endif
			}
			#if DBG
				if(CHECK == rank)
				{
					printf("\n");
				}
			#endif
		}
		//set edges to exchage with neighbours
		for(int i = 0; i < HEIGHT; i++)
	    {
	    	edge_left_column_to_send[i] = current_partOfMatr[i][1];
	    	edge_right_column_to_send[i] = current_partOfMatr[i][loc_wdths[rank] - 2];
	    }
	    #if DBG
		    if(CHECK == rank)
		    {
		    	//printf("before async send and recv [%d]\n", rank);
		    	for(int i = 0; i < HEIGHT; i++)
				{
					for(int j = 0; j < loc_wdths[rank]; j++)
					{
						printf("%.0f ", current_partOfMatr[i][j]);
					}
					printf("  |%d   ", rank);

					for(int j = 0; j < loc_wdths[rank]; j++)
					{
						printf("%.0f ", last_partOfMatr[i][j]);
					}
					printf("  |%d\n", rank);
				}
		    }
	    #endif
	    // reduce norma of whole matrix and check loop condition to break
	    MPI_Allreduce (&local_max, &global_max, 1, MPI_DOUBLE, MPI_MAX, COMMW);
	    if(global_max < EPS)
	    {
	    	if(CHECK == rank){	printf("%lf global_max, %lf EPS\n", global_max, EPS);	}
	    	//gathering result matrix
	    	for(int i = 0; i < HEIGHT; i++)
	    	{
   				MPI_Gatherv(current_partOfMatr[i], loc_wdths[rank], MPI_DOUBLE, matr[i], loc_wdths, starts_h, MPI_DOUBLE, root, COMMW);
   			}
	    	break;
	    }
	    //async send own edges and async receiving neighbour edges
	    if((size - 1) != rank)
	    {
	    	MPI_Isend(edge_right_column_to_send, (HEIGHT - 0), MPI_DOUBLE, rank + 1, 169, COMMW, &request_send_right);
			MPI_Irecv(ear_right_column_to_recv, (HEIGHT - 0), MPI_DOUBLE, rank + 1, 170, COMMW, &request_recv_right);
	    }
	    if(0 != rank)
	    {
	    	MPI_Isend(edge_left_column_to_send, (HEIGHT - 0), MPI_DOUBLE, rank - 1, 170, COMMW, &request_send_left);
			MPI_Irecv(ear_left_column_to_recv, (HEIGHT - 0), MPI_DOUBLE, rank - 1, 169, COMMW, &request_recv_left);
	    }
	    #if DBG
	    	if(CHECK == rank) {	printf("%lf, %lf\n", local_max, global_max);	}
	    #endif
	    //for next iteration current part of matrix becomes old 
	    //and new part will be constructed in last_partOfMatr
		swap = current_partOfMatr;
		current_partOfMatr = last_partOfMatr;
		last_partOfMatr = swap;
	} while(1);
    //grid-algorithm finished

	#if PRINT_RESULT
		if(root == rank)
		{
			for(int i = 0; i < HEIGHT; i++)
			{
				for(int j = 0; j < WIDTH; j++)
				{
					printf("%.0f ", matr[i][j]);
				}
				printf("\n");
			}
		}
	#endif

    MPI_Finalize();

	//free all allocated memory
	for(int i = 0; i < HEIGHT; i++){
    	free(matr[i]);
    	free(current_partOfMatr[i]);
      	free(last_partOfMatr[i]);
    }
    free(matr);
    free(current_partOfMatr);
    free(last_partOfMatr);
    free(loc_wdths);
    free(starts_h);
    free(edge_left_column_to_send);
    free(ear_left_column_to_recv);
    free(edge_right_column_to_send);
    free(ear_right_column_to_recv);

	return 0;
}