#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>

#define COMMW MPI_COMM_WORLD

int main(int argc, char* argv[]){
	int job_index = 0;
	int requset = 0;
	int rank = 0;
	int size = 0;
	int root = 0;
	int NUMBER_OF_JOBS = 17;
	int JOB_REQUEST = 123;
	int NO_JOBS = -1;
	int * jobArray = NULL;
	MPI_Status status;
	

	MPI_Init(&argc, &argv);
    MPI_Comm_size(COMMW, &size);
    MPI_Comm_rank(COMMW, &rank);

    if(root == rank) {
    	jobArray = (int *) malloc(sizeof(int) * NUMBER_OF_JOBS);

    	for(int i = 0; i < NUMBER_OF_JOBS; i++) {
    		jobArray[i] = (0 == i % 2) ? 1 + i * 2 : abs(30 - i * 3);
    	}

    	while(job_index < NUMBER_OF_JOBS) {
    		MPI_Recv(&requset, 1, MPI_INT, MPI_ANY_SOURCE, 169, COMMW, &status);
    		if(JOB_REQUEST == requset) {
    			MPI_Send(&jobArray[job_index], 1, MPI_INT, status.MPI_SOURCE, 169, COMMW);
    			//MPI_Send(&job_index, 1, MPI_INT, status.MPI_SOURCE, 169, COMMW);
    			job_index++;
    		}
    	}
    	for(int i = 1; i < size; i++) {
    		MPI_Recv(&requset, 1, MPI_INT, MPI_ANY_SOURCE, 169, COMMW, &status);
    		if(JOB_REQUEST == requset) {
    			MPI_Send(&NO_JOBS, 1, MPI_INT, status.MPI_SOURCE, 169, COMMW);
    		}	
    	}
    	free(jobArray);
    }
    if(root != rank) {
    	while(1) {
	    	MPI_Send(&JOB_REQUEST, 1, MPI_INT, root, 169, COMMW);
	    	MPI_Recv(&requset, 1, MPI_INT, root, 169, COMMW, &status);
	    	if(NO_JOBS != requset) {
	    		printf("[%d] is doing {%d} job\n", rank, requset);
	    		sleep(requset);
	    	} else {
	    		printf("All jobs are done.[%d]\n", rank);
	    		break;
	    	}
	    }
    }

    MPI_Finalize();

return 0;
}