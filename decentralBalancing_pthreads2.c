#include <mpi.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>

#define COMMW MPI_COMM_WORLD
#define BUFF_MAX_SIZE 256
#define REST_OF_JOBS 2
#define SHIFT_IN_ARRAY 3
#define TRUE 1

int job_index = 0;
int NUMBER_OF_JOBS = 6;
int * jobArray = NULL;
pthread_mutex_t * mutex = NULL;
pthread_cond_t check_rest = PTHREAD_COND_INITIALIZER;
int JOB_REQUEST = 123;
int NO_JOBS = -1;
int DONE = -2;

void * doJob(void * param) {
    int current_job = 0;
    int done_jobs = 0;
    int rank = jobArray[2];
    while(TRUE) {
        pthread_mutex_lock(mutex);
	        if(job_index < NUMBER_OF_JOBS) {
                current_job = jobArray[job_index + SHIFT_IN_ARRAY];
                job_index++;
                done_jobs++;
	        } else {
	            pthread_mutex_unlock(mutex);
	            pthread_cond_signal(&check_rest);
	            //MPI_Send(&DONE, 1, MPI_INT, rank, 169, COMMW);
	            //something else?
	            break;
	        }
			if(REST_OF_JOBS ==(NUMBER_OF_JOBS - (job_index + 1))) {
			    	pthread_cond_signal(&check_rest);	
			    }
	    pthread_mutex_unlock(mutex);
	        
//        printf("[%d]: the {%d} is current job now. local_index == %d\n", rank, current_job, job_index);
        sleep(current_job);

    }
param = (void *) done_jobs; 
return param;
}

void * seacrchJob(void * param) {
	int no_jobs_counter = 0;
	int size = jobArray[1];
	int rank = jobArray[2];
	int response = 0;
	int rank_for_request = (rank + 1) % size;
	int rank_for_send = (rank + 1) % size;
	int * not_to_ask = (int *) malloc(sizeof(int) * size);
	memset(not_to_ask, 0, sizeof(int) * size);
	not_to_ask[rank] = 1;
	MPI_Status status;
	int iter_counter = 0;
	int rest = 0;
	int gotJob = 0;
	//int ACTIVE_R = jobArray[1] - 1; //size - 1;		

    while(no_jobs_counter < (jobArray[1] - 1)) {	//size - 1;
		pthread_mutex_lock(mutex);
			while(job_index < NUMBER_OF_JOBS - 1 - REST_OF_JOBS) {
			//while(REST_OF_JOBS != (NUMBER_OF_JOBS - (job_index + 1))) {
				pthread_cond_wait(&check_rest, mutex);
			}
		pthread_mutex_unlock(mutex);	    		

		while(!gotJob && ((size - 1) != no_jobs_counter)) {
    	//for(int j = 0; j < ACTIVE_R; j++) {
    		if(3 == rank) printf("[%d] before request extra job[%d]\n", rank, rank_for_request);
        	MPI_Send(&JOB_REQUEST, 1, MPI_INT, rank_for_request, 169, COMMW);
	    	MPI_Recv(&response, 1, MPI_INT, MPI_ANY_SOURCE, 170, COMMW, &status);
	    	if(3 == rank) printf("[%d] after request extra job[%d]\n", rank, rank_for_request);
	    	
	    	if(NO_JOBS == response) {
	    		no_jobs_counter++;
	    		not_to_ask[rank_for_request] = 1;
	    		if(3 == rank) printf("[%d] no_jobs_counter - %d (from %d)\n", rank, no_jobs_counter, rank_for_request);
	    		//other threads have no jobs -> current thread is done.
	    	}
	    	
	    	//get next rank for job request
			rank_for_request = (rank_for_request + 1) % size;
		    for(int i = 0; i < size; i++) {
			    if (not_to_ask[rank_for_request]) {
			    	rank_for_request = (rank_for_request + 1) % size;
			    }
	    	}

	    	if(NO_JOBS != response) {
	    		gotJob = 1;
	    		if(3 == rank) printf("[%d] has got extra job {%d}\n", rank, response);
	    		//add job to jobArray, perhaps expand it (and change NUMBER_OF_JOBS)
	    		pthread_mutex_lock(mutex);
	    			NUMBER_OF_JOBS++;
	    			jobArray[NUMBER_OF_JOBS + SHIFT_IN_ARRAY - 1] = response;
	    			//better use realloc
	    		pthread_mutex_unlock(mutex);	    
	    		break;		
	    	} 
	    		    
		}
		gotJob = 0;
	}
	printf("[%d] out of searching cycle\n", rank);
	free(not_to_ask);
}

void * responseForRequests(void * param) {
	int rank = jobArray[2];
	int current_job = 0;
	int requset = 0;
	int size = jobArray[1];
	MPI_Status status;
   while(TRUE) {
   		printf("in responseThread [%d]\n", rank);
    	MPI_Recv(&requset, 1, MPI_INT, MPI_ANY_SOURCE, 169, COMMW, &status);
    	printf("got requset: %d [%d]\n", requset, rank);
    	if(JOB_REQUEST == requset) {
	        pthread_mutex_lock(mutex);
	        	if(NUMBER_OF_JOBS == job_index) {
	        		MPI_Send(&NO_JOBS, 1, MPI_INT, status.MPI_SOURCE, 170, COMMW);
	        		pthread_mutex_unlock(mutex);
	        		break;
	        	}
	            if((NUMBER_OF_JOBS - REST_OF_JOBS) > job_index) {
	                current_job = jobArray[job_index + SHIFT_IN_ARRAY];
	                job_index++;
	            	MPI_Send(&current_job, 1, MPI_INT, status.MPI_SOURCE, 170, COMMW);
	                if((NUMBER_OF_JOBS - REST_OF_JOBS) == job_index) {
	                	pthread_cond_signal(&check_rest);
	                }
	            } else {
					MPI_Send(&NO_JOBS, 1, MPI_INT, status.MPI_SOURCE, 170, COMMW);
					pthread_mutex_unlock(mutex);
					break;
	            }
	        pthread_mutex_unlock(mutex);
		}
		if (DONE == requset) {
			break;
		}
	}
	printf("out of BIGResponse [%d]\n", rank);
	//for ervery process except self and one from main cycle
	for(int i = 2; i < size; i++) {	
		MPI_Recv(&requset, 1, MPI_INT, MPI_ANY_SOURCE, 169, COMMW, &status);
		MPI_Send(&NO_JOBS, 1, MPI_INT, status.MPI_SOURCE, 170, COMMW);
	}
	printf("out of Response [%d]\n", rank);
}

int main(int argc, char* argv[]){
	int rank = 0;
	int size = 0;
	int done_jobs = 0;
	pthread_t jobThread;
	pthread_t responseThread;
	pthread_t searchThread;
    int code = 0;
    int i = 0;        
    int provided = 0;

	MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    MPI_Comm_size(COMMW, &size);
    MPI_Comm_rank(COMMW, &rank);

	    mutex = (pthread_mutex_t *) malloc(sizeof(pthread_mutex_t));
		jobArray = (int *) malloc(sizeof(int) * NUMBER_OF_JOBS * size);
		memset(jobArray, 0, sizeof(int) * NUMBER_OF_JOBS * size);
	    job_index = 0;

	    jobArray[1] = size;
	    jobArray[2] = rank;
	    for(int i = 0; i < NUMBER_OF_JOBS; i++) {
	    	jobArray[i + SHIFT_IN_ARRAY] = rank + 1;
	    }

	    code = pthread_create(&jobThread, NULL, doJob, NULL);        
	    if (code != 0) {
	        printf("error occurs\n");
	        exit(1);
	    }

	    code = pthread_create(&responseThread, NULL, responseForRequests, NULL);        
	    if (code != 0) {
	        printf("error occurs\n");
	        exit(1);
	    }
	    code = pthread_create(&searchThread, NULL, seacrchJob, NULL);
	    if (code!=0) {
	        printf("error occurs\n");
	        exit(1);
	    }

	    pthread_join (jobThread, (void *)&done_jobs);
	    pthread_join (searchThread, NULL);
	    printf("[%d]waiting responseThread end\n", rank);
	    pthread_join (responseThread, NULL);
	    printf("%d has done %d jobs\n",rank, done_jobs);
	    MPI_Barrier(COMMW);
	    int ret_value = pthread_mutex_destroy(mutex);
	    if(ret_value) {
	        //my_perror(ret_value,argc,argv);
	        printf("error occurs\n");
	        exit(EXIT_FAILURE);
	    }
	     
	    ret_value = pthread_cond_destroy(&check_rest);
	    if(ret_value) {
	        //my_perror(ret_value,argc,argv);
	        printf("error occurs\n");
	        exit(EXIT_FAILURE);
	    }

	    free(jobArray);
	    free(mutex);
	    pthread_exit(EXIT_SUCCESS);
    MPI_Finalize();

return 0;
}