#include <mpi.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>

#define COMMW MPI_COMM_WORLD
#define BUFF_SIZE 1024
#define REST_OF_JOBS 2
#define SHIFT_IN_ARRAY 3

int job_index = 0;
int NUMBER_OF_JOBS = 6;
int * jobArray = NULL;
int JOB_REQUEST = 123;
int NO_JOBS = -1;
pthread_mutex_t * mutex = NULL;
pthread_cond_t check_rest = PTHREAD_COND_INITIALIZER;

void * doJob(void * param) {
    int current_job = 0;
    int done_jobs = 0;
    int rank = jobArray[2];
    while(1) {
        pthread_mutex_lock(mutex);
	        if(job_index < NUMBER_OF_JOBS) {
                current_job = jobArray[job_index + SHIFT_IN_ARRAY];
                job_index++;
                done_jobs++;
	        } else {
	            pthread_cond_signal(&check_rest);
	            pthread_mutex_unlock(mutex);
	            break;
	        }
			if(REST_OF_JOBS == (NUMBER_OF_JOBS - (job_index + 1))) {
			    	pthread_cond_signal(&check_rest);	
			    }
	    pthread_mutex_unlock(mutex);
	        
        printf("[%d]: the {%d} is current job now. job_index == %d\n", rank, current_job, job_index);
        sleep(current_job);
    }
	pthread_mutex_lock(mutex);
		jobArray[0] = done_jobs;
	pthread_mutex_unlock(mutex);
return param;
}

void * seacrchJob(void * param) {
	int no_jobs_counter = 0;
	pthread_mutex_lock(mutex);
		int size = jobArray[1];
		int rank = jobArray[2];
	pthread_mutex_unlock(mutex);
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

    while(no_jobs_counter < (jobArray[1] - 1)) {	//size - 1;
		pthread_mutex_lock(mutex);
			//while(job_index < NUMBER_OF_JOBS - 1 - REST_OF_JOBS) {
			while(REST_OF_JOBS < (NUMBER_OF_JOBS - (job_index + 1))) {
				pthread_cond_wait(&check_rest, mutex);
			}
		pthread_mutex_unlock(mutex);	    		

		while(!gotJob && ((size - 1) != no_jobs_counter)) {
        	MPI_Send(&JOB_REQUEST, 1, MPI_INT, rank_for_request, 169, COMMW);
	    	MPI_Recv(&response, 1, MPI_INT, MPI_ANY_SOURCE, 170, COMMW, &status);
	    	
	    	if(NO_JOBS == response) {
	    		no_jobs_counter++;
	    		not_to_ask[rank_for_request] = 1;
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
	    		//printf("[%d] has got extra job {%d}\n", rank, response);
	    		//add job to jobArray, perhaps expand it (and change NUMBER_OF_JOBS)
	    		pthread_mutex_lock(mutex);
	    			NUMBER_OF_JOBS++;
	    			jobArray[NUMBER_OF_JOBS + SHIFT_IN_ARRAY - 1] = response;
	    		pthread_mutex_unlock(mutex);	    
	    		break;		
	    	} 
		}
		gotJob = 0;
	}
	//printf("[%d] out of searching cycle\n", rank);
	free(not_to_ask);
}

void * responseForRequests(void * param) {
	pthread_mutex_lock(mutex);
		int rank = jobArray[2];
	pthread_mutex_unlock(mutex);
	int current_job = 0;
	int requset = 0;
	int size = jobArray[1];
	MPI_Status status;
   while(1) {
    	MPI_Recv(&requset, 1, MPI_INT, MPI_ANY_SOURCE, 169, COMMW, &status);
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
	}
	//for ervery process except self and one from main cycle
	for(int i = 2; i < size; i++) {	
		MPI_Recv(&requset, 1, MPI_INT, MPI_ANY_SOURCE, 169, COMMW, &status);
		MPI_Send(&NO_JOBS, 1, MPI_INT, status.MPI_SOURCE, 170, COMMW);
	}
	//printf("[%d] out of Response\n", rank);
}

int main(int argc, char* argv[]){
	int rank = 0;
	int size = 0;
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
	    code = pthread_mutex_init(mutex, NULL);
		jobArray = (int *) malloc(sizeof(int) * BUFF_SIZE);
		memset(jobArray, 0, sizeof(int) * BUFF_SIZE);
	    job_index = 0;
	    jobArray[1] = size;
	    jobArray[2] = rank;
	    for(int i = 0; i < NUMBER_OF_JOBS; i++) {
	    	jobArray[i + SHIFT_IN_ARRAY] = rank + 1;
	    }

	    code = pthread_create(&jobThread, NULL, doJob, NULL);        
	    if (code != 0) {
	        printf("[%d] error occurs on creating {jobThread} ERROR: %d\n", rank, code);
	        exit(1);
	    }
	    code = pthread_create(&responseThread, NULL, responseForRequests, NULL);        
	    if (code != 0) {
	        printf("[%d] error occurs on creating {responseThread} ERROR: %d\n", rank, code);
	        exit(1);
	    }
	    code = pthread_create(&searchThread, NULL, seacrchJob, NULL);
	    if (code!=0) {
	        printf("[%d] error occurs on creating {searchThread} ERROR: %d\n", rank, code);
	        exit(1);
	    }

	    pthread_join (searchThread, NULL);
	    pthread_join (jobThread, NULL);
	    pthread_join (responseThread, NULL);
	    printf("%d has done %d jobs\n", rank, jobArray[0]);
	    MPI_Barrier(COMMW);
	    int ret_value = pthread_mutex_destroy(mutex);
	    if(ret_value) {
	        printf("[%d] error occurs on destroying {mutex} ERROR: %d\n", rank, ret_value);
	        exit(EXIT_FAILURE);
	    }     
	    ret_value = pthread_cond_destroy(&check_rest);
	    if(ret_value) {
	        printf("[%d] error occurs on destroying {condition variable} ERROR: %d\n", rank,  ret_value);
	        exit(EXIT_FAILURE);
	    }
	    free(jobArray);
	    free(mutex);
	    pthread_exit(EXIT_SUCCESS);
    MPI_Finalize();

return 0;
}