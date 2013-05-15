#include <mpi.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>

#define COMMW MPI_COMM_WORLD
#define BUFF_MAX_SIZE 256

int NUMBER_OF_JOBS = 17;
int * jobArray = NULL;
int job_index = 0;
pthread_mutex_t mutex;
int * done_jobs = NULL;

void my_perror(int code, int argc, char **argv)
{
    char buf[BUFF_MAX_SIZE];
    strerror_r(code, buf, sizeof buf);
    fprintf(stderr, "%s: pthread library error: %s\n", argv[0], buf);
}

void * thread_body(void * param) {
    int local_index = 0;
    int current_job = 0;
    do {
        pthread_mutex_lock(&mutex);
        local_index = job_index;
        if(local_index < NUMBER_OF_JOBS) {
            if(jobArray != NULL) {
                current_job = jobArray[local_index];
                done_jobs[0]++;
                job_index++;    
            }
        } else {
            pthread_mutex_unlock(&mutex);
            break;
        }
        pthread_mutex_unlock(&mutex);
        printf("Master works too {%d}\n", current_job);
        sleep(current_job);
    } while(job_index < NUMBER_OF_JOBS);
}

int main(int argc, char* argv[]){	
	int requset = 0;
	int rank = 0;
	int size = 0;
	int root = 0;
    int local_index = 0;
    int current_job = 0;
	int JOB_REQUEST = 123;
	int NO_JOBS = -1;
	MPI_Status status;
	pthread_t thread;
    int code = 0;
    int i = 0;        
    int provided = 0;
    int priv_done_jobs = 0;

	MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    MPI_Comm_size(COMMW, &size);
    MPI_Comm_rank(COMMW, &rank);

    if(root == rank) {
    	jobArray = (int *) malloc(sizeof(int) * NUMBER_OF_JOBS);
        done_jobs = (int *) malloc(sizeof(int) * size);
        //memset(done_jobs, 0, sizeof(int) * size);
        done_jobs[rank] = 0;

    	for(int i = 0; i < NUMBER_OF_JOBS; i++) {
    		jobArray[i] = (0 == i % 2) ? 1 + i * 2 : abs(30 - i * 3);
    	}
        //===== 
        code = pthread_create(&thread, NULL, thread_body, NULL);
        
        if (code!=0) {
            char buf[256];
            strerror_r(code, buf, sizeof buf);
            fprintf(stderr, "%s: creating thread: %s\n", argv[0], buf);
            exit(1);
        }
        //=====

    	
        do {
            pthread_mutex_lock(&mutex);
            local_index = job_index;
            if(local_index < NUMBER_OF_JOBS) {
                if(jobArray != NULL) {
                    current_job = jobArray[local_index];
                    job_index++;    
                }
            } else {
                pthread_mutex_unlock(&mutex);
                break;
            }
            pthread_mutex_unlock(&mutex);
    		MPI_Recv(&requset, 1, MPI_INT, MPI_ANY_SOURCE, 169, COMMW, &status);
    		if(JOB_REQUEST == requset) {
    			MPI_Send(&current_job, 1, MPI_INT, status.MPI_SOURCE, 169, COMMW);
    		}
    	} while (job_index < NUMBER_OF_JOBS);

    	for(int i = 1; i < size; i++) {
    		MPI_Recv(&requset, 1, MPI_INT, MPI_ANY_SOURCE, 169, COMMW, &status);
    		if(JOB_REQUEST == requset) {
    			MPI_Send(&NO_JOBS, 1, MPI_INT, status.MPI_SOURCE, 169, COMMW);
    		}	
    	}
        printf("Master has done %d\n", done_jobs[rank]);
    	
    }
    if(root != rank) {
    	while(1) {
	    	MPI_Send(&JOB_REQUEST, 1, MPI_INT, root, 169, COMMW);
	    	MPI_Recv(&requset, 1, MPI_INT, root, 169, COMMW, &status);
	    	if(NO_JOBS != requset) {
	    		printf("[%d] is doing {%d} job\n", rank, requset);
                //done_jobs[rank]++;
                priv_done_jobs++;
	    		sleep(requset);
	    	} else {
	    		printf("%d jobs are done.[%d]\n", priv_done_jobs, rank);
	    		break;
	    	}
	    }
    }

    if(root == rank) {
        pthread_join(thread,NULL);
        int ret_value = pthread_mutex_destroy(&mutex);
        if(ret_value)
        {
           my_perror(ret_value,argc,argv);
           printf("error occurs\n");
           exit(EXIT_FAILURE);
        }
        
        free(jobArray);
        pthread_exit(EXIT_SUCCESS);
    }

    MPI_Finalize();

return 0;
}