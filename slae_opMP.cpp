#include <omp.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <iostream>

#define WIDTH 100
#define HEIGHT 100
#define NUM_OF_THREADS 4
#define T 0.01
#define EPS 0.001
#define TYPE 1

int main(int argc, char* argv[])
{
	
	omp_set_dynamic(0);  
	omp_set_num_threads(NUM_OF_THREADS);

	double local_max = 0.0;
	double **matr  = (double **) malloc(sizeof(double *) * HEIGHT);
	double *vector = (double *) malloc(sizeof(double) * WIDTH);
    double *result = (double *) malloc(sizeof(double) * HEIGHT);
    double *b_vector = (double *) malloc(sizeof(double) * WIDTH);
	for(int i = 0; i < HEIGHT; i++)	{
		matr[i] = (double *) malloc(sizeof(double) * WIDTH);
	}

	int iter = 0;
	
	for(int i = 0; i < HEIGHT; i++){
		for(int j = 0; j < WIDTH; j++){
    			matr[i][j] = (double) (i == j) ? 2.0 : 1.0;
		}
    	b_vector[i] = (WIDTH + 1);
    	vector[i] = b_vector[i];
		result[i] = 0.0;
    }
	double norma_b = 0.0;
	#pragma omp parallel for
	for(int i = 0; i < HEIGHT; i++){
		norma_b += b_vector[i] * b_vector[i];
	}
	double prescision = norma_b * EPS * EPS;

	#if (1 == TYPE)
		while(true){

			#pragma omp parallel for
			 for(int i = 0; i < HEIGHT; i++) {
			 	//std::cout << omp_get_thread_num() << " ";//<< std::endl;
			 	result[i] = 0.0;
			 	for(int j = 0; j < WIDTH; j++) {
			 		result[i] += matr[i][j] * vector[j];
			 	}
			 	result[i] = T * (result[i] - b_vector[i]); // T * (Axk - b)
			 }
			 
			local_max = 0.0;
			#pragma omp parallel for reduction(+ : local_max)
			 for(int i = 0; i < HEIGHT; i++){
				 local_max += result[i] * result[i];
			}
			
			 if(sqrt(local_max) < prescision){
			 		break;
			 	} else {
			 		#pragma omp parallel for
			 		for(int i = 0; i < HEIGHT; i++){
			 			//std::cout << omp_get_thread_num() << std::endl;
			 			vector[i] = vector[i] - result[i];
			 		}
			 	}
		}
	#endif
	
	#if (2 == TYPE)
		#pragma omp parallel
		{
			while(true){
				#pragma omp for
				 for(int i = 0; i < HEIGHT; i++) {
			 		//std::cout << omp_get_thread_num() << " ";//<< std::endl;
			 		result[i] = 0.0;
			 		for(int j = 0; j < WIDTH; j++) {
			 			result[i] += matr[i][j] * vector[j];
			 		}
			 		result[i] = T * (result[i] - b_vector[i]); // T * (Axk - b)
				 }
			 
				local_max = 0.0;
				#pragma omp for reduction(+ : local_max)
				 for(int i = 0; i < HEIGHT; i++){
					 local_max += result[i] * result[i];
				}
			
				 if(sqrt(local_max) < prescision){
			 			break;
			 		} else {
			 			#pragma omp for
			 			for(int i = 0; i < HEIGHT; i++){
			 				//std::cout << omp_get_thread_num() << std::endl;
			 				vector[i] = vector[i] - result[i];
			 			}
			 		}

			}
	}
	#endif

	for(int i = 0; i < HEIGHT; i++){
		std::cout << vector[i] <<std::endl;
	}

	free(vector);
	free(b_vector);
	free(result);
	for(int i = 0; i < HEIGHT; i++)	{
		free(matr[i]);
	}
	free(matr);
	
	int i;	std::cin >> i;

return 0;
}
