#include <omp.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <iostream>

#define WIDTH 100
#define HEIGHT 100
#define NUM_OF_THREADS 4
#define TYPE 1

int main()
{
	omp_set_dynamic(0);  
	omp_set_num_threads(NUM_OF_THREADS);

	double *vector = (double *) malloc(sizeof(double) * WIDTH);
	double **matr  = (double **) malloc(sizeof(double *) * HEIGHT);
	double *result = (double *) malloc(sizeof(double) * HEIGHT);
    int *starts_v = (int *) malloc(sizeof(int) * (NUM_OF_THREADS + 1)); 
	int div_h = HEIGHT / NUM_OF_THREADS;
    int mod_h = HEIGHT % NUM_OF_THREADS;

	for(int i = 0; i < HEIGHT; i++)	{
		matr[i] = (double *) malloc(sizeof(double) * WIDTH);
	}

	for(int i = 0; i < HEIGHT; i++)	{
		for(int j = 0; j < WIDTH; j++) {
			matr[i][j] = (i == j) ? 2.0 : 1.0;
		}
	}

	for(int i = 0; i < WIDTH; i++)	{
		vector[i] = 1.0;
	}

	memset(result, 0, sizeof(double) * HEIGHT);

	#if 1 == TYPE

	#pragma omp parallel
	{
		#pragma omp for
			for(int i = 0; i < HEIGHT; i++) {
				//std::cout << omp_get_thread_num() << std::endl;
				for(int j = 0; j < WIDTH; j++) {
					result[i] += matr[i][j] * vector[j];
				}
			}
	}

	#endif

	#if 2 == TYPE

		#pragma omp parallel for
			for(int i = 0; i < HEIGHT; i++) {
				//std::cout << omp_get_thread_num() << std::endl;
				for(int j = 0; j < WIDTH; j++) {
					result[i] += matr[i][j] * vector[j];
				}
			}

	#endif

	#if 3 == TYPE

		#pragma omp parallel
		{
			int rank = omp_get_thread_num();
			int size = omp_get_num_threads();

			for(int i = 0; i < size; i++){
				starts_v[i] = (i <= mod_h) ? (div_h + 1) * i : (div_h + 1) * mod_h + div_h * (i - mod_h);
			}
			starts_v[size] = HEIGHT;

			for(int i = starts_v[rank]; i < starts_v[rank + 1]; i++){
				//std::cout << rank << std::endl;
				for(int j = 0; j < WIDTH; j++) {
						result[i] += matr[i][j] * vector[j];
					}
			}
		}

	#endif

	#if 4 == TYPE
		#pragma omp parallel
		{
			#pragma  omp sections
			{
				#pragma omp section
				{
					//std::cout << omp_get_thread_num() << std::endl;
					for(int i = 0; i < HEIGHT / 2; i++){
						for(int j = 0; j < WIDTH; j++) {
							result[i] += matr[i][j] * vector[j];
						}
					}
				}
				#pragma omp section
				{
					//std::cout << omp_get_thread_num() << std::endl;
					for(int i = HEIGHT / 2; i < HEIGHT; i++){
						for(int j = 0; j < WIDTH; j++) {
							result[i] += matr[i][j] * vector[j];
						}
					}
				}
			}
		}
	#endif

	for(int i = 0; i < HEIGHT; i++)	{
		std::cout << result[i] << " ";
	}

	free(vector);
	free(result);
	free(starts_v);
	for(int i = 0; i < HEIGHT; i++)	{
		free(matr[i]);
	}
	free(matr);

	int i;	std::cin >> i;

return 0;
}