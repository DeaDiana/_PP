#include <mpi.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define COMMW MPI_COMM_WORLD
#define NUM_OF_ARGS 3

int main(int argc, char* argv[]) {
	int rank = 0, size = 0;
	int left_range = -1, right_range = 1, interval = right_range - left_range + 1;
	int swap = 0, i = 0, iterations = 0;
	const int ROOT = 0;
	const int FREE_COEF = 5;
	const int STEP = 1;
	const int NUM_OF_THREADS = 4;
    int div_intrv = 0, mod_intrv = 0;
    double loc_result = 0.0, result = 0.0; int x = 0;
    int *num_points = NULL, *starts_intrv = NULL, *ends_intrv = NULL;
    double *func_value = NULL;

	if(NUM_OF_ARGS == argc) {
		left_range = atoi(argv[1]);
		right_range = atoi(argv[2]);
		if(left_range > right_range) {
			swap = left_range;
			left_range = right_range;
			right_range = swap;
		}
		interval = right_range - left_range + 1;
	}

    omp_set_dynamic(0);  
    omp_set_num_threads(NUM_OF_THREADS);

    MPI_Init(&argc, &argv);
    MPI_Comm_size(COMMW, &size);
    MPI_Comm_rank(COMMW, &rank);

	div_intrv = interval / size;
	mod_intrv = interval % size;

	num_points = (int *) malloc(sizeof(int) * size);
	starts_intrv = (int *) malloc(sizeof(int) * size);
	ends_intrv = (int *) malloc(sizeof(int) * size);

	if((NULL == num_points) || (NULL == starts_intrv) || (NULL == ends_intrv)) {
        fprintf(stderr, "memory  wasn't allocated\n");
        exit(1);
    }

    #pragma omp parallel for private(i)
	for(i = 0; i < size; i++) {
		num_points[i] = (i < mod_intrv) ? (div_intrv + 1) : div_intrv;
	}
	starts_intrv[0] = left_range;
	ends_intrv[0] = starts_intrv[0] + num_points[0] - 1;
	for (i = 1; i < size; i++) {
		starts_intrv[i] = ends_intrv[i - 1];
		ends_intrv[i] = starts_intrv[i] + num_points[i];
	}
	#pragma omp parallel for private(i)
	for(i = 1; i < size; i++) {
		num_points[i]++;
	}
	func_value = (double *) malloc(sizeof(double) * num_points[rank]);
	
	if(NULL == func_value) {
        fprintf(stderr, "memory  wasn't allocated\n");
        exit(1);
    }

    #pragma omp parallel for private(x)
		for(x = starts_intrv[rank]; x <= ends_intrv[rank]; x++) {
			func_value[x - starts_intrv[rank]] = x * sin(x) + FREE_COEF;
	}

	#pragma omp parallel for private(i) reduction(+ : loc_result)
	for(i = 0; i < (num_points[rank] - 1); i++) {
		loc_result += ((func_value[i] + func_value[i + 1]) / 2) * STEP;
	}

	MPI_Reduce(&loc_result, &result, 1, MPI_DOUBLE, MPI_SUM, ROOT, COMMW);

	if(ROOT == rank) {
		printf("%f\n", result);
	}

    MPI_Finalize();

	free(num_points);
	free(starts_intrv);
	free(func_value);

	return 0;
}