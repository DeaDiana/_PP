// grid_slae.cpp : Defines the entry point for the console application.
//

#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>

#define HEIGHT 9
#define WIDTH 21
#define EPS 0.0001

using namespace std;

void calc_intervals()
{
}

void function_value_F()
{
}

int main()
{
	double **current_matrix = NULL;
	double **last_matrix = NULL;
	double hx_step = 1.0;
	double hy_step = 1.0;
	double differ = 0.0;
	double max = 0.0;
	double **swap = NULL;
	ifstream input_for_task;

	current_matrix = (double **) malloc(sizeof(double *) * HEIGHT);
	last_matrix = (double **) malloc(sizeof(double *) * HEIGHT);

	if((NULL == current_matrix) || (NULL == last_matrix))
	{
		printf("memory was not allocated\n");
		exit(1);
	}

	for(int i = 0; i < HEIGHT; i++)
	{
		current_matrix[i] = (double *) malloc(sizeof(double) * WIDTH);
		last_matrix[i] = (double *) malloc(sizeof(double) * WIDTH);
		if((NULL == current_matrix[i]) || (NULL == last_matrix[i]))
		{
			printf("memory was not allocated\n");
			exit(1);
		}
	}

	input_for_task.open("matrix.txt");

	for(int i = 0; i < HEIGHT; i++) 
	{
		for(int j = 0; j < WIDTH; j++)
		{
			input_for_task >> last_matrix[i][j];
			std::cout << last_matrix[i][j] << " ";
		}
		std::cout << std::endl;
	}

	for(int i = 0; i < HEIGHT; i++) 
	{
		for(int j = 0; j < WIDTH; j++)
		{
			current_matrix[i][j] = last_matrix[i][j];
		}
	}

	do
	{
		max = 0.0;
		for(int i = 1; i < (HEIGHT - 1); i++)
		{
			for(int j = 1; j < (WIDTH - 1); j++)
			{
				current_matrix[i][j] = (1 / (2 / (hx_step * hx_step) + 2 / (hy_step * hy_step))) * ((current_matrix[i][j - 1] + last_matrix[i][j + 1]) / (hx_step * hx_step) + (current_matrix[i - 1][j] + last_matrix[i + 1][j]) / (hy_step * hy_step));
				//current_matrix[i][j] = 0.25 * ((current_matrix[i][j - 1] + last_matrix[i][j + 1]) / (hx_step * hx_step) + (current_matrix[i - 1][j] + last_matrix[i + 1][j]) / (hy_step * hy_step));
				//current_matrix[i][j] = 0.25 * ((last_matrix[i][j - 1] + last_matrix[i][j + 1]) / (hx_step * hx_step) + (last_matrix[i - 1][j] + last_matrix[i + 1][j]) / (hy_step * hy_step));
				differ = abs(current_matrix[i][j] - last_matrix[i][j]);
				if(differ > max)
				{
					max = differ;
				}
			}
		}
		std::cout << max << std::endl;
		swap = current_matrix;
		current_matrix = last_matrix;
		last_matrix = swap;
	} while(max > EPS);


	for(int i = 0; i < HEIGHT; i++) 
	{
		for(int j = 0; j < WIDTH; j++)
		{
			//cout << current_matrix[i][j] << " ";
			printf("%.0lf ", current_matrix[i][j]);
		}
		cout << endl;
	}

	for(int i = 0; i < HEIGHT; i++)
	{
		free(current_matrix[i]);
		free(last_matrix[i]);
	}

	free(current_matrix);
	free(last_matrix);

return 0;
}