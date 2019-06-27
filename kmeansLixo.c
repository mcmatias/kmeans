#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <omp.h>
#define DIM 3
int main(void) {
	
	clock_t t2Inicio, t2Fim, t2Decorrido;
	
	int i, j, k, n, c;
	
	double *x, *mean, *sum;
	int *cluster, *count;
	int flips;
	scanf("%d", &k);
	scanf("%d", &n);
	x = (double *)malloc(sizeof(double)*DIM*n);
	mean = (double *)malloc(sizeof(double)*DIM*k);
	sum= (double *)malloc(sizeof(double)*DIM*k);
	cluster = (int *)malloc(sizeof(int)*n);
	count = (int *)malloc(sizeof(int)*k);
	
	t2Inicio = clock();
	for (i = 0; i<n; i++) 
		cluster[i] = 0;
	for (i = 0; i<k; i++)
		scanf("%lf %lf %lf", mean+i*DIM, mean+i*DIM+1, mean+i*DIM+2);
	for (i = 0; i<n; i++)
		scanf("%lf %lf %lf", x+i*DIM, x+i*DIM+1, x+i*DIM+2);
	flips = n;
	while (flips>0) {
		flips = 0;
		for (j = 0; j < k; j++) {
			count[j] = 0; 
			for (i = 0; i < DIM; i++) 
				sum[j*DIM+i] = 0.0;
		}
		

			//aqui eh a duvida
			//se um processo pode interferir no outro
			//como devo fazer?
			#pragma omp parallel private(color, dmin, dx, flips1)
			{	
				double dmin, dx;
				int color;
				//int color1 = color;
				//double dmin1 = dmin;
				//int clusteri = cluster[i];
				int flips1 = flips;
				#pragma omp for
				{
					for (i = 0; i < n; i++) {
					dmin = -1; color = cluster[i];
							for (c = 0; c < k; c++) {
								dx = 0.0;
								for (j = 0; j < DIM; j++) 
									dx +=  (x[i*DIM+j] - mean[c*DIM+j])*(x[i*DIM+j] - mean[c*DIM+j]);
								if (dx < dmin || dmin == -1) {
									color = c;
									dmin = dx;
								}
							}
					
					if (cluster[i] != color) {
						flips1++;
						cluster[i] = color;
			      	}
				}//loop mais externo
			}//omp for
			#pragma omp critical
			{
				flips = flips1;
			}
		}//fim da regiao paralela

	    for (i = 0; i < n; i++) {
			count[cluster[i]]++;
			for (j = 0; j < DIM; j++) 
				sum[cluster[i]*DIM+j] += x[i*DIM+j];
		}
		for (i = 0; i < k; i++) {
			for (j = 0; j < DIM; j++) {
				mean[i*DIM+j] = sum[i*DIM+j]/count[i];
  			}
		}
	}
	for (i = 0; i < k; i++) {
		for (j = 0; j < DIM; j++)
			printf("%5.2f ", mean[i*DIM+j]);
		printf("\n");
	}
	#ifdef DEBUG
	for (i = 0; i < n; i++) {
		for (j = 0; j < DIM; j++)
			printf("%5.2f ", x[i*DIM+j]);
		printf("%d\n", cluster[i]);
	}
	#endif
	t2Fim = clock();
	t2Decorrido = ((t2Fim - t2Inicio)/ (CLOCKS_PER_SEC / 1000));
	printf("tempo2 %Lf\n", (long double) t2Decorrido);
	return(0);
}
