#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#define DIM 3
int main(void) {
	
	clock_t t1Inicio, t1Fim, t1Decorrido;
	clock_t t2Inicio, t2Fim, t2Decorrido;
	clock_t t3Inicio, t3Fim, t3Decorrido;
	clock_t t4Inicio, t4Fim, t4Decorrido;
	clock_t t5Inicio, t5Fim, t5Decorrido;
	clock_t t6Inicio, t6Fim, t6Decorrido;
	int i, j, k, n, c;
	double dmin, dx;
	double *x, *mean, *sum;
	int *cluster, *count, color;
	int flips;
	scanf("%d", &k);
	scanf("%d", &n);
	x = (double *)malloc(sizeof(double)*DIM*n);
	mean = (double *)malloc(sizeof(double)*DIM*k);
	sum= (double *)malloc(sizeof(double)*DIM*k);
	cluster = (int *)malloc(sizeof(int)*n);
	count = (int *)malloc(sizeof(int)*k);
	for (i = 0; i<n; i++) 
		cluster[i] = 0;
	for (i = 0; i<k; i++)
		scanf("%lf %lf %lf", mean+i*DIM, mean+i*DIM+1, mean+i*DIM+2);
	for (i = 0; i<n; i++)
		scanf("%lf %lf %lf", x+i*DIM, x+i*DIM+1, x+i*DIM+2);
	flips = n;
	while (flips>0) {
		flips = 0;
		t1Inicio = clock();
		for (j = 0; j < k; j++) {
			count[j] = 0; 
			for (i = 0; i < DIM; i++) 
				sum[j*DIM+i] = 0.0;
		}
		t1Fim = clock();
		
		t2Inicio = clock();
		#pragma omp parallel private(color, dmin) shared(flips, dx) reduction(+:dx)
		{
			#pragma omp for 
			//shared flips, dx, 
			//private color, dmin
			for (i = 0; i < n; i++) {
				dmin = -1; color = cluster[i];
				t3Inicio = clock();
				for (c = 0; c < k; c++) {
					dx = 0.0;
					t4Inicio = clock();
					for (j = 0; j < DIM; j++)
						//reduction?
						dx +=  (x[i*DIM+j] - mean[c*DIM+j])*(x[i*DIM+j] - mean[c*DIM+j]);
					if (dx < dmin || dmin == -1) {
						color = c;
						dmin = dx;
					}
					t4Fim = clock();
				}
				t3Fim = clock();
				if (cluster[i] != color) {
					flips++;
					cluster[i] = color;
		      	}
			}
		}
		t2Fim = clock();

		t5Inicio = clock();
	    for (i = 0; i < n; i++) {
			count[cluster[i]]++;
			for (j = 0; j < DIM; j++) 
				sum[cluster[i]*DIM+j] += x[i*DIM+j];
		}
		t5Fim = clock();
		t6Inicio = clock();
		for (i = 0; i < k; i++) {
			for (j = 0; j < DIM; j++) {
				mean[i*DIM+j] = sum[i*DIM+j]/count[i];
  			}
		}
		t6Fim = clock();
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
	t1Decorrido = ((t1Inicio - t1Fim)/ (CLOCKS_PER_SEC / 1000));
	printf("tempo1 %Lf\n", (long double) t1Decorrido);
	t2Decorrido = ((t2Inicio - t2Fim)/ (CLOCKS_PER_SEC / 1000));
	printf("tempo2 %Lf\n", (long double) t2Decorrido);
	t3Decorrido = ((t3Inicio - t3Fim)/ (CLOCKS_PER_SEC / 1000));
	printf("tempo3 %Lf\n", (long double) t3Decorrido);
	t4Decorrido = ((t4Inicio - t4Fim)/ (CLOCKS_PER_SEC / 1000));
	printf("tempo4 %Lf\n", (long double) t4Decorrido);
	t5Decorrido = ((t5Inicio - t5Fim)/ (CLOCKS_PER_SEC / 1000));
	printf("tempo5 %Lf\n", (long double) t5Decorrido);
	t6Decorrido = ((t6Inicio - t6Fim)/ (CLOCKS_PER_SEC / 1000));
	printf("tempo6 %Lf\n", (long double) t6Decorrido);
	return(0);
}
