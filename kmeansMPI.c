#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#define DIM 3
int main(void) {
	
	MPI_Status status;
	int i, j, k, n, c;
	double dmin, dx;
	double *x, *mean, *sum;
	int *cluster, *count, color;
	int flips;
	
	int myrank, p;
	int tag = 0;
	MPI_Init (NULL, NULL);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &p);

	printf("inicio\n");
	if(myrank == 0){
		
		scanf("%d", &k);
		scanf("%d", &n);
		printf("inicio processamento myrank0 k=%d, n=%d \n", k, n);
		mean = (double *)malloc(sizeof(double)*DIM*k);
		sum= (double *)malloc(sizeof(double)*DIM*k);
		x = (double *)malloc(sizeof(double)*DIM*n);
		for (i = 0; i<k; i++)
			scanf("%lf %lf %lf", mean+i*DIM, mean+i*DIM+1, mean+i*DIM+2);
		for (i = 0; i<n; i++)
			scanf("%lf %lf %lf", x+i*DIM, x+i*DIM+1, x+i*DIM+2);
			
		

		printf("apos malloc myrank0\n");
		//trouxe o cluster p dentro pq ele precisa do n e o count precisa do k
		cluster = (int *)malloc(sizeof(int)*n);
		count = (int *)malloc(sizeof(int)*k);
		for (i = 0; i<n; i++) {
			cluster[i] = 0;
	}	 
		
		
		for( i=1; i<p; i++ ) {
			MPI_Send( &mean[i*k/p], k/p, MPI_DOUBLE,
					i,tag,MPI_COMM_WORLD );
			MPI_Send( &sum[i*k/p], k/p, MPI_DOUBLE,
					i, tag, MPI_COMM_WORLD );
			MPI_Send( &x[i*n/p], n/p, MPI_DOUBLE,
					i, tag, MPI_COMM_WORLD );
			MPI_Send( &k, k, MPI_INT,
					i, tag, MPI_COMM_WORLD );
			printf("enviou k=%d\n", k);
			MPI_Send( &n, n/p, MPI_INT,
					i, tag, MPI_COMM_WORLD );
		}
		    
		
		//processamento do myrank0
		flips = n;
		while (flips>0) {
			flips = 0;
			for (j = 0; j < k; j++) {
				count[j] = 0; 
				for (i = 0; i < DIM; i++) 
					sum[j*DIM+i] = 0.0;
			}
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
					flips++;
					cluster[i] = color;
		      	}
			}
	
		    for (i = 0; i < n; i++) {
				count[cluster[i]]++;
				for (j = 0; j < DIM; j++) 
					sum[cluster[i]*DIM+j] += x[i*DIM+j];
			}
			for (i = 0; i < (k/p); i++) {
				for (j = 0; j < DIM; j++) {
					mean[i*DIM+j] = sum[i*DIM+j]/count[i];
					printf("mean[i*DIM+j] %lf\n", mean[i*DIM+j]);
	  			}
			}
		}
		printf("terminou myrank0 k=%d\n", k);
		//processamento do myrank0
	} else {
		printf("inicio processamento myranki\n");
		MPI_Recv(&k, k, MPI_INT,
			0, tag, MPI_COMM_WORLD, &status );
		printf("recebeu k, k=%d\n", k);
		MPI_Recv(&n, n/p, MPI_INT,
			0, tag, MPI_COMM_WORLD, &status );
		printf("recebeu n, n=%d\n", n);
		MPI_Recv(mean, k/p, MPI_DOUBLE,
			0, tag, MPI_COMM_WORLD, &status );
		MPI_Recv(sum, k/p, MPI_DOUBLE,
			0, tag, MPI_COMM_WORLD, &status );
		MPI_Recv(x, n/p, MPI_DOUBLE,
			0, tag, MPI_COMM_WORLD, &status );
		
		printf("apos recv myranki k=%d, n=%d \n", k, n);
		
		for (i=0; i<n/p;i++){
			printf("cluster %d\n", cluster[i]);	
			}
		printf("final do for\n");
		//processamento do myranki
		flips = n;//?sera que aqui continua n?
		while (flips>0) {
			flips = 0;
			for (j = 0; j < k; j++) {
				count[j] = 0; 
				for (i = 0; i < DIM; i++) 
					sum[j*DIM+i] = 0.0;
			}
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
					flips++;
					cluster[i] = color;
		      	}
			}
	
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
		printf("terminou myranki\n");
		//processamento do myranki
			
	}
	
	/* Coleta do Resultado */
	if (myrank != 0) {
		printf("coleta myrank0\n");
		MPI_Send(mean, 1,MPI_DOUBLE,
				0,tag, MPI_COMM_WORLD);
	} else {
		for (i=1; i<p; i++){
			MPI_Recv(mean, 1, MPI_DOUBLE,
					i, tag, MPI_COMM_WORLD, &status);
		}
		
	}

    MPI_Finalize();

	//parte sequencial
	
	for (i = 0; i < k; i++) {
		dist = sqrt( pow (mean[i*DIM] - 0.0, 2) + pow (mean[i*DIM+1] - 0.0, 2) + pow (mean[i*DIM+2] - 0.0, 2) );
	}
	pontoMedio = [(x1+x2)/2, (y1+y2)/2, (z1+z2)/2]; 

	//resultado final fica no mean
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
	return(0);
}
