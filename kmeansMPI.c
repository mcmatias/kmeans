#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#define DIM 3
int main(void) {
	
	MPI_Status status;
	int i, j, k, n, c;
	double dmin, dx;
	double *x, *mean, *sum, *meanTemp;
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
			MPI_Send( &k, 1, MPI_INT,
					i, tag, MPI_COMM_WORLD );
			printf("enviou k=%d\n", k);
			MPI_Send( &n, 1, MPI_INT,
					i, tag, MPI_COMM_WORLD );
			//MPI_Send( &mean, DIM*k, MPI_DOUBLE,
					//i,tag,MPI_COMM_WORLD );
			MPI_Send( &sum, DIM*k, MPI_DOUBLE,
					i, tag, MPI_COMM_WORLD );
			MPI_Send( &x[i*DIM*n/p], DIM*n/p, MPI_DOUBLE,
					i, tag, MPI_COMM_WORLD );
		}
		    
		
		//processamento do myrank0
		flips = n/p;
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
			printf("flips = %d\n", flips);
			for (i = 0; i < (k/p); i++) {
				for (j = 0; j < DIM; j++) {
					mean[i*DIM+j] = sum[i*DIM+j]/count[i];
					printf("mean[i*DIM+j] %lf, i= %d j = %d, myid=%d\n", mean[i*DIM+j], i, j, myrank);
	  			}
			}
		}
		printf("terminou myrank0 k=%d\n", k);
		//processamento do myrank0
	} else {

		printf("inicio processamento myranki\n");
		MPI_Recv(&k, 1, MPI_INT,
			0, tag, MPI_COMM_WORLD, &status );
		printf("recebeu k, k=%d\n", k);
		MPI_Recv(&n, 1, MPI_INT,
			0, tag, MPI_COMM_WORLD, &status );
		printf("recebeu n, n=%d\n", n);
		
		meanTemp = (double *)malloc(sizeof(double)*DIM*k);
		sum= (double *)malloc(sizeof(double)*DIM*k);
		x = (double *)malloc(sizeof(double)*DIM*n);
		
		//MPI_Recv(mean, DIM*k, MPI_DOUBLE,
			//0, tag, MPI_COMM_WORLD, &status );
		printf("recebeu mean\n");
		MPI_Recv(sum, DIM*k, MPI_DOUBLE,
			0, tag, MPI_COMM_WORLD, &status );
		printf("recebeu sum\n");
		MPI_Recv(x, DIM*n/p, MPI_DOUBLE,
			0, tag, MPI_COMM_WORLD, &status );
		
		printf("apos recv myranki k=%d, n=%d \n", k, n);
		
		cluster = (int *)malloc(sizeof(int)*n);
		count = (int *)malloc(sizeof(int)*k);
		for (i = 0; i<n; i++) {
			cluster[i] = 0;
		}	 
		
		for (i=0; i<n/p;i++){
			printf("cluster %d\n", cluster[i]);	
			}
		printf("final do for\n");
		//processamento do myranki
		flips = n/p;
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
						dx +=  (x[i*DIM+j] - meanTemp[c*DIM+j])*(x[i*DIM+j] - meanTemp[c*DIM+j]);
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
					meanTemp[i*DIM+j] = sum[i*DIM+j]/count[i];
	  			}
			}
		}
		printf("terminou myranki\n");
		//processamento do myranki
			
	}
	
	/* Coleta do Resultado */
	if (myrank != 0) {
		printf("coleta myrank i\n");
		MPI_Send(meanTemp, 1,MPI_DOUBLE,
				0,tag, MPI_COMM_WORLD);
		printf("enviou meanTemp\n");
	} else {
		for (i=1; i<p; i++){
			printf("entrou no for i=%d \n", i);
			MPI_Recv(meanTemp, 1, MPI_DOUBLE,
					i, tag, MPI_COMM_WORLD, &status);
			printf("recebeu meanTemp no for i=%d \n", i);
		}
		printf("saiu do for\n");
	}

    MPI_Finalize();
    
    printf("\n");
	printf("++++++++++++++++++++++++++++\n");
	printf("\n");
	printf("MPI FINALIZE\n");
	printf("\n");
	printf("++++++++++++++++++++++++++++\n");
	printf("\n");
	//parte sequencial
	
	int *flagMeanTemp;
	flagMeanTemp = (int *)malloc(sizeof(int)*k);
	printf("malloc flagMeanTemp\n");
	for (i = 0; i<k; i++) {
		flagMeanTemp[i] = 0;
		printf("flagMeanTemp[%d] = %d", i, flagMeanTemp[i]);
	}	
	
	for (int f = 0; f < k; f++) {//mean
		double menorDist = 1000.00;
		int indiceMeanTempMenorDist;
		double dist;
		for (int g = 0; g< k; g++) {//meanTemp
			if (flagMeanTemp[g] == 0) {//so calcula dist de meanTemp[g] se esse ponto nao tiver sido utilizado ainda
				dist = sqrt( pow (mean[f*DIM] - meanTemp[g*DIM], 2) + pow (mean[f*DIM+1] - meanTemp[g*DIM+1], 2) + pow (mean[f*DIM+2] - meanTemp[g*DIM+2], 2) );
				if (dist<menorDist) {
					menorDist = dist;
					indiceMeanTempMenorDist = g;
				}
			}
		}
		//ao encontrar os dois pontos de menor distancia
		mean[f*DIM] = (mean[f*DIM] + mean[indiceMeanTempMenorDist*DIM])/2;
		mean[f*DIM+1] = (mean[f*DIM+1] + mean[indiceMeanTempMenorDist*DIM+1])/2;
		mean[f*DIM+2] = (mean[f*DIM+2] + mean[indiceMeanTempMenorDist*DIM+2])/2;
		flagMeanTemp[indiceMeanTempMenorDist] = 1;
	}

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
