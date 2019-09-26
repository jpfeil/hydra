#include<stdio.h>
#include<string.h>
#include "rand_index.h"

#define NELEMS(x) (sizeof(x) / sizeof((x)[0]))
#define MIN(a,b) (((a)<(b))?(a):(b))

long max(long *x, int size){
    //prlongf("Length: %ld\n", size);
	long _max = 0;
	for (int i=0; i < size; i++){
	    //printf("i: %d\n", i);
	    //printf("elem: %d\n", x[i]);
		if(*(x + i) > _max){
			_max = x[i];
		}
	}
	return _max;
}

void print_iarray(long *x, int size){
    for (int i=0; i < size; i++){
        printf("%ld\n", x[i]);
    }
}

double binomial(int n, int k){
    double buff[k + 1];
    memset(buff, 0.0, sizeof(double) * (k + 1));
    buff[0] = 1.0;
    for (int i=1; i < n + 1; i++){
        for (int j=MIN(i, k); j>0; j--){
            buff[j] = buff[j] + buff[j-1];
        }
    }
    return buff[k];
}


double crandi(long *x,
              int maxx,
              long *y,
              int maxy,
              int size){

	long matrix[maxx][maxy];
	double S1=0.0;
	double S2=0.0;
	double S3=0.0;

	// Set matrix values to zero
	memset(matrix, 0, sizeof(long) * maxx * maxy);

	// Fill matrix
	for (int i = 0; i < size; i++){
	    matrix[x[i]][y[i]] += 1;
	}

	// Rowsum
	long ai[maxx];
	memset(ai, 0, sizeof(long) * maxx);
	long bj[maxy];
	memset(bj, 0, sizeof(long) * maxy);

    // Sum Columns and Calculate S1
	for (int i=0; i < maxx; i++){
	    for (int j=0; j < maxy; j++){
	        ai[i] += matrix[i][j];
	        bj[j] += matrix[i][j];

	        S1 += binomial(matrix[i][j], 2);
	    }
	}

	for (int i=0; i < maxx; i++){
	    S2 += binomial(ai[i], 2);
	}

	for (int j=0; j < maxy; j++){
	    S3 += binomial(bj[j], 2);
	}

	double bit = (S2 * S3) / binomial(size, 2);
	return (S1 - bit) / (0.5 * (S2 + S3) - bit);
}


int main(){
    printf("Test Binomial:\n");
    printf("%ld\n", binomial(3, 2));
    printf("%ld\n", binomial(3, 3));
    printf("%ld\n", binomial(1000, 2));



	long X[] = {0, 0, 0, 1, 1, 1, 2, 2, 2};
	long Y[] = {0, 0, 1, 1, 1, 1, 2, 2, 3};
	int size = NELEMS(X);
	int maxx = max(X, size) + 1;
	int maxy = max(Y, size) + 1;
	double res = crandi(X, maxx, Y, maxy, size);
	printf("%f\n", res);
}
