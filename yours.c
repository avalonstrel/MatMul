#include <stdlib.h>
#include <sys/time.h>
#include <printf.h>
#include "handout.h"
#include "time.h"

//these are the implemented methods in 'handout.o' :

//double *generate(int n)
//bool check(double *Yours, double *A, double *B, int n)
//void printMatrix(double *M, int n)
//double *Mult(double *A, double *B, int n)

//the timer should be implemented in the YoursBlocked & YoursRecursive function and printed out in a format like "TIME: 0.000000 seconds"


double *YoursBlocked(int n, double *A, double *B) {
    double *a;
    a = (double *) malloc(n * n * sizeof(double));
    time_t start = clock();
    unsigned int block_size = 32;
    double sum = 0.0;

    for(int j=0; j < n/block_size; j++){
        int j_all = j*block_size;
        for(int k=0; k < n/block_size; k++){
            int k_all = k*block_size;
            for(int i=0; i < n/block_size; i++){
                int i_all = i*block_size;
                for(int l=0; l < block_size; l++){
                    int B_all_ind_ = j_all+l + k_all*n;
                    for(int t=0; t < block_size; t++){
                        sum = 0.0;
                        int A_all_ind = (i_all+t)*n + k_all;
                        int m = 0;
                        int B_all_ind =  B_all_ind_;
                        while(m < block_size){
                            sum += A[A_all_ind] * B[B_all_ind];
                            B_all_ind += n;
                            A_all_ind += 1;
                            m++;
                        }
/*                        for(int m=0; m < block_size; m++){

                            sum += A[A_all_ind+m] * B[(k_all+m)*n+j_all+l];
                        }*/
                        a[(i_all+t)*n+j_all+l] += sum;
                    }
                }
            }
        }
    }
    
    time_t end = clock();
    printf("Time %f\n", (double)(end - start)/CLOCKS_PER_SEC);
// fill your code here, a is your output matrix
    return a;
}
void YoursRecursive(double* a, double *A, double *B, int n, int i, int j, int k, int stride) {
    // fill your code here, a is your output matrix
    if( n > 16){
        int new_n = n/2;
        YoursRecursive(a, A, B, new_n, i, j, k, stride);
        YoursRecursive(a, A, B, new_n, i, j, k + new_n, stride);
        YoursRecursive(a, A, B, new_n, i + new_n, j, k, stride);
        YoursRecursive(a, A, B, new_n, i + new_n, j, k + new_n, stride);
        YoursRecursive(a, A, B, new_n, i, j + new_n, k, stride);
        YoursRecursive(a, A, B, new_n, i, j + new_n, k + new_n, stride);
        YoursRecursive(a, A, B, new_n, i + new_n, j + new_n, k, stride);
        YoursRecursive(a, A, B, new_n, i + new_n, j + new_n, k + new_n, stride);
    }else{
        double sum = 0.0;
        for(int i_ = i; i_ < i + n ; i_++){
            for(int j_ = j; j_ < j + n; j_++){
                sum = 0.0;
                for(int k_ = k; k_ < k + n; k_++){
                    sum += A[i_*stride + k_] * B[k_*stride + j_];
                }
                a[i_*stride + j_] += sum;
            }
        }
    }

}

int main(int argc, char *argv[]) {
    srand((unsigned int) time(NULL));
    int n = atoi(argv[1]);
    double *A, *B;
    A = generate(n);
    B = generate(n);
//    printf("A\n");
//    printMatrix(A, n);
//    printf("B\n");
//    printMatrix(B, n);
    double *Y;
    Y = (double *) malloc(n * n * sizeof(double));
    Y = generate(n);
    Y = YoursBlocked(n,A,B);
//    Y = Mult(A, B, n);
//    printf("Y\n");
//    printMatrix(Y, n);
    if (check(Y, A, B, n))
        printf("B TRUE%d\n", 1);
    else
        printf("B FALSE%d\n", 0);
    double *a;
    a = (double *) malloc(n * n * sizeof(double));
    int block_inds[] = {0, n, 0, n}; 
    time_t start = clock();   
    YoursRecursive(a,A,B, n, 0, 0, 0, n);
    time_t end = clock();
    printf("Time %f\n", (double)(end - start)/CLOCKS_PER_SEC);
    if (check(Y, A, B, n))
        printf("R TRUE%d\n", 1);
    else
        printf("R FALSE%d\n", 0);
    free(A);
    free(B);
    free(Y);
}