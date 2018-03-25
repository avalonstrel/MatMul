#include <stdlib.h>
#include <sys/time.h>
#include <printf.h>
#include "handout.h"
#include "time.h"


#define BLOCK_SIZE 256

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
    unsigned int block_size = BLOCK_SIZE;
    
    if(n < block_size){
        block_size = n;
    }
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
void StrassenRecursiveImpl(double *a, double* A, double*B, int n, int i_A, int j_A, int i_B, int j_B, int A_stride, int B_stride, int O_stride);
//n n/2 
void p1(double *p1, double *A, double *B, int n, int i_A, int j_A, int i_B, int j_B,  int A_stride, int B_stride){
    double * tmp = (double *) malloc(n * n * sizeof(double));
    for(int i=0; i < n; i++){
        for(int j=0; j < n; j++){
            tmp[i*n+j] = B[(i_B+i)*B_stride+j_B+j+n] - B[(i_B+n+i)*B_stride+j_B+j+n];
        }
    }

    StrassenRecursiveImpl(p1, A, tmp, n, i_A, j_A, 0, 0, A_stride, n, n);
}

void p2(double *p2, double *A, double *B, int n, int i_A, int j_A, int i_B, int j_B,  int A_stride, int B_stride){
    double * tmp = (double *) malloc(n * n * sizeof(double));
    for(int i=0; i < n; i++){
        for(int j=0; j < n; j++){
            tmp[i*n+j] = A[(i_A+i)*A_stride+j_A+j] + A[(i_A+i)*A_stride+j_A+j+n];
        }
    }
    StrassenRecursiveImpl(p2, tmp, B, n, 0, 0, i_B + n, j_B + n, n, B_stride, n);
}
void p3(double *p3, double *A, double *B, int n, int i_A, int j_A, int i_B, int j_B,  int A_stride, int B_stride){
    double * tmp = (double *) malloc(n * n * sizeof(double));
    for(int i=0; i < n; i++){
        for(int j=0; j < n; j++){
            tmp[i*n+j] = A[(i_A+n+i)*A_stride+j_A+j] + A[(i_A+n+i)*A_stride+j_A+j+n];
        }
    }
    StrassenRecursiveImpl(p3, tmp, B, n, 0, 0, i_B, j_B, n, B_stride, n);    
}
void p4(double *p4, double *A, double *B, int n, int i_A, int j_A, int i_B, int j_B,  int A_stride, int B_stride){
    double * tmp = (double *) malloc(n * n * sizeof(double));
    for(int i=0; i < n; i++){
        for(int j=0; j < n; j++){
            tmp[i*n+j] = B[(i_B+n+i)*B_stride+j_B+j] - B[(i_B+i)*B_stride+j_B+j];
        }
    }

    StrassenRecursiveImpl(p4, A, tmp, n, i_A + n, j_A + n, 0, 0, A_stride, n, n);
}
void p5(double *p5, double *A, double *B, int n, int i_A, int j_A, int i_B, int j_B,  int A_stride, int B_stride){
    double * tmp = (double *) malloc(n * n * sizeof(double));
    for(int i=0; i < n; i++){
        for(int j=0; j < n; j++){
            tmp[i*n+j] = A[(i_A+i)*A_stride+j_A+j] + A[(i_A+n+i)*A_stride+j_A+n+j];
        }
    }
    double * tmp2 = (double *) malloc(n * n * sizeof(double));
    for(int i=0; i < n; i++){
        for(int j=0; j < n; j++){
            tmp2[i*n+j] = B[(i_B+i)*B_stride+j_B+j] + B[(i_B+n+i)*B_stride+j_B+n+j];
        }
    }    
    StrassenRecursiveImpl(p5, tmp, tmp2, n, 0, 0, 0, 0, n, n, n);
}
void p6(double *p6, double *A, double *B, int n, int i_A, int j_A, int i_B, int j_B,  int A_stride, int B_stride){
    double * tmp = (double *) malloc(n * n * sizeof(double));
    for(int i=0; i < n; i++){
        for(int j=0; j < n; j++){
            tmp[i*n+j] = A[(i_A+i)*A_stride+j_A+n+j] - A[(i_A+n+i)*A_stride+j_A+n+j];
        }
    }
    double * tmp2 = (double *) malloc(n * n * sizeof(double));
    for(int i=0; i < n; i++){
        for(int j=0; j < n; j++){
            tmp2[i*n+j] = B[(i_B+n+i)*B_stride+j_B+j] + B[(i_B+n+i)*B_stride+j_B+n+j];
        }
    }    
    StrassenRecursiveImpl(p6, tmp, tmp2, n, 0, 0, 0, 0, n, n, n);
}
void p7(double *p7, double *A, double *B, int n, int i_A, int j_A, int i_B, int j_B,  int A_stride, int B_stride){
    double * tmp = (double *) malloc(n * n * sizeof(double));
    for(int i=0; i < n; i++){
        for(int j=0; j < n; j++){
            tmp[i*n+j] = A[(i_A+i)*A_stride+j_A+j] - A[(i_A+n+i)*A_stride+j_A+j];
        }
    }
    double * tmp2 = (double *) malloc(n * n * sizeof(double));
    for(int i=0; i < n; i++){
        for(int j=0; j < n; j++){
            tmp2[i*n+j] = B[(i_B+i)*B_stride+j_B+j] + B[(i_B+i)*B_stride+j_B+n+j];
        }
    }    
    StrassenRecursiveImpl(p7, tmp, tmp2, n, 0, 0, 0, 0, n, n, n);
}
void MatAdd(double* A, double *B, int n, int i_A, int j_A, int i_B, int j_B, int A_stride, int B_stride){
    for(int i=0; i < n; i++){
        for(int j=0; j < n; j++){
            A[(i_A+i)*A_stride+j_A+j] += B[(i_B+i)*B_stride+j_B+j];
        }
    }
}
void MatSub(double* A, double *B, int n, int i_A, int j_A, int i_B, int j_B, int A_stride, int B_stride){
    for(int i=0; i < n; i++){
        for(int j=0; j < n; j++){
            A[(i_A+i)*A_stride+j_A+j] -= B[(i_B+i)*B_stride+j_B+j];
        }
    }
}

double * InitMatrix(int n){
    double * p = (double *) malloc(n * n * sizeof(double));

    for(int i=0;i < n*n; i++){
        p[i] = 0.0;
    }
    return p;
}

void StrassenRecursiveImpl(double *a, double* A, double*B, int n, int i_A, int j_A, int i_B, int j_B, int A_stride, int B_stride, int O_stride){
    if(n > BLOCK_SIZE){
        int new_n = n/2;
        double *p1_v = InitMatrix(new_n);
        p1(p1_v, A, B, new_n, i_A, j_A, i_B, j_B, A_stride, B_stride);
        double *p2_v = InitMatrix(new_n);
        p2(p2_v, A, B,  new_n, i_A, j_A, i_B, j_B, A_stride, B_stride);
        double *p3_v = InitMatrix( new_n);
        p3(p3_v, A, B,  new_n, i_A, j_A, i_B, j_B, A_stride, B_stride);
        double *p4_v = InitMatrix( new_n);
        p4(p4_v, A, B,  new_n, i_A, j_A, i_B, j_B, A_stride, B_stride);
        double *p5_v = InitMatrix( new_n);
        p5(p5_v, A, B,  new_n, i_A, j_A, i_B, j_B, A_stride, B_stride);
        double *p6_v = InitMatrix(new_n);
        p6(p6_v, A, B, new_n, i_A, j_A, i_B, j_B, A_stride, B_stride);
        double *p7_v = InitMatrix(new_n);
        p7(p7_v, A, B, new_n, i_A, j_A, i_B, j_B, A_stride, B_stride);
        MatAdd(a, p5_v, new_n, 0, 0, 0, 0, O_stride, new_n);
        MatAdd(a, p4_v, new_n, 0, 0, 0, 0, O_stride, new_n);
        MatAdd(a, p6_v, new_n, 0, 0, 0, 0, O_stride, new_n);
        MatSub(a, p2_v, new_n, 0, 0, 0, 0, O_stride, new_n);
        MatAdd(a, p1_v, new_n, 0, new_n, 0, 0, O_stride, new_n);
        MatAdd(a, p2_v, new_n, 0, new_n, 0, 0, O_stride, new_n);
        MatAdd(a, p3_v, new_n, new_n, 0, 0, 0, O_stride, new_n);
        MatAdd(a, p4_v, new_n, new_n, 0, 0, 0, O_stride, new_n);
        MatAdd(a, p5_v, new_n, new_n, new_n, 0, 0, O_stride, new_n);
        MatAdd(a, p1_v, new_n, new_n, new_n, 0, 0, O_stride, new_n);
        MatSub(a, p3_v, new_n, new_n, new_n, 0, 0, O_stride, new_n);
        MatSub(a, p7_v, new_n, new_n, new_n, 0, 0, O_stride, new_n);
    }else{
        double sum = 0.0;
        int i_stride = (i_A)*A_stride;
        int k_stride = (i_B )*B_stride;
        int o_stride = 0;
        for(int i = 0; i < n ; i++){
            for(int j = 0; j < n; j++){
                sum = 0.0;
                k_stride = (i_B )*B_stride;
                for(int k = 0; k < n; k++){
                    sum += A[i_stride + j_A + k] * B[k_stride + j_B + j];
                    k_stride += B_stride;
                }
                a[o_stride + j] += sum;
            }
            i_stride += A_stride;
            o_stride += O_stride;
        }
    }
}
double *YoursStrassenRecursive(int n, double *A, double *B){
    double *a;
    a = InitMatrix(n);
    time_t start = clock();
    StrassenRecursiveImpl(a, A,B, n, 0, 0, 0, 0, n, n, n);
    time_t end = clock();
    printf("Time %f\n", (double)(end - start)/CLOCKS_PER_SEC);
    return a;    
}
void YoursRecursiveImpl(double* a, double *A, double *B, int n, int i, int j, int k, int stride) {
    // fill your code here, a is your output matrix
    if( n > BLOCK_SIZE){
        int new_n = n/2;
        YoursRecursiveImpl(a, A, B, new_n, i, j, k, stride);
        YoursRecursiveImpl(a, A, B, new_n, i, j, k + new_n, stride);
        YoursRecursiveImpl(a, A, B, new_n, i + new_n, j, k, stride);
        YoursRecursiveImpl(a, A, B, new_n, i + new_n, j, k + new_n, stride);
        YoursRecursiveImpl(a, A, B, new_n, i, j + new_n, k, stride);
        YoursRecursiveImpl(a, A, B, new_n, i, j + new_n, k + new_n, stride);
        YoursRecursiveImpl(a, A, B, new_n, i + new_n, j + new_n, k, stride);
        YoursRecursiveImpl(a, A, B, new_n, i + new_n, j + new_n, k + new_n, stride);
    }else{
        double sum = 0.0;
        int i_stride = i*stride;
        for(int i_ = i; i_ < i + n ; i_++){
            
            for(int j_ = j; j_ < j + n; j_++){
                sum = 0.0;
                for(int k_ = k; k_ < k + n; k_++){
                    sum += A[i_stride + k_] * B[k_*stride + j_];
                }
                
                a[i_stride + j_] += sum;
            }
            i_stride += stride;
        }
    }

}
double * YoursRecursive(int n, double* A, double *B){
    double *a;
    a = (double *) malloc(n * n * sizeof(double));
    time_t start = clock();
    YoursRecursiveImpl(a , A, B, n, 0, 0, 0, n);
    time_t end = clock();
    printf("Time %f\n", (double)(end - start)/CLOCKS_PER_SEC);
    return a;
}
int main(int argc, char *argv[]) {
    srand((unsigned int) time(NULL));
    int n = atoi(argv[1]);
    double *A, *B;
    A = generate(n);
    B = generate(n);

    double *Y;
    Y = (double *) malloc(n * n * sizeof(double));
    Y = generate(n);
    Y = YoursBlocked(n,A,B);

    if (check(Y, A, B, n))
        printf("B TRUE%d\n", 1);
    else
        printf("B FALSE%d\n", 0);

    Y = YoursRecursive(n, A, B);
    if (check(Y, A, B, n))
        printf("R TRUE%d\n", 1);
    else
        printf("R FALSE%d\n", 0);

    Y = YoursStrassenRecursive(n, A, B);
    if (check(Y, A, B, n))
        printf("SR TRUE%d\n", 1);
    else
        printf("SR FALSE%d\n", 0);




    free(A);
    free(B);
    free(Y);
}