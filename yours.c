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
void setSubMatrix(double *p, double *A, int n, int i, int j, int stride){
    int i_n = 0;
    int i_stride = i*stride;
    for(int i_=i; i_ < i + n; i_++){
        i_n += n;
        i_stride += stride;
        for(int j_=j; j_ < j + n; j_++){
            p[i_n+j_] = A[i_stride+j_];
        }
    }
}
void StrassenRecursive(double *a, double* A, double*B, int n, int i_A, int j_A, int i_B, int j_B, int A_stride, int B_stride, int O_stride);
//n n/2 
void p1(double *p1, double *A, double *B, int n, int i_A, int j_A, int i_B, int j_B,  int A_stride, int B_stride){
    double * tmp = (double *) malloc(n * n * sizeof(double));
    for(int i=0; i < n; i++){
        for(int j=0; j < n; j++){
            tmp[i*n+j] = B[(i_B+i)*B_stride+j_B+j+n] - B[(i_B+n+i)*B_stride+j_B+j+n];
        }
    }

    StrassenRecursive(p1, A, tmp, n/2, i_A, j_A, 0, 0, A_stride, n, n);
}

void p2(double *p2, double *A, double *B, int n, int i_A, int j_A, int i_B, int j_B,  int A_stride, int B_stride){
    double * tmp = (double *) malloc(n * n * sizeof(double));
    for(int i=0; i < n; i++){
        for(int j=0; j < n; j++){
            tmp[i*n+j] = A[(i_A+i)*A_stride+j_A+j] + A[(i_A+i)*A_stride+j_A+j+n];
        }
    }
    StrassenRecursive(p2, tmp, B, n/2, 0, 0, i_B + n, j_B + n, n, B_stride, n);
}
void p3(double *p3, double *A, double *B, int n, int i_A, int j_A, int i_B, int j_B,  int A_stride, int B_stride){
    double * tmp = (double *) malloc(n * n * sizeof(double));
    for(int i=0; i < n; i++){
        for(int j=0; j < n; j++){
            tmp[i*n+j] = A[(i_A+n+i)*A_stride+j_A+j] + A[(i_A+n+i)*A_stride+j_A+j+n];
        }
    }
    StrassenRecursive(p3, tmp, B, n/2, 0, 0, i_B, j_B, n, B_stride, n);    
}
void p4(double *p4, double *A, double *B, int n, int i_A, int j_A, int i_B, int j_B,  int A_stride, int B_stride){
    double * tmp = (double *) malloc(n * n * sizeof(double));
    for(int i=0; i < n; i++){
        for(int j=0; j < n; j++){
            tmp[i*n+j] = B[(i_B+n+i)*B_stride+j_B+j] - B[(i_B+i)*B_stride+j_B+j];
        }
    }

    StrassenRecursive(p4, A, tmp, n/2, i_A + n, j_A + n, 0, 0, A_stride, n, n);
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
    StrassenRecursive(p5, tmp, tmp2, n/2, 0, 0, 0, 0, n, n, n);
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
    StrassenRecursive(p6, tmp, tmp2, n/2, 0, 0, 0, 0, n, n, n);
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
            tmp2[i*n+j] = B[(i_B+i)*B_stride+j_B+j] - B[(i_B+i)*B_stride+j_B+n+j];
        }
    }    
    StrassenRecursive(p7, tmp, tmp2, n/2, 0, 0, 0, 0, n, n, n);
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
void StrassenRecursive(double *a, double* A, double*B, int n, int i_A, int j_A, int i_B, int j_B, int A_stride, int B_stride, int O_stride){
    if(n > 16){
        double *p1_v = (double *) malloc(n * n * sizeof(double));
        p1(p1_v, A, B, n, i_A, j_A, i_B, j_B, A_stride, B_stride);
        double *p2_v = (double *) malloc(n * n * sizeof(double));
        p2(p2_v, A, B, n, i_A, j_A, i_B, j_B, A_stride, B_stride);
        double *p3_v = (double *) malloc(n * n * sizeof(double));
        p3(p3_v, A, B, n, i_A, j_A, i_B, j_B, A_stride, B_stride);
        double *p4_v = (double *) malloc(n * n * sizeof(double));
        p4(p4_v, A, B, n, i_A, j_A, i_B, j_B, A_stride, B_stride);
        double *p5_v = (double *) malloc(n * n * sizeof(double));
        p5(p5_v, A, B, n, i_A, j_A, i_B, j_B, A_stride, B_stride);
        double *p6_v = (double *) malloc(n * n * sizeof(double));
        p6(p6_v, A, B, n, i_A, j_A, i_B, j_B, A_stride, B_stride);
        double *p7_v = (double *) malloc(n * n * sizeof(double));
        p7(p7_v, A, B, n, i_A, j_A, i_B, j_B, A_stride, B_stride);
        MatAdd(a, p5_v, n, 0, 0, 0, 0, O_stride, n);
        MatAdd(a, p4_v, n, 0, 0, 0, 0, O_stride, n);
        MatAdd(a, p6_v, n, 0, 0, 0, 0, O_stride, n);
        MatSub(a, p2_v, n, 0, 0, 0, 0, O_stride, n);
        MatAdd(a, p1_v, n, 0, n, 0, 0, O_stride, n);
        MatAdd(a, p2_v, n, 0, n, 0, 0, O_stride, n);
        //MatAdd(a, p2_v, n, 0, n, 0, 0, O_stride, n);
        MatAdd(a, p3_v, n, n, 0, 0, 0, O_stride, n);
        MatAdd(a, p4_v, n, n, 0, 0, 0, O_stride, n);
        MatAdd(a, p5_v, n, n, n, 0, 0, O_stride, n);
        MatAdd(a, p1_v, n, n, n, 0, 0, O_stride, n);
        MatAdd(a, p3_v, n, n, n, 0, 0, O_stride, n);
        MatAdd(a, p7_v, n, n, n, 0, 0, O_stride, n);
    }else{
        double sum = 0.0;
        for(int i = 0; i < n ; i++){
            for(int j = 0; j < n; j++){
                sum = 0.0;
                for(int k = 0; k < n; k++){
                    sum += A[(i_A+i)*A_stride + j_A + k] * B[(k + i_B )*B_stride + j_B + j];
                }
                a[i*n + j] += sum;
            }
        }
    }
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
        int i_stride = i*stride;

        for(int i_ = i; i_ < i + n ; i_++){
            i_stride += stride;
            for(int j_ = j; j_ < j + n; j_++){
                sum = 0.0;
                for(int k_ = k; k_ < k + n; k_++){
                    sum += A[i_stride + k_] * B[k_*stride + j_];
                    printf("i,k,j:%d,%d,%d %d\n", i_,k_, j_ ,i_stride);
                }
                a[i_stride + j_] += sum;
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
    Y = (double *) malloc(n * n * sizeof(double));
    //double* a = (double *) malloc(n * n * sizeof(double));
    //int block_inds[] = {0, n, 0, n};

    time_t start = clock();   
    YoursRecursive(Y , A, B, n, 0, 0, 0, n);
    time_t end = clock();

    printf("Time %f\n", (double)(end - start)/CLOCKS_PER_SEC);
    if (check(Y, A, B, n))
        printf("S TRUE%d\n", 1);
    else
        printf("S FALSE%d\n", 0);

    // Y = (double *) malloc(n * n * sizeof(double));
    
    // start = clock();   
    // StrassenRecursive(Y, A,B, n/2, 0, 0, 0, 0, n, n, n);
    // end = clock();
    
    // printf("Time %f\n", (double)(end - start)/CLOCKS_PER_SEC);
    // if (check(Y, A, B, n))
    //     printf("R TRUE%d\n", 1);
    // else
    //     printf("R FALSE%d\n", 0);
    //     
    free(A);
    free(B);
    free(Y);
}