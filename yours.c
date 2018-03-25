#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include <printf.h>
#include "handout.h"
#include "time.h"


#define BLOCK_SIZE 2

//these are the implemented methods in 'handout.o' :

//double *generate(int n)
//bool check(double *Yours, double *A, double *B, int n)
//void printMatrix(double *M, int n)
//double *Mult(double *A, double *B, int n)

//the timer should be implemented in the YoursBlocked & YoursRecursive function and printed out in a format like "TIME: 0.000000 seconds"
double * InitMatrix(int n){
    double * p = (double *) malloc(n * n * sizeof(double));

    for(int i=0;i < n*n; i++){
        p[i] = 0.0;
    }
    return p;
}
double *Naive(int n, double* A, double *B){
    double *a;
    a = InitMatrix(n);

    time_t start = clock();
    unsigned int block_size = BLOCK_SIZE;
    double sum = 0.0;
    for(int i=0; i < n; i++){
        for(int j=0; j < n; j++){
            sum = 0.0;
            for(int k=0;k < n; k++){
                sum += A[i*n+k] * B[k*n+j];
            }
            a[i*n+j] += sum;
        }
    }
    time_t end = clock();
    printf("Time %f\n", (double)(end - start)/CLOCKS_PER_SEC);
    return a;        
}

double *YoursBlocked(int n, double *A, double *B) {
    double *a;
    a = InitMatrix(n);

    time_t start = clock();
    unsigned int i_block_size = BLOCK_SIZE;
    unsigned int j_block_size = BLOCK_SIZE;
    unsigned int k_block_size = BLOCK_SIZE;
    if(n < i_block_size){
        i_block_size = n;
        j_block_size = n;
        k_block_size = n;
    }
    double sum = 0.0;
    int i_block_num = (n)/i_block_size + 1;
    int j_block_num = (n)/j_block_size + 1;
    int k_block_num = (n)/k_block_size + 1;
    int left_size = n - (j_block_num - 1) * BLOCK_SIZE;
    for(int j=0; j < j_block_num; j++){
        int j_all = j*j_block_size;
        j_block_size = j < j_block_num - 1 ? BLOCK_SIZE :left_size;
        for(int k=0; k < k_block_num; k++){
            int k_all = k*k_block_size;
            k_block_size = k < k_block_num - 1 ? BLOCK_SIZE :left_size;
            for(int i=0; i < i_block_num; i++){
                int i_all = i*i_block_size;
                i_block_size = i < i_block_num - 1 ? BLOCK_SIZE :left_size;
                for(int l=0; l < j_block_size; l++){
                    int B_all_ind_ = j_all+l + k_all*n;
                    for(int t=0; t < i_block_size; t++){
                        sum = 0.0;
                        int A_all_ind = (i_all+t)*n + k_all;
                        int m = 0;
                        int B_all_ind =  B_all_ind_;
                        while(m < k_block_size){
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
void StrassenRecursiveImpl(double *a, double* A, double*B, int n,int pad_i_A, int pad_j_A, int pad_i_B, int pad_j_B, int i_A, int j_A, int i_B, int j_B, int A_stride, int B_stride, int O_stride);
//n n/2 
//
void p1(double *p1, double *A, double *B, int n, int pad, int i_A, int j_A, int i_B, int j_B,  int A_stride, int B_stride){
    double * tmp = InitMatrix(n);

    double f,h;
    int true_n = n-pad;
    for(int i=0; i < n; i++){
        for(int j=0; j < n; j++){
            f = j < true_n ? B[(i_B+i)*B_stride+j_B+j+n]:0;
            h = j < true_n && i < true_n ? B[(i_B+n+i)*B_stride+j_B+j+n]:0;
            tmp[i*n+j] = f - h;
        }
    }
    StrassenRecursiveImpl(p1, A, tmp, n, 0, 0, 0, 0, i_A, j_A, 0, 0, A_stride, n, n);
    

    
}

void p2(double *p2, double *A, double *B, int n, int pad, int i_A, int j_A, int i_B, int j_B,  int A_stride, int B_stride){
    double * tmp = InitMatrix(n);
    double b;
    int true_n = n-pad;
    for(int i=0; i < n; i++){
        for(int j=0; j < n; j++){
            b = j < true_n ? A[(i_A+i)*A_stride+j_A+j+n]:0;
            tmp[i*n+j] = A[(i_A+i)*A_stride+j_A+j] + b;
        }
    }
    StrassenRecursiveImpl(p2, tmp, B, n, 0, 0, pad, pad,  0, 0, i_B + n, j_B + n, n, B_stride, n);
}
void p3(double *p3, double *A, double *B, int n, int pad, int i_A, int j_A, int i_B, int j_B,  int A_stride, int B_stride){
    double * tmp = InitMatrix(n);
    double c, d;
    int true_n = n-pad;
    for(int i=0; i < n; i++){
        for(int j=0; j < n; j++){
            c = i < true_n ? A[(i_A+n+i)*A_stride+j_A+j]:0;
            d = i < true_n && j < true_n ? A[(i_A+n+i)*A_stride+j_A+j+n]:0;
            tmp[i*n+j] = c + d;
        }
    }
    StrassenRecursiveImpl(p3, tmp, B, n, 0, 0, 0, 0, 0, 0,  i_B, j_B, n, B_stride, n);    
}
void p4(double *p4, double *A, double *B, int n, int pad, int i_A, int j_A, int i_B, int j_B,  int A_stride, int B_stride){
    double * tmp = InitMatrix(n);
    double g;
    int true_n = n-pad;
    for(int i=0; i < n; i++){
        for(int j=0; j < n; j++){
            g = i < true_n ? B[(i_B+n+i)*B_stride+j_B+j]:0;
            //e = i < true_n && j < true_n ? B[(i_B+i)*B_stride+j_B+j]:0;
            tmp[i*n+j] = g - B[(i_B+i)*B_stride+j_B+j];
        }
    }

    StrassenRecursiveImpl(p4, A, tmp, n, pad, pad, 0, 0, i_A + n, j_A + n, 0, 0, A_stride, n, n);
}
void p5(double *p5, double *A, double *B, int n, int pad, int i_A, int j_A, int i_B, int j_B,  int A_stride, int B_stride){
    double * tmp = InitMatrix(n);
    double d;
    int true_n = n-pad;
    for(int i=0; i < n; i++){
        for(int j=0; j < n; j++){
            d = i < true_n && j < true_n ? A[(i_A+n+i)*A_stride+j_A+n+j]:0;
            tmp[i*n+j] = A[(i_A+i)*A_stride+j_A+j] + d;
        }
    }
    double * tmp2 = InitMatrix(n);
    double h;
    for(int i=0; i < n; i++){
        for(int j=0; j < n; j++){
            h = i < true_n && j < true_n ? B[(i_B+n+i)*B_stride+j_B+n+j]:0;
            tmp2[i*n+j] = B[(i_B+i)*B_stride+j_B+j] + h;
        }
    }    
    StrassenRecursiveImpl(p5, tmp, tmp2, n, 0, 0, 0, 0, 0, 0, 0, 0, n, n, n);
}
void p6(double *p6, double *A, double *B, int n, int pad, int i_A, int j_A, int i_B, int j_B,  int A_stride, int B_stride){
    double * tmp = InitMatrix(n);
    int true_n = n-pad;
    double b, d;
    for(int i=0; i < n; i++){
        for(int j=0; j < n; j++){
            b = j < true_n ? A[(i_A+i)*A_stride+j_A+n+j]:0;
            d = i < true_n && j < true_n ? A[(i_A+n+i)*A_stride+j_A+n+j]:0;
            tmp[i*n+j] = b - d;
        }
    }
    double * tmp2 = InitMatrix(n);
    double g,h;
    for(int i=0; i < n; i++){
        for(int j=0; j < n; j++){
            g = i < true_n ? B[(i_B+n+i)*B_stride+j_B+j]:0;
            h = i < true_n && j < true_n ? B[(i_B+n+i)*B_stride+j_B+n+j]:0;
            tmp2[i*n+j] = g + h;
        }
    }    
    StrassenRecursiveImpl(p6, tmp, tmp2, n, 0, 0, 0, 0, 0, 0, 0, 0, n, n, n);
}
void p7(double *p7, double *A, double *B, int n, int pad, int i_A, int j_A, int i_B, int j_B,  int A_stride, int B_stride){
    double * tmp = InitMatrix(n);
    int true_n = n-pad;
    double c;
    for(int i=0; i < n; i++){
        for(int j=0; j < n; j++){
            c = i < true_n ? A[(i_A+n+i)*A_stride+j_A+j]:0;
            tmp[i*n+j] = A[(i_A+i)*A_stride+j_A+j] - c;
        }
    }
    double * tmp2 = InitMatrix(n);
    double  f;
    for(int i=0; i < n; i++){
        for(int j=0; j < n; j++){
            f = j < true_n ? B[(i_B+i)*B_stride+j_B+n+j]:0;
            tmp2[i*n+j] = B[(i_B+i)*B_stride+j_B+j] + f;
        }
    }    
    StrassenRecursiveImpl(p7, tmp, tmp2, n, 0, 0, 0, 0, 0, 0, 0, 0, n, n, n);
}
void MatAdd(double* A, double *B, int n_i, int n_j, int i_A, int j_A, int i_B, int j_B, int A_stride, int B_stride){
    for(int i=0; i < n_i; i++){
        for(int j=0; j < n_j; j++){
            A[(i_A+i)*A_stride+j_A+j] += B[(i_B+i)*B_stride+j_B+j];
        }
    }
}
void MatSub(double* A, double *B, int n_i, int n_j, int i_A, int j_A, int i_B, int j_B, int A_stride, int B_stride){
    for(int i=0; i < n_i; i++){
        for(int j=0; j < n_j; j++){
            A[(i_A+i)*A_stride+j_A+j] -= B[(i_B+i)*B_stride+j_B+j];
        }
    }
}



void StrassenRecursiveImpl(double *O, double* A, double*B, int n, int pad_i_A, int pad_j_A, int pad_i_B, int pad_j_B, int i_A, int j_A, int i_B, int j_B, int A_stride, int B_stride, int O_stride){
    if(n > BLOCK_SIZE){

        int new_n = (n+1)/2;
        int pad = new_n - n/2;
        double *p1_v = InitMatrix(new_n);
        p1(p1_v, A, B, new_n, pad, i_A, j_A, i_B, j_B, A_stride, B_stride);
        double *p2_v = InitMatrix(new_n);
        p2(p2_v, A, B,  new_n, pad, i_A, j_A, i_B, j_B, A_stride, B_stride);
        double *p3_v = InitMatrix( new_n);
        p3(p3_v, A, B,  new_n, pad, i_A, j_A, i_B, j_B, A_stride, B_stride);
        double *p4_v = InitMatrix( new_n);
        p4(p4_v, A, B,  new_n, pad, i_A, j_A, i_B, j_B, A_stride, B_stride);
        double *p5_v = InitMatrix( new_n);
        p5(p5_v, A, B,  new_n, pad, i_A, j_A, i_B, j_B, A_stride, B_stride);
        double *p6_v = InitMatrix(new_n);
        p6(p6_v, A, B, new_n, pad, i_A, j_A, i_B, j_B, A_stride, B_stride);
        double *p7_v = InitMatrix(new_n);
        p7(p7_v, A, B, new_n, pad, i_A, j_A, i_B, j_B, A_stride, B_stride);
        MatAdd(O, p5_v, new_n, new_n, 0, 0, 0, 0, O_stride, new_n);
        MatAdd(O, p4_v, new_n, new_n, 0, 0, 0, 0, O_stride, new_n);
        MatAdd(O, p6_v, new_n, new_n, 0, 0, 0, 0, O_stride, new_n);
        MatSub(O, p2_v, new_n, new_n, 0, 0, 0, 0, O_stride, new_n);
        MatAdd(O, p1_v, new_n, new_n-pad, 0, new_n, 0, 0, O_stride, new_n);
        MatAdd(O, p2_v, new_n, new_n-pad, 0, new_n, 0, 0, O_stride, new_n);
        MatAdd(O, p3_v, new_n-pad, new_n,  new_n, 0, 0, 0, O_stride, new_n);
        MatAdd(O, p4_v, new_n-pad, new_n, new_n, 0, 0, 0, O_stride, new_n);
        MatAdd(O, p5_v, new_n-pad, new_n-pad, new_n, new_n, 0, 0, O_stride, new_n);
        MatAdd(O, p1_v, new_n-pad, new_n-pad, new_n, new_n, 0, 0, O_stride, new_n);
        MatSub(O, p3_v, new_n-pad, new_n-pad, new_n, new_n, 0, 0, O_stride, new_n);
        MatSub(O, p7_v, new_n-pad, new_n-pad, new_n, new_n, 0, 0, O_stride, new_n);
    }else{
        double sum = 0.0;
        int i_stride = (i_A)*A_stride;
        int k_stride = (i_B )*B_stride;
        int o_stride = 0;
        int n_pad_i_A = n - pad_i_A;
        int n_pad_j_A = n - pad_j_A;
        int n_pad_i_B = n - pad_i_B;
        int n_pad_j_B = n - pad_i_B;
        double a,b;
        for(int i = 0; i < n ; i++){
            for(int j = 0; j < n; j++){
                sum = 0.0;
                k_stride = (i_B )*B_stride;
                for(int k = 0; k < n; k++){
                    a = i < n_pad_i_A && k < n_pad_j_A ? A[i_stride + j_A + k]:0;
                    b = j < n_pad_j_B && k < n_pad_i_B ? B[k_stride + j_B + j]:0;
                    sum += a * b;
                    k_stride += B_stride;
                    //printf("i,j,k,l:%d,%d,%d,%d",i_A+i, j_A+k, i_B+k, j_B+j);
                }
                
                    O[o_stride + j] += sum;
                
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
    int pad = (n+1)/2 - n/2;
    StrassenRecursiveImpl(a, A, B, n, pad, pad, pad, pad, 0, 0, 0, 0, n, n, n);
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
    a = InitMatrix(n);

    time_t start = clock();
    YoursRecursiveImpl(a, A, B, n, 0, 0, 0, n);
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
    Y = Naive(n,A,B);
    printf("N\n");
    for(int i=0; i<n ;i++){
        printf("\n");
        for(int j=0; j<n; j++)
        printf("%5f ", Y[i*n+j]);
    }
    printf("\n");
    printf("\n");
    if (check(Y, A, B, n))
        printf("N TRUE%d\n", 1);
    else
        printf("N FALSE%d\n", 0);

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
    printf("sr\n");
    for(int i=0; i<n ;i++){
        printf("\n");
        for(int j=0; j<n; j++)
        printf("%5f ", Y[i*n+j]);
    }
    printf("\n");
    if (check(Y, A, B, n))
        printf("SR TRUE%d\n", 1);
    else
        printf("SR FALSE%d\n", 0);




    free(A);
    free(B);
    free(Y);
}