#include <stdlib.h>
#include <stdio.h>
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
double * InitMatrix(int n){
    double * p = (double *) malloc(n * n * sizeof(double));

    for(int i=0;i < n*n; i++){
        p[i] = 0.0;
    }
    return p;
}
void printMat(double *x , int n){
    for(int i=0;i < n;i++){
        for(int j=0;j<n;j++){
            printf("%f ", x[i*n+j]);
        }
        printf("\n");
    }
    printf("\n");
}
int getPadLen(int n, int Q) {
    int cnt = 0;
    int tmp = n;
    while(tmp > Q) {
        cnt++;
        tmp /= 2;
    }
    cnt++;
    // result should be smallest value such that:
    // result >= actual_size AND
    // result % (1<<cnt) == 0

    if (n % (1<<cnt) == 0) {
        return n;
    } else {
        return n + (1<<cnt) - n % (1<<cnt);
    }
}
double *Naive(int n, double* A, double *B){
    double *a;
    a = InitMatrix(n);

    struct timespec time_start={0, 0}, time_end={0, 0};
    clock_gettime(CLOCK_REALTIME, &time_start);
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
    clock_gettime(CLOCK_REALTIME, &time_end);
    printf("Time %llu.%09llu s \t",(long long unsigned)time_end.tv_sec-time_start.tv_sec, (long long unsigned)time_end.tv_nsec-time_start.tv_nsec);
    return a;        
}

double *YoursBlocked(int n, double *A, double *B, int BLOCK_SIZE) {
    double *a;
    a = InitMatrix(n);

    struct timespec time_start={0, 0}, time_end={0, 0};
    clock_gettime(CLOCK_REALTIME, &time_start);
    if(n < BLOCK_SIZE){
        BLOCK_SIZE = n;
    }
    unsigned int i_block_size = BLOCK_SIZE;
    unsigned int j_block_size = BLOCK_SIZE;
    unsigned int k_block_size = BLOCK_SIZE;
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
    
    clock_gettime(CLOCK_REALTIME, &time_end);
    printf("Time :%llu.%09llu ns \t",time_end.tv_sec-time_start.tv_sec, time_end.tv_nsec-time_start.tv_nsec);
// fill your code here, a is your output matrix
    return a;
}
void StrassenRecursiveImpl(double *a, double* A, double*B, int n, int i_A, int j_A, int i_B, int j_B, int A_stride, int B_stride, int O_stride, int b);
//n n/2 
//
void p1(double *p1, double *A, double *B, int n, int i_A, int j_A, int i_B, int j_B,  int A_stride, int B_stride, int b){
    double * tmp = InitMatrix(n);
    int f_i_stride = (i_B)*B_stride;
    int h_i_stride = (i_B+n)*B_stride;
    for(int i=0; i < n; i++){
        for(int j=0; j < n; j++){
            tmp[i*n+j] = B[f_i_stride+j_B+j+n] - B[h_i_stride+j_B+j+n];
        }
        f_i_stride += B_stride;
        h_i_stride += B_stride;
    }
    StrassenRecursiveImpl(p1, A, tmp, n,  i_A, j_A, 0, 0, A_stride, n, n, b);
    
}

void p2(double *p2, double *A, double *B, int n,  int i_A, int j_A, int i_B, int j_B,  int A_stride, int B_stride, int b){
    double * tmp = InitMatrix(n);
    int a_stride = (i_A)*A_stride;
    //int B_stride = (i_A)*A_stride
    for(int i=0; i < n; i++){
        for(int j=0; j < n; j++){
            tmp[i*n+j] = A[a_stride+j_A+j] + A[a_stride+j_A+j+n];
        }
        a_stride += A_stride;
    }
    StrassenRecursiveImpl(p2, tmp, B, n,   0, 0, i_B + n, j_B + n, n, B_stride, n, b);
}
void p3(double *p3, double *A, double *B, int n,  int i_A, int j_A, int i_B, int j_B,  int A_stride, int B_stride, int b){
    double * tmp = InitMatrix(n);
    int b_stride = (i_A+n)*A_stride;
    for(int i=0; i < n; i++){
        for(int j=0; j < n; j++){
            tmp[i*n+j] = A[b_stride+j_A+j] + A[b_stride+j_A+j+n];
        }
        b_stride += A_stride;
    }
    StrassenRecursiveImpl(p3, tmp, B, n, 0, 0,  i_B, j_B, n, B_stride, n, b);    
}
void p4(double *p4, double *A, double *B, int n,  int i_A, int j_A, int i_B, int j_B,  int A_stride, int B_stride, int b){
    double * tmp = InitMatrix(n);
    int g_i_stride = (i_B+n)*B_stride;
    int e_i_stride = (i_B)*B_stride;    
    for(int i=0; i < n; i++){
        for(int j=0; j < n; j++){
            tmp[i*n+j] = B[g_i_stride+j_B+j] - B[e_i_stride+j_B+j];
        }
        g_i_stride += B_stride;
        e_i_stride += B_stride;
    }

    StrassenRecursiveImpl(p4, A, tmp, n,  i_A + n, j_A + n, 0, 0, A_stride, n, n, b);
}
void p5(double *p5, double *A, double *B, int n,  int i_A, int j_A, int i_B, int j_B,  int A_stride, int B_stride, int b){
    double * tmp = InitMatrix(n);
    int a_i_stride = (i_A)*A_stride;
    int d_i_stride = (i_A+n)*A_stride;
    for(int i=0; i < n; i++){
        for(int j=0; j < n; j++){
            
            tmp[i*n+j] = A[a_i_stride+j_A+j] + A[d_i_stride+j_A+n+j];
        }
        a_i_stride += A_stride;
        d_i_stride += A_stride;
    }
    double * tmp2 = InitMatrix(n);
    int e_i_stride = (i_B)*B_stride;
    int h_i_stride = ((i_B+n)*B_stride);
    for(int i=0; i < n; i++){
        for(int j=0; j < n; j++){
            
            tmp2[i*n+j] = B[e_i_stride+j_B+j] + B[h_i_stride+j_B+n+j];
        }
        e_i_stride += B_stride;
        h_i_stride += B_stride;
    }    
    StrassenRecursiveImpl(p5, tmp, tmp2, n,  0, 0, 0, 0, n, n, n, b);
}
void p6(double *p6, double *A, double *B, int n,  int i_A, int j_A, int i_B, int j_B,  int A_stride, int B_stride, int b){
    double * tmp = InitMatrix(n);
    int b_i_stride = (i_A)*A_stride;
    int d_i_stride = (i_A+n)*A_stride;
    for(int i=0; i < n; i++){
        for(int j=0; j < n; j++){

            tmp[i*n+j] = A[b_i_stride+j_A+n+j] - A[d_i_stride+j_A+n+j];
        }
        b_i_stride += A_stride;
        d_i_stride += A_stride;
    }
    double * tmp2 = InitMatrix(n);
    int g_i_stride = (i_B+n)*B_stride;
    
    for(int i=0; i < n; i++){
        for(int j=0; j < n; j++){

            tmp2[i*n+j] = B[g_i_stride+j_B+j] + B[g_i_stride+j_B+n+j];
        }
        g_i_stride += B_stride;
    }    
    StrassenRecursiveImpl(p6, tmp, tmp2, n,  0, 0, 0, 0, n, n, n, b);
}
void p7(double *p7, double *A, double *B, int n,  int i_A, int j_A, int i_B, int j_B,  int A_stride, int B_stride, int b){
    double * tmp = InitMatrix(n);
    int a_i_stride = (i_A)*A_stride;
    int c_i_stride = (i_A+n)*A_stride;
    for(int i=0; i < n; i++){
        for(int j=0; j < n; j++){
            
            tmp[i*n+j] = A[a_i_stride+j_A+j] - A[c_i_stride+j_A+j];
        }
        a_i_stride += A_stride;
        c_i_stride += A_stride;
    }
    double * tmp2 = InitMatrix(n);
    int e_i_stride = (i_B)*B_stride;
    for(int i=0; i < n; i++){
        for(int j=0; j < n; j++){
            tmp2[i*n+j] = B[e_i_stride+j_B+j] + B[e_i_stride+j_B+n+j];
        }
        e_i_stride += B_stride;
    }    
    StrassenRecursiveImpl(p7, tmp, tmp2, n, 0, 0, 0, 0, n, n, n, b);
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



void StrassenRecursiveImpl(double *O, double* A, double*B, int n,  int i_A, int j_A, int i_B, int j_B, int A_stride, int B_stride, int O_stride, int b){
    if(n > b){

        int new_n = n/2;
        double *p1_v = InitMatrix(new_n);
        p1(p1_v, A, B, new_n, i_A, j_A, i_B, j_B, A_stride, B_stride, b);
        double *p2_v = InitMatrix(new_n);
        p2(p2_v, A, B,  new_n, i_A, j_A, i_B, j_B, A_stride, B_stride, b);
        double *p3_v = InitMatrix( new_n);
        p3(p3_v, A, B,  new_n, i_A, j_A, i_B, j_B, A_stride, B_stride, b);
        double *p4_v = InitMatrix( new_n);
        p4(p4_v, A, B,  new_n, i_A, j_A, i_B, j_B, A_stride, B_stride, b);
        double *p5_v = InitMatrix( new_n);
        p5(p5_v, A, B,  new_n, i_A, j_A, i_B, j_B, A_stride, B_stride, b);
        double *p6_v = InitMatrix(new_n);
        p6(p6_v, A, B, new_n, i_A, j_A, i_B, j_B, A_stride, B_stride, b);
        double *p7_v = InitMatrix(new_n);
        p7(p7_v, A, B, new_n, i_A, j_A, i_B, j_B, A_stride, B_stride, b);

        MatAdd(O, p5_v, new_n,  0, 0, 0, 0, O_stride, new_n);
        MatAdd(O, p4_v, new_n,  0, 0, 0, 0, O_stride, new_n);
        MatAdd(O, p6_v, new_n,  0, 0, 0, 0, O_stride, new_n);
        MatSub(O, p2_v, new_n,  0, 0, 0, 0, O_stride, new_n);
        MatAdd(O, p1_v, new_n,  0, new_n, 0, 0, O_stride, new_n);
        MatAdd(O, p2_v, new_n,  0, new_n, 0, 0, O_stride, new_n);
        MatAdd(O, p3_v, new_n,  new_n, 0, 0, 0, O_stride, new_n);
        MatAdd(O, p4_v, new_n, new_n, 0, 0, 0, O_stride, new_n);
        MatAdd(O, p5_v, new_n,  new_n, new_n, 0, 0, O_stride, new_n);
        MatAdd(O, p1_v, new_n,  new_n, new_n, 0, 0, O_stride, new_n);
        MatSub(O, p3_v, new_n,  new_n, new_n, 0, 0, O_stride, new_n);
        MatSub(O, p7_v, new_n,  new_n, new_n, 0, 0, O_stride, new_n);
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
                O[o_stride + j] += sum;
                
            }
            i_stride += A_stride;
            o_stride += O_stride;
        }
    }
}
double *PadMat(double *A, int n, int pad){
    double * r = InitMatrix(n+pad);
    for(int i=0; i < n; i++){
        for(int j=0; j < n; j++){
            r[i*(n+pad)+j] = A[i*n+j];
        }
    }
    return r;
}
double *YoursRecursive(int n, double *A, double *B, int b){
    double *a;
    int padLen = getPadLen(n, b);
    a = InitMatrix(padLen);
    //printf("pad %d", padLen);

    struct timespec time_start={0, 0}, time_end={0, 0};
    clock_gettime(CLOCK_REALTIME, &time_start);
   
   
    
    double *padA = PadMat(A, n, padLen-n);
    double *padB = PadMat(B, n, padLen-n);
    StrassenRecursiveImpl(a, padA, padB, padLen,  0, 0, 0, 0, padLen, padLen, padLen, b);
    double *r = InitMatrix(n);
    for(int i=0;i<n;i++){
        for(int j=0;j <n;j++){
            r[i*n+j] = a[i*(padLen)+j];
        }
    }
    time_t end = clock();

    clock_gettime(CLOCK_REALTIME, &time_end);
    printf("Time %llu.%09llu s \t",time_end.tv_sec-time_start.tv_sec, time_end.tv_nsec-time_start.tv_nsec);
    return r;    
}
int test(int block_size, int n){

    double *A, *B;
    A = generate(n);
    B = generate(n);

    double *Y;
    Y = (double *) malloc(n * n * sizeof(double));
    Y = generate(n);
    Y = Naive(n,A,B);

    if (check(Y, A, B, n))
        printf("N TRUE%d\n", 1);
    else
        printf("N FALSE%d\n", 0);

    Y = YoursBlocked(n,A,B, block_size);

    if (check(Y, A, B, n))
        printf("B TRUE%d\n", 1);
    else
        printf("B FALSE%d\n", 0);

    Y = YoursRecursive(n, A, B, block_size);

    if (check(Y, A, B, n))
        printf("R TRUE%d\n", 1);
    else
        printf("R FALSE%d\n", 0);




    free(A);
    free(B);
    free(Y);    
    return 0;

}
int main(int argc, char *argv[]) {
    for(int i = 0; i < 256; i+=16){
        for(int j = 1000; j < 1099; j += 5){
            test(i, j);
        }
    }
    return 0;
}

// int main(int argc, char *argv[]) {
//     int BLOCK_SIZE = 128;
//     srand((unsigned int) time(NULL));
//     int n = atoi(argv[1]);
//     double *A, *B;
//     A = generate(n);
//     B = generate(n);

//     double *Y;
//     Y = (double *) malloc(n * n * sizeof(double));
//     Y = generate(n);
//     // Y = Naive(n,A,B);

//     // if (check(Y, A, B, n))
//     //     printf("N TRUE%d\n", 1);
//     // else
//     //     printf("N FALSE%d\n", 0);

//     Y = YoursBlocked(n, A, B, BLOCK_SIZE);

//     if (check(Y, A, B, n))
//         printf("B TRUE%d\n", 1);
//     else
//         printf("B FALSE%d\n", 0);

//     // Y = YoursRecursive(n, A, B);
//     // if (check(Y, A, B, n))
//     //     printf("R TRUE%d\n", 1);
//     // else
//     //     printf("R FALSE%d\n", 0);

//     Y = YoursRecursive(n, A, B, BLOCK_SIZE);

//     if (check(Y, A, B, n))
//         printf("R TRUE%d\n", 1);
//     else
//         printf("R FALSE%d\n", 0);




//     free(A);
//     free(B);
//     free(Y);
// }