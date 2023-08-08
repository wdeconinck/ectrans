//
// Wrapper for cublasDgemm function. 
//
// Alan Gray, NVIDIA
//

#include <stdio.h>
#include "cublas_v2.h" 


bool alreadyAllocated_dgemm=false;
bool alreadyAllocated_dgemm_handle=false;

double **d_Aarray;
double **d_Barray;
double **d_Carray;

double **Aarray;
double **Barray;
double **Carray;

cublasHandle_t handle_dgemm;	

extern "C" void cublasDgemmBatched_wrapper (char transa, char transb, int m, int n,int k, double alpha, const double *A, int lda, int tda, const double *B, int ldb, int tdb, double beta, double *C, int ldc, int tdc, int batchCount)
{


  // printf("CUBLAS m=%d,n=%d,k=%d,batchcount=%d\n",m,n,k,batchCount);
    cublasStatus_t stat;

 
  cublasOperation_t op_t1=CUBLAS_OP_N, op_t2=CUBLAS_OP_N;

  if (transa=='T' || transa=='t')	
    op_t1=CUBLAS_OP_T;

  if (transb=='T' || transb=='t')
    op_t2=CUBLAS_OP_T;


  //double **Aarray = (double**) malloc(batchCount*sizeof(double*));
  //double **Barray = (double**) malloc(batchCount*sizeof(double*));
  //double **Carray = (double**) malloc(batchCount*sizeof(double*));



  if (!alreadyAllocated_dgemm_handle){
     stat = cublasCreate(&handle_dgemm);
     if (stat != CUBLAS_STATUS_SUCCESS) {
        printf ("CUBLAS initialization failed\n");
        //return EXIT_FAILURE;
    }
  }
  alreadyAllocated_dgemm_handle=true;

  if (!alreadyAllocated_dgemm){
    cudaError_t errcm1 = cudaMallocHost(&Aarray,batchCount*sizeof(double*));
    cudaError_t errcm2 = cudaMallocHost(&Barray,batchCount*sizeof(double*));
    cudaError_t errcm3 = cudaMallocHost(&Carray,batchCount*sizeof(double*));
        
    cudaError_t errcm4 = cudaMalloc(&d_Aarray,batchCount*sizeof(double*));
    cudaError_t errcm5 = cudaMalloc(&d_Barray,batchCount*sizeof(double*));
    cudaError_t errcm6 = cudaMalloc(&d_Carray,batchCount*sizeof(double*));
   }
  alreadyAllocated_dgemm=true;

  int i;
  for(i=0;i<batchCount;i++){
    Aarray[i]=(double*) &(A[i*lda*tda]);
    Barray[i]=(double*) &(B[i*ldb*tdb]);
    Carray[i]=(double*) &(C[i*ldc*tdc]);
  }

  cudaError_t err1 = cudaMemcpy(d_Aarray,Aarray,batchCount*sizeof(double*),cudaMemcpyHostToDevice);
  cudaError_t err2 = cudaMemcpy(d_Barray,Barray,batchCount*sizeof(double*),cudaMemcpyHostToDevice);
  cudaError_t err3 = cudaMemcpy(d_Carray,Carray,batchCount*sizeof(double*),cudaMemcpyHostToDevice);
  cudaDeviceSynchronize();


  cublasDgemmBatched(handle_dgemm,op_t1,op_t2,m,n,k,&alpha,(const double**) d_Aarray,lda, (const double**) d_Barray,ldb,&beta,(double**) d_Carray,ldc,batchCount);

  cudaDeviceSynchronize();
  
  //cudaFree(Aarray);
  //cudaFree(Barray);
  //cudaFree(Carray);
  
  //cudaFree(d_Aarray);
  //cudaFree(d_Barray);
  //cudaFree(d_Carray);
  //cublasDestroy(handle_dgemm);
  
  
}

extern "C" void cublasDgemmStridedBatched_wrapper (char transa, char transb, int m, int n,int k, double alpha, const double *A, int lda, long long tda, const double *B, int ldb, long long tdb, double beta, double *C, int ldc, long long tdc, int batchCount)
{


  // printf("CUBLAS m=%d,n=%d,k=%d,batchcount=%d\n",m,n,k,batchCount);

 
  cublasOperation_t op_t1=CUBLAS_OP_N, op_t2=CUBLAS_OP_N;

  if (transa=='T' || transa=='t')       
    op_t1=CUBLAS_OP_T;

  if (transb=='T' || transb=='t')       
    op_t2=CUBLAS_OP_T;

  if (!alreadyAllocated_dgemm_handle){
    cublasCreate(&handle_dgemm);
    alreadyAllocated_dgemm_handle=true;
  }
  cublasDgemmStridedBatched(handle_dgemm,op_t1,op_t2,m,n,k,&alpha,(const double *) A,lda,tda, (const double *) B,ldb,tdb,&beta,(double *) C,ldc,tdc,batchCount);

}

extern "C" void cublasDgemmBatched_finalize ()
{



  if (alreadyAllocated_dgemm){
  
    cudaFree(Aarray);
    cudaFree(Barray);
    cudaFree(Carray);
    
    cudaFree(d_Aarray);
    cudaFree(d_Barray);
    cudaFree(d_Carray);
    if (alreadyAllocated_dgemm_handle){
      cublasDestroy(handle_dgemm);
    }
    alreadyAllocated_dgemm_handle=false;

  }
  alreadyAllocated_dgemm=false;
  
}
