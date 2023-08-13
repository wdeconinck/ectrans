#define cufftSafeCall(err) __fftSafeCall(err, __FILE__, __LINE__)
//#include "cuda/cuda_runtime.h"
#include "cufft.h"
#include "stdio.h"
#include "execute_plan_ffth.cuda.h"

#ifdef TRANS_SINGLE
typedef cufftComplex CUDA_DATA_TYPE_COMPLEX;
typedef cufftReal CUDA_DATA_TYPE_REAL;
#else
typedef cufftDoubleComplex CUDA_DATA_TYPE_COMPLEX;
typedef cufftDoubleReal CUDA_DATA_TYPE_REAL;
#endif


    static const char *_cudaGetErrorEnum(cufftResult error)
    {
    switch (error)
    {
    case CUFFT_SUCCESS:
    return "CUFFT_SUCCESS";

    case CUFFT_INVALID_PLAN:
    return "CUFFT_INVALID_PLAN";

    case CUFFT_ALLOC_FAILED:
    return "CUFFT_ALLOC_FAILED";

    case CUFFT_INVALID_TYPE:
    return "CUFFT_INVALID_TYPE";

    case CUFFT_INVALID_VALUE:
    return "CUFFT_INVALID_VALUE";

    case CUFFT_INTERNAL_ERROR:
    return "CUFFT_INTERNAL_ERROR";

    case CUFFT_EXEC_FAILED:
    return "CUFFT_EXEC_FAILED";

    case CUFFT_SETUP_FAILED:
    return "CUFFT_SETUP_FAILED";

    case CUFFT_INVALID_SIZE:
    return "CUFFT_INVALID_SIZE";

    case CUFFT_UNALIGNED_DATA:
    return "CUFFT_UNALIGNED_DATA";

    case CUFFT_INCOMPLETE_PARAMETER_LIST:
    return "CUFFT_INCOMPLETE_PARAMETER_LIST";

    case CUFFT_INVALID_DEVICE:
    return "CUFFT_INVALID_DEVICE";

    case CUFFT_PARSE_ERROR:
    return "CUFFT_PARSE_ERROR";

    case CUFFT_NO_WORKSPACE:
    return "CUFFT_NO_WORKSPACE";

    case CUFFT_NOT_IMPLEMENTED:
    return "CUFFT_NOT_IMPLEMENTED";

    case CUFFT_NOT_SUPPORTED:
    return "CUFFT_NOT_SUPPORTED";
    }

    return "<unknown>";
    }

    inline void __fftSafeCall(cufftResult err, const char *file, const int line)
    {
    if( CUFFT_SUCCESS != err) {
    fprintf(stderr, "CUFFT error at 1\n");
    fprintf(stderr, "CUFFT error in file '%s'\n",__FILE__);
    fprintf(stderr, "CUFFT error at 2\n");
    fprintf(stderr, "CUFFT error at 3\n");
    fprintf(stderr, "CUFFT error %d: %s\nterminating!\n",err,_cudaGetErrorEnum(err)); \
    cudaDeviceReset(); return; \
    }
    }

__global__ void debug(int varId, int N, CUDA_DATA_TYPE_COMPLEX *x) {
    //printf("Hello from GPU\n");
    for (int i = 0; i < N; i++)
    {
        CUDA_DATA_TYPE_COMPLEX a = x[i];
        double b = (double)a.x;
        double c = (double)a.y;
        if (varId == 0) printf("GPU: input[%d]=(%2.4f,%2.4f)\n",i+1,b,c);
        if (varId == 1) printf("GPU: output[%d]=(%2.4f,%2.4f)\n",i+1,b,c);
    }}

__global__ void debugFloat(int varId, int N, CUDA_DATA_TYPE_REAL *x) {
    //printf("Hello from GPU\n");
    for (int i = 0; i < N; i++)
    {
        double a = (double)x[i];
        if (varId == 0) printf("GPU: input[%d]=%2.4f\n",i+1,a);
        if (varId == 1) printf("GPU: output[%d]=%2.4f\n",i+1,a);
    }}


void cudafunction(int ISIGNp, int N, DATA_TYPE *data_in_host, DATA_TYPE *data_out_host, long *iplan)
{
CUDA_DATA_TYPE_COMPLEX *data_in = reinterpret_cast<CUDA_DATA_TYPE_COMPLEX*>(data_in_host);
CUDA_DATA_TYPE_COMPLEX *data_out = reinterpret_cast<CUDA_DATA_TYPE_COMPLEX*>(data_out_host);
cufftHandle* PLANp = reinterpret_cast<cufftHandle*>(iplan);
cufftHandle plan = *PLANp;
int ISIGN = ISIGNp;

//fprintf(stderr, "execute_plan_ffth_c_: transform=%d plan-address = %p\n",ISIGN, plan);

// Check variables on the GPU:
/*int device_count = 0;
cudaGetDeviceCount(&device_count);
for (int i = 0; i < device_count; ++i) {
    cudaSetDevice(i);
    cudaLaunchKernelGGL(debug, dim3(1), dim3(1), 0, 0, 0, N, data_in);
    cudaDeviceSynchronize();
}*/

/*if (cudaDeviceSynchronize() != cudaSuccess){
	fprintf(stderr, "Cuda error: Failed to synchronize\n");
	return;	
}*/

if( ISIGN== -1 ){
  //fprintf(stderr, "Cuda: calling R2C\n");
#ifdef TRANS_SINGLE
  cufftSafeCall(cufftExecR2C(plan, (CUDA_DATA_TYPE_REAL*)data_in, data_out));
#else
  cufftSafeCall(cufftExecD2Z(plan, (CUDA_DATA_TYPE_REAL*)data_in, data_out));
#endif
}
else if( ISIGN== 1){
  //fprintf(stderr, "Cuda: calling C2R\n");
#ifdef TRANS_SINGLE
  cufftSafeCall(cufftExecC2R(plan, data_in, (CUDA_DATA_TYPE_REAL*)data_out));
#else
  cufftSafeCall(cufftExecZ2D(plan, data_in, (CUDA_DATA_TYPE_REAL*)data_out));
#endif
}
else {
  abort();
}

cudaDeviceSynchronize();

/*for (int i = 0; i < device_count; ++i) {
    cudaSetDevice(i);
    cudaLaunchKernelGGL(debugFloat, dim3(1), dim3(1), 0, 0, 1, N, (CUDA_DATA_TYPE_REAL*)data_out);
    cudaDeviceSynchronize();
}*/

//if (cudaDeviceSynchronize() != cudaSuccess){
//	fprintf(stderr, "Cuda error: Failed to synchronize\n");
//	return;	
//}

}
