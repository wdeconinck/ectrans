#include <stdio.h>
#include "execute_plan_ffth.cuda.h"

extern "C" {

void execute_plan_ffth_c_(int ISIGNp, int N, DATA_TYPE *data_in, DATA_TYPE *data_out, long *iplan)
{
    cudafunction(ISIGNp,N,data_in,data_out,iplan);
}}
