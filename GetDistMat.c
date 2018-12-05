#include "mex.h"
//[tarMeasMat]=GetDistMat(tarMeas)
void mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[])
{
	double *inData;
	double *outData;
	int M,N;
	int i,j;
    int cnt;
	inData=mxGetPr(prhs[0]);//�����ǵõ���������ָ�룬��MATLAB�����ַ�����Լ����������������ֵ���иĶ�������ʱֵҲ����
	M=mxGetM(prhs[0]);      //�õ����������к���  
	N=mxGetN(prhs[0]);
	plhs[0]=mxCreateDoubleMatrix(N*(N-1)/2,3,mxREAL);   //Ϊ�����������ַ  
    outData=mxGetPr(plhs[0]);     
    cnt=0;
	for(i = 0;i < N;i++)
		for(j = i+1;j < N;j++)
		{
            outData[2*(N*(N-1)/2)+cnt]=(inData[i*2+0]-inData[j*2+0])*(inData[i*2+0]-inData[j*2+0])
                                                        +(inData[i*2+1]-inData[j*2+1])*(inData[i*2+1]-inData[j*2+1]); 
            outData[cnt]=i+1;
            outData[N*(N-1)/2+cnt]=j+1;
            cnt++;
		}
}