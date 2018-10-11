#include<iostream>
#include<fstream>
#include<complex>
using namespace std;

int filterNum = 24;
int sampleRate = 48000;
 
#define Win_Time 0.01//把10ms里的所有点作为一个点分析 
#define Hop_Time 0.005//每隔5ms分一次帧 
#define Pi 3.1415927

//1.预加重
double* pre_emphasizing(double *sample, int len, double factor)
{
	double *Sample = new double[len];
	Sample[0] = sample[0];
	for(int i = 1; i < len; i++)
	{
		//预加重过程
		Sample[i] = sample[i] - factor * sample[i - 1]; 
	}

	ofstream filesample("./txt/sample.txt");
    ofstream filePreemphasized("./txt/Preemphasized.txt");
	for(int i = 0; i < len; i++)
	{
		filesample << sample[i] << endl; 
		filePreemphasized << Sample[i] << endl; 
	}

	return Sample;
} 

void Hamming( double *hamWin, int hamWinSize )
{
	for (int i = 0; i < hamWinSize; i++)
	{
		hamWin[i] = (double)(0.54 - 0.46 * cos(2 * Pi * (double)i / ((double)hamWinSize - 1) ));
	}
}

static int count[14] = {1,2,4,8,16,32,64,128,256,512,1024,2048,4096,8192};
//计算每一帧的功率谱 
void mfccFFT(double *frameSample, double *FFTSample, int frameSize, int pos)
{
    //对分帧加窗后的各帧信号进行FFT变换得到各帧的频谱
	//并对语音信号的频谱取模平方得到语音信号的功率谱
    double *dataR;
    double *dataI;
	dataR = new double[frameSize];
	dataI = new double[frameSize];
	for(int i = 0; i < frameSize; i++)
	{
		dataR[i] = frameSample[i + pos];
		dataI[i] = 0.0f;
	}
	int n = log10(frameSize)/log10(2);

	int x[n], xx;
	int i, j, k, b, p, L;
	float TR, TI, temp;


	/********** following code invert sequence ************/
	for(i=0;i<count[n];i++) 
	{ 
		xx=0; 
		for(j=0;j<n;j++) 
		x[j]=0; 
			for(j=0;j<n;j++) 
			{
			x[j]=(i/count[j])&0x01;
			} 
			for(j=0;j<n;j++) 
			{
			xx=xx+x[j]*count[n-j-1];
			} 
		dataI[xx]=dataR[i]; 
	}

 	for(i = 0; i < frameSize; i++)
	{
		dataR[i] = dataI[i]; 
		dataI[i] = 0;
		//cout << i << "  " <<dataI[i] << "   "<< i << "  "<< dataR[i] << endl;
	}
 
	/************** following code FFT *******************/
	for(L = 1; L <= n; L++)
	{ /* for(1) */
		b = 1; i = L - 1;
		while(i > 0) 
		{
			b = b * 2; i--;
		} /* b= 2^(L-1) */
		for(j = 0; j <= b-1; j++) /* for (2) */
		{
			p = 1; i = 9 - L;
			while(i > 0) /* p=pow(2,7-L)*j; */
			{
				p = p * 2; i--;
			}
			p = p * j;
			for(k = j; k < frameSize; k = k + 2*b) /* for (3) */
			{
				TR = dataR[k]; TI = dataI[k]; temp = dataR[k + b];
				dataR[k] = dataR[k] + dataR[k + b] * cos(2 * Pi * p / frameSize) + dataI[k + b] * sin(2 * Pi * p / frameSize);
				dataI[k] = dataI[k] - dataR[k + b] * sin(2 * Pi * p / frameSize) + dataI[k + b] * cos(2 * Pi * p / frameSize);
				dataR[k + b] = TR - dataR[k + b] * cos(2 * Pi * p / frameSize) - dataI[k + b] * sin(2 * Pi * p / frameSize);
				dataI[k + b] = TI + temp * sin(2 * Pi * p / frameSize) - dataI[k + b] * cos(2 * Pi * p / frameSize);
			} /* END for (3) */
		} /* END for (2) */
	} /* END for (1) */
	for(i = 0; i < frameSize / 2; i++)
	{ 
	    FFTSample[i + pos] = (dataR[i] * dataR[i] + dataI[i] * dataI[i]); 
	}	
	delete[] dataR;
	delete[] dataI;
}
 
//参数说明：frameSample为处理之后的数组，Sample为对样本进行预加重之后的数组
//          len为Sample的长度，frameSize为每帧的样本点个数，frameSampleLen为处理之后的长度 
double* mfccFrame(double *frameSample, double *Sample, int *len, int frameSize, int &frameSampleLen)
{
	double *hamWin;
	int hamWinSize = sampleRate * Win_Time;
	hamWin = new double[hamWinSize];
	Hamming(hamWin, hamWinSize);//计算hamming窗
	
	int hopStep = Hop_Time * sampleRate;
	int frameNum = ceil(double(*len) / hopStep);//计算一共会有多少帧
	frameSampleLen = frameNum * frameSize;//经过处理之后的长度 
	frameSample = new double[frameSampleLen];
	for(int i = 0; i < frameSampleLen; i++)
	    frameSample[i] = 0;
	
	double *FFTSample = new double[frameSampleLen];
	for(int i = 0; i < frameSampleLen; i++)
	    FFTSample[i] = 0;
	
	for(int i = 0; i * hopStep < *len; i++)//分帧 
	{
		for(int j = 0; j < frameSize; j++)
		{
			if(j < hamWinSize && i * hopStep + j < *len)
			    frameSample[i * frameSize + j] = Sample[i * hopStep + j] * hamWin[j];
		    else
		        frameSample[i * frameSize + j] = 0;//补0 
		}
		mfccFFT(frameSample, FFTSample, frameSize, i * frameSize);
	}
	
	ofstream fileFrame("./txt/Frame.txt");
	for(int j = 0; j < frameNum; j++)       //modify
	{
		for(int i = 0; i < frameSize; i++)
	    	fileFrame << frameSample[j * frameSize + i] << endl;  // del 100 *
	}
	delete []hamWin;
	return FFTSample; 
}
 
void DCT(double mel[400][40], double c[400][40], int frameNum)  
{  
    for(int k = 0; k < frameNum; k++)
    {
        for(int i = 0; i < 13; i++)  
        {  
        	for(int j = 0; j < filterNum; j++)  
        	{    
                c[k][i] += mel[k][j] * cos(Pi * i / (2 * filterNum) *  (2 * j + 1));  
                //if(k == 0 && i ==0)
                    //cout << c[0][0] << endl;
        	}	  
    	}
    }
    //cout << "c[0][0] = " << c[0][0] << endl;
}  
 
void computeMel(double mel[400][40], int sampleRate, double *FFTSample, int frameNum, int frameSize)
{
	double freMax = sampleRate / 2;//实际最大频率 
	double freMin = 0;//实际最小频率 
	double melFremax = 1125 * log(1 + freMax / 700);//将实际频率转换成梅尔频率 
	double melFremin = 1125 * log(1 + freMin / 700);
	double melFilters[filterNum][3];
	double k = (melFremax - melFremin) / (filterNum + 1);
	
	double *m = new double[filterNum + 2];
	double *h = new double[filterNum + 2];
	double *f = new double[filterNum + 2];
	
	for(int i = 0; i < filterNum + 2; i++)
	{
		m[i] = melFremin + k * i;
		h[i] = 700 * (exp(m[i] / 1125) - 1);//将梅尔频率转换成实际频率 
		f[i] = floor((frameSize + 1) * h[i] / sampleRate);
	}		
 
    delete[] m;  
    delete[] h;  
    //delete[] f;
	
	for(int k = 0; k < frameNum; k++)
	{
		for(int j = 0; j < filterNum; j++)
		    mel[k][j] = 0;
	}
	
	//计算出每个三角滤波器的输出: 对每一帧进行处理 	
	for(int i = 0; i < frameNum; i++)
	{
		for(int j = 1; j <= filterNum; j++)
		{
			double temp = 0;
			for(int z = 0; z < frameSize; z++)
			{
				if(z < f[j - 1])
				    temp = 0;
				else if(z >= f[j - 1] && z <= f[j])
				    temp = (z - f[j - 1]) / (f[j] - f[j - 1]);
				else if(z >= f[j] && z <= f[j + 1])
				    temp = (f[j + 1] - z) / (f[j + 1] - f[j]);
				else if(z > f[j + 1])
				    temp = 0;
				mel[i][j - 1] += FFTSample[i * frameSize + z] * temp; 
			}
		}
    }
	
	ofstream fileMel("./txt/Mel.txt");
	for(int j = 0; j < frameNum; j++)       //modify
	{
		for(int i = 0; i < filterNum; i++)
			fileMel << mel[j][i] << endl;
	}
	//取对数 
	for(int i = 0; i < frameNum; i++)
	{
		for(int j = 0; j < filterNum; j++)
		{
			if(mel[i][j] <= 0.00000000001 || mel[i][j] >= 0.00000000001)
				mel[i][j] = log(mel[i][j]);
		}
	} 
}
 
void writeToFile(int frameNum, int frameSize, int hopStep, double DCT[400][40], double *sample, double *Sample, double *frameSample, double *FFTSample, double mel[400][40])
{
	ofstream fileDCT("./txt/DCT.txt");

	ofstream fileFft("./txt/Fft.txt");
	ofstream fileLogMel("./txt/LogMel.txt");

	for(int j = 0; j < frameNum; j++)       //modify
	{
		for(int i = 0; i < frameSize; i++)
		{
			fileFft << FFTSample[j * frameSize + i] << endl; 
		}
	
		for(int i = 0; i < filterNum; i++)
			fileLogMel << mel[j][i] << endl;
	}
	//归一化
	/*for(int i = 0; i < 13; i++) 
	{
		double sum = 0.0f;
		double Mrecording = 0.0f;
		for(int j = 0; j < frameNum; j++)
		{
			sum = sum + DCT[j][i];
		}
		Mrecording = sum / frameNum;
		cout << Mrecording << endl;
		for(int j = 0; j < frameNum; j++)
		{
			DCT[j][i] = abs(DCT[j][i] - Mrecording);
		}
	} */
	
	for (int i = 0; i < 13; i++)//write DCT
	{
		for (int j = 0; j < frameNum; j++)
		    fileDCT << DCT[j][i] << " ";
		fileDCT << endl;
	}
}
 
//MFCC
void MFCC(double *sample, int len)
{
	double factor = 0.99;//预加重参数
	double *Sample;

	//1.预加重 
	Sample = pre_emphasizing(sample, len, factor);
	
	//计算出每帧有多少个点，然后算出最接近点的个数的2^k，使得每帧的点的个数为2^k，以便进行补0 
	int frameSize = (int)pow(2, ceil( log(Win_Time * sampleRate) / log(2.0) ));
	double *frameSample = NULL, *FFTSample = NULL;
	int frameSampleLen;
	
	//分帧、加窗、功率谱 
	FFTSample = mfccFrame(frameSample, Sample, &len, frameSize, frameSampleLen);
	
	int hopStep = Hop_Time * sampleRate;//隔hopStep个样本点分一次帧 
	int frameNum = ceil(double(len) / hopStep);//计算所有样本点一共有多少帧 
 
    double mel[400][40] = {0};
    computeMel(mel, sampleRate, FFTSample, frameNum, frameSize);
    
    double c[400][40] = {0};

    DCT(mel, c, frameNum);
    
    writeToFile(frameNum, frameSize, hopStep, c, sample, Sample, frameSample, FFTSample, mel);
}
 
int main()
{
	ifstream filedata("./data.txt");
	int len = 0;
	//读取wav文件的数值 
	double *sample = new double[200000];
	while(!filedata.eof())
	{
	    filedata >> sample[len];
	    sample[len] = sample[len] * 1000;
		len++;
	}
	cout << len-1 << endl;
	MFCC(sample, len - 1);
	delete [] sample;
	
	
	return 0;
} 

