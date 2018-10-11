#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <math.h>
#include "encode_wav/encode_wav.h"

const double pi = 3.1415926535;  //Π
const int  MAX_SILENCE = 8; //最大静音帧数
const int MIN_LEN = 5; //声音最小长度帧

/*端点检测数据*/
typedef struct
{
	int silence;//静音帧数
	int count;//非静音帧数

	int start;//起始帧
	int end;//终止帧
}VAD_DATA;

/* 阈值 */
const double T1 = 0.0010, T2 = 0.045;
/* 函数声明 */
double sample_normalized(WAVE wav); //归一化
void hanning(int n, double *w);
int vad_param(WAVE wav, double *w, int wlen, int inc);

WAVE res_wav;

/* 主函数 */
int main(int argv, char *argc[])
{
	double *win = NULL; //窗函数
	WAVE wav;

	int wlen=200,inc=80;         //设置帧长和帧移 
	
    //	printf("8字节对齐：\nwav:%ld\n",sizeof(wav));
	if(argv == 2)
	{
		unencode_wav(argc[1], &wav);
	}
	else
		unencode_wav("./wav_file/vad_faqing1.wav", &wav);

    sample_normalized(wav); //normalized sample

 	win = (double *)malloc(wlen * sizeof(double));
	if( win == NULL)
		return -1;
    hanning(wlen, win);  //得到汉宁窗

	vad_param(wav, win, wlen, inc);//端点检测
  //  encode_wav("./out.wav",wav);
	free_wav(wav);
	if(win != NULL)
		free(win);

	return 0;
}

/*
*功能：音频数据归一化
*参数：wav数据
*返回值：无
*/
double sample_normalized(WAVE wav)
{
	int i, len;
	double max = 0.0;
	len = (wav.ChunkSize -44 + 8) / (2 * wav.NumChannels);
	for(i = 0; i < len; i++)   //获取极值
	{
		if(max < fabs(*(wav.sample + i)))
			max = fabs(*(wav.sample + i));
	}
    
	for(i = 0; i < len; i++)   //归一化
	{
		*(wav.sample + i) /= max;
		//printf("%lf\n", *(wav.sample + i));
	}
 //	printf("count :%d \n",len);
	return max;
}

/*
*功能：计算汉宁窗系数
*参数：w窗  wlen帧长
*返回值：无
*/
void hanning(int wlen, double *w)    // 帧长  窗
{
	int i = 0;
	int m;
	if(0 == wlen%2)
	    m = wlen/2;
	else
		m = (wlen+1)/2;

	for(i = 0; i < m; i++)
	{
		*(w + i) = 0.5 * (1 - cos(2 * pi * (i + 1) / (wlen + 1)) );
		//printf("%lf\n",*(w + i));
	}
	for(; i <wlen; i++)
	{
		*(w + i) = *(w + wlen - i - 1);
		//printf("%lf\n",*(w + i));
	}
}
/*
*功能：端点检测
*参数：wav数据  w窗  wlen帧长  inc帧移  f帧数据指针
*返回值：无
*/
int vad_param(WAVE wav, double *w, int wlen, int inc) 
{
	char out_file_name[26];//输出文件名
	int i =0, j = 0;
	int wav_len,  win_len; //样本个数、 帧（窗）长

	int nf; //帧数
	int *indf = NULL; //每帧在数据中的位置
	double *etemp = NULL;//每帧能量
	double max_etemp = 0;
	double *f = NULL; //帧数据  （帧数 * 帧长）

    VAD_DATA vad_data[2000];//预定义 2000 段语音
	int status = 0; //状态标志
	int xn = 0;//记录语言段
	
/*分帧*/
{
	wav_len = (wav.ChunkSize -44 + 8) / (2 * wav.NumChannels); //数据个数
	win_len = wlen; //窗长 = 帧长
	nf = (int)((wav_len - win_len + inc) / inc);
	//printf("wav_len :%d  , nf = %d\n",wav_len,nf);
	/* 申请空间 */
	f = (double *)malloc(nf * wlen * sizeof(double));
	if(f == NULL)
	{
		printf("f:malloc error!\n");
		return -1;
	}
	indf = (int *) malloc(nf * sizeof(int));
	if(indf == NULL)
	{
		printf("f:malloc error!\n");
		return -1;
	}

	/* 数据清零 */
	memset(f, 0, nf * wlen * sizeof(double));

	/* 获取每帧在样本中的位置 */
	for(i = 0; i < nf; i++)
	{
		*(indf + i) = i * inc;
		//printf("indf:%d \n",*(indf + i));
	}

	/* 加窗函数 */
	for(i = 0; i < nf; i++) 
	{
		for(j = 0; j < win_len; j++)
		{
			*(f + *(indf +i ) + j) = *(wav.sample + *(indf +i ) + j) * (*(w + j));
		//	printf("%lf \n", *(f + *(indf +i ) + j));
		}
	}
}
/* 计算短时能量 */
{
	etemp = (double * )malloc(nf * sizeof(double));
	memset(etemp, 0, nf * sizeof(double));

	for(i = 0; i < nf; i++) //求短时平均能量
	{
		for(j = 0; j < win_len; j++)
		{
			*(etemp + i) += *(f + *(indf +i ) + j) * (*(f + *(indf +i ) + j) );
		}
		if(max_etemp < *(etemp + i))
		 	max_etemp = *(etemp + i);
		//printf("%lf\n",*(etemp + i));
	}
	for(i = 0; i < nf; i++)  //能量归一化
	{
		*(etemp + i) /= max_etemp;
		//printf("%lf\n",*(etemp + i));
	}
}
/*端点检测*/
{
	status = 0; // 0 静音  1可能开始  2确定开始
	vad_data[0].count = 0;
	vad_data[0].silence = 0;
	vad_data[0].start = 0;
	vad_data[0].end = 0;

	for(i = 0; i < nf; i++) //开始端点检测
	{
		switch( status )
		{
			case 0: 
			case 1:
				if( *(etemp + i) > T2 )  //确定进入语言段
				{
					status = 2;
					vad_data[xn].start = i - vad_data[xn].count;
					vad_data[xn].silence = 0;
					vad_data[xn].count += 1;
				}
				else if( *(etemp + i) > T1 )
				{
					status = 1;
					vad_data[xn].count += 1;
				}
				else 
				{
					status = 0;
					vad_data[xn].count = 0;
					vad_data[xn].start = 0;
					vad_data[xn].end = 0;
				}
				break;
			case 2:
				if( *(etemp + i) > T1) //保持在语言段
				{
				//	printf("i: %d,count: %d ++\n",i, vad_data[xn].count);
					vad_data[xn].count += 1;
					if(i == nf - 1)
					{
						status = 3;
						vad_data[xn].end = vad_data[xn].start + vad_data[xn].count;
						res_wav = wav; // 使用原始数据初始化

						res_wav.data_sample = NULL;// 不手动分配内存  使用源音频 wav 的内存空

						res_wav.data_sample = (DATA_SAMPLE *)(wav.data_sample + *(indf + vad_data[xn].start));  //起始点
						//res_wav.Subchunk2_Size = vad_data[xn].count * 2 * inc;
						res_wav.ChunkSize = vad_data[xn].count * 2 * inc + WAVEHEAD_SIZE - 8;
						sprintf((char*)out_file_name, "./vad_out_%d.wav",xn);
						printf("1count: %d, start: %d\n",vad_data[xn].count, vad_data[xn].start);
						printf("1Out file name: %s , Size: %d \n", out_file_name, res_wav.ChunkSize + 8);
						encode_wav(out_file_name, res_wav);
					}
				}
				else //语音将结束
				{	
					if((vad_data[xn].silence < MAX_SILENCE) && (i != nf-1))// 静音长度还不够长
					{
					//	printf("silence %d ++\n",vad_data[xn].count);
						vad_data[xn].count += 1;
						vad_data[xn].silence += 1;
					}
					else if(vad_data[xn].count < MIN_LEN) // 语言太短 视为噪音
					{
						printf("噪音\n");
						status = 0;
						vad_data[xn].count = 0;
						vad_data[xn].silence = 0;
						vad_data[xn].start = 0;
					    vad_data[xn].end = 0;
					}
					else//语音结束   保存当前语音
					{
						status = 3;
						vad_data[xn].end = vad_data[xn].start + vad_data[xn].count;
						res_wav = wav; // 使用原始数据初始化
				    					// 不分配内存  res_wav使用源音频 wav 的内存空间
						res_wav.data_sample = NULL;

						res_wav.data_sample = (wav.data_sample + *(indf + vad_data[xn].start));  //起始点
						//res_wav.Subchunk2_Size = vad_data[xn].count * 2 * inc;
						res_wav.ChunkSize = vad_data[xn].count * 2 * inc + WAVEHEAD_SIZE - 8;
						//printf("subsize: %d size: %d\n",res_wav.Subchunk2_Size, res_wav.ChunkSize);
						sprintf((char*)out_file_name, "./vad_out_%d.wav",xn);
						printf("2count: %d, start: %d\n",vad_data[xn].count, vad_data[xn].start);
						printf("2Out file name: %s , Size: %d \n", out_file_name, res_wav.ChunkSize + 8);
						encode_wav(out_file_name, res_wav);
					}
				}
				break;
			case 3: //为下一段语音做准备
						status = 0;
						xn += 1;
						vad_data[xn].count = 0;
						vad_data[xn].silence = 0;
						vad_data[xn].start = 0;
					    vad_data[xn].end = 0;
				break;
			default :
				break;
		}		
	}
}

	free(f);
	free(indf);
	free(etemp);
	f = NULL;
	indf = NULL;
	etemp = NULL;
	return 0;
}


