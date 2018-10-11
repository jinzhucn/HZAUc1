#include <stdio.h>
#include <stdlib.h>
#include "encode_wav.h"

/*wav 编码*/
int encode_wav(char *aim_file_name, WAVE wav)
{ 
	FILE *fp;
	fp = fopen(aim_file_name, "w");
	if(NULL == fp)
	{
	    printf("open %s error !\n",aim_file_name);
	    return (-1);
	}

	/* 声明为单声道 */
{
    wav.Subchunk1_Size = 16; //
	wav.AudioFormat = 1;  //
	wav.NumChannels = 1;//单声道
	wav.BitsPerSample = 16;
	wav.BlockAlign = wav.NumChannels * wav.BitsPerSample / 8;
	wav.ByteRate = wav.BlockAlign * wav.SampleRate;
	wav.Subchunk2_Size = wav.ChunkSize - WAVEHEAD_SIZE + 8 ; //空余4字节保证对齐
}

	/* 写文件 */
	fseek(fp, 0, SEEK_SET);
	fwrite(&wav, WAVEHEAD_SIZE, 1, fp);
	fwrite(wav.data_sample, 2, (wav.ChunkSize - WAVEHEAD_SIZE + 8 ) / 2, fp);

	fclose(fp);
	return (0);
}

/*wav 解码*/
int unencode_wav(const char *exist_wav_file, WAVE *wav)
{
	int i = 0;
	FILE *fp;
	fp = fopen(exist_wav_file, "r");
	if(NULL == fp)
	{
	    printf("open %s error !\n", exist_wav_file);
	    return (-1);
	}

	/*读取文件头 （pcm无损格式，无fact块，文件头 WAVEHEAD_SIZE = 44 字节）*/
	fseek(fp, 0, SEEK_SET);
	fread(wav, WAVEHEAD_SIZE, 1, fp);

	/*打印文件头信息*/
printf("\
ChunkSize: %d\n\
Subchunk1_Size: %d\n\
AudioFormat: %d\n\
NumChannels: %d\n\
SampleRate: %d\n\
ByteRate: %d\n\
BlockAlign: %d\n\
BitsPerSample: %d\n\
Subchunk2_Size: %d\n",
	 wav->ChunkSize,
	 wav->Subchunk1_Size,
	 wav->AudioFormat,
	 wav->NumChannels,
	 wav->SampleRate,
	 wav->ByteRate,
	 wav->BlockAlign,
	 wav->BitsPerSample,
	 wav->Subchunk2_Size);


	/*申请堆空间 存储音频数据，双声道转换为单声道数据存储*/
	wav->data_sample = (DATA_SAMPLE *)malloc( (wav->ChunkSize - WAVEHEAD_SIZE + 8 )/wav->NumChannels);//单声道原始数据（字节数据）
	wav->sample = (double *)malloc((wav->ChunkSize - WAVEHEAD_SIZE + 8 )/(2 * wav->NumChannels) * sizeof(double));//合成后数据（真实数据）
	if(wav->data_sample == NULL || wav->sample == NULL)
		{
			printf("malloc error \n");
			return (-1);
		}
	if(2 == wav->NumChannels)//如果是双声道 则申请双声道数据空间
	{
		wav->data_sample_2 = (DATA_SAMPLE_2 *)malloc((wav->ChunkSize - WAVEHEAD_SIZE + 8 ) );//双声道原始数据
		if(wav->data_sample_2 == NULL)
		{
			printf("malloc error! \n");
			return -1;
		}
	}
	else
	{
		wav->data_sample_2 = NULL;
	}

	/*读取音频数据（若是双声道  则提取单声道）*/
	fseek(fp, WAVEHEAD_SIZE, SEEK_SET);
	if( 1 == wav->NumChannels)//单声道
	{
		while(!feof(fp))
		{
			fread(wav->data_sample + i, 2, 1, fp);
			*(wav->sample + i) = ((wav->data_sample + i)->l) | ((wav->data_sample + i)->h) << 8;
			//printf("%lf\n",*(wav->sample + i));
			i++;
		}
	}
	else if(2 == wav->NumChannels)//双声道
	{
		while(!feof(fp))
		{
			fread(wav->data_sample_2 + i, 4, 1, fp);
			(wav->data_sample + i)->l = (wav->data_sample_2 + i)->ll;//提取左声道数据
			(wav->data_sample + i)->h = (wav->data_sample_2 + i)->lh;
			*(wav->sample + i) = ((wav->data_sample + i)->l) | ((wav->data_sample + i)->h) << 8;
			//printf("%lf\n",*(wav->sample + i));
			i++;
		}
		free(wav->data_sample_2);//释放双声道空间
		wav->data_sample_2 = NULL;
	}

	fclose(fp);
	return (0);
}

/*释放 wav 内存*/
int free_wav(WAVE wav)
{
	if(wav.data_sample == NULL)
		return -1;
	free(wav.data_sample);
	wav.data_sample = NULL;

	if(wav.sample == NULL)
	return -1;
	free(wav.sample);
	wav.sample = NULL;

	return 0;
}