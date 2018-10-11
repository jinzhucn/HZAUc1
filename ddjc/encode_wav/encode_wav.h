#ifndef _ENCODE_WAV_H_
#define _ENCODE_WAV_H_
#include <stdio.h>

#define DATA_BEGIN (36)
#define FMT_BEGIN (12)
#define WAVEHEAD_SIZE (44)
typedef struct
{
	char l; // 样本低字节
	char h; // 样本高字节
}DATA_SAMPLE; //2 Byte

typedef struct
{
	char ll;
	char lh; 
	char rl;
	char rh;
}DATA_SAMPLE_2; //4 Byte

typedef struct
{
	char ChunkID[4]; // 'R''I''F''F'  
	int ChunkSize;   //8
	char Format[4];  //'W''A''V''E' 

	char Subchunk1_ID[4];//'f''m''t'' ' 8
	int Subchunk1_Size; //
	short AudioFormat;  //
	short NumChannels;  //8
	int SampleRate;     //
	int ByteRate;       //8
	short BlockAlign;   //
	short BitsPerSample;//

	char Subchunk2_ID[4];//'d''a''t''a' 8
	int Subchunk2_Size; //空余4字节保证对齐
						//8
	double *sample;  //8	 16bit 样本数据
	DATA_SAMPLE *data_sample;//样本数据  高低字节
	DATA_SAMPLE_2 *data_sample_2;//样本数据  高低字节

}WAVE; 

int encode_wav(char *aim_file_name, WAVE wav);
int unencode_wav(const char *exist_wav_file, WAVE *wav);
int free_wav(WAVE wav);

#endif