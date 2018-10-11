#ifndef  MFCC_H
#define  MFCC_H
#include <stdio.h>
 
typedef struct AudioInfo
{
    int sample_rate;
    int frame_len;
    int frame_shift;
    float* data_in;
    int data_len;  
    AudioInfo():sample_rate(0),frame_len(0),frame_shift(0),data_in(NULL),data_len(0){};//构造函数注释掉，变成纯C版本
};
typedef struct MelBankInfo
{
    float **filter;
    int     nfilters;
    int     nfft;
    int     low;
    int     high;
    MelBankInfo():filter(NULL),nfilters(0),nfft(0),low(0),high(0){};
};
 
typedef struct DctInfo
{
    float  **coeff;
    int     dctlen;
    DctInfo():coeff(NULL),dctlen(0){};
};
 
typedef enum MFCC_TYPE
{
    MFCC_STD,
    MFCC_DIFF_1,
    MFCC_DIFF_2,
};
 
typedef struct MfccInfo
{    
    MelBankInfo   melbank;
    DctInfo       dct;
 
    int     nframes;
    int     out_nframes;//可输出的特征数
    float *frame_data;
    float *data_in;
    float  *window;
    float  *lift_window;   
    int     frame_len;
    int     frame_shift;
    float  *pre1;
    float  *pre2;
    float  *cur;
    float  *next1;
    float  *next2;
 
    float  *diff_pre1;
    float  *diff_pre2;
    float  *diff_cur;
    float  *diff_next1;
    float  *diff_next2;
    MFCC_TYPE m_type;
    MfccInfo():nframes(0),out_nframes(0),window(NULL),lift_window(NULL),cur(NULL),frame_len(0),pre1(NULL),pre2(NULL),next1(NULL),next2(NULL),m_type(MFCC_STD),frame_data(NULL),frame_shift(0){};
};
 
MfccInfo*  MfccInit(AudioInfo audioctx, int nfft, int low, int high, int nfilters, int ndcts, MFCC_TYPE type);
int        Mfcc_Frame_std(MfccInfo *p, int  iframe, float *out, int len);
int        Mfcc_Frame_diff1(MfccInfo *p, int iframe, float *out, int len);
int        Mfcc_Frame_diff2(MfccInfo *p, int iframe, float *out, int len);
void       MfccDestroy(MfccInfo *data);
 
#endif
