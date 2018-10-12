#include<stdio.h>
/*
	用户头文件
*/
#include "fftw3.h"   //fftw库头文件
#include "math.h"
#include "stdlib.h"
//#pragma comment(lib, "libfftw3-3.lib")  
//#pragma comment(lib, "libfftw3f-3.lib")  
//#pragma comment(lib, "libfftw3l-3.lib") 
/*
	宏定义,处理必须的宏
*/
#define LENGTH		184432	//实例音频文件的长度，根据声音样本的长度决定
#define FRMSZ		20		//帧长20ms
#define OVLP		50		//帧重叠50%
#define NOISEFR		6		//噪声帧6帧
#define FLOOR		0.002	
#define FS			24000	//采样率，根据声音样本决定
#define NBAND		5		//子带数，4-8之间
#define PI			3.14159265358979323846
/*
	函数声明
*/
void mband(double *in, double *out);
void berouti(double * SNR, double * a, int nbands, int nframes);
void noiseupt(double * x_magsm, double * n_spect, int nframes, int fftlen);
void hamming(int N, double *win);
void frame(double *sdata, int ndata, double *windows, int nwind, int frmshift,	int offset, int trunc, double *fdata);
void scale_mtx(double *fdata, double *tdata, int ros, int col, double *windows);
/*
	主函数
*/
int main(void)
{
	int frmelen = (int)floor(FRMSZ*FS / 1000);
	int ovlplen = (int)floor(frmelen*OVLP / 100);
	int nframes = (int)(ceil((double)LENGTH / ovlplen));
	int outlen = (frmelen / 2) * (nframes + 1);
	double *pure_speech = (double *)malloc(outlen * sizeof(double));
	if (pure_speech == NULL)
	{
		printf("申请失败\n");
		return 0;
	}
	double *noise_speech = (double *)malloc(LENGTH * sizeof(double));
	if (noise_speech == NULL)
	{
		printf("申请失败\n");
		return 0;
	}
	FILE *fp = fopen("C://Users/35336/Desktop/data.txt", "rb");
	if (fp == NULL)
	{
		printf("打开失败\n");
		return 0;
	}
	for (int i = 0; i < LENGTH; i++)
	{
		int nRes = fscanf(fp, "%lf", &noise_speech[i]);
		if (nRes == -1)
		{
			fclose(fp);
			printf("-----读取异常-----\n");
		}
	}
	fclose(fp);	
	// 调用处理
	mband(noise_speech, pure_speech);


	fp = fopen("C://Users/35336/Desktop/data_new.txt", "w");
	for (int i = 0; i < outlen; i++)
	{
		fprintf(fp, "%lf ", pure_speech[i]);
	}
	fclose(fp);	

	free(pure_speech);
	free(noise_speech);

	return 0;
}

/*
	函数体
*/
void mband(double *in, double *out) 
{
	int frmelen = (int)floor(FRMSZ*FS / 1000);		//帧长
	int ovlplen = (int)floor(frmelen*OVLP / 100);
	int	cmmnlen = frmelen - ovlplen;
	int fftlen = 2;
	while (fftlen < frmelen)
		fftlen *= 2;
	int lobin[NBAND], hibin[NBAND], bandsz[NBAND];
	bandsz[1] = (int)floor(fftlen / (2 * NBAND));
	for (int i = 0; i < NBAND; i++)
	{
		lobin[i] = (i)*bandsz[1] + 1;
		hibin[i] = lobin[i] + bandsz[1] - 1;
		bandsz[i] = bandsz[1];
	}
	//创建hamming窗
	double *win = (double *)malloc(frmelen * sizeof(win));
	if (win == NULL)
	{
		printf("申请失败\n");
		return 0;
	}
	hamming(frmelen, win);
	//噪声估计
	int nframes = (int)(ceil((double)LENGTH / ovlplen));
	double *n_spect = (double *)malloc(fftlen * nframes * sizeof(double));
	if (n_spect == NULL)
	{
		printf("申请失败\n");
		return 0;
	}
	for (int i = 0; i < fftlen * nframes; i++)
		n_spect[i] = 0;
	double *temp = (double *)malloc(fftlen * sizeof(double));
	if (temp == NULL)
	{
		printf("申请失败\n");
		return 0;
	}
	double *n_magsq = (double *)malloc(fftlen * sizeof(double));
	if (n_magsq == NULL)
	{
		printf("申请失败\n");
		return 0;
	}
	double *n_ph = (double *)malloc(fftlen * sizeof(double));
	if (n_ph == NULL)
	{
		printf("申请失败\n");
		return 0;
	}
	fftw_complex *n_fft = (fftw_complex *)fftw_malloc(fftlen * sizeof(fftw_complex));
	fftw_plan p;
	int i = 0;
	for (int j = 0; j < fftlen; j++)
		n_magsq[j] = 0;
	for (int k = 0; k < NOISEFR; k++)
	{
		for (int j = 0; j < fftlen; j++)
		{
			if (j < frmelen)
				temp[j] = in[j + i] * win[j];
			else
				temp[j] = 0;
		}
		p = fftw_plan_dft_r2c_1d(fftlen, temp, n_fft, FFTW_ESTIMATE);
		fftw_execute(p);
		i += frmelen;
		for (int j = 0; j < fftlen; j++)
			n_magsq[j] += pow(sqrt(pow(n_fft[j][0], 2) + pow(n_fft[j][1], 2)), 2);
	}
	for (int j = 0; j < fftlen; j++)
	{
		if (j <= fftlen / 2)
			n_spect[j*nframes + 0] = sqrt(n_magsq[j] / NOISEFR);
		else
			n_spect[j*nframes + 0] = n_spect[(fftlen - j)*nframes + 0];
	}
	//	for (int j = 0; j < fftlen; j++)
	//		for (int i = 1; i < nframes; i++)
	//			n_spect[j*nframes + i] = n_spect[j*nframes + i - 1];

	fftw_destroy_plan(p);
	fftw_free(n_fft);
	free(temp);
	free(n_magsq);
	free(n_ph);
	
	//分帧
	double *frame_data = (double *)malloc(nframes * fftlen * sizeof(double));
	if (frame_data == NULL)
	{
		printf("申请失败\n");
		return 0;
	}
	for (int i = 0; i < nframes * fftlen; i++)
		frame_data[i] = 0;
	frame(in, LENGTH, win, frmelen, ovlplen, 0, 0, frame_data);

	//数据处理
	fftw_complex *x_fft = (fftw_complex *)fftw_malloc(fftlen * nframes * sizeof(fftw_complex));
	fftw_complex *x_fft_temp = (fftw_complex *)fftw_malloc(fftlen * 1 * sizeof(fftw_complex));
	double *frame_data_temp = (double *)malloc(fftlen * 1 * sizeof(double));
	if (frame_data_temp == NULL)
	{
		printf("申请失败\n");
		return 0;
	}
	for (int j = 0; j < nframes; j++)
	{
		for (int i = 0; i < fftlen; i++)
			frame_data_temp[i] = frame_data[i*nframes + j];
		p = fftw_plan_dft_r2c_1d(fftlen, frame_data_temp, x_fft_temp, FFTW_ESTIMATE);
		fftw_execute(p);
		for (int i = 0; i < fftlen; i++)
		{
			if (i < fftlen / 2)
			{
				x_fft[i*nframes + j][0] = x_fft_temp[i][0];
				x_fft[i*nframes + j][1] = x_fft_temp[i][1];
			}
			else
			{
				x_fft[i*nframes + j][0] = x_fft_temp[fftlen - i][0];
				x_fft[i*nframes + j][1] = x_fft_temp[fftlen - i][1];
			}
		}
	}
	free(frame_data_temp);
	fftw_destroy_plan(p);
	//幅值谱和相位谱
	double *x_mag = (double *)malloc(fftlen * nframes * sizeof(double));
	if (x_mag == NULL)
	{
		printf("申请失败\n");
		return 0;
	}
	double *x_ph = (double *)malloc(fftlen * nframes * sizeof(double));
	if (x_ph == NULL)
	{
		printf("申请失败\n");
		return 0;
	}
	for (int i = 0; i < fftlen; i++)
	{
		for (int j = 0; j < nframes; j++)
		{
			x_mag[i*nframes + j] = sqrt(pow(x_fft[i*nframes + j][0], 2) + pow(x_fft[i*nframes + j][1], 2));
			x_ph[i*nframes + j] = atan2(x_fft[i*nframes + j][1], x_fft[i*nframes + j][0]);
			if (i > fftlen / 2)
				x_ph[i*nframes + j] = -x_ph[i*nframes + j];
		}
	}
	//噪声更新
	noiseupt(x_mag, n_spect, nframes, fftlen);

	//统计每个子带的SNR
	double *SNR_x = (double *)malloc(nframes * NBAND * sizeof(double));
	if (SNR_x == NULL)
	{
		printf("申请失败\n");
		return 0;
	}
	double *beta_x = (double *)malloc(nframes * NBAND * sizeof(double));
	if (beta_x == NULL)
	{
		printf("申请失败\n");
		return 0;
	}
	int start, stop;
	double tnorm1, tnorm2;

	for (int i = 0; i < NBAND - 1; i++)
	{
		start = lobin[i];
		stop = hibin[i];
		for (int j = 0; j < nframes; j++)
		{
			tnorm1 = 0; tnorm2 = 0;
			for (int m = start - 1; m < stop; m++)
			{
				tnorm1 += pow(x_mag[m*nframes + j], 2);
				tnorm2 += pow(n_spect[m*nframes + j], 2);
			}
			SNR_x[i*nframes + j] = 10 * log10(tnorm1 / tnorm2);
		}
	}
	//最后一带
	start = lobin[NBAND - 1];
	stop = hibin[NBAND - 1];
	for (int j = 0; j < nframes; j++)
	{
		tnorm1 = 0; tnorm2 = 0;
		for (int m = start - 1; m < fftlen / 2 + 1; m++)
		{
			tnorm1 += pow(x_mag[m*nframes + j], 2);
			tnorm2 += pow(n_spect[m*nframes + j], 2);
		}
		SNR_x[(NBAND - 1)*nframes + j] = 10 * log10(tnorm1 / tnorm2);
	}
	berouti(SNR_x, beta_x, NBAND, nframes);

	//开始谱减过程
	double *sub_speech_x = (double *)malloc(nframes * fftlen * sizeof(double));
	if (sub_speech_x == NULL)
	{
		printf("申请失败\n");
		return 0;
	}
	for (int i = 0; i < nframes * fftlen; i++)
		sub_speech_x[i] = 0;
	for (int i = 0; i < NBAND - 1; i++)
	{
		double *sub_speech = (double *)malloc(nframes * bandsz[i] * sizeof(double));
		if (sub_speech == NULL)
		{
			printf("申请失败\n");
			return 0;
		}
		double *n_spec_sq = (double *)malloc(1 * bandsz[i] * sizeof(double));
		if (n_spec_sq == NULL)
		{
			printf("申请失败\n");
			return 0;
		}
		start = lobin[i];
		stop = hibin[i];
		switch (i)
		{
		case 0:
			for (int j = 0; j < nframes; j++)
			{
				for (int k = 0; k < bandsz[i]; k++)
				{
					n_spec_sq[k] = pow(n_spect[(start - 1 + k)*nframes + j], 2);
					sub_speech[k*nframes + j] = pow(x_mag[(start - 1 + k)*nframes + j], 2) - beta_x[i*nframes + j] * n_spec_sq[k];
				}
			}
			break;
		default:
			for (int j = 0; j < nframes; j++)
			{
				for (int k = 0; k < bandsz[i]; k++)
				{
					n_spec_sq[k] = pow(n_spect[(start - 1 + k)*nframes + j], 2);
					sub_speech[k*nframes + j] = pow(x_mag[(start - 1 + k)*nframes + j], 2) - beta_x[i*nframes + j] * n_spec_sq[k] * 2.5;
				}
			}
			break;
		}
		for (int j = 0; j < nframes; j++)
		{
			for (int k = 0; k < bandsz[i]; k++)
			{
				if (sub_speech[k*nframes + j] < 0)
					sub_speech[k*nframes + j] = FLOOR * pow(x_mag[(start - 1 + k)*nframes + j], 2);
				sub_speech[k*nframes + j] += 0.05*pow(x_mag[(start - 1 + k)*nframes + j], 2);
				sub_speech_x[(start - 1 + k)*nframes + j] += sub_speech[k*nframes + j];
			}
		}
		free(n_spec_sq);
		free(sub_speech);
	}
	//最后一带
	start = lobin[NBAND - 1];
	stop = fftlen / 2 + 1;
	double *sub_speech = (double *)malloc(nframes * (fftlen / 2 + 1) * sizeof(double));
	if (sub_speech == NULL)
	{
		printf("申请失败\n");
		return 0;
	}
	double *n_spec_sq = (double *)malloc(1 * (fftlen / 2 + 1) * sizeof(double));
	if (n_spec_sq == NULL)
	{
		printf("申请失败\n");
		return 0;
	}
	for (int j = 0; j < nframes; j++)
	{
		for (int k = 0; k < (stop - start + 1); k++)
		{
			n_spec_sq[k] = pow(n_spect[(start - 1 + k)*nframes + j], 2);
			sub_speech[k*nframes + j] = pow(x_mag[(start - 1 + k)*nframes + j], 2) - beta_x[(NBAND - 1)*nframes + j] * n_spec_sq[k] * 1.5;
		}
	}
	for (int j = 0; j < nframes; j++)
	{
		for (int k = 0; k < (stop - start + 1); k++)
		{
			if (sub_speech[k*nframes + j] < 0)
				sub_speech[k*nframes + j] = FLOOR * pow(x_mag[(start - 1 + k)*nframes + j], 2);
			sub_speech[k*nframes + j] += 0.01*pow(x_mag[(start - 1 + k)*nframes + j], 2);
			sub_speech_x[(start - 1 + k)*nframes + j] += sub_speech[k*nframes + j];
		}
	}
	free(n_spec_sq);
	free(sub_speech);
	for (int i = (fftlen / 2 + 1); i < fftlen; i++)
	{
		for (int j = 0; j < nframes; j++)
		{
			sub_speech_x[i*nframes + j] = sub_speech_x[(fftlen - i)*nframes + j];
		}
	}
	//重建
	fftw_complex *y_fft_temp = (fftw_complex *)fftw_malloc(fftlen * 1 * sizeof(fftw_complex));
	double *y_ifft_temp = (double *)malloc(fftlen * 1 * sizeof(double));
	if (y_ifft_temp == NULL)
	{
		printf("申请失败\n");
		return 0;
	}
	double *y_ifft = (double *)malloc(fftlen * nframes * sizeof(double));
	if (y_ifft == NULL)
	{
		printf("申请失败\n");
		return 0;
	}
	for (int j = 0; j < nframes; j++)
	{
		for (int i = 0; i < fftlen; i++)
		{
			y_fft_temp[i][0] = sqrt(sub_speech_x[i*nframes + j])*cos(x_ph[i*nframes + j]);
			y_fft_temp[i][1] = sqrt(sub_speech_x[i*nframes + j])*sin(x_ph[i*nframes + j]);
		}
		p = fftw_plan_dft_c2r_1d(fftlen, y_fft_temp, y_ifft_temp, FFTW_ESTIMATE);
		fftw_execute(p);
		for (int i = 0; i < fftlen; i++)
		{
			y_ifft[i*nframes + j] = y_ifft_temp[i] / fftlen;
		}
	}
	int outlen = (frmelen / 2) * (nframes + 1);

	for (int i = 0; i < outlen; i++)
		out[i] = 0;
	//反重叠相加
	for (int i = 0; i < frmelen; i++)
		out[i] = y_ifft[i*nframes + 0];
	start = frmelen - ovlplen + 1;
	int mid = start + ovlplen - 1;
	stop = start + frmelen - 1;
	for (int i = 1; i < nframes; i++)
	{
		for (int j = 0; j < (stop - start) + 1; j++)
		{
			if (j < mid)
				out[start - 1 + j] += y_ifft[(j)*nframes + i];
			else
				out[start - 1 + j] = y_ifft[(j)*nframes + i];
		}
		start = mid + 1;
		mid = start + ovlplen - 1;
		stop = start + frmelen - 1;
	}
	fftw_destroy_plan(p);
	free(win);
	free(n_spect);
	free(frame_data);
	fftw_free(x_fft);
	fftw_free(x_fft_temp);
	free(x_mag);
	free(x_ph);
	free(SNR_x);
	free(beta_x);
	free(sub_speech_x);
	fftw_free(y_fft_temp);
	free(y_ifft_temp);
	free(y_ifft);
}
void berouti(double * SNR, double * a, int nbands, int nframes)
{
	int i, j;	
	for (i = 0; i < nbands; i++)
	{
		for (j = 0; j < nframes; j++)
		{
			if (SNR[i*nframes + j] >= -5.0 && SNR[i*nframes + j] <= 20)
				a[i*nframes + j] = 4 - SNR[i*nframes + j] * 3 / 20;
			else if (SNR[i*nframes + j] < -5.0)
				a[i*nframes + j] = 4.75;
			else
				a[i*nframes + j] = 1;
		}
	}
}

void noiseupt(double * x_magsm, double * n_spect,
				 int nframes, int fftlen)
{
	int *state = (int *)malloc(nframes * sizeof(int));	//分配状态标记数组
	double x_var, n_var, judgevalue, rti = 0;
	int i = 0, j;
	/* 计算第一列数据即第一帧数据 */
	for (j = 0; j < fftlen; j++)
	{
		x_var = pow(x_magsm[j* nframes + i], 2);
		n_var = pow(n_spect[j* nframes + i], 2);
		rti += x_var / n_var - log10(x_var / n_var) - 1;
	}
	judgevalue = rti / fftlen;
	if (judgevalue <= 0.4)
	{
		for (j = 0; j < fftlen; j++)
		{
			n_spect[j* nframes + i] = sqrt(0.9*pow(n_spect[j* nframes + i], 2) + (1 - 0.9)*pow(x_magsm[j* nframes + i], 2));
		}
	}
	for (i = 1; i < nframes; i++)
	{
		rti = 0;
		for (j = 0; j < fftlen; j++)
		{
			x_var = pow(x_magsm[j* nframes + i], 2);
			n_var = pow(n_spect[j* nframes + i - 1], 2);
			rti += x_var / n_var - log10(x_var / n_var) - 1;
		}
		judgevalue = rti / fftlen;
		if (judgevalue <= 0.45)
		{
			for (j = 0; j < fftlen; j++)
			{
				n_spect[j* nframes + i] = sqrt(0.9*pow(n_spect[j* nframes + i - 1], 2) + (1 - 0.9)*pow(x_magsm[j* nframes + i], 2));
			}
		}
		else
		{
			for (j = 0; j < fftlen; j++)
			{
				n_spect[j* nframes + i] = n_spect[j* nframes + i - 1];
			}
		}
	}
	free(state);	//释放
}

void hamming(int N, double *win)
{
	int i;
	for (i = 0; i < N; i++)
	{
		win[i] = 0.54 - 0.46*cos(2 * PI *(double)i / (N - 1));
		win[i] = sqrt(win[i]);
	}
}

void frame(double *sdata, int ndata, double *windows, int nwind, int frmshift,
	int offset, int trunc, double *fdata)
{
	int nframes = (int)(ceil((double)ndata / frmshift));
	double *tdata = (double *)malloc(nwind * nframes * sizeof(double));
	int frm, ixlen, ixend, i;
	int ixstrt = 0;
	for (int i = 0; i < nwind*nframes; i++)
		tdata[i] = 0;
	for (frm = 0; frm < nframes; frm++)
	{
		ixend = min(ndata, ixstrt + nwind - 1);
		ixlen = ixend - ixstrt + 1;
		for (i = 0; i < ixlen; i++)
		{
			if (ixstrt + i < ndata)
				tdata[i*nframes + frm] = sdata[ixstrt + i];
			else
				tdata[i*nframes + frm] = 0;
		}
		ixstrt += frmshift;
	}
	scale_mtx(fdata, tdata, nframes, nwind, windows);
	free(tdata);
}

void scale_mtx(double *fdata, double *tdata, int ros, int col, double *windows)
{
	int i, j;
	for (j = 0; j < ros; j++)
	{
		for (i = 0; i < col; i++)
		{
			fdata[i*ros + j] = tdata[i*ros + j] * windows[i];
		}
	}
}