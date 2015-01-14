/*
spectrum - quick an dirty spectrum analyzer
Written in 2012 by <Ahmet Inan> <xdsopl@googlemail.com>
To the extent possible under law, the author(s) have dedicated all copyright and related and neighboring rights to this software to the public domain worldwide. This software is distributed without any warranty.
You should have received a copy of the CC0 Public Domain Dedication along with this software. If not, see <http://creativecommons.org/publicdomain/zero/1.0/>.
*/

#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <fftw3.h>
#include "window.h"
#include "trans.h"

struct stft {
	struct trans base;
	float *in, *tmp, *win;
	complex float *out;
	fftwf_plan plan;
	int bins, samples;
};

void slide_stft(struct trans *trans, float *in, int N)
{
	struct stft *stft = (struct stft *)trans->data;
	memmove(stft->in + N, stft->in, sizeof(float) * (stft->samples - N));
	memcpy(stft->in, in, sizeof(float) * N);
	for (int i = 0; i < stft->samples; i++)
		stft->tmp[i] = stft->win[i] * stft->in[i];
	fftwf_execute(stft->plan);
}

void get_stft(struct trans *trans, float *out)
{
	struct stft *stft = (struct stft *)trans->data;
	for (int i = 0; i < stft->bins; i++)
		out[i] = cabsf(stft->out[i+1]);
}

void free_stft(struct trans *trans)
{
	struct stft *stft = (struct stft *)trans->data;
	fftwf_destroy_plan(stft->plan);
	free(stft->out);
	free(stft->tmp);
	free(stft->in);
	free(stft->win);
	free(stft);
}

struct trans *create_stft(int bins)
{
	struct stft *stft = (struct stft *)malloc(sizeof(struct stft));
	stft->base.free = free_stft;
	stft->base.slide = slide_stft;
	stft->base.get = get_stft;
	stft->base.data = (void *)stft;
	stft->bins = bins;
	int samples = bins * 2;
	stft->samples = samples;
	stft->in = (float *)malloc(sizeof(float) * samples);
	stft->tmp = (float *)malloc(sizeof(float) * samples);
	memset(stft->in, 0, sizeof(float) * samples);
	stft->out = (complex float *)malloc(sizeof(complex float) * (bins + 1));
	stft->plan = fftwf_plan_dft_r2c_1d(samples, stft->tmp, stft->out, FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);
	stft->win = (float *)malloc(sizeof(float) * samples);
	float sum = 0.0f;
	for (int i = 0; i < samples; i++) {
		stft->win[i] = gauss(i, samples|1, 0.2f);
		sum += stft->win[i];
	}
	for (int i = 0; i < samples; i++)
		stft->win[i] /= sum;
	return (struct trans *)&(stft->base);
}

