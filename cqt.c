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

struct cqt {
	struct trans base;
	float *in;
	complex float *out, *kern;
	int *first, *length, *kern_off;
	fftwf_plan plan;
	int bins, samples;
};

void slide_cqt(struct trans *trans, float *in, int N)
{
	struct cqt *cqt = (struct cqt *)trans->data;
	memmove(cqt->in + N, cqt->in, sizeof(float) * (cqt->samples - N));
	memcpy(cqt->in, in, sizeof(float) * N);
	fftwf_execute(cqt->plan);
}

void get_cqt(struct trans *trans, float *out)
{
	struct cqt *cqt = (struct cqt *)trans->data;
	for (int i = 0; i < cqt->bins; i++) {
		complex float sum = 0.0f;
		for (int n = 0; n < cqt->length[i]; n++)
			sum += cqt->out[cqt->first[i] + n] * cqt->kern[cqt->kern_off[i] + n];
		out[i] = cabsf(sum);
	}
}

void free_cqt(struct trans *trans)
{
	struct cqt *cqt = (struct cqt *)trans->data;
	fftwf_destroy_plan(cqt->plan);
	free(cqt->out);
	free(cqt->in);
	free(cqt->kern);
	free(cqt->first);
	free(cqt->length);
	free(cqt->kern_off);
	free(cqt);
}

struct trans *create_cqt(int bins, float rate, float f_min, float f_max)
{
	struct cqt *cqt = (struct cqt *)malloc(sizeof(struct cqt));
	cqt->base.free = free_cqt;
	cqt->base.slide = slide_cqt;
	cqt->base.get = get_cqt;
	cqt->base.data = (void *)cqt;
	cqt->bins = bins;
	float ratio = expf(logf(f_max / f_min) / (bins - 1.0f));
	int samples = 1 << (int)ceil(log2(rate / (f_min * (ratio - 1.0f))));
	cqt->samples = samples;
#if 0
	fprintf(stderr, "\n");
	for (int i = 0; i < bins; i++) {
		float f_lin = f_max * i / (bins - 1.0f);
		float f_log = f_min * powf(ratio, i);
		int N = rate / (f_log * (ratio - 1.0f));
		fprintf(stderr, "% 4d % 10.3f % 10.3f % 5d\n", i, f_lin, f_log, N);
	}
	fprintf(stderr, "f_min = %g f_max = %g ratio = %g samples = %d\n", f_min, f_max, ratio, samples);
	fprintf(stderr, "latency = %g seconds\n", 0.5f * samples / rate);
#endif
	cqt->first = (int *)malloc(sizeof(int) * bins);
	cqt->length = (int *)malloc(sizeof(int) * bins);
	cqt->kern_off = (int *)malloc(sizeof(int) * bins);
	int kern_size = 4 * samples;
	cqt->kern = (complex float *)malloc(sizeof(complex float) * kern_size);

	complex float *kern_in = (complex float *)malloc(sizeof(complex float) * samples);
	complex float *kern_out = (complex float *)malloc(sizeof(complex float) * samples);
	fftwf_plan kern_plan = fftwf_plan_dft_1d(samples, kern_in, kern_out, FFTW_FORWARD, FFTW_ESTIMATE | FFTW_DESTROY_INPUT);

	for (int i = 0; i < bins; i++) {
		memset(kern_in, 0, sizeof(complex float) * samples);
		float fn = f_min * powf(ratio, i);
		int N = fminf(samples-1, rate / (fn * (ratio - 1.0f)));
		N |= 1;
		float sum = 0.0f;
		for (int n = 0; n < N; n++) {
			float w = hamming(n, N, 0.0f);
			sum += w;
			complex float z = cexpf(I * 2.0f * M_PI * n * fn / rate);
			kern_in[samples / 2 - (N-1) / 2 + n] = w * z;
		}
		for (int n = 0; n < N; n++)
			kern_in[samples / 2 - (N-1) / 2 + n] /= sum;
		fftwf_execute(kern_plan);
#if 0
		if (i == 0) {
			for (int n = 0; n < samples; n++)
				printf("%d %g %g %g\n", n, 10.0f * log10f(cabsf(kern_out[n])), crealf(kern_out[n]), cimagf(kern_out[n]));
			return 0;
		}
#endif
		int nmax = 0;
		for (int n = 0; n < samples; n++)
			nmax = cabsf(kern_out[nmax]) < cabsf(kern_out[n]) ? n : nmax;
		int n0 = nmax, n1 = nmax;
		while (0 < n0 && cabsf(kern_out[n0-1]) <= cabsf(kern_out[n0]))
			n0--;
		while (n1 < samples-1 && cabsf(kern_out[n1+1]) <= cabsf(kern_out[n1]))
			n1++;

		cqt->first[i] = n0;
		n1 = fminf(samples / 2 - 1, n1);
		cqt->length[i] = n1 - n0 + 1;

		if (!i)
			cqt->kern_off[0] = 0;
		else
			cqt->kern_off[i] = cqt->kern_off[i-1] + cqt->length[i-1];

		if (cqt->kern_off[i] + cqt->length[i] > kern_size) {
			kern_size *= 2;
			cqt->kern = (complex float *)realloc(cqt->kern, sizeof(complex float) * kern_size);
		}

		for (int n = 0; n < cqt->length[i]; n++)
			cqt->kern[cqt->kern_off[i] + n] = conjf(kern_out[cqt->first[i] + n]) / samples;
	}
	// fprintf(stderr, "%d coefs\n", cqt->kern_off[bins-1] + cqt->length[bins-1]);

	fftwf_destroy_plan(kern_plan);
	free(kern_in);
	free(kern_out);

	cqt->in = (float *)malloc(sizeof(float) * samples);
	memset(cqt->in, 0, sizeof(float) * samples);
	cqt->out = (complex float *)malloc(sizeof(complex float) * samples / 2);
	cqt->plan = fftwf_plan_dft_r2c_1d(samples, cqt->in, cqt->out, FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);
	return (struct trans *)&(cqt->base);
}

