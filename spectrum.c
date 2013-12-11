/*
spectrum - quick an dirty spectrum analyzer
Written in 2012 by <Ahmet Inan> <xdsopl@googlemail.com>
To the extent possible under law, the author(s) have dedicated all copyright and related and neighboring rights to this software to the public domain worldwide. This software is distributed without any warranty.
You should have received a copy of the CC0 Public Domain Dedication along with this software. If not, see <http://creativecommons.org/publicdomain/zero/1.0/>.
*/

#include <math.h>
#include <complex.h>
#include <stdlib.h>
#include <stdio.h>
#include <SDL.h>
#include <fftw3.h>
#include "pcm.h"
#include "window.h"

#define BINS (512)
#define STEP (8)
static int pause;
static int rainbow = 1;
static int rms_comp = 1;
static int norm_vis = 1;
static int log_vis = 1;

void handle_events()
{
	SDL_Event event;
	while (SDL_PollEvent(&event)) {
		switch (event.type) {
			case SDL_KEYDOWN:
				switch (event.key.keysym.sym) {
					case SDLK_q:
						exit(0);
						break;
					case SDLK_c:
						rms_comp ^= 1;
						break;
					case SDLK_l:
						log_vis ^= 1;
						break;
					case SDLK_n:
						norm_vis ^= 1;
						break;
					case SDLK_r:
						rainbow ^= 1;
						break;
					case SDLK_SPACE:
						pause ^= 1;
						break;
					case SDLK_ESCAPE:
						exit(0);
						break;
					default:
						break;
				}
				break;
			case SDL_QUIT:
				exit(1);
				break;
			default:
				break;
		}
	}

}

uint32_t srgb(float r, float g, float b)
{
	r = fminf(fmaxf(r, 0.0f), 1.0f);
	g = fminf(fmaxf(g, 0.0f), 1.0f);
	b = fminf(fmaxf(b, 0.0f), 1.0f);
#if 1
	float K0 = 0.03928f;
	float a = 0.055f;
	float phi = 12.92f;
	float gamma = 2.4f;
	r = r <= K0 / phi ? r * phi : (1.0f + a) * powf(r, 1.0f / gamma) - a;
	g = g <= K0 / phi ? g * phi : (1.0f + a) * powf(g, 1.0f / gamma) - a;
	b = b <= K0 / phi ? b * phi : (1.0f + a) * powf(b, 1.0f / gamma) - a;
#endif
	return (int)(255.0f * r) << 16 |
		(int)(255.0f * g) << 8 |
		(int)(255.0f * b);
}

uint32_t val_rgb(float v)
{
	uint32_t rgb = 0;
	if (rainbow) {
		float r = 4.0f * v - 2.0f;
		float g = 2.0f - 4.0f * fabsf(v - 0.5f);
		float b = 2.0f - 4.0f * v;
		rgb = srgb(r, g, b);
	} else {
		rgb = srgb(v, v, v);
	}
	return rgb;
}

int main(int argc, char **argv)
{
	struct pcm *pcm;
	char *pcm_name = "default";
	if (argc != 1)
		pcm_name = argv[1];

	if (!open_pcm_read(&pcm, pcm_name))
		return 1;

	info_pcm(pcm);

	float rate = rate_pcm(pcm);
	int channels = channels_pcm(pcm);
	if (channels > 1)
		fprintf(stderr, "using first of %d channels\n", channels);

	short *buff = (short *)malloc(sizeof(short) * channels * STEP);

	atexit(SDL_Quit);
	SDL_Init(SDL_INIT_VIDEO);
	SDL_Surface *screen = SDL_SetVideoMode(1024, BINS / 2, 32, SDL_HWSURFACE|SDL_DOUBLEBUF);
	if (NULL == screen)
		exit(1);
	if (screen->format->BytesPerPixel != 4)
		exit(1);
	uint32_t *fbp = (uint32_t *)screen->pixels;
	int w = screen->w;
	int h = screen->h;
	memset(fbp, 0, sizeof(uint32_t) * w * h);

	for (int i = 0; i < h; i++) {
		for (int k = 0; k < w / 2; k++) {
			fbp[w * i + w / 4 + k] = val_rgb((float)i / (float)(h - 1));
		}
	}

	SDL_WM_SetCaption("Spectrum", "spectrum");
	SDL_EnableKeyRepeat(SDL_DEFAULT_REPEAT_DELAY, SDL_DEFAULT_REPEAT_INTERVAL);

	float *inp = (float *)malloc(sizeof(float) * BINS);
	float *win = (float *)malloc(sizeof(float) * BINS);
	float *tmp = (float *)malloc(sizeof(float) * BINS);
	memset(tmp, 0, sizeof(float) * BINS);
	complex float *out = (complex float *)malloc(sizeof(complex float) * BINS / 2);
	fftwf_plan plan = fftwf_plan_dft_r2c_1d(BINS, inp, out, FFTW_ESTIMATE | FFTW_DESTROY_INPUT);

	for (int i = 0; i < BINS; i++)
//		win[i] = rect(i, BINS|1, 2.0);
		win[i] = gauss(i, BINS|1, 0.1);
//		win[i] = lanczos(i, BINS|1, 0.0);
//		win[i] = hann(i, BINS|1, 0.0);
//		win[i] = hamming(i, BINS|1, 0.0);
//		win[i] = kaiser(i, BINS|1, 3.0);

	float rms_hist_seconds = 0.5;
	int rms_hist_size = (rms_hist_seconds * rate) / STEP;
	float *rms_hist = (float *)malloc(sizeof(float) * rms_hist_size);
	int rms_hist_last = 0;
	for (int i = 0; i < rms_hist_size; i++)
		rms_hist[i] = 1.0;

	float max_hist_seconds = 0.5;
	int max_hist_size = (max_hist_seconds * rate) / STEP;
	float *max_hist = (float *)malloc(sizeof(float) * max_hist_size);
	int max_hist_last = 0;
	for (int i = 0; i < max_hist_size; i++)
		max_hist[i] = 0.0;

	uint64_t cnt = 0;
	for (;;) {
		if (cnt++ > rate / (STEP * 50)) {
			cnt = 0;
			SDL_Flip(screen);
			SDL_Delay(10);
		}

		handle_events();

		if (pause)
			continue;

		for (int i = 0; i < h; i++)
			memmove(fbp + w * i, fbp + w * i + 1, sizeof(uint32_t) * (w - 1));

		if (!read_pcm(pcm, buff, STEP)) {
			close_pcm(pcm);
			exit(0);
		}

		memmove(tmp + STEP, tmp, sizeof(float) * (BINS - STEP));

		for (int i = 0; i < STEP; i++)
			tmp[i] = (float)buff[((STEP-1)-i) * channels] / 32768.0;

		rms_hist[rms_hist_last] = 0;
		for (int i = 0; i < STEP; i++)
			rms_hist[rms_hist_last] += tmp[i] * tmp[i];
		rms_hist_last = (rms_hist_last + 1) % rms_hist_size;

		float rms_sum = 0.0;
		for (int i = 0; i < rms_hist_size; i++)
			rms_sum += rms_hist[i];

		float rms_tmp = sqrtf(rms_sum / (float)(rms_hist_size * STEP));
		float rec_rms = 1.0 / (rms_tmp ? rms_tmp : 1.0);

		for (int i = 0; i < BINS; i++)
			inp[i] = tmp[i] * win[i];

		if (rms_comp) {
			for (int i = 0; i < BINS; i++)
				inp[i] *= rec_rms;
		}

		fftwf_execute(plan);

		float max_tmp = 0.0;
		for (int i = 0; i < max_hist_size; i++)
			max_tmp = max_tmp > max_hist[i] ? max_tmp : max_hist[i];
		float rec_max = 1.0 / (max_tmp ? max_tmp : 1.0);

		float max_power = 0.0;
		for (int i = 0; i < (BINS / 2); i++) {
			float amp = cabsf(out[i]) / (float)BINS;
			float power = amp * amp;
			max_power = max_power < power ? power : max_power;
			if (norm_vis)
				power *= rec_max;
			float decibel = 10.0 * log10f(power);
			float val = 0;
			if (log_vis)
				val = 1.0 + fminf(fmaxf(decibel, -100.0), 0.0) / 100.0;
			else
				val = power;
			fbp[w * i + w - 1] = val_rgb(val);
		}
		max_hist[max_hist_last] = max_power;
		max_hist_last = (max_hist_last + 1) % max_hist_size;
	}
	return 0;
}

