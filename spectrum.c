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
static int rms_comp = 1;

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

uint32_t val_rgb(float v)
{
#if 1
	int R = 255.0 * fminf(fmaxf(4.0 * v - 2.0, 0.0), 1.0);
	int G = 255.0 * fminf(fmaxf(2.0 - 4.0 * fabsf(v - 0.5), 0.0), 1.0);
	int B = 255.0 * fminf(fmaxf(2.0 - 4.0 * v, 0.0), 1.0);
	return (R << 16) | (G << 8) | B;
#else
	return (int)(255.0 * fminf(fmaxf(v, 0.0), 1.0)) * 0x00010101;
#endif
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
//		win[i] = rect(i, BINS, 2.0);
		win[i] = gauss(i, BINS, 0.1);
//		win[i] = lanczos(i, BINS, 0.0);
//		win[i] = hann(i, BINS, 0.0);
//		win[i] = hamming(i, BINS, 0.0);
//		win[i] = kaiser(i, BINS, 3.0);

	float rms_hist_seconds = 0.5;
	int rms_hist_size = rms_hist_seconds * rate;
	float *rms_hist = (float *)malloc(sizeof(float) * rms_hist_size);
	int rms_hist_last = 0;
	for (int i = 0; i < rms_hist_size; i++)
		rms_hist[i] = 1.0;

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

		for (int i = 0; i < STEP; i++, rms_hist_last = (rms_hist_last + 1) % rms_hist_size)
			rms_hist[rms_hist_last] = tmp[i] * tmp[i];

		float rms_sum = 0.0;
		for (int i = 0; i < rms_hist_size; i++)
			rms_sum += rms_hist[i];

		float rms = sqrtf(rms_sum / (float)rms_hist_size);
		float rec_rms = 1.0 / fmax(0.01, rms);

		for (int i = 0; i < BINS; i++)
			inp[i] = tmp[i] * win[i];

		if (rms_comp) {
			for (int i = 0; i < BINS; i++)
				inp[i] *= rec_rms;
		}

		fftwf_execute(plan);

		for (int i = 0; i < (BINS / 2); i++) {
			float amp = cabsf(out[i]) / (float)BINS;
			float power = powf(amp, 2.0);
			float decibel = 10.0 * log10f(power);
			float val = 1.0 + fminf(fmaxf(decibel, -100.0), 0.0) / 100.0;
			fbp[w * i + w - 1] = val_rgb(val);
		}
	}
	return 0;
}

