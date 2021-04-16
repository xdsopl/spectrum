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
#include "pcm.h"
#include "trans.h"

#define BINS (512)
#define STEP (64)
static int pause;
static int rainbow = 1;
static int log_power = 1;
static int log_freq = 1;
static int norm_vis = 0;

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
					case SDLK_l:
						log_power ^= 1;
						break;
					case SDLK_o:
						log_freq ^= 1;
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
	float K0 = 0.04045f;
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
		float r = v * fminf(fmaxf(4.0f * v - 2.0f, 0.0f), 1.0f);
		float g = v * fminf(fmaxf(1.0f - 4.0f * fabsf(v - 0.5f), 0.0f), 1.0f);
		float b = v * fminf(fmaxf(2.0f - 4.0f * v, 0.0f), 1.0f);
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
	float *in = (float *)malloc(sizeof(float) * STEP);
	float *out = (float *)malloc(sizeof(float) * BINS);

	fprintf(stderr, "creating stft ..");
	struct trans *stft = create_stft(BINS);
	fprintf(stderr, " done.\ncreating cqt, please be patient ..");
	struct trans *cqt = create_cqt(BINS, rate, 120.0f, rate / 2.0f);
	fprintf(stderr, " done.\n");

	SDL_Init(SDL_INIT_VIDEO);
	SDL_Surface *screen = SDL_SetVideoMode(1024, BINS, 32, SDL_HWSURFACE|SDL_DOUBLEBUF);
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

		if (!read_pcm(pcm, buff, STEP))
			memset(buff, 0, sizeof(short) * channels * STEP);

		for (int i = 0; i < STEP; i++)
			in[i] = (float)buff[((STEP-1)-i) * channels] / 32768.0f;


		struct trans *trans = stft;
		if (log_freq)
			trans = cqt;

#if 0
		slide_trans(trans, in, STEP);
#else
		slide_trans(stft, in, STEP);
		slide_trans(cqt, in, STEP);
#endif
		get_trans(trans, out);

		float max_tmp = 0.0;
		for (int i = 0; i < max_hist_size; i++)
			max_tmp = max_tmp > max_hist[i] ? max_tmp : max_hist[i];
		float rec_max = 1.0 / (max_tmp ? max_tmp : 1.0);

		float max_power = 0.0;
		for (int i = 0; i < BINS; i++) {
			float power = out[i] * out[i];
			max_power = max_power < power ? power : max_power;
			if (norm_vis)
				power *= rec_max;
			float decibel = 10.0f * log10f(fmaxf(0.000001f, power));
			float val = power;
			if (log_power)
				val = 1.0f + fminf(fmaxf(decibel, -60.0f), 0.0f) / 60.0f;
			fbp[w * i + w - 1] = val_rgb(val);
		}
		max_hist[max_hist_last] = max_power;
		max_hist_last = (max_hist_last + 1) % max_hist_size;
	}
	return 0;
}

