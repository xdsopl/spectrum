
CFLAGS = -g -D_GNU_SOURCE=1 -W -Wall -O3 -std=c99 -fno-math-errno -ffinite-math-only -fno-rounding-math -fno-signaling-nans -fno-trapping-math -fcx-limited-range -fsingle-precision-constant $(shell sdl-config --cflags) $(shell pkg-config fftw3f --cflags)
LDLIBS = -lm -lasound $(shell sdl-config --libs) $(shell pkg-config fftw3f --libs)

all: spectrum

clean:
	rm -f spectrum *.o

spectrum: spectrum.o mmap_file.o pcm.o wav.o alsa.o window.o stft.o cqt.o trans.o

