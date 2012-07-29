
CFLAGS = -g -D_GNU_SOURCE=1 -W -Wall -O3 -std=c99 -ffast-math $(shell sdl-config --cflags) $(shell pkg-config fftw3f --cflags)
LDFLAGS = -lm -lasound $(shell sdl-config --libs) $(shell pkg-config fftw3f --libs)

all: spectrum

clean:
	rm -f spectrum *.o

spectrum: spectrum.o mmap_file.o pcm.o wav.o alsa.o window.o

