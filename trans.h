/*
spectrum - quick an dirty spectrum analyzer
Written in 2012 by <Ahmet Inan> <xdsopl@googlemail.com>
To the extent possible under law, the author(s) have dedicated all copyright and related and neighboring rights to this software to the public domain worldwide. This software is distributed without any warranty.
You should have received a copy of the CC0 Public Domain Dedication along with this software. If not, see <http://creativecommons.org/publicdomain/zero/1.0/>.
*/


#ifndef TRANS_H
#define TRANS_H

struct trans {
	void (*free)(struct trans *);
	void (*slide)(struct trans *, float *, int);
	void (*get)(struct trans *, float *);
	void *data;
};

void free_trans(struct trans *);
void slide_trans(struct trans *, float *, int);
void get_trans(struct trans *, float *);
struct trans *create_stft(int);
struct trans *create_cqt(int, float, float, float);

#endif

