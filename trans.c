/*
spectrum - quick an dirty spectrum analyzer
Written in 2012 by <Ahmet Inan> <xdsopl@googlemail.com>
To the extent possible under law, the author(s) have dedicated all copyright and related and neighboring rights to this software to the public domain worldwide. This software is distributed without any warranty.
You should have received a copy of the CC0 Public Domain Dedication along with this software. If not, see <http://creativecommons.org/publicdomain/zero/1.0/>.
*/

#include "trans.h"

void free_trans(struct trans *trans)
{
	trans->free(trans);
}

void slide_trans(struct trans *trans, float *in, int N)
{
	trans->slide(trans, in, N);
}

void get_trans(struct trans *trans, float *out)
{
	trans->get(trans, out);
}

