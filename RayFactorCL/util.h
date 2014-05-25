/*
Copyright 2010-2011, D. E. Shaw Research.
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

* Redistributions of source code must retain the above copyright
  notice, this list of conditions, and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions, and the following disclaimer in the
  documentation and/or other materials provided with the distribution.

* Neither the name of D. E. Shaw Research nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
#ifndef UTIL_H__
#define UTIL_H__
#include <string.h>
#include "Random123/features/compilerfeatures.h"

#if defined(_MSC_VER)
#include <Windows.h>
R123_STATIC_INLINE double now(){
    LARGE_INTEGER f; // ticks per second
    LARGE_INTEGER t;
    QueryPerformanceFrequency(&f);
    QueryPerformanceCounter(&t); // ticks since epoch
    return ((double)t.QuadPart)/((double)f.QuadPart);
}
#else // _MSC_VER
#include <sys/time.h>
R123_STATIC_INLINE double now(){
    struct timeval tv; 
    gettimeofday(&tv, 0); 
    return 1.e-6*tv.tv_usec + tv.tv_sec;
}
#endif // _MSC_VER

char *nameclean(char *s);


/* timer returns difference between current time and *d, also updates *d with current time. */
R123_STATIC_INLINE double timer(double *d) {
    double dold = *d;
    *d = now();
    return *d - dold;
}

#define WHITESPACE " \t\f\v\n"
char *nameclean(char *s)
{
    char *cp = s, *cp2, *cpend = s + strlen(s);
    size_t i;

    cp2 = s;
    while ((cp2 = strstr(cp2, "(R)")) != NULL) {
	*cp2++ = ' ';
	*cp2++ = ' ';
	*cp2++ = ' ';
    }
    cp2 = s;
    while ((cp2 = strstr(cp2, "CPU")) != NULL) {
	*cp2++ = ' ';
	*cp2++ = ' ';
	*cp2++ = ' ';
    }
    cp2 = s;
    while ((cp2 = strchr(cp2, '@')) != NULL) {
	*cp2++ = ' ';
    }
    while ((i = strcspn(cp, WHITESPACE)) > 0) {
	cp += i;
	i = strspn(cp, WHITESPACE);
	if (i > 0) {
	    cp2 = cp + i;
	    *cp++ = ' ';
	    if (cp2 > cp) {
		memmove(cp, cp2, cpend - cp);
		cpend -= cp2 - cp;
	    }
	}
    }
    return s;
}

#undef WHITESPACE

#endif /* UTIL_H__ */
