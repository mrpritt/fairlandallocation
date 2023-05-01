/*
* A genetic algorithm for fair land allocation
* Copyright (c) 2017 Alex Gliesch, Marcus Ritt, Mayron C. O. Moreira
*
* Permission is hereby granted, free of charge, to any person (the "Person")
* obtaining a copy of this software and associated documentation files (the
* "Software"), to deal in the Software, including the rights to use, copy, modify,
* merge, publish, distribute the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following conditions:
*
* 1. The above copyright notice and this permission notice shall be included in
*    all copies or substantial portions of the Software.
* 2. Under no circumstances shall the Person be permitted, allowed or authorized
*    to commercially exploit the Software.
* 3. Changes made to the original Software shall be labeled, demarcated or
*    otherwise identified and attributed to the Person.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
* FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
* COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
* IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
* CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/
#pragma once
#include <cmath>
#include <cstdarg>
#include <numeric>
#include "defines.h"
#include "fmt/printf.h"

static const int dr[] = {-1, 1, 0, 0, -1, -1, 1, 1};
static const int dc[] = {0, 0, -1, 1, -1, 1, -1, 1};

inline bool eps_eq(double a, double b) { return std::abs(a - b) < EPS; }

inline void empty_func(const char*, ...) {}

#ifdef NDEBUG
#define pr empty_func
#define prf empty_func
#else
#define pr fmt::print
#define prf fmt::printf
#endif

using fmt::format;

template <typename T>
T sum(const vector<T>& v) {
	return accumulate(v.begin(), v.end(), (T)0);
}

template <typename T>
T min(const vector<T>& v) {
	return v.size() ? *min_element(v.begin(), v.end()) : 0;
}

template <typename T>
size_t min_index(const vector<T>& v) {
	return v.size() ? min_element(v.begin(), v.end()) - v.begin() : -1;
}

template <typename T>
size_t max_index(const vector<T>& v) {
	return v.size() ? max_element(v.begin(), v.end()) - v.begin() : -1;
}

template <typename T>
T max(const vector<T>& v) {
	return v.size() ? *max_element(v.begin(), v.end()) : 0;
}

template <typename T>
double avg(const vector<T>& v) {
	return v.size() ? double(sum(v)) / double(v.size()) : 0;
}

template <typename T>
double std_dev(const vector<T>& v) {
	if (v.size() == 0)
		return 0;
	double a = avg(v), s = 0;
	for (auto t : v)
		s += (t - a) * (t - a);
	return sqrt(s / (double)v.size());
}

template <typename T>
void pr_vec(const string& intro, const vector<T>& v, const string& sep) {
	pr("{} ", intro);
	for (const auto& i : v)
		pr("{}{}", i, sep);
	if (sep != "\n")
		pr("\n");
}

// credits to http://web.mit.edu/storborg/Public/hsvtorgb.c
inline void hsv_to_rgb(uchar h, uchar s, uchar v, uchar& r, uchar& g, uchar& b) {
	uchar region, fpart, p, q, t;
	if (s == 0) {
		/* color is grayscale */
		r = g = b = v;
		return;
	}
	/* make hue 0-5 */
	region = h / 43;
	/* find remainder part, make it from 0-255 */
	fpart = (h - (region * 43)) * 6;
	/* calculate temp vars, doing integer multiplication */
	p = (v * (255 - s)) >> 8;
	q = (v * (255 - ((s * fpart) >> 8))) >> 8;
	t = (v * (255 - ((s * (255 - fpart)) >> 8))) >> 8;
	/* assign temp vars based on color cone region */
	switch (region) {
		case 0:
			r = v;
			g = t;
			b = p;
			break;
		case 1:
			r = q;
			g = v;
			b = p;
			break;
		case 2:
			r = p;
			g = v;
			b = t;
			break;
		case 3:
			r = p;
			g = q;
			b = v;
			break;
		case 4:
			r = t;
			g = p;
			b = v;
			break;
		default:
			r = v;
			g = p;
			b = q;
			break;
	}
}