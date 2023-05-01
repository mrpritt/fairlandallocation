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
#include "objective_function.h"
#include <cassert>
#include <iostream>
#include <numeric>
#include "constructive.h"
#include "proterra.h"
#include "solution.h"
#include "util.h"

void objective_function::init() {
	val.assign(pt::lots, 0);
	value = sum_xi = sum_xi_sq = 0;
}

void objective_function::calc_swap(int from, int to, int c, candidate &sp) {
	int xj = pt::val[c];
	sp.sum_xi = sum_xi;
	sp.sum_xi_sq = sum_xi_sq;
	if (from != -1) {
		ll vf = val[from];
		sp.sum_xi_sq += -(vf * vf) + (vf - xj) * (vf - xj);
		sp.sum_xi -= xj;
	}
	if (to != -1) {
		ll vt = val[to];
		sp.sum_xi_sq += -(vt * vt) + (vt + xj) * (vt + xj);
		sp.sum_xi += xj;
	}
	sp.value = sp.sum_xi_sq - (sp.sum_xi * sp.sum_xi) / pt::lots;
#ifdef HARD_DEBUG
	if (from != -1)
		val[from] -= xj;
	if (to != -1)
		val[to] += xj;
	assert(sp.sum_xi_sq >= 0);
	assert_value_acceptable(sp.value);
	if (from != -1)
		val[from] += xj;
	if (to != -1)
		val[to] -= xj;
#endif
}

void objective_function::do_swap(int from, int to, int c, candidate &sp) {
	if (from != -1)
		val[from] -= pt::val[c];
	if (to != -1)
		val[to] += pt::val[c];
	assert(sum_xi != -1);
	assert(sum_xi_sq != -1);
	sum_xi = sp.sum_xi;
	sum_xi_sq = sp.sum_xi_sq;
	value = sp.value;
	assert(sum_xi != -1);
	assert(sum_xi_sq != -1);
	assert_value_acceptable(sp.value);
}

void objective_function::do_swaps(vector<::candidate> &cands) {
	if ((int)cands.size() < pt::lots) {
		// 	if (false) {
		for (int i = 0; i < (int)cands.size(); ++i) {
			if (cands[i].is_invalid())
				continue;
			ll xj = pt::val[cands[i].cell];
			assert(cands[i].lot != -1);
			ll vt = val[cands[i].lot];
			sum_xi_sq += -(vt * vt) + (vt + xj) * (vt + xj);
			sum_xi += xj;
			val[cands[i].lot] += xj;
		}
	} else {
		for (int i = 0; i < (int)cands.size(); ++i) {
			if (cands[i].is_invalid())
				continue;
			val[cands[i].lot] += pt::val[cands[i].cell];
		}
		sum_xi = sum_xi_sq = 0;
		for (int i = 0; i < pt::lots; ++i) {
			ll v = val[i];
			sum_xi += v;
			sum_xi_sq += v * v;
		}
	}
	value = sum_xi_sq - (sum_xi * sum_xi) / pt::lots;
	assert_value_acceptable();
}

ll objective_function::compute_value_brute_force() const {
	double sum = accumulate(val.begin(), val.end(), double(0));
	double mean = sum / (double)pt::lots;
	double var = 0;
	for (int i = 0; i < pt::lots; ++i) {
		double x = (double)val[i] - mean;
		var += x * x;
	}
	return (ll)(var);
}

void objective_function::populate(const solution &s) {
	val.assign(pt::lots, 0);
	for (int i = 0; i < pt::nland; ++i)
		if (s.assigned[i] != -1)
			val[s.assigned[i]] += pt::val[i];
	sum_xi = sum_xi_sq = 0;
	for (int i = 0; i < pt::lots; ++i) {
		ll v = val[i];
		sum_xi += v;
		sum_xi_sq += v * v;
	}
	value = sum_xi_sq - (sum_xi * sum_xi / pt::lots);
#ifdef HARD_DEBUG
	assert_value_acceptable();
#endif
}
