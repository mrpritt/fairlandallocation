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
#include <cassert>
#include <ciso646>
#include <string>
#include <vector>
#include "convexity.h"
#include "defines.h"
#include "objective_function.h"
#include "proterra.h"
#include "rivers_constraint.h"
#include "util.h"

struct candidate;

struct solution {
	solution() = default;
	solution(const solution &) = default;
	solution(solution &&) = default;
	solution &operator=(const solution &) = default;
	solution &operator=(solution &&) = default;

	void init();

	void populate(vi assigned);

	bool is_assigned(int c) const { return assigned[c] != -1; }

	void write_to_png(const string &filename) const;

	bool operator<(const solution &s) const {
		int r = river_value() - s.river_value();
		if (r == 0) {
			double sr = size_ratio() - s.size_ratio();
			if (eps_eq(sr, 0.0))
				return value() < s.value();
			return sr < 0;
		}
		return r < 0;
	}

	int do_swaps(vector<candidate> &candidates);

	ll value() const;

	int river_value() const;

	vi area;

	inline int small_lot() const {
#ifdef HARD_DEBUG
		assert(area[rc.sa[0]] == min(area));
#endif
		return rc.sa[0];
	}

	inline int big_lot() const {
#ifdef HARD_DEBUG
		assert(area[rc.sa[pt::lots - 1]] == max(area));
#endif
		return rc.sa[pt::lots - 1];
	}

	double size_ratio() const;

	vi num_river;

	vi assigned;

	int num_assigned = 0;

	bool check_border_brute_force(int c) const;

	objective_function of;

	rivers_constraint rc;
};