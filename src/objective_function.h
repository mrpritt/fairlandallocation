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
#include "assert.h"
#include "defines.h"
#include <cmath>
#include <vector>

struct solution;
struct candidate;

struct objective_function {
  objective_function() = default;
  objective_function(objective_function&&) = default;
  objective_function& operator=(const objective_function&) = default;
  objective_function& operator=(objective_function&&) = default;

  struct candidate {
    ll value = -1, sum_xi = -1, sum_xi_sq = -1;
    int last_updated = -1;
    bool operator<(const candidate& sp) const { return value < sp.value; }
    bool operator==(const candidate& sp) const { return value == sp.value; }
  };

  void init();

  void calc_swap(int from, int to, int c, candidate& sp);

  void do_swap(int from, int to, int c, candidate& sp);

  void do_swaps(vector<::candidate>& candidates);

  ll compute_value_brute_force() const;

  void assert_value_acceptable() const {
    return assert_value_acceptable(value);
  }
  void assert_value_acceptable(ll value) const {
    (void)value;
    assert(std::abs(value - compute_value_brute_force()) < 2);
  }

  void populate(const solution& s);

  ll value = 0;
  ll sum_xi_sq = 0;
  ll sum_xi = 0;
  vi val;
};
