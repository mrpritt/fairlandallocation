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
#include "defines.h"
#include <vector>

struct solution;
struct candidate;

struct rivers_constraint {
  rivers_constraint() = default;
  rivers_constraint(rivers_constraint&&) = default;
  rivers_constraint& operator=(const rivers_constraint&) = default;
  rivers_constraint& operator=(rivers_constraint&&) = default;

  struct candidate {
    friend rivers_constraint;
    int vio = -1, value = -1;
    int last_updated = -1;
    bool operator<(const candidate& sp) const { return value < sp.value; }
    bool operator==(const candidate& sp) const { return value == sp.value; }

  private:
    int nx = -1, s1 = -1, s2 = -1, s3 = -1, s4 = -1;
  };

  void init();

  void populate(const solution& s);

  void calc_swap(int from, int to, int c, candidate& sp, solution& s);

  void do_swaps(vector<::candidate>& cands, const solution& s);

  int get_value_brute_force(const solution& s) const {
    int nx, vio, cost;
    compute_value_brute_force(nx, vio, cost, s);
    return cost;
  }

  void compute_value_brute_force(int& nx, int& vio, int& cost,
                                 const solution& s) const;

  void compute_value_brute_force(const solution& s) {
    compute_value_brute_force(x, violations, value, s);
  }

  int x;

  int violations;

  int value;

  vi sa, index_sa;
};
