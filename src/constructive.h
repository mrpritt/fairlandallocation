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
#include "objective_function.h"
#include "parameters.h"
#include "proterra.h"
#include "rivers_constraint.h"
#include "solution.h"

struct candidate {
  candidate(int lot = -1, int cell = -1) : lot(lot), cell(cell) {}
  int lot = -1, cell = -1;

  bool is_invalid() const { return cell == -1; }
  void set_invalid() { cell = -1; }

  ll value(const solution& s) const {
    if (of.last_updated < s.num_assigned) {
      auto& t = const_cast<candidate&>(*this);
      const_cast<solution&>(s).of.calc_swap(-1, lot, cell, t.of);
      t.of.last_updated = s.num_assigned;
    }
    assert(of.value != -1);
    return of.value;
  }

  int river_value(const solution& s) const {
    if (rc.last_updated < s.num_assigned) {
      auto& t = const_cast<candidate&>(*this);
      auto& ss = const_cast<solution&>(s);
      ss.rc.calc_swap(-1, lot, cell, t.rc, ss);
      t.rc.last_updated = s.num_assigned;
    }
    assert(rc.value != -1);
    return rc.value;
  }

  int size_ratio(const solution& s) const {
    int ba = (s.area[lot] == s.area[s.big_lot()] ? s.area[lot] + 1
                                                 : s.area[s.big_lot()]);
    int sa = s.area[s.small_lot()];
    if (lot == s.small_lot() and pt::lots > 1 and
        s.area[s.rc.sa[1]] > s.area[lot])
      ++sa;
    if (sa == 0) return numeric_limits<int>::max();
    int sr = (ba * 1000) / sa;
    return sr > 1000 * prm::maximum_size_ratio ? sr : 0;
  }

  objective_function::candidate of;
  rivers_constraint::candidate rc;
};

struct Constructive {

  static void construct_from_seeds(solution& s, const vi& initial_pos,
                                   bool use_alpha = false);

  static void construct(solution& s, bool use_alpha = false);

  static bool compare_candidates(const candidate& c1, const candidate& c2,
                                 const solution& s);
};
