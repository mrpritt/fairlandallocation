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
#include "constructive.h"
#include "random.h"
#include "statistics.h"
#include "util.h"
#include <cassert>
#include <functional>

Constructive cons;

void Constructive::construct_from_seeds(solution& s, const vi& initial_pos,
                                        bool use_alpha /*= false*/) {
#ifdef HARD_DEBUG
  assert((int)initial_pos.size() == pt::lots);
  for (int i = 0; i < pt::lots; ++i)
    for (int j = 0; j < pt::lots; ++j)
      assert(i == j or initial_pos[i] != initial_pos[j]);
#endif
  static vector<candidate> initial_cands;
  initial_cands.clear();
  for (int i = 0; i < pt::lots; ++i)
    initial_cands.emplace_back(i, initial_pos[i]);
  s.do_swaps(initial_cands);
  construct(s, use_alpha);
}

void Constructive::construct(solution& s, bool use_alpha /*= false*/) {
  assert(prm::batch_size >= 1);
  static vector<candidate> cands;
  cands.clear();

  if (s.num_assigned != pt::nland) {
    for (int c = 0; c < pt::nland; ++c)
      if (s.assigned[c] == -1) {
        for (int nb : pt::neighbours[c])
          if (s.assigned[nb] != -1) {
            cands.emplace_back(s.assigned[nb], c);
            break;
          }
      }
    assert(cands.size() > 0);
  }

  int iterations = 0, num_assigned = 0;
  const int bs = prm::batch_size;
  int cands_avg = 0, max_cand = nl<int>::min(), num_constructions = 0;
  while (cands.size()) {
    max_cand = max(max_cand, (int)cands.size());
    cands_avg += cands.size();
    ++num_constructions;
    ++iterations;
    if (stats::time_limit_exceeded()) exit(EXIT_SUCCESS);

    int last_mult = (use_alpha ? prm::mutation_greedy_alpha : 1.0);
    int last = max(0, int(cands.size() - bs * last_mult));
    assert(last >= 0 and last < (int)cands.size());
    if (last > 0) {
      nth_element(cands.begin(), cands.begin() + last, cands.end(),
                  [&](const candidate& c1, const candidate& c2) {
                    return compare_candidates(c2, c1, s);
                  });
#ifdef HARD_DEBUG
      for (const auto& c : cands)
        assert(not c.is_invalid());
#endif
      if (use_alpha) {
        int start = cands.size() - bs;
        for (int a = start; a >= last; --a) {
          int j = rng.rand_int(a, cands.size() - 1);
          if (j >= start) swap(cands[j], cands[a]);
        }
      }
    }

    static vector<candidate> accepted_cands;
    accepted_cands.clear();
    for (int i = max(0, (int)cands.size() - bs); i < (int)cands.size(); ++i) {
      if (not s.is_assigned(cands[i].cell)) accepted_cands.push_back(cands[i]);
    }

    cands.resize(max(0, (int)cands.size() - bs));
    num_assigned += s.do_swaps(accepted_cands);
    for (auto& c : accepted_cands) {
      if (c.is_invalid()) continue;
      assert(c.lot != -1);
      for (int nb : pt::neighbours[c.cell])
        if (not s.is_assigned(nb)) cands.emplace_back(c.lot, nb);
    }
  }
}

bool Constructive::compare_candidates(const candidate& c1, const candidate& c2,
                                      const solution& s) {
  assert(not c1.is_invalid() and not c2.is_invalid());

  if (s.is_assigned(c1.cell) and s.is_assigned(c2.cell))
    return c1.cell < c2.cell;
  else if (s.is_assigned(c2.cell))
    return false;
  else if (s.is_assigned(c1.cell))
    return true;

  int r = c1.river_value(s) - c2.river_value(s);
  if (r == 0) {
    int sr = c1.size_ratio(s) - c2.size_ratio(s);
    if (sr == 0) return c1.value(s) < c2.value(s);
    return sr < 0;
  }
  return r < 0;
}
