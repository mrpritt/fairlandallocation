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
#include "rivers_constraint.h"
#include "constructive.h"
#include "proterra.h"
#include "solution.h"
#include "util.h"
#include <cassert>
#include <iostream>
#include <memory>

void rivers_constraint::init() {
  sa.resize(pt::lots);
  index_sa.resize(pt::lots);
  for (int i = 0; i < pt::lots; ++i) {
    sa[i] = index_sa[i] = i;
  }
  x = violations = value = 0;
}

void rivers_constraint::compute_value_brute_force(int& nx, int& vio, int& cost,
                                                  const solution& s) const {
  nx = -1;
  cost = vio = 0;
  for (int i = 0; i < pt::lots; ++i)
    if (s.num_river[sa[i]] == 0) {
      nx = sa[i];
      break;
    }
  if (nx == -1) return;
  assert(s.num_river[nx] == 0);
  for (int i = 0; i < pt::lots; ++i) {
    assert(sa[index_sa[i]] == i);
    if (s.num_river[i] > 0 and s.area[i] > s.area[nx]) {
      cost += s.area[i] - s.area[nx];
      ++vio;
    }
  }
  assert(cost >= 0);
}

void rivers_constraint::calc_swap(int from, int to, int c, candidate& sp,
                                  solution& s) {
  sp.s1 = -1, sp.s2 = -1, sp.s3 = -1, sp.s4 = -1;
  int cost = value, vio = violations, nx = x;
  auto &area = s.area, &num_river = s.num_river;
  assert(value != -1);

  if (from != -1) {
    assert(s.assigned[c] == from);
    int i = index_sa[from], j;

    for (j = i - 1; j >= 0; --j) {
      if (area[sa[j]] < area[from]) break;
      if (from == x and num_river[sa[j]] > 0) ++vio;
    }

    if (from == x) {
      for (int j = i + 1; j < pt::lots; ++j) {
        if (area[sa[j]] > area[from]) break;
        if (num_river[sa[j]] > 0) {
          ++vio;
        }
      }
    }

#ifdef HARD_DEBUG
    if (j >= 0) assert(area[sa[j]] < area[from]);
    if (j != pt::lots - 1) assert(area[sa[j + 1]] == area[from]);
#endif

    sp.s1 = i, sp.s2 = j + 1;
    assert(sa[sp.s1] == from);
    std::swap(sa[sp.s1], sa[sp.s2]);
    index_sa[sa[sp.s1]] = sp.s1;
    index_sa[sa[sp.s2]] = sp.s2;
    if (pt::nx_river[c]) {
      assert(num_river[from] > 0);
      --num_river[from];
    }
    --area[from];

    if (num_river[from] == 0) {
      assert(from != x or (from == nx and nx == x));
      if (nx == -1 or area[from] < area[nx]) {
        nx = from;
      }
      if (nx != x) {
        assert(nx == from);
        assert(x == -1 or index_sa[nx] < index_sa[x]);
        int ax;
        if (x != -1) {
          ax = area[x];
          if (area[nx] + 1 > ax) {
            --vio;
            cost -= (area[nx] + 1 - area[x]);
          }
          cost += (ax - area[nx]) * vio;
        } else
          ax = nl<int>::max();
        for (j = index_sa[nx] + 1; j < pt::lots; ++j) {
          if (area[sa[j]] > ax) break;
          if (num_river[sa[j]] > 0 and area[sa[j]] > area[nx]) {
            cost += area[sa[j]] - area[nx];
            ++vio;
          }
        }
      } else if (pt::nx_river[c] and x != -1) {
        cost -= (area[from] + 1 - area[x]);
        --vio;
      } else if (from == x) {
        cost += vio;
      }
    } else if (nx != -1) {
      if (area[from] >= area[nx]) {
        --cost;
        if (area[from] == area[nx]) --vio;
      }
    }
  }

  if (to != -1) {
    assert(s.assigned[c] != to);
    int i = index_sa[to], j;

    if (to == nx) {
      for (j = i - 1; j >= 0 and area[sa[j]] == area[to]; --j) {
        if (num_river[sa[j]] == 0) {
          nx = sa[j];
          break;
        }
      }
    }

    for (j = i + 1; j < pt::lots; ++j) {
      if (area[sa[j]] > area[to]) break;
      if (to == nx) {
        if (num_river[sa[j]] == 0) nx = sa[j];
      }
    }

#ifdef HARD_DEBUG
    if (j != pt::lots) assert(area[sa[j]] > area[to]);
    if (j != 0) assert(area[sa[j - 1]] == area[to]);
#endif
    sp.s3 = i, sp.s4 = j - 1;
    assert(sa[sp.s3] == to);
    assert(index_sa[sa[j - 1]] == j - 1);
    std::swap(sa[sp.s3], sa[sp.s4]);
    index_sa[sa[sp.s3]] = sp.s3;
    index_sa[sa[sp.s4]] = sp.s4;
    if (pt::nx_river[c]) ++num_river[to];
    ++area[to];

    if (nx != -1 and nx != to and area[nx] < area[to] and num_river[to] > 0) {
      if (num_river[to] == 1 and pt::nx_river[c]) {
        cost += area[to] - area[nx];
        ++vio;
      } else {
        ++cost;
        if (area[to] == area[nx] + 1) ++vio;
      }
    } else if (nx == to) {
      if (num_river[to] > 0) {
        nx = -1;
        for (j = index_sa[to] + 1; j < pt::lots; ++j) {
          if (nx != -1 and area[sa[j]] != area[nx]) {
            assert(area[sa[j]] > area[nx]);
            break;
          }
          if (num_river[sa[j]] == 0) {
            if (nx == -1) nx = sa[j];
          } else {
            --vio;
            cost -= area[sa[j]] - (area[to] - 1);
          }
        }
        if (nx == -1) {
          cost = vio = 0;
        } else {
          cost -= vio * (area[nx] - (area[to] - 1));
        }
      } else {
        assert(not pt::nx_river[c]);
        cost -= vio;
        for (j = index_sa[nx] + 1; j < pt::lots; ++j) {
          if (area[sa[j]] != area[nx]) break;
          if (num_river[sa[j]] > 0) --vio;
        }
      }
    }
  }

  if (nx == -1) {
    assert(cost == 0);
    assert(vio == 0);
  }

#ifdef HARD_DEBUG
  int nx2, cost2, vio2;
  compute_value_brute_force(nx2, vio2, cost2, s);
  assert(cost == cost2);
  assert(vio == vio2);
#endif

  if (to != -1) {
    assert(sp.s3 != -1 and sp.s4 != -1);
    --area[to];
    if (pt::nx_river[c]) --num_river[to];
    swap(sa[sp.s3], sa[sp.s4]);
    index_sa[sa[sp.s3]] = sp.s3;
    index_sa[sa[sp.s4]] = sp.s4;
  }
  if (from != -1) {
    assert(sp.s1 != -1 and sp.s2 != -1);
    ++area[from];
    if (pt::nx_river[c]) ++num_river[from];
    swap(sa[sp.s1], sa[sp.s2]);
    index_sa[sa[sp.s1]] = sp.s1;
    index_sa[sa[sp.s2]] = sp.s2;
  }
  sp.value = cost;
  sp.nx = nx;
  sp.vio = vio;
}

void rivers_constraint::do_swaps(vector<::candidate>& cands,
                                 const solution& s) {
  if ((int)cands.size() > 0) populate(s);
}

void rivers_constraint::populate(const solution& s) {
  static vii areas;
  areas.resize(pt::lots);

  for (int i = 0; i < pt::lots; ++i) {
    areas[i].ff = s.area[i];
    areas[i].ss = i;
  }

  sort(areas.begin(), areas.end());
  sa.resize(pt::lots);
  index_sa.resize(pt::lots);

  for (int i = 0; i < pt::lots; ++i) {
    sa[i] = areas[i].ss;
    index_sa[sa[i]] = i;
  }
  compute_value_brute_force(s);
}
