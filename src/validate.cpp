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
#include "validate.h"
#include "proterra.h"
#include "random.h"
#include "util.h"
#include <algorithm>
#include <cassert>
#include <queue>

void validate_solution(const solution& s) {
  (void)s;
#ifdef HARD_DEBUG
  static vi any_cell, area, vis, river, border_size;
  static vvi borders;
  int numAssigned = 0;
  any_cell.assign(pt::lots, -1);
  area.assign(pt::lots, 0);
  river.assign(pt::lots, 0);
  borders.assign(pt::lots, vi());
  vis.assign(pt::nland, 0);

  for (int i = 0; i < pt::nland; ++i) {
    if (s.assigned[i] != -1) {
      ++numAssigned;
      ++area[s.assigned[i]];
      if (rng.rand_double(0.0, 1.0) < 1.0 / double(area[s.assigned[i]]))
        any_cell[s.assigned[i]] = i;
      if (pt::nx_river[i]) ++river[s.assigned[i]];
      if (s.check_border_brute_force(i)) {
        borders[s.assigned[i]].push_back(i);
      }
    }
  }

  int vis_total = 0;
  for (int i = 0; i < pt::lots; ++i) {
    assert(area[i] > 0);
    assert(area[i] == s.area[i]);
    assert(any_cell[i] != -1);
    assert(s.assigned[any_cell[i]] == i);
    queue<int> q;
    q.push(any_cell[i]);
    vis[any_cell[i]] = 1;
    int visLot = 0;
    while (q.size()) {
      int c = q.front();
      q.pop();
      ++visLot;
      for (int nb : pt::neighbours[c])
        if (s.assigned[nb] == s.assigned[c] and vis[nb] == 0) {
          assert(s.assigned[nb] == i);
          vis[nb] = 1;
          q.push(nb);
        }
    }
    if (visLot != area[i]) {
      s.write_to_png("error.png");
      assert(visLot == area[i]);
    }
    vis_total += visLot;
  }
  assert(pt::nland == numAssigned);
  assert(vis_total == numAssigned);
  assert(s.num_assigned == numAssigned);
  s.of.assert_value_acceptable();
  assert(s.rc.value == s.rc.get_value_brute_force(s));
#endif
}
