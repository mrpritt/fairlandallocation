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
#include "initial_positions.h"
#include "objective_function.h"
#include "proterra.h"
#include "random.h"
#include "solution.h"
#include "util.h"
#include "voronoi.h"
#include <cassert>
#include <iostream>
#include <numeric>
#include <queue>
#include <stdexcept>
#include <tuple>

vvi initial_positions::cc;
vi initial_positions::cc_num_lots;

void save_initial_positions_as_png(const vi& v, string filename) {
  solution s;
  s.assigned.assign(pt::nland, -1);
  for (int i = 0; i < (int)v.size(); ++i) {
    int u = v[i], r, c;
    tie(r, c) = pt::rc_from_index[u];
    s.assigned[u] = i;
    int dr[] = {0, 0, 1, -1, 1, 1, -1, -1}, dc[] = {1, -1, 0, 0, 1, -1, 1, -1};
    for (int j = 0; j < 8; ++j) {
      int rr = r + dr[j], cc = c + dc[j];
      if (rr >= 0 and rr < pt::r_size and cc >= 0 and cc < pt::c_size) {
        int v = pt::index_from_rc[rr][cc];
        if (v != -1) s.assigned[v] = i;
      }
    }
  }
  s.write_to_png(filename);
}

vi initial_positions::generate_initial_positions() {
  vi v = gen_rand_initial_positions();
  v = k_means_step(v);
  return move(v);
}

vi initial_positions::k_means_step(const vi& v) {
  const int k = pt::lots;
  vvi parcel_cells(k);
  vi assigned = voronoi().construct(v);
  assert((int)assigned.size() == pt::nland);
  for (int i = 0; i < pt::nland; ++i)
    if (assigned[i] != -1) parcel_cells[assigned[i]].push_back(i);
  vi ans;
  vvb centroid(pt::r_size, vb(pt::c_size, false));
  for (int i = 0; i < k; ++i) {
    double r = 0, c = 0;
    for (auto j : parcel_cells[i]) {
      r += pt::rc_from_index[j].ff;
      c += pt::rc_from_index[j].ss;
    }
    r = round(r / (double)parcel_cells[i].size());
    c = round(c / (double)parcel_cells[i].size());
    int ri = r, ci = c;
    if (pt::cell_type[ri][ci] != pt::land or centroid[ri][ci]) {
      queue<ii> q;
      vvb visited(pt::r_size, vb(pt::c_size, false));
      visited[ri][ci] = true;
      q.push(ii(ri, ci));
      while (q.size()) {
        ii v = q.front();
        q.pop();
        if (pt::cell_type[v.ff][v.ss] == pt::land and
            not centroid[v.ff][v.ss]) {
          ri = v.ff;
          ci = v.ss;
          break;
        }
        for (int i = 0; i < 4; ++i) {
          int nr = v.ff + dr[i], nc = v.ss + dc[i];
          if (nr < 0 or nr >= pt::r_size or nc < 0 or nc >= pt::c_size)
            continue;
          if (not visited[nr][nc]) {
            visited[nr][nc] = true;
            q.push(ii(nr, nc));
          }
        }
      }
    }
    centroid[ri][ci] = true;
    ans.push_back(pt::index_from_rc[ri][ci]);
  }
  return move(ans);
}

vi initial_positions::gen_rand_initial_positions() {
  vi v, u;
  for (int i = 0; i < (int)cc.size(); ++i) {
    if (cc_num_lots[i] == 0) continue;
    choice(cc[i], u, cc_num_lots[i]);
    assert((int)u.size() == cc_num_lots[i]);
    for (int i : u)
      v.push_back(i);
  }
  assert((int)v.size() == pt::lots);
  return move(v);
}
