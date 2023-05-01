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
#include "proterra.h"
#include "initial_positions.h"
#include "parameters.h"
#include "random.h"
#include "solution.h"
#include "util.h"
#include <cassert>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <queue>
#include <tuple>

random_number_generator rng;
int pt::river_pct, pt::num_apt_classes, pt::lots, pt::r_size, pt::c_size,
    pt::nland, pt::nriver;
vi pt::in_cell, pt::val;
vb pt::nx_river, pt::global_border;
vvi pt::neighbours, pt::index_from_rc;
vii pt::rc_from_index;
std::vector<std::vector<pt::cell_type_enum>> pt::cell_type;

void pt::init() {
  ifstream f(prm::input_filename);
  if (f.fail())
    throw runtime_error(
        fmt::format("could not open input file {}", prm::input_filename));

  int num_parcels = -1;
  f >> r_size >> c_size >> num_parcels >> num_apt_classes >> river_pct;
  if (num_parcels > 0 and lots <= 0) lots = num_parcels;

  index_from_rc.assign(r_size, vi(c_size, -1));
  cell_type.assign(r_size, vector<cell_type_enum>(c_size, land));
  nland = nriver = 0;
  for (int r = 0; r < r_size; ++r) {
    for (int c = 0; c < c_size; ++c) {
      int x;
      if (not(bool)(f >> x))
        throw runtime_error("Error: could not read cell input from instance.");
      in_cell.push_back(x);
      if (x == -1) {
        cell_type[r][c] = river;
        ++nriver;
      } else if (x == 0 or x == -2) {
        cell_type[r][c] = preserve;
      } else {
        cell_type[r][c] = land;
        ++nland;
      }
    }
  }
  vvii cc;
  vvb visited(r_size, vb(c_size, false));
  for (int r = 0; r < r_size; ++r)
    for (int c = 0; c < c_size; ++c) {
      if (cell_type[r][c] == land and not visited[r][c]) {
        cc.push_back(vii());
        visited[r][c] = true;
        queue<ii> q;
        q.push(mp(r, c));
        while (q.size()) {
          ii p = q.front();
          q.pop();
          cc.back().push_back(p);
          for (int i = 0; i < prm::neighborhood_size; ++i) {
            int nr = p.ff + dr[i], nc = p.ss + dc[i];
            if (inside_boundaries(nr, nc) and cell_type[nr][nc] == land and
                not visited[nr][nc]) {
              visited[nr][nc] = true;
              q.push(mp(nr, nc));
            }
          }
        }
      }
    }

  vi cc_num_lots(cc.size());
  int sum = 0;
  pr("found {} connected components\n", cc.size());
  for (int i = 0; i < (int)cc.size(); ++i) {
    cc_num_lots[i] = int((double)cc[i].size() * lots / (double)nland);
    pr("{} ", cc_num_lots[i]);
    sum += cc_num_lots[i];
    if (cc_num_lots[i] == 0)
      for (auto& p : cc[i]) {
        in_cell[p.ff * c_size + p.ss] = -2;
        cell_type[p.ff][p.ss] = preserve;
      }
  }
  pr("\n");
  pr("proterra::lots: {}, sum: {}\n", pt::lots, sum);
  for (auto& i : cc_num_lots)
    if (i > 1) {
      i += (pt::lots - sum);
      break;
    }
  for (auto& i : cc_num_lots)
    pr("{} ", i);
  pr("\n");

  assert((int)accumulate(cc_num_lots.begin(), cc_num_lots.end(), 0) == lots);
  initial_positions::cc_num_lots = move(cc_num_lots);
  nland = 0;
  int num_preserve = 0;
  for (int i = 0; i < (int)in_cell.size(); ++i) {
    int r = i / c_size, c = i % c_size;
    if (cell_type[r][c] == land) {
      assert(in_cell[i] > 0);
      bool is_next_to_river = false;
      bool is_border = false;
      for (int j = 0; j < prm::neighborhood_size; ++j) {
        int nr = r + dr[j], nc = c + dc[j];
        if (inside_boundaries(nr, nc)) {
          if (cell_type[nr][nc] == river) {
            assert(in_cell[nr * c_size + nc] == -1);
            is_next_to_river = true;
          } else if (cell_type[nr][nc] == preserve) {
            assert(in_cell[nr * c_size + nc] == 0 or
                   in_cell[nr * c_size + nc] == -2);
            is_border = true;
          }
        } else {
          is_border = true;
        }
      }
      index_from_rc[r][c] = nland;
      rc_from_index.emplace_back(r, c);
      val.push_back(in_cell[i]);
      nx_river.push_back(is_next_to_river);
      global_border.push_back(is_border);
      ++nland;
    } else if (cell_type[r][c] == preserve) {
      ++num_preserve;
    }
  }

  assert(nland == (int)val.size());
  neighbours.resize(nland);
  for (int i = 0; i < nland; ++i) {
    int r, c;
    tie(r, c) = rc_from_index[i];
    neighbours[i] = vi();
    for (int j = 0; j < prm::neighborhood_size; ++j) {
      int nr = r + dr[j], nc = c + dc[j];
      if (inside_boundaries(nr, nc) and cell_type[nr][nc] == land) {
        assert(index_from_rc[nr][nc] != -1);
        neighbours[i].push_back(index_from_rc[nr][nc]);
      }
    }
  }
  initial_positions::cc.assign(cc.size(), vi());
  for (int i = 0; i < (int)cc.size(); ++i) {
    if (initial_positions::cc_num_lots[i] == 0) continue;
    for (auto& p : cc[i]) {
      assert(cell_type[p.ff][p.ss] == land);
      initial_positions::cc[i].push_back(index_from_rc[p.ff][p.ss]);
    }
  }
}

void pt::read_solution(solution& s, const string& solution_input_filename) {
  ifstream f(solution_input_filename);
  if (not f or f.fail()) {
    throw runtime_error(
        fmt::format("could not open input file {}", solution_input_filename));
  }
  vi assigned(pt::nland, -1);
  map<int, int> ss;
  int lot_num = 0;
  for (int i = 0; i < r_size; ++i)
    for (int j = 0; j < c_size; ++j) {
      int x;
      f >> x;
      int index = index_from_rc[i][j];
      if (x == -1 or x == 0 or index < 0 or index >= pt::nland) continue;
      if (ss.find(x) == ss.end()) ss[x] = lot_num++;
      assigned[index] = ss[x];
    }
  s.init();
  s.populate(assigned);
}
