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
#include "solution.h"
#include "constructive.h"
#include "lodepng/lodepng.h"
#include "parameters.h"
#include "random.h"
#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>

void solution::init() {
  area.assign(pt::lots, 0);
  num_river.assign(pt::lots, 0);
  assigned.assign(pt::nland, -1);
  of.init();
  rc.init();
  num_assigned = 0;
}

void solution::populate(vi a) {
  assert((int)a.size() == pt::nland);
  init();
  swap(a, assigned);
  for (int i = 0; i < pt::nland; ++i) {
    int lot = assigned[i];
    if (lot != -1) {
      ++num_assigned;
      ++area[lot];
      if (pt::nx_river[i]) ++num_river[lot];
    }
  }
  rc.populate(*this);
  of.populate(*this);
}

int solution::do_swaps(vector<candidate>& cands) {
  int done = 0;
  for (auto& c : cands) {
    if (is_assigned(c.cell)) {
      c.set_invalid();
      continue;
    }
    assert(assigned[c.cell] == -1);
    assigned[c.cell] = c.lot;
    ++num_assigned;
    ++area[c.lot];
    if (pt::nx_river[c.cell]) ++num_river[c.lot];
    ++done;
  }
  rc.do_swaps(cands, *this);
  of.do_swaps(cands);
  return done;
}

bool solution::check_border_brute_force(int c) const {
  for (auto nb : pt::neighbours[c])
    if (assigned[nb] != assigned[c]) return true;
  return false;
}

ll solution::value() const {
#ifdef HARD_DEBUG
  of.assert_value_acceptable();
#endif
  return of.value;
}

int solution::river_value() const {
#ifdef HARD_DEBUG
  if (rc.value != rc.get_value_brute_force(*this))
    pr("rc.value: {}, rc.getValueBruteForce(): {}\n", rc.value,
       rc.get_value_brute_force(*this));
  assert(rc.value == rc.get_value_brute_force(*this));
#endif
  return rc.value;
}

double solution::size_ratio() const {
  int sa = area[small_lot()];
  if (sa == 0) return (1 << 28);
  double sr = double(area[big_lot()]) / double(sa);
  return sr > prm::maximum_size_ratio ? sr : 0;
}

void solution::write_to_png(const string& filename) const {
  static vbyte rcolor(pt::lots), gcolor(pt::lots), bcolor(pt::lots);
  static bool colors_chosen = false;
  const double golden_ratio_conjugate = 0.618033988749895;
  if (not colors_chosen) {
    double h = rng.rand_double(0.0, 1.0), intpart;
    for (int i = 0; i < pt::lots; ++i) {
      h += golden_ratio_conjugate;
      h = modf(h, &intpart);
      hsv_to_rgb(h * 255, 125, 240, rcolor[i], gcolor[i], bcolor[i]);
    }
    colors_chosen = true;
  }

  vbyte data(pt::r_size * pt::c_size * 4);
  auto set_color = [&](int r, int c, int red, int green, int blue) {
    data[4 * (r * pt::c_size + c) + 0] = red;
    data[4 * (r * pt::c_size + c) + 1] = green;
    data[4 * (r * pt::c_size + c) + 2] = blue;
    data[4 * (r * pt::c_size + c) + 3] = 255;
  };

  for (int r = 0; r < pt::r_size; ++r)
    for (int c = 0; c < pt::c_size; ++c) {
      if (pt::cell_type[r][c] == pt::river) {
        set_color(r, c, 0, 0, 255);
      } else if (pt::in_cell[r * pt::c_size + c] == -2) {
        set_color(r, c, 0, 0, 0);
      } else if (pt::cell_type[r][c] == pt::preserve) {
        set_color(r, c, 0, 0, 0);
      } else {
        int i = pt::index_from_rc[r][c];
        if (not is_assigned(i)) {
          set_color(r, c, 255, 255, 255);
        } else {
          set_color(r, c, rcolor[assigned[i]], gcolor[assigned[i]],
                    bcolor[assigned[i]]);
        }
      }
    }
  vbyte png;
  uint error = lodepng::encode(png, data, pt::c_size, pt::r_size);
  if (error) {
    fmt::print("lodepng error: %s\n", lodepng_error_text(error));
    exit(EXIT_FAILURE);
  }
  lodepng::save_file(png, filename);
}
