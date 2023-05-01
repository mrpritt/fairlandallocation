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
#include "convexity.h"
#include "defines.h"
#include "geometry.h"
#include "proterra.h"
#include "solution.h"
#include "util.h"
#include <cassert>
#include <vector>

double recalc_ch_area(int lot, const solution& s) {
  vector<point> pts;
  for (int i = 0; i < pt::nland; ++i) {
    if (s.assigned[i] == lot) {
      int r = pt::rc_from_index[i].ff, c = pt::rc_from_index[i].ss;
      pts.emplace_back(r - 1, c - 1);
      pts.emplace_back(r, c - 1);
      pts.emplace_back(r - 1, c - 1);
      pts.emplace_back(r, c);
    }
  }
  if (pts.size() == 0) return 0;
  pts = convex_hull(pts);
  return polygon_area(pts);
}

double avg_ch_ratio(const solution& s) {
  vd ch_values;
  for (int i = 0; i < pt::lots; ++i) {
    double ch = recalc_ch_area(i, s);
    if (ch and s.area[i]) ch_values.push_back(ch / double(s.area[i]));
  }
  return ch_values.size() ? avg(ch_values) : 0;
}
