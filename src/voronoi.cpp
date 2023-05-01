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
#include "voronoi.h"
#include "proterra.h"
#include "solution.h"
#include "statistics.h"
#include "util.h"
#include <algorithm>
#include <boost/heap/fibonacci_heap.hpp>
#include <cassert>
#include <ciso646>
#include <iostream>

vi voronoi::construct(const vi& pos) {
  vi assigned(pt::nland, -1);
  vi dist(pt::nland, -1);
  vector<boost::heap::fibonacci_heap<iii>::handle_type> handle(pt::nland);
  boost::heap::fibonacci_heap<iii> pq;

  assert((int)pos.size() == pt::lots);
  for (int i = 0; i < pt::lots; ++i) {
    int c = pos[i];
    dist[c] = pt::val[c];
    handle[c] = pq.push(iii(-dist[c], c, i));
  }

  int d, c, p, num_expansions = 0;
  while (pq.size()) {
    tie(d, c, p) = pq.top();
    pq.pop();
    d = -d;
    if (dist[c] != d or assigned[c] != -1) continue;
    assigned[c] = p;
    ++num_expansions;

    for (auto nb : pt::neighbours[c]) {
      int cd = d + pt::val[nb];
      if (dist[nb] == -1) {
        handle[nb] = pq.push(iii(-cd, nb, p));
        dist[nb] = cd;
      } else if (dist[nb] > cd) {
        pq.update(handle[nb], iii(-cd, nb, p));
        dist[nb] = cd;
      }
    }
  }
  return move(assigned);
}
