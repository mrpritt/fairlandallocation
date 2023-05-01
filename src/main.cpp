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
#include "ga.h"
#include "initial_positions.h"
#include "parameters.h"
#include "proterra.h"
#include "random.h"
#include "statistics.h"
#include "util.h"
#include "voronoi.h"
#include <boost/program_options.hpp>
#include <csignal>
#include <iostream>
#include <queue>

bool wrote_stats = false;

void write_stats(int) {
  stats::write_stats();
  wrote_stats = true;
  exit(EXIT_SUCCESS);
}

void write_stats_void() {
  if (not wrote_stats) write_stats(0);
}

solution naive() {
  static vi seeds;
  static vi assigned;
  seeds.resize(pt::lots);
  assigned.assign(pt::nland, -1);

  static queue<int> q;
  for (int i = 0; i < pt::lots; ++i) {
    while (true) {
      int r = rng.rand_int(0, pt::nland - 1);
      if (assigned[r] == -1) {
        seeds[i] = r;
        assigned[r] = i;
        q.push(r);
        break;
      }
    }
  }

  while (q.size()) {
    int c = q.front();
    q.pop();
    for (int nb : pt::neighbours[c]) {
      if (assigned[nb] == -1) {
        assigned[nb] = assigned[c];
        q.push(nb);
      }
    }
  }

  solution s;
  s.populate(assigned);
  return s;
}

int main(int argc, char** argv) {
  prm::parse_cmd_line(argc, argv);
  signal(SIGABRT, write_stats);
  signal(SIGTERM, write_stats);
  signal(SIGINT, write_stats);
  atexit(write_stats_void);

  try {
    pt::init();
    if (prm::solution_input.size()) {
      solution s;
      pt::read_solution(s, prm::solution_input);
      stats::add_stats(s);
    } else if (prm::naive) {
      int repl = 0;
      while (not stats::time_limit_exceeded() and repl < prm::max_generations) {
        solution s = naive();
        stats::add_stats(s);
        ++repl;
      }
    } else {
      ga().run();
    }
  } catch (std::exception& e) {
    pr("\nexception: {}\n", e.what());
    return EXIT_FAILURE;
  }
  write_stats(0);
  return EXIT_SUCCESS;
}
