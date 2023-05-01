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
#include "ga.h"
#include "constructive.h"
#include "initial_positions.h"
#include "matching.h"
#include "parameters.h"
#include "proterra.h"
#include "random.h"
#include "statistics.h"
#include "util.h"
#include "validate.h"
#include <algorithm>
#include <cassert>
#include <queue>
#include <tuple>

bool compare_solutions(const uptr<solution>& a, const uptr<solution>& b) {
  assert(a != nullptr and b != nullptr);
  return *a < *b;
}

void ga::run() {
  pr("running genetic algorithm\n");
  const int pop_size = prm::pop_size;
  const double rates_sum =
      (prm::crossover_ratio + prm::new_rate + prm::keep_ratio);
  const int crossover_size = pop_size * (prm::crossover_ratio / rates_sum);
  const int new_size = pop_size * (prm::new_rate / rates_sum);
  const int keep_size = pop_size - crossover_size - new_size;
  pr("--population size {}\n", pop_size);
  pr("--crossover_size {}\n", crossover_size);
  pr("--keep_size {}\n", keep_size);
  pr("--new_size {}\n", new_size);
  pop.resize(pop_size);
  pop2.resize(pop_size);

  pr("\ngenerating initial population...\n");
  for (int i = 0; i < pop_size; ++i) {
    pr("{}{}", i, (i == pop_size - 1 ? "" : ", "));
    pop[i].reset(new solution());
    pop[i]->init();
    auto pos = initial_positions::generate_initial_positions();
    ++stats::num_new_solutions;
    Constructive::construct_from_seeds(*pop[i], pos);
    if (stats::global_best.num_assigned == 0 or *pop[i] < stats::global_best)
      stats::global_best = *pop[i];
    validate_solution(*pop[i]);
  }

  pr("\nstarting generations...\n");
  int best_since = 0;
  auto best_stats =
      make_tuple(numeric_limits<int>::max(), numeric_limits<double>::max(),
                 numeric_limits<double>::max());
  while (true) {
    int best = 0;
    for (int i = 0; i < pop_size; ++i) {
      if (compare_solutions(pop[i], pop[best])) {
        best = i;
        auto stats = make_tuple(pop[best]->river_value(),
                                pop[best]->size_ratio(), pop[best]->value());
        if (stats < best_stats) {
          best_since = stats::num_generations;
          best_stats = stats;
          pr("best_since = {}\n", best_since);
        }
      }
      if (pop2[i] == nullptr) pop2[i].reset(new solution());
    }

    pr("\nGENERATION #{}:\n", stats::num_generations);
    stats::add_stats(*pop[best]);

    if (stats::time_limit_exceeded() or
        stats::num_generations >= prm::max_generations)
      exit(EXIT_SUCCESS);
    ++stats::num_generations;

    bool restart = prm::restart == -1
                       ? false
                       : (stats::num_generations - best_since) >= prm::restart;
    if (restart) {
      pr("restarting...\n");
      best_since = stats::num_generations;
      pr("best_since = {}\n", best_since);
      best_stats =
          make_tuple(numeric_limits<int>::max(), numeric_limits<double>::max(),
                     numeric_limits<double>::max());
    } else {
      for (int i = 0; i < crossover_size; ++i) {
        int p1, p2;
        select_parents_by_tournament(&p1, &p2);
        if (prm::do_crossover) {
          ++stats::num_crossovers;
          crossover(*pop[p1], *pop[p2], *pop2[i]);
          validate_solution(*pop2[i]);
        } else {
          *pop2[i] = compare_solutions(pop[p1], pop[p2]) ? *pop[p1] : *pop[p2];
        }

        if (prm::do_mutation) {
          ++stats::num_mutations;
          mutation(*pop2[i]);
          validate_solution(*pop2[i]);
        }
      }
    }
    for (int i = restart ? 0 : crossover_size;
         i < (restart ? pop_size : crossover_size + new_size); ++i) {
      pop2[i]->init();
      ++stats::num_new_solutions;
      Constructive::construct_from_seeds(
          *pop2[i], initial_positions::generate_initial_positions());
      validate_solution(*pop2[i]);
      if (restart) {
        auto stats = make_tuple(pop2[i]->river_value(), pop2[i]->size_ratio(),
                                pop2[i]->value());
        if (stats < best_stats) {
          best_stats = stats;
        }
      }
    }
    if (not restart) {
      nth_element(pop.begin(), pop.begin() + keep_size, pop.end(),
                  compare_solutions);
      for (int i = crossover_size + new_size; i < pop_size; ++i) {
        swap(pop2[i], pop[i - crossover_size - new_size]);
      }
    }
    swap(pop2, pop);
  }
}

void ga::crossover(const solution& p1, const solution& p2, solution& child) {
  static vvi cost;
  static vi lmate, rmate;
  cost.resize(pt::lots);
  lmate.resize(pt::lots);
  rmate.resize(pt::lots);
  for (auto& c : cost)
    c.assign(pt::lots, 0);

  for (int c = 0; c < pt::nland; ++c) {
    int lot1 = p1.assigned[c], lot2 = p2.assigned[c];
    if (lot1 < 0 or lot1 >= pt::lots or lot2 < 0 or lot2 >= pt::lots) {
      pr("lp1: {}, lp2: {}\n", lot1, lot2);
      p1.write_to_png("error_p1.png");
      p2.write_to_png("error_p2.png");
    }
    assert(lot1 >= 0 and lot1 < pt::lots);
    assert(lot2 >= 0 and lot2 < pt::lots);
    --cost[lot1][lot2];
  }

  min_cost_bipartite_matching(cost, lmate, rmate);

  static vi assigned, cc_num, cc_largest, cc_start;
  assigned.assign(pt::nland, -1);
  cc_num.assign(pt::lots, 0);
  cc_largest.assign(pt::lots, -1);
  cc_start.assign(pt::lots, -1);

  auto bfs = [&](int start, int lot_assign, int lotp1, int lotp2) {
    static queue<int> q;
    assert(q.empty());
    q.push(start);
    assigned[start] = lot_assign;
    int size = 1;
    while (q.size()) {
      int c = q.front();
      q.pop();
      for (int nb : pt::neighbours[c])
        if (p1.assigned[nb] == lotp1 and p2.assigned[nb] == lotp2 and
            assigned[nb] == -1) {
          q.push(nb);
          assigned[nb] = lot_assign;
          ++size;
        }
    }
    return size;
  };

  for (int c = 0; c < pt::nland; ++c) {
    int lot1 = p1.assigned[c], lot2 = p2.assigned[c];
    assert(lot1 != -1 and lot2 != -1);
    if (assigned[c] == -1 and lmate[lot1] == lot2) {
      assert(rmate[lot2] == lot1);
      ++cc_num[lot1];
      int size = bfs(c, lot1, lot1, lot2);
      if (size > cc_largest[lot1]) {
        cc_largest[lot1] = size;
        cc_start[lot1] = c;
      }
    }
  }

  int lot_num = 0;
  assigned.assign(pt::nland, -1);

  for (int l = 0; l < pt::lots; ++l) {
    if (cc_num[l] > 1) {
      assert(cc_largest[l] > 0);
      ++stats::num_disconnected_crossover_lots;
    }
    if (cc_largest[l] > 0) {
      assert(cc_num[l] > 0);
      int st = cc_start[l], lot1 = p1.assigned[st], lot2 = p2.assigned[st];
      assert(lmate[lot1] == lot2 and rmate[lot2] == lot1);
      bfs(st, lot_num, lot1, lot2);
      ++lot_num;
    }
  }

  if (lot_num != pt::lots) {
    ++stats::num_empty_crossover_lots;
    assert(lot_num < pt::lots);
    vi empty_lots(pt::lots - lot_num);
    iota(empty_lots.begin(), empty_lots.end(), lot_num);
    empty_lots_fix(empty_lots, assigned);
  }
  child.init();
  swap(assigned, child.assigned);
  child.populate(move(child.assigned));
  Constructive::construct(child);
}

void ga::mutation(solution& s) {
  static vi dist, empty_lots, any, any_seen, assigned_new;
  vi& assigned = s.assigned;
  vi& area = s.area;
  dist.assign(pt::nland, -1);
  empty_lots.clear();
  any.assign(pt::lots, -1);
  any_seen.assign(pt::lots, 0);
  assigned_new.assign(pt::nland, -1);

  static queue<int> q;
  assert(q.empty());
  for (int c = 0; c < pt::nland; ++c) {
    assert(assigned[c] != -1);
    any[assigned[c]] = c;
    if (s.check_border_brute_force(c)) {
      dist[c] = 0;
      q.push(c);
    }
  }

  while (q.size()) {
    int c = q.front();
    q.pop();
    assert(assigned[c] != -1);
    assert(area[assigned[c]] > 0);

    if (area[assigned[c]] == 1) {
      ++any_seen[assigned[c]];
      any[assigned[c]] = c;
      continue;
    }

    assigned[c] = -1;
    if (dist[c] + 1 < prm::mutation_brush_size) {
      for (int nb : pt::neighbours[c])
        if (dist[nb] == -1) {
          dist[nb] = dist[c] + 1;
          q.push(nb);
        }
    } else if (dist[c] + 1 == prm::mutation_brush_size) {
      for (int nb : pt::neighbours[c])
        if (assigned[nb] != -1 and dist[nb] == -1) {
          ++any_seen[assigned[nb]];
          if (rng.rand_double(0.0, 1.0) < 1.0 / double(any_seen[assigned[nb]]))
            any[assigned[nb]] = nb;
        }
    }
  }

  assert(q.empty());
  dist.assign(pt::nland, -1);

  for (int i = 0; i < pt::lots; ++i) {
    assert((any[i] != -1 and any_seen[i] >= 0) or
           (any[i] == -1 and any_seen[i] == 0));
    if (any[i] != -1) {
      q.push(any[i]);
      dist[any[i]] = i;
    }
  }

  while (q.size()) {
    int c = q.front();
    q.pop();
    assigned_new[c] = dist[c];
    for (int nb : pt::neighbours[c])
      if (assigned[nb] != -1 and dist[nb] == -1) {
        dist[nb] = dist[c];
        q.push(nb);
      }
  }

  swap(assigned, assigned_new);

  s.populate(move(assigned));
  Constructive().construct(s, true);
}

void ga::select_parents_by_tournament(int* p1, int* p2) {
  const int sz = prm::tournament_size;
  if (sz == 3) {
    int a = rng.rand_int(0, prm::pop_size - 1), b, c;
    do {
      b = rng.rand_int(0, prm::pop_size - 1);
    } while (b == a);
    do {
      c = rng.rand_int(0, prm::pop_size - 1);
    } while (c == a or c == b);
    bool bab = compare_solutions(pop[a], pop[b]),
         bac = compare_solutions(pop[a], pop[c]),
         bbc = compare_solutions(pop[b], pop[c]);
    if (bab) {
      *p1 = a;
      *p2 = bbc ? b : c;
    } else {
      *p1 = b;
      *p2 = bac ? a : c;
    }
  } else {
    assert(sz >= 4);
    static vi nums, chosen;
    if (nums.empty()) {
      nums.resize(prm::pop_size);
      chosen.resize(sz);
      iota(nums.begin(), nums.end(), 0);
    }
    choice(nums, chosen, sz);
    *p1 = chosen[0], *p2 = chosen[1];
    if (compare_solutions(pop[*p2], pop[*p1])) swap(*p1, *p2);
    for (int i = 2; i < sz; ++i) {
      int ci = chosen[i];
      if (compare_solutions(pop[ci], pop[*p2])) {
        if (compare_solutions(pop[ci], pop[*p1])) {
          *p2 = *p1;
          *p1 = ci;
        } else {
          *p2 = ci;
        }
      }
    }
  }
}

void ga::empty_lots_fix(const vi& empty_lots, vi& assigned) {
  if (empty_lots.empty()) return;
  static vi chosen;
  chosen.resize(empty_lots.size());
  int j = 0;
  for (int i = 0; i < pt::nland and j < (int)empty_lots.size(); ++i)
    if (assigned[i] == -1) chosen[j++] = i;
  assert(j == (int)empty_lots.size());
  for (int i = chosen.back() + 1; i < pt::nland; ++i) {
    if (assigned[i] == -1) {
      ++j;
      int k = rng.rand_int(0, j);
      if (k < (int)empty_lots.size()) chosen[k] = i;
    }
  }
  for (int i = 0; i < (int)chosen.size(); ++i)
    assigned[chosen[i]] = empty_lots[i];
}
