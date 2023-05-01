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
#include "statistics.h"
#include "initial_positions.h"
#include "parameters.h"
#include "proterra.h"
#include "solution.h"
#include "util.h"
#include "validate.h"
#include <fstream>
#include <iostream>

timer<> stats::time;
solution stats::global_best;
vi stats::river_violations, stats::river_values;
vd stats::size_ratios, stats::values, stats::cum_times, stats::times;
int stats::nrepl = 0, stats::num_mutations = 0, stats::num_crossovers = 0,
    stats::num_disconnected_crossover_lots = 0,
    stats::num_empty_crossover_lots = 0, stats::num_new_solutions = 0,
    stats::num_generations = 0;

template <typename T>
void pr_min_avg_max(const T& v, const string& name, int best, int worst) {
  pr("--{} min {:.2f} avg {:.2f} max {:.2f}\n", name.c_str(), (double)v[best],
     (double)avg(v), (double)v[worst]);
}

void stats::write_stats() {
  if (global_best.num_assigned != 0 and nrepl == 0) add_stats(global_best);

  viii obj;
  for (int i = 0; i < nrepl; ++i)
    obj.emplace_back(river_values[i], size_ratios[i], values[i]);
  int best = min_index(obj), worst = max_index(obj);

  pr("FINISHED, printing statistics...\n\n", nrepl);
  pr("--instance {}\n", prm::instance_name.c_str());
  pr("--nland {}\n", pt::nland);
  pr("--lots {}\n", pt::lots);
  pr("--batch-size {}\n", prm::batch_size);
  pr("--time {:.2f}\n", time.seconds());
  pr("--repl {}\n", nrepl);
  pr_min_avg_max(values, "value", best, worst);
  pr_min_avg_max(river_violations, "river_violations", best, worst);
  pr_min_avg_max(river_values, "river_value", best, worst);
  pr_min_avg_max(times, "time", best, worst);
  pr_min_avg_max(size_ratios, "size_ratio", best, worst);
  pr("\n");

  if (prm::irace) {
    ld val;
    if (best == -1 or worst == -1)
      val = numeric_limits<double>::max();
    else
      val = ld(1LL << 30) * (ld)river_violations[best] +
            ld(1LL << 20) * (ld)size_ratios[best] + (ld)values[best];
    fmt::print("{}", val);
    return;
  }

  if (prm::naive) num_generations = nrepl;

  fmt::print("{} ", prm::instance_name);
  fmt::print("{} ", pt::c_size);
  fmt::print("{} ", pt::r_size);
  fmt::print("{} ", pt::lots);
  fmt::print("{} ", pt::nland);
  fmt::print("{} ", pt::nriver);
  fmt::print("{} ", pt::r_size * pt::c_size - pt::nriver - pt::nland);
  fmt::print("{} ", pt::num_apt_classes);
  fmt::print("{} ", initial_positions::cc.size());
  fmt::print("{} ", prm::batch_size);
  fmt::print("{} ", prm::crossover_ratio);
  fmt::print("{} ", prm::keep_ratio);
  fmt::print("{} ", prm::new_rate);
  fmt::print("{:.2f} ", time.seconds());
  fmt::print("{} ", num_generations);
  fmt::print("{} ", prm::random_seed);

  assert(nrepl == (int)river_values.size());

  if (best == -1 or worst == -1) {
    for (int i = 0; i < 32; ++i) {
      fmt::print("NA ");
    }
  } else {
    fmt::print("{} ", num_new_solutions);
    fmt::print("{} ", num_crossovers);
    fmt::print("{} ", num_mutations);
    fmt::print("{} ", num_crossovers + num_mutations + num_new_solutions);

    fmt::print("{} ", best);
    fmt::print("{:.2f} ", cum_times[best]);
    fmt::print("{:.2f} ", avg(times));
    fmt::print("{:.2f} ", std_dev(times));
    fmt::print("{:.2f} ", times[0]);

    fmt::print("{:.2f} ", values[best]);
    fmt::print("{:.2f} ", values[worst]);
    fmt::print("{:.2f} ", avg(values));
    fmt::print("{:.2f} ", std_dev(values));
    fmt::print("{:.2f} ", values[0]);

    fmt::print("{} ", 0);
    fmt::print("{} ", 0);
    fmt::print("{:.2f} ", 0.0);
    fmt::print("{:.2f} ", 0.0);
    fmt::print("{} ", 0);

    fmt::print("{} ", 0);
    fmt::print("{} ", 0);
    fmt::print("{:.2f} ", 0.0);
    fmt::print("{:.2f} ", 0.0);
    fmt::print("{} ", 0);

    fmt::print("{} ", river_values[best]);
    fmt::print("{} ", river_values[worst]);
    fmt::print("{:.2f} ", avg(river_values));
    fmt::print("{:.2f} ", std_dev(river_values));
    fmt::print("{} ", river_values[0]);

    fmt::print("{} ", river_violations[best]);
    fmt::print("{} ", river_violations[worst]);
    fmt::print("{:.2f} ", avg(river_violations));
    fmt::print("{:.2f} ", std_dev(river_violations));
    fmt::print("{} ", river_violations[0]);

    fmt::print("{} ", size_ratios[best]
                          ? fmt::format("{:.2f}", size_ratios[best])
                          : "NA");
    fmt::print("{} ", size_ratios[worst]
                          ? fmt::format("{:.2f}", size_ratios[worst])
                          : "NA");
    fmt::print("{} ", avg(size_ratios) ? fmt::format("{:.2f}", avg(size_ratios))
                                       : "NA");
    fmt::print("{:.2f} ", std_dev(size_ratios));
    fmt::print("{} ",
               size_ratios[0] ? fmt::format("{:.2f}", size_ratios[0]) : "NA");

    fmt::print("{} ", 0);
    fmt::print("{} ", 0);
    fmt::print("{:.2f} ", 0.0);
    fmt::print("{:.2f} ", 0.0);
    fmt::print("{} ", 0);

    fmt::print("{} ", 0);
    fmt::print("{} ", 0);
    fmt::print("{:.2f} ", 0.0);
    fmt::print("{:.2f} ", 0.0);
  }
  fmt::print("\n");
  if (global_best.num_assigned > 0 and prm::png.size()) {
    global_best.write_to_png(prm::png);
  }
}

void stats::add_stats(const solution& s) {
  validate_solution(s);
  int i = nrepl++;
  assert(abs(s.of.value - (ll)s.of.compute_value_brute_force()) < 2);
  values.push_back(sqrt(s.of.value / (double)pt::lots));
  river_violations.push_back(s.rc.violations);
  river_values.push_back(s.rc.value);
  int ba = max(s.area), sa = min(s.area);
  size_ratios.push_back(sa ? (double)ba / (double)sa : 0);
  cum_times.push_back(time.seconds());
  times.push_back(i ? time.seconds() - cum_times[i - 1] : time.seconds());
  pr("value {:.2f}\nrivers {}\nsize_ratio {:.2f}\n\ntime {:.2f}\n", values[i],
     river_values[i], size_ratios[i], cum_times[i]);
  pr("\n");

  if (global_best.num_assigned == 0 or s < global_best) global_best = s;
}
