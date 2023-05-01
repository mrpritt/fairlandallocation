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
#pragma once
#include "defines.h"
#include "parameters.h"
#include "solution.h"
#include "timer.h"
#include <string>
#include <vector>

struct solution;

struct stats {
  static void write_stats();
  static void add_stats(const solution& ds);
  inline static bool time_limit_exceeded() {
    return time.seconds() >= prm::time_limit_seconds;
  }
  static timer<> time;
  static solution global_best;
  static int num_generations;
  static int num_new_solutions;
  static int num_empty_crossover_lots;
  static int num_disconnected_crossover_lots;
  static int num_crossovers;
  static int num_mutations;
  static int nrepl;
  static vd values;
  static vi river_violations;
  static vi river_values;
  static vd size_ratios;
  static vd times, cum_times;
};
