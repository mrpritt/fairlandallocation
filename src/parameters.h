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
#include <cstdint>
#include <string>

struct prm {
  static void parse_cmd_line(int argc, char** argv);

  static double new_rate;

  static string desc_nbs;
  static int neighborhood_size;

  static string desc_bat;
  static int batch_size;

  static string desc_if;
  static string input_filename;
  static string instance_name;

  static string desc_si;
  static string solution_input;

  static string desc_png;
  static string png;

  static string desc_naive;
  static bool naive;

  static string desc_mgen;
  static int max_generations;

  static string desc_seed;
  static ll random_seed;

  static string desc_tl;
  static int time_limit_seconds;

  static string desc_msr;
  static int maximum_size_ratio;

  static string desc_dc;
  static bool do_crossover;

  static string desc_dm;
  static bool do_mutation;

  static string desc_ps;
  static int pop_size;

  static string desc_kr;
  static double keep_ratio;

  static string desc_cr;
  static double crossover_ratio;

  static string desc_gra;
  static double mutation_greedy_alpha;

  static string desc_rest;
  static int restart;

  static string desc_brush;
  static int mutation_brush_size;

  static string desc_tourn;
  static int tournament_size;

  static string desc_irace;
  static bool irace;
};
