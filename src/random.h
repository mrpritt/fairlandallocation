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
#include <climits>
#include <cstdint>
#include <ctime>
#include <random>
#include <unistd.h>

struct random_number_generator {
  std::mt19937 engine;

  void seed(int64_t seed) { engine.seed(seed); }

  int64_t seed_unique() {
    int64_t a = (int64_t)clock(), b = (int64_t)time(nullptr),
            c = (int64_t)getpid();
    a = (a - b - c) ^ (c >> 13);
    b = (b - c - a) ^ (a << 8);
    c = (c - a - b) ^ (b >> 13);
    a = (a - b - c) ^ (c >> 12);
    b = (b - c - a) ^ (a << 16);
    c = (c - a - b) ^ (b >> 5);
    a = (a - b - c) ^ (c >> 3);
    b = (b - c - a) ^ (a << 10);
    c = (c - a - b) ^ (b >> 15);
    seed(c);
    return c;
  }

  double rand_double(double from, double to) {
    static std::uniform_real_distribution<double> d;
    return d(engine, decltype(d)::param_type{from, to});
  }

  int rand_int(int from, int to) {
    static std::uniform_int_distribution<int> d;
    return d(engine, decltype(d)::param_type{from, to});
  }

  int rand() { return rand_int(0, INT_MAX); }
};

extern random_number_generator rng;

template <typename T>
inline void choice(const std::vector<T>& v, std::vector<T>& result, int k) {
  int n = v.size();
  result.resize(k);
  for (int i = 0; i < k; ++i)
    result[i] = v[i];
  for (int i = k + 1; i < n; ++i) {
    int j = rng.rand_int(0, i);
    if (j < k) result[j] = v[i];
  }
}
