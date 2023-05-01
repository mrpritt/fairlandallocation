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
#include "../../src/random.h"
#include "../../src/util.h"
#include "FastNoise.h"
#include "hsv_rgb.h"
#include "lodepng.h"
#include <algorithm>
#include <boost/program_options.hpp>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <queue>
#include <random>
#include <string>
#include <vector>

#define INVALID_RIVER (1 << 20)

using namespace std;
using uchar = unsigned char;

random_number_generator rng;
FastNoise noise;
vvi nbs;

int apt, w, h, value_seed, river_seed, river_pct, lots;
bool png;
string output_filename;
string visualize_filename;

int river_soil_width = 3;
int padding = 10;
int min_river_width = 2, max_river_width = 5;

void parse_cmd_line(int argc, char** argv) {
  namespace po = boost::program_options;
  po::options_description desc("Proterra instance generator.");
  desc.add_options()("help", "")("apt", po::value<int>(&apt)->default_value(5),
                                 "number of aptitude classes (soil qualities)")(
      "w", po::value<int>(&w)->default_value(100),
      "map width")("h", po::value<int>(&h)->default_value(100), "map height")(
      "value_seed", po::value<int>(&value_seed)->default_value(0),
      "random value_seed")("river_seed",
                           po::value<int>(&river_seed)->default_value(0),
                           "random river_seed")(
      "river_pct", po::value<int>(&river_pct)->default_value(10),
      "percentage of cells that should be rivers")(
      "lots", po::value<int>(&lots)->default_value(10), "number of lots")(
      "out", po::value<string>(&output_filename)->default_value(""),
      "output filename; if unset, a default name will be used")(
      "png",
      "if this option is set, will generate a png image of the instance")(
      "visualize", po::value<string>(&visualize_filename)->default_value(""),
      "use this option to generate the png image of an existing instance");

  po::variables_map vm;
  try {
    po::store(po::parse_command_line(argc, argv, desc), vm);
    if (vm.count("help")) {
      cout << desc << endl;
      exit(EXIT_SUCCESS);
    } else {
      po::notify(vm);
    }
    png = vm.count("png");

    rng.seed_unique();
    if (value_seed == 0) {
      value_seed = rng.rand_int(1000, 100000);
    }
    if (river_seed == 0) {
      river_seed = rng.rand_int(1000, 100000);
    }

    pr("value_seed: {}\n", value_seed);
    pr("river_seed: {}\n", river_seed);
    if (output_filename.empty() and visualize_filename.empty()) {
      output_filename =
          str(w) + "x" + str(h) + "k" + str(lots) + "s" + str(value_seed);
    }
  } catch (po::error& e) {
    fmt::print("error parsing command line: {}\n", e.what());
    cout << desc << endl;
    exit(EXIT_FAILURE);
  }
}

double interp(double val, double min1, double max1, double min2, double max2) {
  return min2 + (max2 - min2) * (val - min1) / (max1 - min1);
}

void save_png(const vd& values, const vi& rivers, int nw, int nh, int padding) {
  double max_val = -1, min_val = (1 << 28);
  for (auto i : values)
    if (i != -1 and i != -2) {
      max_val = max(max_val, i);
      min_val = min(min_val, i);
    }

  vbyte img(nw * nh * 4);
  for (int i = 0; i < nh; ++i) {
    for (int j = 0; j < nw; ++j) {
      int index = (i + padding) * w + j + padding;
      hsv color;
      if (rivers[index] >= 0) {
        color.h = 215.0;
        color.s = 0.8627;
        color.v = 0.8745;
      } else if (values[index] == -1) {
        color.h = color.s = color.v = 0; /* black */
      } else if (values[index] == -2) {
        color.h = color.s = 0;
        color.v = 1.0; /* white */
      } else {
        color.h = 33.0;
        color.s = 0.68;
        color.v = interp(values[index], min_val, max_val, 0.2, 0.8);
        assert(color.v >= 0.19 and color.v <= 0.81);
      }
      rgb color2 = hsv2rgb(color);
      img[i * nw * 4 + j * 4 + 0] = 255.0 * color2.r;
      img[i * nw * 4 + j * 4 + 1] = 255.0 * color2.g;
      img[i * nw * 4 + j * 4 + 2] = 255.0 * color2.b;
      img[i * nw * 4 + j * 4 + 3] = 255;
    }
  }
  vbyte png;
  unsigned error = lodepng::encode(png, img, nw, nh);
  if (error) {
    fmt::print("enconder error {}: {}\n", error, lodepng_error_text(error));
    exit(EXIT_FAILURE);
  }
  if (output_filename.find(".png") == string::npos) output_filename += ".png";
  lodepng::save_file(png, output_filename);
}

void visualize() {
  ifstream f(visualize_filename);
  if (f.fail())
    throw runtime_error(
        fmt::format("could not open input file {}", visualize_filename));
  int d1;
  f >> w >> h >> d1;
  vd values(w * h, -1);
  vi rivers(w * h, -1);
  for (int i = 0; i < h; ++i)
    for (int j = 0; j < w; ++j) {
      int x;
      f >> x;
      if (x == -1)
        rivers[i * w + j] = 1;
      else
        values[i * w + j] = (x == 0 ? -1 : x);
    }

  if (output_filename.empty()) {
    size_t i = visualize_filename.find_last_of(".");
    if (i != string::npos)
      visualize_filename.erase(visualize_filename.begin() + i,
                               visualize_filename.end());
    i = visualize_filename.find_last_of("/");
    if (i != string::npos)
      visualize_filename.erase(visualize_filename.begin(),
                               next(visualize_filename.begin() + i));
    output_filename = visualize_filename;
  }
  save_png(values, rivers, w, h, 0);
}

void save_instance(const vd& values, const vi& rivers) {
  int nw = w - 2 * padding, nh = h - 2 * padding;

  ofstream f(output_filename + ".input");
  f << nh << " " << nw << " " << lots << " " << apt << " " << river_pct << endl;
  for (int i = 0; i < nh; ++i) {
    for (int j = 0; j < nw; ++j) {
      int index = (i + padding) * w + j + padding;
      if (rivers[index] != -1)
        f << -1 << " ";
      else
        f << (int)values[index] << " ";
    }
    f << endl;
  }

  if (png) save_png(values, rivers, nw, nh, padding);
}

vi make_rivers() {
  double lo = 0.0, hi = 1.0;
  int octaves = 10.0;
  noise.SetSeed(river_seed);
  rng.seed(river_seed);

  auto generate_rivers = [&](double frequency) {
    noise.SetFractalOctaves(octaves);
    noise.SetFrequency(frequency);
    vd m(w * h);
    for (int i = 0; i < w * h; ++i)
      m[i] = noise.GetValueFractal(i / w, i % w);

    vb is_river(w * h, false);
    vi rivers;
    for (int i = 0; i < w * h; ++i)
      if (not is_river[i])
        for (int j = 0; j < (int)nbs[i].size(); ++j)
          if ((m[i] > 0) != (m[nbs[i][j]] > 0)) {
            is_river[i] = is_river[nbs[i][j]] = true;
            rivers.push_back(i);
            rivers.push_back(nbs[i][j]);
            break;
          }

    vi cc(w * h, -1);
    vector<vi> cc_cells;
    int num_cc = 0;
    for (int i = 0; i < (int)rivers.size(); ++i) {
      assert(is_river[rivers[i]]);
      if (cc[rivers[i]] == -1) {
        cc[rivers[i]] = ++num_cc;
        cc_cells.push_back(vi());
        queue<int> q;
        q.push(rivers[i]);
        while (q.size()) {
          int v = q.front();
          cc_cells[cc_cells.size() - 1].push_back(v);
          q.pop();
          for (int u : nbs[v])
            if (is_river[u] and cc[u] == -1) {
              cc[u] = num_cc;
              q.push(u);
            }
        }
      }
    }
    for (int c = 0; c < num_cc; ++c) {
      int width = rng.rand_int(min_river_width, max_river_width);
      vi dist(w * h, 0);
      queue<int> q;
      for (auto u : cc_cells[c]) {
        dist[u] = width;
        q.push(u);
      }
      while (q.size()) {
        int v = q.front();
        q.pop();
        cc[v] = c;
        if (dist[v] == 0) continue;
        for (auto u : nbs[v])
          if (dist[u] < dist[v] - 1) {
            dist[u] = dist[v] - 1;
            q.push(u);
          }
      }
    }

    int num_cells_with_river = 0;
    for (int i = padding; i < h - padding; ++i)
      for (int j = padding; j < w - padding; ++j)
        num_cells_with_river += cc[i * w + j] >= 0;
    double pct = 100.0 * double(num_cells_with_river) /
                 double((w - 2 * padding) * (h - 2 * padding));
    return mp(pct, move(cc));
  };

  pair<double, vi> v;
  while (lo <= hi) {
    double mid = lo + (hi - lo) / 2.0;
    if (mid == 0) break;
    v = generate_rivers(mid);
    pr("value: {}; frequency: {}\n", v.first, mid);
    if (abs(v.first - river_pct) < 1.0)
      break;
    else if (v.first < river_pct)
      lo = mid + 0.0001;
    else
      hi = mid - 0.0001;
  }

  auto& rivers = v.second;
  return move(rivers);
}

vd make_values(const vi& rivers) {
  double frequency = 0.0010;
  int octaves = 16;
  noise.SetFrequency(frequency);
  noise.SetFractalOctaves(octaves);
  noise.SetSeed(value_seed);
  rng.seed(value_seed);

  vd apt(::apt);
  for (int i = 0; i < ::apt; ++i)
    apt[i] = 30.0 + 70.0 * rng.rand_double(0.0, 1.0);
  sort(apt.begin(), apt.end());

  pr_vec("apt_classes", apt, " ");

  vd values(w * h);
  double min_val = 2, max_val = -2;
  for (int i = 0; i < w * h; ++i) {
    values[i] = noise.GetValueFractal(i / w, i % w);
    min_val = min(min_val, values[i]);
    max_val = max(max_val, values[i]);
  }

  pr("min {}, max {}\n", min_val, max_val);
  for (int i = 0; i < w * h; ++i) {
    int index =
        int((values[i] - min_val) * (::apt - 1.0) / (max_val - min_val));
    index = min(index, ::apt - 2);
    values[i] = apt[index];
  }

  vi dist(w * h, 0);
  for (int i = 0; i < w * h; ++i) {
    if (rivers[i] != -1 and rivers[i] != INVALID_RIVER) {
      int width = rng.rand_int(min_river_width, max_river_width);
      queue<int> q;
      q.push(i);
      dist[i] = width;
      while (q.size()) {
        int v = q.front();
        q.pop();
        if (rivers[v] == -1) values[v] = apt[::apt - 1];
        if (dist[v] <= 1) continue;
        for (auto u : nbs[v])
          if (dist[u] < dist[v] - 1) {
            dist[u] = (rivers[u] == -1 ? dist[v] - 1 : width);
            q.push(u);
          }
      }
    }
  }
  return move(values);
}

int main(int argc, char** argv) {
  parse_cmd_line(argc, argv);
  if (not visualize_filename.empty()) {
    visualize();
    exit(EXIT_SUCCESS);
  }

  w += 2 * padding;
  h += 2 * padding;

  nbs.resize(w * h);
  int di[] = {0, 0, -1, 1, -1, 1, -1, 1}, dj[] = {-1, 1, 0, 0, -1, -1, 1, 1};
  for (int i = 0; i < h; ++i)
    for (int j = 0; j < w; ++j)
      for (int k = 0; k < 8; ++k) {
        int i2 = i + di[k], j2 = j + dj[k];
        if (i2 < 0 or i2 >= h or j2 < 0 or j2 >= w) continue;
        nbs[i * w + j].push_back(i2 * w + j2);
      }

  noise.SetFractalLacunarity(2.0);
  noise.SetFractalGain(0.5);
  noise.SetFractalType(FastNoise::FractalType::FBM);
  noise.SetNoiseType(FastNoise::NoiseType::SimplexFractal);
  noise.SetInterp(FastNoise::Interp::Quintic);
  noise.SetSeed(value_seed);

  auto rivers = make_rivers();
  auto values = make_values(rivers);
  save_instance(values, rivers);
}
