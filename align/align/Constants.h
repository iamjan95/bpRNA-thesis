#pragma once
#include <map>
#include <vector>
#include <string>
#include <cmath>
#include <string>
#include <algorithm>
#include <iostream>
#include <chrono>
using namespace std::chrono;

#define GAP -3
#define EXTEND -6

using namespace std;

enum Direction
{
    UP,
    LEFT
};




static map<pair<char, char>, double> gap_score = {
    {{'-', 'M'}, -2.9},
    {{'-', 'I'}, -2.8},
    {{'-', 'R'}, -2.9},
    {{'-', 'L'}, -2.9},
    {{'-', 'B'}, -2.9},
    {{'-', 'E'}, -0.5},
    {{'-', 'X'}, -1.7},
    {{'-', 'H'}, -1.8}
};

static map<pair<char, char>, double> extend_score = {
      {{'-', 'M'}, -5.7},
      {{'-', 'I'}, -5.6},
      {{'-', 'R'}, -5.8},
      {{'-', 'L'}, -5.8},
      {{'-', 'B'}, -5.8},
      {{'-', 'E'}, -1.0},
      {{'-', 'X'}, -3.4},
      {{'-', 'H'}, -3.6}
};


static map<pair<char, char>, double> score = {
    {{'M','M'}, 4}, {{'M','I'}, -4}, {{'M','R'}, -4}, {{'M','L'}, -4},
    {{'M','B'}, -4}, {{'M','E'}, -2}, {{'M','X'}, -2}, {{'M','H'}, -4},
    {{'I','I'}, 2}, {{'I','R'}, -4}, {{'I','L'}, -4}, {{'I','B'}, 0},
    {{'I','E'}, -2}, {{'I','X'}, -2}, {{'I','H'}, -4}, {{'R','R'}, 6},
    {{'R','L'}, -8}, {{'R','B'}, -4}, {{'R','E'}, -4}, {{'R','X'}, -2},
    {{'R','H'}, -4}, {{'L','L'}, 6}, {{'L','B'}, -4}, {{'L','E'}, -4},
    {{'L','X'}, -2}, {{'L','H'}, -4}, {{'B','B'}, 2}, {{'B','E'}, -2},
    {{'B','X'}, -2}, {{'B','H'}, -4}, {{'E','E'}, 2}, {{'E','X'}, -2},
    {{'E','H'}, -4}, {{'X','X'}, 2}, {{'X','H'}, -4}, {{'H','H'}, 4}
};
