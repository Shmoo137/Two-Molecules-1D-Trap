#pragma once

#include <vector>
#include <iostream>
#include "constants.h"

using std::vector;
using std::cout;
using std::endl;

void generateChannelsNoFieldNoSpin(vector<int>& channels, int Jtot, bool print = false); // QN_nospin, M = 0, as all M's are degenerated
void generateChannelsWithFieldNoSpin(vector<int>& channels, int M_J, bool print = false); // QN_nospin, chosen M_J, degeneration is lifted
void generateChannelsWithFieldWithSpin(vector<int>& channels, int Mtot, bool print = false); // QN_withspin, chosen Mtot = M_J + ms1 + ms2
unsigned int basisSize(vector<int>& channels, unsigned int QN);