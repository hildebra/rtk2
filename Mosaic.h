#pragma once

#include "options.h"

// Discover likely split/mosaic gene clusters from MGS membership, minimap2 PAF
// relatedness, and an abundance/coverage matrix. The mode is report-only.
void runMosaicMode(options* opts);
