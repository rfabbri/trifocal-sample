// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2020, Viktor Larsson
// Copyright (c) 2020, Pierre MOULON

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/multiview/solver_resection_p2pt_fabbri.hpp"
#include "openMVG/numeric/numeric.h"

using namespace openMVG;

namespace openMVG {
namespace euclidean_resection {


void UP2PSolver_Kukelova::Solve
(
  const Mat & x,
  const Mat & X,
  std::vector<Mat34> * models
)
{
  // XXX call Ariel's code
}

} // namespace euclidean_resection
} // namespace openMVG
