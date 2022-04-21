// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.
//
//:\file
//\author Gabriel ANDRADE Rio de Janeiro State U.
//\author Ricardo Fabbri Rio de Janeiro State U. (rfabbri.github.io) 
//\date Tue Jun  1 11:55:58 -03 2021
//\author Pierre MOULON
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_MULTIVIEW_SOLVER_TRIFOCAL_THREE_POINT_HPP
#define OPENMVG_MULTIVIEW_SOLVER_TRIFOCAL_THREE_POINT_HPP

#include <vector>
#include "openMVG/multiview/trifocal/trifocal_model.hpp"

namespace openMVG {
namespace trifocal {
  
struct Trifocal3PointPositionTangentialSolver {
  enum { MINIMUM_SAMPLES = 3 };
  enum { MAX_MODELS = 312 };

  static void Solve(
      const Mat &datum_0,
      const Mat &datum_1,
      const Mat &datum_2,
      std::vector<trifocal_model_t> *trifocal_tensor);
};

} // namespace trifocal
} // namespace OpenMVG

#define  // OPENMVG_MULTIVIEW_SOLVER_TRIFOCAL_THREE_POINT_HPP
