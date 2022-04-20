// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.
//
//:\file
//\author Ricardo Fabbri Rio de Janeiro State U. (rfabbri.github.io) 
//\author Gabriel ANDRADE Rio de Janeiro State U.
//\author Pierre MOULON
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_MULTIVIEW_TRIFOCAL_MODEL_HPP
#define OPENMVG_MULTIVIEW_TRIFOCAL_MODEL_HPP

#include "openMVG/numeric/eigen_alias_definition.hpp"

namespace openMVG {
namespace trifocal {
  
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION_INITIALIZER_LIST(std::array<openMVG::Mat34,3>)
using trifocal_model_t = std::array<Mat34, 3>;

} // namespace trifocal
} // namespace OpenMVG

#endif  // OPENMVG_GEOMETRY_TRIFOCAL_MODEL_HPP
