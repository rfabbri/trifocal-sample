//:\file
//\author Ricardo Fabbri Rio de Janeiro State U. (rfabbri.github.io) 
//\date Tue Jun  1 11:55:58 -03 2021
//\author Gabriel ANDRADE Rio de Janeiro State U.
//\author Pierre MOULON
#ifndef OPENMVG_MULTIVIEW_SOLVER_TRIFOCAL_THREE_POINT_HPP
#define OPENMVG_MULTIVIEW_SOLVER_TRIFOCAL_THREE_POINT_HPP

#include <iostream>
#include <array>
#include <vector>
#include "openMVG/numeric/extract_columns.hpp"
#include "openMVG/numeric/eigen_alias_definition.hpp"
#include "openMVG/multiview/trifocal_model.hpp"


namespace openMVG {
namespace trifocal {
  
//------------------------------------------------------------------------------
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
