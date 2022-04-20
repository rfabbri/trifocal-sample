//:\file
//\author Ricardo Fabbri Rio de Janeiro State U. (rfabbri.github.io) 
//\date Tue Apr 19 20:25:06 -03 2022
//\author Gabriel ANDRADE Rio de Janeiro State U.
//\author Pierre MOULON
#ifndef OPENMVG_MULTIVIEW_TRIFOCAL_METRICS_HPP
#define OPENMVG_MULTIVIEW_TRIFOCAL_METRICS_HPP

#include "openMVG/multiview/trifocal_model.hpp"

namespace openMVG {
namespace trifocal {
  
struct NormalizedSquaredPointReprojectionOneViewError {
    static double Error(
      const trifocal_model_t &tt,
      const Vec &bearing_0, // x,y,tangentialx,tangentialy
      const Vec &bearing_1,
      const Vec &bearing_2);
};

} // namespace trifocal
} // namespace OpenMVG
#endif  // OPENMVG_MULTIVIEW_TRIFOCAL_METRICS_HPP
