// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2012, 2013 Pierre MOULON.
//
// Author: Ariel Kovaljski and Ricardo Fabbri

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/multiview/solver_resection_p2pt_fabbri.hpp"
#include "testing/testing.h"

using namespace openMVG;

/// XXX Ariel: incluir rotacao, translacao veerdadeiros

TEST(P2Pt_Fabbri_PAMI20, Multiview) {

  // Solve the problem and check that fitted value are good enough
  const Mat x = d._x[nResectionCameraIndex];
  const Mat bearing_vectors = (d._K[0].inverse() * x.colwise().homogeneous()).colwise().normalized();
  const Mat X = d._X;
  
  openMVG::euclidean_resection::PoseResectionKernel_P3P_Nordberg kernel(bearing_vectors, X);

  std::vector<Mat34> Ps;
  kernel.Fit({0,1}, &Ps); // 3 points sample are required, lets take the first three

  bool bFound = false;
  size_t index = -1;
  for (size_t i = 0; i < Ps.size(); ++i)  {
    Mat34 GT_ProjectionMatrix = d.P(nResectionCameraIndex).array()
    / d.P(nResectionCameraIndex).norm();
    Mat34 COMPUTED_ProjectionMatrix = Ps[i].array() / Ps[i].norm();
    if ( NormLInfinity(GT_ProjectionMatrix - COMPUTED_ProjectionMatrix) < 1e-8 )
    {
      bFound = true;
      index = i;
    }
  }
  EXPECT_TRUE(bFound);

  // Check that for the found matrix the residual is small
  for (Mat::Index i = 0; i < x.cols(); ++i) {
    EXPECT_NEAR(0.0, kernel.Error(i, Ps[index]), 1e-8);
  }
}

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
