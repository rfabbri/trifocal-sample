// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/multiview/solver_resection_kernel.hpp"
#include "openMVG/multiview/solver_resection_p3p.hpp"
#include "openMVG/multiview/solver_resection_up2p_kukelova.hpp"
#include "openMVG/multiview/solver_resection_p2pt_fabbri.hpp"
#include "openMVG/multiview/test_data_sets.hpp"

#include "testing/testing.h"

using namespace openMVG;

TEST(Resection_Kernel_DLT, Multiview) {

  const int nViews = 3;
  const int nbPoints = 10;
  const NViewDataSet d =
    NRealisticCamerasRing(nViews, nbPoints,
                          nViewDatasetConfigurator(1,1,0,0,5,0)); // Suppose a camera with Unit matrix as K

  // Solve the problem and check that fitted value are good enough
  for (int nResectionCameraIndex = 0; nResectionCameraIndex < nViews; ++nResectionCameraIndex)
  {
    const Mat x = d._x[nResectionCameraIndex];
    const Mat X = d._X;
    openMVG::resection::kernel::PoseResectionKernel kernel(x, X);

    std::vector<Mat34> Ps;
    kernel.Fit({0,1,2,3,4,5}, &Ps);
    for (Mat::Index i = 0; i < x.cols(); ++i) {
      EXPECT_NEAR(0.0, kernel.Error(i, Ps[0]), 1e-8);
    }

    CHECK_EQUAL(1, Ps.size());

    // Check that Projection matrix is near to the GT:
    Mat34 GT_ProjectionMatrix = d.P(nResectionCameraIndex).array()
                                / d.P(nResectionCameraIndex).norm();
    Mat34 COMPUTED_ProjectionMatrix = Ps[0].array() / Ps[0].norm();
    EXPECT_MATRIX_NEAR(GT_ProjectionMatrix, COMPUTED_ProjectionMatrix, 1e-8);
  }
}

TEST(P3P_Kneip_CVPR11, Multiview) {

  const int nViews = 3;
  const int nbPoints = 3;
  const NViewDataSet d =
    NRealisticCamerasRing(nViews, nbPoints,
                          nViewDatasetConfigurator(1,1,0,0,5,0)); // Suppose a camera with Unit matrix as K

  // Solve the problem and check that fitted value are good enough
  for (int nResectionCameraIndex = 0; nResectionCameraIndex < nViews; ++nResectionCameraIndex)
  {
    const Mat x = d._x[nResectionCameraIndex];
    const Mat bearing_vectors = (d._K[0].inverse() * x.colwise().homogeneous()).colwise().normalized();
    const Mat X = d._X;
    openMVG::euclidean_resection::PoseResectionKernel_P3P_Kneip kernel(bearing_vectors, X);

    std::vector<Mat34> Ps;
    kernel.Fit({0,1,2}, &Ps); // 3 points sample are required, lets take the first three

    bool bFound = false;
    char index = -1;
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
}

TEST(P3P_Ke_CVPR17, Multiview) {

  const int nViews = 3;
  const int nbPoints = 3;
  const NViewDataSet d =
    NRealisticCamerasRing(nViews, nbPoints,
                          nViewDatasetConfigurator(1,1,0,0,5,0)); // Suppose a camera with Unit matrix as K



  // Solve the problem and check that fitted value are good enough
  for (int nResectionCameraIndex = 0; nResectionCameraIndex < nViews; ++nResectionCameraIndex)
  {
    const Mat x = d._x[nResectionCameraIndex];
    const Mat bearing_vectors = (d._K[0].inverse() * x.colwise().homogeneous()).colwise().normalized();
    const Mat X = d._X;
    openMVG::euclidean_resection::PoseResectionKernel_P3P_Ke kernel(bearing_vectors, X);
    std::vector<Mat34> Ps;
    kernel.Fit({0,1,2}, &Ps); // 3 points sample are required, lets take the first three

    bool bFound = false;
    char index = -1;
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
}


TEST(P3P_Nordberg_ECCV18, Multiview) {

  const int nViews = 3;
  const int nbPoints = 3;
  const NViewDataSet d =
    NRealisticCamerasRing(nViews, nbPoints,
                          nViewDatasetConfigurator(1,1,0,0,5,0)); // Suppose a camera with Unit matrix as K



  // Solve the problem and check that fitted value are good enough
  for (int nResectionCameraIndex = 0; nResectionCameraIndex < nViews; ++nResectionCameraIndex)
  {
    const Mat x = d._x[nResectionCameraIndex];
    const Mat bearing_vectors = (d._K[0].inverse() * x.colwise().homogeneous()).colwise().normalized();
    const Mat X = d._X;
    openMVG::euclidean_resection::PoseResectionKernel_P3P_Nordberg kernel(bearing_vectors, X);

    std::vector<Mat34> Ps;
    kernel.Fit({0,1,2}, &Ps); // 3 points sample are required, lets take the first three

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
}

TEST(UP2PSolver_Kukelova, Multiview) {

  const int nViews = 6;
  const int nbPoints = 2;
  const NViewDataSet d =
    NRealisticCamerasRing(nViews, nbPoints,
                          nViewDatasetConfigurator(1,1,0,0,5,0)); // Suppose a camera with Unit matrix as K

  // Solve the problem and check that fitted value are good enough
  for (int nResectionCameraIndex = 0; nResectionCameraIndex < nViews; ++nResectionCameraIndex)
  {
    const Mat x = d._x[nResectionCameraIndex];
    const Mat bearing_vectors = (d._K[0].inverse() * x.colwise().homogeneous()).colwise().normalized();
    const Mat X = d._X;
    openMVG::euclidean_resection::PoseResectionKernel_UP2P_Kukelova kernel(bearing_vectors, X);

    std::vector<Mat34> Ps;
    kernel.Fit({0,1}, &Ps); // 2 points sample are required, lets take the first three

    bool bFound = false;
    size_t index = -1;
      for (size_t i = 0; i < Ps.size(); ++i)  {
      const Mat34 GT_ProjectionMatrix = d.P(nResectionCameraIndex).array()
        / d.P(nResectionCameraIndex).norm();
      const Mat34 COMPUTED_ProjectionMatrix = Ps[i].array() / Ps[i].norm();
      std::cout << "GT:\n " << GT_ProjectionMatrix << std::endl;
      std::cout << "Computed:\n " << COMPUTED_ProjectionMatrix << std::endl;
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
}

TEST(P2Pt_Fabbri_PAMI20, Multiview) {
  static constexpr double sample_gama1[3] = {0.51515532818894982, 0.1487011661471217, 1};
  static constexpr double sample_tgt1[3]  = {0.527886693031222, 0.84931480578202578, 0};
  static constexpr double sample_gama2[3] = {0.16081537437895527, -0.48875444114985156, 1};
  static constexpr double sample_tgt2[3]  = {-0.27224516854045233, -0.96222791905368288, 0};
  static constexpr double sample_Gama1[3] = {-0.43359202230568356, 3.5783969397257605, -1.3498869401565212};
  static constexpr double sample_Tgt1[3]  = {0.70708731473372055, 0.69669519863759266, -0.12100962580713076};
  static constexpr double sample_Gama2[3] = {0.34262446653864992, 2.7694370298848772, 3.0349234663318545};
  static constexpr double sample_Tgt2[3]  = {-0.041895437077508819, -0.13618523302227314, 0.98979712803117059};

  static constexpr double sample_Rot[3][3] = {
    -0.917851347078651,    0.362140358287978,  -0.162490816863487,
     0.0846221494628621,  -0.221430418224084,  -0.971497638548544,
    -0.387798912435553,   -0.905440738416481,   0.172595112125579
  };
  static constexpr double sample_Transl[3] = {
    0.862173320368137,  0.318765239858977, 8.6923117036948,
  };

  Mat bearing_vectors(3,2);
  Mat tangent_vectors(3,2);
  Mat X(3,2);
  Mat T(3,2);
  
  for (unsigned i=0; i < 3; ++i) {
    bearing_vectors(i,0) = sample_gama1[i];
    bearing_vectors(i,1) = sample_gama2[i];
    tangent_vectors(i,0) = sample_tgt1[i];
    tangent_vectors(i,1) = sample_tgt2[i];
    X(i,0) = sample_Gama1[i];
    X(i,1) = sample_Gama2[i];
    T(i,0) = sample_Tgt1[i];
    T(i,1) = sample_Tgt2[i];
  }
  
  std::vector<Mat34> Ps;
  openMVG::euclidean_resection::P2PtSolver_Fabbri::Solve(
      bearing_vectors, tangent_vectors, X, T, &Ps);
  
  Mat34 GT_ProjectionMatrix;
  for (unsigned j = 0 ; j < 3; ++j)
    for (unsigned k = 0 ; k < 3; ++k)
      GT_ProjectionMatrix(j,k) = sample_Rot[j][k];
  for (unsigned j = 0 ; j < 3; ++j)
      GT_ProjectionMatrix(j,3) = sample_Transl[j];

  GT_ProjectionMatrix /= GT_ProjectionMatrix.norm();
  
  bool bFound = false;
  size_t index = -1;
  for (size_t i = 0; i < Ps.size(); ++i)  {
    Mat34 COMPUTED_ProjectionMatrix = Ps[i].array() / Ps[i].norm();
    if ( NormLInfinity(GT_ProjectionMatrix - COMPUTED_ProjectionMatrix) < 1e-3 ) {
      bFound = true;
      index = i;
    }
  }
  EXPECT_TRUE(bFound);
}

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
