// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2018 Michael Persson
// Adapted to openMVG by Romain Janvier and Pierre Moulon

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/multiview/solver_resection_p3p_nordberg.hpp"
#include "openMVG/multiview/projection.hpp"

#include <array>

namespace openMVG
{
namespace euclidean_resection
{

// p2pt.hxx here

void P2PtSolver_Fabbri::Solve(
    const Mat &bearing_vectors,
    const Mat &tangent_vectors,
    const Mat &X, // 3D points
    const Mat &T, // 3D tangents
    std::vector<Mat34> *models)
{
  assert(2 == bearing_vectors.rows());
  assert(2 == X.rows());
  assert(bearing_vectors.cols() == X.cols());

  unsigned nsols;
  double degen;
	double rotation_translation_solutions[RT_MAX_LEN][4][3];

  p2pt<double>::pose_from_point_tangents(
    bearing_vectors.col(0).data(), tangent_vectors.col(0).data(),
    bearing_vectors.col(1).data(), tangent_vectors.col(1).data(),
    X.col(0).data(), T.col(0).data(),
    X.col(1).data(), T.col(1).data(),
    &rotation_translation_solutions, &nsols, &degen
  );

	for (unsigned i = 0; i < nsols; ++i) {
    Mat34 P;
    for (unsigned j = 0 ; j < 3; ++i)
      for (unsigned k = 0 ; k < 3; ++k)
        P(j,k) = rotation_translation_solutions[i][j][k];

    for (unsigned k = 0 ; k < 3; ++k)
      P(k,3) = rotation_translation_solutions[i][3][k];
    models->push_back(P);
  }
};

} // namespace euclidean_resection
} // namespace openMVG
