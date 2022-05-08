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
    std::vector<Mat34> *models)
{
  assert(3 == bearing_vectors.rows());
  assert(3 == X.rows());
  assert(bearing_vectors.cols() == X.cols());
  std::vector<std::tuple<Mat3, Vec3>> rotation_translation_solutions;

  unsigned nsols;
  double degen;
	double rotation_translation_solutions[RT_MAX_LEN][4][3];

  p2pt<double>::pose_from_point_tangents(
    sample_gama1, sample_tgt1,
    sample_gama2, sample_tgt2,
    sample_Gama1, sample_Tgt1,
    sample_Gama2, sample_Tgt2,
    &rotation_translation_solutions, &nsols, &degen
  );

	for (unsigned i = 0; i < nsols; ++i) {
    Mat34 P;
    rotation_translation_solutions[i]
    P_From_KRt(Mat3::Identity(),           // intrinsics
                std::get<0>(rot_trans_it), // rotation
                std::get<1>(rot_trans_it), // translation
                &P);
    models->push_back(P);
  }
};

} // namespace euclidean_resection
} // namespace openMVG
