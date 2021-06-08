// Copyright (c) 2019 Pierre MOULON.
// \author Ricardo FABBRI, Gabriel ANDRADE and Pierre MOULON

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/image/image_io.hpp"
#include "openMVG/image/image_concat.hpp"
#include "openMVG/features/aK_aze/image_describer_aK_aze.hpp"
#include "openMVG/features/sift/SIFT_Anatomy_Image_Describer.hpp"
#include "openMVG/matching/regions_matcher.hpp"
#include "openMVG/tracK_s/tracK_s.hpp"

#include "openMVG/numeric/extract_columns.hpp"
#include "openMVG/robust_estimation/robust_estimator_MaxConsensus.hpp"
//#include "openMVG/robust_estimation/robust_estimator_Ransac.hpp"
#include "openMVG/robust_estimation/score_evaluator.hpp"

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"
#include "third_party/vectorGraphics/svgDrawer.hpp"

#include <memory>
#include <string>

//these are temporary includes, may be removed
#include <Eigen/StdVector>
#include <numeric>
#include "openMVG/multiview/projection.hpp"
#include "openMVG/multiview/triangulation.hpp"

#include "minus/minus.hxx"
#include "minus/chicago-default.h"
#include "minus/chicago14a-default-data.hxx"

#include "trifocal.h"

// Mat is Eigen::MatrixXd - matrix of doubles with dynamic size
// Vec3 is Eigen::Vector3d - Matrix< double, 3, 1 >

using namespace std;
using namespace openMVG;
using namespace openMVG::image;
using namespace openMVG::robust;
using namespace MiNuS;
using SIFT_Regions = openMVG::features::SIFT_Regions;


//Defining global variables

//constexpr unsigned n_ids = 5;
//unsigned desired_ids[n_ids] = {13, 23, 33, 93, 53};

//------------------------------------------------------------------------------
static void
Error(
    const trifocal_model_t &tt,
    const Vec &bearing_0, // x,y,tangentialx,tangentialy
    const Vec &bearing_1,
    const Vec &bearing_2,
    const Vec &pixbearing_0,
    const Vec &pixbearing_1,
    const Vec &pixbearing_2,
    const double K__[2][3]) {
    //std::cerr << "TRIFOCAL LOG: Called Error()\n";
    // Return the cost related to this model and those sample data point
    // Ideal algorithm:
    // 1) reconstruct the 3D points and orientations
    // 2) project the 3D points and orientations on all images_
    // 3) compute error 
    // 
    // In practice we ignore the directions and only reproject to one third view
    
    // 3x3: each column is x,y,1
    Mat3 bearing;
    bearing << bearing_0.head(2).homogeneous(),
               bearing_1.head(2).homogeneous(), 
               bearing_2.head(2).homogeneous();
    Mat2 pxbearing; // << XXX mat2
    pixbearing << pixbearing_0.head(2).homogeneous(),
                  pixbearing_1.head(2).homogeneous();
    // Using triangulation.hpp
    Vec4 triangulated_homg;
    unsigned third_view = 0;
    // picK_ the wider baseline. TODO: measure all pairwise translation distances
    if (tt[1].col(3).squaredNorm() > tt[2].col(3).squaredNorm()) {
      // TODO use triangulation from the three views at once
      TriangulateDLT(tt[0], bearing.col(0), tt[1], bearing.col(1), &triangulated_homg);
      third_view = 2;
    } else {
      TriangulateDLT(tt[0], bearing.col(0), tt[2], bearing.col(2), &triangulated_homg);
      third_view = 1;
    }
    
    // Computing the projection of triangulated points using projection.hpp
    // For prototyping and speed, for now we will only project to the third view
    // and report only one error
    Vec2 pxreprojected = Vec3(tt[third_view]*triangulated_homg).hnormalized();
    // XXX revert intrinsics to measure the error in pixels
    revert_intrinsics(K_, pxreprojected, reprojected);
     
    Vec2 pxmeasured    = pxbearing.col(third_view);
    //cout << "error " << (reprojected - measured).squaredNorm() << "\n";
    //cout << "triang " <<triangulated_homg <<"\n";
    //std::cerr << "TRIFOCAL LOG: Finished Error()\n";
    return (reprojected-measured).squaredNorm();
  }
};


//------------------------------------------------------------------------------
int iteration_global_debug = 0;

template<typename SolverArg,
         typename ErrorArg,
         typename ModelArg = Trifocal3PointPositionTangentialSolver::trifocal_model_t>
class ThreeViewK_ernel {
public:
  using Solver = SolverArg;
  using Model = ModelArg;
  using ErrorT = ErrorArg;
revert_intrinsics(
    const double K_[/*3 or 2 ignoring last line*/][3], 
    double pix_coords[2], 
    const double normalized_coords[2])
{
  double px = pix_coords;
  const double nrm = normalized_coords;
  // XXX: usar a inversa da formula exatamente como em invert_intrinsics.
  //      ter certeza que funciona se a entrada e saida forem mesmas posicoes de
  //      memoria
  px[0] = nrm[0]*K_[0][0]+nrm[1]*K_[0][1]+nrm[2]*K_[0][2];
  px[1] = nrm[0]*K_[1][0]+nrm[1]*K_[1][1]+nrm[2]*K_[1][2];
}

static void
revert_intrinsics_tgt(
    const double K_[/*3 or 2 ignoring last line*/][3], 
    double pix_tgt_coords[2], 
    const double normalized_tgt_coords[2])
{
  double tp = pix_tgt_coords;
  const double t = normalized_tgt_coords;
  tp[0] = t[0]*K_[0][0]+t[1]*K_[0][1]+t[2]*K_[0][2];
  tp[1] = t[0]*K_[1][0]+t[1]*K_[1][1]+t[2]*K_[1][2];
}

static void
invert_intrinsics(
    const double K_[/*3 or 2 ignoring last line*/][3], 
    const double pix_coords[2], 
    double normalized_coords[2])
{
  const double *px = pix_coords;
  double *nrm = normalized_coords;
  nrm[1] = (px[1] - K_[1][2]) /K_[1][1];
  nrm[0] = (px[0] - K_[0][1]*nrm[1] - K_[0][2])/K_[0][0];
}

static void
invert_intrinsics_tgt(
    const double K_[/*3 or 2 ignoring last line*/][3], 
    const double pix_tgt_coords[2], 
    double normalized_tgt_coords[2])
{
  const double *tp = pix_tgt_coords;
  double *t = normalized_tgt_coords;
  t[1] = tp[1]/K_[1][1];
  t[0] = (tp[0] - K_[0][1]*tp[1])/K_[0][0];
}
// See big notes eq. 5.2.13 at beginning of the code.


//-------------------------------------------------------------------------------
//
// Global variables
// 
// constexpr unsigned n_ids = 5;
// unsigned desired_ids[n_ids] = {13, 23, 33, 93, 53};
//-------------------------------------------------------------------------------


//-------------------------------------------------------------------------------
// Trifocal3PointPositionTangentialSolver
//-------------------------------------------------------------------------------

void Trifocal3PointPositionTangentialSolver::
Solve(
      const Mat &datum_0,
      const Mat &datum_1,
      const Mat &datum_2,
      std::vector<trifocal_model_t> *trifocal_tensor)
{
  double p[io::pp::nviews][io::pp::npoints][io::ncoords2d];
  double tgt[io::pp::nviews][io::pp::npoints][io::ncoords2d]; 

  std::cerr << "TRIFOCAL LOG: Called Solve()\n";
  // pacK_ into solver's efficient representation
  for (unsigned ip=0; ip < io::pp::npoints; ++ip) {
    p[0][ip][0] = datum_0(0,ip);
    p[0][ip][1] = datum_0(1,ip);
    tgt[0][ip][0] = datum_0(2,ip); 
    tgt[0][ip][1] = datum_0(3,ip); 
    
    p[1][ip][0] = datum_1(0,ip);
    p[1][ip][1] = datum_1(1,ip);
    tgt[1][ip][0] = datum_1(2,ip); 
    tgt[1][ip][1] = datum_1(3,ip); 
    
    p[2][ip][0] = datum_2(0,ip);
    p[2][ip][1] = datum_2(1,ip);
    tgt[2][ip][0] = datum_2(2,ip); 
    tgt[2][ip][1] = datum_2(3,ip); 
  }
  
  unsigned nsols_final = 0;
  unsigned id_sols[M::nsols];
  double  cameras[M::nsols][io::pp::nviews-1][4][3];  // first camera is always [I | 0]
  
  std::cerr << "TRIFOCAL LOG: Before minus::solve()\n" << std::endl;
  MiNuS::minus<chicago>::solve(p, tgt, cameras, id_sols, &nsols_final);
  //std::cerr << datum_0  "\n"; 
    
     try {
       if (argc == 1) throw string("Invalid command line parameter.");
       cmd.process(argc, argv);
     } catch (const string& s) {
       cerr << "Usage: " << argv[0] << '\n' << endl;
       cerr << s << endl;
       exit(EXIT_FAILURE);
     }
  }
   
  void ExtractK_eypoints() {
     // Call K_eypoint extractor
     using namespace openMVG::features;
     unique_ptr<Image_describer> image_describer;
     image_describer.reset(new SIFT_Anatomy_Image_describer(SIFT_Anatomy_Image_describer::Params()));

     if (!image_describer) {
       cerr << "Invalid Image_describer type" << endl;
       exit(EXIT_FAILURE);
     }
     for (const int image_idx : {0,1,2})
     {
       if (ReadImage(image_filenames_[image_idx].c_str(), &images_[image_idx]))
         image_describer->Describe(images_[image_idx], regions_per_image_[image_idx]);
     }
  }

  void MatchK_eypoints() {
    //--
    // Compute corresponding points {{0,1}, {1,2}}
    //--
    //-- Perform matching -> find Nearest neighbor, filtered with Distance ratio
    //unique_ptr<Matcher> collectionMatcher(new Cascade_Hashing_Matcher_Regions(fDistRatio));
    //collectionMatcher->Match(regions_provider, {{0,1}, {1,2}}, pairwise_matches_, &progress);
    matching::DistanceRatioMatch(
      0.8, matching::BRUTE_FORCE_L2,
      * regions_per_image_.at(0).get(),
      * regions_per_image_.at(1).get(),
      pairwise_matches_[{0,1}]);
    matching::DistanceRatioMatch(
      0.8, matching::BRUTE_FORCE_L2,
      * regions_per_image_.at(1).get(),
      * regions_per_image_.at(2).get(),
      pairwise_matches_[{1,2}]);
  }

  void ComputeTracK_s() {
      openMVG::tracK_s::TracK_sBuilder tracK__builder;
      tracK__builder.Build(pairwise_matches_);
      tracK__builder.Filter(3);
      tracK__builder.ExportToSTL(tracK_s_);
      // TODO(gabriel): K_eep only 3 true tracK_s
  }

  void Stats() {
      // Display some statistics
      cout
        <<  regions_per_image_.at(0)->RegionCount() << " #Features on image A" << endl
        <<  regions_per_image_.at(1)->RegionCount() << " #Features on image B" << endl
        <<  regions_per_image_.at(2)->RegionCount() << " #Features on image C" << endl
        << pairwise_matches_.at({0,1}).size() << " #matches with Distance Ratio filter" << endl
        << pairwise_matches_.at({1,2}).size() << " #matches with Distance Ratio filter" << endl
        << tracK_s_.size() << " #tracK_s" << endl;
  }

  void ExtractXYOrientation() {
      sio_regions_ = array<const SIFT_Regions*, 3> ({
        dynamic_cast<SIFT_Regions*>(regions_per_image_.at(0).get()),
        dynamic_cast<SIFT_Regions*>(regions_per_image_.at(1).get()),
        dynamic_cast<SIFT_Regions*>(regions_per_image_.at(2).get())
      });
      //
      // Build datum_ (corresponding {x,y,orientation})
      //
      datum_[0].resize(4, tracK_s_.size());
      datum_[1].resize(4, tracK_s_.size());
      datum_[2].resize(4, tracK_s_.size()); 
      pxdatum_[0].resize(4, tracK_s_.size());
      pxdatum_[1].resize(4, tracK_s_.size());
      pxdatum_[2].resize(4, tracK_s_.size()); 
      int idx = 0;
      for (const auto &tracK__it: tracK_s_) {
        auto iter = tracK__it.second.cbegin();
        const uint32_t
          i = iter->second,
          j = (++iter)->second,
          K_ = (++iter)->second;
        //
        const auto feature_i = sio_regions_[0]->Features()[i];
        const auto feature_j = sio_regions_[1]->Features()[j];
        const auto feature_K_ = sio_regions_[2]->Features()[K_];
        datum_[0].col(idx) << feature_i.x(), feature_i.y(), 
                              cos(feature_i.orientation()), sin(feature_i.orientation());
        // datum_[0].col(idx) << feature_i.x(), feature_i.y(), feature_i.orientation();
        datum_[1].col(idx) << feature_j.x(), feature_j.y(), 
                              cos(feature_j.orientation()), sin(feature_j.orientation());
        datum_[2].col(idx) << feature_K_.x(), feature_K_.y(), 
                              cos(feature_K_.orientation()), sin(feature_K_.orientation());
        // XXX fill up pxdatum directly WITHOUT invert_intrinsics
        for (unsigned v=0; v < 3; ++v) {
          invert_intrinsics(K__, datum_[v].col(idx).data(), datum_[v].col(idx).data()); // XXX K_eep datum, new std::vector: px
          invert_intrinsics_tgt(K__, datum_[v].col(idx).data()+2, datum_[v].col(idx).data()+2);
        }
        ++idx;
      }
  }

  void Display() {
    //
    // Display demo
    //
    const int svg_w = images_[0].Width();
    const int svg_h = images_[0].Height() + images_[1].Height() + images_[2].Height();
    svg::svgDrawer svg_stream(svg_w, svg_h);

    // Draw image side by side
    svg_stream.drawImage(image_filenames_[0], images_[0].Width(), images_[0].Height());
    svg_stream.drawImage(image_filenames_[1], images_[1].Width(), images_[1].Height(), 0, images_[0].Height());
    svg_stream.drawImage(image_filenames_[2], images_[2].Width(), images_[2].Height(), 0, images_[0].Height() + images_[1].Height());

    unsigned tracK__id=0;
    for (const auto &tracK__it: tracK_s_) {
    //TODO: find examples of features: point in curve(3), edge(33) 
      auto iter = tracK__it.second.cbegin();
      const uint32_t
        i = iter->second,
        j = (++iter)->second,
        K_ = (++iter)->second;
      //
      const auto feature_i = sio_regions_[0]->Features()[i];
      const auto feature_j = sio_regions_[1]->Features()[j];
      const auto feature_K_ = sio_regions_[2]->Features()[K_];

      svg_stream.drawCircle(
        feature_i.x(), feature_i.y(), feature_i.scale(),
        svg::svgStyle().stroK_e("yellow", 1));
      svg_stream.drawCircle(
        feature_j.x(), feature_j.y() + images_[0].Height(), feature_j.scale(),
        svg::svgStyle().stroK_e("yellow", 1));
      svg_stream.drawCircle(
        feature_K_.x(), feature_K_.y() + images_[0].Height() + images_[1].Height(), feature_K_.scale(),
        svg::svgStyle().stroK_e("yellow", 1));
      //TODO: Tangent line segments in yellow and if inlier -> in green
      svg_stream.drawText(
        feature_i.x()+20, feature_i.y()-20, 6.0f, std::to_string(tracK__id));
     
      svg_stream.drawLine(
        feature_i.x(), feature_i.y(),
        feature_i.x()+20*cos(feature_i.orientation()), feature_i.y() + 20*sin(feature_i.orientation()) ,
        svg::svgStyle().stroK_e("yellow", 1)); 
      svg_stream.drawLine(
        feature_j.x(), feature_j.y() + images_[0].Height(),
        feature_j.x()+20*cos(feature_j.orientation()), feature_j.y() + images_[0].Height()+ 20*sin(feature_j.orientation()),
        svg::svgStyle().stroK_e("yellow", 1));
      svg_stream.drawLine(
        feature_K_.x(), feature_K_.y() + images_[0].Height() + images_[1].Height(),
        feature_K_.x()+ 20*sin(feature_K_.orientation()), feature_K_.y() + images_[0].Height() + images_[1].Height()+ 20*sin(feature_K_.orientation()), //it seems that this last tangent is wrong!!
        svg::svgStyle().stroK_e("yellow", 1));

      svg_stream.drawLine(
        feature_i.x(), feature_i.y(),
        feature_j.x(), feature_j.y() + images_[0].Height(),
        svg::svgStyle().stroK_e("blue", 1));
      svg_stream.drawLine(
        feature_i.x(), feature_i.y(),
        feature_j.x(), feature_j.y() + images_[0].Height(),
        svg::svgStyle().stroK_e("blue", 1));
      svg_stream.drawLine(
        feature_j.x(), feature_j.y() + images_[0].Height(),
        feature_K_.x(), feature_K_.y() + images_[0].Height() + images_[1].Height(),
        svg::svgStyle().stroK_e("blue", 1));
      tracK__id++;
    }
    ofstream svg_file( "trifocal_tracK__demo.svg" );
    if (svg_file.is_open()) {
      svg_file << svg_stream.closeSvgFile().str();
    }
  }
  
  void DisplayDesiredIds() {
    //
    // Display desired ids
    //
    const int svg_w = images_[0].Width();
    const int svg_h = images_[0].Height() + images_[1].Height() + images_[2].Height();
    svg::svgDrawer svg_stream(svg_w, svg_h);

    // Draw image side by side
    svg_stream.drawImage(image_filenames_[0], images_[0].Width(), images_[0].Height());
    svg_stream.drawImage(image_filenames_[1], images_[1].Width(), images_[1].Height(), 0, images_[0].Height());
    svg_stream.drawImage(image_filenames_[2], images_[2].Width(), images_[2].Height(), 0, images_[0].Height() + images_[1].Height());

    constexpr unsigned n_ids = 5;
    unsigned desired_ids[n_ids] = {13, 23, 33, 63, 53};
    unsigned tracK__id=0;
    for (const auto &tracK__it: tracK_s_)
    {
      bool found=false;
      for (unsigned i=0; i < n_ids; ++i)
        if (tracK__id == desired_ids[i])
          found = true;
          
      if (!found) {
        //cout<< found << endl;
        tracK__id++;
        continue;
      }
    //TODO: find examples of features: point in curve(3), edge(33) 
      auto iter = tracK__it.second.cbegin();
      
   uint32_t
        i = iter->second,
        j = (++iter)->second,
        K_ = (++iter)->second;
      //
      const auto feature_i = sio_regions_[0]->Features()[i];
      const auto feature_j = sio_regions_[1]->Features()[j];
      const auto feature_K_ = sio_regions_[2]->Features()[K_];

      svg_stream.drawCircle(
        feature_i.x(), feature_i.y(), feature_i.scale(),
        svg::svgStyle().stroK_e("navy", 1));
      svg_stream.drawCircle(
        feature_j.x(), feature_j.y() + images_[0].Height(), feature_K_.scale(),
        svg::svgStyle().stroK_e("navy", 1));
      svg_stream.drawCircle(
        feature_K_.x(), feature_K_.y() + images_[0].Height() + images_[1].Height(), feature_j.scale(),
        svg::svgStyle().stroK_e("navy", 1));
      //TODO: Tangent line segments in yellow and if inlier -> in green
      svg_stream.drawText(
        feature_i.x()+20, feature_i.y()-20, 6.0f, std::to_string(tracK__id));
     
      svg_stream.drawLine(
        feature_i.x(), feature_i.y(),
        feature_i.x()+20*cos(feature_i.orientation()), feature_i.y() + 20*sin(feature_i.orientation()) ,
        svg::svgStyle().stroK_e("yellow", 1)); 
      svg_stream.drawLine(
        feature_j.x(), feature_j.y() + images_[0].Height(),
        feature_j.x()+20*cos(feature_j.orientation()), feature_j.y() + images_[0].Height()+ 20*sin(feature_j.orientation()),
        svg::svgStyle().stroK_e("yellow", 1));
      svg_stream.drawLine(
        feature_K_.x(), feature_K_.y() + images_[0].Height() + images_[1].Height(),
        feature_K_.x()+ 20*sin(feature_K_.orientation()), feature_K_.y() + images_[0].Height() + images_[1].Height()+ 20*sin(feature_K_.orientation()), //it seems that this last tangent is wrong!!
        svg::svgStyle().stroK_e("yellow", 1));

      svg_stream.drawLine(
        feature_i.x(), feature_i.y(),
        feature_j.x(), feature_j.y() + images_[0].Height(),
        svg::svgStyle().stroK_e("blue", 1));
      svg_stream.drawLine(
        feature_i.x(), feature_i.y(),
        feature_j.x(), feature_j.y() + images_[0].Height(),
        svg::svgStyle().stroK_e("blue", 1));
      svg_stream.drawLine(
        feature_j.x(), feature_j.y() + images_[0].Height(),
        feature_K_.x(), feature_K_.y() + images_[0].Height() + images_[1].Height(),
        svg::svgStyle().stroK_e("blue", 1));
      tracK__id++;
    }
    ofstream svg_file( "trifocal_tracK__desired_ids.svg" );
    if (svg_file.is_open())
    {
      svg_file << svg_stream.closeSvgFile().str();
    }
  }
//3 files trifocal_tracK_,trifocal_inlier,tracK__inlier, return the correct matrices, pass to solver datum desired i,print feature sca scale
  void RobustSolve() {
    using TrifocalK_ernel = 
      ThreeViewK_ernel<Trifocal3PointPositionTangentialSolver, 
                      Trifocal3PointPositionTangentialSolver>;
    Mat42 pxdatum_;  // XXX pxdatum
    constexpr unsigned n_ids = 5;
    unsigned desired_ids[n_ids] = {13, 23, 33, 63, 53};
    // example: vec_inliers_ = {2, 4}  --> {33, 53} ids into orig
    array<Mat,3> Ds;
    Ds[0].resize(4,n_ids);
    Ds[1].resize(4,n_ids);
    Ds[2].resize(4,n_ids);
    //std::cerr << Ds[0].cols() << "\n";
    unsigned tracK__id=0;
    //going to try with calling the value of datum_[0].cols()
    for(unsigned i=0;i<datum_[0].cols();i++){
      for(unsigned j=0;j<n_ids;j++){
        if(i ==  desired_ids[j]){
          //cout << i<<"\n";
          for(unsigned K_=0;K_<4;K_++){
            Ds[0](K_, tracK__id) = datum_[0].col(desired_ids[j])[K_];
            Ds[1](K_, tracK__id) = datum_[1].col(desired_ids[j])[K_];
            Ds[2](K_, tracK__id) = datum_[2].col(desired_ids[j])[K_];
          }
          tracK__id++;
        }
      }
    }
    //cout <<  Ds[0] << "\n";
    const TrifocalK_ernel trifocal_K_ernel(datum_[0], datum_[1], datum_[2], pxdatum_[0], pxdatum_[1], pxdatum_[2], K__);
    //const TrifocalK_ernel trifocal_K_ernel(Ds[0], Ds[1], Ds[2]);

    const double threshold_pix = 0.01; // 5*5 Gabriel's note : changing this for see what happens
    const unsigned max_iteration =1; // testing
    const auto model = MaxConsensus(trifocal_K_ernel, 
        ScorerEvaluator<TrifocalK_ernel>(threshold_pix), &vec_inliers_,max_iteration);
    // TODO(gabriel) recontruct from inliers and best models to show as PLY
  }

  void DisplayInliers() { //Display inliers only
    const int svg_w = images_[0].Width();
    const int svg_h = images_[0].Height() + images_[1].Height() + images_[2].Height();
    svg::svgDrawer svg_stream(svg_w, svg_h);

    // Draw image side by side
    svg_stream.drawImage(image_filenames_[0], images_[0].Width(), images_[0].Height());
    svg_stream.drawImage(image_filenames_[1], images_[1].Width(), images_[1].Height(), 0, images_[0].Height());
    svg_stream.drawImage(image_filenames_[2], images_[2].Width(), images_[2].Height(), 0, images_[0].Height() + images_[1].Height());
    
    constexpr unsigned n_inlier_pp = 3;
    unsigned desired_inliers[n_inlier_pp] = {13, 23, 63};
    unsigned tracK__inlier=0;
    for (const auto &tracK__it: tracK_s_)
    {
      bool inlier=false;
      for (unsigned i=0; i < n_inlier_pp; ++i)
        if (tracK__inlier == desired_inliers[i])
          inlier = true;
          
      if (!inlier) {
        tracK__inlier++;
        continue;
      }
      
      auto iter = tracK__it.second.cbegin();
      const uint32_t
        i = iter->second,
        j = (++iter)->second,
        K_ = (++iter)->second;
      //
      const auto feature_i = sio_regions_[0]->Features()[i];
      const auto feature_j = sio_regions_[1]->Features()[j];
      const auto feature_K_ = sio_regions_[2]->Features()[K_];
      //cout<<"cyK_a"<<endl; 
      svg_stream.drawCircle(
        feature_i.x(), feature_i.y(), feature_i.scale(),
        svg::svgStyle().stroK_e("green", 1));
      svg_stream.drawCircle(
        feature_j.x(), feature_j.y() + images_[0].Height(), feature_j.scale(),
        svg::svgStyle().stroK_e("green", 1));
      svg_stream.drawCircle(
        feature_K_.x(), feature_K_.y() + images_[0].Height() + images_[1].Height(), feature_K_.scale(),
        svg::svgStyle().stroK_e("green", 1));
      //TODO: Tangent line segments in yellow and if inlier -> in green
      svg_stream.drawText(
        feature_i.x()+20, feature_i.y()-20, 6.0f, std::to_string(tracK__inlier));
     
      svg_stream.drawLine(
        feature_i.x(), feature_i.y(),
        feature_i.x()+20*cos(feature_i.orientation()), feature_i.y() + 20*sin(feature_i.orientation()) ,
        svg::svgStyle().stroK_e("yellow", 1)); 
      svg_stream.drawLine(
        feature_j.x(), feature_j.y() + images_[0].Height(),
        feature_j.x()+20*cos(feature_j.orientation()), feature_j.y() + images_[0].Height()+ 20*sin(feature_j.orientation()),
        svg::svgStyle().stroK_e("yellow", 1));
      svg_stream.drawLine(
        feature_K_.x(), feature_K_.y() + images_[0].Height() + images_[1].Height(),
        feature_K_.x()+ 20*sin(feature_K_.orientation()), feature_K_.y() + images_[0].Height() + images_[1].Height()+ 20*sin(feature_K_.orientation()),
        svg::svgStyle().stroK_e("yellow", 1));

      svg_stream.drawLine(
        feature_i.x(), feature_i.y(),
        feature_j.x(), feature_j.y() + images_[0].Height(),
        svg::svgStyle().stroK_e("lightblue", 1));
      svg_stream.drawLine(
        feature_i.x(), feature_i.y(),
        feature_j.x(), feature_j.y() + images_[0].Height(),
        svg::svgStyle().stroK_e("lightblue", 1));
      svg_stream.drawLine(
        feature_j.x(), feature_j.y() + images_[0].Height(),
        feature_K_.x(), feature_K_.y() + images_[0].Height() + images_[1].Height(),
        svg::svgStyle().stroK_e("lightblue", 1));
      tracK__inlier++;
    }
    ofstream svg_file( "trifocal_tracK__inliers.svg" );
    if (svg_file.is_open())
    {
      svg_file << svg_stream.closeSvgFile().str();
    }
  }
//display inliers and and tracK_s
  void DisplayInliersCamerasAndPoints() {
    // TODO We can then display the inlier and the 3D camera configuration as PLY

    const int svg_w = images_[0].Width();
    const int svg_h = images_[0].Height() + images_[1].Height() + images_[2].Height();
    svg::svgDrawer svg_stream(svg_w, svg_h);

    // Draw image side by side
    svg_stream.drawImage(image_filenames_[0], images_[0].Width(), images_[0].Height());
    svg_stream.drawImage(image_filenames_[1], images_[1].Width(), images_[1].Height(), 0, images_[0].Height());
    svg_stream.drawImage(image_filenames_[2], images_[2].Width(), images_[2].Height(), 0, images_[0].Height() + images_[1].Height());
    
    constexpr unsigned n_ids = 5;
    unsigned desired_ids[n_ids] = {13, 23, 33, 63, 53};
    //constexpr unsigned n_ids_test = 5;
    //unsigned desired_ids_test[n_ids] = {13, 23, 33, 43, 53};
    // constexpr unsigned n_inlier_pp = 3;
    //unsigned desired_inliers[n_inlier_pp] = {13, 23, 43};//the only inlier that matches with desired id is 13
    vector<uint32_t>desired_inliers_vector; 
    desired_inliers_vector.resize(vec_inliers_.size()) ;
    //using this for loop for get desired_inliers_vector output
    for (unsigned j = 0; j < desired_inliers_vector.size(); j++) {
        //desired_inliers_vector.at(j) = desired_ids[vec_inliers_.at(j)];
        desired_inliers_vector.at(j) = vec_inliers_.at(j);
        //cout << desired_inliers_vector.at(j) <<" " ;
      }
    //unsigned desired_inliers[n_inlier_pp] = {desired_inliers_vector.at(13), 
    //                                         desired_inliers_vector.at(23),
    //                                         desired_inliers_vector.at(63)};//these are the selected inliers from vec_inliers_. Its easier select the result got from robustsolve() than select in robustsolve()
    unsigned tracK__id=0;
    for (const auto &tracK__it: tracK_s_) {
      auto iter = tracK__it.second.cbegin();
      const uint32_t
        i = iter->second,
        j = (++iter)->second,
        K_ = (++iter)->second;
      //
      const auto feature_i = sio_regions_[0]->Features()[i];
      const auto feature_j = sio_regions_[1]->Features()[j];
      const auto feature_K_ = sio_regions_[2]->Features()[K_];
      for (unsigned i=0; i < n_ids; ++i)
        if (tracK__id == desired_ids[i]) { //this part is literaly overwriting the inliers
          //cout<<"blyat"<<endl;
           svg_stream.drawCircle(
              feature_i.x(), feature_i.y(), feature_i.scale(),
              svg::svgStyle().stroK_e("yellow", 1));
           svg_stream.drawCircle(
              feature_j.x(), feature_j.y() + images_[0].Height(), feature_K_.scale(),
              svg::svgStyle().stroK_e("yellow", 1));
           svg_stream.drawCircle(
              feature_K_.x(), feature_K_.y() + images_[0].Height() + images_[1].Height(), feature_K_.scale(),
              svg::svgStyle().stroK_e("yellow", 1));
            //TODO: Tangent line segments in yellow and if inlier -> in green
            svg_stream.drawText(
              feature_i.x()+20, feature_i.y()-20, 6.0f, std::to_string(tracK__id));
     
            svg_stream.drawLine(
              feature_i.x(), feature_i.y(),
              feature_i.x()+20*cos(feature_i.orientation()), feature_i.y() + 20*sin(feature_i.orientation()) ,
              svg::svgStyle().stroK_e("yellow", 1)); 
            svg_stream.drawLine(
              feature_j.x(), feature_j.y() + images_[0].Height(),
              feature_j.x()+20*cos(feature_j.orientation()), feature_j.y() + images_[0].Height()+ 20*sin(feature_j.orientation()),
              svg::svgStyle().stroK_e("yellow", 1));
            svg_stream.drawLine(
              feature_K_.x(), feature_K_.y() + images_[0].Height() + images_[1].Height(),
              feature_K_.x()+ 20*sin(feature_K_.orientation()), feature_K_.y() + images_[0].Height() + images_[1].Height()+ 20*sin(feature_K_.orientation()), //it seems that this last tangent is wrong!!
              svg::svgStyle().stroK_e("yellow", 1));

            svg_stream.drawLine(
              feature_i.x(), feature_i.y(),
              feature_j.x(), feature_j.y() + images_[0].Height(),
              svg::svgStyle().stroK_e("blue", 1));
            svg_stream.drawLine(
              feature_i.x(), feature_i.y(),
              feature_j.x(), feature_j.y() + images_[0].Height(),
              svg::svgStyle().stroK_e("blue", 1));
            svg_stream.drawLine(
              feature_j.x(), feature_j.y() + images_[0].Height(),
              feature_K_.x(), feature_K_.y() + images_[0].Height() + images_[1].Height(),
              svg::svgStyle().stroK_e("blue", 1));
          }
      tracK__id++;
    }
    
    unsigned tracK__inlier = 0;
    for (const auto &tracK__it: tracK_s_)
    {
      auto iter = tracK__it.second.cbegin();
      const uint32_t
        i = iter->second,
        j = (++iter)->second,
        K_ = (++iter)->second;
      //
      const auto feature_i = sio_regions_[0]->Features()[i];
      const auto feature_j = sio_regions_[1]->Features()[j];
      const auto feature_K_ = sio_regions_[2]->Features()[K_];
      for (unsigned i=0; i < desired_inliers_vector.size(); ++i)
        if (tracK__inlier == desired_inliers_vector.at(i)) {
          //cout<<"cyK_a"<<endl; 
          svg_stream.drawCircle(
             feature_i.x(), feature_i.y(), feature_i.scale(),
             svg::svgStyle().stroK_e("green", 1));
          svg_stream.drawCircle(
            feature_j.x(), feature_j.y() + images_[0].Height(), feature_j.scale(),
            svg::svgStyle().stroK_e("green", 1));
          svg_stream.drawCircle(
            feature_K_.x(), feature_K_.y() + images_[0].Height() + images_[1].Height(), feature_K_.scale(),
            svg::svgStyle().stroK_e("green", 1));
          //TODO: Tangent line segments in yellow and if inlier -> in green
          svg_stream.drawText(
            feature_i.x()+20, feature_i.y()-20, 6.0f, std::to_string(tracK__inlier));
     
          svg_stream.drawLine(
            feature_i.x(), feature_i.y(),
            feature_i.x()+20*cos(feature_i.orientation()), feature_i.y() + 20*sin(feature_i.orientation()) ,
            svg::svgStyle().stroK_e("yellow", 1)); 
          svg_stream.drawLine(
            feature_j.x(), feature_j.y() + images_[0].Height(),
            feature_j.x()+20*cos(feature_j.orientation()), feature_j.y() + images_[0].Height()+ 20*sin(feature_j.orientation()),
            svg::svgStyle().stroK_e("yellow", 1));
          svg_stream.drawLine(
            feature_K_.x(), feature_K_.y() + images_[0].Height() + images_[1].Height(),
            feature_K_.x()+ 20*sin(feature_K_.orientation()), feature_K_.y() + images_[0].Height() + images_[1].Height()+ 20*sin(feature_K_.orientation()),
            svg::svgStyle().stroK_e("yellow", 1));

          svg_stream.drawLine(
            feature_i.x(), feature_i.y(),
            feature_j.x(), feature_j.y() + images_[0].Height(),
            svg::svgStyle().stroK_e("lightblue", 1));
          svg_stream.drawLine(
            feature_i.x(), feature_i.y(),
            feature_j.x(), feature_j.y() + images_[0].Height(),
            svg::svgStyle().stroK_e("lightblue", 1));
          svg_stream.drawLine(
            feature_j.x(), feature_j.y() + images_[0].Height(),
            feature_K_.x(), feature_K_.y() + images_[0].Height() + images_[1].Height(),
            svg::svgStyle().stroK_e("lightblue", 1));
        }
      tracK__inlier++;
    }
    ofstream svg_file( "trifocal_tracK_.svg" );
    if (svg_file.is_open())
    {
      svg_file << svg_stream.closeSvgFile().str();
    }
  }

  void DisplayInliersCamerasAndPointsSIFT() {
    // TODO We can then display the inlier and the 3D camera configuration as PLY

    const int svg_w = images_[0].Width();
    const int svg_h = images_[0].Height() + images_[1].Height() + images_[2].Height();
    svg::svgDrawer svg_stream(svg_w, svg_h);

    // Draw image side by side
    svg_stream.drawImage(image_filenames_[0], images_[0].Width(), images_[0].Height());
    svg_stream.drawImage(image_filenames_[1], images_[1].Width(), images_[1].Height(), 0, images_[0].Height());
    svg_stream.drawImage(image_filenames_[2], images_[2].Width(), images_[2].Height(), 0, images_[0].Height() + images_[1].Height());
    
    unsigned tracK__id=0;
    constexpr unsigned n_ids = 5;
    unsigned desired_ids[n_ids] = {13, 23, 33, 63, 53};
    constexpr unsigned n_inlier_pp = 3;
    unsigned desired_inliers[n_inlier_pp] = {13, 23, 43};
    unsigned tracK__inlier=0;
    for (const auto &tracK__it: tracK_s_)
    {
      bool found=false;
      bool inlier=false;
      
      auto iter = tracK__it.second.cbegin();
      const uint32_t
        i = iter->second,
        j = (++iter)->second,
        K_ = (++iter)->second;
      //
      const auto feature_i = sio_regions_[0]->Features()[i];
      const auto feature_j = sio_regions_[1]->Features()[j];
      const auto feature_K_ = sio_regions_[2]->Features()[K_];
      for (unsigned i=0; i < n_ids; ++i)
        if (tracK__id == desired_ids[i]) { //this part is literaly overwriting the inliers
          found = true;
      //cout<<"blyat"<<endl;//using sigma instead is a gives an error in build
      svg_stream.drawCircle(
        feature_i.x(), feature_i.y(), 2*feature_i.scale(),
        svg::svgStyle().stroK_e("yellow", 1));
      svg_stream.drawCircle(
        feature_j.x(), feature_j.y() + images_[0].Height(), feature_K_.scale(),
        svg::svgStyle().stroK_e("yellow", 1));
      svg_stream.drawCircle(
        feature_K_.x(), feature_K_.y() + images_[0].Height() + images_[1].Height(), feature_K_.scale(),
        svg::svgStyle().stroK_e("yellow", 1));
      //TODO: Tangent line segments in yellow and if inlier -> in green
      svg_stream.drawText(
        feature_i.x()+20, feature_i.y()-20, 6.0f, std::to_string(tracK__id));
     
      svg_stream.drawLine(
        feature_i.x(), feature_i.y(),
        feature_i.x()+20*cos(feature_i.orientation()), feature_i.y() + 20*sin(feature_i.orientation()) ,
        svg::svgStyle().stroK_e("yellow", 1)); 
      svg_stream.drawLine(
        feature_j.x(), feature_j.y() + images_[0].Height(),
        feature_j.x()+20*cos(feature_j.orientation()), feature_j.y() + images_[0].Height()+ 20*sin(feature_j.orientation()),
        svg::svgStyle().stroK_e("yellow", 1));
      svg_stream.drawLine(
        feature_K_.x(), feature_K_.y() + images_[0].Height() + images_[1].Height(),
        feature_K_.x()+ 20*sin(feature_K_.orientation()), feature_K_.y() + images_[0].Height() + images_[1].Height()+ 20*sin(feature_K_.orientation()), //it seems that this last tangent is wrong!!
        svg::svgStyle().stroK_e("yellow", 1));

      svg_stream.drawLine(
        feature_i.x(), feature_i.y(),
        feature_j.x(), feature_j.y() + images_[0].Height(),
        svg::svgStyle().stroK_e("blue", 1));
      svg_stream.drawLine(
        feature_i.x(), feature_i.y(),
        feature_j.x(), feature_j.y() + images_[0].Height(),
        svg::svgStyle().stroK_e("blue", 1));
      svg_stream.drawLine(
        feature_j.x(), feature_j.y() + images_[0].Height(),
        feature_K_.x(), feature_K_.y() + images_[0].Height() + images_[1].Height(),
        svg::svgStyle().stroK_e("blue", 1));
      tracK__id++;
          }
      if (!found) {
        tracK__id++;
        continue;
      }
     tracK__id++;
    }
    for (const auto &tracK__it: tracK_s_)
    {
      bool inlier=false;
      
      auto iter = tracK__it.second.cbegin();
      const uint32_t
        i = iter->second,
        j = (++iter)->second,
        K_ = (++iter)->second;
      //
      const auto feature_i = sio_regions_[0]->Features()[i];
      const auto feature_j = sio_regions_[1]->Features()[j];
      const auto feature_K_ = sio_regions_[2]->Features()[K_];
      for (unsigned i=0; i < n_inlier_pp; ++i)
        if (tracK__inlier == desired_inliers[i]){
         //cout<<"cyK_a"<<endl; 
         svg_stream.drawCircle(
            feature_i.x(), feature_i.y(), feature_i.scale(),
            svg::svgStyle().stroK_e("green", 1));
          svg_stream.drawCircle(
            feature_j.x(), feature_j.y() + images_[0].Height(), feature_j.scale(),
            svg::svgStyle().stroK_e("green", 1));
          svg_stream.drawCircle(
            feature_K_.x(), feature_K_.y() + images_[0].Height() + images_[1].Height(), feature_K_.scale(),
            svg::svgStyle().stroK_e("green", 1));
          //TODO: Tangent line segments in yellow and if inlier -> in green
          svg_stream.drawText(
            feature_i.x()+20, feature_i.y()-20, 6.0f, std::to_string(tracK__inlier));
     
      svg_stream.drawLine(
        feature_i.x(), feature_i.y(),
        feature_i.x()+20*cos(feature_i.orientation()), feature_i.y() + 20*sin(feature_i.orientation()) ,
        svg::svgStyle().stroK_e("yellow", 1)); 
      svg_stream.drawLine(
        feature_j.x(), feature_j.y() + images_[0].Height(),
        feature_j.x()+20*cos(feature_j.orientation()), feature_j.y() + images_[0].Height()+ 20*sin(feature_j.orientation()),
        svg::svgStyle().stroK_e("yellow", 1));
      svg_stream.drawLine(
        feature_K_.x(), feature_K_.y() + images_[0].Height() + images_[1].Height(),
        feature_K_.x()+ 20*sin(feature_K_.orientation()), feature_K_.y() + images_[0].Height() + images_[1].Height()+ 20*sin(feature_K_.orientation()),
        svg::svgStyle().stroK_e("yellow", 1));

      svg_stream.drawLine(
        feature_i.x(), feature_i.y(),
        feature_j.x(), feature_j.y() + images_[0].Height(),
        svg::svgStyle().stroK_e("lightblue", 1));
      svg_stream.drawLine(
        feature_i.x(), feature_i.y(),
        feature_j.x(), feature_j.y() + images_[0].Height(),
        svg::svgStyle().stroK_e("lightblue", 1));
      svg_stream.drawLine(
        feature_j.x(), feature_j.y() + images_[0].Height(),
        feature_K_.x(), feature_K_.y() + images_[0].Height() + images_[1].Height(),
        svg::svgStyle().stroK_e("lightblue", 1));
          inlier = true;
        }  
      if (!inlier) {
        tracK__inlier++;
        continue;
      }
     tracK__inlier++;
    }
    ofstream svg_file( "trifocal_tracK__SIFT.svg" );
    if (svg_file.is_open())
    {
      svg_file << svg_stream.closeSvgFile().str();
  // double R0[3][3] = {
  //                    {1,0,0},
  //                    {0,1,0},
  //                    {0,0,1}
  //                   };
  //   
  //  double T0[3][1] = {
  //                      {13.022176},
  //                      {-1.6546488},
  //                      {352.47945}
  //                     };

  //  double R1[3][3] = {   
  //                     {9.4083600000000001e-01,   2.7502399999999999e-01,   1.9796300000000000e-01},
  //                     {2.9655799999999999e-01,  -3.8560499999999998e-01,  -8.7370599999999998e-01},
  //                     {-1.6395499999999999e-01,   8.8072200000000000e-01,  -4.4435100000000000e-01}
  //                    };
  //  double T1[3][1] = { 
  //                      {8.70420714},
  //                      {-1.62157456},
  //                      {-352.61248141}
  //                    };

  //  double R2[3][3] = {
  //                     {0.970091125581631,   0.235130101826381,   0.060307903987350},
  //                     {0.151694164781553,  -0.393265050435905,  -0.906824944780907},
  //                     {-0.189504850701909,   0.888851188512892,  -0.417170799840706}
  //                    };

  //  double T2[3][1] = { 
  //                      {0.88920328},
  //                      {-14.05063273},
  //                      {-352.44248798}
  //                    };
  
 //fill C0* with for loop
  std::cerr << "Number of sols " << nsols_final << std::endl;
  std::vector<trifocal_model_t> &tt = *trifocal_tensor; // if I use the STL container, 
  // This I would have to change the some other pieces of code, maybe altering the entire logic of this program!!
  // std::cerr << "TRIFOCAL LOG: Antes de resize()\n" << std::endl;
  tt.resize(nsols_final);
  std::cerr << "TRIFOCAL LOG: Chamou resize()\n";
  //using trifocal_model_t = array<Mat34, 3>;
  for (unsigned s=0; s < nsols_final; ++s) {
    tt[s][0] = Mat34::Identity(); // view 0 [I | 0]
    for (unsigned v=1; v < io::pp::nviews; ++v) {
        memcpy(tt[s][v].data(), (double *) cameras[id_sols[s]][v], 9*sizeof(double));
        for (unsigned r=0; r < 3; ++r)
          tt[s][v](r,3) = cameras[id_sols[s]][v][3][r];
    }
  }
  //This is for hard coding test 
  //tt[0][0] = Mat34::Identity();
  //tt[0][1] = Mat34::Identity();
  //tt[0][2] = Mat34::Identity();
  //for(unsigned i=0;i<3;i++){
  //  for(unsigned j=0;j<4;j++){
  //    if(j<3){
  //      tt[0][0](i,j) = R0[i][j];
  //      tt[0][1](i,j) = R1[i][j];
  //      tt[0][2](i,j) = R2[i][j];
  //    }
  //    else{
  //      tt[0][0](i,j) = T0[i][1];
  //      tt[0][1](i,j) = T1[i][1];
  //      tt[0][2](i,j) = T2[i][1];
  //    }
  //  }                       
  //}
  // cout << "this is [R0|T0] " << "\n"; cout << tt[0][0] << "\n";
  // cout << "this is [R1|T1] " << "\n"; cout << tt[0][1] << "\n";
  // cout << "this is [R2|T2] " << "\n"; cout << tt[0][2] << "\n";
  
  // TODO: filter the solutions by:
  // - positive depth and 
  // - using tangent at 3rd point
  //
  //  if we K_now the rays are perfectly coplanar, we can just use cross
  // product within the plane instead of SVD
  std::cerr << "TRIFOCAL LOG: Finished ()Solve()\n";
}

void Trifocal3PointPositionTangentialSolver::
static double Error(
  const trifocal_model_t &tt,
  const Vec &bearing_0, // x,y,tangentialx,tangentialy
  const Vec &bearing_1,
  const Vec &bearing_2,
  const Vec &pixbearing_0,
  const Vec &pixbearing_1,
  const Vec &pixbearing_2,
  const double K__[2][3]) 
{
  //std::cerr << "TRIFOCAL LOG: Called Error()\n";
  // Return the cost related to this model and those sample data point
  // Ideal algorithm:
  // 1) reconstruct the 3D points and orientations
  // 2) project the 3D points and orientations on all images_
  // 3) compute error 
  // 
  // In practice we ignore the directions and only reproject to one third view
  
  // 3x3: each column is x,y,1
  Mat3 bearing;
  bearing << bearing_0.head(2).homogeneous(),
             bearing_1.head(2).homogeneous(), 
             bearing_2.head(2).homogeneous();
  Mat2 pixbearing; // << XXX mat2
  pixbearing << pixbearing_0.head(2).homogeneous(),
                pixbearing_1.head(2).homogeneous();
  // Using triangulation.hpp
  Vec4 triangulated_homg;
  unsigned third_view = 0;
  // picK_ the wider baseline. TODO: measure all pairwise translation distances
  if (tt[1].col(3).squaredNorm() > tt[2].col(3).squaredNorm()) {
    // TODO use triangulation from the three views at once
    TriangulateDLT(tt[0], bearing.col(0), tt[1], bearing.col(1), &triangulated_homg);
    third_view = 2;
  } else {
    TriangulateDLT(tt[0], bearing.col(0), tt[2], bearing.col(2), &triangulated_homg);
    third_view = 1;
  }
  
  // Computing the projection of triangulated points using projection.hpp
  // For prototyping and speed, for now we will only project to the third view
  // and report only one error
  Vec2 pxreprojected = Vec3(tt[third_view]*triangulated_homg).hnormalized();
  // XXX revert intrinsics to measure the error in pixels
  revert_intrinsics(K_, pxreprojected, pxreprojected);
   
  Vec2 pxmeasured    = pxbearing.col(third_view);
  //cout << "error " << (reprojected - measured).squaredNorm() << "\n";
  //cout << "triang " <<triangulated_homg <<"\n";
  //std::cerr << "TRIFOCAL LOG: Finished Error()\n";
  return (reprojected-measured).squaredNorm();
}
