#include "estimators/line.h"

#include "base/camera_models.h"
#include "base/cost_functions.h"
#include "base/pose.h"
#include "base/projection.h"

namespace {

using namespace colmap;  // needed for macro

/**
 * Given a point2d, calculate its location on the normalized projection plane in
 * world space
 */
template <typename CameraModel>
Eigen::Vector3d ImagePointToWorld(const double* const params, const Image& img,
                                  const Eigen::Vector2d& point2D) {
  Eigen::Vector3d point3D(0, 0, 1);
  CameraModel::ImageToWorld(params, point2D.x(), point2D.y(), &point3D.x(),
                            &point3D.y());
  point3D -= img.Tvec();
  point3D = QuaternionRotatePoint(InvertQuaternion(img.Qvec()), point3D);
  return point3D;
}

Eigen::Vector3d ImagePointToWorld(const Camera& cam, const Image& img,
                                  const Eigen::Vector2d& point2D) {
  switch (cam.ModelId()) {
#define CAMERA_MODEL_CASE(CameraModel) \
  case CameraModel::kModelId:          \
    return ImagePointToWorld<CameraModel>(cam.ParamsData(), img, point2D);

    CAMERA_MODEL_SWITCH_CASES

#undef CAMERA_MODEL_CASE
  }
}

template <typename CameraModel>
Eigen::Vector2d WorldPointToImage(const double* const params, const Image& img,
                                  Eigen::Vector3d point3D) {
  Eigen::Vector2d point2D;
  point3D = QuaternionRotatePoint(img.Qvec(), point3D);
  point3D += img.Tvec();
  CameraModel::WorldToImage(params, point3D.x(), point3D.y(), &point2D.x(),
                            &point2D.y());
  return point2D;
}

Eigen::Vector2d WorldPointToImage(const Camera& cam, const Image& img,
                                  const Eigen::Vector3d& point3D) {
  switch (cam.ModelId()) {
#define CAMERA_MODEL_CASE(CameraModel) \
  case CameraModel::kModelId:          \
    return WorldPointToImage<CameraModel>(cam.ParamsData(), img, point3D);

    CAMERA_MODEL_SWITCH_CASES

#undef CAMERA_MODEL_CASE
  }
}

/**
 * Calculates intersection point of a plane with a line
 */
Eigen::Vector3d IntersectLineWithPlane(const Eigen::Vector3d& line_point,
                                       const Eigen::Vector3d& line_direction,
                                       const Eigen::Vector3d& plane_point,
                                       const Eigen::Vector3d& plane_normal) {
  const double distance = plane_normal.dot(line_point - plane_point);
  return line_point - line_direction.normalized() * distance /
                          plane_normal.dot(line_direction.normalized());
}

/**
 * Returns the closest point to the target line on the source line.
 */
Eigen::Vector3d ClosestPointToLineOnLine(
    const Eigen::Vector3d& source_line_point,
    const Eigen::Vector3d& source_line_direction,
    const Eigen::Vector3d& target_line_point,
    const Eigen::Vector3d& target_line_direction) {
  const Eigen::Vector3d normal =
      target_line_direction
          .cross(target_line_direction.cross(source_line_direction))
          .normalized();
  return IntersectLineWithPlane(source_line_point, source_line_direction,
                                target_line_point, normal);
}

}  // namespace

namespace colmap {

// Projects the two endpoints of the 2d line onto the 3d line
// Returns their position as distance along the line
std::pair<double, double> ProjectOnLine(const Camera& cam, const Image& img, const Line2D& line2d, const Line3D& line3d) {
    const Eigen::Vector3d line_dir = line3d.XYZ2() - line3d.XYZ1();
    const Eigen::Vector3d img_pt1 = ImagePointToWorld(cam, img, line2d.XY1());
    const Eigen::Vector3d img_pt2 = ImagePointToWorld(cam, img, line2d.XY2());
    const Eigen::Vector3d img_dir1 = img_pt1 - img.ProjectionCenter();
    const Eigen::Vector3d img_dir2 = img_pt2 - img.ProjectionCenter();
    const Eigen::Vector3d line_pt1 = ClosestPointToLineOnLine(line3d.XYZ1(), line_dir, img_pt1, img_dir1);
    const Eigen::Vector3d line_pt2 = ClosestPointToLineOnLine(line3d.XYZ1(), line_dir, img_pt2, img_dir2);
    std::pair<double, double> res{(line_pt1 - line3d.XYZ1()).dot(line_dir.normalized()), (line_pt2 - line3d.XYZ1()).dot(line_dir.normalized())};
    if (res.first > res.second) {
      std::swap(res.first, res.second);
    }
    return res;
}

bool CheckLineOverlap(const Camera& cam1, const Image& img1, const Line2D& line2d1, const Camera& cam2, const Image& img2, const Line2D& line2d2, const Line3D& line3d) {
  const double min_line_overlap = 0.7;

  const std::pair<double, double> line_points =
      ProjectOnLine(cam1, img1, line2d1, line3d);
  const std::pair<double, double> track_points =
      ProjectOnLine(cam2, img2, line2d2, line3d);
  const double overlap = (std::min(line_points.second, track_points.second) -
                          std::max(line_points.first, track_points.first)) /
                         std::max(line_points.second - line_points.first,
                                  track_points.second - track_points.first);
  //std::cerr << "got1 " << line_points.first << " -> " << line_points.second << "\n";
  //std::cerr << "got2 " << track_points.first << " -> " << track_points.second << "\n";
  //std::cerr <<"overlap is " << overlap << "\n";
  CHECK_LT(overlap, 1.0);
  return overlap > min_line_overlap;
}

bool CheckLineOverlap(const Camera& cam, const Image& img, const Line2D& line2d, const Line3D& line3d, const Reconstruction& reconstruction) {

    for (const auto& track_el : line3d.Track().Elements()) {
        if (track_el.image_id == img.ImageId()) continue;
      const class Image& image = reconstruction.Image(track_el.image_id);
      const class Camera& camera = reconstruction.Camera(image.CameraId());
      const Line2D& track_line = image.Line2D(track_el.point2D_idx);
      if (CheckLineOverlap(cam, img, line2d, camera, image, track_line, line3d)) {
          return true;
      }
    }

    return false;
}

Eigen::Vector2d LineReprojectionCost(const Camera& cam, const Image& img,
                                     const Line2D& line2D,
                                     const Line3D& line3D) {
  Eigen::Vector2d result;
  switch (cam.ModelId()) {
#define CAMERA_MODEL_CASE(CameraModel)                                        \
  case CameraModel::kModelId:                                                 \
    LineBundleAdjustmentCostFunction<CameraModel>(line2D.XY1(), line2D.XY2()) \
        .                                                                     \
        operator()<double>(img.Qvec().data(), img.Tvec().data(),              \
                           line3D.XYZ1().data(), line3D.XYZ2().data(),        \
                           cam.ParamsData(), result.data());                  \
    break;

    CAMERA_MODEL_SWITCH_CASES

#undef CAMERA_MODEL_CASE
  }
  return result;
}

Line3D EstimateLine3D(const Camera& cam1, const Image& img1,
                      const Line2D& line1, const Camera& cam2,
                      const Image& img2, const Line2D& line2) {
  const Eigen::Vector3d img1_pt1 = ImagePointToWorld(cam1, img1, line1.XY1());
  const Eigen::Vector3d img1_pt2 = ImagePointToWorld(cam1, img1, line1.XY2());
  const Eigen::Vector3d img2_pt1 = ImagePointToWorld(cam2, img2, line2.XY1());
  const Eigen::Vector3d img2_pt2 = ImagePointToWorld(cam2, img2, line2.XY2());
  const Eigen::Vector3d img1_dir1 = img1_pt1 - img1.ProjectionCenter();
  const Eigen::Vector3d img1_dir2 = img1_pt2 - img1.ProjectionCenter();
  const Eigen::Vector3d img2_dir1 = img2_pt1 - img2.ProjectionCenter();
  const Eigen::Vector3d img2_dir2 = img2_pt2 - img2.ProjectionCenter();
  const Eigen::Vector3d plane_normal1 = img1_dir1.cross(img1_dir2).normalized();
  const Eigen::Vector3d plane_normal2 = img2_dir1.cross(img2_dir2).normalized();

  if (plane_normal1.dot(plane_normal2) > 0.999) {
    // Ill-formed problem, planes are nearly parallel
    return {};
  }

  // Find extremal points
  std::array<Eigen::Vector3d, 4> line3D_pts;
  line3D_pts[0] =
      IntersectLineWithPlane(img1.ProjectionCenter(), img1_dir1,
                             img2.ProjectionCenter(), plane_normal2);
  line3D_pts[1] =
      IntersectLineWithPlane(img1.ProjectionCenter(), img1_dir2,
                             img2.ProjectionCenter(), plane_normal2);
  line3D_pts[2] =
      IntersectLineWithPlane(img2.ProjectionCenter(), img2_dir1,
                             img1.ProjectionCenter(), plane_normal1);
  line3D_pts[3] =
      IntersectLineWithPlane(img2.ProjectionCenter(), img2_dir2,
                             img1.ProjectionCenter(), plane_normal1);
  for (size_t i = 0; i < 4; i++) {
    if (!HasPointPositiveDepth(img1.ProjectionMatrix(), line3D_pts[i])) {
      return {};
    }
    if (!HasPointPositiveDepth(img2.ProjectionMatrix(), line3D_pts[i])) {
      return {};
    }
  }
  const Eigen::Vector3d line3D_dir =
      plane_normal1.cross(plane_normal2).normalized();

  std::sort(line3D_pts.begin(), line3D_pts.end(),
            [&](const auto& l, const auto& r) {
              return line3D_dir.dot(l - line3D_pts[0]) <
                     line3D_dir.dot(r - line3D_pts[0]);
            });
  Line3D line3d(std::move(line3D_pts[0]), std::move(line3D_pts[3]));
  if (!CheckLineOverlap(cam1, img1, line1, cam2, img2, line2, line3d)) {
      return {};
  }
  return line3d;
}

void RecalculateEndpoints(const Reconstruction& reconstruction,
                          Line3D* const linePtr) {
  Line3D& line = *linePtr;
  CHECK(line.Track().Length() >= 2);

  const Eigen::Vector3d line_dir = line.XYZ2() - line.XYZ1();
  std::vector<Eigen::Vector3d> points3d;
  for (const auto& track_el : line.Track().Elements()) {
    const Image& img = reconstruction.Image(track_el.image_id);
    const Camera& cam = reconstruction.Camera(img.CameraId());
    const Line2D& line2d = img.Line2D(track_el.point2D_idx);
    const Eigen::Vector3d img_pt1 = ImagePointToWorld(cam, img, line2d.XY1());
    const Eigen::Vector3d img_pt2 = ImagePointToWorld(cam, img, line2d.XY2());
    const Eigen::Vector3d img_dir1 = img_pt1 - img.ProjectionCenter();
    const Eigen::Vector3d img_dir2 = img_pt2 - img.ProjectionCenter();
    points3d.push_back(ClosestPointToLineOnLine(line.XYZ1(), line_dir, img_pt1, img_dir1));
    points3d.push_back(ClosestPointToLineOnLine(line.XYZ1(), line_dir, img_pt2, img_dir2));
  }

  std::sort(
      points3d.begin(), points3d.end(), [&](const auto& l, const auto& r) {
        return line_dir.dot(l - line.XYZ1()) < line_dir.dot(r - line.XYZ1());
      });
/*  std::cerr << "Updated line from \n"
            << line.XYZ1() << " -> \n"
            << line.XYZ2() << " to \n"
            << points3d[0] << " -> \n"
            << points3d.back() << "\n";*/
  line.SetXYZ(points3d[0], points3d.back());
}

line2D_t MatchLine(const Camera& cam, const Image& img, const Line3D& line3D, const Reconstruction& reconstruction) {
  const double max_reproj_err = 4;
  double min_reproj_err = std::numeric_limits<double>::max();
  line2D_t best_line2D_idx = kInvalidLine2DIdx;

  for (line2D_t line2D_idx = 0; line2D_idx < img.NumLines2D(); ++line2D_idx) {
    const Line2D& test_line = img.Line2D(line2D_idx);

    if (!HasPointPositiveDepth(img.ProjectionMatrix(), line3D.XYZ1())) {
      continue;
    }
    if (!HasPointPositiveDepth(img.ProjectionMatrix(), line3D.XYZ2())) {
      continue;
    }
    if (test_line.HasLine3D()) {
      continue;  // alternatively check if new match would be better
    }
    if (!CheckLineOverlap(cam, img, test_line, line3D, reconstruction)) {
        continue;
    }
    Eigen::Vector2d reprojErr =
        LineReprojectionCost(cam, img, test_line, line3D);
    if (reprojErr.x() < max_reproj_err && reprojErr.y() < max_reproj_err &&
        reprojErr.norm() < min_reproj_err) {
      min_reproj_err = reprojErr.norm();
      best_line2D_idx = line2D_idx;
    }
  }
  return best_line2D_idx;
}

std::vector<Line3D> EstimateLines(const Camera& cam1, const Image& img1,
                                  const Camera& cam2, const Image& img2,
                                  const Camera& test_camera,
                                  const Image& test_image) {
  std::cerr << "estimating lines " << img1.NumLines2D() << " "
            << img2.NumLines2D() << " " << test_image.NumLines2D() << "\n";
  std::list<Line3D> result;

  // Observed lines for each 2D line, since each 2D line can only be associated
  // to one 3d line, this is used to select the best candidate. There's one map
  // of candidates for each of the three images
  std::vector<
      std::map<line2D_t, std::map<double, std::reference_wrapper<Line3D>>>>
      trackCandidates(3);
  for (line2D_t line2D_idx1 = 0; line2D_idx1 < img1.NumLines2D();
       ++line2D_idx1) {
    for (line2D_t line2D_idx2 = 0; line2D_idx2 < img2.NumLines2D();
         ++line2D_idx2) {
      const Line2D& line1 = img1.Line2D(line2D_idx1);
      const Line2D& line2 = img2.Line2D(line2D_idx2);

      Line3D line3D = EstimateLine3D(cam1, img1, line1, cam2, img2, line2);
      if (line3D.Length() < 1e-6) {
        continue;
      }

      const double max_reproj_err = 4;

      for (line2D_t line2D_idx3 = 0; line2D_idx3 < test_image.NumLines2D();
           ++line2D_idx3) {
        const Line2D& test_line = test_image.Line2D(line2D_idx3);
        if (!CheckLineOverlap(test_camera, test_image, test_line, cam1, img1, line1, line3D) && !CheckLineOverlap(test_camera, test_image, test_line, cam2, img2, line2, line3D)) {
            continue;
        }
        Eigen::Vector2d reprojErr1 =
            LineReprojectionCost(cam1, img1, line1, line3D);
        CHECK_LT(reprojErr1.norm(), 1e-4);
        Eigen::Vector2d reprojErr2 =
            LineReprojectionCost(cam2, img2, line2, line3D);
        CHECK_LT(reprojErr2.norm(), 1e-4);
        Eigen::Vector2d reprojErr =
            LineReprojectionCost(test_camera, test_image, test_line, line3D);
        // std::cerr << "err: \n" << reprojErr << "\n";
        if (reprojErr.x() < max_reproj_err && reprojErr.y() < max_reproj_err) {
          result.push_back(std::move(line3D));
          result.back().SetError(reprojErr.norm());
          result.back().Track().AddElement(img1.ImageId(), line2D_idx1);
          result.back().Track().AddElement(img2.ImageId(), line2D_idx2);
          result.back().Track().AddElement(test_image.ImageId(), line2D_idx3);
          /*trackCandidates[0][line2D_idx1].emplace(reprojErr.norm(),
                                                  result.back());
          trackCandidates[1][line2D_idx2].emplace(reprojErr.norm(),
                                                  result.back());
          trackCandidates[2][line2D_idx3].emplace(reprojErr.norm(),
          result.back());*/
          break;
        }
      }
    }
  }
  /*
  for (const auto& candidate : trackCandidates[0]) {
    if (img1.Line2D(candidate.first).HasLine3D())
      continue;  // alternatively check if new match would be better
    candidate.second.begin()->second.get().Track().AddElement(img1.ImageId(),
                                                              candidate.first);
  }
  for (const auto& candidate : trackCandidates[1]) {
    if (img2.Line2D(candidate.first).HasLine3D())
      continue;  // alternatively check if new match would be better
    candidate.second.begin()->second.get().Track().AddElement(img2.ImageId(),
                                                              candidate.first);
  }
  for (const auto& candidate : trackCandidates[2]) {
    if (test_image.Line2D(candidate.first).HasLine3D())
      continue;  // alternatively check if new match would be better
    candidate.second.begin()->second.get().Track().AddElement(
        test_image.ImageId(), candidate.first);
  }
  */
  std::vector<Line3D> lines3D;
  const auto isConstrained = [](const Line3D& line) {
    return line.Track().Length() > 1;
  };
  std::copy_if(std::make_move_iterator(result.begin()),
               std::make_move_iterator(result.end()),
               std::back_inserter(lines3D), isConstrained);
  return lines3D;
}

}  // namespace colmap
