#ifndef COLMAP_SRC_ESTIMATORS_LINE_H_
#define COLMAP_SRC_ESTIMATORS_LINE_H_

#include "base/camera.h"
#include "base/image.h"
#include "base/line2d.h"
#include "base/line3d.h"

namespace colmap {

Line3D EstimateLine3D(const Camera& cam1, const Image& img1,
                      const Line2D& line1, const Camera& cam2,
                      const Image& img2, const Line2D& line2);

std::vector<Line3D> EstimateLines(const Camera& cam1, const Image& img1,
                                  const Camera& cam2, const Image& img2,
                                  const Camera& test_camera,
                                  const Image& test_image);

point2D_t MatchLine(const Camera& cam, const Image& img, const Line3D& line);

}  // namespace colmap



#endif
