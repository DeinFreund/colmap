// Copyright (c) 2018, ETH Zurich and UNC Chapel Hill.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//
//     * Neither the name of ETH Zurich and UNC Chapel Hill nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// Author: Johannes L. Schoenberger (jsch-at-demuc-dot-de)

#include "base/line_track.h"

namespace colmap {

LineTrack::LineTrack() {}

LineTrackElement::LineTrackElement()
    : image_id(kInvalidImageId),
      line2D_idx(kInvalidLine2DIdx),
      start_parameter(std::numeric_limits<double>::quiet_NaN()),
      end_parameter(std::numeric_limits<double>::quiet_NaN()),
      fixed_start(false),
      fixed_end(false) {}

LineTrackElement::LineTrackElement(const image_t image_id_,
                                   const line2D_t line2D_idx_,
                                   double start_parameter_,
                                   double end_parameter_, bool fixed_start_,
                                   bool fixed_end_)
    : image_id(image_id_),
      line2D_idx(line2D_idx_),
      start_parameter(start_parameter_),
      end_parameter(end_parameter_),
      fixed_start(fixed_start_),
      fixed_end(fixed_end_) {}

void LineTrack::DeleteElement(const LineTrackElement& el) {
  elements_.erase(std::remove_if(elements_.begin(), elements_.end(),
                                 [&el](const LineTrackElement& element) {
                                   return element == el;
                                 }),
                  elements_.end());
}

}  // namespace colmap
