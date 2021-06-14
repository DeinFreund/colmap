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

#ifndef COLMAP_SRC_BASE_LINETRACK_H_
#define COLMAP_SRC_BASE_LINETRACK_H_

#include <vector>

#include "util/logging.h"
#include "util/types.h"

namespace colmap {

// LineTrack class stores all observations of a 3D line.
struct LineTrackElement {
  LineTrackElement();
  LineTrackElement(const image_t image_id, const line2D_t line2D_idx,
                   double start_parameter, double end_parameter,
                   bool fixed_start, bool fixed_end);
  // The image in which the lineTrack element is observed.
  image_t image_id;
  // The line in the image that the lineTrack element is observed.
  line2D_t line2D_idx;

  double start_parameter, end_parameter;
  bool fixed_start, fixed_end;

  bool operator==(const LineTrackElement& o) const {
    return o.image_id == image_id && o.line2D_idx == line2D_idx;
  }
  bool operator!=(const LineTrackElement& o) const { return !(o == *this); }
};

class LineTrack {
 public:
  LineTrack();

  // The number of lineTrack elements.
  inline size_t Length() const;

  // Access all elements.
  inline const std::vector<LineTrackElement>& Elements() const;
  inline std::vector<LineTrackElement>& Elements();
  inline void SetElements(const std::vector<LineTrackElement>& elements);

  // Access specific elements.
  inline const LineTrackElement& Element(const size_t idx) const;
  inline LineTrackElement& Element(const size_t idx);
  inline void SetElement(const size_t idx, const LineTrackElement& element);

  // Append new elements.
  inline void AddElement(LineTrackElement element);

  template <class... Args>
  inline void AddElement(Args&&... args) {
    AddElement(LineTrackElement(std::forward<Args>(args)...));
  }

  inline void AddElements(const std::vector<LineTrackElement>& elements);

  // Delete existing element.
  inline void DeleteElement(const size_t idx);
  void DeleteElement(const LineTrackElement& el);
  void DeleteElement(const image_t image_id, const line2D_t line2D_idx) {
    DeleteElement(LineTrackElement(image_id, line2D_idx, 0, 0, false, false));
  }

  // Requests that the lineTrack capacity be at least enough to contain the
  // specified number of elements.
  inline void Reserve(const size_t num_elements);

  // Shrink the capacity of lineTrack vector to fit its size to save memory.
  inline void Compress();

 private:
  std::vector<LineTrackElement> elements_;
};

////////////////////////////////////////////////////////////////////////////////
// Implementation
////////////////////////////////////////////////////////////////////////////////

size_t LineTrack::Length() const { return elements_.size(); }

const std::vector<LineTrackElement>& LineTrack::Elements() const {
  return elements_;
}

std::vector<LineTrackElement>& LineTrack::Elements() { return elements_; }

void LineTrack::SetElements(const std::vector<LineTrackElement>& elements) {
  elements_ = elements;
}

// Access specific elements.
const LineTrackElement& LineTrack::Element(const size_t idx) const {
  return elements_.at(idx);
}

LineTrackElement& LineTrack::Element(const size_t idx) {
  return elements_.at(idx);
}

void LineTrack::SetElement(const size_t idx, const LineTrackElement& element) {
  elements_.at(idx) = element;
}

void LineTrack::AddElement(LineTrackElement element) {
  elements_.emplace_back(std::move(element));
}

void LineTrack::AddElements(const std::vector<LineTrackElement>& elements) {
  elements_.insert(elements_.end(), elements.begin(), elements.end());
}

void LineTrack::DeleteElement(const size_t idx) {
  CHECK_LT(idx, elements_.size());
  elements_.erase(elements_.begin() + idx);
}

void LineTrack::Reserve(const size_t num_elements) {
  elements_.reserve(num_elements);
}

void LineTrack::Compress() { elements_.shrink_to_fit(); }

}  // namespace colmap

#endif  // COLMAP_SRC_BASE_LINETRACK_H_
