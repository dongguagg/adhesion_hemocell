/*
This file is part of the HemoCell library

HemoCell is developed and maintained by the Computational Science Lab
in the University of Amsterdam. Any questions or remarks regarding this library
can be sent to: info@hemocell.eu

When using the HemoCell library in scientific work please cite the
corresponding paper: https://doi.org/10.3389/fphys.2017.00563

The HemoCell library is free software: you can redistribute it and/or
modify it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.
*/
#ifndef HEMOCELLSURFACEFORCE_H
#define HEMOCELLSURFACEFORCE_H

#include "hemocell.h"

#include <map>
#include <string>
#include <vector>

namespace hemo {

enum class CellSurfaceFace {
  XPositive,
  XNegative,
  YPositive,
  YNegative,
  ZPositive,
  ZNegative
};

enum class SurfacePatchMode { GeodesicRadius, OneRing };

enum class SurfaceForceWeighting { NodalArea, Uniform };

struct SurfaceForceVertex {
  plint vertexId = -1;
  T weight = 0.0;
  hemo::Array<T, 3> initialPosition = {0.0, 0.0, 0.0};
  T alignment = 0.0;
  T geodesicDistance = 0.0;
  T nodalArea = 0.0;
};

CellSurfaceFace parseCellSurfaceFace(const std::string &value);
SurfacePatchMode parseSurfacePatchMode(const std::string &value);
SurfaceForceWeighting parseSurfaceForceWeighting(const std::string &value);

const char *cellSurfaceFaceName(CellSurfaceFace value);
const char *surfacePatchModeName(SurfacePatchMode value);
const char *surfaceForceWeightingName(SurfaceForceWeighting value);

/**
 * Apply a configured total force to a fixed material patch of a cell membrane.
 *
 * The patch is selected once from the initial geometry. A face-center seed is
 * found from the direction of each vertex relative to the cell center, after
 * which the patch is expanded over the triangle graph. The force is added to
 * sv.force on owner particles and therefore remains fully coupled to the fluid.
 */
class HemoCellSurfaceForce {
public:
  HemoCellSurfaceForce(
      HemoCell &hemocell, const std::string &cellType,
      plint targetBaseCellId, CellSurfaceFace face,
      const hemo::Array<T, 3> &forcePicoNewton,
      SurfacePatchMode patchMode, T patchGeodesicRadiusMicrometer,
      SurfaceForceWeighting weighting);

  void initializeFromCurrentGeometry();
  void restoreGripState(const std::string &filename);
  void saveGripState(const std::string &filename) const;
  void saveSelectedVerticesCsv(const std::string &filename) const;
  void applyForce(T loadScale);

  const std::vector<SurfaceForceVertex> &selectedVertices() const;
  hemo::Array<T, 3> forcePicoNewton() const;
  hemo::Array<T, 3> forceLbm() const;
  hemo::Array<T, 3> initialCellCenter() const;
  hemo::Array<T, 3> currentWeightedPosition() const;
  plint seedVertexId() const;
  plint targetBaseCellId() const;
  CellSurfaceFace face() const;
  T sumWeights() const;

private:
  struct VertexGeometry {
    hemo::Array<T, 3> position = {0.0, 0.0, 0.0};
  };

  class GatherTargetGeometry;
  class AddSurfaceForceFunctional;

  std::map<int, VertexGeometry> gatherTargetGeometry() const;
  void validateTimescales() const;
  void validateSelectedState(bool validateCurrentTarget) const;

  HemoCell &hemocell_;
  HemoCellField &cellField_;
  std::string cellType_;
  plint targetBaseCellId_;
  CellSurfaceFace face_;
  hemo::Array<T, 3> forcePicoNewton_;
  hemo::Array<T, 3> forceLbm_;
  SurfacePatchMode patchMode_;
  T patchGeodesicRadiusMicrometer_;
  SurfaceForceWeighting weighting_;
  plint seedVertexId_ = -1;
  hemo::Array<T, 3> initialCellCenter_ = {0.0, 0.0, 0.0};
  std::vector<SurfaceForceVertex> selectedVertices_;
  std::vector<T> vertexWeights_;
  bool initialized_ = false;
};

} // namespace hemo

#endif
