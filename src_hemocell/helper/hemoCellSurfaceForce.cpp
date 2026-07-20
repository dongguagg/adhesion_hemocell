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
#include "hemoCellSurfaceForce.h"

#include "hemoCellFunctional.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <queue>
#include <set>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace hemo {

namespace {

hemo::Array<T, 3> faceDirection(CellSurfaceFace face) {
  switch (face) {
  case CellSurfaceFace::XPositive:
    return {1.0, 0.0, 0.0};
  case CellSurfaceFace::XNegative:
    return {-1.0, 0.0, 0.0};
  case CellSurfaceFace::YPositive:
    return {0.0, 1.0, 0.0};
  case CellSurfaceFace::YNegative:
    return {0.0, -1.0, 0.0};
  case CellSurfaceFace::ZPositive:
    return {0.0, 0.0, 1.0};
  case CellSurfaceFace::ZNegative:
    return {0.0, 0.0, -1.0};
  }
  throw std::invalid_argument("Unknown cell surface face enum value");
}

T dot(const hemo::Array<T, 3> &a, const hemo::Array<T, 3> &b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

T norm(const hemo::Array<T, 3> &a) { return std::sqrt(dot(a, a)); }

hemo::Array<T, 3> subtract(const hemo::Array<T, 3> &a,
                           const hemo::Array<T, 3> &b) {
  return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}

hemo::Array<T, 3> cross(const hemo::Array<T, 3> &a,
                        const hemo::Array<T, 3> &b) {
  return {a[1] * b[2] - a[2] * b[1],
          a[2] * b[0] - a[0] * b[2],
          a[0] * b[1] - a[1] * b[0]};
}

bool finiteArray(const hemo::Array<T, 3> &value) {
  return std::isfinite(value[0]) && std::isfinite(value[1]) &&
         std::isfinite(value[2]);
}

bool nearlyEqual(T a, T b, T tolerance = 1.0e-12) {
  const T scale = std::max<T>(1.0, std::max(std::fabs(a), std::fabs(b)));
  return std::fabs(a - b) <= tolerance * scale;
}

void requireKey(std::istream &stream, const char *expected) {
  std::string key;
  if (!(stream >> key) || key != expected) {
    throw std::runtime_error(std::string("Invalid surface-force state: expected '") +
                             expected + "'");
  }
}

} // namespace

CellSurfaceFace parseCellSurfaceFace(const std::string &value) {
  if (value == "xPositive") {
    return CellSurfaceFace::XPositive;
  }
  if (value == "xNegative") {
    return CellSurfaceFace::XNegative;
  }
  if (value == "yPositive") {
    return CellSurfaceFace::YPositive;
  }
  if (value == "yNegative") {
    return CellSurfaceFace::YNegative;
  }
  if (value == "zPositive") {
    return CellSurfaceFace::ZPositive;
  }
  if (value == "zNegative") {
    return CellSurfaceFace::ZNegative;
  }
  throw std::invalid_argument(
      "surfaceForce/face must be one of xPositive, xNegative, yPositive, "
      "yNegative, zPositive, or zNegative; got '" +
      value + "'");
}

SurfacePatchMode parseSurfacePatchMode(const std::string &value) {
  if (value == "geodesicRadius") {
    return SurfacePatchMode::GeodesicRadius;
  }
  if (value == "oneRing") {
    return SurfacePatchMode::OneRing;
  }
  throw std::invalid_argument(
      "surfaceForce/patchMode must be geodesicRadius or oneRing; got '" +
      value + "'");
}

SurfaceForceWeighting parseSurfaceForceWeighting(const std::string &value) {
  if (value == "nodalArea") {
    return SurfaceForceWeighting::NodalArea;
  }
  if (value == "uniform") {
    return SurfaceForceWeighting::Uniform;
  }
  throw std::invalid_argument(
      "surfaceForce/weighting must be nodalArea or uniform; got '" + value +
      "'");
}

const char *cellSurfaceFaceName(CellSurfaceFace value) {
  switch (value) {
  case CellSurfaceFace::XPositive:
    return "xPositive";
  case CellSurfaceFace::XNegative:
    return "xNegative";
  case CellSurfaceFace::YPositive:
    return "yPositive";
  case CellSurfaceFace::YNegative:
    return "yNegative";
  case CellSurfaceFace::ZPositive:
    return "zPositive";
  case CellSurfaceFace::ZNegative:
    return "zNegative";
  }
  return "unknown";
}

const char *surfacePatchModeName(SurfacePatchMode value) {
  switch (value) {
  case SurfacePatchMode::GeodesicRadius:
    return "geodesicRadius";
  case SurfacePatchMode::OneRing:
    return "oneRing";
  }
  return "unknown";
}

const char *surfaceForceWeightingName(SurfaceForceWeighting value) {
  switch (value) {
  case SurfaceForceWeighting::NodalArea:
    return "nodalArea";
  case SurfaceForceWeighting::Uniform:
    return "uniform";
  }
  return "unknown";
}

class HemoCellSurfaceForce::GatherTargetGeometry
    : public HemoCellFunctional {
public:
  GatherTargetGeometry(std::map<int, VertexGeometry> &vertices,
                       plint targetBaseCellId, unsigned char cellType)
      : vertices_(vertices), targetBaseCellId_(targetBaseCellId),
        cellType_(cellType) {}

  void processGenericBlocks(plb::Box3D,
                            std::vector<plb::AtomicBlock3D *> blocks) {
    HemoCellParticleField *particleField =
        dynamic_cast<HemoCellParticleField *>(blocks[0]);
    if (!particleField || !particleField->cellFields ||
        particleField->cellFields->number_of_cells <= 0) {
      return;
    }

    for (const HemoCellParticle &particle : particleField->particles) {
      if (!particleField->isContainedABS(particle.sv.position,
                                         particleField->localDomain)) {
        continue;
      }
      if (particle.sv.celltype != cellType_ ||
          particleField->cellFields->base_cell_id(particle.sv.cellId) !=
              targetBaseCellId_) {
        continue;
      }
      VertexGeometry geometry;
      geometry.position = particle.sv.position;
      vertices_[particle.sv.vertexId] = geometry;
    }
  }

  GatherTargetGeometry *clone() const {
    return new GatherTargetGeometry(*this);
  }

private:
  std::map<int, VertexGeometry> &vertices_;
  plint targetBaseCellId_;
  unsigned char cellType_;
};

class HemoCellSurfaceForce::AddSurfaceForceFunctional
    : public HemoCellFunctional {
public:
  AddSurfaceForceFunctional(plint targetBaseCellId, unsigned char cellType,
                            const std::vector<T> &weights,
                            const hemo::Array<T, 3> &totalForce)
      : targetBaseCellId_(targetBaseCellId), cellType_(cellType),
        weights_(weights), totalForce_(totalForce) {}

  void processGenericBlocks(plb::Box3D,
                            std::vector<plb::AtomicBlock3D *> blocks) {
    HemoCellParticleField *particleField =
        dynamic_cast<HemoCellParticleField *>(blocks[0]);
    if (!particleField || !particleField->cellFields ||
        particleField->cellFields->number_of_cells <= 0) {
      return;
    }

    for (HemoCellParticle &particle : particleField->particles) {
      if (!particleField->isContainedABS(particle.sv.position,
                                         particleField->localDomain)) {
        continue;
      }
      if (particle.sv.celltype != cellType_ ||
          particleField->cellFields->base_cell_id(particle.sv.cellId) !=
              targetBaseCellId_ ||
          particle.sv.vertexId >= weights_.size()) {
        continue;
      }

      const T weight = weights_[particle.sv.vertexId];
      if (weight <= 0.0) {
        continue;
      }
      for (int d = 0; d < 3; ++d) {
        particle.sv.force[d] += weight * totalForce_[d];
      }
    }
  }

  AddSurfaceForceFunctional *clone() const {
    return new AddSurfaceForceFunctional(*this);
  }

private:
  plint targetBaseCellId_;
  unsigned char cellType_;
  std::vector<T> weights_;
  hemo::Array<T, 3> totalForce_;
};

HemoCellSurfaceForce::HemoCellSurfaceForce(
    HemoCell &hemocell, const std::string &cellType, plint targetBaseCellId,
    CellSurfaceFace face, const hemo::Array<T, 3> &forcePicoNewton,
    SurfacePatchMode patchMode, T patchGeodesicRadiusMicrometer,
    SurfaceForceWeighting weighting)
    : hemocell_(hemocell), cellField_(*hemocell.cellfields->operator[](cellType)),
      cellType_(cellType), targetBaseCellId_(targetBaseCellId), face_(face),
      forcePicoNewton_(forcePicoNewton), forceLbm_({0.0, 0.0, 0.0}),
      patchMode_(patchMode),
      patchGeodesicRadiusMicrometer_(patchGeodesicRadiusMicrometer),
      weighting_(weighting) {
  if (targetBaseCellId_ < 0) {
    throw std::invalid_argument("surfaceForce/targetCellId must be non-negative");
  }
  if (!finiteArray(forcePicoNewton_)) {
    throw std::invalid_argument("surfaceForce force components must be finite");
  }
  if (!std::isfinite(patchGeodesicRadiusMicrometer_) ||
      (patchMode_ == SurfacePatchMode::GeodesicRadius &&
       patchGeodesicRadiusMicrometer_ <= 0.0)) {
    throw std::invalid_argument(
        "surfaceForce/patchGeodesicRadius must be finite and positive for "
        "geodesicRadius mode");
  }
  if (!(param::df > 0.0) || !std::isfinite(param::df)) {
    throw std::runtime_error(
        "Cannot convert surface force because param::df is not positive");
  }

  for (int d = 0; d < 3; ++d) {
    forceLbm_[d] = forcePicoNewton_[d] * 1.0e-12 / param::df;
  }
  if (!finiteArray(forceLbm_)) {
    throw std::invalid_argument(
        "surfaceForce force is too large to represent in lattice units");
  }
  validateTimescales();

  const hemo::Array<T, 3> direction = faceDirection(face_);
  const T forceMagnitude = norm(forcePicoNewton_);
  pcout << "(HemoCellSurfaceForce) target base cell " << targetBaseCellId_
        << ", face " << cellSurfaceFaceName(face_) << ", total force [pN] = ("
        << forcePicoNewton_[0] << ", " << forcePicoNewton_[1] << ", "
        << forcePicoNewton_[2] << ")" << std::endl;
  if (forceMagnitude > 0.0) {
    T cosine = dot(forcePicoNewton_, direction) / forceMagnitude;
    cosine = std::max<T>(-1.0, std::min<T>(1.0, cosine));
    const T angleDegrees = std::acos(cosine) * 180.0 / PI;
    pcout << "(HemoCellSurfaceForce) angle between selected face direction "
             "and force vector: "
          << angleDegrees << " degrees" << std::endl;
  } else {
    pcout << "(HemoCellSurfaceForce) zero-force regression mode enabled"
          << std::endl;
  }
}

void HemoCellSurfaceForce::validateTimescales() const {
  if (cellField_.timescale != 1) {
    throw std::runtime_error(
        "HemoCellSurfaceForce requires target material timescale equal to 1");
  }
  if (hemocell_.cellfields->particleVelocityUpdateTimescale != 1) {
    throw std::runtime_error(
        "HemoCellSurfaceForce requires particle velocity timescale equal to 1");
  }
  if (hemocell_.cellfields->adhesionTimescale != 1) {
    throw std::runtime_error(
        "HemoCellSurfaceForce requires adhesion timescale equal to 1");
  }
}

std::map<int, HemoCellSurfaceForce::VertexGeometry>
HemoCellSurfaceForce::gatherTargetGeometry() const {
  std::map<int, VertexGeometry> vertices;
  std::vector<plb::MultiBlock3D *> blocks;
  blocks.push_back(cellField_.getParticleField3D());
  applyProcessingFunctional(
      new GatherTargetGeometry(vertices, targetBaseCellId_, cellField_.ctype),
      cellField_.getParticleField3D()->getBoundingBox(), blocks);
  HemoCellGatheringFunctional<VertexGeometry>::gather(vertices);
  return vertices;
}

void HemoCellSurfaceForce::initializeFromCurrentGeometry() {
  if (initialized_) {
    throw std::runtime_error(
        "HemoCellSurfaceForce has already been initialized");
  }
  validateTimescales();
  hemocell_.cellfields->syncEnvelopes();

  const std::map<int, VertexGeometry> gathered = gatherTargetGeometry();
  if (gathered.size() != static_cast<std::size_t>(cellField_.numVertex)) {
    std::ostringstream message;
    message << "Target base cell " << targetBaseCellId_ << " has "
            << gathered.size() << " owner vertices after MPI gather, expected "
            << cellField_.numVertex;
    throw std::runtime_error(message.str());
  }

  std::vector<hemo::Array<T, 3>> positions(cellField_.numVertex);
  initialCellCenter_ = {0.0, 0.0, 0.0};
  for (int vertexId = 0; vertexId < cellField_.numVertex; ++vertexId) {
    const std::map<int, VertexGeometry>::const_iterator vertex =
        gathered.find(vertexId);
    if (vertex == gathered.end() || !finiteArray(vertex->second.position)) {
      throw std::runtime_error(
          "Target cell does not contain a finite, contiguous vertexId set");
    }
    positions[vertexId] = vertex->second.position;
    initialCellCenter_ += positions[vertexId];
  }
  initialCellCenter_ /= static_cast<T>(cellField_.numVertex);

  const auto findSeed = [&](const hemo::Array<T, 3> &direction) -> plint {
    const T tieTolerance = 1.0e-12;
    T bestAlignment = -std::numeric_limits<T>::infinity();
    plint bestVertex = -1;
    for (int vertexId = 0; vertexId < cellField_.numVertex; ++vertexId) {
      const hemo::Array<T, 3> radial =
          subtract(positions[vertexId], initialCellCenter_);
      const T radialNorm = norm(radial);
      const T projection = dot(radial, direction);
      if (!(radialNorm > 0.0) || !(projection > 0.0)) {
        continue;
      }
      const T alignment = projection / radialNorm;
      if (alignment > bestAlignment + tieTolerance ||
          (std::fabs(alignment - bestAlignment) <= tieTolerance &&
           (bestVertex < 0 || vertexId < bestVertex))) {
        bestAlignment = alignment;
        bestVertex = vertexId;
      }
    }
    if (bestVertex < 0) {
      throw std::runtime_error(
          "Could not find a membrane vertex on the requested cell face");
    }
    return bestVertex;
  };

  const hemo::Array<T, 3> direction = faceDirection(face_);
  seedVertexId_ = findSeed(direction);
  const hemo::Array<T, 3> oppositeDirection =
      {-direction[0], -direction[1], -direction[2]};
  const plint oppositeSeed = findSeed(oppositeDirection);

  std::vector<std::map<int, T>> graph(cellField_.numVertex);
  std::vector<T> nodalAreas(cellField_.numVertex, 0.0);
  const auto addEdge = [&](int a, int b) {
    const T length = norm(subtract(positions[a], positions[b]));
    if (!(length > 0.0) || !std::isfinite(length)) {
      throw std::runtime_error("Membrane triangle graph contains a zero edge");
    }
    const std::map<int, T>::iterator old = graph[a].find(b);
    if (old == graph[a].end() || length < old->second) {
      graph[a][b] = length;
      graph[b][a] = length;
    }
  };

  for (const hemo::Array<plint, 3> &triangle : cellField_.triangle_list) {
    const int a = static_cast<int>(triangle[0]);
    const int b = static_cast<int>(triangle[1]);
    const int c = static_cast<int>(triangle[2]);
    if (a < 0 || b < 0 || c < 0 || a >= cellField_.numVertex ||
        b >= cellField_.numVertex || c >= cellField_.numVertex || a == b ||
        b == c || c == a) {
      throw std::runtime_error("Membrane triangle list contains an invalid ID");
    }
    addEdge(a, b);
    addEdge(b, c);
    addEdge(c, a);

    const T area = 0.5 * norm(cross(subtract(positions[b], positions[a]),
                                    subtract(positions[c], positions[a])));
    if (!(area > 0.0) || !std::isfinite(area)) {
      throw std::runtime_error("Membrane triangle has invalid initial area");
    }
    nodalAreas[a] += area / 3.0;
    nodalAreas[b] += area / 3.0;
    nodalAreas[c] += area / 3.0;
  }

  std::vector<T> distances(cellField_.numVertex,
                           std::numeric_limits<T>::infinity());
  distances[seedVertexId_] = 0.0;
  typedef std::pair<T, int> QueueEntry;
  std::priority_queue<QueueEntry, std::vector<QueueEntry>,
                      std::greater<QueueEntry>>
      queue;
  queue.push(QueueEntry(0.0, static_cast<int>(seedVertexId_)));
  while (!queue.empty()) {
    const T distance = queue.top().first;
    const int vertexId = queue.top().second;
    queue.pop();
    if (distance > distances[vertexId]) {
      continue;
    }
    for (const std::pair<const int, T> &edge : graph[vertexId]) {
      const T candidate = distance + edge.second;
      if (candidate < distances[edge.first]) {
        distances[edge.first] = candidate;
        queue.push(QueueEntry(candidate, edge.first));
      }
    }
  }
  for (int vertexId = 0; vertexId < cellField_.numVertex; ++vertexId) {
    if (graph[vertexId].empty() || !std::isfinite(distances[vertexId])) {
      throw std::runtime_error(
          "Membrane triangle graph is disconnected or incomplete");
    }
  }

  std::vector<bool> selected(cellField_.numVertex, false);
  if (patchMode_ == SurfacePatchMode::OneRing) {
    selected[seedVertexId_] = true;
    for (const std::pair<const int, T> &edge : graph[seedVertexId_]) {
      selected[edge.first] = true;
    }
  } else {
    const T radiusLbm =
        patchGeodesicRadiusMicrometer_ * 1.0e-6 / param::dx;
    for (int vertexId = 0; vertexId < cellField_.numVertex; ++vertexId) {
      selected[vertexId] = distances[vertexId] <= radiusLbm;
    }
  }

  if (selected[oppositeSeed]) {
    throw std::runtime_error(
        "Selected membrane patch reaches the opposite face seed; reduce the "
        "geodesic radius");
  }

  selectedVertices_.clear();
  T normalizer = 0.0;
  for (int vertexId = 0; vertexId < cellField_.numVertex; ++vertexId) {
    if (!selected[vertexId]) {
      continue;
    }
    SurfaceForceVertex vertex;
    vertex.vertexId = vertexId;
    vertex.initialPosition = positions[vertexId];
    const hemo::Array<T, 3> radial =
        subtract(positions[vertexId], initialCellCenter_);
    vertex.alignment = dot(radial, direction) / norm(radial);
    vertex.geodesicDistance = distances[vertexId];
    vertex.nodalArea = nodalAreas[vertexId];
    selectedVertices_.push_back(vertex);
    normalizer += weighting_ == SurfaceForceWeighting::NodalArea
                      ? vertex.nodalArea
                      : 1.0;
  }

  if (selectedVertices_.empty() ||
      selectedVertices_.size() >= static_cast<std::size_t>(cellField_.numVertex) ||
      !(normalizer > 0.0) || !std::isfinite(normalizer)) {
    throw std::runtime_error("Surface-force patch is empty or invalid");
  }

  vertexWeights_.assign(cellField_.numVertex, 0.0);
  for (SurfaceForceVertex &vertex : selectedVertices_) {
    vertex.weight = weighting_ == SurfaceForceWeighting::NodalArea
                        ? vertex.nodalArea / normalizer
                        : 1.0 / normalizer;
    vertexWeights_[vertex.vertexId] = vertex.weight;
  }
  initialized_ = true;
  validateSelectedState(false);

  pcout << "(HemoCellSurfaceForce) selected seed vertex " << seedVertexId_
        << " and " << selectedVertices_.size() << " membrane vertices using "
        << surfacePatchModeName(patchMode_) << " / "
        << surfaceForceWeightingName(weighting_) << "; sum(weights) = "
        << sumWeights() << std::endl;
}

void HemoCellSurfaceForce::validateSelectedState(
    bool validateCurrentTarget) const {
  if (!initialized_ || seedVertexId_ < 0 || selectedVertices_.empty() ||
      seedVertexId_ >= cellField_.numVertex ||
      vertexWeights_.size() != static_cast<std::size_t>(cellField_.numVertex) ||
      !finiteArray(initialCellCenter_)) {
    throw std::runtime_error("Surface-force grip state is not initialized");
  }

  std::set<plint> ids;
  T weightSum = 0.0;
  int zeroDistanceCount = 0;
  bool containsSeed = false;
  for (const SurfaceForceVertex &vertex : selectedVertices_) {
    if (vertex.vertexId < 0 || vertex.vertexId >= cellField_.numVertex ||
        !ids.insert(vertex.vertexId).second || !std::isfinite(vertex.weight) ||
        vertex.weight < 0.0 || !finiteArray(vertex.initialPosition) ||
        !std::isfinite(vertex.alignment) ||
        !std::isfinite(vertex.geodesicDistance) ||
        vertex.geodesicDistance < 0.0 || !std::isfinite(vertex.nodalArea) ||
        vertex.nodalArea <= 0.0) {
      throw std::runtime_error(
          "Surface-force grip state contains an invalid vertex record");
    }
    if (!nearlyEqual(vertexWeights_[vertex.vertexId], vertex.weight)) {
      throw std::runtime_error(
          "Surface-force lookup weights do not match selected vertices");
    }
    if (vertex.vertexId == seedVertexId_) {
      containsSeed = true;
    }
    if (nearlyEqual(vertex.geodesicDistance, 0.0)) {
      ++zeroDistanceCount;
    }
    weightSum += vertex.weight;
  }
  if (!containsSeed || zeroDistanceCount != 1) {
    throw std::runtime_error(
        "Surface-force grip must contain exactly one zero-distance seed");
  }
  if (!nearlyEqual(weightSum, 1.0, 1.0e-10)) {
    throw std::runtime_error("Surface-force vertex weights do not sum to one");
  }

  if (validateCurrentTarget) {
    const std::map<int, VertexGeometry> current = gatherTargetGeometry();
    if (current.size() != static_cast<std::size_t>(cellField_.numVertex)) {
      throw std::runtime_error(
          "Checkpoint target cell does not contain the complete vertex set");
    }
    for (const SurfaceForceVertex &vertex : selectedVertices_) {
      if (current.find(vertex.vertexId) == current.end()) {
        throw std::runtime_error(
            "Checkpoint target cell is missing a selected vertexId");
      }
    }
  }
}

void HemoCellSurfaceForce::applyForce(T loadScale) {
  if (!std::isfinite(loadScale) || loadScale < 0.0 || loadScale > 1.0) {
    throw std::invalid_argument(
        "HemoCellSurfaceForce loadScale must be finite and within [0,1]");
  }
  validateTimescales();
  validateSelectedState(false);

  hemo::Array<T, 3> appliedForce = {forceLbm_[0] * loadScale,
                                    forceLbm_[1] * loadScale,
                                    forceLbm_[2] * loadScale};
  std::vector<plb::MultiBlock3D *> blocks;
  blocks.push_back(cellField_.getParticleField3D());
  applyProcessingFunctional(
      new AddSurfaceForceFunctional(targetBaseCellId_, cellField_.ctype,
                                    vertexWeights_, appliedForce),
      cellField_.getParticleField3D()->getBoundingBox(), blocks);
  hemocell_.cellfields->syncEnvelopes();
}

void HemoCellSurfaceForce::saveGripState(const std::string &filename) const {
  validateSelectedState(false);
  if (plb::global::mpi().isMainProcessor()) {
    std::ofstream output(filename.c_str());
    if (!output) {
      throw std::runtime_error("Could not write surface-force state to " +
                               filename);
    }
    output << std::setprecision(17);
    output << "HEMOCELL_SURFACE_FORCE 1\n";
    output << "cellType " << cellType_ << '\n';
    output << "targetBaseCellId " << targetBaseCellId_ << '\n';
    output << "face " << cellSurfaceFaceName(face_) << '\n';
    output << "forcePicoNewton " << forcePicoNewton_[0] << ' '
           << forcePicoNewton_[1] << ' ' << forcePicoNewton_[2] << '\n';
    output << "patchMode " << surfacePatchModeName(patchMode_) << '\n';
    output << "patchGeodesicRadiusMicrometer "
           << patchGeodesicRadiusMicrometer_ << '\n';
    output << "weighting " << surfaceForceWeightingName(weighting_) << '\n';
    output << "seedVertexId " << seedVertexId_ << '\n';
    output << "initialCellCenter " << initialCellCenter_[0] << ' '
           << initialCellCenter_[1] << ' ' << initialCellCenter_[2] << '\n';
    output << "selectedCount " << selectedVertices_.size() << '\n';
    for (const SurfaceForceVertex &vertex : selectedVertices_) {
      output << "vertex " << vertex.vertexId << ' ' << vertex.weight << ' '
             << vertex.initialPosition[0] << ' ' << vertex.initialPosition[1]
             << ' ' << vertex.initialPosition[2] << ' ' << vertex.alignment
             << ' ' << vertex.geodesicDistance << ' ' << vertex.nodalArea
             << '\n';
    }
    output << "end\n";
  }
  plb::global::mpi().barrier();
}

void HemoCellSurfaceForce::restoreGripState(const std::string &filename) {
  if (initialized_) {
    throw std::runtime_error(
        "Cannot restore over an initialized HemoCellSurfaceForce state");
  }
  std::ifstream input(filename.c_str());
  if (!input) {
    throw std::runtime_error(
        "Surface-force checkpoint sidecar is required but missing: " +
        filename);
  }

  requireKey(input, "HEMOCELL_SURFACE_FORCE");
  int version = 0;
  if (!(input >> version) || version != 1) {
    throw std::runtime_error("Unsupported surface-force state version");
  }

  std::string savedCellType;
  plint savedTarget = -1;
  std::string savedFace;
  hemo::Array<T, 3> savedForce = {0.0, 0.0, 0.0};
  std::string savedPatchMode;
  T savedRadius = 0.0;
  std::string savedWeighting;
  std::size_t selectedCount = 0;

  requireKey(input, "cellType");
  input >> savedCellType;
  requireKey(input, "targetBaseCellId");
  input >> savedTarget;
  requireKey(input, "face");
  input >> savedFace;
  requireKey(input, "forcePicoNewton");
  input >> savedForce[0] >> savedForce[1] >> savedForce[2];
  requireKey(input, "patchMode");
  input >> savedPatchMode;
  requireKey(input, "patchGeodesicRadiusMicrometer");
  input >> savedRadius;
  requireKey(input, "weighting");
  input >> savedWeighting;
  requireKey(input, "seedVertexId");
  input >> seedVertexId_;
  requireKey(input, "initialCellCenter");
  input >> initialCellCenter_[0] >> initialCellCenter_[1] >>
      initialCellCenter_[2];
  requireKey(input, "selectedCount");
  input >> selectedCount;
  if (!input) {
    throw std::runtime_error("Surface-force state header is truncated");
  }

  if (savedCellType != cellType_ || savedTarget != targetBaseCellId_ ||
      savedFace != cellSurfaceFaceName(face_) ||
      savedPatchMode != surfacePatchModeName(patchMode_) ||
      savedWeighting != surfaceForceWeightingName(weighting_) ||
      !nearlyEqual(savedRadius, patchGeodesicRadiusMicrometer_)) {
    throw std::runtime_error(
        "Surface-force checkpoint state does not match current XML settings");
  }
  for (int d = 0; d < 3; ++d) {
    if (!nearlyEqual(savedForce[d], forcePicoNewton_[d])) {
      throw std::runtime_error(
          "Checkpoint force vector does not match current XML settings");
    }
  }
  if (selectedCount == 0 ||
      selectedCount >= static_cast<std::size_t>(cellField_.numVertex)) {
    throw std::runtime_error("Checkpoint selected vertex count is invalid");
  }

  selectedVertices_.clear();
  vertexWeights_.assign(cellField_.numVertex, 0.0);
  for (std::size_t i = 0; i < selectedCount; ++i) {
    requireKey(input, "vertex");
    SurfaceForceVertex vertex;
    input >> vertex.vertexId >> vertex.weight >> vertex.initialPosition[0] >>
        vertex.initialPosition[1] >> vertex.initialPosition[2] >>
        vertex.alignment >> vertex.geodesicDistance >> vertex.nodalArea;
    if (!input || vertex.vertexId < 0 ||
        vertex.vertexId >= cellField_.numVertex) {
      throw std::runtime_error(
          "Surface-force state contains a malformed vertex record");
    }
    selectedVertices_.push_back(vertex);
    vertexWeights_[vertex.vertexId] = vertex.weight;
  }
  requireKey(input, "end");
  initialized_ = true;
  validateTimescales();
  validateSelectedState(true);

  pcout << "(HemoCellSurfaceForce) restored seed vertex " << seedVertexId_
        << " and " << selectedVertices_.size() << " selected vertices from "
        << filename << std::endl;
}

void HemoCellSurfaceForce::saveSelectedVerticesCsv(
    const std::string &filename) const {
  validateSelectedState(false);
  if (plb::global::mpi().isMainProcessor()) {
    std::ofstream output(filename.c_str());
    if (!output) {
      throw std::runtime_error("Could not write selected vertices to " +
                               filename);
    }
    output << std::setprecision(17);
    output << "target_cell_id,face,vertex_id,x_um,y_um,z_um,alignment,"
              "geodesic_distance_um,nodal_area_um2,weight\n";
    const T lengthToMicrometer = param::dx * 1.0e6;
    const T areaToSquareMicrometer = lengthToMicrometer * lengthToMicrometer;
    for (const SurfaceForceVertex &vertex : selectedVertices_) {
      output << targetBaseCellId_ << ',' << cellSurfaceFaceName(face_) << ','
             << vertex.vertexId << ','
             << vertex.initialPosition[0] * lengthToMicrometer << ','
             << vertex.initialPosition[1] * lengthToMicrometer << ','
             << vertex.initialPosition[2] * lengthToMicrometer << ','
             << vertex.alignment << ','
             << vertex.geodesicDistance * lengthToMicrometer << ','
             << vertex.nodalArea * areaToSquareMicrometer << ','
             << vertex.weight << '\n';
    }
  }
  plb::global::mpi().barrier();
}

const std::vector<SurfaceForceVertex> &
HemoCellSurfaceForce::selectedVertices() const {
  return selectedVertices_;
}

hemo::Array<T, 3> HemoCellSurfaceForce::forcePicoNewton() const {
  return forcePicoNewton_;
}

hemo::Array<T, 3> HemoCellSurfaceForce::forceLbm() const {
  return forceLbm_;
}

hemo::Array<T, 3> HemoCellSurfaceForce::initialCellCenter() const {
  return initialCellCenter_;
}

hemo::Array<T, 3> HemoCellSurfaceForce::currentWeightedPosition() const {
  validateSelectedState(false);
  const std::map<int, VertexGeometry> current = gatherTargetGeometry();
  hemo::Array<T, 3> weightedPosition = {0.0, 0.0, 0.0};
  for (const SurfaceForceVertex &vertex : selectedVertices_) {
    const std::map<int, VertexGeometry>::const_iterator position =
        current.find(vertex.vertexId);
    if (position == current.end()) {
      throw std::runtime_error(
          "Target cell is missing a selected surface-force vertex");
    }
    for (int d = 0; d < 3; ++d) {
      weightedPosition[d] += vertex.weight * position->second.position[d];
    }
  }
  return weightedPosition;
}

plint HemoCellSurfaceForce::seedVertexId() const { return seedVertexId_; }

plint HemoCellSurfaceForce::targetBaseCellId() const {
  return targetBaseCellId_;
}

CellSurfaceFace HemoCellSurfaceForce::face() const { return face_; }

T HemoCellSurfaceForce::sumWeights() const {
  T total = 0.0;
  for (const SurfaceForceVertex &vertex : selectedVertices_) {
    total += vertex.weight;
  }
  return total;
}

} // namespace hemo
