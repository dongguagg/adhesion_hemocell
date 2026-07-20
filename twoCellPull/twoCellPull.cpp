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

#include "hemocell.h"
#include "rbcHighOrderModel.h"
#include "helper/cellInfo.h"
#include "helper/genericFunctions.h"
#include "helper/hemoCellSurfaceForce.h"

#include "palabos3D.h"
#include "palabos3D.hh"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

using namespace hemo;

namespace {

struct DiagnosticVertex {
  hemo::Array<T, 3> position = {0.0, 0.0, 0.0};
  int baseCellId = -1;
  int vertexId = -1;
};

class GatherTwoCellVertices : public HemoCellFunctional {
public:
  GatherTwoCellVertices(std::map<int, DiagnosticVertex> &vertices,
                        unsigned char cellType, int numVertices)
      : vertices_(vertices), cellType_(cellType), numVertices_(numVertices) {}

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
                                         particleField->localDomain) ||
          particle.sv.celltype != cellType_) {
        continue;
      }
      const int baseCellId =
          particleField->cellFields->base_cell_id(particle.sv.cellId);
      if (baseCellId < 0 || baseCellId > 1 ||
          particle.sv.vertexId >= numVertices_) {
        continue;
      }
      DiagnosticVertex value;
      value.position = particle.sv.position;
      value.baseCellId = baseCellId;
      value.vertexId = particle.sv.vertexId;
      vertices_[baseCellId * numVertices_ + particle.sv.vertexId] = value;
    }
  }

  GatherTwoCellVertices *clone() const {
    return new GatherTwoCellVertices(*this);
  }

private:
  std::map<int, DiagnosticVertex> &vertices_;
  unsigned char cellType_;
  int numVertices_;
};

T forceScale(unsigned int iteration, unsigned int tRelax,
             unsigned int tForceRamp) {
  if (iteration < tRelax) {
    return 0.0;
  }
  if (tForceRamp == 0 || iteration >= tRelax + tForceRamp) {
    return 1.0;
  }
  return static_cast<T>(iteration - tRelax) /
         static_cast<T>(tForceRamp);
}

const char *simulationStage(unsigned int iteration, unsigned int tRelax,
                            unsigned int tForceRamp) {
  if (iteration < tRelax) {
    return "relaxation";
  }
  if (tForceRamp > 0 && iteration < tRelax + tForceRamp) {
    return "ramp";
  }
  return "constant_force";
}

std::map<int, DiagnosticVertex>
gatherDiagnosticVertices(HemoCell &hemocell, HemoCellField &cellField) {
  std::map<int, DiagnosticVertex> vertices;
  std::vector<plb::MultiBlock3D *> blocks;
  blocks.push_back(cellField.getParticleField3D());
  applyProcessingFunctional(
      new GatherTwoCellVertices(vertices, cellField.ctype,
                                cellField.numVertex),
      cellField.getParticleField3D()->getBoundingBox(), blocks);
  HemoCellGatheringFunctional<DiagnosticVertex>::gather(vertices);
  return vertices;
}

std::map<int, CellInformation> gatherCellInformation(HemoCell &hemocell) {
  std::map<int, CellInformation> local;
  CellInformationFunctionals::calculateCellInformation(&hemocell, local);

  std::map<int, CellInformation> owners;
  for (const std::pair<const int, CellInformation> &entry : local) {
    if (entry.second.centerLocal && entry.second.base_cell_id >= 0 &&
        entry.second.base_cell_id <= 1) {
      owners[entry.second.base_cell_id] = entry.second;
    }
  }
  HemoCellGatheringFunctional<CellInformation>::gather(owners);
  return owners;
}

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

T distance(const hemo::Array<T, 3> &a, const hemo::Array<T, 3> &b) {
  const hemo::Array<T, 3> delta = subtract(a, b);
  return std::sqrt(delta[0] * delta[0] + delta[1] * delta[1] +
                   delta[2] * delta[2]);
}

void validateSimulationConfiguration(unsigned int tRelax,
                                     unsigned int tForceRamp,
                                     unsigned int tmax,
                                     unsigned int tmeas,
                                     unsigned int tcheckpoint) {
  if (tRelax > tmax ||
      static_cast<unsigned long long>(tRelax) + tForceRamp > tmax) {
    throw std::invalid_argument(
        "sim requires tRelax <= tmax and tRelax+tForceRamp <= tmax");
  }
  if (tmeas == 0 || tcheckpoint == 0) {
    throw std::invalid_argument(
        "sim/tmeas and sim/tcheckpoint must both be positive");
  }
}

void writeCsvHeader(plb_ofstream &stream) {
  stream
      << "iteration,load_iteration,stage,target_cell_id,face,seed_vertex_id,"
         "load_scale,configured_force_x_pN,configured_force_y_pN,"
         "configured_force_z_pN,applied_force_x_pN,applied_force_y_pN,"
         "applied_force_z_pN,selected_vertex_count,sum_weights,"
         "grip_center_x_um,grip_center_y_um,grip_center_z_um,"
         "target_center_x_um,target_center_y_um,target_center_z_um,"
         "other_center_x_um,other_center_y_um,other_center_z_um,"
         "relative_center_x_um,relative_center_y_um,relative_center_z_um,"
         "minimum_cell_cell_distance_um,minimum_lower_cell_wall_distance_um,"
         "torque_x_pN_um,torque_y_pN_um,torque_z_pN_um,"
         "lower_area_um2,lower_volume_um3,lower_bbox_x_min_um,"
         "lower_bbox_x_max_um,lower_bbox_y_min_um,lower_bbox_y_max_um,"
         "lower_bbox_z_min_um,lower_bbox_z_max_um,upper_area_um2,"
         "upper_volume_um3,upper_bbox_x_min_um,upper_bbox_x_max_um,"
         "upper_bbox_y_min_um,upper_bbox_y_max_um,upper_bbox_z_min_um,"
         "upper_bbox_z_max_um"
      << std::endl;
}

void writePullState(HemoCell &hemocell, HemoCellField &cellField,
                    HemoCellSurfaceForce &surfaceForce, plb_ofstream &stream,
                    unsigned int loadIteration, unsigned int tRelax,
                    unsigned int tForceRamp, T loadScale) {
  const std::map<int, DiagnosticVertex> vertices =
      gatherDiagnosticVertices(hemocell, cellField);
  const std::map<int, CellInformation> cells =
      gatherCellInformation(hemocell);
  const int numVertices = cellField.numVertex;
  if (vertices.size() != static_cast<std::size_t>(2 * numVertices) ||
      cells.find(0) == cells.end() || cells.find(1) == cells.end()) {
    throw std::runtime_error(
        "TwoCellPull diagnostics require two complete RBCs with base IDs 0 and 1");
  }

  T minimumCellCellDistance = std::numeric_limits<T>::infinity();
  for (int lowerVertex = 0; lowerVertex < numVertices; ++lowerVertex) {
    const hemo::Array<T, 3> &lower =
        vertices.at(lowerVertex).position;
    for (int upperVertex = 0; upperVertex < numVertices; ++upperVertex) {
      const hemo::Array<T, 3> &upper =
          vertices.at(numVertices + upperVertex).position;
      minimumCellCellDistance =
          std::min(minimumCellCellDistance, distance(lower, upper));
    }
  }

  T minimumLowerWallDistance = std::numeric_limits<T>::infinity();
  for (int vertexId = 0; vertexId < numVertices; ++vertexId) {
    const hemo::Array<T, 3> &position = vertices.at(vertexId).position;
    const T dxToGrid = position[0] - std::floor(position[0] + 0.5);
    const T dyToGrid = position[1] - std::floor(position[1] + 0.5);
    const T wallDistance =
        std::sqrt(dxToGrid * dxToGrid + dyToGrid * dyToGrid +
                  position[2] * position[2]);
    minimumLowerWallDistance =
        std::min(minimumLowerWallDistance, wallDistance);
  }

  hemo::Array<T, 3> gripCenter = {0.0, 0.0, 0.0};
  const int targetCell = static_cast<int>(surfaceForce.targetBaseCellId());
  const int otherCell = targetCell == 0 ? 1 : 0;
  for (const SurfaceForceVertex &selected :
       surfaceForce.selectedVertices()) {
    const DiagnosticVertex &current =
        vertices.at(targetCell * numVertices + selected.vertexId);
    for (int d = 0; d < 3; ++d) {
      gripCenter[d] += selected.weight * current.position[d];
    }
  }

  const hemo::Array<T, 3> targetCenter = cells.at(targetCell).position;
  const hemo::Array<T, 3> otherCenter = cells.at(otherCell).position;
  const hemo::Array<T, 3> relativeCenter =
      subtract(targetCenter, otherCenter);
  const T lengthToMicrometer = param::dx * 1.0e6;
  const T areaToSquareMicrometer =
      lengthToMicrometer * lengthToMicrometer;
  const T volumeToCubicMicrometer =
      areaToSquareMicrometer * lengthToMicrometer;
  const hemo::Array<T, 3> configuredForce =
      surfaceForce.forcePicoNewton();
  const hemo::Array<T, 3> appliedForce =
      {configuredForce[0] * loadScale, configuredForce[1] * loadScale,
       configuredForce[2] * loadScale};
  hemo::Array<T, 3> leverArm = subtract(gripCenter, targetCenter);
  leverArm *= lengthToMicrometer;
  const hemo::Array<T, 3> torque = cross(leverArm, appliedForce);

  const CellInformation &lower = cells.at(0);
  const CellInformation &upper = cells.at(1);
  stream << std::setprecision(17) << hemocell.iter << ',' << loadIteration
         << ',' << simulationStage(loadIteration, tRelax, tForceRamp) << ','
         << targetCell << ',' << cellSurfaceFaceName(surfaceForce.face()) << ','
         << surfaceForce.seedVertexId() << ',' << loadScale << ','
         << configuredForce[0] << ',' << configuredForce[1] << ','
         << configuredForce[2] << ',' << appliedForce[0] << ','
         << appliedForce[1] << ',' << appliedForce[2] << ','
         << surfaceForce.selectedVertices().size() << ','
         << surfaceForce.sumWeights() << ','
         << gripCenter[0] * lengthToMicrometer << ','
         << gripCenter[1] * lengthToMicrometer << ','
         << gripCenter[2] * lengthToMicrometer << ','
         << targetCenter[0] * lengthToMicrometer << ','
         << targetCenter[1] * lengthToMicrometer << ','
         << targetCenter[2] * lengthToMicrometer << ','
         << otherCenter[0] * lengthToMicrometer << ','
         << otherCenter[1] * lengthToMicrometer << ','
         << otherCenter[2] * lengthToMicrometer << ','
         << relativeCenter[0] * lengthToMicrometer << ','
         << relativeCenter[1] * lengthToMicrometer << ','
         << relativeCenter[2] * lengthToMicrometer << ','
         << minimumCellCellDistance * lengthToMicrometer << ','
         << minimumLowerWallDistance * lengthToMicrometer << ',' << torque[0]
         << ',' << torque[1] << ',' << torque[2] << ','
         << lower.area * areaToSquareMicrometer << ','
         << lower.volume * volumeToCubicMicrometer << ','
         << lower.bbox[0] * lengthToMicrometer << ','
         << lower.bbox[1] * lengthToMicrometer << ','
         << lower.bbox[2] * lengthToMicrometer << ','
         << lower.bbox[3] * lengthToMicrometer << ','
         << lower.bbox[4] * lengthToMicrometer << ','
         << lower.bbox[5] * lengthToMicrometer << ','
         << upper.area * areaToSquareMicrometer << ','
         << upper.volume * volumeToCubicMicrometer << ','
         << upper.bbox[0] * lengthToMicrometer << ','
         << upper.bbox[1] * lengthToMicrometer << ','
         << upper.bbox[2] * lengthToMicrometer << ','
         << upper.bbox[3] * lengthToMicrometer << ','
         << upper.bbox[4] * lengthToMicrometer << ','
         << upper.bbox[5] * lengthToMicrometer << std::endl;
}

} // namespace

int main(int argc, char *argv[]) {
  if (argc < 2) {
    std::cout << "Usage: " << argv[0] << " <configuration.xml>" << std::endl;
    return 1;
  }

  try {
    HemoCell hemocell(argv[1], argc, argv);
    Config *cfg = hemocell.cfg;

    const T dx = (*cfg)["domain"]["dx"].read<T>();
    const T lx = (*cfg)["domain"]["Lx"].read<T>();
    const T ly = (*cfg)["domain"]["Ly"].read<T>();
    const T lz = (*cfg)["domain"]["Lz"].read<T>();
    const plint nx = static_cast<plint>(lx * 1.0e-6 / dx + 0.5);
    const plint ny = static_cast<plint>(ly * 1.0e-6 / dx + 0.5);
    const plint nz = static_cast<plint>(lz * 1.0e-6 / dx + 0.5);
    const plint fluidEnvelope =
        (*cfg)["domain"]["fluidEnvelope"].read<plint>();
    const plint nprocX =
        (*cfg)["decomposition"]["nprocX"].read<plint>();
    const plint nprocY =
        (*cfg)["decomposition"]["nprocY"].read<plint>();
    const plint nprocZ =
        (*cfg)["decomposition"]["nprocZ"].read<plint>();
    if (nx <= 1 || ny <= 1 || nz <= 1 || nprocX <= 0 || nprocY <= 0 ||
        nprocZ <= 0) {
      throw std::invalid_argument("Invalid domain size or MPI decomposition");
    }
    const plint configuredRanks = nprocX * nprocY * nprocZ;
    if (configuredRanks != global::mpi().getSize()) {
      throw std::invalid_argument(
          "MPI rank count must equal decomposition/nprocX*nprocY*nprocZ");
    }

    pcout << "(TwoCellPull) Calculating base LBM parameters" << std::endl;
    param::lbm_base_parameters(*cfg);
    param::printParameters();
    pcout << "(TwoCellPull) Initializing " << nx << 'x' << ny << 'x' << nz
          << " lattice with a " << nprocX << 'x' << nprocY << 'x' << nprocZ
          << " block decomposition" << std::endl;

    const SparseBlockStructure3D blockStructure =
        createRegularDistribution3D(nx, ny, nz, nprocX, nprocY, nprocZ);
    const MultiBlockManagement3D management(
        blockStructure, new OneToOneThreadAttribution, fluidEnvelope);
    hemocell.initializeLattice(management);
    hemocell.lattice->periodicity().toggleAll(false);
    hemocell.lattice->periodicity().toggle(0, true);
    hemocell.lattice->periodicity().toggle(1, true);
    hemocell.lattice->toggleInternalStatistics(false);

    const Box3D bottom(0, nx - 1, 0, ny - 1, 0, 0);
    const Box3D top(0, nx - 1, 0, ny - 1, nz - 1, nz - 1);
    OnLatticeBoundaryCondition3D<T, DESCRIPTOR> *boundaryCondition =
        createLocalBoundaryCondition3D<T, DESCRIPTOR>();
    boundaryCondition->setVelocityConditionOnBlockBoundaries(*hemocell.lattice,
                                                              bottom);
    boundaryCondition->setVelocityConditionOnBlockBoundaries(*hemocell.lattice,
                                                              top);
    setBoundaryVelocity(*hemocell.lattice, bottom,
                        plb::Array<T, 3>(0.0, 0.0, 0.0));
    setBoundaryVelocity(*hemocell.lattice, top,
                        plb::Array<T, 3>(0.0, 0.0, 0.0));

    hemocell.latticeEquilibrium(1.0, plb::Array<T, 3>(0.0, 0.0, 0.0));
    setExternalVector(*hemocell.lattice, hemocell.lattice->getBoundingBox(),
                      DESCRIPTOR<T>::ExternalField::forceBeginsAt,
                      plb::Array<T, DESCRIPTOR<T>::d>(0.0, 0.0, 0.0));
    hemocell.lattice->initialize();
    hemocell.outputInSiUnits = true;

    hemocell.initializeCellfield();
    hemocell.addCellType<RbcHighOrderModel>("RBC", RBC_FROM_SPHERE);
    HemoCellField &rbc = *hemocell.cellfields->operator[]("RBC");
    hemocell.setSystemPeriodicity(0, true);
    hemocell.setSystemPeriodicity(1, true);
    hemocell.setSystemPeriodicity(2, false);

    const unsigned int materialTimescale =
        (*cfg)["ibm"]["stepMaterialEvery"].read<unsigned int>();
    const unsigned int particleTimescale =
        (*cfg)["ibm"]["stepParticleEvery"].read<unsigned int>();
    hemocell.setMaterialTimeScaleSeparation("RBC", materialTimescale);
    hemocell.setParticleVelocityUpdateTimeScaleSeparation(particleTimescale);
    hemocell.setAdhesionTimeScaleSeperation(1);

    hemocell.setAdhesion(
        (*cfg)["cellCellAdhesion"]["r0"].read<T>(),
        (*cfg)["cellCellAdhesion"]["rc"].read<T>(),
        (*cfg)["cellCellAdhesion"]["epsilon"].read<T>(),
        (*cfg)["cellCellAdhesion"]["D0"].read<T>(),
        (*cfg)["cellCellAdhesion"]["alpha"].read<T>());
    hemocell.setBoundaryAdhesion(
        (*cfg)["cellWallAdhesion"]["r0"].read<T>(),
        (*cfg)["cellWallAdhesion"]["rc"].read<T>(),
        (*cfg)["cellWallAdhesion"]["epsilon"].read<T>(),
        (*cfg)["cellWallAdhesion"]["D0"].read<T>(),
        (*cfg)["cellWallAdhesion"]["alpha"].read<T>());

    hemocell.setOutputs(
        "RBC", {OUTPUT_POSITION, OUTPUT_TRIANGLES, OUTPUT_VELOCITY,
                OUTPUT_FORCE, OUTPUT_FORCE_REPULSION, OUTPUT_VERTEX_ID,
                OUTPUT_CELL_ID});
    hemocell.setFluidOutputs({OUTPUT_VELOCITY, OUTPUT_FORCE});

    const unsigned int warmup =
        (*cfg)["parameters"]["warmup"].read<unsigned int>();
    const unsigned int tRelax =
        (*cfg)["sim"]["tRelax"].read<unsigned int>();
    const unsigned int tForceRamp =
        (*cfg)["sim"]["tForceRamp"].read<unsigned int>();
    const unsigned int tmax = (*cfg)["sim"]["tmax"].read<unsigned int>();
    const unsigned int tmeas =
        (*cfg)["sim"]["tmeas"].read<unsigned int>();
    const unsigned int tcheckpoint =
        (*cfg)["sim"]["tcheckpoint"].read<unsigned int>();
    validateSimulationConfiguration(tRelax, tForceRamp, tmax, tmeas,
                                    tcheckpoint);

    const plint targetCellId =
        (*cfg)["surfaceForce"]["targetCellId"].read<plint>();
    const CellSurfaceFace face = parseCellSurfaceFace(
        (*cfg)["surfaceForce"]["face"].read<std::string>());
    const hemo::Array<T, 3> forcePicoNewton = {
        (*cfg)["surfaceForce"]["forceX"].read<T>(),
        (*cfg)["surfaceForce"]["forceY"].read<T>(),
        (*cfg)["surfaceForce"]["forceZ"].read<T>()};
    const SurfacePatchMode patchMode = parseSurfacePatchMode(
        (*cfg)["surfaceForce"]["patchMode"].read<std::string>());
    const T patchRadius =
        (*cfg)["surfaceForce"]["patchGeodesicRadius"].read<T>();
    const SurfaceForceWeighting weighting = parseSurfaceForceWeighting(
        (*cfg)["surfaceForce"]["weighting"].read<std::string>());
    HemoCellSurfaceForce surfaceForce(
        hemocell, "RBC", targetCellId, face, forcePicoNewton, patchMode,
        patchRadius, weighting);

    plb_ofstream pullLog;
    if (cfg->checkpointed) {
      pcout << "(TwoCellPull) Loading checkpoint" << std::endl;
      hemocell.loadCheckPoint();
      if (hemocell.cellfields->number_of_cells != 2) {
        throw std::runtime_error(
            "TwoCellPull checkpoint must contain exactly two cells");
      }
      surfaceForce.restoreGripState(hemo::global.checkpointDirectory +
                                    "grip_state.dat");
      pullLog.open((global::directories().getOutputDir() +
                    "twoCellPull.csv")
                       .c_str(),
                   std::ofstream::out | std::ofstream::app);
    } else {
      pullLog.open((global::directories().getOutputDir() +
                    "twoCellPull.csv")
                       .c_str());
      writeCsvHeader(pullLog);

      pcout << "(TwoCellPull) Warming up the cell-free fluid for " << warmup
            << " iterations" << std::endl;
      for (unsigned int i = 0; i < warmup; ++i) {
        hemocell.lattice->collideAndStream();
      }

      hemocell.loadParticles();
      if (hemocell.cellfields->number_of_cells != 2) {
        throw std::runtime_error(
            "TwoCellPull requires exactly two cells in RBC.pos");
      }
      surfaceForce.initializeFromCurrentGeometry();
      if (global::mpi().isMainProcessor()) {
        mkpath(hemo::global.checkpointDirectory.c_str(), 0777);
      }
      global::mpi().barrier();
      surfaceForce.saveSelectedVerticesCsv(
          global::directories().getOutputDir() + "grip_vertices.csv");
      surfaceForce.saveGripState(hemo::global.checkpointDirectory +
                                 "grip_state.dat");
      hemocell.writeOutput();
      writePullState(hemocell, rbc, surfaceForce, pullLog, hemocell.iter,
                     tRelax, tForceRamp, 0.0);
    }

    pcout << "(TwoCellPull) Static walls; force starts at iteration "
          << tRelax << " with ramp duration " << tForceRamp << std::endl;
    while (hemocell.iter < tmax) {
      const unsigned int loadIteration = hemocell.iter;
      const T scale = forceScale(loadIteration, tRelax, tForceRamp);
      if (scale > 0.0) {
        surfaceForce.applyForce(scale);
      }
      hemocell.iterate();

      if (hemocell.iter % tmeas == 0) {
        hemocell.writeOutput();
        writePullState(hemocell, rbc, surfaceForce, pullLog, loadIteration,
                       tRelax, tForceRamp, scale);
      }
      if (hemocell.iter % tcheckpoint == 0) {
        hemocell.saveCheckPoint();
        surfaceForce.saveGripState(hemo::global.checkpointDirectory +
                                   "grip_state.dat");
      }
    }

    pullLog.close();
    delete boundaryCondition;
    pcout << "(TwoCellPull) Simulation finished" << std::endl;
  } catch (const std::exception &error) {
    plb::pcerr << "(TwoCellPull) Fatal error: " << error.what() << std::endl;
    return 1;
  }
  return 0;
}
