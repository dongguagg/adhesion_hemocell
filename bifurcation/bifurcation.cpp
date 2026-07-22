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

The library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#include "hemocell.h"
#include "rbcHighOrderModel.h"
#include "cellInfo.h"
#include "fluidInfo.h"
#include "preInlet.h"
#include <helper/voxelizeDomain.h>

#include "palabos3D.h"
#include "palabos3D.hh"

#include <algorithm>
#include <array>
#include <cstdlib>
#include <fstream>
#include <queue>
#include <set>
#include <tuple>
#include <utility>
#include <vector>

using namespace hemo;

namespace {

struct FluidMask {
  explicit FluidMask(Box3D const& domain)
      : box(domain), ny(domain.getNy()), nz(domain.getNz()),
        values(static_cast<std::size_t>(domain.nCells()), 0) {}

  std::size_t index(plint x, plint y, plint z) const {
    return static_cast<std::size_t>(
        ((x - box.x0) * ny + (y - box.y0)) * nz + (z - box.z0));
  }

  bool contains(plint x, plint y, plint z) const {
    return x >= box.x0 && x <= box.x1 && y >= box.y0 && y <= box.y1 &&
           z >= box.z0 && z <= box.z1;
  }

  bool isFluid(plint x, plint y, plint z) const {
    return contains(x, y, z) && values[index(x, y, z)] != 0;
  }

  Box3D box;
  plint ny;
  plint nz;
  std::vector<unsigned char> values;
};

struct OutletSection {
  plint x = 0;
  std::vector<Dot3D> points;
  T centerY = 0;
  T centerZ = 0;
};

using GridPointKey = std::tuple<plint, plint, plint>;

struct BoundaryParticleDiagnostics {
  unsigned long long total = 0;
  unsigned long long preInlet = 0;
  unsigned long long mainDomain = 0;
  unsigned long long inletHits = 0;
  unsigned long long outlet1Hits = 0;
  unsigned long long outlet2Hits = 0;
};

GridPointKey gridPointKey(Dot3D const& point) {
  return std::make_tuple(point.x, point.y, point.z);
}

std::set<GridPointKey> gridPointSet(std::vector<Dot3D> const& points) {
  std::set<GridPointKey> result;
  for (Dot3D const& point : points) {
    result.insert(gridPointKey(point));
  }
  return result;
}

BoundaryParticleDiagnostics collectBoundaryParticleDiagnostics(
    HemoCell& hemocell, std::set<GridPointKey> const& inlet,
    std::set<GridPointKey> const& outlet1,
    std::set<GridPointKey> const& outlet2) {
  std::array<unsigned long long, 6> local = {{0, 0, 0, 0, 0, 0}};
  plb::MultiParticleField3D<HemoCellParticleField>& particleField =
      *hemocell.cellfields->immersedParticles;

  for (plint blockId : particleField.getLocalInfo().getBlocks()) {
    HemoCellParticleField& block = particleField.getComponent(blockId);
    Dot3D const location = block.atomicLattice->getLocation();
    for (Dot3D const& localPoint : block.boundaryParticles) {
      Dot3D const point = localPoint + location;
      GridPointKey const key = gridPointKey(point);
      ++local[0];
      ++local[hemocell.partOfpreInlet ? 1 : 2];
      if (inlet.count(key) != 0) {
        ++local[3];
      }
      if (outlet1.count(key) != 0) {
        ++local[4];
      }
      if (outlet2.count(key) != 0) {
        ++local[5];
      }
    }
  }

  std::array<unsigned long long, 6> globalCounts = {{0, 0, 0, 0, 0, 0}};
  MPI_Allreduce(local.data(), globalCounts.data(),
                static_cast<int>(globalCounts.size()),
                MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

  BoundaryParticleDiagnostics result;
  result.total = globalCounts[0];
  result.preInlet = globalCounts[1];
  result.mainDomain = globalCounts[2];
  result.inletHits = globalCounts[3];
  result.outlet1Hits = globalCounts[4];
  result.outlet2Hits = globalCounts[5];
  return result;
}

void printBoundaryParticleDiagnostics(
    char const* selection, BoundaryParticleDiagnostics const& diagnostics) {
  pcout << "(Bifurcation) " << selection
        << " boundary-surface entries: total=" << diagnostics.total
        << ", preInlet sidewall=" << diagnostics.preInlet
        << ", main STL/open boundaries=" << diagnostics.mainDomain
        << ", main-inlet hits=" << diagnostics.inletHits
        << ", outlet-1 hits=" << diagnostics.outlet1Hits
        << ", outlet-2 hits=" << diagnostics.outlet2Hits << endl;
}

bool configureSolidBoundaryRepulsion(
    HemoCell& hemocell, T kRep, T RepCutoff, unsigned int timestep,
    std::vector<Dot3D> const& inletPoints,
    std::vector<OutletSection> const& outlets) {
  std::set<GridPointKey> const inlet = gridPointSet(inletPoints);
  std::set<GridPointKey> const outlet1 = gridPointSet(outlets[0].points);
  std::set<GridPointKey> const outlet2 = gridPointSet(outlets[1].points);

  // Build the legacy list once for an explicit compatibility diagnostic. The
  // following four-argument call replaces it with the list used by the run.
  hemocell.cellfields->populateBoundaryParticles();
  BoundaryParticleDiagnostics const allBoundary =
      collectBoundaryParticleDiagnostics(hemocell, inlet, outlet1, outlet2);

  hemocell.enableBoundaryParticles(
      kRep, RepCutoff, timestep,
      BoundaryParticleSelection::SolidBounceBackOnly);
  BoundaryParticleDiagnostics const solidOnly =
      collectBoundaryParticleDiagnostics(hemocell, inlet, outlet1, outlet2);

  printBoundaryParticleDiagnostics("AllBoundaryDynamics", allBoundary);
  printBoundaryParticleDiagnostics("SolidBounceBackOnly", solidOnly);
  pcout << "(Bifurcation) Solid-only selection excluded "
        << allBoundary.total - solidOnly.total
        << " non-solid/open boundary-surface entries." << endl;

  const bool legacyOpenBoundariesDetected =
      allBoundary.inletHits > 0 && allBoundary.outlet1Hits > 0 &&
      allBoundary.outlet2Hits > 0;
  const bool openBoundariesExcluded =
      solidOnly.inletHits == 0 && solidOnly.outlet1Hits == 0 &&
      solidOnly.outlet2Hits == 0;
  const bool solidWallsRetained =
      solidOnly.preInlet > 0 && solidOnly.mainDomain > 0;

  if (!legacyOpenBoundariesDetected || !openBoundariesExcluded ||
      !solidWallsRetained) {
    pcout << "(Bifurcation) Boundary-particle selection verification failed: "
             "legacy mode must include all three open sections, solid-only "
             "mode must exclude them, and both preInlet/main solid walls "
             "must remain non-empty."
          << endl;
    return false;
  }
  return true;
}

FluidMask collectGlobalFluidMask(MultiScalarField3D<int>& flagMatrix) {
  FluidMask mask(flagMatrix.getBoundingBox());
  MultiBlockManagement3D const& management =
      flagMatrix.getMultiBlockManagement();

  for (plint blockId : flagMatrix.getLocalInfo().getBlocks()) {
    Box3D const bulk = management.getBulk(blockId);
    ScalarField3D<int>& block = flagMatrix.getComponent(blockId);
    Dot3D const location = block.getLocation();

    for (plint x = bulk.x0; x <= bulk.x1; ++x) {
      for (plint y = bulk.y0; y <= bulk.y1; ++y) {
        for (plint z = bulk.z0; z <= bulk.z1; ++z) {
          if (block.get(x - location.x, y - location.y, z - location.z) ==
              1) {
            mask.values[mask.index(x, y, z)] = 1;
          }
        }
      }
    }
  }

  MPI_Allreduce(MPI_IN_PLACE, mask.values.data(),
                static_cast<int>(mask.values.size()), MPI_UNSIGNED_CHAR,
                MPI_MAX, MPI_COMM_WORLD);
  return mask;
}

std::vector<std::vector<Dot3D>> connectedFluidComponents(
    FluidMask const& mask, plint x) {
  const plint ny = mask.box.getNy();
  const plint nz = mask.box.getNz();
  std::vector<unsigned char> visited(
      static_cast<std::size_t>(ny * nz), 0);
  std::vector<std::vector<Dot3D>> components;

  const auto planeIndex = [&](plint y, plint z) {
    return static_cast<std::size_t>((y - mask.box.y0) * nz +
                                    (z - mask.box.z0));
  };

  for (plint y = mask.box.y0; y <= mask.box.y1; ++y) {
    for (plint z = mask.box.z0; z <= mask.box.z1; ++z) {
      if (!mask.isFluid(x, y, z) || visited[planeIndex(y, z)] != 0) {
        continue;
      }

      components.emplace_back();
      std::queue<std::pair<plint, plint>> pending;
      pending.push(std::make_pair(y, z));
      visited[planeIndex(y, z)] = 1;

      while (!pending.empty()) {
        const plint currentY = pending.front().first;
        const plint currentZ = pending.front().second;
        pending.pop();
        components.back().push_back(Dot3D(x, currentY, currentZ));

        for (plint dy = -1; dy <= 1; ++dy) {
          for (plint dz = -1; dz <= 1; ++dz) {
            if (dy == 0 && dz == 0) {
              continue;
            }
            const plint nextY = currentY + dy;
            const plint nextZ = currentZ + dz;
            if (!mask.isFluid(x, nextY, nextZ) ||
                visited[planeIndex(nextY, nextZ)] != 0) {
              continue;
            }
            visited[planeIndex(nextY, nextZ)] = 1;
            pending.push(std::make_pair(nextY, nextZ));
          }
        }
      }
    }
  }
  return components;
}

bool continuesInPositiveX(FluidMask const& mask,
                          std::vector<Dot3D> const& component) {
  if (component.empty() || component.front().x >= mask.box.x1) {
    return false;
  }

  const plint nextX = component.front().x + 1;
  for (Dot3D const& point : component) {
    for (plint dy = -1; dy <= 1; ++dy) {
      for (plint dz = -1; dz <= 1; ++dz) {
        if (mask.isFluid(nextX, point.y + dy, point.z + dz)) {
          return true;
        }
      }
    }
  }
  return false;
}

std::vector<OutletSection> identifyDaughterOutlets(FluidMask const& mask) {
  std::vector<OutletSection> outlets;

  for (plint x = mask.box.x1; x >= mask.box.x0 && outlets.size() < 2;
       --x) {
    std::vector<std::vector<Dot3D>> components =
        connectedFluidComponents(mask, x);

    for (std::vector<Dot3D> const& component : components) {
      if (continuesInPositiveX(mask, component)) {
        continue;
      }

      OutletSection outlet;
      outlet.x = x;
      outlet.points = component;
      for (Dot3D const& point : outlet.points) {
        outlet.centerY += point.y;
        outlet.centerZ += point.z;
      }
      outlet.centerY /= static_cast<T>(outlet.points.size());
      outlet.centerZ /= static_cast<T>(outlet.points.size());
      outlets.push_back(outlet);
    }
  }

  std::sort(outlets.begin(), outlets.end(),
            [](OutletSection const& lhs, OutletSection const& rhs) {
              return lhs.centerY < rhs.centerY;
            });
  return outlets;
}

std::vector<Dot3D> fluidPointsInBox(FluidMask const& mask,
                                    Box3D const& domain) {
  std::vector<Dot3D> points;
  for (plint x = domain.x0; x <= domain.x1; ++x) {
    for (plint y = domain.y0; y <= domain.y1; ++y) {
      for (plint z = domain.z0; z <= domain.z1; ++z) {
        if (mask.isFluid(x, y, z)) {
          points.push_back(Dot3D(x, y, z));
        }
      }
    }
  }
  return points;
}

T axialFlowRate(MultiBlockLattice3D<T, DESCRIPTOR>& lattice,
                std::vector<Dot3D> const& section) {
  T flowRate = 0;
  plb::Array<T, 3> velocity;
  for (Dot3D const& point : section) {
    lattice.get(point.x, point.y, point.z).computeVelocity(velocity);
    flowRate += velocity[0];
  }
  return flowRate;
}

}  // namespace

int main(int argc, char* argv[]) {
  if (argc < 2) {
    cout << "Usage: " << argv[0] << " <configuration.xml>" << endl;
    return EXIT_FAILURE;
  }

  HemoCell hemocell(argv[1], argc, argv);
  Config* cfg = hemocell.cfg;

  const std::array<int, 3> preInletBlocks = {{
      (*cfg)["preInlet"]["parameters"]["pABx"].read<int>(),
      (*cfg)["preInlet"]["parameters"]["pABy"].read<int>(),
      (*cfg)["preInlet"]["parameters"]["pABz"].read<int>()}};
  const std::array<int, 3> mainBlocks = {{
      (*cfg)["domain"]["mABx"].read<int>(),
      (*cfg)["domain"]["mABy"].read<int>(),
      (*cfg)["domain"]["mABz"].read<int>()}};
  if (std::any_of(preInletBlocks.begin(), preInletBlocks.end(),
                  [](int blocks) { return blocks <= 0; }) ||
      std::any_of(mainBlocks.begin(), mainBlocks.end(),
                  [](int blocks) { return blocks <= 0; })) {
    pcout << "(Bifurcation) All pAB and mAB decomposition counts must be "
             "positive."
          << endl;
    return EXIT_FAILURE;
  }

  const int preInletRanks =
      preInletBlocks[0] * preInletBlocks[1] * preInletBlocks[2];
  const int mainRanks = mainBlocks[0] * mainBlocks[1] * mainBlocks[2];
  const int requiredRanks = preInletRanks + mainRanks;
  if (global::mpi().getSize() != requiredRanks) {
    pcout << "(Bifurcation) Configured MPI decomposition requires "
          << requiredRanks << " ranks (preInlet " << preInletRanks
          << " + main domain " << mainRanks << "), but the run has "
          << global::mpi().getSize() << "." << endl;
    return EXIT_FAILURE;
  }
  pcout << "(Bifurcation) MPI decomposition: preInlet="
        << preInletBlocks[0] << "x" << preInletBlocks[1] << "x"
        << preInletBlocks[2] << " (" << preInletRanks
        << " ranks), main domain=" << mainBlocks[0] << "x"
        << mainBlocks[1] << "x" << mainBlocks[2] << " (" << mainRanks
        << " ranks)." << endl;

  hlog << "(Bifurcation) Reading and voxelizing STL file "
       << (*cfg)["domain"]["geometry"].read<string>() << endl;
  MultiScalarField3D<int>* flagMatrix = nullptr;
  VoxelizedDomain3D<double>* voxelizedDomain = nullptr;
  getFlagMatrixFromSTL((*cfg)["domain"]["geometry"].read<string>(),
                       (*cfg)["domain"]["fluidEnvelope"].read<int>(),
                       (*cfg)["domain"]["refDirN"].read<int>(),
                       (*cfg)["domain"]["refDir"].read<int>(),
                       voxelizedDomain, flagMatrix,
                       (*cfg)["domain"]["blockSize"].read<int>(),
                       (*cfg)["domain"]["particleEnvelope"].read<int>());

  param::lbm_base_parameters(*cfg);
  param::printParameters();

  hemocell.preInlet = new hemo::PreInlet(&hemocell, flagMatrix);
  Box3D inletSlice = flagMatrix->getBoundingBox();
  inletSlice.x1 = inletSlice.x0;
  hemocell.preInlet->preInletFromSlice(Direction::Xneg, inletSlice);

  hemocell.initializeLattice(voxelizedDomain->getMultiBlockManagement());
  if (!hemocell.partOfpreInlet) {
    hemocell.lattice->periodicity().toggleAll(false);
  }

  hemocell.preInlet->initializePreInlet();
  boundaryFromFlagMatrix(hemocell.lattice, flagMatrix,
                         hemocell.partOfpreInlet);
  hemocell.preInlet->createBoundary();

  FluidMask const fluidMask = collectGlobalFluidMask(*flagMatrix);
  std::vector<OutletSection> const outlets =
      identifyDaughterOutlets(fluidMask);
  if (outlets.size() != 2) {
    pcout << "(Bifurcation) Expected exactly two +x terminal outlet "
             "components, found "
          << outlets.size() << "." << endl;
    return EXIT_FAILURE;
  }

  const T outletDensity = (*cfg)["flow"]["outletDensity"].read<T>();
  OnLatticeBoundaryCondition3D<T, DESCRIPTOR>* outletBoundary =
      createZouHeBoundaryCondition3D<T, DESCRIPTOR>();
  for (pluint outletId = 0; outletId < outlets.size(); ++outletId) {
    OutletSection const& outlet = outlets[outletId];
    pcout << "(Bifurcation) Outlet " << outletId + 1 << ": x="
          << outlet.x << ", nodes=" << outlet.points.size()
          << ", center(y,z)=(" << outlet.centerY << ","
          << outlet.centerZ << "), rho=" << outletDensity << endl;

    for (Dot3D const& point : outlet.points) {
      Box3D const node(point.x, point.x, point.y, point.y, point.z,
                       point.z);
      outletBoundary->addPressureBoundary0P(
          node, *hemocell.domain_lattice, boundary::density);
      setBoundaryDensity(*hemocell.domain_lattice, node, outletDensity);
    }
  }

  std::vector<Dot3D> const inletPoints =
      fluidPointsInBox(fluidMask, hemocell.preInlet->fluidInlet);
  if (inletPoints.empty()) {
    pcout << "(Bifurcation) Main-domain inlet contains no fluid nodes."
          << endl;
    delete outletBoundary;
    return EXIT_FAILURE;
  }
  pcout << "(Bifurcation) Main inlet: x="
        << hemocell.preInlet->fluidInlet.x0
        << ", nodes=" << inletPoints.size() << endl;

  hemocell.lattice->toggleInternalStatistics(false);
  hemocell.latticeEquilibrium(1., plb::Array<T, 3>(0., 0., 0.));

  // The main bifurcation domain is not body-force driven.
  setExternalVector(*hemocell.domain_lattice,
                    hemocell.domain_lattice->getBoundingBox(),
                    DESCRIPTOR<T>::ExternalField::forceBeginsAt,
                    plb::Array<T, 3>(0., 0., 0.));

  hemocell.preInlet->calculateDrivingForce();
  hemocell.preInlet->setDrivingForce();
  pcout << "(Bifurcation) Periodic pre-inlet driving force = "
        << hemocell.preInlet->drivingForce
        << "; main-domain external force = (0,0,0)." << endl;

  hemocell.lattice->initialize();
  hemocell.initializeCellfield();

  hemocell.addCellType<RbcHighOrderModel>("RBC", RBC_FROM_SPHERE);
  hemocell.setMaterialTimeScaleSeparation(
      "RBC", (*cfg)["ibm"]["stepMaterialEvery"].read<int>());
  hemocell.setInitialMinimumDistanceFromSolid(
      "RBC",
      (*cfg)["ibm"]["initialMinimumDistanceFromSolid"].read<T>());
  hemocell.setParticleVelocityUpdateTimeScaleSeparation(
      (*cfg)["ibm"]["stepParticleEvery"].read<int>());

  const unsigned int interactionEvery =
      (*cfg)["ibm"]["stepInteractionEvery"].read<unsigned int>();
  if (interactionEvery == 0) {
    pcout << "(Bifurcation) stepInteractionEvery must be greater than zero."
          << endl;
    delete outletBoundary;
    return EXIT_FAILURE;
  }

  const T cellCellR0 = (*cfg)["cellCellAdhesion"]["r0"].read<T>();
  const T cellCellRc = (*cfg)["cellCellAdhesion"]["rc"].read<T>();
  const T cellCellEpsilon =
      (*cfg)["cellCellAdhesion"]["epsilon"].read<T>();
  const T cellCellD0 = (*cfg)["cellCellAdhesion"]["D0"].read<T>();
  const T cellCellAlpha =
      (*cfg)["cellCellAdhesion"]["alpha"].read<T>();
  const T kRep = (*cfg)["boundaryRepulsion"]["kRep"].read<T>();
  const T RepCutoff =
      (*cfg)["boundaryRepulsion"]["RepCutoff"].read<T>();

  hemocell.setAdhesion(cellCellR0, cellCellRc, cellCellEpsilon, cellCellD0,
                       cellCellAlpha);
  hemocell.setAdhesionTimeScaleSeperation(interactionEvery);
  pcout << "(Bifurcation) RBC-RBC adhesion and RBC-wall repulsion update every "
        << interactionEvery << " iteration(s)." << endl;

  // OUTPUT_FORCE_REPULSION contains the sum of RBC-RBC adhesion and original
  // solid-wall repulsion.
  vector<int> cellOutputs = {
      OUTPUT_POSITION,      OUTPUT_TRIANGLES,    OUTPUT_VELOCITY,
      OUTPUT_FORCE,         OUTPUT_FORCE_VOLUME, OUTPUT_FORCE_BENDING,
      OUTPUT_FORCE_LINK,    OUTPUT_FORCE_AREA,   OUTPUT_FORCE_VISC,
      OUTPUT_FORCE_REPULSION, OUTPUT_CELL_ID,     OUTPUT_VERTEX_ID,
      OUTPUT_RES_TIME};
  hemocell.setOutputs("RBC", cellOutputs);
  hemocell.setFluidOutputs(
      {OUTPUT_VELOCITY, OUTPUT_DENSITY, OUTPUT_FORCE, OUTPUT_BOUNDARY});

  if (!cfg->checkpointed) {
    std::ifstream rbcPositions("RBC.pos");
    if (!rbcPositions.good()) {
      pcout << "(Bifurcation) Fresh run requires the existing RBC.pos input "
               "in the case directory."
            << endl;
      delete outletBoundary;
      return EXIT_FAILURE;
    }
    pcout << "(Bifurcation) Fresh run uses the existing RBC.pos directly; "
             "the case does not generate, rewrite, or adjust that file."
          << endl;

    const plint warmup = (*cfg)["parameters"]["warmup"].read<plint>();
    pcout << "(Bifurcation) Warming up the cell-free flow for " << warmup
          << " iterations." << endl;
    for (plint warmupIter = 0; warmupIter < warmup; ++warmupIter) {
      hemocell.preInlet->setDrivingForce();
      hemocell.lattice->collideAndStream();
      hemocell.preInlet->applyPreInletVelocityBoundary();
    }

    if (!configureSolidBoundaryRepulsion(
            hemocell, kRep, RepCutoff,
            interactionEvery, inletPoints, outlets)) {
      delete outletBoundary;
      return EXIT_FAILURE;
    }
    hemocell.loadParticles();
    hemocell.writeOutput();
  } else {
    pcout << "(Bifurcation) Checkpoint restart restores particle fields from "
             "the checkpoint; RBC.pos is not read."
          << endl;
    hemocell.loadCheckPoint();
    if (!configureSolidBoundaryRepulsion(
            hemocell, kRep, RepCutoff,
            interactionEvery, inletPoints, outlets)) {
      delete outletBoundary;
      return EXIT_FAILURE;
    }
    hemocell.preInlet->setDrivingForce();
  }

  const unsigned int tmax =
      (*cfg)["sim"]["tmax"].read<unsigned int>();
  const unsigned int tmeas =
      (*cfg)["sim"]["tmeas"].read<unsigned int>();
  const unsigned int tcheckpoint =
      (*cfg)["sim"]["tcheckpoint"].read<unsigned int>();

  pcout << "(Bifurcation) Starting pure-RBC flow with RBC-RBC adhesion and "
           "solid-wall repulsion."
        << endl;
  while (hemocell.iter < tmax) {
    // Install the force before the pre-inlet collide-and-stream performed by
    // HemoCell::iterate(). The main-domain external force remains zero.
    hemocell.preInlet->setDrivingForce();
    hemocell.iterate();
    hemocell.preInlet->applyPreInlet();

    if (tmeas > 0 && hemocell.iter % tmeas == 0) {
      const T inletFlow =
          axialFlowRate(*hemocell.domain_lattice, inletPoints);
      const T outletFlow1 =
          axialFlowRate(*hemocell.domain_lattice, outlets[0].points);
      const T outletFlow2 =
          axialFlowRate(*hemocell.domain_lattice, outlets[1].points);
      const T massResidual = inletFlow - outletFlow1 - outletFlow2;
      FluidStatistics const fluidVelocity =
          FluidInfo::calculateVelocityStatistics(&hemocell);
      const pluint rbcCount =
          CellInformationFunctionals::getNumberOfCellsFromType(&hemocell,
                                                                "RBC");

      pcout << "(Bifurcation) iter=" << hemocell.iter
            << ", time=" << hemocell.iter * param::dt
            << " s, RBC=" << rbcCount << ", Qin=" << inletFlow
            << ", Qout1=" << outletFlow1 << ", Qout2=" << outletFlow2
            << ", flow residual=" << massResidual
            << ", max|u|=" << fluidVelocity.max << endl;
      hemocell.writeOutput();
    }

    if (tcheckpoint > 0 && hemocell.iter % tcheckpoint == 0) {
      hemocell.saveCheckPoint();
    }
  }

  pcout << "(Bifurcation) Simulation finished." << endl;
  delete outletBoundary;
  return EXIT_SUCCESS;
}
