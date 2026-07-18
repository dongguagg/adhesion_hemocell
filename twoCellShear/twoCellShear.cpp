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
#include "helper/cellInfo.h"

#include "palabos3D.h"
#include "palabos3D.hh"

#include <fstream>
#include <map>
#include <string>

using namespace hemo;

namespace {

T topWallVelocity(unsigned int iteration, unsigned int tRelax,
                  unsigned int tRamp, T targetVelocity) {
  if (iteration < tRelax) {
    return 0.0;
  }
  if (tRamp == 0 || iteration >= tRelax + tRamp) {
    return targetVelocity;
  }
  return targetVelocity * static_cast<T>(iteration - tRelax) /
         static_cast<T>(tRamp);
}

const char *simulationStage(unsigned int iteration, unsigned int tRelax,
                            unsigned int tRamp) {
  if (iteration < tRelax) {
    return "relaxation";
  }
  if (tRamp > 0 && iteration < tRelax + tRamp) {
    return "ramp";
  }
  return "shear";
}

void setWallVelocities(HemoCell &hemocell, const Box3D &bottom,
                       const Box3D &top, T topVelocity) {
  setBoundaryVelocity(*hemocell.lattice, bottom,
                      plb::Array<T, 3>(0.0, 0.0, 0.0));
  setBoundaryVelocity(*hemocell.lattice, top,
                      plb::Array<T, 3>(topVelocity, 0.0, 0.0));
}

void writeCellState(HemoCell &hemocell, plb_ofstream &stream,
                    unsigned int tRelax, unsigned int tRamp,
                    T topVelocity) {
  std::map<int, CellInformation> cellInfo;
  CellInformationFunctionals::calculateCellInformation(&hemocell, cellInfo);

  if (!global::mpi().isMainProcessor()) {
    return;
  }

  const T micrometerToLbm = 1.0e-6 / param::dx;
  const T velocityToSI = param::dx / param::dt;
  hemo::Array<T, 3> relativeCenter({0.0, 0.0, 0.0});

  const std::map<int, CellInformation>::const_iterator lower = cellInfo.find(0);
  const std::map<int, CellInformation>::const_iterator upper = cellInfo.find(1);
  if (lower != cellInfo.end() && upper != cellInfo.end()) {
    relativeCenter = (upper->second.position - lower->second.position) /
                     micrometerToLbm;
  }

  for (int cellId = 0; cellId <= 1; ++cellId) {
    const std::map<int, CellInformation>::const_iterator cell =
        cellInfo.find(cellId);
    if (cell == cellInfo.end()) {
      continue;
    }

    const CellInformation &info = cell->second;
    const hemo::Array<T, 3> position = info.position / micrometerToLbm;
    const hemo::Array<T, 3> velocity = info.velocity * velocityToSI;
    const hemo::Array<T, 6> bbox = info.bbox / micrometerToLbm;
    const T area = info.area / (micrometerToLbm * micrometerToLbm);
    const T volume = info.volume /
                     (micrometerToLbm * micrometerToLbm * micrometerToLbm);

    stream << hemocell.iter << ','
           << simulationStage(hemocell.iter, tRelax, tRamp) << ',' << cellId
           << ','
           << (cellId == 0 ? "lower_wall_adhering" : "upper") << ','
           << topVelocity
           << ',' << topVelocity * param::dx / param::dt << ',' << position[0]
           << ',' << position[1] << ',' << position[2] << ',' << velocity[0]
           << ',' << velocity[1] << ',' << velocity[2] << ',' << bbox[0] << ','
           << bbox[1] << ',' << bbox[2] << ',' << bbox[3] << ',' << bbox[4]
           << ',' << bbox[5] << ',' << area << ',' << volume << ','
           << relativeCenter[0] << ',' << relativeCenter[1] << ','
           << relativeCenter[2] << std::endl;
  }
}

} // namespace

int main(int argc, char *argv[]) {
  if (argc < 2) {
    cout << "Usage: " << argv[0] << " <configuration.xml>" << endl;
    return -1;
  }

  HemoCell hemocell(argv[1], argc, argv);
  Config *cfg = hemocell.cfg;

  if (global::mpi().getSize() != 2) {
    pcerr << "(TwoCellShear) This case must be run with exactly two MPI ranks."
          << endl;
    return -1;
  }

  const T dx = (*cfg)["domain"]["dx"].read<T>();
  const T lx = (*cfg)["domain"]["Lx"].read<T>();
  const T ly = (*cfg)["domain"]["Ly"].read<T>();
  const T lz = (*cfg)["domain"]["Lz"].read<T>();
  const plint nx = static_cast<plint>(lx * 1.0e-6 / dx + 0.5);
  const plint ny = static_cast<plint>(ly * 1.0e-6 / dx + 0.5);
  const plint nz = static_cast<plint>(lz * 1.0e-6 / dx + 0.5);
  const plint fluidEnvelope =
      (*cfg)["domain"]["fluidEnvelope"].read<plint>();

  pcout << "(TwoCellShear) Calculating shear-flow parameters" << endl;
  param::lbm_shear_parameters(*cfg, nz);
  param::printParameters();

  pcout << "(TwoCellShear) Initializing " << nx << 'x' << ny << 'x' << nz
        << " lattice with a 1x1x2 block decomposition" << endl;
  const SparseBlockStructure3D blockStructure =
      createRegularDistribution3D(nx, ny, nz, 1, 1, 2);
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
  setWallVelocities(hemocell, bottom, top, 0.0);

  hemocell.latticeEquilibrium(1.0, plb::Array<T, 3>(0.0, 0.0, 0.0));
  setExternalVector(*hemocell.lattice, hemocell.lattice->getBoundingBox(),
                    DESCRIPTOR<T>::ExternalField::forceBeginsAt,
                    plb::Array<T, DESCRIPTOR<T>::d>(0.0, 0.0, 0.0));
  hemocell.lattice->initialize();
  hemocell.outputInSiUnits = true;

  hemocell.initializeCellfield();
  hemocell.addCellType<RbcHighOrderModel>("RBC", RBC_FROM_SPHERE);
  hemocell.setSystemPeriodicity(0, true);
  hemocell.setSystemPeriodicity(1, true);
  hemocell.setSystemPeriodicity(2, false);

  hemocell.setMaterialTimeScaleSeparation(
      "RBC", (*cfg)["ibm"]["stepMaterialEvery"].read<unsigned int>());
  hemocell.setParticleVelocityUpdateTimeScaleSeparation(
      (*cfg)["ibm"]["stepParticleEvery"].read<unsigned int>());
  hemocell.setAdhesionTimeScaleSeperation(1);

  const T cellCellR0 = (*cfg)["cellCellAdhesion"]["r0"].read<T>();
  const T cellCellRc = (*cfg)["cellCellAdhesion"]["rc"].read<T>();
  const T cellCellEpsilon =
      (*cfg)["cellCellAdhesion"]["epsilon"].read<T>();
  const T cellCellD0 = (*cfg)["cellCellAdhesion"]["D0"].read<T>();
  const T cellCellAlpha = (*cfg)["cellCellAdhesion"]["alpha"].read<T>();
  hemocell.setAdhesion(cellCellR0, cellCellRc, cellCellEpsilon, cellCellD0,
                       cellCellAlpha);

  const T cellWallR0 = (*cfg)["cellWallAdhesion"]["r0"].read<T>();
  const T cellWallRc = (*cfg)["cellWallAdhesion"]["rc"].read<T>();
  const T cellWallEpsilon =
      (*cfg)["cellWallAdhesion"]["epsilon"].read<T>();
  const T cellWallD0 = (*cfg)["cellWallAdhesion"]["D0"].read<T>();
  const T cellWallAlpha = (*cfg)["cellWallAdhesion"]["alpha"].read<T>();
  hemocell.setBoundaryAdhesion(cellWallR0, cellWallRc, cellWallEpsilon,
                               cellWallD0, cellWallAlpha);

  const T maximumAttractionRatio =
      cellWallAlpha * cellWallD0 / (cellCellAlpha * cellCellD0);
  pcout << "(TwoCellShear) Cell-wall maximum pair attraction is "
        << maximumAttractionRatio
        << " times the cell-cell maximum pair attraction" << endl;

  hemocell.setOutputs(
      "RBC", {OUTPUT_POSITION, OUTPUT_TRIANGLES, OUTPUT_VELOCITY, OUTPUT_FORCE,
              OUTPUT_FORCE_REPULSION, OUTPUT_FORCE_VOLUME,
              OUTPUT_FORCE_BENDING, OUTPUT_FORCE_LINK, OUTPUT_FORCE_AREA,
              OUTPUT_FORCE_VISC});
  hemocell.setFluidOutputs({OUTPUT_VELOCITY, OUTPUT_FORCE});

  const unsigned int warmup =
      (*cfg)["parameters"]["warmup"].read<unsigned int>();
  const unsigned int tRelax = (*cfg)["sim"]["tRelax"].read<unsigned int>();
  const unsigned int tRamp = (*cfg)["sim"]["tRamp"].read<unsigned int>();
  const unsigned int tmax = (*cfg)["sim"]["tmax"].read<unsigned int>();
  const unsigned int tmeas = (*cfg)["sim"]["tmeas"].read<unsigned int>();
  const unsigned int tcheckpoint =
      (*cfg)["sim"]["tcheckpoint"].read<unsigned int>();

  const T targetTopVelocity = param::shearrate_lbm * (nz - 1);
  pcout << "(TwoCellShear) Target shear rate: "
        << (*cfg)["domain"]["shearrate"].read<T>() << " s^-1" << endl;
  pcout << "(TwoCellShear) Discrete wall distance: "
        << (nz - 1) * param::dx * 1.0e6 << " um; target top-wall velocity: "
        << targetTopVelocity << " lu ("
        << targetTopVelocity * param::dx / param::dt << " m/s)" << endl;

  plb_ofstream cellLog;
  if (cfg->checkpointed) {
    // Loading the checkpoint restores its output-directory configuration.
    pcout << "(TwoCellShear) Loading checkpoint" << endl;
    hemocell.loadCheckPoint();
    const std::string cellLogName =
        global::directories().getOutputDir() + "twoCellShear.csv";
    cellLog.open(cellLogName.c_str(), std::ofstream::out | std::ofstream::app);
  } else {
    const std::string cellLogName =
        global::directories().getOutputDir() + "twoCellShear.csv";
    cellLog.open(cellLogName.c_str());
    cellLog << "iteration,stage,cell_id,role,top_wall_velocity_lbm,"
               "top_wall_velocity_m_per_s,center_x_um,center_y_um,center_z_um,"
               "velocity_x_m_per_s,velocity_y_m_per_s,velocity_z_m_per_s,"
               "bbox_x_min_um,bbox_x_max_um,bbox_y_min_um,bbox_y_max_um,"
               "bbox_z_min_um,bbox_z_max_um,area_um2,volume_um3,"
               "relative_center_x_um,relative_center_y_um,"
               "relative_center_z_um"
            << std::endl;

    pcout << "(TwoCellShear) Warming up the cell-free fluid for " << warmup
          << " iterations" << endl;
    for (unsigned int i = 0; i < warmup; ++i) {
      hemocell.lattice->collideAndStream();
    }

    hemocell.loadParticles();
    hemocell.writeOutput();
    writeCellState(hemocell, cellLog, tRelax, tRamp, 0.0);
  }

  while (hemocell.iter < tmax) {
    const T currentTopVelocity =
        topWallVelocity(hemocell.iter, tRelax, tRamp, targetTopVelocity);
    setWallVelocities(hemocell, bottom, top, currentTopVelocity);
    hemocell.iterate();

    if (hemocell.iter % tmeas == 0) {
      hemocell.writeOutput();
      writeCellState(hemocell, cellLog, tRelax, tRamp, currentTopVelocity);
    }
    if (hemocell.iter % tcheckpoint == 0) {
      hemocell.saveCheckPoint();
    }
  }

  cellLog.close();
  delete boundaryCondition;
  pcout << "(TwoCellShear) Simulation finished" << endl;
  return 0;
}
