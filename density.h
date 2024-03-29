//
// Created by Yiming Tang on 2019-05-12.
//

#ifndef GMX_DENSITY_DENSITY_H
#define GMX_DENSITY_DENSITY_H

#include <gromacs/trajectoryanalysis.h>
#include <iostream>

using namespace std;
using namespace gmx;

class density : public gmx::TrajectoryAnalysisModule{
public:
    density();

    virtual void initOptions(IOptionsContainer *options,
                             TrajectoryAnalysisSettings *settings) override;

    virtual void initAnalysis(const TrajectoryAnalysisSettings &settings,
                              const TopologyInformation &top) override;

    virtual void analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc,
                              TrajectoryAnalysisModuleData *pdata) override;

    virtual void finishAnalysis(int nframes) override;

    virtual void writeOutput() override;

private:

    class ModuleData;

    // File names for output control

    // This file prints the histogram for density.
    // This file cut the simulation box into several pieces.
    std::string fnHistogram_;

    //// These parameters are for use of maximum/minimum density.

    std::string fnRawDensity_;

    // These file print the maximum / minimum density within spheres of radius defined by probeRadius_;
    std::string fnMaxMinDensity_;

    double probeRadius_;
    int probeMesh_;
    double percentage_maxmin_;

    AnalysisNeighborhood nb_;

    AnalysisData dataMaxMinDensity_;

    //// End of parameters for use of maximum/minimum density.

    // Minimum and maximum density and step values for histogram
    double hist_min_density_;
    double hist_max_density_;
    double hist_step_density_;

    // Mesh density in each direction
    int mesh_number_;

    Selection sel_;

    AnalysisData dataHistogram_;

    AnalysisDataAverageModulePointer avemHistogram_;

    unsigned massCount_[10000];
    int densityCount_[10000];

    //AnalysisDataAverageModulePointer avegHistogram_;

    t_topology *top_;
    t_atoms atoms_;
};

#endif //GMX_DENSITY_DENSITY_H
