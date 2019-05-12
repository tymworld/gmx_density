//
// Created by Yiming Tang on 2019-05-12.
//

#include "density.h"
#include <iostream>
#include <math>
#include <gromacs/trajectoryanalysis.h>

using namespace std;
using namespace gmx;

density::density(): hist_min_density_(0.0), hist_max_density_(200.0), mesh_number_(20)
{
    registerAnalysisDataset(&dataHistogram_, "histogram");
}

void density::initOptions(IOptionsContainer *options, TrajectoryAnalysisSettings *settings)
{
    static const char *const desc[] =
            {
            "Density Profile Analyzer by Yiming Tang @ Fudan.\n",
            "Introductions are still on its way. Bite me!"
            };

    settings->setHelpText(desc);

    options->addOption(FileNameOption("hist")
                               .filetype(eftPlot).outputFile()
                               .store(&fnHistogram_).defaultBasename("histogram")
                               .description("Histogram of density profiles"));

    options->addOption(DoubleOption("histmin").store(&hist_min_density_)
                               .defaultValue(0.0)
                               .description("Minimum density for histogram"));

    options->addOption(DoubleOption("histmax").store(&hist_max_density_)
                               .defaultValue(200.0)
                               .description("Maximum density for histogram"));

    options->addOption(DoubleOption("histstep").store(&hist_max_density_)
                               .defaultValue(1.0)
                               .description("Step density for histogram"));

    options->addOption(DoubleOption("mesh").store(&mesh_number_)
                               .defaultValue(20)
                               .description("Number of meshes to use for each direction"));

    options->addOption(SelectionOption("select").store(&sel_).required()
                               .description("Group contains molecules to be analyzed"));

    settings->setFlag(TrajectoryAnalysisSettings::efRequireTop);

}

void density::initAnalysis(const TrajectoryAnalysisSettings &settings, const TopologyInformation &top)
{
    dataHistogram_.setColumnCount(0, int((hist_max_density_ - hist_min_density_) / hist_step_density_) + 1);

    avegHistogram_.reset(new AnalysisDataAverageModule());
    dataHistogram_.addModule(avegHistogram_);

    this->top_ = top.topology();
    this->atoms_ = top.topology()->atoms;

}

void density::analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc, TrajectoryAnalysisModuleData *pdata)
{
    AnalysisDataHandle dhHistogram = pdata->dataHandle(dataHistogram_);
    dhHistogram.startFrame(frnr, fr.time);
    const Selection &sel = sel_;
    int *tempCount = new int[int(pow(mesh_number_, 3))];

    




}

