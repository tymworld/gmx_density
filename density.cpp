//
// Created by Yiming Tang on 2019-05-12.
//

#include "density.h"
#include <iostream>
#include <cmath>
#include <gromacs/trajectoryanalysis.h>
#include <gromacs/fileio/pdbio.h>
#include <gromacs/pbcutil/pbc.h>
#include <gromacs/pbcutil/rmpbc.h>
#include <gromacs/selection/nbsearch.h>
#include <vector>
#include <numeric>

using namespace std;
using namespace gmx;

density::density(): hist_min_density_(0.0), hist_max_density_(200.0), mesh_number_(20)
{
    registerAnalysisDataset(&dataHistogram_, "histogram");
    registerAnalysisDataset(&dataMaxMinDensity_, "maxmin-density");
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
                               .description("Histogram of density profiles, averaged"));

    options->addOption(FileNameOption("maxmin-density")
                               .filetype(eftPlot).outputFile()
                               .store(&fnMaxMinDensity_).defaultBasename("maxmin_density")
                               .description("Maximum/Minimum density"));

    options->addOption(DoubleOption("histmin").store(&hist_min_density_)
                               .defaultValue(0.0)
                               .description("Minimum density for histogram"));

    options->addOption(DoubleOption("histmax").store(&hist_max_density_)
                               .defaultValue(200.0)
                               .description("Maximum density for histogram"));

    options->addOption(DoubleOption("histstep").store(&hist_step_density_)
                               .defaultValue(1.0)
                               .description("Step density for histogram"));

    options->addOption(DoubleOption("rprobe").store(&probeRadius_)
                               .defaultValue(2.0)
                               .description("Radius of density probes for max/min density (nm)"));

    options->addOption(DoubleOption("percentage-max-min").store(&percentage_maxmin_)
                               .defaultValue(0.1)
                               .description("Percentage of the largest/smallest density values to calculate max/min value"));

    options->addOption(IntegerOption("probemesh").store(&probeMesh_)
                               .defaultValue(100)
                               .description("Number of meshes to use for each direction for max/min density"));

    options->addOption(IntegerOption("mesh").store(&mesh_number_)
                               .defaultValue(20)
                               .description("Number of meshes to use for each direction"));

    options->addOption(SelectionOption("select").store(&sel_).required()
                               .description("Group contains molecules to be analyzed"));

    settings->setFlag(TrajectoryAnalysisSettings::efRequireTop);


}

void density::initAnalysis(const TrajectoryAnalysisSettings &settings, const TopologyInformation &top)
{
    //dataHistogram_.setColumnCount(0, int((hist_max_density_ - hist_min_density_) / hist_step_density_) + 1);


    if(!fnHistogram_.empty())
    {
        dataHistogram_.setColumnCount(0, (int)((hist_max_density_ - hist_min_density_)/hist_step_density_ ) + 1);

        this->top_ = top.topology();
        this->atoms_ = top.topology()->atoms;

        avemHistogram_.reset(new AnalysisDataAverageModule());
        dataHistogram_.addModule(avemHistogram_);
    }

    if(!fnMaxMinDensity_.empty())
    {
        dataMaxMinDensity_.setColumnCount(0, 2);

        AnalysisDataPlotModulePointer plotMaxMin(new AnalysisDataPlotModule(settings.plotSettings()));
        plotMaxMin->setFileName(fnMaxMinDensity_);
        plotMaxMin->setTitle("Max/Min Density");
        plotMaxMin->setYLabel("Density (mg/mL)");
        plotMaxMin->setXAxisIsTime();

        plotMaxMin->appendLegend("Maximum Density");
        plotMaxMin->appendLegend("Minimum Density");

        dataMaxMinDensity_.addModule(plotMaxMin);

        nb_.setCutoff(probeRadius_);

    }


    top_ = top.topology();
    atoms_ = top.topology()->atoms;



}

void density::analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc, TrajectoryAnalysisModuleData *pdata)
{
    if(!fnHistogram_.empty()) {
        AnalysisDataHandle dhHistogram = pdata->dataHandle(dataHistogram_);
        dhHistogram.startFrame(frnr, fr.time);
        const Selection &sel = sel_;

        //unsigned *massCount_;
        //int *densityCount_;

        //massCount_ = new unsigned(pow(mesh_number_, 3));
        //densityCount_ = new int((hist_max_density_ - hist_min_density_)/hist_step_density_  + 1);


        for (int i = 0; i < pow(mesh_number_, 3); i++) {
            massCount_[i] = 0;
        }

        for (int i = 0; i < (hist_max_density_ - hist_min_density_) / hist_step_density_ + 1; i++) {
            densityCount_[i] = 0;
        }

        double box_x = pbc->box[0][0];
        double box_y = pbc->box[1][1];
        double box_z = pbc->box[2][2];

        // We now go through all atoms and add their masses to the mass matrix

        gmx::ArrayRef<float const[3]>::iterator iter_coordinate;
        gmx::ArrayRef<const real>::iterator iter_mass;

        //cout << "DEBUG: Begin to go through atoms" << endl;

        for (iter_coordinate = sel.coordinates().begin(), iter_mass = sel.masses().begin();
             iter_coordinate != sel.coordinates().end() && iter_mass != sel.masses().end();
             ++iter_coordinate, ++iter_mass) {

            int x_index = int(floor((*iter_coordinate)[0] / (box_x / mesh_number_)));
            int y_index = int(floor((*iter_coordinate)[1] / (box_y / mesh_number_)));
            int z_index = int(floor((*iter_coordinate)[2] / (box_z / mesh_number_)));

            if (x_index == mesh_number_) x_index -= 1;
            if (y_index == mesh_number_) y_index -= 1;
            if (z_index == mesh_number_) z_index -= 1;

            if (x_index == -1) x_index = mesh_number_ - 1;
            if (y_index == -1) y_index = mesh_number_ - 1;
            if (z_index == -1) z_index = mesh_number_ - 1;

            massCount_[x_index * (mesh_number_ * mesh_number_) + y_index * mesh_number_ +
                       z_index] += (int) (*iter_mass);
        }

        //cout << "DEBUG: Gone through all atoms." << endl;

        // We now go through elements of mass matrix, calculate their density and put them into right statistics.

        double density_temp;
        for (int x_index = 0; x_index < mesh_number_; x_index++) {
            for (int y_index = 0; y_index < mesh_number_; y_index++) {
                for (int z_index = 0; z_index < mesh_number_; z_index++) {

                    //cout << massCount_[x_index * (mesh_number_ * mesh_number_) + y_index * mesh_number_ + z_index] << '\t';
                    density_temp =
                            massCount_[x_index * (mesh_number_ * mesh_number_) + y_index * mesh_number_ + z_index]
                            / (box_x * box_y * box_z / pow(mesh_number_, 3)) * 10 / 6.02;
                    //cout << density_temp << endl;

                    densityCount_[(int) (floor((density_temp - hist_min_density_) / hist_step_density_))] += 1;
                }
            }
        }

        //cout << "DEBUG: Calculated the mass matrix." << endl;

        for (int i = 0; i < (hist_max_density_ - hist_min_density_) / hist_step_density_ + 1; i++) {
            dhHistogram.setPoint(i, densityCount_[i]);
        }

        dhHistogram.finishFrame();
    }




    if(!fnMaxMinDensity_.empty())
    {


        double tempMass;

        real max_density = 0.0;
        real min_density = 0.0;
        real density_temp;

        const Selection &sel = sel_;
        AnalysisDataHandle dhMaxMin = pdata->dataHandle(dataMaxMinDensity_);
        dhMaxMin.startFrame(frnr, fr.time);
        AnalysisNeighborhoodSearch nbsearch = nb_.initSearch(pbc, sel);

        int maxmin_point_number = (int)floor(pow(probeMesh_, 3) * percentage_maxmin_);
        /*
        vector<double> max_vector;
        vector<double> min_vector;
        */
        vector<float> density_vector;

        for(int x_index = 0; x_index < probeMesh_; x_index++)
        {
            for(int y_index = 0; y_index < probeMesh_; y_index++)
            {
                for(int z_index = 0; z_index < probeMesh_; z_index++)
                {
                    real mass_temp = 0;
                    float const searchPoint[3] = {(float)x_index / probeMesh_ * pbc->box[0][0],
                                                  (float)y_index / probeMesh_ * pbc->box[1][1],
                                                  (float)z_index / probeMesh_ * pbc->box[2][2]};
                    AnalysisNeighborhoodPairSearch pairSearch = nbsearch.startPairSearch(searchPoint);
                    AnalysisNeighborhoodPair pair;
                    while(pairSearch.findNextPair(&pair))
                    {
                        //cout << "DEBUG: REFINDEX is " << pair.refIndex() << ", ";
                        //cout << "Mass is " << atoms_.atom[pair.refIndex()].m << endl;
                        mass_temp += atoms_.atom[pair.refIndex()].m;
                    }

                    density_temp = (real)(mass_temp / (4.0/3 * M_PI * pow(probeRadius_,3) ) * 10 / 6.02);
                    /*

                    if(max_vector.size() < maxmin_point_number)
                    {
                        max_vector.insert(max_vector.end(), density_temp);
                        sort(max_vector.begin(), max_vector.end());
                    } else
                    {
                        if(density_temp > max_vector[0])
                        {
                            max_vector.erase(max_vector.begin());
                            max_vector.insert(max_vector.begin(), density_temp);
                            sort(max_vector.begin(), max_vector.end());
                        }
                    }

                    if(min_vector.size() < maxmin_point_number)
                    {
                        min_vector.insert(min_vector.end(), density_temp);
                        sort(min_vector.begin(), min_vector.end());
                    } else
                    {
                        if(density_temp < min_vector[maxmin_point_number - 1])
                        {
                            min_vector.pop_back();
                            min_vector.insert(min_vector.begin(), density_temp);
                            sort(min_vector.begin(), min_vector.end());
                        }
                    }

                    */

                    density_vector.insert(density_vector.end(), density_temp);


                    /*


                    if(x_index == 0 && y_index == 0 && z_index == 0)
                    {
                        max_density = density_temp;
                        min_density = density_temp;
                    } else
                    {
                        max_density = density_temp > max_density ? density_temp : max_density;
                        min_density = density_temp < min_density ? density_temp : min_density;
                    }
                    //cout << density_temp << endl;
                    */

                }
            }
        }
        //cout << "HAHA" << endl;
        /*

        max_density = accumulate(begin(max_vector), end(max_vector), 0.0) / max_vector.size();
        min_density = accumulate(begin(min_vector), end(min_vector), 0.0) / min_vector.size();
        */ 
        sort(density_vector.begin(), density_vector.end());

        long double min_density_sum = 0.0;
        long double max_density_sum = 0.0;

        for(int i = 0; i < maxmin_point_number; i++)
        {
            min_density_sum += density_vector[i];
            max_density_sum += density_vector[pow(probeMesh_, 3) - 1 - i];
        }

        min_density = min_density_sum / maxmin_point_number;
        max_density = max_density_sum / maxmin_point_number;

        dhMaxMin.setPoint(0, max_density);
        dhMaxMin.setPoint(1, min_density);

        dhMaxMin.finishFrame();


    }

}

void density::finishAnalysis(int) {


}

void density::writeOutput()
{
    if(!fnHistogram_.empty())
    {
        if(fnHistogram_.compare(".xvg")){}
        else
        {
            fnHistogram_ += ".dat";
        }

        FILE *fpHistogram = fopen(fnHistogram_.c_str(), "w");
        fprintf(fpHistogram, "# This file was created by gmx_density by Yiming Tang\n");
        fprintf(fpHistogram, "@    title \"Histogram\"\n");
        fprintf(fpHistogram, "@    xaxis  label \"Density(mg/mL)\"\n");
        fprintf(fpHistogram, "@    yaxis  label \"Probability\"\n");
        fprintf(fpHistogram, "@TYPE BAR\n");

        for(int i=0; i<dataHistogram_.columnCount(); i++)
        {
            fprintf(fpHistogram, "%f  %f\n", hist_min_density_ + hist_step_density_ * i,
                    avemHistogram_->average(0,i));
        }
        fclose(fpHistogram);
    }


}

