#include <iostream>
#include "density.h"

int main(int argc, char *argv[]) {
    return gmx::TrajectoryAnalysisCommandLineRunner::runAsMain<density>(argc, argv);
}