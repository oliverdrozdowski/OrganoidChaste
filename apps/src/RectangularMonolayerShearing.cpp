/*

Copyright (c) 2005-2022, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

/**
 * @file
 *
 * This file gives an example of how you can create your own executable
 * in a user project.
 */

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include "FileFinder.hpp"

#include "SimulationTime.hpp"

#include "Exception.hpp"
#include "ExecutableSupport.hpp"
#include "OutputFileHandler.hpp"
#include "PetscException.hpp"
#include "PetscTools.hpp"

#include "CellApicalAreaWriter.hpp"
#include "CellBasalAreaWriter.hpp"
#include "CellCentroidWriter.hpp"
#include "CellLabelWriter.hpp"
#include "CellMutationStatesWriter.hpp"
#include "CellOpeningAngleWriter.hpp"
#include "CellProliferativeTypesWriter.hpp"
#include "CellSymmetryWriter.hpp"
#include "CellThicknessWriter.hpp"
#include "CellVolumesWriter.hpp"
#include "PopulationSurfaceTensionEnergyWriter.hpp"

#include "CellsGenerator.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "FaceTypeWriter.hpp"
#include "FiniteThicknessHoneycombVertexMeshGenerator.hpp"
#include "FiniteThicknessSimulation3d.hpp"
#include "MechanicalMosaicCellMutationState.hpp"
#include "MonolayerVertexBasedCellPopulation.hpp"
#include "NoCellCycleModel.hpp"
#include "OffLatticeSimulation.hpp"
#include "SmallAreaDampingModifier.hpp"
#include "SmartPointers.hpp"

#include "DecayingRandomForce.hpp"
#include "FixNodeDimensionsBoundaryCondition.hpp"
#include "GeneralizedVolumeConservingForce.hpp"
#include "GeometricalTargetVolumeModifier.hpp"
#include "RectangularEdgeSubForce.hpp"
#include "SimpleTargetVolumeModifier.hpp"
#include "SimplySupportedEdgeBoundaryCondition.hpp"
#include "SurfaceEvolverSaveModifier.hpp"
#include "SurfaceTensionSubForce.hpp"

struct paramsMonolayerStretch
{
    double GammaA = 1.0;
    double GammaB = 1.0;
    int NumCellsY = 20;
    int NumCellsX = 20;
    double StressVariable = 0.0;
    double MaxSimulationTime = 10000.0;
    double CenterForce = 0.0;
    double AreaDifferenceCutoff = 1e-5;
    unsigned PrintEveryNSteps = 1000;
    double MaxTimeLoopStepping = 10000.0;
    double FinalDT = 0.1;
    std::string FileComment = "";
};

// Function to parse command line arguments
paramsMonolayerStretch ParseArguments(int argc, char* argv[])
{
    paramsMonolayerStretch params;

    for (int i = 1; i < argc; i++)
    {
        std::string arg = argv[i];

        if (arg == "-h")
        {
            std::cout << "Usage: program [options]" << std::endl;
            std::cout << "Options:" << std::endl;
            std::cout << "  -GammaA <value>   : dim.less apical tension (std: 1.0)" << std::endl;
            std::cout << "  -GammaB <value>   : dim.less basal tension (std: 1.0)" << std::endl;
            std::cout << "  -NumX <value>   : Number hexagonal cells in x direction (std: 20)" << std::endl;
            std::cout << "  -NumY <value>   : Number hexagonal cells in y direction (std: 20)" << std::endl;
            std::cout << "  -StressVariable <value>   : Stress variable for given tensor (std: 0.0)" << std::endl;
            std::cout << "  -CenterForce <value>   : Upward Force per node for center cell in symm. breaking (std: 0.0)" << std::endl;
            std::cout << "  -MaxSimTime <value>   : MaxSimulationTime  (std: 1e4)" << std::endl;
            std::cout << "  -AreaDiffCutoff <value>   : Area difference cutoff to terminate (std: 1e-5)" << std::endl;
            std::cout << "  -PrintEveryNSteps <value>   : How often to save to file (std: 1000)" << std::endl;
            std::cout << "  -LoopStepping <value>   : How large the timestep will be increased if break condition not fulfilled (std: 1e4)" << std::endl;
            std::cout << "  -FinalDT <value>   : dt in final simulation loop (std: 0.1)" << std::endl;
            std::cout << "  -Comment <value>   : Comment to append to filename" << std::endl;
            EXCEPTION("Aborting as help flag was given.");
        }
        else if (arg == "-GammaA" && i + 1 < argc)
        {
            params.GammaA = std::stod(argv[i + 1]);
        }
        else if (arg == "-GammaB" && i + 1 < argc)
        {
            params.GammaB = std::stod(argv[i + 1]);
        }
        else if (arg == "-NumX" && i + 1 < argc)
        {
            params.NumCellsX = std::stod(argv[i + 1]);
        }
        else if (arg == "-NumY" && i + 1 < argc)
        {
            params.NumCellsY = std::stod(argv[i + 1]);
        }
        else if (arg == "-StressVariable" && i + 1 < argc)
        {
            params.StressVariable = std::stod(argv[i + 1]);
        }
        else if (arg == "-CenterForce" && i + 1 < argc)
        {
            params.CenterForce = std::stod(argv[i + 1]);
        }
        else if (arg == "-MaxSimTime" && i + 1 < argc)
        {
            params.MaxSimulationTime = std::stod(argv[i + 1]);
        }
        else if (arg == "-AreaDiffCutoff" && i + 1 < argc)
        {
            params.AreaDifferenceCutoff = std::stod(argv[i + 1]);
        }
        else if (arg == "-PrintEveryNSteps" && i + 1 < argc)
        {
            params.PrintEveryNSteps = std::stoi(argv[i + 1]);
        }
        else if (arg == "-LoopStepping" && i + 1 < argc)
        {
            params.MaxTimeLoopStepping = std::stod(argv[i + 1]);
        }
        else if (arg == "-FinalDT" && i + 1 < argc)
        {
            params.FinalDT = std::stod(argv[i + 1]);
        }
        else if (arg == "-Comment" && i + 1 < argc)
        {
            params.FileComment = argv[i + 1];
        }
        i++;
    }

    return params;
}

/*
 * Create the filepath where to save the files
 */
std::string CreateFilePath(const paramsMonolayerStretch& params)
{
    // Get the current date and time
    auto now = std::chrono::system_clock::now();
    std::time_t now_time = std::chrono::system_clock::to_time_t(now);

    // Convert the current date and time to the desired format
    std::tm* local_time = std::localtime(&now_time);
    std::ostringstream date_time;
    date_time << std::put_time(local_time, "%d%m%Y_%H_%M");

    std::ostringstream filepath;
    filepath << "Sheared_Honeycomb_"
             << "GA" << params.GammaA << "_"
             << "GB" << params.GammaB << "_"
             << "NX" << params.NumCellsX << "_"
             << "NY" << params.NumCellsY << "_"
             << "FN" << params.StressVariable << "_"
             << "ADC" << params.AreaDifferenceCutoff
             << params.FileComment << "_"
             << date_time.str();
    return filepath.str();
}

void SaveConstantsToFile(const std::string& filepath, const paramsMonolayerStretch& params,
                         double numNodesX, double numNodesY, double numNodesTop, double numNodesBottom)
{
    OutputFileHandler handler(filepath, false);
    out_stream output_file = handler.OpenOutputFile("MonolayerBucklingParameters.txt");

    if (output_file->is_open())
    {
        *output_file << "GammaA = " << params.GammaA << std::endl;
        *output_file << "GammaB = " << params.GammaB << std::endl;
        *output_file << "NumCellsY = " << params.NumCellsY << std::endl;
        *output_file << "NumCellsX = " << params.NumCellsX << std::endl;
        *output_file << "StressVariable = " << params.StressVariable << std::endl;
        *output_file << "CenterForce = " << params.CenterForce << std::endl;
        *output_file << "MaxSimulationTime = " << params.MaxSimulationTime << std::endl;
        *output_file << "AreaDifferenceCutoff = " << params.AreaDifferenceCutoff << std::endl;
        *output_file << "PrintEveryNSteps = " << params.PrintEveryNSteps << std::endl;
        *output_file << "MaxTimeLoopStepping = " << params.MaxTimeLoopStepping << std::endl;
        *output_file << "FinalDT = " << params.FinalDT << std::endl;
        *output_file << "NumNodesX = " << numNodesX << std::endl;
        *output_file << "NumNodesY = " << numNodesY << std::endl;
        *output_file << "NumNodesTop = " << numNodesTop << std::endl;
        *output_file << "NumNodesBottom = " << numNodesBottom << std::endl;

        output_file->close();
        ExecutableSupport::Print("\nConstants saved to file successfully.\n");
    }
    else
    {
        EXCEPTION("Error opening output file. Maybe filepath wrong?");
    }
}

double CalculateMaximalAreaDifference(
    MonolayerVertexBasedCellPopulation<3>* pCellPopulation, std::map<int, double>& rMapFaceAreas)
{
    // If map is empty, we initialize the map and return a negative value
    double max_area_diff = -1.0;
    bool initialize = rMapFaceAreas.empty();

    MonolayerVertexMesh<3, 3>& r_mesh = pCellPopulation->rGetMesh();
    // Loop over faces, because if we look at cell's faces, we go over the same faces multiple times
    unsigned num_faces = r_mesh.GetNumFaces();
    for (unsigned index_face = 0; index_face < num_faces; ++index_face)
    {
        MonolayerVertexElement<2, 3>* p_face = r_mesh.GetFace(index_face);
        double face_area_new = r_mesh.CalculateAreaOfFace(p_face);

        if (!initialize && abs(face_area_new - rMapFaceAreas[index_face]) > max_area_diff)
        {
            max_area_diff = abs(face_area_new - rMapFaceAreas[index_face]);
        }
        rMapFaceAreas[index_face] = face_area_new;
    }
    return max_area_diff;
}

int main(int argc, char* argv[])
{
    // This sets up PETSc and prints out copyright information, etc.
    ExecutableSupport::StandardStartup(&argc, &argv);

    int exit_code = ExecutableSupport::EXIT_OK;

    // You should put all the main code within a try-catch, to ensure that
    // you clean up PETSc before quitting.
    try
    {
        // Instantiate Simulation Time
        SimulationTime::Instance()->SetStartTime(0.0);

        // Parse the input and read the config file
        paramsMonolayerStretch simulation_parameters = ParseArguments(argc, argv);

        double MAX_TIME = simulation_parameters.MaxSimulationTime;
        double AREA_DIFFERENCE_CUTOFF = simulation_parameters.AreaDifferenceCutoff;

        double gamma_a = simulation_parameters.GammaA;
        double gamma_b = simulation_parameters.GammaB;

        double height = cbrt(2.0 / sqrt(3.0)) * cbrt((gamma_a + gamma_b) * (gamma_a + gamma_b));
        unsigned num_cells_x = simulation_parameters.NumCellsX;
        unsigned num_cells_y = simulation_parameters.NumCellsY;
        unsigned num_cells = num_cells_x * num_cells_y;

        // length at which T1s become favorable: l^*=0.66/3^(2/3)/gamma^(1/3)
        double t1_length = 0.66 / cbrt(3. * 3. * (gamma_a + gamma_b) / 2.);
        t1_length = 0.0;
        FiniteThicknessHoneycombVertexMeshGenerator generator(num_cells_x, num_cells_y, false,
                                                              t1_length, 0.001, height, 1.0 / height);
        // generator.CurveMonolayer(-2.0 * num_cells_x);
        MutableMonolayerVertexMesh<3, 3>* p_mesh = generator.GetMesh();
        num_cells = p_mesh->GetNumElements();
        // We use protorosettes as intermediate states
        p_mesh->SetProtorosetteFormationProbability(1.0);
        p_mesh->SetProtorosetteResolutionProbabilityPerTimestep(0.1);

        // Get the edge nodes
        c_vector<double, 3> top_force = zero_vector<double>(3);
        c_vector<double, 3> bottom_force = zero_vector<double>(3);

        // We generate the cell objects
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state_wt);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);

        for (unsigned cell_index = 0; cell_index < num_cells; cell_index++)
        {
            c_vector<double, 3> centr = p_mesh->GetCentroidOfElement(cell_index);
            NoCellCycleModel* p_model = new NoCellCycleModel;
            p_model->SetDimension(3);

            CellPtr p_cell;

            // Randomly assign cell state (i.e. tensions)
            p_cell = CellPtr(new Cell(p_state_wt, p_model));
            p_cell->SetCellProliferativeType(p_differentiated_type);
            double birth_time = -RandomNumberGenerator::Instance()->ranf() * 5.0;
            p_cell->SetBirthTime(birth_time);
            cells.push_back(p_cell);
        }

        MonolayerVertexBasedCellPopulation<3> cell_population(*p_mesh, cells);
        cell_population.SetOutputCellRearrangementLocations(false);
        cell_population.SetRestrictVertexMovementBoolean(false);
        cell_population.SetDoInitialVolumeRelaxation(true);

        cell_population.AddCellWriter<CellVolumesWriter>();
        cell_population.AddCellWriter<CellThicknessWriter>();
        cell_population.AddCellWriter<CellLabelWriter>();
        cell_population.AddCellWriter<CellMutationStatesWriter>();
        cell_population.AddCellWriter<CellCentroidWriter>();
        cell_population.AddCellWriter<CellApicalAreaWriter>();
        cell_population.AddCellWriter<CellBasalAreaWriter>();
        cell_population.AddCellWriter<CellOpeningAngleWriter>();
        cell_population.AddFaceWriter<FaceTypeWriter>();

        FiniteThicknessSimulation3d simulator(cell_population);
        simulator.SetOutputDirectory(CreateFilePath(simulation_parameters));
        simulator.SetEndTime(40.0);
        simulator.SetSamplingTimestepMultiple(6000);
        simulator.SetDt(0.003);

        // Make a dictionary for mosaic cell surface tension parameters
        std::map<boost::shared_ptr<AbstractCellMutationState>, std::array<double, 3> > dictionary_surface_tensions;
        std::array<double, 3> wild_type_surface_tensions = { gamma_a, gamma_b, 1.0 }; // apical, basal, lateral

        dictionary_surface_tensions[p_state_wt] = wild_type_surface_tensions;

        // Add the forces
        MAKE_PTR(GeneralizedVolumeConservingForce<3>, p_force);
        MAKE_PTR(SurfaceTensionSubForce<3>, p_tension_subforce);

        p_force->AddSubForce(p_tension_subforce);
        p_force->SetSimulationInstance(&simulator);
        p_tension_subforce->SetSurfaceTensionParametersByMutation(dictionary_surface_tensions);
        // p_tension_subforce->SetSimulatedAnnealingParameters(0.1, 50.0, 0.0); // temp zero anyways
        p_tension_subforce->SetSimulationInstance(&simulator);
        simulator.AddForce(p_force);

        MAKE_PTR(PopulationSurfaceTensionEnergyWriter<3>, p_energy_writer);
        p_energy_writer->SetTensionForcePointer(p_tension_subforce);
        cell_population.AddMonolayerPopulationWriter(p_energy_writer);

        // Make the map from surface tensions to colors in surface evolver
        std::map<double, int> color_map;
        color_map[1.0] = 2;
        color_map[0.6] = 3;
        color_map[0.7] = 4;
        color_map[0.8] = 5;
        color_map[0.9] = 6;
        color_map[1.1] = 7;
        color_map[1.2] = 8;
        color_map[1.3] = 9;

        double avg_cell_volume = 1.0;

        /*
        MAKE_PTR(SurfaceEvolverSaveModifier<3>, p_surface_evolver_modifier);
        p_surface_evolver_modifier->SetSaveAtEnd();
        p_surface_evolver_modifier->SetMapTensionToColor(color_map);
        p_surface_evolver_modifier->SetSurfaceTensionForce(p_force);
        simulator.AddSimulationModifier(p_surface_evolver_modifier);
        */

        MAKE_PTR_ARGS(GeometricalTargetVolumeModifier<3>, p_growth_modifier, (&cell_population));
        p_growth_modifier->SetGrowthDuration(0.0);
        p_growth_modifier->SetT1AdaptationDuration(0.100);
        p_growth_modifier->SetReferenceTargetVolume(avg_cell_volume);
        simulator.AddSimulationModifier(p_growth_modifier);

        // Add simply supported boundary condition at edges
        c_vector<double, 3> mid_zero_point = zero_vector<double>(3);
        c_vector<double, 3> mid_normal_vector = zero_vector<double>(3);
        mid_normal_vector[2] = 1.0;
        mid_zero_point[2] = height / 2.0;
        MAKE_PTR_ARGS(SimplySupportedEdgeBoundaryCondition<3>,
                      p_supported_bc,
                      (&cell_population, mid_zero_point, mid_normal_vector));
        simulator.AddCellPopulationBoundaryCondition(p_supported_bc);

        simulator.Solve();

        simulator.SetEndTime(400.0);
        simulator.SetSamplingTimestepMultiple(400);
        simulator.SetDt(0.1);

        simulator.Solve();

        // Fix the boundary nodes for shearing
        double s_edge = 1.0 / cbrt(3.0 * 3.0 * (gamma_a + gamma_b) / 2.0);
        double l_x = (num_cells_x)*sqrt(3) * s_edge;
        double l_y = (num_cells_y) * 1.5 * s_edge;

        unsigned num_nodes = p_mesh->GetNumNodes();
        std::vector<Node<3>*> node_vector_bottom;
        std::vector<Node<3>*> node_vector_top;
        for (unsigned index_node = 0; index_node < num_nodes; ++index_node)
        {
            if (p_mesh->GetNode(index_node)->rGetLocation()[1] < s_edge / 2.0 + 1e-3)
            {
                node_vector_bottom.push_back(p_mesh->GetNode(index_node));
            }
            else if (p_mesh->GetNode(index_node)->rGetLocation()[1] > l_y - 1e-3)
            {
                node_vector_top.push_back(p_mesh->GetNode(index_node));
            }
        }
        unsigned num_top_nodes = node_vector_top.size();
        unsigned num_bottom_nodes = node_vector_bottom.size();

        // Write the simulation parameters into file in folder
        SaveConstantsToFile(simulator.GetOutputDirectory(), simulation_parameters, 0, 0,
                            num_top_nodes, num_bottom_nodes);

        MAKE_PTR_ARGS(FixNodeDimensionsBoundaryCondition<3>,
                      p_lower_boundary_bc,
                      (&cell_population, node_vector_bottom, false, true, false));
        simulator.AddCellPopulationBoundaryCondition(p_lower_boundary_bc);

        MAKE_PTR_ARGS(FixNodeDimensionsBoundaryCondition<3>,
                      p_upper_boundary_bc,
                      (&cell_population, node_vector_top, false, true, false));
        simulator.AddCellPopulationBoundaryCondition(p_upper_boundary_bc);

        // We assume a stress tensor (sigma_xx, sigma_xy // sigma_xy, sigma_yy)
        // with sigma_xx = StressVariable, sigma_yy = 0.5*StressVariable, sigma_xy = 2.0*StressVariable
        // Stress is force_per_node * num_nodes / length
        double sigma_yx = simulation_parameters.StressVariable * l_x / num_top_nodes;

        top_force[0] = sigma_yx;
        bottom_force[0] = -sigma_yx;

        MAKE_PTR(RectangularEdgeSubForce<3>, p_force_edge);
        p_force_edge->SetTopEdgeNodes(node_vector_top);
        p_force_edge->SetBottomEdgeNodes(node_vector_bottom);

        p_force_edge->SetTopEdgeForce(top_force);
        p_force_edge->SetBottomEdgeForce(bottom_force);

        // p_force_edge->UseInternalCoordinates();
        // p_force_edge->ParallelizeOppositeEdgeNormals();
        //  simulator.AddForce(p_force_edge);
        p_force->AddSubForce(p_force_edge);

        simulator.SetEndTime(1000.0);
        simulator.SetSamplingTimestepMultiple(400);
        simulator.SetDt(0.1);

        simulator.Solve();

        // Now set up convergence monitoring with the maximum area change
        // Then iterate until the area change over the intervall of 2000
        // is negligible or we reach a maximum time
        std::map<int, double> map_faces_to_areas;
        CalculateMaximalAreaDifference(&cell_population, map_faces_to_areas);
        double max_area_diff = 1e10;
        double max_time = 4000.0;
        while (max_area_diff > AREA_DIFFERENCE_CUTOFF && max_time <= MAX_TIME)
        {
            simulator.SetEndTime(max_time);
            simulator.SetSamplingTimestepMultiple(simulation_parameters.PrintEveryNSteps);
            simulator.SetDt(simulation_parameters.FinalDT);
            simulator.Solve();

            max_time += simulation_parameters.MaxTimeLoopStepping;
            max_area_diff = CalculateMaximalAreaDifference(&cell_population, map_faces_to_areas);
        }
        ExecutableSupport::Print("\nDone with simulation loop.\n");
    }
    catch (const Exception& e)
    {
        ExecutableSupport::PrintError(e.GetMessage());
        exit_code = ExecutableSupport::EXIT_ERROR;
    }

    // Optional - write the machine info to file.
    ExecutableSupport::WriteMachineInfoFile("machine_info");

    // End by finalizing PETSc, and returning a suitable exit code.
    // 0 means 'no error'
    ExecutableSupport::FinalizePetsc();
    return exit_code;
}
