#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>

#include "FileFinder.hpp"

#include "SimulationTime.hpp"

#include "Exception.hpp"
#include "ExecutableSupport.hpp"
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
#include "FiniteThicknessRandomizedSphereMeshGenerator.hpp"
#include "FiniteThicknessSimulation3d.hpp"
#include "MechanicalMosaicCellMutationState.hpp"
#include "MonolayerVertexBasedCellPopulation.hpp"
#include "NoCellCycleModel.hpp"
#include "OffLatticeSimulation.hpp"
#include "SmallAreaDampingModifier.hpp"
#include "SmartPointers.hpp"

#include "DecayingRandomForce.hpp"
#include "GeneralizedVolumeConservingForce.hpp"
#include "GeometricalTargetVolumeModifier.hpp"
#include "SimpleTargetVolumeModifier.hpp"
#include "SurfaceEvolverSaveModifier.hpp"
#include "SurfaceTensionSubForce.hpp"

struct paramsMosaicSphere
{
    unsigned Seed = 0;
    unsigned NumCells = 100;
    unsigned PrintEveryNSteps = 1000;
    double GammaA = 1.2;
    double GammaB = 0.95;
    double GammaL = 1.0;
    double MaxSimulationTime = 1e4;
    double MaxTimeLoopStepping = 10000.0;
    double FinalDT = 0.1;
    double AreaDifferenceCutoff = 1e-5;
    double SmallAreaVolumeDampingCoefficient = 0.1;
    double SmallAreaDampingThreshold = 0.04;
    double SmallAreaForceDampingCoefficient = 0.1;
    std::string FileComment = "";
};

paramsMosaicSphere ParseArguments(int argc, char* argv[])
{
    paramsMosaicSphere params;

    for (int i = 1; i < argc; i++)
    {
        std::string arg = argv[i];

        if (arg == "-h")
        {
            std::cout << "Usage: program [options]" << std::endl;
            std::cout << "Options:" << std::endl;
            std::cout << "  -NumCells <value>  : How many cells to use." << std::endl;
            std::cout << "  -Seed <value>  : Seed to use in random number generator (std:0)" << std::endl;
            std::cout << "  -GammaA <value>  : Apical tension (std: 1.2)" << std::endl;
            std::cout << "  -GammaB <value>  : Basal tension (std: 0.95)." << std::endl;
            std::cout << "  -GammaL <value>  : Lateral tension (std: 1.0)." << std::endl;
            std::cout << "  -SmallAreaForceDampingCoefficient <value>   : How much forces are damped for small faces (std: 0.1)" << std::endl;
            std::cout << "  -SmallAreaVolumeDampingCoefficient <value>   : How much volume corrections are damped for cells with small faces (std: 0.1)" << std::endl;
            std::cout << "  -SmallAreaDampingThreshold <value>   : Below which face area forces are damped (std: 0.04)" << std::endl;
            std::cout << "  -MaxSimTime <value>   : MaxSimulationTime  (std: 1e4)" << std::endl;
            std::cout << "  -AreaDiffCutoff <value>   : Area difference cutoff to terminate (std: 1e-5)" << std::endl;
            std::cout << "  -PrintEveryNSteps <value>   : How often to save to file (std: 1000)" << std::endl;
            std::cout << "  -LoopStepping <value>   : How large the timestep will be increased if break condition not fulfilled (std: 1e4)" << std::endl;
            std::cout << "  -FinalDT <value>   : dt in final simulation loop (std: 0.1)" << std::endl;
            std::cout << "  -Comment <value>   : Comment to append to filename" << std::endl;
            EXCEPTION("Aborting as help flag was given.");
        }
        else if (arg == "-NumCells" && i + 1 < argc)
        {
            params.NumCells = std::stod(argv[i + 1]);
        }
        else if (arg == "-Seed" && i + 1 < argc)
        {
            params.Seed = std::stoi(argv[i + 1]);
        }
        else if (arg == "-GammaA" && i + 1 < argc)
        {
            params.GammaA = std::stod(argv[i + 1]);
        }
        else if (arg == "-GammaB" && i + 1 < argc)
        {
            params.GammaB = std::stod(argv[i + 1]);
        }
        else if (arg == "-GammaL" && i + 1 < argc)
        {
            params.GammaL = std::stod(argv[i + 1]);
        }
        else if (arg == "-SmallAreaForceDampingCoefficient" && i + 1 < argc)
        {
            params.SmallAreaForceDampingCoefficient = std::stod(argv[i + 1]);
        }
        else if (arg == "-SmallAreaVolumeDampingCoefficient" && i + 1 < argc)
        {
            params.SmallAreaVolumeDampingCoefficient = std::stod(argv[i + 1]);
        }
        else if (arg == "-SmallAreaDampingThreshold" && i + 1 < argc)
        {
            params.SmallAreaDampingThreshold = std::stod(argv[i + 1]);
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
std::string CreateFilePath(const paramsMosaicSphere& params)
{
    // Get the current date and time
    auto now = std::chrono::system_clock::now();
    std::time_t now_time = std::chrono::system_clock::to_time_t(now);

    // Convert the current date and time to the desired format
    std::tm* local_time = std::localtime(&now_time);
    std::ostringstream date_time;
    date_time << std::put_time(local_time, "%d%m%Y_%H_%M");

    std::ostringstream filepath;
    filepath << "RandomizedSphere_"
             << "NUM" << params.NumCells << "_"
             << "GA" << params.GammaA << "_"
             << "GB" << params.GammaB << "_"
             << "GL" << params.GammaL << "_"
             << "VDMP" << params.SmallAreaVolumeDampingCoefficient << "_"
             << "FDMP" << params.SmallAreaForceDampingCoefficient << "_"
             << "DMPTHR" << params.SmallAreaDampingThreshold << "_"
             << "DT" << params.FinalDT << "_"
             << "ADC" << params.AreaDifferenceCutoff << "_"
             << "SEED" << params.Seed
             << params.FileComment << "_"
             << date_time.str();
    return filepath.str();
}

void SaveConstantsToFile(const std::string& filepath, const paramsMosaicSphere& params)
{
    OutputFileHandler handler(filepath, false);
    out_stream output_file = handler.OpenOutputFile("SphereParameters.txt");

    if (output_file->is_open())
    {
        *output_file << "NumCells = " << params.NumCells << std::endl;
        *output_file << "Seed = " << params.Seed << std::endl;
        *output_file << "GammaA = " << params.GammaA << std::endl;
        *output_file << "GammaB = " << params.GammaB << std::endl;
        *output_file << "GammaL = " << params.GammaL << std::endl;
        *output_file << "SmallAreaDampingThreshold = " << params.SmallAreaDampingThreshold << std::endl;
        *output_file << "SmallAreaForceDampingCoefficient = " << params.SmallAreaForceDampingCoefficient << std::endl;
        *output_file << "SmallAreaVolumeDampingCoefficient = " << params.SmallAreaVolumeDampingCoefficient << std::endl;
        *output_file << "MaxSimulationTime = " << params.MaxSimulationTime << std::endl;
        *output_file << "AreaDifferenceCutoff = " << params.AreaDifferenceCutoff << std::endl;
        *output_file << "PrintEveryNSteps = " << params.PrintEveryNSteps << std::endl;
        *output_file << "MaxTimeLoopStepping = " << params.MaxTimeLoopStepping << std::endl;
        *output_file << "FinalDT = " << params.FinalDT << std::endl;

        output_file->close();
        ExecutableSupport::Print("\nConstants saved to file successfully.\n");
    }
    else
    {
        EXCEPTION("Error opening output file. Maybe filepath wrong?");
    }
}

double calculateMaximalAreaDifference(MonolayerVertexBasedCellPopulation<3>* pCellPopulation, std::map<int, double>& rMapFaceAreas)
{
    // If map is empty, we initialize the map and return a negative value
    double max_area_diff = -1.0;
    bool initialize = rMapFaceAreas.empty();

    MonolayerVertexMesh<3, 3>& r_mesh = pCellPopulation->rGetMesh();
    // Loop over faces, because if we look at cell's faces, we go over the same faces multiple times
    unsigned num_faces = r_mesh.GetNumFaces();
    for (unsigned index_face = 0; index_face < num_faces; ++index_face)
    {
        MonolayerVertexElement<3 - 1, 3>* p_face = r_mesh.GetFace(index_face);
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

        // Parse the input
        paramsMosaicSphere simulation_parameters = ParseArguments(argc, argv);

        if (simulation_parameters.Seed > 0)
        {
            RandomNumberGenerator::Instance()->Reseed(simulation_parameters.Seed);
        }

        double MAX_TIME = simulation_parameters.MaxSimulationTime;
        double AREA_DIFFERENCE_CUTOFF = simulation_parameters.AreaDifferenceCutoff;

        // Set up the shell according to the parameters

        unsigned num_cells = simulation_parameters.NumCells;

        // alpha = apical/lateral surface tension
        double alpha = simulation_parameters.GammaA / simulation_parameters.GammaL;
        // beta = basal/lateral surface tension
        double beta = simulation_parameters.GammaB / simulation_parameters.GammaL;

        double height = cbrt(2.0 / sqrt(3.0)) * cbrt((alpha + beta) * (alpha + beta));
        double radius = sqrt(num_cells / M_PI / height) / 2.0;
        double hex_mid_edge = cbrt(2.0 / 3.0 / 3.0 / (alpha + beta));
        std::cout << "\nHeight, radius:" << height << " , " << radius << std::flush;

        // length at which T1s become favorable: l^*=0.66/3^(2/3)/gamma^(1/3)
        double t1_length = hex_mid_edge * 0.66 / 4.0;
        FiniteThicknessRandomizedSphereMeshGenerator generator(num_cells, t1_length, 0.001, height, radius - height / 2.0, true);

        MutableMonolayerVertexMesh<3, 3>* p_mesh = generator.GetMesh();
        // We use protorosettes as intermediate states
        p_mesh->SetProtorosetteFormationProbability(1.0);
        p_mesh->SetProtorosetteResolutionProbabilityPerTimestep(0.1);

        // We generate the cell objects
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state_wt);
        MAKE_PTR(MechanicalMosaicCellMutationState, p_state_mosaic_1);
        MAKE_PTR(MechanicalMosaicCellMutationState, p_state_mosaic_2);
        MAKE_PTR(MechanicalMosaicCellMutationState, p_state_mosaic_3);
        MAKE_PTR(MechanicalMosaicCellMutationState, p_state_mosaic_4);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);

        for (unsigned cell_index = 0; cell_index < num_cells; cell_index++)
        {
            NoCellCycleModel* p_model = new NoCellCycleModel;
            p_model->SetDimension(3);

            CellPtr p_cell = nullptr;
            // unsigned rand_int = RandomNumberGenerator::Instance()->randMod(5);

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
        cell_population.AddCellWriter<CellSymmetryWriter>();
        cell_population.AddCellWriter<CellCentroidWriter>();
        cell_population.AddCellWriter<CellApicalAreaWriter>();
        cell_population.AddCellWriter<CellBasalAreaWriter>();
        cell_population.AddCellWriter<CellOpeningAngleWriter>();
        cell_population.AddFaceWriter<FaceTypeWriter>();

        FiniteThicknessSimulation3d simulator(cell_population);
        simulator.SetOutputDirectory(CreateFilePath(simulation_parameters));

        // Write the simulation parameters into file in folder
        SaveConstantsToFile(simulator.GetOutputDirectory(), simulation_parameters);

        // Make a dictionary for mosaic cell surface tension parameters
        std::map<boost::shared_ptr<AbstractCellMutationState>, std::array<double, 3> > dictionary_surface_tensions;

        std::array<double, 3> wild_type_surface_tensions = { simulation_parameters.GammaA,
                                                             simulation_parameters.GammaB, simulation_parameters.GammaL }; // apical, basal, lateral

        dictionary_surface_tensions[p_state_wt] = wild_type_surface_tensions;
        // Add the forces
        MAKE_PTR(GeneralizedVolumeConservingForce<3>, p_force);
        MAKE_PTR(SurfaceTensionSubForce<3>, p_tension_subforce);

        p_force->AddSubForce(p_tension_subforce);
        p_force->SetSimulationInstance(&simulator);
        p_tension_subforce->SetSurfaceTensionParametersByMutation(dictionary_surface_tensions);
        // p_tension_subforce->SetSimulatedAnnealingParameters(0.1, 50.0, 0.0);
        p_tension_subforce->SetSimulationInstance(&simulator);
        simulator.AddForce(p_force);

        MAKE_PTR(PopulationSurfaceTensionEnergyWriter<3>, p_energy_writer);
        p_energy_writer->SetTensionForcePointer(p_tension_subforce);
        cell_population.AddMonolayerPopulationWriter(p_energy_writer);

        // Make the map from surface tensions to colors in surface evolver
        std::map<double, int> color_map;
        color_map[simulation_parameters.GammaA] = 1;
        color_map[simulation_parameters.GammaB] = 2;
        color_map[simulation_parameters.GammaL] = 3;

        MAKE_PTR(SurfaceEvolverSaveModifier<3>, p_surface_evolver_modifier);
        p_surface_evolver_modifier->SetSaveAtEnd();
        p_surface_evolver_modifier->SetMapTensionToColor(color_map);
        p_surface_evolver_modifier->SetSurfaceTensionSubForce(p_tension_subforce);
        p_surface_evolver_modifier->SetUseRandomizedVolumes(false, 0.0);
        p_surface_evolver_modifier->SetWriteFaceTypeIntoFile(true);
        p_surface_evolver_modifier->SetConsiderLumenAsCell(true);
        simulator.AddSimulationModifier(p_surface_evolver_modifier);

        // Introduce damping to allow for area->zero for small faces
        MAKE_PTR_ARGS(SmallAreaDampingModifier<3>, p_small_area_modifier, (&cell_population, &(*p_force)));
        p_small_area_modifier->SetSmallAreaThreshold(simulation_parameters.SmallAreaDampingThreshold);
        p_small_area_modifier->SetVolumeDampingRatio(simulation_parameters.SmallAreaVolumeDampingCoefficient);
        simulator.AddSimulationModifier(p_small_area_modifier);

        // Target volume modifier
        MAKE_PTR_ARGS(GeometricalTargetVolumeModifier<3>, p_growth_modifier, (&cell_population));
        p_growth_modifier->SetGrowthDuration(0.0);
        p_growth_modifier->SetT1AdaptationDuration(0.100);
        p_growth_modifier->SetReferenceTargetVolume(1.0);
        simulator.AddSimulationModifier(p_growth_modifier);

        // Initial short-time simulation solve
        simulator.SetEndTime(1.20);
        simulator.SetSamplingTimestepMultiple(199);
        // simulator.SetSamplingTimestepMultiple();
        simulator.SetDt(0.003);
        p_small_area_modifier->SetForceDampingRatio(simulation_parameters.SmallAreaForceDampingCoefficient * 0.003);
        simulator.Solve();

        // Initial long-time simulation solve
        cell_population.SetDoInitialVolumeRelaxation(false);
        simulator.SetEndTime(300.0);
        simulator.SetSamplingTimestepMultiple(999);
        p_growth_modifier->SetT1AdaptationDuration(0.400);
        simulator.SetDt(0.1);
        p_small_area_modifier->SetForceDampingRatio(simulation_parameters.SmallAreaForceDampingCoefficient * 0.1);

        // p_tension_subforce->SetSimulatedAnnealingParameters(0.1, 200.0, 0.0);
        // p_tension_subforce->SetT1TransitionParameters(1.0, true);
        p_tension_subforce->SetSimulationInstance(&simulator);
        // p_tension_subforce->SetPerformActiveT1Swaps();
        simulator.Solve();

        // Now set up convergence monitoring with the maximum area change
        // Then iterate until the area change over the intervall of 1000
        // is negligible or we reach a maximum time

        p_mesh->SetCellRearrangementThreshold(0.66 * hex_mid_edge / 2.0);
        // p_tension_subforce->SetSimulatedAnnealingParameters(0.0, 200.0, 0.0);
        // p_tension_subforce->SetPerformActiveT1Swaps(false);

        std::map<int, double> map_faces_to_areas;
        calculateMaximalAreaDifference(&cell_population, map_faces_to_areas);
        double max_area_diff = 1e10;
        double max_time = 1000.0;
        while (max_area_diff > AREA_DIFFERENCE_CUTOFF && max_time <= MAX_TIME)
        {
            simulator.SetEndTime(max_time);
            simulator.SetSamplingTimestepMultiple(simulation_parameters.PrintEveryNSteps);
            simulator.SetDt(simulation_parameters.FinalDT);
            p_small_area_modifier->SetForceDampingRatio(simulation_parameters.SmallAreaForceDampingCoefficient * simulation_parameters.FinalDT);
            simulator.Solve();

            max_time += simulation_parameters.MaxTimeLoopStepping;
            max_area_diff = calculateMaximalAreaDifference(&cell_population, map_faces_to_areas);
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
