/**
 * Perform a two-sided compression of a shell centered around the origin. The initial mesh is a MutableMonolayerVertexMesh read from a given
 * directory and which was written with MonolayerVertexMeshWriter. The compression is done symmetrically from both sides.
 */

#include <algorithm>
#include <cmath>
#include <ctime> // for current timestamp
#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>

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
#include "CellsGenerator.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "FaceTypeWriter.hpp"
#include "FiniteThicknessSimulation3d.hpp"
#include "FiniteThicknessVertexMeshFromReaderGenerator.hpp"
#include "MechanicalMosaicCellMutationState.hpp"
#include "MonolayerVertexBasedCellPopulation.hpp"
#include "NoCellCycleModel.hpp"
#include "OffLatticeSimulation.hpp"
#include "PopulationSurfaceTensionEnergyWriter.hpp"
#include "PopulationT1CounterWriter.hpp"
#include "SmallAreaDampingModifier.hpp"
#include "SmartPointers.hpp"

#include "BoundaryHeightWriter.hpp"
#include "PopulationBoundaryForceWriter.hpp"
#include "PopulationEdgeLengthWriter.hpp"
#include "PopulationLumenVolumeWriter.hpp"

#include "ActiveT1ProbabilityModifier.hpp"
#include "DecayingRandomForce.hpp"
#include "GeneralizedVolumeConservingForce.hpp"
#include "GeometricalTargetVolumeModifier.hpp"
#include "LumenPressureSubForce.hpp"
#include "MovingBoundarySubForce.hpp"
#include "RandomSubForce.hpp"
#include "SimpleTargetVolumeModifier.hpp"
#include "SurfaceEvolverSaveModifier.hpp"
#include "SurfaceTensionSubForce.hpp"

#include "AbstractMovingBoundary.hpp"
#include "OneStepProfile.hpp"
#include "PlateMovingBoundary.hpp"
#include "SphereMovingBoundary.hpp"
#include "Temperature.hpp"

struct paramsMonolayer
{
    double GammaA = 1.0;
    double GammaB = 1.0;
    double MaxSimulationTime = 10000.0;
    double AreaDifferenceCutoff = 1e-5;
    unsigned PrintEveryNStepsCompression = 100;
    unsigned PrintEveryNStepsLoopRelaxation = 200;
    double DTCompression = 0.01;
    double DTRelaxation = 0.01;
    double DTLoopRelaxation = 0.1;
    double TimeForInitialRelaxation = 100;
    double T1Threshold = 0.3;
    unsigned FixLumen = 0;
    double IndentationPercent = 0.7;
    double BoundaryVelocity = 0.4;
    double DampingConstantRelax = 1.0;
    double DampingConstantStep = 1.0;
    std::string InputMeshPath;
    double RandomForceStrength = 0.0;
    double ActiveT1ProbRelaxation = 0.0;
    unsigned ConstantActiveT1Prob = 0;
    double CellRearrangementRatio = 1.5;
    unsigned UseHardPotential = 1;
    double SphereIndenterRadius = 0.0;
    double LumenPressure = 0.0;
    unsigned Seed = 0;
    std::string FileComment = "";
};

// Function to parse command line arguments
paramsMonolayer ParseArguments(int argc, char* argv[])
{
    paramsMonolayer params;

    for (int i = 1; i < argc; i++)
    {
        std::string arg = argv[i];

        if (arg == "-GammaA" && i + 1 < argc)
        {
            params.GammaA = std::stod(argv[i + 1]);
        }
        else if (arg == "-GammaB" && i + 1 < argc)
        {
            params.GammaB = std::stod(argv[i + 1]);
        }
        else if (arg == "-MaxSimTime" && i + 1 < argc)
        {
            params.MaxSimulationTime = std::stod(argv[i + 1]);
        }
        else if (arg == "-AreaDiffCutoff" && i + 1 < argc)
        {
            params.AreaDifferenceCutoff = std::stod(argv[i + 1]);
        }
        else if (arg == "-PrintEveryNStepsC" && i + 1 < argc)
        {
            params.PrintEveryNStepsCompression = std::stoi(argv[i + 1]);
        }
        else if (arg == "-PrintEveryNStepsR" && i + 1 < argc)
        {
            params.PrintEveryNStepsLoopRelaxation = std::stoi(argv[i + 1]);
        }
        else if (arg == "-DTC" && i + 1 < argc)
        {
            params.DTCompression = std::stod(argv[i + 1]);
        }
        else if (arg == "-DTR" && i + 1 < argc)
        {
            params.DTRelaxation = std::stod(argv[i + 1]);
        }
        else if (arg == "-DTLR" && i + 1 < argc)
        {
            params.DTLoopRelaxation = std::stod(argv[i + 1]);
        }
        else if (arg == "-T1Threshold" && i + 1 < argc)
        {
            params.T1Threshold = std::stod(argv[i + 1]);
        }
        else if (arg == "-Lumen" && i + 1 < argc)
        {
            params.FixLumen = std::stod(argv[i + 1]);
        }
        else if (arg == "-Indent" && i + 1 < argc)
        {
            params.IndentationPercent = std::stod(argv[i + 1]);
        }
        else if (arg == "-BoundaryVel" && i + 1 < argc)
        {
            params.BoundaryVelocity = std::stod(argv[i + 1]);
        }
        else if (arg == "-DCRelax" && i + 1 < argc)
        {
            params.DampingConstantRelax = std::stod(argv[i + 1]);
        }
        else if (arg == "-DCStep" && i + 1 < argc)
        {
            params.DampingConstantStep = std::stod(argv[i + 1]);
        }
        else if (arg == "-InputMeshPath" && i + 1 < argc)
        {
            params.InputMeshPath = argv[i + 1];
        }
        else if (arg == "-RandomStrength" && i + 1 < argc)
        {
            params.RandomForceStrength = std::stod(argv[i + 1]);
        }
        else if (arg == "-ActiveT1Prob" && i + 1 < argc)
        {
            params.ActiveT1ProbRelaxation = std::stod(argv[i + 1]);
        }
        else if (arg == "-ConstActiveT1Prob" && i + 1 < argc)
        {
            params.ConstantActiveT1Prob = std::stod(argv[i + 1]);
        }
        else if (arg == "-CellRearrangementRatio" && i + 1 < argc)
        {
            params.CellRearrangementRatio = std::stod(argv[i + 1]);
        }
        else if (arg == "-HardPotential" && i + 1 < argc)
        {
            params.UseHardPotential = std::stod(argv[i + 1]);
        }
        else if (arg == "-TimeForInitialRelaxation" && i + 1 < argc)
        {
            params.TimeForInitialRelaxation = std::stod(argv[i + 1]);
        }
        else if (arg == "-SphereIndenterRadius" && i + 1 < argc)
        {
            params.SphereIndenterRadius = std::stod(argv[i + 1]);
        }
        else if (arg == "-LumenPressure" && i + 1 < argc)
        {
            params.LumenPressure = std::stod(argv[i + 1]);
        }
        else if (arg == "-Seed" && i + 1 < argc)
        {
            params.Seed = std::stod(argv[i + 1]);
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
std::string CreateFilePath(const paramsMonolayer& params, unsigned numCells)
{
    // Get the current date and time
    auto now = std::chrono::system_clock::now();
    std::time_t now_time = std::chrono::system_clock::to_time_t(now);

    // Convert the current date and time to the desired format
    std::tm* local_time = std::localtime(&now_time);
    std::ostringstream date_time;
    date_time << std::put_time(local_time, "%Y%m%d%H%M%S");

    std::ostringstream filepath;
    filepath << date_time.str() << "_"
             << "Compress_Sphere_"
             << "GA" << params.GammaA << "_"
             << "GB" << params.GammaB << "_"
             << "N" << numCells << "_"
             << "ADC" << params.AreaDifferenceCutoff << "_"
             << "IP" << params.IndentationPercent << "_"
             << "HP" << params.UseHardPotential << "_"
             << params.FileComment;
    return filepath.str();
}

void SaveConstantsToFile(const std::string& filepath, const paramsMonolayer& params, double numCells)
{
    OutputFileHandler handler(filepath, false);
    out_stream output_file = handler.OpenOutputFile("CompressSphereParameters.txt");

    if (output_file->is_open())
    {
        *output_file << "GammaA = " << params.GammaA << std::endl;
        *output_file << "GammaB = " << params.GammaB << std::endl;
        *output_file << "MaxSimulationTime = " << params.MaxSimulationTime << std::endl;
        *output_file << "AreaDifferenceCutoff = " << params.AreaDifferenceCutoff << std::endl;
        *output_file << "PrintEveryNStepsCompression = " << params.PrintEveryNStepsCompression << std::endl;
        *output_file << "PrintEveryNStepsLoopRelaxation = " << params.PrintEveryNStepsLoopRelaxation << std::endl;
        *output_file << "DTCompression = " << params.DTCompression << std::endl;
        *output_file << "DTRelaxation = " << params.DTRelaxation << std::endl;
        *output_file << "DTLoopRelaxation = " << params.DTLoopRelaxation << std::endl;
        *output_file << "T1Threshold = " << params.T1Threshold << std::endl;
        *output_file << "FixLumen = " << params.FixLumen << std::endl;
        *output_file << "IndentationPercent = " << params.IndentationPercent << std::endl;
        *output_file << "BoundaryVelocity = " << params.BoundaryVelocity << std::endl;
        *output_file << "DampingConstantRelax = " << params.DampingConstantRelax << std::endl;
        *output_file << "DampingConstantStep = " << params.DampingConstantStep << std::endl;
        *output_file << "InputMeshPath = " << params.InputMeshPath << std::endl;
        *output_file << "ActiveT1Prob = " << params.ActiveT1ProbRelaxation << std::endl;
        *output_file << "ConstantActiveT1Prob = " << params.ConstantActiveT1Prob << std::endl;
        *output_file << "CellRearrangementRatio = " << params.CellRearrangementRatio << std::endl;
        *output_file << "HardPotential = " << params.UseHardPotential << std::endl;
        *output_file << "TimeForInitialRelaxation = " << params.TimeForInitialRelaxation << std::endl;
        *output_file << "SphereIndenterRadius = " << params.SphereIndenterRadius << std::endl;
        *output_file << "LumenPressure = " << params.LumenPressure << std::endl;
        *output_file << "Seed = " << params.Seed << std::endl;
        output_file->close();
        ExecutableSupport::Print("\nConstants saved to file successfully.\n");
    }
    else
    {
        EXCEPTION("Error opening output file. Maybe filepath wrong?");
    }
}

std::unordered_map<std::string, double> parseParametersFromFile(const std::string& filename)
{
    std::unordered_map<std::string, double> params;
    std::ifstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "Error opening file: " << filename << std::endl;
        return params;
    }

    std::string line;
    while (std::getline(file, line))
    {
        std::istringstream iss(line);
        std::string key;
        double value;

        if (iss >> key >> value)
        {
            params[key] = value;
            std::cout << "Added: Key: " << key << ", Value: " << value << std::endl;
        }
        else
        {
            std::cout << "Parsed: Key: '" << key << "', Value: " << value << std::endl;

            std::cout << "Parsing error for line: " << line << std::endl;
        }
    }

    file.close();
    return params;
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

        paramsMonolayer parameters = ParseArguments(argc, argv);

        std::cout << "Input String for Mesh : " << parameters.InputMeshPath << std::endl;

        if (parameters.Seed > 0)
            RandomNumberGenerator::Instance()->Reseed(parameters.Seed);

        // Read mesh and configure it
        MonolayerVertexMeshReader<3, 3> mesh_reader(parameters.InputMeshPath + "final_mesh/evolved_random_sphere");
        std::cout << "Read Mesh : " << parameters.InputMeshPath << std::endl;
        FiniteThicknessVertexMeshFromReaderGenerator generator(mesh_reader, 0.0, 0.001);
        MutableMonolayerVertexMesh<3, 3>* p_mesh = generator.GetMesh();
        // We use protorosettes as intermediate states
        p_mesh->SetProtorosetteFormationProbability(1.0);
        p_mesh->SetProtorosetteResolutionProbabilityPerTimestep(0.1);
        p_mesh->SetCellRearrangementThreshold(parameters.T1Threshold);
        // p_mesh -> SetCheckForInternalIntersections(true);
        std::string outputPath = CreateFilePath(parameters, p_mesh->GetNumElements());
        SaveConstantsToFile(outputPath, parameters, p_mesh->GetNumElements());

        // Generate cells
        std::vector<CellPtr> cells;
        CellsGenerator<NoCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumElements());

        // Create population
        MonolayerVertexBasedCellPopulation<3> cell_population(*p_mesh, cells);
        cell_population.SetOutputCellRearrangementLocations(false);
        cell_population.SetRestrictVertexMovementBoolean(false);

        // Add writers
        cell_population.AddCellWriter<CellThicknessWriter>();
        cell_population.AddCellWriter<CellSymmetryWriter>();
        cell_population.AddCellWriter<CellCentroidWriter>();
        cell_population.AddMonolayerPopulationWriter<PopulationLumenVolumeWriter>();
        cell_population.AddMonolayerPopulationWriter<BoundaryHeightWriter>();
        cell_population.AddMonolayerPopulationWriter<PopulationBoundaryForceWriter>();
        cell_population.AddMonolayerPopulationWriter<PopulationT1CounterWriter>();
        cell_population.AddMonolayerPopulationWriter<PopulationEdgeLengthWriter>();

        // Create simulator
        FiniteThicknessSimulation3d simulator(cell_population);
        simulator.SetOutputDirectory(outputPath);

        // Create Forces
        MAKE_PTR(GeneralizedVolumeConservingForce<3>, p_force);
        MAKE_PTR(SurfaceTensionSubForce<3>, p_tension_subforce);
        MAKE_PTR(MovingBoundarySubForce<3>, p_moving_boundary_subforce);
        p_force->AddSubForce(p_tension_subforce);
        p_force->AddSubForce(p_moving_boundary_subforce);
        p_force->SetSimulationInstance(&simulator);
        p_tension_subforce->SetSimulationInstance(&simulator);
        p_tension_subforce->CreateSurfaceTensionParametersForCells(parameters.GammaA, parameters.GammaB, 1.0, p_mesh); // apical, basal, lateral

        if (parameters.LumenPressure > 0.0)
        {
            MAKE_PTR_ARGS(LumenPressureSubForce<3>, p_lumen_pressure_subforce, (parameters.LumenPressure));
            p_force->AddSubForce(p_lumen_pressure_subforce);
        }

        MAKE_PTR(PopulationSurfaceTensionEnergyWriter<3>, p_energy_writer);
        p_energy_writer->SetTensionForcePointer(p_tension_subforce);
        cell_population.AddMonolayerPopulationWriter(p_energy_writer);

        simulator.AddForce(p_force);

        // Add modifier
        MAKE_PTR_ARGS(GeometricalTargetVolumeModifier<3>, p_growth_modifier, (&cell_population));
        p_growth_modifier->SetGrowthDuration(0.0);
        p_growth_modifier->SetT1AdaptationDuration(0.100);
        double avg_cell_volume = 1.0;
        p_growth_modifier->SetReferenceTargetVolume(avg_cell_volume);
        simulator.AddSimulationModifier(p_growth_modifier);

        // Add Boundaries
        double boundary_force_constant = 0.05;

        if (!parameters.UseHardPotential)
        {
            boundary_force_constant = cell_population.GetDampingConstantNormal() / parameters.DTCompression;
        }

        // find positions
        double max_z_position = cell_population.GetMaxZPosition();
        double min_z_position = cell_population.GetMinZPosition();
        double starting_position = (abs(max_z_position) >= abs(min_z_position)) ? abs(max_z_position) : abs(min_z_position);
        // start a bit away from vertices
        starting_position = starting_position * 1.1;

        double distance_per_step = (starting_position * (1 - parameters.IndentationPercent));
        double time_for_one_step = distance_per_step / parameters.BoundaryVelocity;
        double time_after_step = time_for_one_step;

        MAKE_PTR_ARGS(OneStepProfile<3>, linear_move_profile_top, (0.0, time_after_step, starting_position, (starting_position - distance_per_step)));
        MAKE_PTR_ARGS(OneStepProfile<3>, linear_move_profile_bottom, (0.0, time_after_step, -1.0 * starting_position, ((-1.0 * starting_position) + distance_per_step)));

        if (parameters.SphereIndenterRadius > 0.0)
        {
            MAKE_PTR_ARGS(SphereMovingBoundary<3>, boundary_top, (boundary_force_constant, linear_move_profile_top, true, parameters.SphereIndenterRadius));
            MAKE_PTR_ARGS(SphereMovingBoundary<3>, boundary_bottom, (boundary_force_constant, linear_move_profile_bottom, false, parameters.SphereIndenterRadius));
            cell_population.AddMovingBoundary(boundary_top);
            cell_population.AddMovingBoundary(boundary_bottom);
        }
        else
        {
            MAKE_PTR_ARGS(PlateMovingBoundary<3>, boundary_top, (boundary_force_constant, linear_move_profile_top, true, parameters.UseHardPotential));
            MAKE_PTR_ARGS(PlateMovingBoundary<3>, boundary_bottom, (boundary_force_constant, linear_move_profile_bottom, false, parameters.UseHardPotential));
            cell_population.AddMovingBoundary(boundary_top);
            cell_population.AddMovingBoundary(boundary_bottom);
        }

        if (parameters.FixLumen)
        {
            double target_vol = p_mesh->GetLumenVolume();
            std::cout << "Target Lumen Volume: " << target_vol << std::endl;
            cell_population.SetTargetLumenVolume(target_vol);
            p_force->SetConstantLumenVolumeFlag(true);
        }

        cell_population.SetDampingConstantNormal(parameters.DampingConstantStep);

        simulator.SetEndTime(time_after_step);
        simulator.SetDt(parameters.DTCompression);
        simulator.SetOutputDirectory(outputPath);
        simulator.SetSamplingTimestepMultiple(parameters.PrintEveryNStepsCompression);

        std::cout << "Start step" << std::endl;
        simulator.Solve();

        // keep the boundary at the current position
        linear_move_profile_top->Stop();
        linear_move_profile_bottom->Stop();

        // Add constant temperature
        MAKE_PTR(Temperature, temperature);
        temperature->SetTemperatureToConstant();

        // Add constant random noise if wanted
        MAKE_PTR_ARGS(RandomSubForce<3>, p_random_subforce, (parameters.RandomForceStrength, temperature));
        p_force->AddSubForce(p_random_subforce);

        // Add active T1 modifier
        MAKE_PTR(ActiveT1ProbabilityModifier<3>, p_activeT1_modifier);
        p_activeT1_modifier->SetTemperature(temperature);

        if (parameters.ConstantActiveT1Prob)
        {
            p_activeT1_modifier->SetActiveT1Rate(parameters.ActiveT1ProbRelaxation);
        }
        else
        {
            p_activeT1_modifier->SetActiveT1BoltzmannParameter(parameters.ActiveT1ProbRelaxation);
        }
        simulator.AddSimulationModifier(p_activeT1_modifier);

        // Set passive T1 threshold
        p_mesh->SetCellRearrangementRatio(parameters.CellRearrangementRatio);

        cell_population.SetDampingConstantNormal(parameters.DampingConstantRelax);

        simulator.SetOutputDirectory(outputPath);
        simulator.SetEndTime(time_after_step + parameters.TimeForInitialRelaxation);
        simulator.SetDt(parameters.DTRelaxation);
        simulator.SetSamplingTimestepMultiple(parameters.PrintEveryNStepsLoopRelaxation);
        std::cout << "Start Relaxation" << std::endl;
        simulator.Solve();

        // Remove random force
        p_force->RemoveSubForce(p_random_subforce);

        simulator.SetEndTime(8000);
        simulator.SetSamplingTimestepMultiple(parameters.PrintEveryNStepsLoopRelaxation);
        simulator.SetDt(parameters.DTLoopRelaxation);
        simulator.SetOutputDirectory(outputPath);
        simulator.Solve();

        unsigned loop_counter_final_relaxation = 0;
        p_activeT1_modifier->SetActiveT1BoltzmannParameter(0.0);
        p_activeT1_modifier->SetActiveT1Rate(0.0);

        std::map<int, double> map_faces_to_areas;
        CalculateMaximalAreaDifference(&cell_population, map_faces_to_areas);
        double max_area_diff = 1e10;
        double max_time = 9000.0;
        while (max_area_diff > parameters.AreaDifferenceCutoff && max_time <= parameters.MaxSimulationTime)
        {
            std::cout << "Start final relaxation loop " << loop_counter_final_relaxation++ << std::endl;
            simulator.SetEndTime(max_time);
            simulator.SetSamplingTimestepMultiple(parameters.PrintEveryNStepsLoopRelaxation);
            simulator.SetDt(parameters.DTLoopRelaxation);
            simulator.SetOutputDirectory(outputPath);
            simulator.Solve();

            max_time += 1000.0;
            max_area_diff = CalculateMaximalAreaDifference(&cell_population, map_faces_to_areas);
        }
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
