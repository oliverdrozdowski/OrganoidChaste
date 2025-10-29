/**
 * Generate a random spherical shell and let it evolve over time. The final configuration of the mesh is written to a file with
 * MonolayerVertexMeshWriter.
 */

#include <algorithm>
#include <cmath>
#include <ctime> // for current timestamp
#include <fstream>
#include <iostream>
#include <string>

#include "AbstractSubForce.hpp"
#include "ActiveT1ProbabilityModifier.hpp"
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
#include "DecayingRandomForce.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "Exception.hpp"
#include "ExecutableSupport.hpp"
#include "FaceTypeWriter.hpp"
#include "FileFinder.hpp"
#include "FiniteThicknessRandomizedSphereMeshGenerator.hpp"
#include "FiniteThicknessSimulation3d.hpp"
#include "GeneralizedVolumeConservingForce.hpp"
#include "GeometricalTargetVolumeModifier.hpp"
#include "LumenPressureSubForce.hpp"
#include "MechanicalMosaicCellMutationState.hpp"
#include "MonolayerVertexBasedCellPopulation.hpp"
#include "MonolayerVertexMeshWriter.hpp"
#include "NoCellCycleModel.hpp"
#include "OffLatticeSimulation.hpp"
#include "OneStepProfile.hpp"
#include "PetscException.hpp"
#include "PetscTools.hpp"
#include "PopulationEdgeLengthWriter.hpp"
#include "PopulationLumenVolumeWriter.hpp"
#include "PopulationSurfaceTensionEnergyWriter.hpp"
#include "PopulationT1CounterWriter.hpp"
#include "RandomSubForce.hpp"
#include "SimpleTargetVolumeModifier.hpp"
#include "SimulationTime.hpp"
#include "SmallAreaDampingModifier.hpp"
#include "SmartPointers.hpp"
#include "SurfaceEvolverSaveModifier.hpp"
#include "SurfaceTensionSubForce.hpp"
#include "Temperature.hpp"

struct paramsMonolayer
{
    unsigned Seed = 0;
    double GammaA = 1.0;
    double GammaB = 1.0;
    int NumCells = 500;
    double MaxSimulationTime = 10000.0;
    double AreaDifferenceCutoff = 1e-5;
    unsigned PrintEveryNSteps = 100;
    unsigned PrintEveryNStepsFinal = 100;
    double DT1 = 0.003;
    double DT2 = 0.1;
    double DT3 = 0.1;
    double TimeForInitialRelaxation = 10.0;
    double TimeForRelaxationWithRandomForceConstant = 200.0;
    double TimeForRelaxationWithRandomForceDecay = 200.0;
    double RandomForceStrength = 100.0;
    double T1Threshold = 0.3;
    double ActiveT1Prob = 0.5;
    unsigned ConstantActiveT1Prob = 0;
    double LumenPressure = 0.0;
    double ConstantLumenFactor = 0;
    unsigned DoRSA = 1;
    std::string FileComment = "";
};

// Function to parse command line arguments
paramsMonolayer ParseArguments(int argc, char* argv[])
{
    paramsMonolayer params;

    for (int i = 1; i < argc; i++)
    {
        std::string arg = argv[i];

        if (arg == "-seed" && i + 1 < argc)
        {
            params.Seed = std::stod(argv[i + 1]);
        }
        else if (arg == "-GammaA" && i + 1 < argc)
        {
            params.GammaA = std::stod(argv[i + 1]);
        }
        else if (arg == "-GammaB" && i + 1 < argc)
        {
            params.GammaB = std::stod(argv[i + 1]);
        }
        else if (arg == "-NumCells" && i + 1 < argc)
        {
            params.NumCells = std::stod(argv[i + 1]);
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
        else if (arg == "-PrintEveryNStepsFinal" && i + 1 < argc)
        {
            params.PrintEveryNStepsFinal = std::stoi(argv[i + 1]);
        }
        else if (arg == "-DT1" && i + 1 < argc)
        {
            params.DT1 = std::stod(argv[i + 1]);
        }
        else if (arg == "-DT2" && i + 1 < argc)
        {
            params.DT2 = std::stod(argv[i + 1]);
        }
        else if (arg == "-DT3" && i + 1 < argc)
        {
            params.DT3 = std::stod(argv[i + 1]);
        }
        else if (arg == "-Time1" && i + 1 < argc)
        {
            params.TimeForInitialRelaxation = std::stod(argv[i + 1]);
        }
        else if (arg == "-Time2" && i + 1 < argc)
        {
            params.TimeForRelaxationWithRandomForceConstant = std::stod(argv[i + 1]);
        }
        else if (arg == "-Time3" && i + 1 < argc)
        {
            params.TimeForRelaxationWithRandomForceDecay = std::stod(argv[i + 1]);
        }
        else if (arg == "-RandomStrength" && i + 1 < argc)
        {
            params.RandomForceStrength = std::stod(argv[i + 1]);
        }
        else if (arg == "-T1Threshold" && i + 1 < argc)
        {
            params.T1Threshold = std::stod(argv[i + 1]);
        }
        else if (arg == "-ActiveT1Prob" && i + 1 < argc)
        {
            params.ActiveT1Prob = std::stod(argv[i + 1]);
        }
        else if (arg == "-ConstActiveT1Prob" && i + 1 < argc)
        {
            params.ConstantActiveT1Prob = std::stod(argv[i + 1]);
        }
        else if (arg == "-LumenPressure" && i + 1 < argc)
        {
            params.LumenPressure = std::stod(argv[i + 1]);
        }
        else if (arg == "-ConstantLumenFactor" && i + 1 < argc)
        {
            params.ConstantLumenFactor = std::stod(argv[i + 1]);
        }
        else if (arg == "-DoRSA" && i + 1 < argc)
        {
            params.DoRSA = std::stod(argv[i + 1]);
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
std::string CreateFilePath(const paramsMonolayer& params)
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
             << "Evolve_Sphere_"
             << "GA" << params.GammaA << "_"
             << "GB" << params.GammaB << "_"
             << "N" << params.NumCells << "_"
             << "T1" << params.T1Threshold << "_"
             << params.FileComment;
    return filepath.str();
}

void SaveConstantsToFile(const std::string& filepath, const paramsMonolayer& params)
{
    OutputFileHandler handler(filepath, false);
    out_stream output_file = handler.OpenOutputFile("parameters.txt");

    if (output_file->is_open())
    {
        *output_file << "Seed = " << params.Seed << std::endl;
        *output_file << "GammaA = " << params.GammaA << std::endl;
        *output_file << "GammaB = " << params.GammaB << std::endl;
        *output_file << "NumCells = " << params.NumCells << std::endl;
        *output_file << "MaxSimulationTime = " << params.MaxSimulationTime << std::endl;
        *output_file << "AreaDifferenceCutoff = " << params.AreaDifferenceCutoff << std::endl;
        *output_file << "PrintEveryNSteps = " << params.PrintEveryNSteps << std::endl;
        *output_file << "PrintEveryNStepsFinal = " << params.PrintEveryNStepsFinal << std::endl;
        *output_file << "DT1 = " << params.DT1 << std::endl;
        *output_file << "DT2 = " << params.DT2 << std::endl;
        *output_file << "DT3 = " << params.DT3 << std::endl;
        *output_file << "Time1 = " << params.TimeForInitialRelaxation << std::endl;
        *output_file << "Time2 = " << params.TimeForRelaxationWithRandomForceConstant << std::endl;
        *output_file << "Time3 = " << params.TimeForRelaxationWithRandomForceDecay << std::endl;
        *output_file << "RandomStrength = " << params.RandomForceStrength << std::endl;
        *output_file << "T1Threshold = " << params.T1Threshold << std::endl;
        *output_file << "ActiveT1Prob = " << params.ActiveT1Prob << std::endl;
        *output_file << "ConstActiveT1Prob = " << params.ConstantActiveT1Prob << std::endl;
        *output_file << "LumenPressure = " << params.LumenPressure << std::endl;
        *output_file << "ConstantLumenFactor = " << params.ConstantLumenFactor << std::endl;
        *output_file << "DoRSA = " << params.DoRSA << std::endl;

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

        paramsMonolayer parameters = ParseArguments(argc, argv);

        std::string outputPath = CreateFilePath(parameters);

        SaveConstantsToFile(outputPath, parameters);

        double height = cbrt(4.0 * 4.0) * cbrt(3.0 * sqrt(3.0) / 2.0) / cbrt(6.0 * 6.0) * cbrt((parameters.GammaA + parameters.GammaB) * (parameters.GammaA + parameters.GammaB));
        double avg_cell_volume = 1.0;

        // for this h and N. make a bit larger because the random distribution of
        // points would fail otherwise
        double inner_radius = sqrt(parameters.NumCells / 4.0 / M_PI / height) - height / 2.0;

        //  length at which T1s become favorable: l^*=0.66/3^(2/3)/gamma^(1/3)
        double t1_length = parameters.T1Threshold * 0.5;
        // double t1_length = 0.0;

        std::cout << "T1 Length: " << t1_length << std::endl;

        // Set the seed of the random number generator
        if (parameters.Seed > 0)
            RandomNumberGenerator::Instance()->Reseed(parameters.Seed);

        // Generate and configure Mesh
        FiniteThicknessRandomizedSphereMeshGenerator generator(
            parameters.NumCells, t1_length, 0.001, height, inner_radius, parameters.DoRSA);
        MutableMonolayerVertexMesh<3, 3>* p_mesh = generator.GetMesh();
        // We use protorosettes as intermediate states
        p_mesh->SetProtorosetteFormationProbability(1.0);
        p_mesh->SetProtorosetteResolutionProbabilityPerTimestep(0.1);

        // Generate Cells
        std::vector<CellPtr> cells;
        CellsGenerator<NoCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumElements());

        // Create Population
        MonolayerVertexBasedCellPopulation<3> cell_population(*p_mesh, cells);
        cell_population.SetOutputCellRearrangementLocations(false);
        cell_population.SetRestrictVertexMovementBoolean(false);
        cell_population.SetDoInitialVolumeRelaxation(true);

        // Add writers
        cell_population.AddCellWriter<CellThicknessWriter>();
        cell_population.AddCellWriter<CellSymmetryWriter>();
        cell_population.AddCellWriter<CellCentroidWriter>();
        cell_population.AddCellWriter<CellApicalAreaWriter>();

        cell_population.AddMonolayerPopulationWriter<PopulationLumenVolumeWriter>();
        cell_population.AddMonolayerPopulationWriter<PopulationT1CounterWriter>();
        cell_population.AddMonolayerPopulationWriter<PopulationEdgeLengthWriter>();

        // Create Simulator
        FiniteThicknessSimulation3d simulator(cell_population);
        simulator.SetOutputDirectory(outputPath);
        simulator.SetSamplingTimestepMultiple(parameters.PrintEveryNSteps);
        simulator.SetDt(parameters.DT1);

        // Create Forces
        MAKE_PTR(GeneralizedVolumeConservingForce<3>, p_force);
        MAKE_PTR(SurfaceTensionSubForce<3>, p_tension_subforce);
        p_force->AddSubForce(p_tension_subforce);
        p_force->SetSimulationInstance(&simulator);
        p_tension_subforce->SetSimulationInstance(&simulator);
        p_tension_subforce->CreateSurfaceTensionParametersForCells(parameters.GammaA, parameters.GammaB, 1.0, p_mesh); // apical, basal, lateral
        simulator.AddForce(p_force);

        MAKE_PTR(PopulationSurfaceTensionEnergyWriter<3>, p_energy_writer);
        p_energy_writer->SetTensionForcePointer(p_tension_subforce);
        cell_population.AddMonolayerPopulationWriter(p_energy_writer);

        // Add Modifier
        MAKE_PTR_ARGS(GeometricalTargetVolumeModifier<3>, p_growth_modifier, (&cell_population));
        p_growth_modifier->SetGrowthDuration(0.0);
        p_growth_modifier->SetT1AdaptationDuration(0.100);
        p_growth_modifier->SetReferenceTargetVolume(avg_cell_volume);
        simulator.AddSimulationModifier(p_growth_modifier);

        double target_vol = 4.0 / 3.0 * M_PI * pow(inner_radius, 3);
        std::cout << "Target Lumen Volume: " << target_vol << std::endl;
        cell_population.SetTargetLumenVolume(target_vol * 1.2);
        p_force->SetConstantLumenVolumeFlag(true);

        simulator.SetEndTime(parameters.TimeForInitialRelaxation);
        std::cout << "Start first relaxation" << std::endl;
        simulator.Solve();

        p_force->SetConstantLumenVolumeFlag(false);
        if (parameters.ConstantLumenFactor > 0)
        {
            cell_population.SetTargetLumenVolume(target_vol * parameters.ConstantLumenFactor);
            p_force->SetConstantLumenVolumeFlag(true);
        }
        else if (parameters.LumenPressure > 0)
        {
            MAKE_PTR_ARGS(LumenPressureSubForce<3>, p_lumen_pressure_subforce, (parameters.LumenPressure));
            p_force->AddSubForce(p_lumen_pressure_subforce);
        }

        // System Temperature
        MAKE_PTR(Temperature, temperature);

        // Evolve with random force
        MAKE_PTR_ARGS(RandomSubForce<3>, p_random_subforce, (parameters.RandomForceStrength, temperature));
        p_force->AddSubForce(p_random_subforce);

        // Add active T1 modifier
        MAKE_PTR(ActiveT1ProbabilityModifier<3>, p_activeT1_modifier);
        p_activeT1_modifier->SetTemperature(temperature);

        if (parameters.ConstantActiveT1Prob)
        {
            p_activeT1_modifier->SetActiveT1Rate(parameters.ConstantActiveT1Prob / parameters.DT1);
        }
        else
        {
            p_activeT1_modifier->SetActiveT1BoltzmannParameter(parameters.ActiveT1Prob);
        }
        simulator.AddSimulationModifier(p_activeT1_modifier);

        // Set passive T1 threshold
        p_mesh->SetCellRearrangementThreshold(parameters.T1Threshold);

        cell_population.SetDoInitialVolumeRelaxation(false);

        simulator.SetEndTime(parameters.TimeForInitialRelaxation + parameters.TimeForRelaxationWithRandomForceConstant);
        simulator.SetOutputDirectory(outputPath);
        simulator.SetSamplingTimestepMultiple(parameters.PrintEveryNSteps);
        simulator.SetDt(parameters.DT2);

        std::cout << "Start constant relaxation with random force" << std::endl;
        simulator.Solve();

        // Let the temperature decay
        temperature->SetStartTimeToNow();
        temperature->SetRelativeDecayTime(parameters.TimeForRelaxationWithRandomForceDecay);

        simulator.SetEndTime(parameters.TimeForInitialRelaxation + parameters.TimeForRelaxationWithRandomForceConstant + parameters.TimeForRelaxationWithRandomForceDecay);
        simulator.SetOutputDirectory(outputPath);
        simulator.SetSamplingTimestepMultiple(parameters.PrintEveryNSteps);
        simulator.SetDt(parameters.DT2);

        std::cout << "Start relaxation with random force decay" << std::endl;
        simulator.Solve();

        std::cout << "Remove random force" << std::endl;
        p_force->RemoveSubForce(p_random_subforce);
        p_activeT1_modifier->SetActiveT1Rate(0.0);
        p_activeT1_modifier->SetActiveT1BoltzmannParameter(0.0);
        p_mesh->SetCellRearrangementThreshold(0.001);
        p_mesh->SetCellRearrangementRatio(parameters.T1Threshold / 0.001);

        // p_force->SetConstantLumenVolumeFlag(false);

        // Now set up convergence monitoring with the maximum area change
        // Then iterate until the area change over the interval of 2000
        // is negligible or we reach a maximum time
        std::map<int, double> map_faces_to_areas;
        CalculateMaximalAreaDifference(&cell_population, map_faces_to_areas);
        double max_area_diff = 1e10;
        double max_time = 2000.0;

        unsigned loop_counter = 0;

        while (max_area_diff > parameters.AreaDifferenceCutoff && max_time <= parameters.MaxSimulationTime)
        {
            std::cout << "Start loop number " << loop_counter++ << std::endl;
            simulator.SetEndTime(max_time);
            simulator.SetSamplingTimestepMultiple(parameters.PrintEveryNStepsFinal);
            simulator.SetDt(parameters.DT3);
            simulator.SetOutputDirectory(outputPath);
            simulator.Solve();

            max_time += 1000.0;
            max_area_diff = CalculateMaximalAreaDifference(&cell_population, map_faces_to_areas);
        }

        try
        {
            // Save Mesh to File
            std::ostringstream meshPath;
            meshPath << outputPath << "/final_mesh/";

            MonolayerVertexMeshWriter<3, 3> mesh_writer(meshPath.str(), "evolved_random_sphere", false);

            mesh_writer.WriteFilesUsingMesh(*p_mesh);
        }
        catch (const Exception& e)
        {
            ExecutableSupport::PrintError(e.GetMessage());
            exit_code = ExecutableSupport::EXIT_ERROR;
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
