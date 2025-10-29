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
#include "FaceIndexWriter.hpp"
#include "FaceTensionWriter.hpp"
#include "FaceTypeWriter.hpp"
#include "PopulationSurfaceTensionEnergyWriter.hpp"

#include "CellProliferativeTypesWriter.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "StemCellProliferativeType.hpp"

#include "CellsGenerator.hpp"
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
#include "GeometricalGrowthTargetVolumeModifier.hpp"
#include "SimpleTargetVolumeModifier.hpp"
#include "SurfaceEvolverSaveModifier.hpp"
#include "SurfaceTensionSubForce.hpp"

#include "ActiveT1ProbabilityModifier.hpp"
#include "Temperature.hpp"

struct paramsGrowingSphere
{
    unsigned NumCells = 12;
    unsigned NumTransitGenerations = 1;
    unsigned PrintEveryNSteps = 1000;
    double MaxSimulationTime = 300.0;
    double G1Duration = 25.0;
    double G2Duration = 25.0;
    double FinalDT = 0.006;
    double T1Rate = 2.0;
    bool T1RateGrows = false;
    double T1GrowthRate = 0.1;
    std::string FileComment = "";
};

paramsGrowingSphere ParseArguments(int argc, char* argv[])
{
    paramsGrowingSphere params;

    for (int i = 1; i < argc; i++)
    {
        std::string arg = argv[i];

        if (arg == "-h")
        {
            std::cout << "Usage: program [options]" << std::endl;
            std::cout << "Options:" << std::endl;
            std::cout << "  -NumCells <value>  : How many cells to use." << std::endl;
            std::cout << "  -NumTransitGenerations <value>  : How many dividing transit generations." << std::endl;
            std::cout << "  -MaxSimTime <value>   : MaxSimulationTime  (std: 300.0)" << std::endl;
            std::cout << "  -PrintEveryNSteps <value>   : How often to save to file (std: 1000)" << std::endl;
            std::cout << "  -G1Duration <value>   : Duration of G1 phase (growth) (std: 25.0)" << std::endl;
            std::cout << "  -G2Duration <value>   : Duration of G2 phase (post-division) (std: 25.0)" << std::endl;
            std::cout << "  -FinalDT <value>   : dt in final simulation loop (std: 0.006)" << std::endl;
            std::cout << "  -T1Rate <value>   : T1 transition rate (std: 2.0)" << std::endl;
            std::cout << "  -GrowingT1Rate <value>   : Sets that T1Rate grows linearly with growth rate (std: 0.1)" << std::endl;
            std::cout << "  -Comment <value>   : Comment to append to filename" << std::endl;
            EXCEPTION("Aborting as help flag was given.");
        }
        else if (arg == "-NumCells" && i + 1 < argc)
        {
            params.NumCells = std::stoi(argv[i + 1]);
        }
        else if (arg == "-NumTransitGenerations" && i + 1 < argc)
        {
            params.NumTransitGenerations = std::stoi(argv[i + 1]);
        }
        else if (arg == "-MaxSimTime" && i + 1 < argc)
        {
            params.MaxSimulationTime = std::stod(argv[i + 1]);
        }
        else if (arg == "-PrintEveryNSteps" && i + 1 < argc)
        {
            params.PrintEveryNSteps = std::stoi(argv[i + 1]);
        }
        else if (arg == "-G1Duration" && i + 1 < argc)
        {
            params.G1Duration = std::stod(argv[i + 1]);
        }
        else if (arg == "-G2Duration" && i + 1 < argc)
        {
            params.G2Duration = std::stod(argv[i + 1]);
        }
        else if (arg == "-T1Rate" && i + 1 < argc)
        {
            params.T1Rate = std::stod(argv[i + 1]);
        }
        else if (arg == "-FinalDT" && i + 1 < argc)
        {
            params.FinalDT = std::stod(argv[i + 1]);
        }
        else if (arg == "-Comment" && i + 1 < argc)
        {
            params.FileComment = argv[i + 1];
        }
        else if (arg == "-GrowingT1Rate" && i + 1 < argc)
        {
            params.T1RateGrows = true;
            params.T1GrowthRate = std::stod(argv[i + 1]);
        }
        i++;
    }

    return params;
}

/*
 * Create the filepath where to save the files
 */
std::string CreateFilePath(const paramsGrowingSphere& params)
{
    // Get the current date and time
    auto now = std::chrono::system_clock::now();
    std::time_t now_time = std::chrono::system_clock::to_time_t(now);

    // Convert the current date and time to the desired format
    std::tm* local_time = std::localtime(&now_time);
    std::ostringstream date_time;
    date_time << std::put_time(local_time, "%d%m%Y_%H_%M");

    std::ostringstream filepath;
    filepath << "LinearGrowingT1Sphere_"
             << "NUM" << params.NumCells << "_"
             << "TGEN" << params.NumTransitGenerations << "_"
             << "G1D" << params.G1Duration << "_"
             << "G2D" << params.G2Duration << "_"
             << "DT" << params.FinalDT << "_"
             << "T1R" << params.T1Rate
             << params.FileComment << "_"
             << date_time.str();
    return filepath.str();
}

void SaveConstantsToFile(const std::string& filepath, const paramsGrowingSphere& params)
{
    OutputFileHandler handler(filepath, false);
    out_stream output_file = handler.OpenOutputFile("GrowingSphereParameters.txt");

    if (output_file->is_open())
    {
        *output_file << "NumCells = " << params.NumCells << std::endl;
        *output_file << "NumTransitGenerations = " << params.NumTransitGenerations << std::endl;
        *output_file << "G1Duration = " << params.G1Duration << std::endl;
        *output_file << "G2Duration = " << params.G2Duration << std::endl;
        *output_file << "T1Rate = " << params.T1Rate << std::endl;
        *output_file << "MaxSimulationTime = " << params.MaxSimulationTime << std::endl;
        *output_file << "PrintEveryNSteps = " << params.PrintEveryNSteps << std::endl;
        *output_file << "T1RateGrows = " << params.T1RateGrows << std::endl;
        *output_file << "T1GrowthRate = " << params.T1GrowthRate << std::endl;
        *output_file << "FinalDT = " << params.FinalDT << std::endl;

        output_file->close();
        ExecutableSupport::Print("\nConstants saved to file successfully.\n");
    }
    else
    {
        EXCEPTION("Error opening output file. Maybe filepath wrong?");
    }
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
        paramsGrowingSphere simulation_parameters = ParseArguments(argc, argv);

        // Set up the shell according to the parameters

        unsigned num_cells = simulation_parameters.NumCells;

        double gamma_apical = 0.9;
        double gamma_basal = 0.9;
        double gamma_lateral = 1.0;

        // alpha = apical/lateral surface tension
        // beta = basal/lateral surface tension
        double alpha_avg = gamma_apical / gamma_lateral;
        double beta_avg = gamma_basal / gamma_lateral;

        double height = cbrt(2.0 / sqrt(3.0)) * cbrt((alpha_avg + beta_avg) * (alpha_avg + beta_avg));
        double radius = sqrt(num_cells / M_PI / height) / 2.0;
        std::cout << "\nHeight, radius:" << height << " , " << radius << std::flush;

        // length at which T1s become favorable: l^*=0.66/3^(2/3)/gamma^(1/3)
        double hex_mid_edge = cbrt(2.0 / 3.0 / 3.0 / (alpha_avg + beta_avg));
        double t1_length = hex_mid_edge * 0.66 / 4.0;
        FiniteThicknessRandomizedSphereMeshGenerator generator(num_cells, t1_length, 0.001, height, radius - height / 2.0, true);

        MutableMonolayerVertexMesh<3, 3>* p_mesh = generator.GetMesh();
        // We use protorosettes as intermediate states
        p_mesh->SetProtorosetteFormationProbability(1.0);
        p_mesh->SetProtorosetteResolutionProbabilityPerTimestep(0.1);

        // We generate the cell objects
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        MAKE_PTR(StemCellProliferativeType, p_stem_type);

        for (unsigned cell_index = 0; cell_index < num_cells; cell_index++)
        {
            FixedG1GenerationalCellCycleModel* p_model = new FixedG1GenerationalCellCycleModel;
            p_model->SetDimension(3);
            CellPtr p_cell(new Cell(p_state, p_model));
            p_cell->SetCellProliferativeType(p_stem_type);
            p_model->SetStemCellG1Duration(simulation_parameters.G1Duration);
            p_model->SetTransitCellG1Duration(simulation_parameters.G1Duration);
            p_model->SetMaxTransitGenerations(simulation_parameters.NumTransitGenerations);
            p_model->SetSDuration(0.001);
            p_model->SetG2Duration(simulation_parameters.G2Duration);
            p_model->SetMDuration(0.001);
            double birth_time = -RandomNumberGenerator::Instance()->ranf() * 5.0;
            p_cell->SetBirthTime(birth_time);
            cells.push_back(p_cell);
        }

        // Make a dictionary for mosaic cell surface tension parameters
        std::map<boost::shared_ptr<AbstractCellMutationState>, std::array<double, 3> > dictionary_surface_tensions;
        std::array<double, 3> wild_type_surface_tensions = { gamma_apical, gamma_basal, gamma_lateral }; // apical, basal, lateral
        dictionary_surface_tensions[p_state] = wild_type_surface_tensions;

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
        cell_population.AddCellWriter<CellProliferativeTypesWriter>();
        cell_population.AddFaceWriter<FaceTypeWriter>();
        cell_population.AddFaceWriter<FaceIndexWriter>();

        FiniteThicknessSimulation3d simulator(cell_population);
        simulator.SetOutputDirectory(CreateFilePath(simulation_parameters));

        // Write the simulation parameters into file in folder
        SaveConstantsToFile(simulator.GetOutputDirectory(), simulation_parameters);

        // Add the forces
        MAKE_PTR(GeneralizedVolumeConservingForce<3>, p_force);
        MAKE_PTR(SurfaceTensionSubForce<3>, p_tension_subforce);

        p_force->AddSubForce(p_tension_subforce);
        p_force->SetSimulationInstance(&simulator);
        // p_tension_subforce->SetSurfaceTensionParametersByMutation(dictionary_surface_tensions);
        //  p_tension_subforce->SetSimulatedAnnealingParameters(0.1, 50.0, 0.0);
        p_tension_subforce->SetSurfaceTensionParametersByMutation(dictionary_surface_tensions);
        p_tension_subforce->SetPerformActiveT1Swaps(false);
        // p_tension_subforce->SetSimulatedAnnealingParameters(0.003, 19000000.0, 1.0);
        // p_tension_subforce->SetT1TransitionParameters(simulation_parameters.T1Rate, simulation_parameters.T1RateGrows);

        p_tension_subforce->SetSimulationInstance(&simulator);
        simulator.AddForce(p_force);

        // Add additional writers
        // Cannot use macro to create smart pointer
        boost::shared_ptr<FaceTensionWriter<3, 3> > p_tension_writer(new FaceTensionWriter<3, 3>(&(*p_tension_subforce)));
        cell_population.AddFaceWriter(p_tension_writer);

        MAKE_PTR(PopulationSurfaceTensionEnergyWriter<3>, p_energy_writer);
        p_energy_writer->SetTensionForcePointer(p_tension_subforce);
        cell_population.AddMonolayerPopulationWriter(p_energy_writer);

        // Growth target volume modifier
        MAKE_PTR_ARGS(GeometricalGrowthTargetVolumeModifier<3>, p_growth_modifier, (&cell_population));
        p_growth_modifier->SetT1AdaptationDuration(0.100);
        p_growth_modifier->SetReferenceTargetVolume(1.0);
        // p_growth_modifier->SetT1AdaptationDuration(0.0);
        simulator.AddSimulationModifier(p_growth_modifier);

        // Initial short-time simulation solve
        double initial_time_no_t1s = 2.0 * (simulation_parameters.G1Duration + simulation_parameters.G2Duration);
        simulator.SetEndTime(initial_time_no_t1s);
        simulator.SetSamplingTimestepMultiple(400);
        // simulator.SetSamplingTimestepMultiple(5);
        simulator.SetDt(0.003);
        simulator.Solve();

        // Initial long-time simulation solve
        p_tension_subforce->SetPerformActiveT1Swaps(true);
        // p_tension_subforce->SetT1TransitionParameters(simulation_parameters.T1Rate, simulation_parameters.T1RateGrows);
        // p_mesh->SetActiveT1SwapRate(simulation_parameters.T1Rate);

        // Add active T1 modifier
        MAKE_PTR(ActiveT1ProbabilityModifier<3>, p_activeT1_modifier);
        p_activeT1_modifier->SetActiveT1Rate(simulation_parameters.T1Rate);
        simulator.AddSimulationModifier(p_activeT1_modifier);

        // use annealing for growing t1 rate
        if (simulation_parameters.T1RateGrows)
        {
            // p_tension_subforce->SetSimulatedAnnealingParameters(
            //     0.003,
            //     (1.0 / simulation_parameters.T1GrowthRate) * simulation_parameters.T1Rate,
            //     (simulation_parameters.T1Rate + simulation_parameters.T1GrowthRate * initial_time_no_t1s) / simulation_parameters.T1Rate);

            // Use annealing for growing T1 rate with temperature object
            // System Temperature
            MAKE_PTR(Temperature, temperature);
            temperature->SetStartTime(0.0);
            temperature->SetAbsoluteDecayTime((1.0 / simulation_parameters.T1GrowthRate));

            // Make T-dependent T1 rate (modifier)
            p_activeT1_modifier->SetTemperature(temperature);
        }

        cell_population.SetDoInitialVolumeRelaxation(false);
        simulator.SetEndTime(simulation_parameters.MaxSimulationTime);
        simulator.SetSamplingTimestepMultiple(simulation_parameters.PrintEveryNSteps);
        simulator.SetDt(simulation_parameters.FinalDT);

        p_tension_subforce->SetSimulationInstance(&simulator);
        p_tension_subforce->SetPerformActiveT1Swaps();
        simulator.Solve();

        ExecutableSupport::Print("\nDone with simulation.\n");
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
