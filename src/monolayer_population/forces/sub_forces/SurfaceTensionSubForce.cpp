#include "SurfaceTensionSubForce.hpp"
#include "LinearSystem.hpp"
#include "MonolayerVertexElement.hpp"
#include "PetscVecTools.hpp"
#include "RandomNumberGenerator.hpp"
#include "SimulationTime.hpp"

template <unsigned DIM>
SurfaceTensionSubForce<DIM>::SurfaceTensionSubForce()
        : AbstractSubForce<DIM>(),
          mSurfaceTensionParameters(), // Initialize as empty vector, then use 1.0 as standard value
          mPerformActiveT1Swaps(false),
          mAnnealingDecayTime(1.0),
          mInitialTemperature(1e-10),
          mApicalStandardTension(0.0),
          mBasalStandardTension(0.0),
          mLateralStandardTension(0.0),
          mT1AreTemperatureDependent(true),
          mT1TransitionRate(100.0),
          mModeOfSurfaceTensionCalculation(0),
          mpSimulation(nullptr)
{
}

template <unsigned DIM>
SurfaceTensionSubForce<DIM>::~SurfaceTensionSubForce()
{
}

template <unsigned DIM>
NodeSubForceMap<DIM> SurfaceTensionSubForce<DIM>::GetForceContributions(
    MonolayerVertexBasedCellPopulation<DIM>& rCellPopulation,
    AbstractCellBasedSimulation<DIM, DIM>* pSimulation)
{
    // Define some helper variables
    MonolayerVertexBasedCellPopulation<DIM>* p_cell_population = static_cast<MonolayerVertexBasedCellPopulation<DIM>*>(&rCellPopulation);
    MutableMonolayerVertexMesh<DIM, DIM>& r_mesh = p_cell_population->rGetMesh();
    unsigned num_nodes = p_cell_population->GetNumNodes();
    unsigned num_elements = p_cell_population->GetNumElements();
    unsigned num_faces = p_cell_population->rGetMesh().GetNumFaces();

    // First we update the surface tensions, which is necessary after cell division or remeshing
    // with different mutation states
    UpdateSurfaceTensions(&(rCellPopulation));

    // If we do not have any Surface Tensions use a standard value
    if (mSurfaceTensionParameters.empty())
    {
        mSurfaceTensionParameters = std::vector<double>(num_faces, 1.0);
    }

    // Begin by computing the areas of each element in the mesh, to avoid having to do this multiple times
    std::map<MonolayerVertexElement<DIM - 1, DIM>*, double> faces_areas_old;
    std::vector<double> cell_types(num_elements);

    // Map of nodes to forces which will be returned
    NodeSubForceMap<DIM> map_forces_on_nodes;

    // NOTE: NOW IN MUTABLE MESH

    // Set up everything for the simulated annealing
    // SimulationTime* p_simulation_time = SimulationTime::Instance();
    // double current_time = p_simulation_time->GetTime();
    // double temperature;
    // if(current_time < mAnnealingDecayTime)
    //{
    //   temperature = mInitialTemperature * (1.0 -current_time/mAnnealingDecayTime);
    // }
    // else
    //{
    //   temperature=1.0e-10;
    // }
    //
    //// Set the probability for Active T1s
    // if(DoPerformActiveT1Swaps())
    //{
    //   double time_step = p_simulation_time->GetTimeStep();
    //   double probability = 0.0;
    //   if(IsT1TransitonRateTemperatureDependent())
    //     probability = temperature*time_step*mT1TransitionRate / r_mesh.GetNumFaces();
    //   else
    //     probability = time_step*mT1TransitionRate / r_mesh.GetNumFaces();
    //   r_mesh.SetActiveT1SwapProbability(probability);
    // }
    // else
    //{
    //   r_mesh.SetActiveT1SwapProbability(0.0);
    // }
    //
    // TODO: negative annealing decay time must be implemented elsewhere!

    // SimulationTime* p_simulation_time = SimulationTime::Instance();
    // double current_time = p_simulation_time->GetTime();
    // double temperature;
    // if (current_time < mAnnealingDecayTime)
    //{
    //     temperature = mInitialTemperature * (1.0 - current_time / mAnnealingDecayTime);
    // }
    // else
    //{
    //     // Negative AnnealingDecayTime allows for growth
    //     if (mAnnealingDecayTime < 0.0)
    //         temperature = mInitialTemperature - current_time / mAnnealingDecayTime;
    //     else
    //         temperature = 1.0e-10;
    // }
    //
    //  Set the probability for Active T1s
    // if (DoPerformActiveT1Swaps())
    //{
    //     double time_step = p_simulation_time->GetTimeStep();
    //     double probability = 0.0;
    //     if (IsT1TransitonRateTemperatureDependent())
    //         probability = temperature * time_step * mT1TransitionRate / r_mesh.GetNumFaces();
    //     else
    //         probability = time_step * mT1TransitionRate / r_mesh.GetNumFaces();
    //     r_mesh.SetActiveT1SwapProbability(probability);
    // }
    // else
    //{
    //     r_mesh.SetActiveT1SwapProbability(0.0);
    // }

    // Iterate over vertices in the cell population
    for (unsigned node_index = 0; node_index < num_nodes; node_index++)
    {
        Node<DIM>* p_this_node = p_cell_population->GetNode(node_index);

        /*
         * The force on this Node is given by the gradient of the total free
         * energy of the CellPopulation, evaluated at the position of the vertex.
         *
         * Note that since the movement of this Node only affects the free energy
         * of the CellPtrs containing it, we can just consider the contributions
         * to the free energy gradient from each of these CellPtrs.
         */
        c_vector<double, DIM> surface_tension_contribution = zero_vector<double>(DIM);
        std::vector<c_vector<double, DIM> > vector_volume_gradients;

        // Find the indices of the elements owned by this node
        std::set<unsigned> containing_elem_indices = p_cell_population->GetNode(node_index)->rGetContainingElementIndices();

        // Iterate over these elements
        for (std::set<unsigned>::iterator iter = containing_elem_indices.begin();
             iter != containing_elem_indices.end();
             ++iter)
        {
            // Get this element, its index and its number of nodes
            MonolayerVertexElement<DIM, DIM>* p_element = p_cell_population->GetElement(*iter);
            unsigned num_faces_elem = p_element->GetNumFaces();

            // For each face in the elements which contain the node/vertex
            // we calculate the term from the area gradients
            for (unsigned index_face = 0; index_face < num_faces_elem; ++index_face)
            {
                MonolayerVertexElement<DIM - 1, DIM>* p_face = p_element->GetFace(index_face);
                unsigned face_index = p_face->GetIndex();
                unsigned local_node_index_in_face = p_face->GetNodeLocalIndex(node_index);
                // If the node is not even in the face, we just continue to the next face
                if (local_node_index_in_face == UINT_MAX)
                {
                    continue;
                }

                double surface_tension_parameter = mSurfaceTensionParameters[face_index];

                // Note that the 1/2 from counting faces twice is part of the surface tension parameter
                surface_tension_contribution -= surface_tension_parameter * r_mesh.GetAreaGradientOfFaceAtNode(p_face, local_node_index_in_face);
            }
        }
        map_forces_on_nodes[p_this_node] = surface_tension_contribution;
    }
    return map_forces_on_nodes;
}

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
NodeSubForceMap<1> SurfaceTensionSubForce<1>::GetForceContributions(
    MonolayerVertexBasedCellPopulation<1>& rCellPopulation,
    AbstractCellBasedSimulation<1, 1>* pSimulation)
{
    EXCEPTION("SurfaceTensionSubForce not possible in 1D.");
}
/// \endcond Get Doxygen to ignore, since it's confused by these templates

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
void SurfaceTensionSubForce<1>::CreateSurfaceTensionParametersForCells(double apical_surface_tension,
                                                                       double basal_surface_tension,
                                                                       double lateral_surface_tension,
                                                                       MonolayerVertexMesh<1, 1>* p_mesh)
{
    EXCEPTION("CreateSurfaceTensionParametersForCells not implemented in one dimension.");
}
/// \endcond Get Doxygen to ignore, since it's confused by these templates

template <unsigned DIM>
void SurfaceTensionSubForce<DIM>::CreateSurfaceTensionParametersForCells(double apical_surface_tension,
                                                                         double basal_surface_tension,
                                                                         double lateral_surface_tension,
                                                                         MonolayerVertexMesh<DIM, DIM>* p_mesh)
{
    // Set the mode to uniform
    mModeOfSurfaceTensionCalculation = 1;

    unsigned num_faces = p_mesh->GetNumFaces();
    if (mSurfaceTensionParameters.size() != num_faces)
    {
        mSurfaceTensionParameters = std::vector<double>(num_faces, 1.0);
    }

    for (unsigned faceIndex = 0; faceIndex < num_faces; ++faceIndex)
    {
        MonolayerVertexElement<DIM - 1, DIM>* p_face = p_mesh->GetFace(faceIndex);
        bool boundary_face = p_face->IsBoundaryFace();
        double boundary_factor = boundary_face ? 2.0 : 1.0;
        if (p_face->GetFaceType() == MonolayerVertexElementType::Apical)
        {
            mSurfaceTensionParameters[faceIndex] = apical_surface_tension;
        }
        else if (p_face->GetFaceType() == MonolayerVertexElementType::Basal)
        {
            mSurfaceTensionParameters[faceIndex] = basal_surface_tension;
        }
        else if (p_face->GetFaceType() == MonolayerVertexElementType::Lateral)
        {
            mSurfaceTensionParameters[faceIndex] = boundary_factor * lateral_surface_tension / 2.0;
        }
        else
        {
            mSurfaceTensionParameters[faceIndex] = boundary_factor * 0.0;
        }
    }
    // Save into internal variables
    mApicalStandardTension = apical_surface_tension;
    mBasalStandardTension = basal_surface_tension;
    mLateralStandardTension = lateral_surface_tension;
}

template <unsigned DIM>
double SurfaceTensionSubForce<DIM>::GetApicalStandardTension()
{
    return mApicalStandardTension;
}

template <unsigned DIM>
double SurfaceTensionSubForce<DIM>::GetBasalStandardTension()
{
    return mBasalStandardTension;
}

template <unsigned DIM>
double SurfaceTensionSubForce<DIM>::GetLateralStandardTension()
{
    return mLateralStandardTension;
}

template <unsigned DIM>
double SurfaceTensionSubForce<DIM>::GetSurfaceTensionParameter(unsigned face)
{
    if (face >= mSurfaceTensionParameters.size())
    {
        EXCEPTION("Error! Index to large in GetSurfaceTensionParameter");
    }
    return mSurfaceTensionParameters[face];
}

template <unsigned DIM>
std::vector<double> SurfaceTensionSubForce<DIM>::GetSurfaceTensionParameters()
{
    return mSurfaceTensionParameters;
}

template <unsigned DIM>
bool SurfaceTensionSubForce<DIM>::DoPerformActiveT1Swaps()
{
    return mPerformActiveT1Swaps;
}

template <unsigned DIM>
void SurfaceTensionSubForce<DIM>::SetSurfaceTensionParameter(unsigned face, double surfaceTensionParameter)
{
    // Set mode to matrix
    mModeOfSurfaceTensionCalculation = 2;

    if (face >= mSurfaceTensionParameters.size())
    {
        EXCEPTION("Error! Index to large in GetSurfaceTensionParameter");
    }
    mSurfaceTensionParameters[face] = surfaceTensionParameter;
}

template <unsigned DIM>
void SurfaceTensionSubForce<DIM>::SetSurfaceTensionParametersByMutation(std::map<boost::shared_ptr<AbstractCellMutationState>, std::array<double, 3> > mappingMutationTension)
{
    // Set mode to dictionary
    mModeOfSurfaceTensionCalculation = 3;

    mMappingMutationTension = mappingMutationTension;
}

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
void SurfaceTensionSubForce<1>::UpdateSurfaceTensions(AbstractCellPopulation<1>* pCellPopulation)
{
    EXCEPTION("SurfaceTensionSubForce not possible in 1D.");
}
/// \endcond Get Doxygen to ignore, since it's confused by these templates

template <unsigned DIM>
void SurfaceTensionSubForce<DIM>::UpdateSurfaceTensions(AbstractCellPopulation<DIM>* pCellPopulation)
{
    if (mModeOfSurfaceTensionCalculation == 1)
    {
        MonolayerVertexBasedCellPopulation<DIM>* p_cell_population = static_cast<MonolayerVertexBasedCellPopulation<DIM>*>(pCellPopulation);

        // Update with the old values.
        CreateSurfaceTensionParametersForCells(mApicalStandardTension, mBasalStandardTension,
                                               mLateralStandardTension, &(p_cell_population->rGetMesh()));
    }
    else if (mModeOfSurfaceTensionCalculation == 3)
    {
        MonolayerVertexBasedCellPopulation<DIM>* p_cell_population = static_cast<MonolayerVertexBasedCellPopulation<DIM>*>(pCellPopulation);

        // Update with the old mapping to tensions
        UpdateSurfaceTensionsByMutation(p_cell_population);
    }
    // In the matrix case we cannot consistently update.
}

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
void SurfaceTensionSubForce<1>::UpdateSurfaceTensionsByMutation(AbstractCellPopulation<1>* pCellPopulation)
{
    EXCEPTION("SurfaceTensionSubForce not possible in 1D.");
}
/// \endcond Get Doxygen to ignore, since it's confused by these templates

template <unsigned DIM>
void SurfaceTensionSubForce<DIM>::UpdateSurfaceTensionsByMutation(AbstractCellPopulation<DIM>* pCellPopulation)
{
    assert(mModeOfSurfaceTensionCalculation == 3);

    MonolayerVertexBasedCellPopulation<DIM>* p_cell_population = static_cast<MonolayerVertexBasedCellPopulation<DIM>*>(pCellPopulation);

    // We initialize an empty surface tension parameter vector
    unsigned num_faces = p_cell_population->rGetMesh().GetNumFaces();
    mSurfaceTensionParameters = std::vector<double>(num_faces, 0.0);

    for (typename MonolayerVertexMesh<DIM, DIM>::MonolayerVertexElementIterator elem_iter = p_cell_population->rGetMesh().GetElementIteratorBegin();
         elem_iter != p_cell_population->rGetMesh().GetElementIteratorEnd();
         ++elem_iter)
    {
        // Get the tensions for the cell with the correpsonding mutation
        CellPtr p_cell = p_cell_population->GetCellUsingLocationIndex(elem_iter->GetIndex());
        boost::shared_ptr<AbstractCellMutationState> p_mutation = p_cell->GetMutationState();
        std::array<double, 3> surface_tensions_cell = mMappingMutationTension[p_mutation];
        double apical_surface_tension = surface_tensions_cell[0];
        double basal_surface_tension = surface_tensions_cell[1];
        double lateral_surface_tension = surface_tensions_cell[2];

        unsigned num_faces_in_element = elem_iter->GetNumFaces();
        for (unsigned face_index = 0; face_index < num_faces_in_element; ++face_index)
        {
            MonolayerVertexElement<DIM - 1, DIM>* p_face = elem_iter->GetFace(face_index);
            bool boundary_face = p_face->IsBoundaryFace();
            unsigned global_face_index = p_face->GetIndex();
            assert(global_face_index < num_faces);

            // If it is a shared face, we take the average of the cell-specific face values
            // In the sum of forces we also go twice over this face. Therefore multiply by 1/4.
            // For boundaries we remultiply by 2, because this face is only shared by ONE element!
            // Therefore we need to save it as 1.0*x and don't average
            double boundary_factor = boundary_face ? 2.0 : 1.0;

            if (p_face->GetFaceType() == MonolayerVertexElementType::Apical)
            {
                mSurfaceTensionParameters[global_face_index] += apical_surface_tension;
            }
            else if (p_face->GetFaceType() == MonolayerVertexElementType::Basal)
            {
                mSurfaceTensionParameters[global_face_index] += basal_surface_tension;
            }
            else if (p_face->GetFaceType() == MonolayerVertexElementType::Lateral)
            {
                // We sum over each face twice and divide by 2.0 to get the average.
                mSurfaceTensionParameters[global_face_index] += boundary_factor * lateral_surface_tension / 4.0;
            }
            else
            {
                mSurfaceTensionParameters[global_face_index] += boundary_factor * 0.0;
            }
        }
    }
}

template <unsigned DIM>
void SurfaceTensionSubForce<DIM>::SetSurfaceTensionParameters(std::vector<double> surfaceTensionParameters)
{
    mSurfaceTensionParameters = surfaceTensionParameters;
}

template <unsigned DIM>
void SurfaceTensionSubForce<DIM>::SetT1TransitionParameters(double T1TransitionRate, bool T1AreTemperatureDependent)
{
    mT1TransitionRate = T1TransitionRate;
    mT1AreTemperatureDependent = T1AreTemperatureDependent;
}

template <unsigned DIM>
double SurfaceTensionSubForce<DIM>::GetT1TransitionRate()
{
    return mT1TransitionRate;
}

template <unsigned DIM>
bool SurfaceTensionSubForce<DIM>::IsT1TransitonRateTemperatureDependent()
{
    return mT1AreTemperatureDependent;
}

template <unsigned DIM>
void SurfaceTensionSubForce<DIM>::SetPerformActiveT1Swaps(bool doActiveT1Swaps)
{
    mPerformActiveT1Swaps = doActiveT1Swaps;
}

template <>
double SurfaceTensionSubForce<1>::CalculateTotalEnergy(AbstractCellPopulation<1>& rCellPopulation)
{
    EXCEPTION("SurfaceTensionSubForce not possible in 1D.");
}

template <unsigned DIM>
double SurfaceTensionSubForce<DIM>::CalculateTotalEnergy(AbstractCellPopulation<DIM>& rCellPopulation)
{
    double total_energy = 0.0;

    // First we update the surface tensions, which is necessary after cell division or remeshing
    // with different mutation states
    UpdateSurfaceTensions(&(rCellPopulation));

    MonolayerVertexBasedCellPopulation<DIM>* p_cell_population = static_cast<MonolayerVertexBasedCellPopulation<DIM>*>(&rCellPopulation);
    unsigned num_faces = p_cell_population->rGetMesh().GetNumFaces();

    // If we do not have any Surface Tensions use a standard value
    if (mSurfaceTensionParameters.empty())
    {
        mSurfaceTensionParameters = std::vector<double>(num_faces, 1.0);
    }

    for (typename MonolayerVertexMesh<DIM, DIM>::MonolayerVertexElementIterator elem_iter = p_cell_population->rGetMesh().GetElementIteratorBegin();
         elem_iter != p_cell_population->rGetMesh().GetElementIteratorEnd();
         ++elem_iter)
    {
        unsigned num_faces_in_element = elem_iter->GetNumFaces();
        for (unsigned face_index = 0; face_index < num_faces_in_element; ++face_index)
        {
            MonolayerVertexElement<DIM - 1, DIM>* p_face = elem_iter->GetFace(face_index);
            unsigned global_face_index = p_face->GetIndex();
            double face_area = p_cell_population->rGetMesh().CalculateAreaOfFace(p_face);
            total_energy += mSurfaceTensionParameters[global_face_index] * face_area;
        }
    }

    return total_energy;
}

template <unsigned DIM>
void SurfaceTensionSubForce<DIM>::SetSimulatedAnnealingParameters(double relativeNoiseStrength,
                                                                  double annealingDecayTime,
                                                                  double initialTemperature)
{
    mRelativeNoiseStrength = relativeNoiseStrength;
    mAnnealingDecayTime = annealingDecayTime;
    mInitialTemperature = initialTemperature;
    EXCEPTION("SimulatedAnnealing has been removed from SurfaceTensionSubForce. Use ActiveT1ProbabilityModifier instead!");
}

template <unsigned DIM>
void SurfaceTensionSubForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    if (mSurfaceTensionParameters.empty())
    {
        *rParamsFile << "\t\t\t<SurfaceTensionParameterZero>" << 1.0 << "</SurfaceTensionParameterZero>\n";
    }
    else
    {
        *rParamsFile << "\t\t\t<SurfaceTensionParameterZero>" << mSurfaceTensionParameters[0] << "</SurfaceTensionParameterZero>\n";
    }
    *rParamsFile << "\t\t\t<ModeOfSurfaceTensionCalculation>" << mModeOfSurfaceTensionCalculation << "</ModeOfSurfaceTensionCalculation>\n";

    *rParamsFile << "\t\t\t<StandardApicalSurfaceTension>" << mApicalStandardTension << "</StandardApicalSurfaceTension>\n";
    *rParamsFile << "\t\t\t<StandardBasalSurfaceTension>" << mBasalStandardTension << "</StandardBasalSurfaceTension>\n";
    *rParamsFile << "\t\t\t<StandardLateralSurfaceTension>" << mLateralStandardTension << "</StandardLateralSurfaceTension>\n";

    *rParamsFile << "\t\t\t<MappingMutationTension>\n";
    for (auto it = mMappingMutationTension.begin(); it != mMappingMutationTension.end(); ++it)
    {
        *rParamsFile << "\t\t\t\t<MutationCellType>" << it->first->GetColour() << "</MutationCellType>\n";
        *rParamsFile << "\t\t\t\t<ApicalSurfaceTension>" << it->second[0] << "</ApicalSurfaceTension>\n";
        *rParamsFile << "\t\t\t\t<BasalSurfaceTension>" << it->second[1] << "</BasalSurfaceTension>\n";
        *rParamsFile << "\t\t\t\t<LateralSurfaceTension>" << it->second[2] << "</LateralSurfaceTension>\n";
    }
    *rParamsFile << "\t\t\t<\nMappingMutationTension>\n";

    *rParamsFile << "\t\t\t<T1AreTemperatureDependent>" << mT1AreTemperatureDependent << "</T1AreTemperatureDependent>\n";
    *rParamsFile << "\t\t\t<T1TransitionRate>" << mT1TransitionRate << "</T1TransitionRate>\n";

    *rParamsFile << "\t\t\t<SimulatedAnnealingNoiseStrength>" << mRelativeNoiseStrength << "</SimulatedAnnealingNoiseStrength>\n";
    *rParamsFile << "\t\t\t<SimulatedAnnealingDecayTime>" << mAnnealingDecayTime << "</SimulatedAnnealingDecayTime>\n";
    *rParamsFile << "\t\t\t<SimulatedAnnealingInitialTemperature>" << mInitialTemperature << "</SimulatedAnnealingInitialTemperature>\n";

    *rParamsFile << "\t\t\t<ActiveT1Swaps>" << mPerformActiveT1Swaps << "</ActiveT1Swaps>\n";

    // Call method on direct parent class
    AbstractForce<DIM>::OutputForceParameters(rParamsFile);
}

// Explicit instantiation
template class SurfaceTensionSubForce<1>;
template class SurfaceTensionSubForce<2>;
template class SurfaceTensionSubForce<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"