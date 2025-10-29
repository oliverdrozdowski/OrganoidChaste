/*

Copyright (c) 2005-2019, University of Oxford.
Copyright (c) 2025, Oliver M. Drozdowski and Ulrich S. Schwarz (Heidelberg University)

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

#ifndef SMALLAREADAMPINGMODIFIER_HPP_
#define SMALLAREADAMPINGMODIFIER_HPP_

#include <boost/serialization/base_object.hpp>
#include "AbstractCellBasedSimulationModifier.hpp"
#include "ChasteSerialization.hpp"
#include "GeneralizedVolumeConservingForce.hpp"
#include "MonolayerVertexBasedCellPopulation.hpp"
#include "SurfaceTensionForce.hpp"

// Forward declaration prevents circular include chain
template <unsigned SPACE_DIM>
class MonolayerVertexBasedCellPopulation;

/**
 * A struct saving both damping coefficients: The one thresholding the
 * volume correction step and the one thresholding the force
 */
struct DampingCoefficientPair
{
    double volumeCorrectionDamping;
    double forceDamping;
};

/**
 * A modifier that dynamically decreases the mobility of nodes that make
 * up a polygonal face that is very small in area. Effectively we increase
 * the friction for the node in the overdamped limit.
 * With this we can decrease the collapse speed of faces which go to zero
 * area without lowering the time step to zero. If we have vanishing face
 * areas and do not correct with this modifier, the simulation is not able
 * to fully reach area zero.
 *
 * This function needs both pointers to the
 * SurfaceTensionForce/GeneralizedVolumeConservingForce and to the population
 * as the damping is directly given to the force for calculation.
 */
template <unsigned DIM>
class SmallAreaDampingModifier : public AbstractCellBasedSimulationModifier<DIM>
{
    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Boost Serialization method for archiving/checkpointing.
     * Archives the object and its member variables.
     *
     * @param archive  The boost archive.
     * @param version  The current version of this class.
     */
    template <class Archive>
    void serialize(Archive& archive, const unsigned int version)
    {
        archive& boost::serialization::base_object<AbstractCellBasedSimulationModifier<DIM> >(*this);
        archive & mSmallAreaThreshold;
        archive & mForceDampingRatio;
        archive & mVolumeDampingRatio;
        archive & mpPopulation;
        archive & mpForce;
    }

    /**
     * The threshold below which we dampen the vertex movement.
     */
    double mSmallAreaThreshold;

    /**
     * The damping ratio with which we set an upper limit for the node movement through force
     * The force is thresholded below mForceDampingRatio * characteristic length,
     * where the characteristic length is the square root of the small face area.
     * Note that the damping ratio should be smaller for small time steps!
     */
    double mForceDampingRatio;

    /**
     * The damping ratio with which we set an upper limit for the node movement in
     * the volume correction step. The force is thresholded below
     * mVolumeDampingRatio * characteristic length,
     * where the characteristic length is the square root of the small face area.
     * This ratio is independent of time-steps.
     */
    double mVolumeDampingRatio;

    /**
     * Pointer to the Mutable Monolayer Vertex Mesh, to calculate face areas
     */
    MonolayerVertexBasedCellPopulation<DIM>* mpPopulation;

    /**
     * Pointer to the Force to which we give the damping ratios.
     * For now only implemented for Surface Tension Force
     */
    AbstractForce<DIM>* mpForce;

public:
    /**
     * Default constructor.
     *
     * @param pPopulation pointer to MonolayerVertexBasedCellPopulation
     * @param pForce pointer to force
     */
    SmallAreaDampingModifier(MonolayerVertexBasedCellPopulation<DIM>* pPopulation,
                             AbstractForce<DIM>* pForce);

    /**
     * Empty constructor for archiving.
     */
    SmallAreaDampingModifier();

    /**
     * Destructor.
     */
    virtual ~SmallAreaDampingModifier();

    /**
     * Overridden UpdateAtEndOfTimeStep() method.
     *
     * Specify what to do in the simulation at the end of each time step.
     *
     * @param rCellPopulation reference to the cell population
     */
    virtual void UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM, DIM>& rCellPopulation);

    /**
     * Overridden SetupSolve() method.
     *
     * Specify what to do in the simulation before the start of the time loop.
     *
     * @param rCellPopulation reference to the cell population
     * @param outputDirectory the output directory, relative to where Chaste output is stored
     */
    virtual void SetupSolve(AbstractCellPopulation<DIM, DIM>& rCellPopulation, std::string outputDirectory);

    /**
     * Setter for mSmallAreaThreshold
     *
     * @param smallAreaThreshold area threshold
     */
    void SetSmallAreaThreshold(double smallAreaThreshold);

    /**
     * Setter for mForceDampingRatio
     *
     * @param forceDampingRatio
     */
    void SetForceDampingRatio(double forceDampingRatio);

    /**
     * Setter for mForceDampingRatio
     *
     * @param volumeDampingRatio
     */
    void SetVolumeDampingRatio(double volumeDampingRatio);

    /**
     * Getter for mSmallAreaThreshold
     *
     * @return smallAreaThreshold area threshold
     */
    double GetSmallAreaThreshold();

    /**
     * Getter for mForceDampingRatio
     *
     * @return forceDampingRatio
     */
    double GetForceDampingRatio();

    /**
     * Getter for mVolumeDampingRatio
     *
     * @return dampingRatio
     */
    double GetVolumeDampingRatio();

    /**
     * The function that determines the node which need to be damped additionally,
     * then calculates the damping threshold and then gives a map to the force
     * of these nodes with the corresponding thresholds
     *
     * @param pCellPopulation pointer to MonolayerVertexBasedCellPopulation
     */
    void GiveMapOfDampedNodesToForce(MonolayerVertexBasedCellPopulation<DIM>* pCellPopulation);

    /**
     * Overridden OutputSimulationModifierParameters() method.
     * Output any simulation modifier parameters to file.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputSimulationModifierParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(SmallAreaDampingModifier)

#endif /*SMALLAREADAMPINGMODIFIER_HPP_*/
