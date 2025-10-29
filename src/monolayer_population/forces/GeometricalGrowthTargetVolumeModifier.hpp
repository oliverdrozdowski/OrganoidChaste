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

#ifndef GEOMETRICALGROWTHTARGETVOLUMEMODIFIER_HPP_
#define GEOMETRICALGROWTHTARGETVOLUMEMODIFIER_HPP_

#include <boost/serialization/base_object.hpp>
#include "AbstractTargetVolumeModifier.hpp"
#include "ChasteSerialization.hpp"
#include "MonolayerVertexBasedCellPopulation.hpp"

// Forward declaration prevents circular include chain
template <unsigned SPACE_DIM>
class MonolayerVertexBasedCellPopulation;

/**
 * A target area modifier class in which the target area of a cell grows linearly,
 * starting from mReferenceTargetArea, over a prescribed duration.
 *
 * If used with a phase-based cell-cycle model (such as FixedG1GenerationalCellCycleModel),
 * the target area of a cell doubles from the value mReferenceTargetArea over
 * the course of the cell's G2 phase. This rule is based on the description of
 * cell growth and division proposed by Farhadifar et al
 * (The influence of cell mechanics, cell-cell interactions,
 * and proliferation on epithelial packing, Current Biology 17(24):2095-2104, 2007,
 * http://dx.doi.org/10.1016/j.cub.2007.11.049).
 *
 * If used with a non-phase-based cell-cycle model, the target area of a cell
 * increases linearly at a rate mGrowthRate as soon as the cell's age exceeds mAgeToStartGrowing.
 *
 * Here mReferenceTargetArea, mAgeToStartGrowing and mGrowthRate are settable member variables.
 * The default value of mReferenceTargetArea is 1.0 and the default value of mAgeToStartGrowing
 * and mGrowthRate is DOUBLE_UNSET.
 *
 * If mAgeToStartGrowing is set by the user, then mGrowthRate must also be set;
 * note that in this case, these values are used to prescribe target area growth
 * as described earlier, regardless of whether a phase-based cell-cycle model is present.
 *
 * Addtionally we slowly adapt the cell volume after volume-changing T1 transitions to minimize
 * the occurence of intersections. Here mT1AdaptationDuration determines the time it takes
 * for the new changed volume to be corrected linearly to the target area.
 *
 */
template <unsigned DIM>
class GeometricalGrowthTargetVolumeModifier : public AbstractTargetVolumeModifier<DIM>
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
        archive& boost::serialization::base_object<AbstractTargetVolumeModifier<DIM> >(*this);
        archive & mAgeToStartGrowing;
        archive & mGrowthRate;
    }

    /**
     * The age of a proliferating cell at which its target area should start growing.
     * Defaults to DOUBLE_UNSET. If this variable is set using SetAgeToStartGrowing(),
     * then it is used regardless of whether a phase-based cell-cycle model is used.
     */
    double mAgeToStartGrowing;

    /**
     * The growth rate of a proliferating cell's target area, when it is growing.
     * Defaults to DOUBLE_UNSET. This variable must be set if mAgeToStartGrowing is set,
     * otherwise an exception is thrown.
     */
    double mGrowthRate;

    /**
     * The duration that a cell's target volume needs to increase from the volume
     * after T1 transitions to mReferenceVolumeArea  Defaults to zero, i.e. instanteneously.
     */
    double mT1AdaptationDuration;

    /**
     * Pointer to the Mutable Monolayer Vertex Mesh, to get access to protorosette elements
     */
    MonolayerVertexBasedCellPopulation<DIM>* mpPopulation;

public:
    /**
     * Default constructor.
     */
    GeometricalGrowthTargetVolumeModifier(MonolayerVertexBasedCellPopulation<DIM>* pPopulation);

    /**
     * Empty constructor for archiving.
     */
    GeometricalGrowthTargetVolumeModifier();

    /**
     * Destructor.
     */
    virtual ~GeometricalGrowthTargetVolumeModifier();

    /**
     * Helper method to update the target volume property of an individual cell.
     *
     * @param pCell pointer to a cell
     */
    void UpdateTargetVolumeOfCell(const CellPtr pCell);

    /**
     * @return #mAgeToStartGrowing.
     */
    double GetAgeToStartGrowing();

    /**
     * Set #mAgeToStartGrowing.
     *
     * @param ageToStartGrowing the new value of #mAgeToStartGrowing
     */
    void SetAgeToStartGrowing(double ageToStartGrowing);

    /**
     * @return #mGrowthRate
     */
    double GetGrowthRate();

    /**
     * Set #mGrowthRate.
     *
     * @param growthRate the new value of #mGrowthRate
     */
    void SetGrowthRate(double growthRate);

    /**
     * @return #mT1AdaptationDuration
     */
    double GetT1AdaptationDuration();

    /**
     * Set #mGrowthDuration.
     *
     * @param growthDuration the new value of #mT1AdaptationDuration
     */
    void SetT1AdaptationDuration(double t1AdaptationDuration);

    /**
     * Overridden OutputSimulationModifierParameters() method.
     * Output any simulation modifier parameters to file.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputSimulationModifierParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(GeometricalGrowthTargetVolumeModifier)

#endif /*GEOMETRICALGROWTHTARGETVOLUMEMODIFIER_HPP_*/