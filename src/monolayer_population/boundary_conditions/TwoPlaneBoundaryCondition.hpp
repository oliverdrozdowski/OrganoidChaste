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

#ifndef TWOPLANEBOUNDARYCONDITION_HPP_
#define TWOPLANEBOUNDARYCONDITION_HPP_

#include "AbstractCellPopulationBoundaryCondition.hpp"

#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>
#include "ChasteSerialization.hpp"

/**
 * A simply supported cell population boundary condition class, which fixes the movement
 * of apical and basal boundary nodes to specified planes in the domain. It thus implements
 * a simply supported boundary condition, where we, however, do not allow for a free height.
 * The height of the lateral face will be fixed by the inclination angle of the cell sheet.
 */
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM = ELEMENT_DIM>
class TwoPlaneBoundaryCondition : public AbstractCellPopulationBoundaryCondition<ELEMENT_DIM, SPACE_DIM>
{
private:
    /**
     * A point on the basal boundary plane.
     */
    c_vector<double, SPACE_DIM> mBasalPointOnSubstrate;

    /**
     * The outward-facing unit normal vector to the basal boundary plane.
     */
    c_vector<double, SPACE_DIM> mBasalNormalToSubstrate;

    /**
     * A point on the apical boundary plane.
     */
    c_vector<double, SPACE_DIM> mApicalPointOnSubstrate;

    /**
     * The outward-facing unit normal vector to the apical boundary plane.
     */
    c_vector<double, SPACE_DIM> mApicalNormalToSubstrate;

    /**
     * Whether to jiggle the cells on the basal bottom surface, initialised to false
     * in the constructor.
     */
    bool mUseJiggledNodesOnSubstrate;

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Serialize the object and its member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template <class Archive>
    void serialize(Archive& archive, const unsigned int version)
    {
        archive& boost::serialization::base_object<AbstractCellPopulationBoundaryCondition<ELEMENT_DIM, SPACE_DIM> >(*this);
        // archive & mUseJiggledNodesOnSubstrate;
    }

public:
    /**
     * Constructor.
     *
     * @param pCellPopulation pointer to the cell population
     * @param point a point on the boundary plane
     * @param normal the outward-facing unit normal vector to the boundary plane
     */
    TwoPlaneBoundaryCondition(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation,
                              c_vector<double, SPACE_DIM> basalPoint,
                              c_vector<double, SPACE_DIM> basalNormal,
                              c_vector<double, SPACE_DIM> apicalPoint,
                              c_vector<double, SPACE_DIM> apicalNormal);

    /**
     * @return #mBasalPointOnSubstrate.
     */
    const c_vector<double, SPACE_DIM>& rGetBasalPointOnSubstrate() const;

    /**
     * @return #mBasalNormalToSubstrate.
     */
    const c_vector<double, SPACE_DIM>& rGetBasalNormalToSubstrate() const;

    /**
     * @return #mApicalPointOnSubstrate.
     */
    const c_vector<double, SPACE_DIM>& rGetApicalPointOnSubstrate() const;

    /**
     * @return #mApicalNormalToSubstrate.
     */
    const c_vector<double, SPACE_DIM>& rGetApicalNormalToSubstrate() const;

    /**
     * Set method for mUseJiggledNodesOnSubstrate
     *
     * @param useJiggledNodesOnSubstrate whether to jiggle the nodes on the surface of the plane, can help stop overcrowding on plane.
     */
    void SetUseJiggledNodesOnSubstrate(bool useJiggledNodesOnSubstrate);

    /** @return #mUseJiggledNodesOnSubstrate. */
    bool GetUseJiggledNodesOnSubstrate();

    /**
     * Overridden ImposeBoundaryCondition() method.
     *
     * Apply the cell population boundary conditions.
     *
     * @param rOldLocations the node locations before any boundary conditions are applied
     */
    void ImposeBoundaryCondition(const std::map<Node<SPACE_DIM>*, c_vector<double, SPACE_DIM> >& rOldLocations);

    /**
     * Overridden VerifyBoundaryCondition() method.
     * Verify the boundary conditions have been applied.
     * This is called after ImposeBoundaryCondition() to ensure the condition is still satisfied.
     *
     * @return whether the boundary conditions are satisfied.
     */
    bool VerifyBoundaryCondition();

    /**
     * Overridden OutputCellPopulationBoundaryConditionParameters() method.
     * Output cell population boundary condition parameters to file.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputCellPopulationBoundaryConditionParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(TwoPlaneBoundaryCondition)

namespace boost
{
namespace serialization
{
    /**
     * Serialize information required to construct a TwoPlaneBoundaryCondition.
     */
    template <class Archive, unsigned ELEMENT_DIM, unsigned SPACE_DIM>
    inline void save_construct_data(
        Archive& ar, const TwoPlaneBoundaryCondition<ELEMENT_DIM, SPACE_DIM>* t, const unsigned int file_version)
    {
        // Save data required to construct instance
        const AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* const p_cell_population = t->GetCellPopulation();
        ar << p_cell_population;

        // Archive c_vectors one component at a time
        c_vector<double, SPACE_DIM> basalPoint = t->rGetBasalPointOnSubstrate();
        for (unsigned i = 0; i < SPACE_DIM; i++)
        {
            ar << basalPoint[i];
        }
        c_vector<double, SPACE_DIM> basalNormal = t->rGetBasalNormalToSubstrate();
        for (unsigned i = 0; i < SPACE_DIM; i++)
        {
            ar << basalNormal[i];
        }

        c_vector<double, SPACE_DIM> apicalPoint = t->rGetApicalPointOnSubstrate();
        for (unsigned i = 0; i < SPACE_DIM; i++)
        {
            ar << apicalPoint[i];
        }
        c_vector<double, SPACE_DIM> apicalNormal = t->rGetApicalNormalToSubstrate();
        for (unsigned i = 0; i < SPACE_DIM; i++)
        {
            ar << apicalNormal[i];
        }
    }

    /**
     * De-serialize constructor parameters and initialize a TwoPlaneBoundaryCondition.
     */
    template <class Archive, unsigned ELEMENT_DIM, unsigned SPACE_DIM>
    inline void load_construct_data(
        Archive& ar, TwoPlaneBoundaryCondition<ELEMENT_DIM, SPACE_DIM>* t, const unsigned int file_version)
    {
        // Retrieve data from archive required to construct new instance
        AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* p_cell_population;
        ar >> p_cell_population;

        // Archive c_vectors one component at a time
        c_vector<double, SPACE_DIM> basalPoint;
        for (unsigned i = 0; i < SPACE_DIM; i++)
        {
            ar >> basalPoint[i];
        }
        c_vector<double, SPACE_DIM> basalNormal;
        for (unsigned i = 0; i < SPACE_DIM; i++)
        {
            ar >> basalNormal[i];
        }

        c_vector<double, SPACE_DIM> apicalPoint;
        for (unsigned i = 0; i < SPACE_DIM; i++)
        {
            ar >> apicalPoint[i];
        }
        c_vector<double, SPACE_DIM> apicalNormal;
        for (unsigned i = 0; i < SPACE_DIM; i++)
        {
            ar >> apicalNormal[i];
        }

        // Invoke inplace constructor to initialise instance
        ::new (t) TwoPlaneBoundaryCondition<ELEMENT_DIM, SPACE_DIM>(p_cell_population, basalPoint, basalNormal,
                                                                    apicalPoint, apicalNormal);
    }
} // namespace serialization
} // namespace boost

#endif /*TWOPLANEBOUNDARYCONDITION_HPP_*/
