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

#ifndef RECTANGULAREDGEFORCE_HPP_
#define RECTANGULAREDGEFORCE_HPP_

#include <boost/serialization/base_object.hpp>
#include "ChasteSerialization.hpp"
#include "Exception.hpp"

#include "AbstractCellBasedSimulation.hpp"
#include "AbstractForce.hpp"
#include "MonolayerVertexBasedCellPopulation.hpp"

#include <iostream>

/**
 * A force class for use in MonolayerVertex-based simulations. This force applies
 * a constant force onto nodes belonging to a certain edge. It is meant to be used
 * in rectangular grids, where left, top, right and bottom edges can be defined
 * but it is universally usable if constant forces are supposed to be supplied to
 * a list of known nodes.
 */

template <unsigned DIM>
class RectangularEdgeForce : public AbstractForce<DIM>
{
private:
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
        archive& boost::serialization::base_object<AbstractForce<DIM> >(*this);
        archive & mLeftEdgeNodes;
        archive & mTopEdgeNodes;
        archive & mRightEdgeNodes;
        archive & mBottomEdgeNodes;
        archive & mLeftEdgeForce;
        archive & mTopEdgeForce;
        archive & mRightEdgeForce;
        archive & mBottomEdgeForce;
        archive & mUseInternalCoordinates;
        archive & mpSimulation;
    }

    double CalculateTotalEnergy(AbstractCellPopulation<DIM>& rCellPopulation);

protected:
    /**
     * Vector containing pointers to the nodes belonging to the left edge.
     */
    std::vector<Node<DIM>*> mLeftEdgeNodes;

    /**
     * Vector containing pointers to the nodes belonging to the upper edge.
     */
    std::vector<Node<DIM>*> mTopEdgeNodes;

    /**
     * Vector containing pointers to the nodes belonging to the right edge.
     */
    std::vector<Node<DIM>*> mRightEdgeNodes;

    /**
     * Vector containing pointers to the nodes belonging to the bottom edge.
     */
    std::vector<Node<DIM>*> mBottomEdgeNodes;

    c_vector<double, DIM> mLeftEdgeForce;
    c_vector<double, DIM> mTopEdgeForce;
    c_vector<double, DIM> mRightEdgeForce;
    c_vector<double, DIM> mBottomEdgeForce;

    bool mUseInternalCoordinates = false;

    AbstractCellBasedSimulation<DIM, DIM>* mpSimulation;

public:
    /**
     * Constructor.
     */
    RectangularEdgeForce();

    /**
     * Destructor.
     */
    virtual ~RectangularEdgeForce();

    /**
     * Overridden AddForceContribution() method.
     *
     * Adds the force on each node in the vertex-based cell population according to the edges
     *
     *
     * @param rCellPopulation reference to the cell population
     */
    virtual void AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation);

    /**
     * @return mLeftEdgeNodes
     */
    std::vector<Node<DIM>*> GetLeftEdgeNodes();

    /**
     * @return mTopEdgeNodes
     */
    std::vector<Node<DIM>*> GetTopEdgeNodes();

    /**
     * @return mRightEdgeNodes
     */
    std::vector<Node<DIM>*> GetRightEdgeNodes();

    /**
     * @return mBottomEdgeNodes
     */
    std::vector<Node<DIM>*> GetBottomEdgeNodes();

    /**
     * @param leftEdgeNodes
     */
    void SetLeftEdgeNodes(std::vector<Node<DIM>*> leftEdgeNodes);

    /**
     * @param topEdgeNodes
     */
    void SetTopEdgeNodes(std::vector<Node<DIM>*> topEdgeNodes);

    /**
     * @param rightEdgeNodes
     */
    void SetRightEdgeNodes(std::vector<Node<DIM>*> rightEdgeNodes);

    /**
     * @param bottomEdgeNodes
     */
    void SetBottomEdgeNodes(std::vector<Node<DIM>*> bottomEdgeNodes);

    /**
     * @return mLeftEdgeForce
     */
    c_vector<double, DIM> GetLeftEdgeForce();

    /**
     * @return mTopEdgeForce
     */
    c_vector<double, DIM> GetTopEdgeForce();

    /**
     * @return mRightEdgeForce
     */
    c_vector<double, DIM> GetRightEdgeForce();

    /**
     * @return mBottomEdgeForce
     */
    c_vector<double, DIM> GetBottomEdgeForce();

    /**
     * @param leftEdgeForce
     */
    void SetLeftEdgeForce(c_vector<double, DIM> leftEdgeForce);

    /**
     * @param topEdgeForce
     */
    void SetTopEdgeForce(c_vector<double, DIM> topEdgeForce);

    /**
     * @param rightEdgeForce
     */
    void SetRightEdgeForce(c_vector<double, DIM> rightEdgeForce);

    /**
     * @param bottomEdgeForce
     */
    void SetBottomEdgeForce(c_vector<double, DIM> bottomEdgeForce);

    /**
     * Whether to use internal coordinates. Then the first component of the force
     * vectors coorresponds to the edge normal (average vector to population center),
     * the second component to the average normal of the apical faces and the third
     * component to the cross-product, which should be tangential to the edge.
     *
     * @param useInternalCoordinates
     */
    void UseInternalCoordinates(bool useInternalCoordinates = true);

    /**
     * Overridden OutputForceParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputForceParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(RectangularEdgeForce)

#endif /*RECTANGULAREDGEFORCE_HPP_*/