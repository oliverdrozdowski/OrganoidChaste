#ifndef FIXNODEDIMENSIONSBOUNDARYCONDITION_HPP_
#define FIXNODEDIMENSIONSBOUNDARYCONDITION_HPP_

#include "AbstractCellPopulationBoundaryCondition.hpp"

#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>
#include "ChasteSerialization.hpp"

/**
 * A boundary condition which fixes the x and/or the y and/or the z position of nodes.
 * A node list is supplied in the beginning and for these nodes the positions in the dimensions will be kept constant
 */
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM = ELEMENT_DIM>
class FixNodeDimensionsBoundaryCondition : public AbstractCellPopulationBoundaryCondition<ELEMENT_DIM, SPACE_DIM>
{
private:
    /**
     * Vector of nodes to consider
     */
    std::vector<Node<SPACE_DIM>*> mNodesToConsider;

    /**
     * Original positions
     */
    std::map<Node<SPACE_DIM>*, c_vector<double, SPACE_DIM> > mOriginalNodePositions;

    bool mKeepConstantX;
    bool mKeepConstantY;
    bool mKeepConstantZ;

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
        archive& mNodesToConsider;
        archive& mOriginalNodePositions;
        archive& mKeepConstantX;
        archive& mKeepConstantY;
        archive& mKeepConstantZ;
    }

public:
    /**
     * Constructor.
     *
     * @param pCellPopulation pointer to population
     * @param nodesToConsider vector of nodes to consider
     * @param keepConstantX whether to keep x constant
     * @param keepConstantY whether to keep y constant
     * @param keepConstantZ whether to keep z constant
     */
    FixNodeDimensionsBoundaryCondition(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation,
                                       std::vector<Node<SPACE_DIM>*> nodesToConsider,
                                       bool keepConstantX,
                                       bool keepConstantY,
                                       bool keepConstantZ);

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

    /**
     * @param keepConstantX
     */
    void SetKeepConstantX(bool keepConstantX)
    {
        mKeepConstantX = keepConstantX;
    }

    /**
     * @param keepConstantY
     */
    void SetKeepConstantY(bool keepConstantY)
    {
        mKeepConstantY = keepConstantY;
    }

    /**
     * @param keepConstantZ
     */
    void SetKeepConstantZ(bool keepConstantZ)
    {
        mKeepConstantZ = keepConstantZ;
    }

    /**
     * @param originalNodePositions
     */
    void SetOriginalNodePositions(std::map<Node<SPACE_DIM>*, c_vector<double, SPACE_DIM> > nodePositions)
    {
        mOriginalNodePositions = nodePositions;
    }

    /**
     * @param nodesToConsider
     */
    void SetNodesToConsider(std::vector<Node<SPACE_DIM>*> nodes)
    {
        mNodesToConsider = nodes;
    }

    /**
     * @return keepConstantX
     */
    bool GetKeepConstantX() const
    {
        return mKeepConstantX;
    }

    /**
     * @return keepConstantY
     */
    bool GetKeepConstantY() const
    {
        return mKeepConstantY;
    }

    /**
     * @return keepConstantZ
     */
    bool GetKeepConstantZ() const
    {
        return mKeepConstantZ;
    }

    /**
     * @return originalNodePositions
     */
    const std::map<Node<SPACE_DIM>*, c_vector<double, SPACE_DIM> > GetOriginalNodePositions() const
    {
        return mOriginalNodePositions;
    }

    /**
     * @return nodesToConsider
     */
    const std::vector<Node<SPACE_DIM>*> GetNodesToConsider() const
    {
        return mNodesToConsider;
    }
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(FixNodeDimensionsBoundaryCondition)

namespace boost
{
namespace serialization
{
    /**
     * Serialize information required to construct a FixNodeDimensionsBoundaryCondition.
     */
    template <class Archive, unsigned ELEMENT_DIM, unsigned SPACE_DIM>
    inline void save_construct_data(
        Archive& ar, const FixNodeDimensionsBoundaryCondition<ELEMENT_DIM, SPACE_DIM>* t, const unsigned int file_version)
    {
        // Save data required to construct instance
        const AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* const p_cell_population = t->GetCellPopulation();
        ar << p_cell_population;

        bool keepX = t->GetKeepConstantX();
        bool keepY = t->GetKeepConstantY();
        bool keepZ = t->GetKeepConstantZ();
        ar << keepX;
        ar << keepY;
        ar << keepZ;

        std::map<Node<SPACE_DIM>*, c_vector<double, SPACE_DIM> > originalNodePos = t->GetOriginalNodePositions();
        ar << originalNodePos;

        std::vector<Node<SPACE_DIM>*> nodesToConsider = t->GetNodesToConsider();
        ar << nodesToConsider;
    }

    /**
     * De-serialize constructor parameters and initialize a FixNodeDimensionsBoundaryCondition.
     */
    template <class Archive, unsigned ELEMENT_DIM, unsigned SPACE_DIM>
    inline void load_construct_data(
        Archive& ar, FixNodeDimensionsBoundaryCondition<ELEMENT_DIM, SPACE_DIM>* t, const unsigned int file_version)
    {
        // Retrieve data from archive required to construct new instance
        AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* p_cell_population;
        ar >> p_cell_population;

        bool keepX, keepY, keepZ;
        ar >> keepX;
        ar >> keepY;
        ar >> keepZ;

        std::map<Node<SPACE_DIM>*, c_vector<double, SPACE_DIM> > originalNodePos;
        ar >> originalNodePos;

        std::vector<Node<SPACE_DIM>*> nodesToConsider;
        ar >> nodesToConsider;

        // Invoke inplace constructor to initialise instance
        ::new (t) FixNodeDimensionsBoundaryCondition<ELEMENT_DIM, SPACE_DIM>(p_cell_population, nodesToConsider,
                                                                             keepX, keepY, keepZ);
    }
} // namespace serialization
} // namespace boost

#endif /*FIXNODEDIMENSIONSBOUNDARYCONDITION_HPP_*/
