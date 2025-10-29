#ifndef MIDPLANESPHERICALBOUNDARYCONDITION_HPP_
#define MIDPLANESPHERICALBOUNDARYCONDITION_HPP_

#include "AbstractCellPopulationBoundaryCondition.hpp"

#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>
#include "ChasteSerialization.hpp"

/**
 * A cylindrical cell population boundary condition class, which fixes the movement
 * of the midpoint of apical and basal boundary nodes to a cylindrical surface.
 * On the surface of the cylinder the center points are, however, free to move.
 */
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM = ELEMENT_DIM>
class MidPlaneSphericalBoundaryCondition : public AbstractCellPopulationBoundaryCondition<ELEMENT_DIM, SPACE_DIM>
{
private:
    /**
     * Center point
     */
    c_vector<double, SPACE_DIM> mMidPoint;

    /**
     * The radius of the sphere
     */
    double mRadius;

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
     * @param mid_point mid point of sphere
     * @param radius
     */
    MidPlaneSphericalBoundaryCondition(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation,
                                       c_vector<double, SPACE_DIM> mid_point,
                                       double radius);

    /**
     * @return #mMidPoint.
     */
    const c_vector<double, SPACE_DIM>& rGetMidPoint() const;

    /**
     * @return #mRadius.
     */
    double GetRadius() const;

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
EXPORT_TEMPLATE_CLASS_ALL_DIMS(MidPlaneSphericalBoundaryCondition)

namespace boost
{
namespace serialization
{
    /**
     * Serialize information required to construct a MidPlaneSphericalBoundaryCondition.
     */
    template <class Archive, unsigned ELEMENT_DIM, unsigned SPACE_DIM>
    inline void save_construct_data(
        Archive& ar, const MidPlaneSphericalBoundaryCondition<ELEMENT_DIM, SPACE_DIM>* t, const unsigned int file_version)
    {
        // Save data required to construct instance
        const AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* const p_cell_population = t->GetCellPopulation();
        ar << p_cell_population;

        // Archive c_vectors one component at a time
        c_vector<double, SPACE_DIM> midPoint = t->rGetMidPoint();
        for (unsigned i = 0; i < SPACE_DIM; i++)
        {
            ar << midPoint[i];
        }
        double radius = t->GetRadius();
        ar << radius;
    }

    /**
     * De-serialize constructor parameters and initialize a MidPlaneSphericalBoundaryCondition.
     */
    template <class Archive, unsigned ELEMENT_DIM, unsigned SPACE_DIM>
    inline void load_construct_data(
        Archive& ar, MidPlaneSphericalBoundaryCondition<ELEMENT_DIM, SPACE_DIM>* t, const unsigned int file_version)
    {
        // Retrieve data from archive required to construct new instance
        AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* p_cell_population;
        ar >> p_cell_population;

        // Archive c_vectors one component at a time
        c_vector<double, SPACE_DIM> midPoint;
        for (unsigned i = 0; i < SPACE_DIM; i++)
        {
            ar >> midPoint[i];
        }
        double radius;
        ar >> radius;

        // Invoke inplace constructor to initialise instance
        ::new (t) MidPlaneSphericalBoundaryCondition<ELEMENT_DIM, SPACE_DIM>(p_cell_population, midPoint, radius);
    }
} // namespace serialization
} // namespace boost

#endif /*MIDPLANESPHERICALBOUNDARYCONDITION_HPP_*/
