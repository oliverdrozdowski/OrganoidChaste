#ifndef PlateMovingBoundary_HPP_
#define PlateMovingBoundary_HPP_

#include "AbstractMoveProfile.hpp"
#include "AbstractMovingBoundary.hpp"

#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>

template <unsigned DIM>
class PlateMovingBoundary : public AbstractMovingBoundary<DIM>
{
private:
    /**
     * Boolean on if to use a hard potential as the boundary. The alternative case is a "soft" material which forces vertices that passed the boundary back to the boundary position
     */
    bool mUseHardPotential;

public:
    /**
     * Constructor for simple plate (can either be on top or bottom of the sphere)
     */
    PlateMovingBoundary(double force_constant,
                        boost::shared_ptr<AbstractMoveProfile<DIM> > move_profile, bool above, bool use_hard_potential = true);

    /**
     * Deconstructor
     */
    virtual ~PlateMovingBoundary();

    /**
     * Calculate the force on a node based on its position. Return 0 when node is
     * not beyond the plate
     */
    virtual c_vector<double, DIM>
    CalculateForceOnNode(c_vector<double, DIM> node_location, double damping_constant, double dt);
};

#endif