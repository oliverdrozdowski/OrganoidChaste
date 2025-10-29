#ifndef SPHEREMOVINGBOUNDARY_HPP_
#define SPHEREMOVINGBOUNDARY_HPP_

#include "AbstractMoveProfile.hpp"
#include "AbstractMovingBoundary.hpp"

#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>

template <unsigned DIM>
class SphereMovingBoundary : public AbstractMovingBoundary<DIM>
{
private:
    double mRadius;

public:
    /**
     * Constructor for simple plate (can either be on top or bottom of the sphere)
     */
    SphereMovingBoundary(double force_constant,
                         boost::shared_ptr<AbstractMoveProfile<DIM> > move_profile, bool above,
                         double radius);

    /**
     * Deconstructor
     */
    virtual ~SphereMovingBoundary();

    /**
     * Calculate the force on a node based on its position. Return 0 when node is
     * not beyond the plate
     */
    virtual c_vector<double, DIM>
    CalculateForceOnNode(c_vector<double, DIM> node_location, double damping_constant, double dt);
};

#endif