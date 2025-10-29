#ifndef ABSTRACTMOVINGBOUNDARY_HPP_
#define ABSTRACTMOVINGBOUNDARY_HPP_

#include "AbstractMoveProfile.hpp"

#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>

/**
 * A class to implement moving boundary conditions in 3D vertex model
 * simulations. The class uses an AbstractMoveProfile to implement different
 * temporal profiles which modify the boundary conditions. The moving boundaries
 * are implemented via a force acting on nodes, if they experience the boundary.
 */
template <unsigned DIM>
class AbstractMovingBoundary
{
protected:
    /**
     * the force constant k in F= -k * dx
     */
    double mForceConstant{};

    /**
     * Move profile z(t) for the moving boundary
     */
    boost::shared_ptr<AbstractMoveProfile<DIM> > mMoveProfile{};

    /**
     * Position of the boundary on respect to the cell population. Can be above or
     * below
     */
    bool mAbove{};

    /**
     * Total force applied by boundary on the population in one step
     */
    c_vector<double, DIM> mForceAppliedOnPopulation;

public:
    /**
     * Constructor for simple plate (can either be on top or bottom of the sphere)
     */
    AbstractMovingBoundary(double force_constant, boost::shared_ptr<AbstractMoveProfile<DIM> > move_profile, bool above);

    /**
     * Deconstructor
     */
    virtual ~AbstractMovingBoundary();

    /**
     * Calculate the force on a node based on its position. Return 0 when node is
     * not beyond the plate
     */
    virtual c_vector<double, DIM>
    CalculateForceOnNode(c_vector<double, DIM> node_location, double damping_constant, double dt) = 0;

    /**
     * Add force contribution
     */
    void AddAppliedForceOnPopulation(c_vector<double, DIM> force);

    /**
     * Reset applied force on boundary to 0
     */
    void ResetAppliedForceOnPopulation();

    /**
     * Get the total applied force from the boundary on the population
     */
    c_vector<double, DIM> GetAppliedForceOnPopulation();

    void SetMoveProfile(boost::shared_ptr<AbstractMoveProfile<DIM> > move_profile);

    void SetForceConstant(double k);

    double GetForceConstant();

    double GetZPosition();
};

#endif // ABSTRACTMOVINGBOUNDARY_HPP