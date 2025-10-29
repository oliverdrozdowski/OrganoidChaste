#ifndef STATICMOVEPROFILE_HPP_
#define STATICMOVEPROFILE_HPP_

#include "AbstractMoveProfile.hpp"

#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>

template <unsigned DIM>
class StaticMoveProfile : public AbstractMoveProfile<DIM>
{
public:
    /**
     * Constructor for simple plate (can either be on top or bottom of the sphere)
     */

    StaticMoveProfile(double position);

    /**
     * Get the z position at the current simulation time
     */
    virtual double GetZPosition();

    /**
     * Update position
     */
    void UpdatePosition(double position);
};

#endif // STATICMOVEVEPROFILE_HPP_