#include "LinearSystem.hpp"

#include "MonolayerVertexBasedCellPopulation.hpp"
#include "MonolayerVertexElement.hpp"
#include "SphereMovingBoundary.hpp"

template <unsigned DIM>
SphereMovingBoundary<DIM>::SphereMovingBoundary(
    double force_constant, boost::shared_ptr<AbstractMoveProfile<DIM> > move_profile, bool above,
    double radius)
        : AbstractMovingBoundary<DIM>(force_constant, move_profile, above),
          mRadius(radius) {}

template <unsigned DIM>
SphereMovingBoundary<DIM>::~SphereMovingBoundary() {}

template <unsigned DIM>
c_vector<double, DIM> SphereMovingBoundary<DIM>::CalculateForceOnNode(
    c_vector<double, DIM> node_location, double damping_constant, double dt)
{
    EXCEPTION("SphereMovingBoundary<DIM>::CalculateForceOnNode(node_location) is "
              "only implemented for DIM = 3");
}

template <>
c_vector<double, 3> SphereMovingBoundary<3>::CalculateForceOnNode(
    c_vector<double, 3> node_location, double damping_constant, double dt)
{

    c_vector<double, 3> force = zero_vector<double>(3);

    double z_direction = this->mAbove ? -1.0 : 1.0;
    double boundary_z_position = this->mMoveProfile->GetZPosition();

    c_vector<double, 3> indenter_center = zero_vector<double>(3);
    indenter_center[2] = boundary_z_position - z_direction * mRadius;

    c_vector<double, 3> center_to_node = node_location - indenter_center;

    double distance = norm_2(center_to_node);

    if (distance < mRadius)
    {
        force = damping_constant / dt * (mRadius / distance - 1.0) * center_to_node;
    }

    return force;
}

// Explicit instantiation
template class SphereMovingBoundary<1>;
template class SphereMovingBoundary<2>;
template class SphereMovingBoundary<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
