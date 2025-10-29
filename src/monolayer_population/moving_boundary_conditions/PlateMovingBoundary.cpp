#include "LinearSystem.hpp"

#include "MonolayerVertexBasedCellPopulation.hpp"
#include "MonolayerVertexElement.hpp"
#include "PlateMovingBoundary.hpp"
#include "SimulationTime.hpp"

template <unsigned DIM>
PlateMovingBoundary<DIM>::PlateMovingBoundary(double force_constant, boost::shared_ptr<AbstractMoveProfile<DIM> > move_profile, bool above, bool use_hard_potential)
        : AbstractMovingBoundary<DIM>(force_constant, move_profile, above), mUseHardPotential(use_hard_potential) {}

template <unsigned DIM>
PlateMovingBoundary<DIM>::~PlateMovingBoundary() {}

template <unsigned DIM>
c_vector<double, DIM> PlateMovingBoundary<DIM>::CalculateForceOnNode(
    c_vector<double, DIM> node_location, double damping_constant, double dt)
{
    // normal of the plane is defined in the direction of the cells
    double z_normal_to_plane = this->mAbove ? -1.0 : 1.0;

    c_vector<double, DIM> normal_to_plane = zero_vector<double>(DIM);
    normal_to_plane[DIM - 1] = z_normal_to_plane;

    c_vector<double, DIM> force = zero_vector<double>(DIM);

    c_vector<double, DIM> z_position_vector = zero_vector<double>(DIM);
    z_position_vector[DIM - 1] = this->mMoveProfile->GetZPosition();

    double signed_distance = inner_prod(node_location - z_position_vector, normal_to_plane);

    if (mUseHardPotential)
    {
        if (signed_distance < 0.0)
        {
            // if the node is on the "wrong" side of the boundary, exert a force on it that puts it on the "correct" side with the exact distance
            force[DIM - 1] = 2.0 * damping_constant / dt * signed_distance * z_normal_to_plane * -1.0;

            //  EXCEPTION("Node is outside of moving boundary");
        }
        else
        {
            force[DIM - 1] = 12.0 * pow(this->mForceConstant, 12) * pow(signed_distance, -13) * z_normal_to_plane;
        }
    }
    else
    {
        if (signed_distance <= 0.0)
        {
            force[DIM - 1] = (-1.0) * damping_constant / dt * signed_distance * z_normal_to_plane;
        }
    }

    return force;
}

// Explicit instantiation
template class PlateMovingBoundary<1>;
template class PlateMovingBoundary<2>;
template class PlateMovingBoundary<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
