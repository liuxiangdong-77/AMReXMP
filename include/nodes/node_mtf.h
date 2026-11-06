// mtf_quantities.h
#ifndef MPM_MTF_QUANTITIES_H
#define MPM_MTF_QUANTITIES_H

#include <Eigen/Dense>

namespace mpm {


struct RiemannInvariants{
  //! Outgoing wave Riemann invariant (R1 = v + σ/(ρC))
  double R1 = 0.0;

  //! Incoming wave Riemann invariant (R2 = v - σ/(ρC))
  double R2 = 0.0;

  RiemannInvariants() = default;

  //! Constructor with values
  //! \param[in] r1 Outgoing wave invariant
  //! \param[in] r2 Incoming wave invariant
  RiemannInvariants(double r1, double r2) : R1(r1), R2(r2) {}

  void setZero() {
    R1 = 0.0;
    R2 = 0.0;
  }

  bool isValid() const {
    return !std::isnan(R1) && !std::isnan(R2);
  }

  //! String representation for debugging
  std::string toString() const {
    return "R1=" + std::to_string(R1) + ", R2=" + std::to_string(R2);
  }
};

//! MTF physical quantities container for extrapolation points
//! \brief Contains interpolated physical quantities at extrapolation points
//! \details This class stores the physical quantities obtained by interpolating
//! from particles at specific time steps to extrapolation points
template <unsigned Tdim>
struct MTFQuantities {
  //! Velocity vector
  Eigen::Matrix<double, Tdim, 1> velocity = Eigen::Matrix<double, Tdim, 1>::Zero();
  
  //! Stress in Voigt notation (6 components)
  Eigen::Matrix<double, 6, 1> stress = Eigen::Matrix<double, 6, 1>::Zero();
  
  //! Density
  double density = 0.0;
  
  //! Coordinates of the extrapolation point
  Eigen::Matrix<double, Tdim, 1> coordinates = Eigen::Matrix<double, Tdim, 1>::Zero();
  
  //! Time step offset for this interpolation
  //! 0 = current time step, -1 = previous time step, -2 = two steps back, etc.
  int time_offset = 0;

  //! Default constructor
  MTFQuantities() = default;
  
  //! Constructor with time offset
  //! \param[in] offset Time step offset
  MTFQuantities(int offset) : time_offset(offset) {}

  //! Constructor with coordinates
  //! \param[in] coords Coordinates of the point
  explicit MTFQuantities(const Eigen::Matrix<double, Tdim, 1>& coords) 
    : coordinates(coords) {}
  
  //! Constructor with full initialization
  //! \param[in] coords Coordinates
  //! \param[in] vel Velocity vector
  //! \param[in] str Stress tensor (Voigt)
  //! \param[in] dens Density
  //! \param[in] offset Time offset
  MTFQuantities(const Eigen::Matrix<double, Tdim, 1>& coords,
                const Eigen::Matrix<double, Tdim, 1>& vel,
                const Eigen::Matrix<double, 6, 1>& str,
                double dens, int offset = 0)
    : velocity(vel), stress(str), density(dens), 
      coordinates(coords), time_offset(offset) {}
  
  
  //! Reset all quantities to zero
  void setZero() {
    velocity.setZero();
    stress.setZero();
    density = 0.0;
    coordinates.setZero();
  }
  
  //! Check if quantities are valid (basic validation)
  bool isValid() const {
    return density > 0.0 && 
           !std::isnan(velocity.norm()) && 
           !std::isnan(stress.norm());
  }
  
  //! Get time identifier string
  std::string timeIdentifier() const {
    if (time_offset == 0) return "t0";
    return "t" + std::to_string(time_offset); // e.g., "t-1", "t-2"
  }
  
  //! String representation for debugging
  std::string toString() const {
    std::stringstream ss;
    ss << "MTFQuantities[" << timeIdentifier() << "]: "
       << "coord=[" << coordinates.transpose() << "], "
       << "vel=[" << velocity.transpose() << "], "
       << "density=" << density << ", "
       << "stress_norm=" << stress.norm();
    return ss.str();
  }
  
};

}  // namespace mpm

#endif  // MPM_MTF_QUANTITIES_H