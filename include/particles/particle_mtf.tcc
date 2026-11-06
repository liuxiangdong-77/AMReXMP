//! Check if particle has MTF history variables
template <unsigned Tdim>
bool mpm::Particle<Tdim>::has_mtf_history() const {
  if (this->material(mpm::ParticlePhase::Solid) != nullptr) {
    auto state_vars = this->material(mpm::ParticlePhase::Solid)->state_variables();
    for (const auto& var : state_vars) {
      if (var.find("_t") != std::string::npos) {
        return true;
      }
    }
  }
  return false;
}

//! Get history variable at specific time step for MTF
template <unsigned Tdim>
double mpm::Particle<Tdim>::get_history_variable(
    const std::string& var_name, 
    unsigned time_step,
    unsigned phase) const {
  
  auto material_ptr = this->material(phase);
  if (material_ptr == nullptr || phase >= state_variables_.size()) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  // 直接使用基类接口，不需要知道具体类型
  return material_ptr->get_history_variable(&state_variables_[phase], 
                                           var_name, time_step);
}

//! Get history velocity for MTF
template <unsigned Tdim>
typename mpm::Particle<Tdim>::VectorDim mpm::Particle<Tdim>::get_history_velocity(
    unsigned time_step,
    unsigned phase) const {
  
  VectorDim history_vel = VectorDim::Zero();
  
  for (unsigned i = 0; i < Tdim; ++i) {
    double vel_component = get_history_variable("velocity_" + std::to_string(i), 
                                               time_step, phase);
    if (!std::isnan(vel_component)) {
      history_vel(i) = vel_component;
    }
  }
  
  return history_vel;
}

// //! Get history stress for MTF
// template <unsigned Tdim>
// Eigen::Matrix<double, 6, 1> mpm::Particle<Tdim>::get_history_stress(
//     unsigned time_step,
//     unsigned phase) const {
//   Eigen::Matrix<double, 6, 1> history_stress = Eigen::Matrix<double, 6, 1>::Zero();
//   for (unsigned i = 0; i < 6; ++i) {
//     double stress_component = get_history_variable("stress_" + std::to_string(i), 
//                                                   time_step, phase);
//     if (!std::isnan(stress_component)) {
//       history_stress(i) = stress_component;
//     }
//   }
//   return history_stress;
// }

// //! Get history density for MTF
// template <unsigned Tdim>
// double mpm::Particle<Tdim>::get_history_density(
//     unsigned time_step,
//     unsigned phase) const {
//   return get_history_variable("density", time_step, phase);
// }


//! Advance MTF history variables (should be called at end of time step)
template <unsigned Tdim>
void mpm::Particle<Tdim>::advance_mtf_history(unsigned long long time_step) {
  for (unsigned phase = 0; phase < material_.size(); ++phase) {
    auto material_ptr = material_[phase];
    if (material_ptr != nullptr && phase < state_variables_.size()) {
      // 存储当前时间步的状态
      console_->debug("Particle {}: Storing history for phase {}", 
                         this->id(), phase);
      // 直接使用基类接口！
      material_ptr->advance_history_variables(&state_variables_[phase], this,
                                              time_step);
      // // 验证存储是否成功
      // double density_after = get_history_density(0, phase);
      // console_->debug("Particle {}: After advance, density_t0 = {}", 
      //                    this->id(), density_after);
    }
  }
}

// Assign mtf_acceleration to the particle
template <unsigned Tdim>
bool mpm::Particle<Tdim>::assign_mtf_acceleration(
    const Eigen::Matrix<double, Tdim, 1>& mtf_acceleration) {
  // Assign acceleration
  mtf_acceleration_ = mtf_acceleration;
  return true;
}