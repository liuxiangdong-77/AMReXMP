//! Read material properties
template <unsigned Tdim>
mpm::LinearElasticMTF<Tdim>::LinearElasticMTF(unsigned id,
                                              const Json& material_properties)
    : LinearElastic<Tdim>(id, material_properties) {
  // 多传一个history_size_的参数
  if (material_properties.find("history_size") != material_properties.end()) {
    history_size_ = material_properties["history_size"].get<unsigned>();//这里得在json里写
    console_->info("Material {}: history_size set to {}", id, history_size_);
  } else {
    history_size_ = 1;
    console_->warn("Material {}: history_size not provided, defaulting to {}", id, history_size_);
  }
}

//! Initialise state variables
//！初始化历史变量
template <unsigned Tdim>
mpm::dense_map mpm::LinearElasticMTF<Tdim>::initialise_state_variables() {
  
  mpm::dense_map state_vars;
  // Initialize history storage for each time step
  // 为每一时间步初始化历史变量容器
  for (unsigned t = 0; t < history_size_; ++t) {
    // // Density history
    // state_vars["density_t" + std::to_string(t)] = 0.0;//state_vars["density_t0"]=0.0;
    
    // Velocity history (initialize to zero for all components)
    for (unsigned i = 0; i < Tdim; ++i) {
      state_vars["velocity_" + std::to_string(i) + "_t" + std::to_string(t)] = 0.0;
    }//state_vars["velocity_0_t0"]=0.0;
    
    // // Stress history (initialize to zero for all components)
    // for (unsigned i = 0; i < 6; ++i) {
    //   state_vars["stress_" + std::to_string(i) + "_t" + std::to_string(t)] = 0.0;
    // }
  }
  return state_vars;
}

//! State variables
//! 存储历史变量
template <unsigned Tdim>
std::vector<std::string> mpm::LinearElasticMTF<Tdim>::state_variables() const {
  std::vector<std::string> state_vars;
  
  // Generate variable names for all history steps
  for (unsigned t = 0; t < history_size_; ++t) {
    // // Density
    // state_vars.push_back("density_t" + std::to_string(t));
    
    // Velocity components
    for (unsigned i = 0; i < Tdim; ++i) {
      state_vars.push_back("velocity_" + std::to_string(i) + "_t" + std::to_string(t));
    }
    
    // // Stress components
    // for (unsigned i = 0; i < 6; ++i) {
    //   state_vars.push_back("stress_" + std::to_string(i) + "_t" + std::to_string(t));
    // }
  }
  return state_vars;
}

//! Advance history variables
//在每个时间步结束时被调用，将历史数据向后移动，并存储当前值
template <unsigned Tdim>
void mpm::LinearElasticMTF<Tdim>::advance_history_variables(
    mpm::dense_map* state_vars, const ParticleBase<Tdim>* particle,
    unsigned long long time_step) {
  
  if (state_vars == nullptr || particle == nullptr) {
    console_->warn("Invalid input for advance_history_variables\n");
    return;
  }
  // console_->info("Advancing history variables for particle {} at time step = {}, history size = {}",
  //                particle->id(), time_step, history_size_);
  try {
    // 获取当前物理量
    // double current_density = particle->mass_density();
    auto current_velocity = particle->velocity();
    // auto current_stress = particle->stress();
    
    // // 验证当前值的有效性
    // if (std::isnan(current_density) || current_density <= 0) {
    //   console_->warn("Invalid current density: {}", current_density);
    //   current_density = 2000.0; // 使用默认密度
    // }
    
    // 统一时序推进：每步先右移历史，再写入当前到 t0
    for (unsigned t = history_size_ - 1; t > 0; --t) {
      for (unsigned i = 0; i < Tdim; ++i) {
        std::string current_key =
            "velocity_" + std::to_string(i) + "_t" + std::to_string(t);
        std::string previous_key =
            "velocity_" + std::to_string(i) + "_t" + std::to_string(t - 1);

        auto prev_it = state_vars->find(previous_key);
        if (prev_it != state_vars->end()) {
          (*state_vars)[current_key] = prev_it->second;
        }
        // 若 previous_key 不存在（理论上不会发生，initialise_state_variables
        // 已保证所有槽被初始化），保留该槽原值，不覆盖
        // else {
        //   (*state_vars)[current_key] = current_velocity(i);
        // }
      }
    }

    // // 平滑启动：在历史未填满阶段，把未初始化段补成当前值（避免启动尖峰）
    // if (time_step < history_size_ - 1) {
    //   const unsigned max_valid_history = static_cast<unsigned>(time_step + 1);
    //   for (unsigned t = max_valid_history; t < history_size_; ++t) {
    //     for (unsigned i = 0; i < Tdim; ++i) {
    //       (*state_vars)["velocity_" + std::to_string(i) + "_t" +
    //                     std::to_string(t)] = current_velocity(i);
    //     }
    //   }
    // }

    for (unsigned i = 0; i < Tdim; ++i) {
      (*state_vars)["velocity_" + std::to_string(i) + "_t0"] = current_velocity(i);
    }
    
    // ============= 调试输出 =============
    if (console_->level() <= spdlog::level::debug) {
      console_->debug("Step {}: History values after advancement:", time_step);
      for (unsigned t = 0; t < history_size_; ++t) {
        // std::string density_key = "density_t" + std::to_string(t);
        // auto it = state_vars->find(density_key);
        // if (it != state_vars->end()) {
        //   console_->debug("  {}: density={}", density_key, it->second);
        // }
        
        // 只检查第一个速度分量作为示例
        std::string vel_key = "velocity_0_t" + std::to_string(t);
        auto it = state_vars->find(vel_key);
        if (it != state_vars->end()) {
          console_->debug("  {}: velocity_0={}", vel_key, it->second);
        }
        
        // // 只检查第一个应力分量作为示例
        // std::string stress_key = "stress_0_t" + std::to_string(t);
        // it = state_vars->find(stress_key);
        // if (it != state_vars->end()) {
        //   console_->debug("  {}: stress_xx={}", stress_key, it->second);
        // }
      }
    }
    
  } catch (std::exception& exception) {
    console_->error("Failed to advance history variables at step {}: {}", 
                   time_step, exception.what());
  }
}

//! Get history variable at specific time step
template <unsigned Tdim>
double mpm::LinearElasticMTF<Tdim>::get_history_variable(
    const mpm::dense_map* state_vars, const std::string& var_name, 
    unsigned time_step) const {
  
  // 虚函数的具体实现
  if (state_vars == nullptr) {
    console_->warn("State variables pointer is null\n");
    return std::numeric_limits<double>::quiet_NaN();
  }
  
  if (time_step >= history_size_) {
    console_->warn("Requested time step {} exceeds history size {}\n", 
                   time_step, history_size_);
    return std::numeric_limits<double>::quiet_NaN();
  }
  
  std::string full_var_name = var_name + "_t" + std::to_string(time_step);

  
  if (state_vars->find(full_var_name) != state_vars->end()) {
    return state_vars->at(full_var_name);
  } else {
    console_->warn("History variable not found: {}\n", full_var_name);
    return std::numeric_limits<double>::quiet_NaN();
  }
}
