/*================================MTF=====================================*/
template <unsigned Tdim>
bool mpm::Mesh<Tdim>::create_nodal_mtf_constraint(
    int set_id,
    const std::shared_ptr<mpm::MTFConstraint>& constraint, unsigned phase, double dt) {
  bool status = true;
  try {
    bool constraint_exists = false;

    if (!mtf_console_) {
        mtf_console_ = spdlog::get("MTF");
        if (!mtf_console_) mtf_console_ = spdlog::stdout_color_mt("MTF");
    }

    for (const auto& existing_constraint : mtf_constraints_) {
      if (existing_constraint->setid() == set_id) {
        constraint_exists = true;
        mtf_console_->debug("节点集{}的MTF约束已存在,跳过重复创建", set_id);
        break;
      }
    }

    if (!constraint_exists && (set_id >= 0 || node_sets_.find(set_id) != node_sets_.end())) {
      // 第一次创建：添加到向量
      mtf_constraints_.emplace_back(constraint);
      mtf_console_->info("创建节点集{}的MTF约束,当前共{}个约束", 
                         set_id, mtf_constraints_.size());
    }

    // this->apply_mtf_constraints(phase, dt);

  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

//! Store particle history for MTF
template <unsigned Tdim>
void mpm::Mesh<Tdim>::store_particle_history_for_mtf(unsigned phase) {
  try {
    // Iterate over all particles to store history
    #pragma omp parallel for schedule(runtime)
    for (auto pitr = particles_.cbegin(); pitr != particles_.cend(); ++pitr) {
      (*pitr)->store_mtf_history(phase);
    }
    
    if (mtf_console_) {
      mtf_console_->debug("Stored MTF history for {} particles", particles_.size());
    }
  } catch (std::exception& exception) {
    if (mtf_console_) {
      mtf_console_->error("Failed to store MTF history: {}", exception.what());
    }
  }
}

//! Apply MTF constraints - 这是核心函数，移植自 mtf_constraint.tcc
template <unsigned Tdim>
void mpm::Mesh<Tdim>::apply_mtf_constraints(unsigned phase, double dt) {
  try {

    if (mtf_constraints_.empty()) {
      mtf_console_->debug("No MTF constraints to apply");//报错点lxd
      return;
    }

    mtf_console_->debug("Applying {} MTF constraints", mtf_constraints_.size());//size测试成功，Applying 1 MTF constraints
    
    // 对每个MTF约束集
    for (size_t i = 0; i < mtf_constraints_.size(); ++i) {
      auto mtf_constraint = mtf_constraints_[i];//成功
      auto set_id = mtf_constraint->setid();//成功
      
      // 获取该节点集的所有节点
      auto& nodes_in_set = (set_id == -1) ? nodes_ : node_sets_.at(set_id);//成功
      
      mtf_console_->debug("Applying MTF to nset {}: {} nodes, pos={}, order={}",
                         set_id, nodes_in_set.size(),
                         static_cast<int>(mtf_constraint->position()),
                         mtf_constraint->order());//测试成功（Applying MTF to nset 0: 1 nodes, dir=0, pos=10, order=3）
      
      // 对节点集中的每个节点应用MTF
      for (size_t j = 0; j < nodes_in_set.size(); ++j) {
        auto node = nodes_in_set[j];
        bool status = this->apply_mtf_to_node(
            node, phase, dt,
            mtf_constraint->position(),
            mtf_constraint->order(),
            mtf_constraint->drift_coeff()
        );
        
        if (!status) {
          mtf_console_->warn("Failed to apply MTF to node {} in nset {}", 
                            node->id(), set_id);
        }
      }
    }
  } catch (std::exception& exception) {
    if (mtf_console_) {
      mtf_console_->error("Failed to apply MTF constraints: {}", exception.what());
    }
  }
}

//! Apply MTF absorbing boundary as external force (Neumann form)
//! 特征分解型吸收边界：将入射特征速度 v_in 代入 LK 阻尼公式，作为结点外力施加
//! F_absorb = -rho*c*A*(v - v_in) = -0.5*(rho*c*A*v + f_int)
//! 与 Dirichlet 型（直接覆盖速度）相比，此版本保留结点 F=ma 的动量一致性。
template <unsigned Tdim>
void mpm::Mesh<Tdim>::apply_mtf_absorbing_forces(unsigned phase) {
  try {
    if (mtf_constraints_.empty()) {
      if (mtf_console_) mtf_console_->debug("No MTF absorbing constraints to apply");
      return;
    }

    // 对每个MTF约束集
    for (size_t i = 0; i < mtf_constraints_.size(); ++i) {
      auto mtf_constraint = mtf_constraints_[i];
      auto set_id = mtf_constraint->setid();

      // 获取该节点集的所有节点
      auto& nodes_in_set = (set_id == -1) ? nodes_ : node_sets_.at(set_id);

      // 对节点集中的每个节点施加吸收外力
      for (size_t j = 0; j < nodes_in_set.size(); ++j) {
        auto node = nodes_in_set[j];
        bool status = this->apply_mtf_absorbing_force(
            node, phase, mtf_constraint->position());

        if (!status && mtf_console_) {
          mtf_console_->warn("Failed to apply MTF absorbing force to node {} in nset {}",
                            node->id(), set_id);
        }
      }
    }
  } catch (std::exception& exception) {
    if (mtf_console_) {
      mtf_console_->error("Failed to apply MTF absorbing forces: {}", exception.what());
    }
  }
}

//! Apply MTF absorbing force to a single node
template <unsigned Tdim>
bool mpm::Mesh<Tdim>::apply_mtf_absorbing_force(
    const std::shared_ptr<mpm::NodeBase<Tdim>>& node,
    unsigned phase, mpm::Edge_Position boundary_position) {

  bool status = false;
  try {
    if (!node) {
      if (mtf_console_) mtf_console_->error("apply_mtf_absorbing_force: Node pointer is null");
      return false;
    }
    if (boundary_position == mpm::Edge_Position::None) return true;
    if (node->mass(phase) < 1e-15) return true;
    if (node->material_ids().empty()) return true;

    // 关键1：拆分角点为单一法向量
    std::vector<VectorDim> single_normals =
        split_corner_boundary_normals(boundary_position);
    if (single_normals.empty()) return true;

    // 关键2：P/S波速与密度
    VectorDim wave_speeds = node->compute_instantaneous_wave_speed(phase);
    double cp = wave_speeds(0);
    double cs = wave_speeds(1);
    if (cp <= 1e-6 || cs <= 1e-6) return true;

    double density = node->compute_mtf_density(phase);
    if (density <= 1e-15) return true;

    // 关键3：边界控制面积 A
    double h_min = this->compute_average_cell_size();
    double area = 0.0;
    if constexpr (Tdim == 2) {
      area = (single_normals.size() > 1) ? 0.5 * h_min : h_min;
    } else if constexpr (Tdim == 3) {
      if (single_normals.size() == 3)      area = 0.25 * h_min * h_min;
      else if (single_normals.size() == 2) area = 0.5 * h_min * h_min;
      else                                 area = h_min * h_min;
    }

    // 关键4：当前结点速度与内力
    VectorDim v = node->velocity(phase);
    VectorDim f_int = node->internal_force(phase);

    // 关键5：对每个法向计算吸收外力 F = -0.5*(rho*c*A*v + f_int)
    const double epsilon = 1e-10;
    VectorDim total_force = VectorDim::Zero();

    for (const auto& n : single_normals) {
      int normal_axis = -1;
      int tangent_axis = -1;
      if constexpr (Tdim == 2) {
        if (std::abs(n(0)) > epsilon && std::abs(n(1)) < epsilon) {
          normal_axis = 0;
          tangent_axis = 1;
        } else if (std::abs(n(1)) > epsilon && std::abs(n(0)) < epsilon) {
          normal_axis = 1;
          tangent_axis = 0;
        }
      }
      if (normal_axis < 0) {
        if (mtf_console_) mtf_console_->warn("Node{}: Abnormal normal, skip absorbing force", node->id());
        continue;
      }

      // P波（法向）：F_n = -0.5*(rho*cp*A*v_n + f_int_n)
      double F_n = -0.5 * (density * cp * area * v(normal_axis) + f_int(normal_axis));
      // S波（切向）：F_t = -0.5*(rho*cs*A*v_t + f_int_t)
      double F_t = -0.5 * (density * cs * area * v(tangent_axis) + f_int(tangent_axis));

      VectorDim force = VectorDim::Zero();
      force(normal_axis) = F_n;
      force(tangent_axis) = F_t;
      total_force += force;
    }

    // 关键6：角点统一做平均
    if (single_normals.size() > 1) total_force *= 0.5;

    // 关键7：作为结点外力施加（Neumann，参与 F=ma）
    node->update_external_force(true, phase, total_force);

    if (mtf_console_) {
      mtf_console_->debug("Node{}: MTF absorbing force norm = {:.4e}",
                          node->id(), total_force.norm());
    }

    status = true;
  } catch (std::exception& exception) {
    if (mtf_console_) {
      mtf_console_->error("apply_mtf_absorbing_force failed: {}", exception.what());
    }
    status = false;
  }
  return status;
}

// //! Apply MTF to a single node
// //！MTF主函数
// template <unsigned Tdim>
// bool mpm::Mesh<Tdim>::apply_mtf_to_node(
//     const std::shared_ptr<mpm::NodeBase<Tdim>>& node,
//     unsigned phase, double dt, mpm::Edge_Position boundary_position,
//     unsigned int order) {
//   bool status = false;
//   try {
//     /*==================================test========================================*/
//     if (!node) {
//       mtf_console_->error("apply_mtf_to_node: Node pointer is null");
//       return false;
//     }
//     // 检查mesh对象是否有效
//     if (!this) {
//       mtf_console_->error("apply_mtf_to_node: Mesh pointer is invalid");
//       return false;
//     }
//     // 检查容器是否有效
//     if (mtf_constraints_.size() == 0) {
//       mtf_console_->warn("No MTF constraints to apply");
//       return true;
//     }
//     // Only apply to boundary nodes with specific positions
//     if (boundary_position == mpm::Edge_Position::None) {
//       return true; // Not a boundary node, skip MTF
//     }
//     mtf_console_->debug("Node{}: Applying MTF constraint at {} with n={} order",
//                     node->id(), 
//                     static_cast<int>(boundary_position),
//                     order);//测试成功（Node65: Applying MTF constraint at 10 with n=3 order）
//     const double mass_tolerance = 1e-15;                
//     if (node->mass(phase) < mass_tolerance) {
//       mtf_console_->warn("Node{}: Skipping MTF - zero mass at phase {}", 
//                         node->id(), phase);//报错
//       return true;
//     }
//     if (node->material_ids().empty()) {
//       mtf_console_->warn("Node{}: Skipping MTF - no material assigned", 
//                         node->id());
//       return true; // 不视为错误，只是跳过 //测试成功
//     }
//     /*=================================test==========================================*/
//     // 关键1：拆分角点为单一法向量（右上结点→右+上，2个法向量）
//       std::vector<VectorDim> single_normals = split_corner_boundary_normals(boundary_position);//成功
//       if (single_normals.empty()) {
//         mtf_console_->warn("Node{}: No valid boundary normals, skip MTF", node->id());
//         return true;
//       }
//       // 关键2：获取P波/S波波速（cp=法向用，cs=切向用）
//       VectorDim wave_speeds = node->compute_instantaneous_wave_speed(phase);//成功
//       double cp = wave_speeds(0);//成功
//       double cs = wave_speeds(1);//成功
//       if (cp <= 1e-6 || cs <= 1e-6) {
//         mtf_console_->warn("Node{}: Invalid P/S wave speed (cp={}, cs={})", node->id(), cp, cs);
//         return true;
//       }
//       // 关键3：累加所有单一法向量对应的4套外推点结果（2个法向量×2套/法向量）
//       VectorDim total_mtf_velocity = VectorDim::Zero();
//       for (const auto& single_normal : single_normals) {
//         // 每套单一法向量对应2套外推点：
//         // 1. P波（cp）沿该法向外推 → 直接求解该法向对应的法向速度
//         VectorDim velocity_p = compute_mtf_velocity_single_wave(
//             node, phase, dt, single_normal, cp, order, true);
//         // 2. S波（cs）沿该法向外推 → 直接求解该法向对应的切向速度
//         VectorDim velocity_s = compute_mtf_velocity_single_wave(
//             node, phase, dt, single_normal, cs, order, false);
//         // 累加该法向量对应的P/S波速度
//         total_mtf_velocity += velocity_s + velocity_p;
//       }
//       if(single_normals.size() > 1){
//         total_mtf_velocity *= 0.5;
//       }
//       // 关键4改：更新结点MTF速度
//       node->update_mtf_velocity(total_mtf_velocity, phase);
//       // 供测试用，2D结点较多，先注释掉3.17
//       // VectorDim node_mtf_velo = node->mtf_velocity(phase);
//       // if (node_mtf_velo.norm() < 1e-15) {
//       //     mtf_console_->error("节点{}: mtf_velocity仍然是0!", node->id());
//       // } else {
//       //     mtf_console_->info("节点{}: mtf_velocity成功设置=[{}, {}]", 
//       //                       node->id(), node_mtf_velo(0), node_mtf_velo(1));//成功
//       // }
//       status = true;
//     } catch (std::exception& exception) {
//       mtf_console_->error("apply_dwa_mtf_to_node failed: {}", exception.what());
//       status = false;
//     }
//     return status;
// }


//! Apply MTF to a single node
//! MTF主函数 — 方案C：MTF外推 + 特征分解校正
template <unsigned Tdim>
bool mpm::Mesh<Tdim>::apply_mtf_to_node(
    const std::shared_ptr<mpm::NodeBase<Tdim>>& node,
    unsigned phase, double dt, mpm::Edge_Position boundary_position,
    unsigned int order, double drift_coeff) {
  
  bool status = false;
  try {

    /*==================================test========================================*/
    if (!node) {
      mtf_console_->error("apply_mtf_to_node: Node pointer is null");
      return false;
    }
    if (!this) {
      mtf_console_->error("apply_mtf_to_node: Mesh pointer is invalid");
      return false;
    }
    if (mtf_constraints_.size() == 0) {
      mtf_console_->warn("No MTF constraints to apply");
      return true;
    }
    if (boundary_position == mpm::Edge_Position::None) {
      return true;
    }
    mtf_console_->debug("Node{}: Applying MTF constraint at {} with n={} order",
                    node->id(), 
                    static_cast<int>(boundary_position),
                    order);
    
    const double mass_tolerance = 1e-15;                
    if (node->mass(phase) < mass_tolerance) {
      mtf_console_->warn("Node{}: Skipping MTF - zero mass at phase {}", 
                        node->id(), phase);
      return true;
    }
    if (node->material_ids().empty()) {
      mtf_console_->warn("Node{}: Skipping MTF - no material assigned", 
                        node->id());
      return true;
    }
    /*=================================test==========================================*/

    // 关键1：拆分角点为单一法向量（右上结点→右+上，2个法向量）
    std::vector<VectorDim> single_normals = split_corner_boundary_normals(boundary_position);
    if (single_normals.empty()) {
      mtf_console_->warn("Node{}: No valid boundary normals, skip MTF", node->id());
      return true;
    }

    // 关键2：获取P波/S波波速（cp=法向用，cs=切向用）
    VectorDim wave_speeds = node->compute_instantaneous_wave_speed(phase);
    double cp = wave_speeds(0);
    double cs = wave_speeds(1);
    if (cp <= 1e-6 || cs <= 1e-6) {
      mtf_console_->warn("Node{}: Invalid P/S wave speed (cp={}, cs={})", node->id(), cp, cs);
      return true;
    }

    // 关键3（新增）：获取特征分解所需的物理量
    double density = node->compute_mtf_density(phase);
    if (density <= 1e-15) {
      mtf_console_->warn("Node{}: Invalid density for characteristic decomposition", node->id());
      return true;
    }

    // 边界控制面积 A（和 Lysmer absorbing_constraint 中一致）
    // 2D: Edge结点 A = h_min, Corner结点 A = 0.5*h_min
    // 3D: Face结点 A = h^2, Edge结点 A = 0.5*h^2, Corner结点 A = 0.25*h^2
    double h_min = this->compute_average_cell_size();
    double area = 0.0;
    if constexpr (Tdim == 2) {
      area = (single_normals.size() > 1) ? 0.5 * h_min : h_min;
    } else if constexpr (Tdim == 3) {
      if (single_normals.size() == 3)      area = 0.25 * h_min * h_min;
      else if (single_normals.size() == 2) area = 0.5 * h_min * h_min;
      else                                 area = h_min * h_min;
    }

    // 当前结点的全速度和内力（特征分解用）
    // 确认：f_int = -σ·∇N·V（带负号），traction = -f_int/A
    VectorDim v_total = node->velocity(phase);
    VectorDim f_int = node->internal_force(phase);

    // ====== 关键4：对每个法向做 MTF外推 + 延时v_in校正 ======
    VectorDim total_mtf_velocity = VectorDim::Zero();
    VectorDim total_force = VectorDim::Zero();

    for (const auto& single_normal : single_normals) {
      // --- 4a. MTF 外推全速度（和原来完全一样，不改） ---
      VectorDim velocity_p = compute_mtf_velocity_single_wave(
          node, phase, dt, single_normal, cp, order, true, drift_coeff);
      VectorDim velocity_s = compute_mtf_velocity_single_wave(
          node, phase, dt, single_normal, cs, order, false, drift_coeff);

      // --- 4b. 特征分解：计算当前步的内行波分量 v_in ---
      int normal_axis = -1;
      int tangent_axis = -1;
      const double epsilon = 1e-10;
      if constexpr (Tdim == 2) {
        if (std::abs(single_normal(0)) > epsilon && 
            std::abs(single_normal(1)) < epsilon) {
          normal_axis = 0;
          tangent_axis = 1;
        } else if (std::abs(single_normal(1)) > epsilon && 
                   std::abs(single_normal(0)) < epsilon) {
          normal_axis = 1;
          tangent_axis = 0;
        } else {
          mtf_console_->warn("Node{}: Abnormal normal in characteristic decomposition", node->id());
          total_mtf_velocity += velocity_s + velocity_p;
          continue;
        }
      }

      // P波内行波：法向分量
      double v_in_normal = 0.5 * (v_total(normal_axis)
                           - f_int(normal_axis) / (density * cp * area));
      // S波内行波：切向分量
      double v_in_tangent = 0.5 * (v_total(tangent_axis)
                            - f_int(tangent_axis) / (density * cs * area));


      //同步衰减
      v_in_normal  *= 1.0/(1.0 + drift_coeff);
      v_in_tangent *= 1.0/(1.0 + drift_coeff);

      // MTF外推 - v_in 校正，在*0.5之前对每个法向都减一次
      velocity_p(normal_axis)  -= v_in_normal;
      velocity_s(tangent_axis) -= v_in_tangent;

      // 累加
      total_mtf_velocity += velocity_s + velocity_p;

      // 该法向对应的吸收外力：t = -ρc·(v_total - v_in) = -ρc·v_out
      VectorDim force = VectorDim::Zero();
      force(normal_axis) =
          -density * cp * area * (v_total(normal_axis) - v_in_normal);
      force(tangent_axis) =
          -density * cs * area * (v_total(tangent_axis) - v_in_tangent);
      total_force += force;

      mtf_console_->debug(
          "Node{}: CharDecomp (n=[{},{}]): "
          "v_tot=[{:.4e},{:.4e}], f_int=[{:.4e},{:.4e}], "
          "rho={:.1f}, cp={:.1f}, cs={:.1f}, A={:.4e}, ",
          node->id(), single_normal(0), single_normal(1),
          v_total(0), v_total(1), f_int(0), f_int(1),
          density, cp, cs, area);
    }

    // 关键5：角点统一做平均（和原来一样，但现在v_in校正也包含在内）
    if (single_normals.size() > 1) {
      total_mtf_velocity *= 0.5;
      total_force *= 0.5;
    }

    // 关键7：更新结点MTF吸收外力（Neumann，参与 F=ma）
    node->update_mtf_velocity(total_mtf_velocity, phase);

    status = true;
  } catch (std::exception& exception) {
    mtf_console_->error("apply_mtf_to_node failed: {}", exception.what());
    status = false;
  }
  return status;
}

// //! Apply MTF to a single node
// //! MTF主函数 — 动态波速 + 方案C特征分解校正
// template <unsigned Tdim>
// bool mpm::Mesh<Tdim>::apply_mtf_to_node(
//     const std::shared_ptr<mpm::NodeBase<Tdim>>& node,
//     unsigned phase, double dt, mpm::Edge_Position boundary_position,
//     unsigned int order, double drift_coeff) {
//   bool status = false;
//   try {
//     /*==================================检查========================================*/
//     if (!node) {
//       mtf_console_->error("apply_mtf_to_node: Node pointer is null");
//       return false;
//     }
//     if (!this) {
//       mtf_console_->error("apply_mtf_to_node: Mesh pointer is invalid");
//       return false;
//     }
//     if (mtf_constraints_.size() == 0) {
//       mtf_console_->warn("No MTF constraints to apply");
//       return true;
//     }
//     if (boundary_position == mpm::Edge_Position::None) {
//       return true;
//     }
//     mtf_console_->debug("Node{}: Applying MTF constraint at {} with n={} order",
//                     node->id(),
//                     static_cast<int>(boundary_position),
//                     order);
//     const double mass_tolerance = 1e-15;
//     if (node->mass(phase) < mass_tolerance) {
//       mtf_console_->warn("Node{}: Skipping MTF - zero mass at phase {}",
//                         node->id(), phase);
//       return true;
//     }
//     if (node->material_ids().empty()) {
//       mtf_console_->warn("Node{}: Skipping MTF - no material assigned",
//                         node->id());
//       return true;
//     }
//     /*==================================检查========================================*/
//     // 关键1：拆分角点为单一法向量
//     std::vector<VectorDim> single_normals = split_corner_boundary_normals(boundary_position);
//     if (single_normals.empty()) {
//       mtf_console_->warn("Node{}: No valid boundary normals, skip MTF", node->id());
//       return true;
//     }
//     // 关键2：获取P波/S波波速
//     VectorDim wave_speeds = node->compute_instantaneous_wave_speed(phase);
//     double cp = wave_speeds(0);
//     double cs = wave_speeds(1);
//     if (cp <= 1e-6 || cs <= 1e-6) {
//       mtf_console_->warn("Node{}: Invalid P/S wave speed (cp={}, cs={})", node->id(), cp, cs);
//       return true;
//     }
//     // 关键3：获取特征分解所需的物理量
//     double density = node->compute_mtf_density(phase);
//     if (density <= 1e-15) {
//       mtf_console_->warn("Node{}: Invalid density", node->id());
//       return true;
//     }
//     // 边界控制面积
//     double h_min = this->compute_average_cell_size();
//     double area = 0.0;
//     if constexpr (Tdim == 2) {
//       area = (single_normals.size() > 1) ? 0.5 * h_min : h_min;
//     } else if constexpr (Tdim == 3) {
//       if (single_normals.size() == 3)      area = 0.25 * h_min * h_min;
//       else if (single_normals.size() == 2) area = 0.5 * h_min * h_min;
//       else                                 area = h_min * h_min;
//     }
//     // 当前结点的全速度和内力
//     VectorDim v_total = node->velocity(phase);
//     VectorDim f_int = node->internal_force(phase);
//     // ====== 关键4：对每个法向做 动态波速估计 + MTF外推 + 特征分解 ======
//     VectorDim total_mtf_velocity = VectorDim::Zero();
//     for (const auto& single_normal : single_normals) {
//       // --- 4a. 判断法向/切向坐标轴 ---
//       int normal_axis = -1;
//       int tangent_axis = -1;
//       const double epsilon = 1e-10;
//       if constexpr (Tdim == 2) {
//         if (std::abs(single_normal(0)) > epsilon &&
//             std::abs(single_normal(1)) < epsilon) {
//           normal_axis = 0;
//           tangent_axis = 1;
//         } else if (std::abs(single_normal(1)) > epsilon &&
//                    std::abs(single_normal(0)) < epsilon) {
//           normal_axis = 1;
//           tangent_axis = 0;
//         } else {
//           mtf_console_->warn("Node{}: Abnormal normal", node->id());
//           continue;
//         }
//       }
//       // --- 4b. 动态波速估计（基于波场回溯匹配） ---
//       // 回溯步数：可调参数，建议 5~20
//       int n_back = 10;
//       // P波视波速估计
//       double theta_p = estimate_transmitting_angle(
//           node, phase, dt, single_normal,
//           normal_axis, tangent_axis, cp, n_back);
//       double cos_theta_p = std::cos(theta_p);
//       cos_theta_p = std::max(std::abs(cos_theta_p), 0.3);  // 限幅防掠入射
//       double cp_app = cp / cos_theta_p;  // P波视波
//       // S波视波速估计
//       double theta_s = estimate_transmitting_angle(
//           node, phase, dt, single_normal,
//           normal_axis, tangent_axis, cs, n_back);
//       double cos_theta_s = std::cos(theta_s);
//       cos_theta_s = std::max(std::abs(cos_theta_s), 0.3);
//       double cs_app = cs / cos_theta_s;  // S波视波速
//       mtf_console_->debug(
//           "Node{}: Dynamic wave speed: theta_p={:.1f}deg cp_app={:.1f}, "
//           "theta_s={:.1f}deg cs_app={:.1f}",
//           node->id(),
//           theta_p * 180.0 / M_PI, cp_app,
//           theta_s * 180.0 / M_PI, cs_app);
//       // --- 4c. MTF 外推（用动态视波速替代固定cp/cs） ---
//       VectorDim velocity_p = compute_mtf_velocity_single_wave(
//           node, phase, dt, single_normal, cp_app, order, true, drift_coeff);
//       VectorDim velocity_s = compute_mtf_velocity_single_wave(
//           node, phase, dt, single_normal, cs_app, order, false, drift_coeff);
//       // --- 4d. 特征分解校正（用动态视波速） ---
//       // v_in = 0.5 * (v - f_int / (rho * c_app * A))
//       double v_in_normal = 0.5 * (v_total(normal_axis)
//                            - f_int(normal_axis) / (density * cp_app * area));
//       double v_in_tangent = 0.5 * (v_total(tangent_axis)
//                             - f_int(tangent_axis) / (density * cs_app * area));
//       // 同步衰减（和MTF消飘系数一致）
//       if (drift_coeff > 0.0) {
//         double drift_scale = 1.0 / (1.0 + drift_coeff);
//         v_in_normal  *= drift_scale;
//         v_in_tangent *= drift_scale;
//       }
//       // MTF外推 - v_in校正
//       // velocity_p(normal_axis)  -= v_in_normal;
//       // velocity_s(tangent_axis) -= v_in_tangent;
//       // 累加
//       total_mtf_velocity += velocity_s + velocity_p;
//     }
//     // 关键5：角点统一做平均
//     if (single_normals.size() > 1) {
//       total_mtf_velocity *= 0.5;
//     }
//     // 关键6：更新结点MTF速度
//     node->update_mtf_velocity(total_mtf_velocity, phase);
//     status = true;
//   } catch (std::exception& exception) {
//     mtf_console_->error("apply_mtf_to_node failed: {}", exception.what());
//     status = false;
//   }
//   return status;
// }


//! Estimate transmitting angle by wave-field image matching
//  基于波场回溯匹配的透射角动态识别（简化版李小军2024方法）
//
//  原理：
//    当前时刻 t，边界结点 B 的速度为 v_B(t)
//    如果波沿角度 θ 传播，那么在 t - n_back*dt 时刻，
//    这个波位于从 B 出发沿 θ 方向往域内 R = c*n_back*dt 距离处
//    在那个位置插值历史速度 v_ref，和 v_B(t) 比较
//    误差最小的 θ 就是透射角
//
template <unsigned Tdim>
double mpm::Mesh<Tdim>::estimate_transmitting_angle(
    const std::shared_ptr<mpm::NodeBase<Tdim>>& node,
    unsigned phase, double dt,
    const VectorDim& single_normal,
    int normal_axis, int tangent_axis,
    double wave_speed, int n_back) {

  double best_theta = 0.0;  // 默认垂直入射

  try {
    VectorDim v_current = node->velocity(phase);

    // 波场方向性不显著时（速度太小），直接返回0°
    if (v_current.norm() < 1e-12) {
      return 0.0;
    }

    VectorDim node_coord = node->coordinates();
    double R = wave_speed * n_back * dt;  // 回溯距离

    // 防止回溯距离太小导致识别不准
    if (R < 1e-10) {
      return 0.0;
    }

    double min_error = std::numeric_limits<double>::max();

    // 法向方向的符号：single_normal 指向域外
    // 往域内回溯 = 减去 single_normal 方向
    // 法向分量：-single_normal(normal_axis) * R * cos(theta)
    // 切向偏移：R * sin(theta)

    // 角度扫描：-75° ~ +75°，间隔5°
    // 不扫到 ±90° 因为掠入射时视波速趋于无穷大，没有实际意义
    for (int deg = -75; deg <= 75; deg += 5) {
      double theta_rad = deg * M_PI / 180.0;

      // 从边界结点出发，沿 theta 方向往域内回溯 R 距离
      // 法向（指向域内）= -single_normal
      // 切向偏移方向：沿 tangent_axis
      VectorDim ref_point = node_coord;
      ref_point(normal_axis) -= single_normal(normal_axis) * R * std::cos(theta_rad);
      ref_point(tangent_axis) += R * std::sin(theta_rad);

      // 插值 t - n_back*dt 时刻在参考点处的速度
      VectorDim v_ref = VectorDim::Zero();
      try {
        v_ref = this->interpolate_particle_mtf_velo_at_point(
            ref_point, phase, -n_back);
      } catch (...) {
        // 参考点在域外或无粒子，跳过这个角度
        continue;
      }

      // 如果插值结果是零（参考点处没有粒子），跳过
      if (v_ref.norm() < 1e-15) {
        continue;
      }

      // 计算误差：当前边界速度 vs 参考点历史速度
      double error = (v_current - v_ref).squaredNorm();

      if (error < min_error) {
        min_error = error;
        best_theta = theta_rad;
      }
    }

    // 安全检查：如果所有角度都失败了，返回0
    if (min_error >= std::numeric_limits<double>::max() * 0.5) {
      best_theta = 0.0;
    }

    mtf_console_->debug(
        "Node{}: Angle scan result: best_theta={:.1f}deg, min_error={:.4e}",
        node->id(), best_theta * 180.0 / M_PI, min_error);

  } catch (std::exception& exception) {
    mtf_console_->error("estimate_transmitting_angle failed: {}", exception.what());
    best_theta = 0.0;
  }

  return best_theta;
}



//MTF公式
template <unsigned Tdim>
Eigen::Matrix<double, Tdim, 1> mpm::Mesh<Tdim>::compute_mtf_velocity_single_wave(
    const std::shared_ptr<NodeBase<Tdim>>& node,
    unsigned phase,
    double dt,
    const VectorDim& single_normal,
    double wave_speed,
    unsigned int order,
    bool is_p_wave,
    double drift_coeff) {
  VectorDim mtf_velocity = VectorDim::Zero();
  try {
    // 关键：外推方向固定为单一法向量（朝向计算域内，不单位化，无切向量）
    VectorDim extrapolate_direction = single_normal;//成功

    // 核心步骤1：判断法向对应的坐标轴（x轴=0，y轴=1）
    int normal_axis = -1; // 法向坐标轴索引（0=x，1=y）
    int tangent_axis = -1; // 切向坐标轴索引（0=x，1=y）
    const double epsilon = 1e-10;
    if constexpr (Tdim == 2) {
      // 判断法向量是否沿x轴（Left/Right边界）
      if (std::abs(single_normal(0)) > epsilon && std::abs(single_normal(1)) < epsilon) {
        normal_axis = 0;
        tangent_axis = 1;
      }
      // 判断法向量是否沿y轴（Bottom/Top边界）
      else if (std::abs(single_normal(1)) > epsilon && std::abs(single_normal(0)) < epsilon) {
        normal_axis = 1;
        tangent_axis = 0;
      }
      // 异常处理：非坐标轴方向的法向量（当前规则域不会出现）
      else {
        mtf_console_->warn("Node{}: Abnormal normal vector (not along x/y axis), skip wave calculation", node->id());
        return mtf_velocity;
      }
    }//成功。Bottom: 法向对应Y轴，normlal_axis = 1;

    // MTF核心公式：直接求解法向（P波）/切向（S波）加速度，无投影
    for (unsigned m = 1; m <= order; ++m) {
      // 1. 二项式系数（复用你的正确实现）
      double binomial_coeff = this->binomial_coefficient(order, m);//成功
      double coeff = std::pow(-1.0, m + 1) * binomial_coeff;//成功

      //零频漂移修正系数(可修改)
      double amend_coeff = 1 / std::pow(1 + drift_coeff, m);//可修改系数4.7

      // 2. 计算外推点坐标（沿法向，距离=m*wave_speed*dt，不单位化保证准确）
      double extrapolate_distance = m * wave_speed * dt;//成功
      VectorDim extrapoint_coord = node->coordinates() - (extrapolate_direction * extrapolate_distance);//成功

      // 3. 时间偏移（对应历史时间步，复用你的逻辑）
      int time_offset = 1 - static_cast<int>(m);

      // 4. 插值外推点速度
      VectorDim velo_extrapolated = this->interpolate_particle_mtf_velo_at_point(
        extrapoint_coord, phase, time_offset);

      // 5. 核心修改：直接提取坐标轴分量，按波类型累加（无误差，贴合规则域）
      if (is_p_wave) {
        // P波：仅累加【法向坐标轴】对应的速度分量
        mtf_velocity(normal_axis) += coeff * velo_extrapolated(normal_axis) * amend_coeff;
      } else {
        // S波：仅累加【切向坐标轴】对应的速度分量
        mtf_velocity(tangent_axis) += coeff * velo_extrapolated(tangent_axis) * amend_coeff;
      }
    }

    // 调试信息
    std::string wave_type = is_p_wave ? "P-wave (normal accel)" : "S-wave (tangent accel)";
    mtf_console_->debug("Node{}: {} (dir=[{},{}], speed={}) = [{},{}]",
                    node->id(), wave_type,
                    extrapolate_direction(0), extrapolate_direction(1),
                    wave_speed, mtf_velocity(0), mtf_velocity(1));
  } catch (std::exception& exception) {
    mtf_console_->error("compute_mtf_velocity_single_wave failed: {}", exception.what());
  }
  return mtf_velocity;
}

// template <unsigned Tdim>
// Eigen::Matrix<double, Tdim, 1> mpm::Mesh<Tdim>::compute_mtf_velocity_single_wave(
//     const std::shared_ptr<NodeBase<Tdim>>& node,
//     unsigned phase, double dt,
//     const VectorDim& single_normal,
//     double wave_speed,
//     unsigned int order,
//     bool is_p_wave) {
//   if (order != 1) {
//     mtf_console_->warn("先固定order=1测试,高阶以后再开");
//     return VectorDim::Zero();
//   }
//   VectorDim mtf_vel = VectorDim::Zero();
//   try {
//     // 法向轴判断（0=x, 1=y）
//     int normal_axis = (std::abs(single_normal(0)) > 1e-10) ? 0 : 1;
//     int tangent_axis = 1 - normal_axis;
//     // 真实网格间距ds（关键！不能只用wave_speed*dt）
//     double ds = this->compute_average_cell_size();   // 或沿法向真实网格尺寸
//     // α系数（解决cΔt ≠ ds的核心）
//     double alpha = (wave_speed * dt - ds) / (wave_speed * dt + ds);
//     // 当前边界节点历史速度 v0^n
//     VectorDim v0_n = node->velocity(phase);
//     // 内侧外推点当前步速度 v1^n
//     VectorDim extrapoint = node->coordinates() - (single_normal * ds);
//     VectorDim v1_n = this->interpolate_particle_mtf_velo_at_point(extrapoint, phase, 0);
//     // **正确一阶MTF**（带α）
//     if (is_p_wave) {
//       mtf_vel(normal_axis) = v1_n(normal_axis) + alpha * (v0_n(normal_axis) - v1_n(normal_axis));
//     } else {
//       mtf_vel(tangent_axis) = v1_n(tangent_axis) + alpha * (v0_n(tangent_axis) - v1_n(tangent_axis));
//     }
//     mtf_console_->debug("Node{} 一阶MTF ({} wave, alpha={:.3f}): vel=[{:.4e},{:.4e}]", 
//                         node->id(), is_p_wave?"P":"S", alpha, mtf_vel(0), mtf_vel(1));
//   } catch (...) {}
//   return mtf_vel;
// }

template <unsigned Tdim>
std::vector<Eigen::Matrix<double, Tdim, 1>> mpm::Mesh<Tdim>::split_corner_boundary_normals(mpm::Edge_Position boundary_position) {
  std::vector<VectorDim> single_normals;
  if constexpr (Tdim == 2) { // 仅适配二维，符合你的场景
    switch (boundary_position) {
      // 普通边界：1个单一法向量（对应2套外推点）
      case Edge_Position::Left:
        single_normals.emplace_back(VectorDim(-1.0, 0.0));
        break;
      case Edge_Position::Right:
        single_normals.emplace_back(VectorDim(1.0, 0.0));
        break;
      case Edge_Position::Bottom:
        single_normals.emplace_back(VectorDim(0.0, -1.0));
        break;
      case Edge_Position::Top:
        single_normals.emplace_back(VectorDim(0.0, 1.0));
        break;

      // 角点：2个单一法向量（对应4套外推点，2*2）
      case Edge_Position::TopRight: // 右上结点：右+上
        single_normals.emplace_back(VectorDim(1.0, 0.0));  // 右法向
        single_normals.emplace_back(VectorDim(0.0, 1.0));  // 上法向
        break;
      case Edge_Position::TopLeft: // 左上结点：左+上
        single_normals.emplace_back(VectorDim(-1.0, 0.0)); // 左法向
        single_normals.emplace_back(VectorDim(0.0, 1.0));  // 上法向
        break;
      case Edge_Position::BottomRight: // 右下结点：右+下
        single_normals.emplace_back(VectorDim(1.0, 0.0));  // 右法向
        single_normals.emplace_back(VectorDim(0.0, -1.0)); // 下法向
        break;
      case Edge_Position::BottomLeft: // 左下结点：左+下
        single_normals.emplace_back(VectorDim(-1.0, 0.0)); // 左法向
        single_normals.emplace_back(VectorDim(0.0, -1.0)); // 下法向
        break;

      // 其他位置：空
      default:
        break;
    }
  }
  return single_normals;
}


//获取粒子速度（当前时间步）
template <unsigned Tdim>
Eigen::Matrix<double, Tdim, 1> mpm::Mesh<Tdim>::interpolate_particle_mtf_velo_at_point(
    const VectorDim& point_coord,
    unsigned phase,
    int time_offset) {
  VectorDim interpolated_velo = VectorDim::Zero();
  try {
    std::shared_ptr<Cell<Tdim>> found_cell = nullptr;
    auto all_cells = this->cells(); //成功
    // 1. 定位包含外推点的cell
    #pragma omp parallel for schedule(runtime)
    for (size_t i = 0; i < all_cells.size(); ++i) {
      auto cell = all_cells[i];
      VectorDim xi;
      if (!found_cell && cell->is_mtfpoint_in_cell(point_coord, &xi)) {
        #pragma omp critical
        {
          if (!found_cell) found_cell = cell;
        }
      }
    }// 成功！！
    // 2. Fallback：找最近cell
    if (!found_cell) {
      // mtf_console_->warn("MTF point [{}, {}] not in any cell, use nearest", point_coord(0), point_coord(1));
      double min_distance = std::numeric_limits<double>::max();
      // 核心修改：用 all_cells.for_each() 替代范围for循环
      all_cells.for_each([&](std::shared_ptr<mpm::Cell<Tdim>> cell) {
        // 以下业务逻辑完全不变，直接复用
        if (!cell) return;  // 注意：for_each 中用 return 替代 continue（效果一致）
        double distance = (cell->centroid() - point_coord).norm();
        if (distance < min_distance) {
          min_distance = distance;
          found_cell = cell;
        }
      }); 
      if (!found_cell) throw std::runtime_error("No valid cell found for extrapolation point");
    }
    // 3. 获取cell+相邻cell的所有粒子ID
    std::vector<Index> all_particle_ids;
    auto main_particle_ids = found_cell->particles(); // 成功
    all_particle_ids.insert(all_particle_ids.end(), main_particle_ids.begin(), main_particle_ids.end());
  for (const auto& neighbour_cell_id : found_cell->neighbours()) {
    auto neighbour_cell = this->cell(neighbour_cell_id);
    if (neighbour_cell) {
      auto neighbour_particles = neighbour_cell->particles();
      all_particle_ids.insert(all_particle_ids.end(), neighbour_particles.begin(), neighbour_particles.end());
    }
  }
  if (all_particle_ids.empty()) throw std::runtime_error("No particles found for extrapolation");


  // 4. 直接收集粒子的mtf_velocity + 新增收集粒子坐标（MLS插值需要）
  const double support_radius = 2.0 * this->compute_average_cell_size();//可改


  const double epsilon = 1.E-15;
  // 【修改1】新增：存储有效粒子坐标（MLS插值需要粒子坐标来构建基函数）
  std::vector<VectorDim> valid_particle_coords;
  std::vector<VectorDim> valid_particle_velos;
  std::vector<double> valid_particle_weights;
  // 【修改2】给坐标向量预留空间（和其他两个向量保持一致）
  valid_particle_coords.reserve(all_particle_ids.size());
  valid_particle_velos.reserve(all_particle_ids.size());
  valid_particle_weights.reserve(all_particle_ids.size());
  
  for (auto particle_id : all_particle_ids) {
    auto particle_ptr = this->particle(particle_id);
    if (!particle_ptr) continue;
    VectorDim particle_coord = particle_ptr->coordinates();//成功
    double distance = (particle_coord - point_coord).norm();
    if (distance > support_radius) continue;
    double weight = this->compute_mls_weight(distance, support_radius);//成功 =0.5686
    if (weight < epsilon) continue;

    //1阶MTF
    // VectorDim particle_velo = particle_ptr->velocity();

    //高阶MTF
    VectorDim particle_velo = VectorDim::Zero();
    Eigen::Matrix<double, 6, 1> dummy_stress;
    double dummy_density;

    this->get_particle_quantities(particle_ptr, phase, time_offset,
                                  particle_velo, dummy_stress, dummy_density);

    if (particle_velo.hasNaN()) {
      mtf_console_->warn("Particle {}: Invalid mtf velocity, skip", particle_id);
      continue;
    }
    // 【修改3】新增：存入有效粒子坐标（MLS必需）
    valid_particle_coords.push_back(particle_coord);
    valid_particle_velos.push_back(particle_velo);
    valid_particle_weights.push_back(weight);
  }

  if (valid_particle_weights.empty()) throw std::runtime_error("No valid particles for MLS velocity interpolation");
  // 5. 【核心修改】替换简单加权平均为 MLS 插值速度
  interpolated_velo = this->mls_interpolate_velo_at_point(
      point_coord,
      valid_particle_coords,   // 有效粒子坐标
      valid_particle_velos,   // 有效粒子加速度
      valid_particle_weights   // 有效粒子MLS权重
  );
  } catch (std::exception& exception) {
    mtf_console_->error("interpolate_particle_mtf_velo_at_point failed: {}", exception.what());
  }
  return interpolated_velo;
  // double min_dist = std::numeric_limits<double>::max();
  // VectorDim nearest_accel = VectorDim::Zero();
  // for (auto particle_id : all_particle_ids) {
  //     auto particle_ptr = this->particle(particle_id);
  //     if (!particle_ptr) continue; 
  //     VectorDim particle_coord = particle_ptr->coordinates();
  //     double distance = (particle_coord - point_coord).norm();//成功 
  //     if (distance < min_dist) {
  //         min_dist = distance;
  //         nearest_accel = particle_ptr->mtf_acceleration();
  //     }
  // }
  // mtf_console_->info("最近粒子距离={}, 加速度=[{}, {}]",
  //                   min_dist, nearest_accel(0), nearest_accel(1));
  //return nearest_accel;
}

// 1D MLS
template <unsigned Tdim>
Eigen::Matrix<double, Tdim, 1> mpm::Mesh<Tdim>::mls_interpolate_velo_at_point(
    const VectorDim& point_coord,  // 外推点坐标（插值目标点）
    const std::vector<VectorDim>& valid_particle_coords,  // 搜索域内有效粒子坐标
    const std::vector<VectorDim>& valid_particle_velos,  // 搜索域内有效粒子加速度
    const std::vector<double>& valid_particle_weights) {  // 搜索域内有效粒子的MLS权重
  // 初始化插值结果（外推点的加速度）
  VectorDim interpolated_velo = VectorDim::Zero();
  try {
    // ==================== 常数基函数核心：仅含常数项[1]，无需矩阵运算 ====================
    const unsigned basis_size = 1;  // 论文定义：常数基函数维度=1（仅[1]）
    double total_weight = 0.0;       // 权重总和（替代矩矩阵，避免求逆）
    // ==================== 仅需计算“权重×加速度”的累加和 ====================
    for (size_t i = 0; i < valid_particle_weights.size(); ++i) {
      double weight = valid_particle_weights[i];
      VectorDim particle_velo = valid_particle_velos[i];
      // 1. 累加权重（常数基函数的矩矩阵为1×1的标量，即权重总和）
      total_weight += weight;
      // 2. 累加“权重×粒子加速度”（对应论文中b(x)的计算）
      //interpolated_displacement += weight * particle_displacement;
      interpolated_velo += weight * particle_velo;
    }
    // ==================== 重构外推点加速度：加权平均（论文§3.1常数基函数逻辑） ====================
    const double epsilon = 1e-15;
    if (total_weight < epsilon) {
      mtf_console_->warn("Total weight of constant basis is zero, return zero velocity");
      //return interpolated_displacement;
      return interpolated_velo;
    }
    // 常数基函数插值结果 = （权重×加速度总和） / 权重总和（本质是加权平均）
    interpolated_velo /= total_weight;
    // 调试信息：输出关键参数，验证常数基函数生效
    mtf_console_->debug("Constant basis MLS: Total weight={:.4e}, Interpolated velocity=[{:.4e}, {:.4e}]",
                    total_weight, interpolated_velo(0), interpolated_velo(1));
  } catch (std::exception& exception) {
    mtf_console_->error("Constant basis MLS interpolate failed: {}", exception.what());
  }
  //return interpolated_displacement;
  return interpolated_velo;
}

// // MLS核心算法
// template <unsigned Tdim>
// Eigen::Matrix<double, Tdim, 1> mpm::Mesh<Tdim>::mls_interpolate_velo_at_point(
//     const VectorDim& point_coord,  // 外推点坐标（插值目标点）
//     const std::vector<VectorDim>& valid_particle_coords,  // 搜索域内有效粒子坐标
//     const std::vector<VectorDim>& valid_particle_velos,  // 搜索域内有效粒子速度
//     const std::vector<double>& valid_particle_weights) {  // 搜索域内有效粒子的MLS权重
//   // 初始化插值结果（外推点的速度）
//   VectorDim interpolated_velo = VectorDim::Zero();
//   try {
//     // ==================== 复用你原有MLS核心：线性多项式基函数（无改动） ====================
//     const unsigned basis_size = Tdim + 1;  // 二维=3，三维=4，和你原有逻辑一致
//     Eigen::Matrix<double, basis_size, basis_size> moment_matrix = Eigen::Matrix<double, basis_size, basis_size>::Zero();//M(x)
//     Eigen::VectorXd rhs_velo = Eigen::VectorXd::Zero(basis_size * Tdim);  // 右端项适配速度（Tdim个分量）3*2=6个分量的列向量
//     // ==================== 复用你原有MLS核心：构建矩矩阵和右端项（仅修改右端项为加速度） ====================
//     for (size_t i = 0; i < valid_particle_weights.size(); ++i) {
//       double weight = valid_particle_weights[i];
//       VectorDim particle_coord = valid_particle_coords[i];  // 成功
//       VectorDim particle_velo = valid_particle_velos[i];  // 
//       // 1. 计算基函数 Φ(x)
//       Eigen::VectorXd basis_at_particle = Eigen::VectorXd::Zero(basis_size);
//       basis_at_particle(0) = 1.0;  // 常数项
//       for (unsigned dim = 0; dim < Tdim; ++dim) {
//         // 基函数的一次项：(x_i - x_target)
//         basis_at_particle(dim + 1) = particle_coord(dim) - point_coord(dim);//P(xi-x)
//       }
//       // 2. 构建矩矩阵 M = Σ(w_i * Φ_i * Φ_i^T)
//       moment_matrix += weight * basis_at_particle * basis_at_particle.transpose();//M(x)，矩阵值很小（病态矩阵！）找到问题！10-5量级
//       // 3. 构建速度右端项 b = Σ(w_i * Φ_i * a_i)（仅此处修改，替换为粒子加速度）
//       for (unsigned dim = 0; dim < Tdim; ++dim) {
//         // 每个速度分量对应一组右端项，和你原有多分量插值逻辑一致
//         rhs_velo.segment(dim * basis_size, basis_size) += 
//             weight * particle_velo(dim) * basis_at_particle;//b(x)
//       }
//     }
//     // ==================== 复用你原有MLS核心：求解矩矩阵+重构目标点加速度（无改动） ====================
//     // 1. 奇异性判断（兜底逻辑，复用你原有代码）
//     Eigen::FullPivLU<Eigen::Matrix<double, basis_size, basis_size>> lu_decomp(moment_matrix);
//     if (!lu_decomp.isInvertible()) {
//       mtf_console_->warn("MLS moment matrix is singular, fallback to weighted average for acceleration");
//       // 兜底：加权平均（和你原有兜底逻辑一致，仅适配速度）
//       return this->weighted_average_velo_fallback(valid_particle_velos, valid_particle_weights);
//     }
//     // 2. 求矩矩阵的逆（复用你原有逻辑）
//     Eigen::Matrix<double, basis_size, basis_size> moment_inv = lu_decomp.inverse();//值很大？
//     // 3. 计算目标点的基函数（完全复用你原有逻辑，无改动）
//     Eigen::VectorXd basis_at_point = Eigen::VectorXd::Zero(basis_size);
//     basis_at_point(0) = 1.0;
//     for (unsigned dim = 0; dim < Tdim; ++dim) {
//       basis_at_point(dim + 1) = point_coord(dim) - point_coord(dim);  // 目标点自身，该项为0，和你原有逻辑一致 //P(z-x)
//     }
//     // 4. 重构外推点速度
//     for (unsigned dim = 0; dim < Tdim; ++dim) {
//       // 复用你原有：v_target = Φ_target^T * M^{-1} * b_dim
//       Eigen::VectorXd coeffs_velo = moment_inv * rhs_velo.segment(dim * basis_size, basis_size);
//       interpolated_velo(dim) = basis_at_point.dot(coeffs_velo);//式(1)
//     }
//   } catch (std::exception& exception) {
//     mtf_console_->error("mls_interpolate_velo_at_point failed: {}", exception.what());
//   }
//   return interpolated_velo;
// }

// // MLS-2D QR
// template <unsigned Tdim>
// Eigen::Matrix<double, Tdim, 1> mpm::Mesh<Tdim>::mls_interpolate_velo_at_point(
//     const VectorDim& point_coord,  
//     const std::vector<VectorDim>& valid_particle_coords,  
//     const std::vector<VectorDim>& valid_particle_velos,  
//     const std::vector<double>& valid_particle_weights) {  
//   VectorDim interpolated_velo = VectorDim::Zero();
//   try {
//     const unsigned basis_size = Tdim + 1;  // 2D=3，3D=4，逻辑不变
//     Eigen::Matrix<double, basis_size, basis_size> moment_matrix = Eigen::Matrix<double, basis_size, basis_size>::Zero();
//     Eigen::VectorXd rhs_velo = Eigen::VectorXd::Zero(basis_size * Tdim);
//     const double support_radius = 2.0 * this->compute_average_cell_size();
//     // 1. 构建矩矩阵和右端项（核心逻辑不变）
//     for (size_t i = 0; i < valid_particle_weights.size(); ++i) {
//       double weight = valid_particle_weights[i];
//       if (std::fabs(weight) < 1e-15) continue;  // 新增：过滤极小权重，避免噪声
//       VectorDim particle_coord = valid_particle_coords[i];
//       VectorDim particle_velo = valid_particle_velos[i];
//       // 计算2D线性基：[1, x_i-x, y_i-y]
//       Eigen::VectorXd basis_at_particle = Eigen::VectorXd::Zero(basis_size);
//       basis_at_particle(0) = 1.0;  
//       for (unsigned dim = 0; dim < Tdim; ++dim) {
//         // basis_at_particle(dim + 1) = particle_coord(dim) - point_coord(dim);
//         // 归一化，防止数值问题
//         basis_at_particle(dim + 1) = (particle_coord(dim) - point_coord(dim)) / support_radius;
//       }
//       // 构建矩矩阵 M = Σw_i * Φ_i * Φ_i^T
//       moment_matrix += weight * basis_at_particle * basis_at_particle.transpose();
//       // 构建右端项 b = Σw_i * Φ_i * v_i
//       for (unsigned dim = 0; dim < Tdim; ++dim) {
//         rhs_velo.segment(dim * basis_size, basis_size) += 
//             weight * particle_velo(dim) * basis_at_particle;
//       }
//     }
//     // 2. 替换LU分解为QR分解（核心修改）
//     Eigen::ColPivHouseholderQR<Eigen::Matrix<double, basis_size, basis_size>> qr(moment_matrix);
//     // 检查矩阵秩（比LU的isInvertible()更鲁棒）
//     if (qr.rank() < basis_size) {
//       mtf_console_->warn("MLS moment matrix rank deficient (rank={}), fallback to weighted average", qr.rank());
//       return this->weighted_average_velo_fallback(valid_particle_velos, valid_particle_weights);
//     }
//     // 3. 目标点基函数（不变）
//     Eigen::VectorXd basis_at_point = Eigen::VectorXd::Zero(basis_size);
//     basis_at_point(0) = 1.0;
//     // 目标点自身的线性基项为0，无需循环赋值（等价于你的代码，更简洁）
//     for (unsigned dim = 0; dim < Tdim; ++dim) {
//       basis_at_point(dim + 1) = 0.0;
//     }
//     // 4. 重构速度（直接解方程，不求逆！）
//     for (unsigned dim = 0; dim < Tdim; ++dim) {
//       // 解 M·c = b_dim （替代 c = M^{-1}·b_dim）
//       Eigen::VectorXd coeffs_velo = qr.solve(rhs_velo.segment(dim * basis_size, basis_size));
//       interpolated_velo(dim) = basis_at_point.dot(coeffs_velo);
//     }
//     mtf_console_->debug("2D Linear MLS: Moment matrix trace={:.4e}, Interpolated velo=[{:.4e}, {:.4e}]",
//                         moment_matrix.trace(), interpolated_velo(0), interpolated_velo(1));
//   } catch (std::exception& exception) {
//     mtf_console_->error("mls_interpolate_velo_at_point failed: {}", exception.what());
//   }
//   return interpolated_velo;
// }

template <unsigned Tdim>
Eigen::Matrix<double, Tdim, 1> mpm::Mesh<Tdim>::weighted_average_velo_fallback(
    const std::vector<VectorDim>& valid_particle_velos,
    const std::vector<double>& valid_particle_weights) {

  VectorDim avg_velo = VectorDim::Zero();
  double total_weight = std::accumulate(valid_particle_weights.begin(), valid_particle_weights.end(), 0.0);
  const double epsilon = 1.E-15;

  if (total_weight < epsilon) {
    mtf_console_->warn("Total weight for velocity fallback is zero, return zero velocity");
    return avg_velo;
  }

  // 复用你原有加权平均逻辑，仅适配速度分量
  for (size_t i = 0; i < valid_particle_weights.size(); ++i) {
    double normalized_weight = valid_particle_weights[i] / total_weight;
    avg_velo += normalized_weight * valid_particle_velos[i];
  }

  return avg_velo;
}

//! Calculate binomial coefficient C_m^N = N! / (m!(N-m)!)
template <unsigned Tdim>
double mpm::Mesh<Tdim>::binomial_coefficient(unsigned N, unsigned m) {
  if (m < 0 || m > N) return 0.0;
  if (m == 0 || m == N) return 1.0;
  
  double result = 1.0;
  for (unsigned i = 1; i <= m; ++i) {
    result *= static_cast<double>(N - i + 1) / static_cast<double>(i);
  }
  return result;//结果正确2025.12.16 （C13 = 3）
}

//! Compute MLS weight function (cubic spline)
//! localized weighting function centered at xi(MP).
template <unsigned Tdim>
double mpm::Mesh<Tdim>::compute_mls_weight(
    double distance, double support_radius) {
  if (distance >= support_radius) return 0.0;
  
  double q = distance / support_radius;
  if (q <= 0.5) {
    return 2.0/3.0 - 4.0 * q * q + 4.0 * q * q * q;
  } else {
    return 4.0/3.0 - 4.0 * q + 4.0 * q * q - 4.0/3.0 * q * q * q;
  }
}//没问题
























//! Interpolate physical quantities at given point coordinates from particles using MLS
//！ 获取外推点所在cell及周围cell中的particles
template <unsigned Tdim>
mpm::MTFQuantities<Tdim> mpm::Mesh<Tdim>::interpolate_at_point(
    const VectorDim& point_coord,
    unsigned phase,
    int time_offset) {
  
  mpm::MTFQuantities<Tdim> quantities(time_offset);
  quantities.coordinates = point_coord;

  try {
    // Step 1: 找到包含外推点的cell
    std::shared_ptr<mpm::Cell<Tdim>> found_cell = nullptr;

    auto all_cells = this->cells();//测试成功
    
    // Use the same logic as locate_particle_cells to find the cell
    // 定位外推点在哪个cell中
    // 方法1：使用MTF专用检查（处理边界点）
    #pragma omp parallel for schedule(runtime)
    for (size_t i = 0; i < all_cells.size(); ++i) {
      auto cell = all_cells[i];
      VectorDim xi; //cell局部坐标系中的坐标
      if (!found_cell && cell->is_mtfpoint_in_cell(point_coord, &xi)) {
        #pragma omp critical
        {
          if (!found_cell) {
            found_cell = cell;
          }
        }
      }
    }

    // 方法2: 如果MTF检查也失败，使用最近cell作为fallback
    if (!found_cell) {
    mtf_console_->warn("MTF point [{}, {}] not assigned to any cell by rules, using nearest",
                      point_coord[0], point_coord[1]);
    
    double min_distance = std::numeric_limits<double>::max();
    std::shared_ptr<mpm::Cell<Tdim>> nearest_cell = nullptr;
    
    #pragma omp parallel
    {
      double local_min_distance = std::numeric_limits<double>::max();
      std::shared_ptr<mpm::Cell<Tdim>> local_nearest_cell = nullptr;
      
      #pragma omp for
      for (size_t i = 0; i < all_cells.size(); ++i) {
        auto cell = all_cells[i];
        if (!cell) continue;
        
        try {
          double distance = (cell->centroid() - point_coord).norm();
          if (distance < local_min_distance) {
            local_min_distance = distance;
            local_nearest_cell = cell;
          }
        } catch (...) {
          // 忽略错误
        }
      }
      
      // 归约操作
      #pragma omp critical
      {
        if (local_min_distance < min_distance) {
          min_distance = local_min_distance;
          nearest_cell = local_nearest_cell;
        }
      }
    }
    
    found_cell = nearest_cell;
    if (found_cell) {
      mtf_console_->debug("Using nearest cell {} (distance={})", 
                        found_cell->id(), min_distance);
    }
  }

    // Step 2: 获取包含外推点的cell及其所有相邻cell中的粒子
    std::vector<mpm::Index> all_particle_ids;
    
    // 添加主cell中的粒子
    auto main_particle_ids = found_cell->particles();//这里出错了2025.12.22 
    //1. 这里调用cell->particles()
    //2. 调用vector的拷贝构造函数
    //3. 在拷贝过程中调用size()函数时崩溃 //2026.1.5
    //这表示found_cell是空指针（大概率），或cell对象已被销毁但指针还在
    //核心原因：外推点在网格边线上，不属于任何一个cell
    all_particle_ids.insert(all_particle_ids.end(), 
                           main_particle_ids.begin(), main_particle_ids.end());
    
    // 添加相邻cell中的粒子
    for (const auto& neighbour_cell_id : found_cell->neighbours()) {

      auto neighbour_cell = this->cell(neighbour_cell_id);

      if (neighbour_cell) {
        auto neighbour_particles = neighbour_cell->particles();//测试成功
        all_particle_ids.insert(all_particle_ids.end(),
                               neighbour_particles.begin(), neighbour_particles.end());
      }
    }//测试成功

    if (all_particle_ids.empty()) {
      throw std::runtime_error("No particles found for extrapolation point");
    }

    // Step 3: 使用MLS进行插值
    quantities = this->mls_interpolate_at_point(point_coord, all_particle_ids, 
                                                phase, time_offset);

  } catch (std::exception& exception) {
    mtf_console_->error("interpolate_at_point failed: {}", exception.what());
    quantities.setZero();
  }
  
  return quantities;
}

//! MLS interpolation at point using particles from neighboring cells
//! 通过周围particles插值得到外推点物理量
template <unsigned Tdim>
mpm::MTFQuantities<Tdim> mpm::Mesh<Tdim>::mls_interpolate_at_point(
    const VectorDim& point_coord,
    const std::vector<mpm::Index>& all_particle_ids,
    unsigned phase,
    int time_offset) {
  
  mpm::MTFQuantities<Tdim> quantities(time_offset);
  quantities.coordinates = point_coord;

  try {
    // Step 1: 计算MLS支持域半径
    const double support_radius = 1.5 * this->compute_average_cell_size();//测试成功，support_radius = 0.15
    const double epsilon = 1.E-15;
    
    // 线性多项式基函数尺寸
    const unsigned basis_size = Tdim + 1;
    
    // 初始化MLS系统
    Eigen::MatrixXd moment_matrix = Eigen::MatrixXd::Zero(basis_size, basis_size);
    Eigen::VectorXd rhs_velocity = Eigen::VectorXd::Zero(basis_size * Tdim);
    Eigen::VectorXd rhs_stress = Eigen::VectorXd::Zero(basis_size * 6);
    Eigen::VectorXd rhs_density = Eigen::VectorXd::Zero(basis_size);
    
    // 收集所有有效粒子的信息
    // 这些变量是局部变量，每个外推点调用时都会重新创建，不存在重复存储的问题
    std::vector<VectorDim> valid_particle_coords;
    std::vector<VectorDim> valid_particle_velocities;
    std::vector<Eigen::Matrix<double, 6, 1>> valid_particle_stresses;
    std::vector<double> valid_particle_densities;
    std::vector<double> valid_particle_weights;

    // 显式调用clear和reserve
    valid_particle_coords.clear();
    valid_particle_velocities.clear();
    valid_particle_stresses.clear();
    valid_particle_densities.clear();
    valid_particle_weights.clear();
    
    // 预留足够空间
    size_t estimated_size = std::min(all_particle_ids.size(), size_t(100));
    valid_particle_coords.reserve(estimated_size);
    valid_particle_velocities.reserve(estimated_size);
    valid_particle_stresses.reserve(estimated_size);
    valid_particle_densities.reserve(estimated_size);
    valid_particle_weights.reserve(estimated_size);
    
    mtf_console_->debug("MLS: Starting with fresh vectors, reserved {} slots", 
                      estimated_size);

    // 第一步：收集粒子和计算权重
    for (auto particle_id : all_particle_ids) {

      //auto particle = this->particle(particle_id);

      // 修改粒子获取逻辑
      std::shared_ptr<mpm::ParticleBase<Tdim>> particle_ptr = nullptr;

      // 对mesh访问加锁
      {
        // 检查mesh指针有效性
        if (!this) {
          mtf_console_->error("MLS: this pointer is null!");
          continue;
        }
        
        // 检查粒子ID有效性
        if (particle_id >= std::numeric_limits<mpm::Index>::max() / 2) {
          mtf_console_->warn("Suspicious particle ID: {}", particle_id);
          continue;
        }
        
        // 安全获取粒子
        try {
          particle_ptr = this->particle(particle_id);
        } catch (const std::exception& e) {
          mtf_console_->error("Failed to get particle {}: {}", particle_id, e.what());
          continue;
        }
      }

      if (!particle_ptr) {
        mtf_console_->warn("Particle {} not found or null", particle_id);
        continue;
      }
      
      VectorDim particle_coord = particle_ptr->coordinates();
      
      // 计算距离和权重
      double distance = (particle_coord - point_coord).norm();//测试成功，distance = 0.0512937
      if (distance > support_radius) continue;
      
      double weight = this->compute_mls_weight(distance, support_radius);//测试成功
      if (weight < epsilon) continue;

      // 获取粒子物理量
      VectorDim particle_velocity;
      Eigen::Matrix<double, 6, 1> particle_stress;
      double particle_density;
      
      this->get_particle_quantities(particle_ptr, phase, time_offset, 
                                   particle_velocity, particle_stress, particle_density);//怎么直接跳过了？？有问题2025.12.17
      
      // 存储有效粒子信息
      valid_particle_coords.push_back(particle_coord);
      valid_particle_velocities.push_back(particle_velocity);
      valid_particle_stresses.push_back(particle_stress);
      valid_particle_densities.push_back(particle_density);
      valid_particle_weights.push_back(weight);//测试成功
    }

    if (valid_particle_weights.empty()) {
      throw std::runtime_error("No valid particles for MLS interpolation");
    }

    // 第二步：构建MLS系统
    for (size_t i = 0; i < valid_particle_weights.size(); ++i) {
      double weight = valid_particle_weights[i];
      VectorDim particle_coord = valid_particle_coords[i];
      
      // 计算基函数在粒子位置的值(相对于插值点)
      Eigen::VectorXd basis_at_particle = Eigen::VectorXd::Zero(basis_size);
      this->compute_mls_basis(point_coord, particle_coord, basis_at_particle);//计算P(xi-x)插值基函数 //有问题
      
      // 更新矩矩阵M(x)
      moment_matrix += weight * basis_at_particle * basis_at_particle.transpose();
      
      // 更新右端项 - 分别处理每个分量  //eq(2)中的b(x)
      // 速度（每个分量独立重构）
      for (unsigned dim = 0; dim < Tdim; ++dim) {
        rhs_velocity.segment(dim * basis_size, basis_size) += 
            weight * valid_particle_velocities[i](dim) * basis_at_particle;
      }
      
      // 应力（每个分量独立重构）
      for (unsigned comp = 0; comp < 6; ++comp) {
        rhs_stress.segment(comp * basis_size, basis_size) += 
            weight * valid_particle_stresses[i](comp) * basis_at_particle;
      }//测试成功，基函数3个分量*应力6个分量=18
      
      // 密度
      rhs_density += weight * valid_particle_densities[i] * basis_at_particle;
    }

    // 第三步：求解MLS系统
    // 检查矩矩阵是否可逆
    Eigen::FullPivLU<Eigen::MatrixXd> lu_decomp(moment_matrix);//通过LU分解检查矩矩阵可逆性
    if (!lu_decomp.isInvertible()) {
      mtf_console_->warn("MLS moment matrix is singular, using fallback interpolation");
      return this->weighted_average_fallback(valid_particle_coords, 
                                           valid_particle_velocities,
                                           valid_particle_stresses,
                                           valid_particle_densities,
                                           valid_particle_weights);
    }
    
    Eigen::MatrixXd moment_inv = lu_decomp.inverse();
    
    // 第四步：在插值点重构函数值
    // 计算插值点处的基函数（相对于自身，所以是[1, 0, 0, ...]）
    Eigen::VectorXd basis_at_point = Eigen::VectorXd::Zero(basis_size);
    basis_at_point(0) = 1.0;
    
    // 重构速度
    quantities.velocity.setZero();
    for (unsigned dim = 0; dim < Tdim; ++dim) {
      //基函数系数c(x)，式(2)
      Eigen::VectorXd coeffs_velocity = moment_inv * rhs_velocity.segment(dim * basis_size, basis_size);
      quantities.velocity(dim) = basis_at_point.dot(coeffs_velocity);//向量点积，式(1)
    }
    
    // 重构应力
    quantities.stress.setZero();
    for (unsigned comp = 0; comp < 6; ++comp) {
      Eigen::VectorXd coeffs_stress = moment_inv * rhs_stress.segment(comp * basis_size, basis_size);
      quantities.stress(comp) = basis_at_point.dot(coeffs_stress);
    }
    
    // 重构密度
    Eigen::VectorXd coeffs_density = moment_inv * rhs_density;
    quantities.density = basis_at_point.dot(coeffs_density);//测试成功

    mtf_console_->debug("MLS interpolation: {} particles, support_radius={}, density={}", 
                   valid_particle_weights.size(), support_radius, quantities.density);//测试成功，但step 2 时 density=0，应该加时间步判断策略

  } catch (std::exception& exception) {
    mtf_console_->error("MLS interpolation failed: {}", exception.what());
    quantities.setZero();
  }
  
  return quantities;
}

//! Get particle quantities (current or historical)
//！高阶MTF需要
template <unsigned Tdim>
void mpm::Mesh<Tdim>::get_particle_quantities(
    const std::shared_ptr<mpm::ParticleBase<Tdim>>& particle,
    unsigned phase, int time_offset,
    VectorDim& velocity, Eigen::Matrix<double, 6, 1>& stress,
    double& density) {
  
  //在这里加时间步判断逻辑
  if (time_offset == 0) {
    // 当前时间步
    velocity = particle->velocity();
    stress = particle->stress();
    density = particle->mass_density();
  } else {
    // 历史时间步
    unsigned history_step = static_cast<unsigned>(-time_offset);

    // 检查粒子是否有足够的MTF历史变量
    if (!particle->has_mtf_history()) {
      mtf_console_->warn("Particle {} has no MTF history", particle->id());
      velocity = particle->velocity();
      stress = particle->stress();
      density = particle->mass_density();
      return;
    }

    // 检查历史步是否在有效范围内（与材料真实history_size一致）
    unsigned max_history_steps = 0;
    auto material_ptr = particle->material(phase);
    if (material_ptr != nullptr && material_ptr->supports_mtf_history()) {
      const unsigned hsize = material_ptr->history_size();
      max_history_steps = (hsize > 0) ? (hsize - 1) : 0;
    }
    if (history_step > max_history_steps) {
      mtf_console_->warn("Particle {}: history_step {} exceeds max {}", 
                        particle->id(), history_step, max_history_steps);
      velocity = particle->velocity();
      stress = particle->stress();
      density = particle->mass_density();
      return;
    }
    
    try {

      velocity = particle->get_history_velocity(history_step, phase);
      // stress = particle->get_history_stress(history_step, phase);
      // density = particle->get_history_density(history_step, phase);
      stress = particle->stress();
      density = particle->mass_density();
      
      // // 检查是否获取到了有效数据
      // if (std::isnan(density) || density <= 0.0) {
      //   mtf_console_->warn("Particle {}: invalid history density at step {}", 
      //                     particle->id(), history_step);
      //   // 使用当前值作为fallback
      //   velocity = particle->velocity();
      //   stress = particle->stress();
      //   density = particle->mass_density();
      // }

      // 检查是否获取到了有效的速度数据
      if (velocity.hasNaN()) {
        mtf_console_->warn("Particle {}: invalid history velocity at step {}", 
                          particle->id(), history_step);
        velocity = particle->velocity(); // fallback 兜底
      }
      
    } catch (const std::exception& e) {
      mtf_console_->error("Failed to get particle {} history at step {}: {}", 
                         particle->id(), history_step, e.what());
      // 使用当前值作为fallback
      velocity = particle->velocity();
      stress = particle->stress();
      density = particle->mass_density();
    }
    
  }
}

//! Compute MLS basis functions
template <unsigned Tdim>
void mpm::Mesh<Tdim>::compute_mls_basis(
    const VectorDim& point_coord, const VectorDim& particle_coord,
    Eigen::VectorXd& basis_values) {
  basis_values(0) = 1.0;  // 常数项
  
  for (unsigned i = 0; i < Tdim; ++i) {
    basis_values(i + 1) = particle_coord(i) - point_coord(i);  // 线性项
  }
}

//! Fallback method: weighted average when MLS fails
template <unsigned Tdim>
mpm::MTFQuantities<Tdim> mpm::Mesh<Tdim>::weighted_average_fallback(
    const std::vector<VectorDim>& coords,
    const std::vector<VectorDim>& velocities,
    const std::vector<Eigen::Matrix<double, 6, 1>>& stresses,
    const std::vector<double>& densities,
    const std::vector<double>& weights) {
  
  mpm::MTFQuantities<Tdim> quantities(0);
  
  double total_weight = std::accumulate(weights.begin(), weights.end(), 0.0);
  if (total_weight < 1.E-15) {
    quantities.setZero();
    return quantities;
  }
  
  // 简单加权平均
  for (size_t i = 0; i < weights.size(); ++i) {
    double normalized_weight = weights[i] / total_weight;
    quantities.velocity += normalized_weight * velocities[i];
    quantities.stress += normalized_weight * stresses[i];
    quantities.density += normalized_weight * densities[i];
  }
  
  return quantities;
}

//! Compute Riemann invariant R1 (outgoing wave)
template <unsigned Tdim>
double mpm::Mesh<Tdim>::compute_riemann_invariant_R1(
    unsigned phase, const mpm::MTFQuantities<Tdim>& quantities,
    const VectorDim& boundary_normal, double wave_speed) {
  
  double R1 = 0.0;
  
  try {
    // Extract normal component of velocity
    double normal_velocity = quantities.velocity.dot(boundary_normal);
    
    // 法向应力是应力张量在边界法向方向上的投影
    double normal_stress = 0.0;

    // 将应力张量从Voigt记号转换为矩阵形式
    Eigen::Matrix<double, Tdim, Tdim> stress_tensor;
    
    if constexpr (Tdim == 1) {
      stress_tensor(0, 0) = quantities.stress[0]; // σ_xx
    } 
    else if constexpr (Tdim == 2) {
      stress_tensor(0, 0) = quantities.stress[0]; // σ_xx
      stress_tensor(1, 1) = quantities.stress[1]; // σ_yy
      stress_tensor(0, 1) = quantities.stress[2]; // σ_xy
      stress_tensor(1, 0) = quantities.stress[2]; // σ_xy
    } 
    else if constexpr (Tdim == 3) {
      stress_tensor(0, 0) = quantities.stress[0]; // σ_xx
      stress_tensor(1, 1) = quantities.stress[1]; // σ_yy
      stress_tensor(2, 2) = quantities.stress[2]; // σ_zz
      stress_tensor(1, 2) = quantities.stress[3]; // σ_yz
      stress_tensor(2, 1) = quantities.stress[3]; // σ_zy
      stress_tensor(0, 2) = quantities.stress[4]; // σ_xz
      stress_tensor(2, 0) = quantities.stress[4]; // σ_zx
      stress_tensor(0, 1) = quantities.stress[5]; // σ_xy
      stress_tensor(1, 0) = quantities.stress[5]; // σ_yx
    }
    
    // 计算法向应力: σ_n = n^T · σ · n
    normal_stress = boundary_normal.transpose() * stress_tensor * boundary_normal;

    // Riemann invariant for outgoing wave: R1 = v + σ/(ρc)
    if (quantities.density > 0.0 && wave_speed > 0.0) {
      R1 = normal_velocity + normal_stress / (quantities.density * wave_speed);
    }
    
  } catch (std::exception& exception) {
    mtf_console_->error("compute_riemann_invariant_R1 failed: {}", exception.what());
    R1 = 0.0;
  }
  
  return R1;
}

/*================================MTF=====================================*/
