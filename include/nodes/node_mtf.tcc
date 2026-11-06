//! Compute boundary normal vector based on boundary position
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
Eigen::Matrix<double, Tdim, 1> mpm::Node<Tdim, Tdof, Tnphases>::compute_boundary_normal(
    mpm::Edge_Position boundary_position) const {
  
  VectorDim normal = VectorDim::Zero();
  
  switch (boundary_position) {
    // 边界面
    case mpm::Edge_Position::Left:
      normal(0) = -1.0;  // 左边界，法向指向左
      break;
    case mpm::Edge_Position::Right:
      normal(0) = 1.0;   // 右边界，法向指向右
      break;
    case mpm::Edge_Position::Bottom:
      if constexpr (Tdim >= 2) normal(1) = -1.0;  // 下边界，法向指向下
      break;
    case mpm::Edge_Position::Top:
      if constexpr (Tdim >= 2) normal(1) = 1.0;   // 上边界，法向指向上
      break;
    case mpm::Edge_Position::Front:
      if constexpr (Tdim >= 3) normal(2) = -1.0;  // 前边界，法向指向前
      break;
    case mpm::Edge_Position::Back:
      if constexpr (Tdim >= 3) normal(2) = 1.0;   // 后边界，法向指向后
      break;
      
    // 2D角点
    case mpm::Edge_Position::BottomLeft:
      if constexpr (Tdim >= 2) normal << -1.0, -1.0;  // 左下角
      break;
    case mpm::Edge_Position::TopLeft:
      if constexpr (Tdim >= 2) normal << -1.0, 1.0;   // 左上角
      break;
    case mpm::Edge_Position::BottomRight:
      if constexpr (Tdim >= 2) normal << 1.0, -1.0;   // 右下角
      break;
    case mpm::Edge_Position::TopRight:
      if constexpr (Tdim >= 2) normal << 1.0, 1.0;    // 右上角
      break;
      
    // 3D角点（边角）
    case mpm::Edge_Position::BottomLeftFront:
      if constexpr (Tdim >= 3) normal << -1.0, -1.0, -1.0;  // 左前下角
      break;
    case mpm::Edge_Position::BottomLeftBack:
      if constexpr (Tdim >= 3) normal << -1.0, -1.0, 1.0;   // 左后下角
      break;
    case mpm::Edge_Position::TopLeftFront:
      if constexpr (Tdim >= 3) normal << -1.0, 1.0, -1.0;   // 左前上角
      break;
    case mpm::Edge_Position::TopLeftBack:
      if constexpr (Tdim >= 3) normal << -1.0, 1.0, 1.0;    // 左后上角
      break;
    case mpm::Edge_Position::BottomRightFront:
      if constexpr (Tdim >= 3) normal << 1.0, -1.0, -1.0;   // 右前下角
      break;
    case mpm::Edge_Position::BottomRightBack:
      if constexpr (Tdim >= 3) normal << 1.0, -1.0, 1.0;    // 右后下角
      break;
    case mpm::Edge_Position::TopRightFront:
      if constexpr (Tdim >= 3) normal << 1.0, 1.0, -1.0;    // 右前上角
      break;
    case mpm::Edge_Position::TopRightBack:
      if constexpr (Tdim >= 3) normal << 1.0, 1.0, 1.0;     // 右后上角
      break;
      
    default:
      normal.setZero();
      break;
  }
  
  return normal.normalized(); // 返回单位化后的法向量
}

//! Compute instantaneous wave speed at current boundary node
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
Eigen::Matrix<double, Tdim, 1> mpm::Node<Tdim, Tdof, Tnphases>::compute_instantaneous_wave_speed(
    unsigned phase) const {
  
  Eigen::Matrix<double, Tdim, 1> wave_speed = VectorDim::Zero();
  const double tolerance = 1.E-16;

  // 调试信息
  if (console_ && console_->level() <= spdlog::level::debug) {
    console_->debug("Node {}: computing wave speed, mass = {}, material_ids size = {}", 
                    id(), mass(phase), material_ids_.size());
  }
  
  if (mass(phase) > tolerance && material_ids_.size() > 0) {
    // Get material properties from the first material
    auto mat_id = material_ids_.begin();//没问题，mat_id = 0;
    
    // Get wave velocities from material properties
    VectorDim wave_velocity = property_handle_->property(
        "wave_velocities", prop_id_, *mat_id, 2); //MatrixXd: 动态大小矩阵类型； 2: 属性的维度 //??波速这里有问题？

    //提取P波和S波属性
    double pwave_v = this->property_handle_->property(
          "wave_velocities", prop_id_, *mat_id, 2)(0);
    double swave_v = this->property_handle_->property(
          "wave_velocities", prop_id_, *mat_id, 2)(1);
    
      pwave_v /= this->mass(*mat_id);
      swave_v /= this->mass(*mat_id);

      if constexpr (Tdim == 2) {
        wave_speed << pwave_v, swave_v;
      } else {
        wave_speed << pwave_v, swave_v, swave_v;
      }
  }
  
  return wave_speed;
}

template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::compute_rayleigh_damp(
    unsigned phase, double alpha) {

  Eigen::Matrix<double, Tdim, 1> rayleigh_force = Eigen::Matrix<double, Tdim, 1>::Zero();

  rayleigh_force = alpha * this->mass(phase) * this->velocity(phase);

  this->update_external_force(true, phase, -rayleigh_force);
}

template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::update_mtf_velocity(
    VectorDim& new_velocity,  unsigned phase) {
      // Apply the MTF-updated velocity
      node_mutex_.lock();
      mtf_velocity_.col(phase) = new_velocity;
      node_mutex_.unlock();
}//lxd

template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::update_displacement_mtf(
    VectorDim& new_displacement, unsigned phase){
      node_mutex_.lock();
      displacement_.col(phase) = new_displacement;
      node_mutex_.unlock();
}

template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::update_mtf_acceleration(
    VectorDim& mtf_accel, unsigned phase){
      node_mutex_.lock();
      mtf_acceleration_.col(phase) = mtf_accel;
      node_mutex_.unlock();
}

template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
double mpm::Node<Tdim, Tdof, Tnphases>::compute_mtf_density(
    unsigned phase) const {
  const double tolerance = 1.E-16;
  if (mass(phase) > tolerance && material_ids_.size() > 0) {
    auto mat_id = material_ids_.begin();
    double density = this->property_handle_->property(
        "density", prop_id_, *mat_id)(0);
    density /= this->mass(*mat_id);
    return density;
  }
  return 0.0;
}
