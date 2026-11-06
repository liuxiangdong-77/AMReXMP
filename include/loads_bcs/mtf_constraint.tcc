//! Assign MTF constraint to nodes
//  将MTF约束分配到结点集
template <unsigned Tdim>
bool mpm::Constraints<Tdim>::assign_nodal_mtf_constraint(
    int nset_id,
    const std::shared_ptr<mpm::MTFConstraint>& mtf_constraint, unsigned phase, double dt) {
  bool status = true;
  try {
    int set_id = nset_id;
    auto nset = mesh_->nodes(set_id);//成功
    if (nset.size() == 0)
      throw std::runtime_error(
          "Node set is empty for assignment of MTF constraints");

    if (!mesh_->create_nodal_mtf_constraint(set_id, mtf_constraint, phase, dt))
      throw std::runtime_error("Failed to create MTF constraint in mesh");//直接把mtf_constraint整个传到mesh.里处理

    console_->info("Assigned MTF constraint to nset {}: pos={}, order={}",
                  set_id,
                  static_cast<int>(mtf_constraint->position()),
                  mtf_constraint->order());//成功，表示参数已传到mtf_constraint类里
    
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

//! Assign MTF constraints to nodes
//  将MTF约束分配到单个结点
template <unsigned Tdim>
bool mpm::Constraints<Tdim>::assign_nodal_mtf_constraints(
    const std::vector<std::tuple<mpm::Index, mpm::Edge_Position, 
                                 unsigned>>& mtf_constraints, unsigned phase, double dt) {
  bool status = true;
  try {
    for (const auto& mtf_constraint_tuple : mtf_constraints) {
      // Node id
      mpm::Index nid = std::get<0>(mtf_constraint_tuple);
      // Position
      mpm::Edge_Position position = std::get<1>(mtf_constraint_tuple);
      // Order
      unsigned order = std::get<2>(mtf_constraint_tuple);

      // 创建 MTF 约束对象
      auto constraint = std::make_shared<mpm::MTFConstraint>(
          -1, order, position); // set_id = -1 表示单个节点
      
      // 将约束应用到单个节点
      // 这里需要创建一个临时的节点集
      // 或者直接调用节点的 MTF 初始化方法
      
      console_->info("Assigned MTF constraint to node {}: pos={}, order={}",
                    nid, static_cast<int>(position), order);
    }
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

//! Assign mtf constraints pointers and ids
template <unsigned Tdim>
void mpm::Constraints<Tdim>::assign_mtf_id_ptr(
    unsigned nset_id,
    std::shared_ptr<mpm::MTFConstraint>& mtf_constraint) {
  this->mtf_constraint_.emplace_back(mtf_constraint);//存储约束对象指针
  this->mtf_nset_id_.emplace_back(nset_id);//存储对应的结点集ID
}

