#ifndef MPM_MTF_CONSTRAINT_H_
#define MPM_MTF_CONSTRAINT_H_

#include "data_types.h"

namespace mpm {

//! MTFConstraint class to store friction constraint on a set
//! \brief MTFConstraint class to store a constraint on a set
//! \details MTFConstraint stores the constraint as a static value
//! 负责接收json输入文件参数的类
class MTFConstraint {
 public:
  // Constructor
  //! \param[in] setid  set id
  //! \param[in] n 透射次数
  //! \param[in] position Nodal position along boundary
  MTFConstraint(int setid, unsigned int order, double drift_coeff = 0.0,
                mpm::Edge_Position position = mpm::Edge_Position::None)
      : setid_{setid},
        order_{order},
        drift_coeff_{drift_coeff},
        position_{position} {};//构造函数初始化列表

  // Set id
  int setid() const { return setid_; }

  //N
  unsigned int order() const { return order_; }

  // 零频消飘修正系数
  double drift_coeff() const { return drift_coeff_; }

  // Return position
  mpm::Edge_Position position() const { return position_; }

 private:
  // ID
  int setid_;
  // N
  unsigned int order_; //透射次数

  // 零频消飘修正系数
  double drift_coeff_;

  // Node position
  mpm::Edge_Position position_;
};
}  // namespace mpm
#endif  // MPM_MTF_CONSTRAINT_H_ //lxd