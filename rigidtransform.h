#include <utility>
#include <Eigen/SVD>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/Geometry>
#include <vector>
#include <cassert>

Eigen::Vector3d RigidTransform(Eigen::Vector3d local, Eigen::Matrix3d rotation, Eigen::Vector3d transform) {
    return (rotation * local + transform);
}

Eigen::Matrix3d Quaterniond2Rotational(Eigen::Quaterniond q) {
    return q.normalized().toRotationMatrix();
}

typedef std::pair<Eigen::Matrix3d, Eigen::Vector3d> TransformType;
typedef std::vector<Eigen::Vector3d>                PointsType;

TransformType computeRigidTransform(const PointsType& src, const PointsType& dst, const unsigned int pair_point_num = 0) {
	assert(src.size() == dst.size());

  auto pair_size = 0;
  if (0 == pair_point_num)
    pair_size = src.size();
  else
    pair_size = pair_point_num;
	auto center_src = Eigen::Vector3d{0, 0, 0}; 
  auto center_dst = Eigen::Vector3d{0, 0, 0};

	for (auto i = 0; i<pair_size; ++ i) {
		center_src += src[i];
		center_dst += dst[i];
	}
  // cal center;
	center_src /= (double)pair_size;
	center_dst /= (double)pair_size;

	auto S = Eigen::MatrixXd(pair_size, 3); 
  auto D = Eigen::MatrixXd(pair_size, 3);

	for (auto i = 0; i < pair_size; ++ i) {
		for (int j = 0; j < 3; ++ j)
			S(i, j) = src[i][j] - center_src[j];
		for (int j = 0; j < 3; ++ j)
			D(i, j) = dst[i][j] - center_dst[j];
	}

	auto Dt = D.transpose();
	auto H = Dt*S;
	Eigen::Matrix3d W, U, V;

	Eigen::JacobiSVD<Eigen::MatrixXd> svd;
	auto H_ = Eigen::MatrixXd(3, 3);
	for (auto i = 0; i < 3; ++ i) 
    for (int j = 0; j < 3; ++ j)
      H_(i, j) = H(i, j);
	svd.compute(H_, Eigen::ComputeThinU | Eigen::ComputeThinV );
	if (!svd.computeU() || !svd.computeV()) {
		//std::cerr << "decomposition error" << endl;
		return std::make_pair(Eigen::Matrix3d::Identity(), Eigen::Vector3d::Zero());
	}
	auto Vt = svd.matrixV().transpose();
	R = svd.matrixU()*Vt;
	Eigen::Vector3d t = center_dst - R*center_src;	
	
	return std::make_pair(R, t);
}
