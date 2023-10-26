#ifndef _AREAEXPANSIONRATE_H_
#define _AREAEXPANSIONRATE_H_

#include <string>
#include "HalfEdgeMesh.h"

class AreaExpansionRate{
public:
	int PreSurfaceMode; // 0: planar,  1: input_pre_file
private:
	HalfEdgeMesh HEM_post;
	HalfEdgeMesh HEM_pre;

	HalfEdgeMesh HEM_post_param;
	HalfEdgeMesh HEM_pre_param;

	HalfEdgeMesh HEM_post_on_pre;
	Eigen::VectorXi CorrespPostOnPre; //(V_post_uvの各頂点がF_preのどれに属するのかの対応)

	std::string InputPostFileName;
	std::string InputPreFileName;
	std::string InputPostFileFormat;
	std::string InputPreFileFormat;

	Eigen::VectorXd AreaExpRatePerFace;
	Eigen::VectorXd AreaExpRatePerVertex;

public:
	AreaExpansionRate();
	double CalcAngle(const Eigen::Vector3d& v, const Eigen::Vector3d& w, const Eigen::Vector3d& axis){ // [-π, π]
		// 参考: https://www.cs.utexas.edu/users/evouga/uploads/4/5/6/8/45689883/turning.pdf
		double theta = 2.0 * atan2((v.cross(w).dot(axis) / axis.norm()), v.dot(w) + v.norm() * w.norm());
		return theta;
	}
	void ReadInput();
	void HarmonicMaps();
	void HarmonicMapPost();
	bool HitTestPointOnFace(Eigen::Vector2d&, Face*);
	void SetCorrespPostOnPre();
	void MapPostToPre(int);
	void CalcAreaExpRatePerFace();
	void CalcAreaExpRatePerVertex();
	void CalcAreaExpRatePerFacePostFromPlane();

	void WriteVTK(HalfEdgeMesh&, std::string);
	void WriteOFF(HalfEdgeMesh&, std::string);

	void WriteAreaExpRate();
};

#endif
