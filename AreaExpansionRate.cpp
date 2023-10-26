#include <iostream>
#include <fstream>
#include <Eigen/Core>
#include <Eigen/LU>
#include "AreaExpansionRate.h"

AreaExpansionRate::AreaExpansionRate(){
	ReadInput();
	std::string fname_post = "input/" + InputPostFileName + "." + InputPostFileFormat;
	std::string fname_pre = "input/" + InputPreFileName + "." + InputPreFileFormat;

	std::cout << "Construct triangle mesh ... " << std::flush;
	HEM_post.ReadOFF(fname_post);
	if(PreSurfaceMode == 1){
		HEM_pre.ReadOFF(fname_pre);
		HEM_post.CopyHalfEdgeMesh(&HEM_post_on_pre);
		HEM_pre.WriteVTK("pre_growth_shape.vtk");
	}
	HEM_post.WriteVTK("post_growth_shape.vtk");

	std::cout << "done." << std::endl;
}

void AreaExpansionRate::ReadInput(){
	std::cout << "ReadInput ... " << std::endl;
	std::ifstream fin("input/input.txt");
	if(!fin){
		std::cout << "Error: cannot open input.txt file." << std::endl;
		exit(1);
	}
	std::string dummy;
	fin >> dummy;	fin >> PreSurfaceMode; 	
	if(PreSurfaceMode == 0){
		std::cout << "PreSurfaceMode: Map post shape on plane" << std::endl;
	}else{
		std::cout << "PreSurfaceMode: Map post shape on pre shape" << std::endl;
	}

	fin >> dummy;	fin >> InputPostFileName;
	fin >> dummy;	fin >> InputPostFileFormat;
	std::cout << "Input post file: " << InputPostFileName + "." + InputPostFileFormat << std::endl;
	fin >> dummy;	fin >> InputPreFileName;
	fin >> dummy;	fin >> InputPreFileFormat;
	if(PreSurfaceMode == 1){
		std::cout << "Input pre file: " << InputPreFileName + "." + InputPreFileFormat << std::endl;
	}

	std::string check;
	fin >> check;
	if(check != "End"){
		std::cout << "Error: invalid input.txt file." << std::endl;
		exit(1);
	}
	std::cout << "done." << std::endl;
}

void AreaExpansionRate::HarmonicMaps(){
	std::cout << "HarmonicMaps ... " << std::flush;
#pragma omp parallel sections
{
#pragma omp sections
{
	HEM_post.HarmonicParam(&HEM_post_param);
}
#pragma omp sections
{
	HEM_pre.HarmonicParam(&HEM_pre_param);
}
}
	HEM_post_param.WriteOFF("post_param.off");
	HEM_pre_param.WriteOFF("pre_param.off");
	//HEM_post_param.WriteVTK("post_param.vtk");
	//HEM_pre_param.WriteVTK("pre_param.vtk");
	std::cout << "done." << std::endl;
}

void AreaExpansionRate::HarmonicMapPost(){
	std::cout << "HarmonicMapPost ... " << std::flush;
	HEM_post.HarmonicParam(&HEM_post_param);
	HEM_post_param.WriteOFF("post_param.off");
	std::cout << "done." << std::endl;
}

bool AreaExpansionRate::HitTestPointOnFace(Eigen::Vector2d& v, Face* f){
	Eigen::Vector2d A(f->p_he->p_v->Pos(0), f->p_he->p_v->Pos(1));
	Eigen::Vector2d B(f->p_he->next->p_v->Pos(0), f->p_he->next->p_v->Pos(1));
	Eigen::Vector2d C(f->p_he->prev->p_v->Pos(0), f->p_he->prev->p_v->Pos(1));
	Eigen::Vector2d P = v;
	Eigen::Vector2d AB = B - A;
	Eigen::Vector2d BC = C - B;
	Eigen::Vector2d CA = A - C;
	Eigen::Vector2d AP = P - A;
	Eigen::Vector2d BP = P - B;
	Eigen::Vector2d CP = P - C;

	double c1 = AB(0) * BP(1) - AB(1) * BP(0);
	double c2 = BC(0) * CP(1) - BC(1) * CP(0);
	double c3 = CA(0) * AP(1) - CA(1) * AP(0);

	if((c1 > 0.0 && c2 > 0.0 && c3 > 0.0) || (c1 < 0.0 && c2 < 0.0 && c3 < 0.0)){
		return true;
	}else{
		return false;
	}
}

void AreaExpansionRate::SetCorrespPostOnPre(){
	std::cout << "SetCorrespPostOnPre ... " << std::flush;
	int n_post_bnd = (int)HEM_post_param.p_bv.size();
	int n_post_v = (int)HEM_post_param.p_v.size();
	int n_pre_f = (int)HEM_pre_param.p_f.size();
	Eigen::MatrixXd V_post_param_uv;
	V_post_param_uv.resize(n_post_v, 2);
	for(int i = 0; i < n_post_v; i++){
		V_post_param_uv(i, 0) = HEM_post_param.p_v[i]->Pos(0);
		V_post_param_uv(i, 1) = HEM_post_param.p_v[i]->Pos(1);
	}

	CorrespPostOnPre.resize(n_post_v);
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for(int i = 0; i < n_post_v; i++){
		if(HEM_post_param.p_v[i]->p_he->p_e->BoundaryFlag == 0){ // まずは内点のみ
			bool hit_flag = false;
			Eigen::Vector2d v = V_post_param_uv.row(i);
			for(int j = 0; j < n_pre_f; j++){
				Face* f = HEM_pre_param.p_f[j];
				if(HitTestPointOnFace(v, f)){
					CorrespPostOnPre(i) = j;
					hit_flag = true;
					MapPostToPre(i);
					break;
				}
			}
			if(!hit_flag){
				std::cout << "Error: vertex " << i << " cannot be assigned on any face." << std::endl;
			}
		}
	}

	// 境界点は境界辺上に
	int n_bnd_pre = (int)HEM_pre_param.p_bv.size();
	std::vector<double> theta;
	theta.resize(n_bnd_pre+1);

	std::vector<Vertex*> p_bv_loop_pre;
	p_bv_loop_pre.push_back(HEM_pre_param.p_bv[0]);
	theta[0] = atan2(p_bv_loop_pre[0]->Pos(1), p_bv_loop_pre[0]->Pos(0));
	for(int i = 0; i < n_bnd_pre - 1; i++){
		Vertex* tmp = p_bv_loop_pre[i]->p_he->next->p_v; // 境界辺のp_heは境界辺に沿ったハーフエッジにしている
		p_bv_loop_pre.push_back(tmp);
		Eigen::Vector3d axis(0.0, 0.0, 1.0);
		double delta_theta = CalcAngle(p_bv_loop_pre[i]->Pos, p_bv_loop_pre[i+1]->Pos, axis);
		theta[i+1] = theta[i] + delta_theta;
	}
	Eigen::Vector3d axis(0.0, 0.0, 1.0);
	double delta_theta = CalcAngle(p_bv_loop_pre[n_bnd_pre-1]->Pos, p_bv_loop_pre[0]->Pos, axis);
	theta[n_bnd_pre] = theta[n_bnd_pre - 1] + delta_theta;

	int n_bnd_post = (int)HEM_post_param.p_bv.size();
	for(int i = 0; i < n_bnd_post; i++){
		bool hit_flag = false;
		int vid = HEM_post_param.p_bv[i]->Id;
		double thetai = atan2(V_post_param_uv(vid, 1), V_post_param_uv(vid, 0));
		if(thetai < theta[0]) thetai += 2.0 * M_PI;
		for(int j = 0; j < n_bnd_pre; j++){
			if(theta[j] <= thetai && thetai <= theta[j+1]){
				CorrespPostOnPre(vid) = p_bv_loop_pre[j]->p_he->p_f->Id;
				double phi = thetai - theta[j];
				Eigen::Vector2d vA(p_bv_loop_pre[j]->p_he->p_v->Pos(0), p_bv_loop_pre[j]->p_he->p_v->Pos(1));
				Eigen::Vector2d vB(p_bv_loop_pre[(j+1)%n_bnd_pre]->p_he->p_v->Pos(0), p_bv_loop_pre[(j+1)%n_bnd_pre]->p_he->p_v->Pos(1));
				double cosphi = cos(phi);
				double sinphi = sin(phi);
				Eigen::Matrix2d A;
				A << vB(0) - vA(0), 
					-vA(0) * cosphi + vA(1) * sinphi,
					vB(1) - vA(1),
					-vA(0) * sinphi - vA(1) * cosphi;
				Eigen::Vector2d b(-vA(0), -vA(1));
				Eigen::Vector2d st;
				st = A.inverse() * b;		
				Eigen::Vector2d newPos = (1.0 - st(0)) * vA + st(0) * vB;

				HEM_post_param.p_v[vid]->Pos(0) = newPos(0);
				HEM_post_param.p_v[vid]->Pos(1) = newPos(1);
				MapPostToPre(vid);

				hit_flag = true;
				break;
			}
		}
		if(!hit_flag){
			std::cout << "Error: vertex " << vid << " cannot be assigned on any face. theta = " << thetai / M_PI * 180.0 << std::endl;
		}
	}

	//HEM_post_on_pre.WriteOFF("test_map.off");
	std::cout << "done." << std::endl;
}
void AreaExpansionRate::MapPostToPre(int vid){
	double s0, s1, s2, s;
	int v_id0 = HEM_pre_param.p_f[CorrespPostOnPre(vid)]->p_he->p_v->Id;
	int v_id1 = HEM_pre_param.p_f[CorrespPostOnPre(vid)]->p_he->next->p_v->Id;
	int v_id2 = HEM_pre_param.p_f[CorrespPostOnPre(vid)]->p_he->prev->p_v->Id;
	Eigen::Vector3d v0 = HEM_pre_param.p_v[v_id0]->Pos - HEM_post_param.p_v[vid]->Pos;
	Eigen::Vector3d v1 = HEM_pre_param.p_v[v_id1]->Pos - HEM_post_param.p_v[vid]->Pos;
	Eigen::Vector3d v2 = HEM_pre_param.p_v[v_id2]->Pos - HEM_post_param.p_v[vid]->Pos;

	s0 = (v1.cross(v2)).norm() * 0.5;
	s1 = (v2.cross(v0)).norm() * 0.5;
	s2 = (v0.cross(v1)).norm() * 0.5;
	s = s0 + s1 + s2;

	Eigen::Vector3d v3d0 = HEM_pre.p_v[v_id0]->Pos;
	Eigen::Vector3d v3d1 = HEM_pre.p_v[v_id1]->Pos;
	Eigen::Vector3d v3d2 = HEM_pre.p_v[v_id2]->Pos;
	Eigen::Vector3d v_mapped = (s0/s) * v3d0 + (s1/s) * v3d1 + (s2/s) * v3d2;
	HEM_post_on_pre.p_v[vid]->Pos = v_mapped;
}

void AreaExpansionRate::CalcAreaExpRatePerFace(){
	std::cout << "CalcAreaExpRatePerFace ... " << std::flush;
	int n_post_f = (int)HEM_post.p_f.size();
	AreaExpRatePerFace.resize(n_post_f);

	HEM_post.CalcFaceAreas();
	HEM_post_on_pre.CalcFaceAreas();
	for(int i = 0; i < n_post_f; i++){
		AreaExpRatePerFace[i] = HEM_post.p_f[i]->Area / HEM_post_on_pre.p_f[i]->Area;
	}

	WriteVTK(HEM_post_on_pre, "post_on_pre");
	WriteOFF(HEM_post_on_pre, "post_on_pre");
	WriteVTK(HEM_post, "post");
	std::cout << "done." << std::endl;
}
void AreaExpansionRate::CalcAreaExpRatePerFacePostFromPlane(){
	std::cout << "CalcAreaExpRatePerFacePostFromPlane ... " << std::flush;
	HEM_post_param.CopyHalfEdgeMesh(&HEM_post_on_pre);
	int n_post_f = (int)HEM_post.p_f.size();
	AreaExpRatePerFace.resize(n_post_f);
	HEM_post.CalcFaceAreas();
	HEM_post_param.CalcFaceAreas();
	for(int i = 0; i < n_post_f; i++){
		AreaExpRatePerFace[i] = HEM_post.p_f[i]->Area / HEM_post_param.p_f[i]->Area;
	}
	WriteVTK(HEM_post_param, "post_on_pre");
	WriteOFF(HEM_post_param, "post_on_pre");
	WriteVTK(HEM_post, "post");
	std::cout << "done." << std::endl;
}

void AreaExpansionRate::CalcAreaExpRatePerVertex(){
	std::cout << "CalcAreaExpRatePerVertex ... " << std::flush;
	int n_post_v = (int)HEM_post.p_v.size();
	AreaExpRatePerVertex.resize(n_post_v);

	HEM_post.CalcVertexMixedAreas();
	HEM_post_on_pre.CalcFaceAreas();
	HEM_post_on_pre.CalcIAngs();
	HEM_post_on_pre.CalcVertexMixedAreas();
	for(int i = 0; i < n_post_v; i++){
		AreaExpRatePerVertex[i] = HEM_post.p_v[i]->MixedArea / HEM_post_on_pre.p_v[i]->MixedArea;
	}
	std::cout << "done." << std::endl;
}

void AreaExpansionRate::WriteVTK(HalfEdgeMesh& hem, std::string name){
	Eigen::MatrixXd V;
	Eigen::MatrixXi F;

	int n_v = (int)hem.p_v.size();
	int n_f = (int)hem.p_f.size();

	V.resize(n_v, 3);
	F.resize(n_f, 4);

	for(int i = 0; i < n_v; i++){
		V.row(i) = hem.p_v[i]->Pos;
	}

	for(int i = 0; i < n_f; i++){
		F(i, 0) = 3;
		F(i, 1) = hem.p_f[i]->p_he->p_v->Id;
		F(i, 2) = hem.p_f[i]->p_he->next->p_v->Id;
		F(i, 3) = hem.p_f[i]->p_he->prev->p_v->Id;
	}

	std::string fname = "result_" + InputPostFileName + "_" + name + ".vtk";
	std::ofstream fout(fname);
	if(!fout){
		std::cout << "Error: cannot open output vtk file." << std::endl;
		exit(1);
	}
	fout << "# vtk DataFile Version 3.0\n";
	fout << "AreaExpansionRate\n";
	fout << "ASCII\n";
	fout << "DATASET UNSTRUCTURED_GRID\n";

	fout << "POINTS " << n_v << " float\n";
	for(int i = 0; i < n_v; i++){
		fout << V.row(i) << "\n";
	}

	fout << "CELLS " << n_f << " " << 4 * n_f << "\n";
	for(int i = 0; i < n_f; i++){
		fout << F.row(i) << "\n";
	}

	fout << "CELL_TYPES " << n_f << "\n";
	for(int i = 0; i < n_f; i++){
		fout << "5\n";
	}

	/*
	fout << "POINT_DATA " << n_v << "\n";
	fout << "SCALARS AreaExpansionRate float\n";
	fout << "LOOKUP_TABLE default\n";
	for(int i = 0; i < n_v; i++){
		fout << AreaExpRatePerVertex[i] << "\n";
	}
	*/
	fout << "CELL_DATA " << n_f << "\n";
	fout << "SCALARS AreaExpansionRate float\n";
	fout << "LOOKUP_TABLE default\n";
	for(int i = 0; i < n_f; i++){
		fout << AreaExpRatePerFace[i] << "\n";
	}

	fout.close();
}
void AreaExpansionRate::WriteOFF(HalfEdgeMesh& hem, std::string name){
	std::string fname = "result_" + InputPostFileName + "_" + name + ".off";
	std::ofstream fout(fname);
	if(!fout){
		std::cout << "Error: cannot open output off file." << std::endl;
		exit(1);
	}
	hem.WriteOFF(fname);
}

void AreaExpansionRate::WriteAreaExpRate(){
	int n_f = (int)HEM_post_on_pre.p_f.size();

	std::string fname = "area_expansion_rate.dat";
	std::ofstream fout(fname);
	if(!fout){
		std::cout << "Error: cannot open output area_expansion_rate file." << std::endl;
		exit(1);
	}
	fout << n_f << "\n" << std::endl;
	fout << AreaExpRatePerFace << std::endl;
	fout.close();
}
