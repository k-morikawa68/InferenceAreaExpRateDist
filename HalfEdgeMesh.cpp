#include <iostream>
#include <fstream>
#include <algorithm>
#include "HalfEdgeMesh.h"

void HalfEdgeMesh::CopyHalfEdgeMesh(HalfEdgeMesh* hem){
	int n_v = (int)p_v.size();
	int n_f = (int)p_f.size();
	int n_e = (int)p_e.size();
	int n_he = (int)p_he.size();
	int n_bnd = (int)p_bv.size();
	hem->p_v.resize(n_v);
	hem->p_f.resize(n_f);
	hem->p_e.resize(n_e);
	hem->p_he.resize(n_he);
	hem->p_bv.resize(n_bnd);

	for(int i = 0; i < n_v; i++){
		hem->p_v[i] = new Vertex(p_v[i]->Pos);
		hem->p_v[i]->Id = i;
		hem->p_v[i]->MixedArea = p_v[i]->MixedArea;
		hem->p_v[i]->NormalVec = p_v[i]->NormalVec;
	}

	for(int i = 0; i < n_f; i++){
		int j = 0;
		HalfEdge* he = p_f[i]->p_he;
		do{
			int he_id = 3*i + j;
			int v_id = he->p_v->Id;
			int e_id = he->p_e->Id;
			hem->p_he[he_id] = new HalfEdge(hem->p_v[v_id]);
			hem->p_he[he_id]->Id = he_id;
			hem->p_he[he_id]->IAng = p_he[he_id]->IAng;
			hem->p_he[he_id]->Cotan = p_he[he_id]->Cotan;
			if(j == 0){
				hem->p_f[i] = new Face(hem->p_he[he_id]);
				hem->p_f[i]->Id = p_f[i]->Id;
				hem->p_f[i]->Area = p_f[i]->Area;
				hem->p_f[i]->NormalVec = p_f[i]->NormalVec;
				hem->p_f[i]->FFMat = p_f[i]->FFMat;
				hem->p_f[i]->SFMat = p_f[i]->SFMat;
			}
			
			he = he->next;
			j++;
		}while(he != p_f[i]->p_he);
	}

	for(int i = 0; i < n_he; i++){
		if(p_he[i]->pair != nullptr) hem->p_he[i]->pair = hem->p_he[p_he[i]->pair->Id];
		hem->p_he[i]->next = hem->p_he[p_he[i]->next->Id];
		hem->p_he[i]->prev = hem->p_he[p_he[i]->prev->Id];

		hem->p_he[i]->p_f = hem->p_f[p_he[i]->p_f->Id];
	}

	for(int i = 0; i < n_e; i++){
		hem->p_e[i] = new Edge(hem->p_he[p_e[i]->p_he_left->Id]);
	}

	for(int i = 0; i < n_v; i++){
		hem->p_v[i]->p_he = hem->p_he[p_v[i]->p_he->Id];
	}

	for(int i = 0; i < n_bnd; i++){
		hem->p_bv[i] = hem->p_v[p_bv[i]->Id];
	}
}
void HalfEdgeMesh::ReadOFF(std::string fname){
	std::ifstream fin(fname);
	if(!fin){
		std::cout << "Error: cannot open input mesh file." << std::endl;
		exit(1);
	}
	std::string str;
	fin >> str;
	int n_v, n_f, n_e;
	fin >> n_v >> n_f >> n_e;

	for(int i = 0 ; i < n_v; i++){
		double x, y, z;
		fin >> x >> y >> z;
		Vertex* tmp_v = new Vertex(x, y, z);
		p_v.push_back(tmp_v);
	}

	for(int i = 0; i < n_f; i++){
		int n;
		int v0, v1, v2;
		fin >> n >> v0 >> v1 >> v2;
		AddFace(p_v[v0], p_v[v1], p_v[v2]);
	}

	SetVertexIds();
	SetFaceIds();

	// 境界付きの場合、頂点周りのハーフエッジを反時計回りに巡回するときに、境界に当たったところで止まってしまう
	// ので、境界頂点のメンバとするハーフエッジは、境界辺のハーフエッジ(反時計回りに巡回させるときの始点)にしておく
	// そのため、頂点周りの巡回をするときは、必ずhei = hei->prev->pair(反時計回り)と回さないといけないので注意
	// さらに、そのときは、whileの条件に hei(->prev->pair)がnullptrでないことを入れるとともに、最後にhei->prevで処理を追加する必要があるので要注意（ミスを多発するので何とかしたい）
	for(int i = 0; i < n_v; i++){
		HalfEdge* hei = p_v[i]->p_he;	
		do{
			if(hei->prev->pair == nullptr){ // 境界点
				HalfEdge* hei2 = p_v[i]->p_he;
				p_bv.push_back(p_v[i]);
				if(hei2->pair != nullptr){
					do{
						hei2 = hei2->pair->next; // 時計回りに巡回
						p_v[i]->p_he = hei2;
					}while(hei2->pair != nullptr);
					break;
				}else{
					break;
				}

			}
			hei = hei->prev->pair; // 反時計回りに巡回
		}while(hei != p_v[i]->p_he);
	}

	ConstructEdge();
	SetEdgeIds();

	SetHalfEdgeIds();
}

void HalfEdgeMesh::WriteOFF(std::string fname){
	Eigen::MatrixXd V;
	Eigen::MatrixXi F;

	int n_v = (int)p_v.size();
	int n_f = (int)p_f.size();

	std::ofstream fout(fname);
	if(!fout){
		std::cout << "Error: cannot open output mesh file." << std::endl;
		exit(1);
	}

	V.resize(n_v, 3);
	F.resize(n_f, 4);

	for(int i = 0; i < n_v; i++){
		V.row(i) = p_v[i]->Pos;
	}

	SetVertexIds();
	for(int i = 0; i < n_f; i++){
		F(i, 0) = 3;
		F(i, 1) = p_f[i]->p_he->p_v->Id;
		F(i, 2) = p_f[i]->p_he->next->p_v->Id;
		F(i, 3) = p_f[i]->p_he->prev->p_v->Id;
	}

	fout << "OFF\n" << n_v << " " << n_f << " 0\n";
	fout << V << "\n";
	fout << F << "\n";
	fout.close();
}

void HalfEdgeMesh::WriteVTK(std::string fname){
	Eigen::MatrixXd V;
	Eigen::MatrixXi F;

	int n_v = (int)p_v.size();
	int n_f = (int)p_f.size();

	V.resize(n_v, 3);
	F.resize(n_f, 4);

	for(int i = 0; i < n_v; i++){
		V.row(i) = p_v[i]->Pos;
	}

	SetVertexIds();
	for(int i = 0; i < n_f; i++){
		F(i, 0) = 3;
		F(i, 1) = p_f[i]->p_he->p_v->Id;
		F(i, 2) = p_f[i]->p_he->next->p_v->Id;
		F(i, 3) = p_f[i]->p_he->prev->p_v->Id;
	}

	std::ofstream fout(fname);
	if(!fout){
		std::cout << "Error: cannot open output vtk file." << std::endl;
		exit(1);
	}
	fout << "# vtk DataFile Version 3.0\n";
	fout << "HalfEdgeMesh\n";
	fout << "ASCII\n";
	fout << "DATASET UNSTRUCTURED_GRID\n";

	fout << "POINTS " << n_v << " float\n";
	fout << V << "\n";

	fout << "CELLS " << n_f << " " << 4 * n_f << "\n";
	fout << F << "\n";

	fout << "CELL_TYPES " << n_f << "\n";
	for(int i = 0; i < n_f; i++){
		fout << "5\n";
	}

	fout.close();
}

void HalfEdgeMesh::AddFace(Vertex* v0, Vertex* v1, Vertex* v2){
	HalfEdge* he0 = new HalfEdge(v0);
	HalfEdge* he1 = new HalfEdge(v1);
	HalfEdge* he2 = new HalfEdge(v2);

	he0->next = he1; he0->prev = he2;
	he1->next = he2; he1->prev = he0;
	he2->next = he0; he2->prev = he1;

	Face* tmp_f = new Face(he0);
	p_f.push_back(tmp_f);
	he0->p_f = tmp_f;
	he1->p_f = tmp_f;
	he2->p_f = tmp_f;

	SetHalfEdgePair(he0);
	SetHalfEdgePair(he1);
	SetHalfEdgePair(he2);

	p_he.push_back(he0);
	p_he.push_back(he1);
	p_he.push_back(he2);
}

void HalfEdgeMesh::SetVertexIds(){
	int n_v = (int)p_v.size();
	for(int i = 0; i < n_v; i++){
		p_v[i]->Id = i;
	}
}
void HalfEdgeMesh::SetFaceIds(){
	int n_f = (int)p_f.size();
	for(int i = 0; i < n_f; i++){
		p_f[i]->Id = i;
	}
}
void HalfEdgeMesh::SetEdgeIds(){
	int n_e = (int)p_e.size();
	for(int i = 0; i < n_e; i++){
		p_e[i]->Id = i;
	}
}
void HalfEdgeMesh::SetHalfEdgeIds(){
	int n_he = (int)p_he.size();
	for(int i = 0; i < n_he; i++){
		p_he[i]->Id = i;
	}
}

void HalfEdgeMesh::SetHalfEdgePair(HalfEdge* he){
	std::vector<Face*>::iterator it_f;
	for(it_f = p_f.begin(); it_f != p_f.end(); it_f++){
		HalfEdge* he_in_face = (*it_f)->p_he;
		do{
			if(he->p_v == he_in_face->next->p_v && he->next->p_v == he_in_face->p_v){
				he->pair = he_in_face;
				he_in_face->pair = he;
				return;
			}
			he_in_face = he_in_face->next;
		}while(he_in_face != (*it_f)->p_he);
	}
}

void HalfEdgeMesh::ConstructEdge(){
	int n_f = p_f.size();
	for(int i = 0; i < n_f; i++){
		Face* f = p_f[i];
		HalfEdge* he = f->p_he;
		do{
			Edge* tmp_e = new Edge(he);
			int flag = 0;
			int currentNumEdge = p_e.size();
			if(tmp_e->BoundaryFlag != 1){
				for(int k = 0; k < currentNumEdge; k++){
					if(tmp_e->p_he_right == p_e[k]->p_he_left){
						flag = 1;
						delete tmp_e;
						break;
					}
				}
			}
			if(flag != 1){
				p_e.push_back(tmp_e);
			}
			he = he->next;
		}while(he != f->p_he);
	}
}

void HalfEdgeMesh::EdgeCollapse(HalfEdge* he){
	HalfEdge* heLT = he->prev->pair; // (8)
	HalfEdge* heRT = he->next->pair; // (7)
	HalfEdge* heLB = he->pair->next->pair; // (9)
	HalfEdge* heRB = he->pair->prev->pair; // (10)

	// 移動後の頂点の座標値を設定
	he->next->p_v->Pos = 0.5 * (he->p_v->Pos + he->next->p_v->Pos);

	// 頂点 → ハーフエッジの再設定
	// 頂点が参照するハーフエッジが削除してしまうかもしれないので設定し直す
	he->next->p_v->p_he = heRB; // (2)のハーフエッジの始点が参照するハーフエッジを(10)に
	he->prev->p_v->p_he = heRT; // (3)のハーフエッジの始点が参照するハーフエッジを(7)に
	he->pair->prev->p_v->p_he = heLB; // (6)のハーフエッジの始点が参照するハーフエッジを(9)に

	// ハーフエッジ → 頂点の更新
	// (1)の始点の頂点周りのハーフエッジを巡回しながら参照頂点を付け替える
	HalfEdge* he1 = he->prev->pair;
	do {
		he1->p_v = he->next->p_v;
		he1 = he1->prev->pair;
	}while(he1 != he);

	// pair 関係の更新
	SetHalfEdgePair(heRT, heLT);
	SetHalfEdgePair(heLB, heRB);

	// 不要になった要素の削除
	DeleteEdge(he->p_e);
	DeleteVertex(he->p_v);
	DeleteFace(he->pair->p_f);
	DeleteFace(he->p_f);
}

void HalfEdgeMesh::EdgeFlip(HalfEdge* he){
	HalfEdge* heLT = he->prev->pair;
	HalfEdge* heRT = he->next->pair;
	HalfEdge* heLB = he->pair->next->pair;
	HalfEdge* heRB = he->pair->prev->pair;

	// 頂点 → ハーフエッジの再設定
	he->next->p_v->p_he = heRB;
	he->prev->p_v->p_he = heRT;
	he->p_v->p_he = heLT;
	he->pair->prev->p_v->p_he = heLB;

	// ハーフエッジ → 頂点の更新
	he->p_v = heLB->p_v;
	he->next->p_v = heRT->p_v;
	he->prev->p_v = heLT->p_v;
	he->pair->p_v = heRT->p_v;
	he->pair->next->p_v = heLB->p_v;
	he->pair->prev->p_v = heRB->p_v;

	// pair 関係の更新
	SetHalfEdgePair(heLT, he->next);
	SetHalfEdgePair(heLB, he->prev);
	SetHalfEdgePair(heRT, he->pair->prev);
	SetHalfEdgePair(heRB, he->pair->next);

	// Edgeの更新も必要
}

void HalfEdgeMesh::EdgeSplit(HalfEdge* he){
	Eigen::Vector3d new_v_pos = 0.5 * (he->p_v->Pos + he->next->p_v->Pos);
	Vertex* new_v = new Vertex(new_v_pos);
	p_v.push_back(new_v);

	he->next->p_v = new_v;
	he->pair->p_v = new_v;

	HalfEdge* he1 = he;
	HalfEdge* he2 = he->next;
	HalfEdge* he3 = he->prev;
	HalfEdge* he4 = he->pair;
	HalfEdge* he5 = he->pair->next;
	HalfEdge* he6 = he->pair->prev;
	HalfEdge* he7 = he->next->pair;
	HalfEdge* he8 = he->prev->pair;
	HalfEdge* he9 = he->pair->next->pair;
	HalfEdge* he10 = he->pair->prev->pair;
	HalfEdge* he11 = new HalfEdge(new_v);
	HalfEdge* he12 = new HalfEdge(he10->p_v);
	HalfEdge* he13 = new HalfEdge(he7->p_v);
	HalfEdge* he14 = new HalfEdge(he10->p_v);
	HalfEdge* he15 = new HalfEdge(new_v);
	HalfEdge* he16 = new HalfEdge(he9->p_v);

	he11->next = he12;	he11->prev = he13;
	he12->next = he13;	he12->prev = he11;
	he13->next = he11;	he13->prev = he12;
	he14->next = he15;	he14->prev = he16;
	he15->next = he16;	he15->prev = he14;
	he16->next = he14;	he16->prev = he15;

	// pair 関係の更新
	he11->pair = he14;
	he12->pair = he7;	he7->pair = he12;
	he13->pair = he2;	he2->pair = he13;
	he14->pair = he11;
	he15->pair = he6;	he6->pair = he15;
	he16->pair = he10;	he10->pair = he16;

	// 頂点 → ハーフエッジの再設定
	new_v->p_he = he2;
	he10->p_v->p_he = he10; // he10->p_vがhe2かhe4だったらいけないので

	Face* tmp_f1 = new Face(he11);
	p_f.push_back(tmp_f1);
	he11->p_f = tmp_f1;
	he12->p_f = tmp_f1;
	he13->p_f = tmp_f1;
	Face* tmp_f2 = new Face(he14);
	p_f.push_back(tmp_f2);
	he14->p_f = tmp_f2;
	he15->p_f = tmp_f2;
	he16->p_f = tmp_f2;

	// Edgeの更新も必要
}

void HalfEdgeMesh::CalcFaceArea(Face* f){
	Eigen::Vector3d v1 = f->p_he->next->p_v->Pos - f->p_he->p_v->Pos;
	Eigen::Vector3d v2 = f->p_he->prev->p_v->Pos - f->p_he->p_v->Pos;
	Eigen::Vector3d cross12 = v1.cross(v2);
	f->Area = 0.5 * cross12.norm();
}
void HalfEdgeMesh::CalcFaceNormal(Face* f){
	Eigen::Vector3d v1 = f->p_he->next->p_v->Pos - f->p_he->p_v->Pos;
	Eigen::Vector3d v2 = f->p_he->prev->p_v->Pos - f->p_he->p_v->Pos;
	Eigen::Vector3d cross12 = v1.cross(v2);
	cross12.normalize();
	f->NormalVec = cross12;
}
void HalfEdgeMesh::CalcDAng(Edge* e){
	if(e->BoundaryFlag == 0){
		HalfEdge* he1 = e->p_he_left;
		HalfEdge* he2 = e->p_he_right;
		Eigen::Vector3d v1 = he1->p_v->Pos;
		Eigen::Vector3d v2 = he2->p_v->Pos;
		Eigen::Vector3d axis = v2 - v1;
		Face* f1 = he1->p_f;
		Face* f2 = he2->p_f;
		double ang = CalcAngle(f1->NormalVec, f2->NormalVec, axis);
		e->DAng = ang;
	}else{
		e->DAng = 0.0; // 境界辺
	}
}
void HalfEdgeMesh::CalcIAng(HalfEdge* he){
	Eigen::Vector3d v1 = he->p_v->Pos - he->prev->p_v->Pos;	
	Eigen::Vector3d v2 = he->next->p_v->Pos - he->prev->p_v->Pos;	
	Eigen::Vector3d axis = v1.cross(v2);
	double ang = CalcAngle(v1, v2, axis);
	he->IAng = ang;
}
void HalfEdgeMesh::CalcCotan(HalfEdge* he){
	Eigen::Vector3d v0 = he->prev->p_v->Pos;
	Eigen::Vector3d v1 = he->p_v->Pos;
	Eigen::Vector3d v2 = he->next->p_v->Pos;
	Eigen::Vector3d e1 = v1 - v0;
	Eigen::Vector3d e2 = v2 - v0;
	double e1norm = e1.norm();
	double e2norm = e2.norm();
	double c = e1.dot(e2) / (e1norm * e2norm);
	Eigen::Vector3d e1crosse2;
	e1crosse2 = e1.cross(e2);
	double s = e1crosse2.norm() / (e1norm * e2norm);
	he->Cotan = c / s;
}
void HalfEdgeMesh::CalcVertexMixedArea(Vertex* v){
	// 参考: http://rodolphe-vaillant.fr/entry/20/compute-harmonic-weights-on-a-triangular-mesh#mixed_voro_area
	double sum = 0.0;	
	HalfEdge* he = v->p_he;
	do{
		double deg90 = 0.5 * M_PI;
		if(he->IAng < deg90 && he->next->IAng < deg90 && he->prev->IAng < deg90){
			// 鋭角三角形なので、Voronoi
			Eigen::Vector3d p = v->Pos;
			Eigen::Vector3d q = he->next->p_v->Pos;
			Eigen::Vector3d r = he->prev->p_v->Pos;
			double ang_q = he->prev->IAng;
			double ang_r = he->IAng;
			double cot_q = cos(ang_q) / sin(ang_q);
			double cot_r = cos(ang_r) / sin(ang_r);
			double tmp_r = (p - r).squaredNorm();
			double tmp_q = (p - q).squaredNorm();
			sum += 0.125 * (tmp_r * cot_q + tmp_q * cot_r);
		}else{
			if(he->next->IAng > deg90){
				// vが鈍角の頂点
				sum += 0.5 * he->p_f->Area;
			}else{
				// vではない頂点が鈍角
				sum += 0.25 * he->p_f->Area;
			}
		}

		HalfEdge* he_tmp = he;
		he = he->prev->pair;
	}while(he != v->p_he && he != nullptr);
	v->MixedArea = sum;
}
void HalfEdgeMesh::CalcFFMat(Face* f){
	Eigen::Vector3d I_entry;
	Eigen::Vector3d e[3];
	Eigen::Vector3d v[3];

	v[0] = f->p_he->p_v->Pos;
	v[1] = f->p_he->next->p_v->Pos;
	v[2] = f->p_he->prev->p_v->Pos;

	for(int i = 0; i < 3; i++){
		e[i] = v[(i+2)%3] - v[(i+1)%3];
		I_entry(i) = e[i].dot(e[i]);
	}

	double E, F, G;
	E = I_entry(0);
	F = -0.5 * (I_entry(0) + I_entry(1) - I_entry(2));
	G = I_entry(1);
	Eigen::Matrix2d result;
	result << E, F, F, G;
	f->FFMat = result;
}
void HalfEdgeMesh::CalcSFMat(Face* f){
	Eigen::Vector3d II_entry;
	Eigen::Vector3d e[3];
	Eigen::Vector3d n[3];
	Eigen::Vector3d v[3];
	Eigen::Vector3d n_opp[3];
	Eigen::Vector3d n_T;

	v[0] = f->p_he->p_v->Pos;
	v[1] = f->p_he->next->p_v->Pos;
	v[2] = f->p_he->prev->p_v->Pos;

	n_T = f->NormalVec;
	if(f->p_he->next->pair != nullptr){
		n_opp[0] = f->p_he->next->pair->p_f->NormalVec;
	}else{
		n_opp[0] = n_T;
	}
	if(f->p_he->prev->pair != nullptr){
		n_opp[1] = f->p_he->prev->pair->p_f->NormalVec;
	}else{
		n_opp[1] = n_T;
	}
	if(f->p_he->pair != nullptr){
		n_opp[2] = f->p_he->pair->p_f->NormalVec;
	}else{
		n_opp[2] = n_T;
	}

	for(int i = 0; i < 3; i++){
		n[i] = 0.5 * (n_T + n_opp[i]);
		e[i] = v[(i+2)%3] - v[(i+1)%3];
	}
	for(int i = 0; i < 3; i++){
		II_entry(i) = 2.0 * (n[(i+1)%3] - n[(i+2)%3]).dot(e[i]);
	}

	double L, M, N;
	L = II_entry(0);
	M = -0.5 * (II_entry(0) + II_entry(1) - II_entry(2));
	N = II_entry(1);
	Eigen::Matrix2d result;
	result << L, M, M, N;
	f->SFMat = result;
}

void HalfEdgeMesh::CalcFaceAreas(){
	int n_f = (int)p_f.size();
	for(int i = 0; i < n_f; i++){
		CalcFaceArea(p_f[i]);
	}
}
void HalfEdgeMesh::CalcFaceNormals(){
	int n_f = (int)p_f.size();
	for(int i = 0; i < n_f; i++){
		CalcFaceNormal(p_f[i]);
	}
}
void HalfEdgeMesh::CalcDAngs(){
	int n_e = (int)p_e.size();
	for(int i = 0; i < n_e; i++){
		CalcDAng(p_e[i]);
	}
}
void HalfEdgeMesh::CalcIAngs(){
	int n_f = (int)p_f.size();
	for(int i = 0; i < n_f; i++){
		HalfEdge* he = p_f[i]->p_he;
		do{
			CalcIAng(he);	
			he = he->next;
		}while(he != p_f[i]->p_he);
	}
}
void HalfEdgeMesh::CalcCotans(){
	int n_f = (int)p_f.size();
	for(int i = 0; i < n_f; i++){
		HalfEdge* he = p_f[i]->p_he;
		do{
			CalcCotan(he);
			he = he->next;
		}while(he != p_f[i]->p_he);
	}
}
void HalfEdgeMesh::CalcVertexMixedAreas(){
	int n_v = (int)p_v.size();
	for(int i = 0; i < n_v; i++){
		CalcVertexMixedArea(p_v[i]);
	}
}
void HalfEdgeMesh::CalcFFMats(){
	int n_f = (int)p_f.size();
	for(int i = 0; i < n_f; i++){
		CalcFFMat(p_f[i]);
	}
}
void HalfEdgeMesh::CalcSFMats(){
	int n_f = (int)p_f.size();
	for(int i = 0; i < n_f; i++){
		CalcSFMat(p_f[i]);
	}
}
void HalfEdgeMesh::CalcCotMat(Eigen::SparseMatrix<double>& cotmat){
	int n_v = (int)p_v.size();
	cotmat.resize(n_v, n_v);
	std::vector< Eigen::Triplet<double> > triplets;

	for(int i = 0; i < n_v; i++){
		HalfEdge* hei = p_v[i]->p_he;
		double value_ii = 0.0;
		do{
			int j_global = hei->next->p_v->Id;
			double wij = hei->Cotan;
			if(hei->pair != nullptr){
				wij += hei->pair->Cotan;
			}
			wij *= 0.5;
			triplets.emplace_back(i, j_global, wij);
			value_ii += wij;

			HalfEdge* hei_tmp = hei;
			hei = hei->prev->pair;
			if(hei == nullptr){
				hei_tmp = hei_tmp->prev;
				j_global = hei_tmp->p_v->Id;
				wij = hei_tmp->Cotan;
				wij *= 0.5;
				triplets.emplace_back(i, j_global, wij);
				value_ii += wij;
			}
		}while(hei != p_v[i]->p_he && hei != nullptr);
		triplets.emplace_back(i, i, -value_ii);
	}
	cotmat.setFromTriplets(triplets.begin(), triplets.end());
}	
void HalfEdgeMesh::CalcCotMat(Eigen::SparseMatrix<double>& cotmat, Eigen::SparseMatrix<double>& cotmat_int, Eigen::SparseMatrix<double>& cotmat_bnd){
	int n_v = (int)p_v.size();
	cotmat.resize(n_v, n_v);
	std::vector< Eigen::Triplet<double> > triplets;

	for(int i = 0; i < n_v; i++){
		HalfEdge* hei = p_v[i]->p_he;
		double value_ii = 0.0;
		do{
			int j_global = hei->next->p_v->Id;
			double wij = hei->Cotan;
			if(hei->pair != nullptr){
				wij += hei->pair->Cotan;
			}
			wij *= 0.5;
			triplets.emplace_back(i, j_global, wij);
			value_ii += wij;

			HalfEdge* hei_tmp = hei;
			hei = hei->prev->pair;
			if(hei == nullptr){
				hei_tmp = hei_tmp->prev;
				j_global = hei_tmp->p_v->Id;
				wij = hei_tmp->Cotan;
				wij *= 0.5;
				triplets.emplace_back(i, j_global, wij);
				value_ii += wij;
			}
		}while(hei != p_v[i]->p_he && hei != nullptr);
		triplets.emplace_back(i, i, -value_ii);
	}
	cotmat.setFromTriplets(triplets.begin(), triplets.end());

	int n_bnd = (int)p_bv.size();
	int n_int = n_v - n_bnd;
	cotmat_int.resize(n_int, n_int);
	cotmat_bnd.resize(n_int, n_bnd);
	std::vector< Eigen::Triplet<double> > triplets_int;
	std::vector< Eigen::Triplet<double> > triplets_bnd;
	std::vector<int> bv;
	for(int i = 0; i < n_bnd; i++){
		bv.push_back(p_bv[i]->Id);
	}
	std::vector<int> bv_sort;
	bv_sort = bv;
	std::sort(bv_sort.begin(), bv_sort.end());

	int ii = 0;
	int jj = 0;
	int bi = 0;
	int bj = 0;
	for(int i = 0; i < n_v; i++){
		if(i != bv_sort[bi]){
			for(int j = 0; j < n_v; j++){
				if(j != bv_sort[bj]){
					triplets_int.emplace_back(ii, jj, cotmat.coeff(i, j));
					jj++;
				}else{
					triplets_bnd.emplace_back(ii, bj, cotmat.coeff(i, j));
					bj++;
				}
			}
			ii++;
			jj = 0;
			bj = 0;
		}else{
			bi++;
		}
	}
	cotmat_int.setFromTriplets(triplets_int.begin(), triplets_int.end());
	cotmat_bnd.setFromTriplets(triplets_bnd.begin(), triplets_bnd.end());
}

void HalfEdgeMesh::MapBoundaryToCircle(HalfEdgeMesh* hem){
	int n_v = (int)hem->p_v.size();
	int n_bnd = (int)hem->p_bv.size();
	std::vector<double> len(n_bnd);
	double sum_len = 0.0;

	std::vector<Vertex*> p_bv_loop;
	p_bv_loop.push_back(hem->p_bv[0]);
	for(int i = 0; i < n_bnd - 1; i++){
		Vertex* tmp = p_bv_loop[i]->p_he->next->p_v; // 境界辺のp_heは境界辺に沿ったハーフエッジにしている
		p_bv_loop.push_back(tmp);
	}

	// Boundaryの最初の点の位置の角度(θ_0)は、もとのメッシュ(の極座標でのφ)に合わせる
	double theta0 = atan2(hem->p_bv[0]->Pos(1), hem->p_bv[0]->Pos(0));
	len[0] = 0.0;
	for(int i = 1; i < n_bnd; i++){
		len[i] = len[i-1] + (p_bv_loop[i-1]->Pos - p_bv_loop[i]->Pos).norm();
	}
	double total_len = len[n_bnd-1] + (p_bv_loop[0]->Pos - p_bv_loop[n_bnd-1]->Pos).norm();
	for(int i = 0; i < n_bnd; i++){
		double frac = 2.0 * M_PI * len[i] / total_len; 
		hem->p_v[p_bv_loop[i]->Id]->Pos(0) = cos(frac + theta0);
		hem->p_v[p_bv_loop[i]->Id]->Pos(1) = sin(frac + theta0);
	}
}
void HalfEdgeMesh::HarmonicParam(HalfEdgeMesh* hem){
	CopyHalfEdgeMesh(hem);
	MapBoundaryToCircle(hem);
	int n_bnd = (int)p_bv.size();
	int n_v = (int)p_v.size();
	int n_int = n_v - n_bnd;

	Eigen::SparseMatrix<double> cotmat, cotmat_int, cotmat_bnd;
	CalcCotans();
	CalcCotMat(cotmat, cotmat_int, cotmat_bnd);

	for(int i = 0; i < n_v; i++){
		hem->p_v[i]->Pos(2) = 0.0;
	}

	// Calc RHS ------------------------------
	std::vector<int> bv;
	for(int i = 0; i < n_bnd; i++){
		bv.push_back(p_bv[i]->Id);
	}
	std::vector<int> bv_sort = bv;
	std::sort(bv_sort.begin(), bv_sort.end());

	Eigen::VectorXd u_bnd, v_bnd;
	u_bnd.resize(n_bnd);
	v_bnd.resize(n_bnd);
	int ii = 0;
	int bi = 0;
	for(int i = 0; i < n_v; i++){
		if(i == bv_sort[bi]){
			u_bnd(ii) = hem->p_v[i]->Pos(0);
			v_bnd(ii) = hem->p_v[i]->Pos(1);
			ii++;
			bi++;
		}
	}

	Eigen::VectorXd ubar, vbar;
	ubar.resize(n_int);
	vbar.resize(n_int);

	ubar = -cotmat_bnd * u_bnd;
	vbar = -cotmat_bnd * v_bnd;
	// ---------------------------------------

	Eigen::VectorXd u_int, v_int;
	u_int.resize(n_int);
	v_int.resize(n_int);

	// Solve ---------------------------------
	//Eigen::SparseLU< Eigen::SparseMatrix<double> > solver;
	Eigen::SimplicialLDLT< Eigen::SparseMatrix<double> > solver;
	solver.compute(cotmat_int);
	if(solver.info() != Eigen::Success){
		std::cout << "Decomposition failed" << std::endl;
		exit(1);
	}

#pragma omp parallel sections
{
#pragma omp sections
{
	u_int = solver.solve(ubar);
	if(solver.info() != Eigen::Success){
		std::cout << "Failed solver" << std::endl;
		exit(1);
	}
}
#pragma omp sections
{
	v_int = solver.solve(vbar);
	if(solver.info() != Eigen::Success){
		std::cout << "Failed solver" << std::endl;
		exit(1);
	}
}
}
	// ---------------------------------------
	
	ii = 0;
	bi = 0;
	for(int i = 0; i < n_v; i++){
		if(i != bv_sort[bi]){
			hem->p_v[i]->Pos(0) = u_int(ii);
			hem->p_v[i]->Pos(1) = v_int(ii);
			ii++;
		}else{
			bi++;
		}
	}
}
