#ifndef _HALFEDGEMESH_H_
#define _HALFEDGEMESH_H_

// EdgeCollapse, EdgeSplitで新たに定められる頂点位置は、
// ちゃんと立体に沿う位置にすべきなので、要修正

// EdgeCollapse, EdgeFlipはできない場合があるので、
// そういう場合はfalseを返すような感じで書くべき？

// 参考: https://mitani.cs.tsukuba.ac.jp/lecture/jikken/polygon_operation.pdf 

#include <vector>
#include <string>
#include <algorithm>
#include <Eigen/Dense>
#include <Eigen/Sparse>

class Vertex;
class HalfEdge;
class Face;
class Edge;

class Vertex{
public:
	int Id; // トポロジー処理の際に更新しないので、使用するときは改めてSetVertexIds()が必要
	Eigen::Vector3d Pos;
	HalfEdge* p_he;

	Vertex(double x, double y, double z){
		Pos(0) = x;
		Pos(1) = y;
		Pos(2) = z;
		p_he = nullptr;
	}
	Vertex(Eigen::Vector3d v){
		Pos = v;
		p_he = nullptr;
	}

	double MixedArea;
	Eigen::Vector3d NormalVec;
};

class HalfEdge{
public:
	int Id;
	Vertex* p_v; // half-edgeの始点
	Face* p_f;
	HalfEdge* pair;
	HalfEdge* next;
	HalfEdge* prev;

	Edge* p_e;

	HalfEdge(Vertex* v){
		p_v = v;
		if(v->p_he == nullptr){
			v->p_he = this;
		}
		pair = nullptr;
	}

	double IAng;
	double Cotan;
};

class Face{
public:
	int Id;
	HalfEdge* p_he;

	Face(HalfEdge* he){
		p_he = he;
	}

	double Area;
	Eigen::Vector3d NormalVec;
	Eigen::Matrix2d FFMat;
	Eigen::Matrix2d SFMat;
};

class Edge{
public:
	int Id;
	Vertex* p_v[2];
	HalfEdge* p_he_left;
	HalfEdge* p_he_right;
	int BoundaryFlag;

	Edge(HalfEdge* he){
		p_he_left = he;
		he->p_e = this;
		BoundaryFlag = 0;
		if(p_he_left->pair != nullptr){
			p_he_right = p_he_left->pair;
			p_he_left->pair->p_e = this;
		}else{
			p_he_right = nullptr;
			BoundaryFlag = 1;
		}
		p_v[0] = he->p_v;
		p_v[1] = he->next->p_v;
	}

	double DAng;
};


class HalfEdgeMesh{
public:
	std::vector<Vertex*> p_v;
	std::vector<Face*> p_f;
	std::vector<Edge*> p_e;

	std::vector<HalfEdge*> p_he;

	std::vector<Vertex*> p_bv;

private:
	void SetHalfEdgePair(HalfEdge*);
	void SetHalfEdgePair(HalfEdge* he0, HalfEdge* he1){
		he0->pair = he1;
		he1->pair = he0;
	}
	void DeleteVertex(Vertex* v){
		p_v.erase(std::remove(std::begin(p_v), std::end(p_v), v), std::cend(p_v));
		delete v;
	}
	void DeleteFace(Face* f){
		p_f.erase(std::remove(std::begin(p_f), std::end(p_f), f), std::cend(p_f));
		delete f->p_he->next;
		delete f->p_he->prev;
		delete f->p_he;
		delete f;
	}
	void DeleteEdge(Edge* e){
		p_e.erase(std::remove(std::begin(p_e), std::end(p_e), e), std::cend(p_e));
		delete e;
	}
public:
	void CopyHalfEdgeMesh(HalfEdgeMesh* hem); // hemは空のHalfEdgeMeshでないといけない
	void ReadOFF(std::string);
	void WriteOFF(std::string);
	void WriteVTK(std::string);
	void AddFace(Vertex*, Vertex*, Vertex*);
	void ConstructEdge();
	void SetVertexIds();
	void SetFaceIds();
	void SetEdgeIds();
	void SetHalfEdgeIds();

	void EdgeSplit(HalfEdge*);
	void EdgeFlip(HalfEdge*);
	void EdgeCollapse(HalfEdge*);

	double CalcAngle(const Eigen::Vector3d& v, const Eigen::Vector3d& w, const Eigen::Vector3d& axis){ // [-π, π]
		// 参考: https://www.cs.utexas.edu/users/evouga/uploads/4/5/6/8/45689883/turning.pdf
		double theta = 2.0 * atan2((v.cross(w).dot(axis) / axis.norm()), v.dot(w) + v.norm() * w.norm());
		return theta;
	}
	void CalcFaceArea(Face*);
	void CalcFaceNormal(Face*);
	void CalcDAng(Edge*); // NormalVecを使う
	void CalcIAng(HalfEdge*);
	void CalcCotan(HalfEdge*);
	void CalcVertexMixedArea(Vertex*); // AreaとIAngを使う
	void CalcFFMat(Face*);
	void CalcSFMat(Face*); // NormalVecを使う

	void CalcFaceAreas();
	void CalcFaceNormals();
	void CalcDAngs();
	void CalcIAngs();
	void CalcCotans();
	void CalcVertexMixedAreas();
	void CalcFFMats();
	void CalcSFMats();
	void CalcCotMat(Eigen::SparseMatrix<double>& cotmat);
	void CalcCotMat(Eigen::SparseMatrix<double>& cotmat, Eigen::SparseMatrix<double>& cotmat_int, Eigen::SparseMatrix<double>& cotmat_bnd);

	void MapBoundaryToCircle(HalfEdgeMesh* hem);
	void HarmonicParam(HalfEdgeMesh* hem);
};

#endif
