#pragma once
#include <math_lib\math_lib.h>
#include <vector>
#include <map>
#include <assimp\scene.h>

struct Face
{
	SVec3f V1;
	SVec3f V2;
	SVec3f V3;
	
	bool operator==(const Face &f) const{
		if (V1 == f.V1 && V2 == f.V2 && V3 == f.V3)
			return true;
		if (V1 == f.V3 && V2 == f.V1 && V3 == f.V2)
			return true;
		if (V1 == f.V2 && V2 == f.V3 && V3 == f.V1)
			return true;
		return false;
	};
	
	bool operator!=(const Face &f) const{
		return !(*this == f);
	};
	
	Face& operator=(const Face& f){
		if (this == &f)
			return *this;

		V1 = f.V1;
		V2 = f.V2;
		V3 = f.V3;
		return *this;
	}
};


class Edge
{
public:
	Edge(void);
	~Edge(void);
	
	Edge(SVec3f V1, SVec3f V2);

	Edge(SVec3f V1, SVec3f V2, 
		unsigned int i);
	
	Edge(SVec3f V1, SVec3f V2, 
		Face* lFace, Face* rFace, 
		unsigned int i);

	bool contains(const SVec3f &V1, const SVec3f &V2);
	bool contains(const SVec3f &V1, const SVec3f &V2, bool &orientation);
	
	void inverse();
	
	void setV1edge(Edge* edge);
	void setV2edge(Edge* edge);

	Face* getLeftFace(){
		return leftFace_;
	}

	Face* getRightFace(){
		return rightFace_;
	}

	Edge* getV1edge(){
		return nextV1edge_;
	}
	
	Edge* getV2edge(){
		return nextV2edge_;
	}

	SVec3f getV1(){
		return v1_;
	}

	SVec3f getV2(){
		return v2_;
	}

	void setIndex(unsigned int i){
		index_ = i;
		isIndexSetted_ = false;
	}

	void setV1(SVec3f V1){
		v1_ = V1;
	}

	void setV2(SVec3f V2){
		v2_ = V2;
	}

	void setLeftFace(Face* lFace){
		if (lFace != NULL){
			leftFace_ = new Face(*lFace);
		}
		else
			leftFace_ = NULL;
		isLeftFaceSetted_ = true;
	}

	void setRightFace(Face* rFace){
		if (rFace != NULL){
		rightFace_ = new Face(*rFace);
		}
		else
			rightFace_ = NULL;
		isRightFaceSetted_ = true;
	}

	bool isIndexSetted(){
		return isIndexSetted_;
	}

	bool isLeftFaceSetted(){
		return isLeftFaceSetted_;
	}

	bool isRightFaceSetted(){
		return isRightFaceSetted_;
	}

	bool isV1edgeSetted(){
		return isV1edgeSetted_;
	}

	bool isV2edgeSetted(){
		return isV2edgeSetted_;
	}

private:
	SVec3f v1_;
	SVec3f v2_;
	Face* leftFace_;
	Face* rightFace_;
	unsigned int index_;
	Edge* nextV1edge_;
	Edge* nextV2edge_;

	bool isIndexSetted_;
	bool isLeftFaceSetted_;
	bool isRightFaceSetted_;
	bool isV1edgeSetted_;
	bool isV2edgeSetted_;
};


class Mesh
{
public:
	Mesh(void);
	~Mesh(void);
	
	bool loadMesh(aiMesh* mesh);

	void saveToVRML();

private:
	unsigned int verticesCount_;
	unsigned int facesCount_;
	unsigned int edgesCount_;

	std::vector<Face> faces_;
	std::vector<SVec3f> vertices_;
	std::map<SVec3f, Edge> edges_;

	void halfEdgeCollapse(Edge* edge);

	void loadVertices(aiMesh* mesh);
	void loadFaces(aiMesh* mesh);
	void loadEdges();
	
	void initializeEdge(Edge* edge);
	int leftFace(const SVec3f &V1, const SVec3f &V2);
	int rightFace(const SVec3f &V1, const SVec3f &V2);
	
	void nextV1edge(Edge *currEdge, Edge* edge);
	void nextV2edge(Edge *currEdge, Edge* edge);
};

