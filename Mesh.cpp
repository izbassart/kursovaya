#include "StdAfx.h"
#include "Mesh.h"

#include <algorithm>


Mesh::Mesh(void)
{
}


Mesh::~Mesh(void)
{
}

void Mesh::halfEdgeCollapse(Edge* edge){
	std::vector<SVec3f>edgeKeys_;
	std::vector<int>faceKeys_;
	
	Edge* e_ = edge->getV1edge();
	int f_ = leftFace(edge->getV1(), edge->getV2());

	while ((*e_) == (*edge)){
		SVec3f key_ = e_->getV1() / 2.0 + e_->getV2() / 2.0;
		edgeKeys_.push_back(key_);
		
		faceKeys_.push_back(f_);

		e_ = e_->getV1edge();
		f_ = leftFace(edge->getV1(), edge->getV2());
	}
	
	int i = 0;
	for ( ; i < edgeKeys_.size(); i++){
		SVec3f currKey = edgeKeys_.at(i);
		edges_.erase(currKey);
	}
	
}
bool Mesh::loadMesh(aiMesh* mesh){

	verticesCount_ = mesh->mNumVertices;
	facesCount_ = mesh->mNumFaces;

	loadVertices(mesh);
	loadFaces(mesh);
	loadEdges();
	return true;
}

void Mesh::loadEdges(){

	edgesCount_ = 0;
	unsigned int i = 0;

	for ( ; i < faces_.size(); i++){
		Face f_ = faces_.at(i);
		Edge* edge_ = new Edge(f_.V1, f_.V2);
		Edge* nextV1edge_ = new Edge();
		Edge* nextV2edge_ = new Edge();
	
		initializeEdge(edge_);
		nextV1edge(edge_, nextV1edge_);
		nextV2edge(edge_, nextV2edge_);

		edge_->setV1edge(nextV1edge_);
		edge_->setV2edge(nextV2edge_);

		SVec3f key = f_.V1 / 2.0 + f_.V2 / 2.0;
		edges_[key] = *edge_;

		edge_ = new Edge(f_.V2, f_.V3);
		nextV1edge_ = new Edge();
		nextV2edge_ = new Edge();
	
		initializeEdge(edge_);
		nextV1edge(edge_, nextV1edge_);
		nextV2edge(edge_, nextV2edge_);

		edge_->setV1edge(nextV1edge_);
		edge_->setV2edge(nextV2edge_);

		key = f_.V2 / 2.0 + f_.V3 / 2.0;
		edges_[key] = *edge_;

		edge_ = new Edge(f_.V3, f_.V1);
		nextV1edge_ = new Edge();
		nextV2edge_ = new Edge();
	
		initializeEdge(edge_);
		nextV1edge(edge_, nextV1edge_);
		nextV2edge(edge_, nextV2edge_);

		edge_->setV1edge(nextV1edge_);
		edge_->setV2edge(nextV2edge_);

		key = f_.V3 / 2.0 + f_.V1 / 2.0;
		edges_[key] = *edge_;
	}
}

void Mesh::initializeEdge(Edge* edge){

	int fNum_ = leftFace(edge->getV1(), edge->getV2());
	if (fNum_ == -1)
		edge->setLeftFace(NULL);
	else
		edge->setLeftFace(&faces_.at(fNum_));

	fNum_ = rightFace(edge->getV1(), edge->getV2());
	if (fNum_ == -1)
		edge->setRightFace(NULL);
	else
		edge->setRightFace(&faces_.at(fNum_));
}

int Mesh::leftFace(const SVec3f &V1, const SVec3f &V2){
	
	int i = 0;
	for ( ; i < faces_.size(); i++)
	{
		Face* currFace = new Face();
		*currFace = faces_.at(i);
		if (currFace->V1 == V1 && currFace->V2 == V2)
			return i;
		if (currFace->V1 == V2 && currFace->V3 == V1)
			return i;
		if (currFace->V2 == V1 && currFace->V3 == V2)
			return i;
	}
	return -1;
}

int Mesh::rightFace(const SVec3f &V1, const SVec3f &V2){

	int i = 0;
	for ( ; i < faces_.size(); i++)
	{
		Face* currFace = new Face();
		*currFace = faces_.at(i);
		if (currFace->V1 == V1 && currFace->V3 == V2)
			return i;
		if (currFace->V2 == V2 && currFace->V3 == V1)
			return i;
		if (currFace->V1 == V2 && currFace->V2 == V1)
			return i;
	}
	return -1;
}

void Mesh::nextV1edge(Edge *currEdge, Edge* edge){

	const Face* lFace = currEdge->getLeftFace();
	if (lFace != NULL){

		if (lFace->V1 == currEdge->getV1() && lFace->V2 == currEdge->getV2())
			edge->setV2(lFace->V3);
		else if (lFace->V1 == currEdge->getV2() && lFace->V3 == currEdge->getV1())
			edge->setV2(lFace->V2);
		else if (lFace->V2 == currEdge->getV1() && lFace->V3 == currEdge->getV2())
			edge->setV2(lFace->V1);
		edge->setV1(currEdge->getV1());
		
		int fNum_ = leftFace(edge->getV1(), edge->getV2());
		if (fNum_ == - 1)
			edge->setLeftFace(NULL);
		else
			edge->setLeftFace(&faces_.at(fNum_));
		
		fNum_ = rightFace(edge->getV1(), edge->getV2());
		if (fNum_ == -1)
			edge->setRightFace(NULL);
		else
			edge->setRightFace(&faces_.at(fNum_));
	}
	else
		edge = NULL;
}

void Mesh::nextV2edge(Edge* currEdge, Edge* edge)
{
	const Face* rFace = currEdge->getRightFace();
	if (rFace != NULL){
	
		if (rFace->V1 == currEdge->getV2() && rFace->V2 == currEdge->getV1())
			edge->setV1(rFace->V3);
		else if (rFace->V1 == currEdge->getV1() && rFace->V3 == currEdge->getV2())
			edge->setV1(rFace->V2);
		else if (rFace->V2 == currEdge->getV2() && rFace->V3 == currEdge->getV1())
			edge->setV1(rFace->V1);
		edge->setV2(currEdge->getV2());
		
		int fNum_ = leftFace(edge->getV1(), edge->getV2());
		if (fNum_ != -1)
			edge->setLeftFace(&faces_.at(fNum_));
		else
			edge->setLeftFace(NULL);

		fNum_ = rightFace(edge->getV1(), edge->getV2());
		if (fNum_ != -1)
			edge->setRightFace(&faces_.at(fNum_));
		else
			edge->setRightFace(NULL);
	}
	else
		edge = NULL;
}

void Mesh::loadVertices(aiMesh* mesh){

	for (int i = 0; mesh->mNumVertices; i++){
		
		SVec3f currVec;
		currVec.x = mesh->mVertices[i].x;
		currVec.y = mesh->mVertices[i].y;
		currVec.z = mesh->mVertices[i].z;
		vertices_.push_back(currVec);
	}
}

void Mesh::loadFaces(aiMesh* mesh){
	
	for (int i = 0; i < mesh->mNumFaces; i++){
		
		Face currFace;
		aiFace currAiFace = mesh->mFaces[i];
		int v1num = currAiFace.mIndices[0];
		int v2num = currAiFace.mIndices[1];
		int v3num = currAiFace.mIndices[2];
		currFace.V1.x = mesh->mVertices[v1num].x;
		currFace.V1.y = mesh->mVertices[v1num].y;
		currFace.V1.z = mesh->mVertices[v1num].z;
		currFace.V2.x = mesh->mVertices[v2num].x;
		currFace.V2.y = mesh->mVertices[v2num].y;
		currFace.V2.z = mesh->mVertices[v2num].z;
		currFace.V3.x = mesh->mVertices[v3num].x;
		currFace.V3.y = mesh->mVertices[v3num].y;
		currFace.V3.z = mesh->mVertices[v3num].z;
		faces_.push_back(currFace);
	}
}

void Mesh::saveToVRML()
{

}

Edge::~Edge(void)
{
	delete leftFace_;
	delete rightFace_;
}

Edge::Edge(void)
{
	isIndexSetted_ = false;
	isLeftFaceSetted_ = false;
	isRightFaceSetted_ = false;
	isV1edgeSetted_ = false;
	isV2edgeSetted_ = false;

	index_ = -1;
	nextV1edge_ = NULL;
	nextV2edge_ = NULL;
	leftFace_ = NULL;
	rightFace_ = NULL;
}

Edge::Edge(SVec3f V1, SVec3f V2)
{
	isIndexSetted_ = false;
	isLeftFaceSetted_ = false;
	isRightFaceSetted_ = false;
	isV1edgeSetted_ = false;
	isV2edgeSetted_ = false;

	index_ = -1;
	
	nextV1edge_ = NULL;
	nextV2edge_ = NULL;
	leftFace_ = NULL;
	rightFace_ = NULL;

	v1_ = V1;
	v2_ = V2;
}

Edge::Edge(SVec3f V1, SVec3f V2, unsigned int i)
{
	isIndexSetted_ = true;
	isLeftFaceSetted_ = false;
	isRightFaceSetted_ = false;
	isV1edgeSetted_ = false;
	isV2edgeSetted_ = false;

	index_ = i;
	nextV1edge_ = NULL;
	nextV2edge_ = NULL;
	leftFace_ = NULL;
	rightFace_ = NULL;

	v1_ = V1;
	v2_ = V2;
}

Edge::Edge(SVec3f V1, SVec3f V2, Face* lFace, Face* rFace, unsigned int i)
{
	v1_ = V1;
	v2_ = V2;

	leftFace_ = new Face();
	*leftFace_ = *lFace;
	rightFace_ = new Face();
	*leftFace_ = *lFace;

	index_ = i;

	isLeftFaceSetted_ = true;
	isRightFaceSetted_ = true;
	isIndexSetted_ = true;

	isV1edgeSetted_ = false;
	isV2edgeSetted_ = false;

	nextV1edge_ = NULL;
	nextV2edge_ = NULL;
}

bool Edge::contains(const SVec3f &V1, const SVec3f &V2)
{
	if (v1_ == V1 && v2_ == V2)
		return true;
	else if (v2_ == V1 && v1_ == V2)
		return true;
	return false;
}

bool Edge::contains(const SVec3f &V1, const SVec3f &V2, bool &orientation)
{
	if (v1_ == V1 && v2_ == V2)
	{
		orientation = true;
		return true;
	}
	else if (v2_ == V1 && v1_ == V2)
	{
		orientation = false;
		return true;
	}
	return false;
}

void Edge::inverse()
{
	SVec3f tempV_ = v1_;
	v1_ = v2_;
	v2_ = tempV_;

	if (isV1edgeSetted_ && isV2edgeSetted_){
		Edge* tempE_ = new Edge(nextV1edge_->getV1(), nextV2edge_->getV2());
		
		if (nextV1edge_->isLeftFaceSetted())
			tempE_->setLeftFace(nextV1edge_->getLeftFace());
		if (nextV1edge_->isRightFaceSetted())
			tempE_->setRightFace(nextV1edge_->getRightFace());
		if (nextV1edge_->isV1edgeSetted())
			tempE_->setV1edge(nextV1edge_->getV1edge());
		if (nextV1edge_->isV2edgeSetted())
			tempE_->setV2edge(nextV1edge_->getV2edge());

		nextV1edge_->setV1(nextV2edge_->getV1());
		nextV1edge_->setV2(nextV2edge_->getV2());

		if (nextV2edge_->isLeftFaceSetted())
			nextV1edge_->setLeftFace(nextV2edge_->getLeftFace());
		if (nextV2edge_->isRightFaceSetted())
			nextV1edge_->setRightFace(nextV2edge_->getRightFace());
		if (nextV2edge_->isV1edgeSetted())
			nextV1edge_->setV1edge(nextV2edge_->getV1edge());
		if (nextV2edge_->isV2edgeSetted())
			nextV1edge_->setV2edge(nextV2edge_->getV2edge());

		nextV2edge_->setV1(tempE_->getV1());
		nextV2edge_->setV2(tempE_->getV2());

		if (tempE_->isLeftFaceSetted())
			nextV2edge_->setLeftFace(tempE_->getLeftFace());
		if (tempE_->isRightFaceSetted())
			nextV2edge_->setRightFace(tempE_->getRightFace());
		if (tempE_->isV1edgeSetted())
			nextV2edge_->setV1edge(tempE_->getV1edge());
		if (tempE_->isV2edgeSetted())
			nextV2edge_->setV2edge(tempE_->getV2edge());

		delete tempE_;
	}

	if (isLeftFaceSetted_ && isRightFaceSetted_){
		Face tempF_ = *leftFace_;
		*leftFace_ = *rightFace_;
		*rightFace_ = tempF_;
	}
}

void Edge::setV1edge(Edge* edge){

	if (edge != NULL){

		nextV1edge_ = new Edge(edge->getV1(), edge->getV2());
		
		if (edge->isLeftFaceSetted()){
			nextV1edge_->setLeftFace(edge->getLeftFace());
		}
		
		if (edge->isRightFaceSetted()){
			nextV1edge_->setRightFace(edge->getRightFace());
		}
		
		if (edge->isV1edgeSetted()){
			nextV1edge_->setV1edge(edge->getV1edge());
		}
		
		if (edge->isV2edgeSetted()){
			nextV1edge_->setV2edge(edge->getV2edge());
		}
	}
	else
		nextV1edge_ = NULL;

	isV1edgeSetted_ = true;
}

void Edge::setV2edge(Edge* edge){
		
	if (edge != NULL){
		
		nextV2edge_ = new Edge(edge->getV1(), edge->getV2());
		
		if (edge->isLeftFaceSetted()){
			nextV2edge_->setLeftFace(edge->getLeftFace());
		}
		
		if (edge->isRightFaceSetted()){
			nextV2edge_->setRightFace(edge->getRightFace());
		}
		
		if (edge->isV1edgeSetted()){
			nextV2edge_->setV1edge(edge->getV1edge());
		}
		
		if (edge->isV2edgeSetted()){
			nextV2edge_->setV2edge(edge->getV2edge());
		}
	}
	else
		nextV2edge_ = NULL;
	
	isV2edgeSetted_ = true;
}

bool operator == (Edge e1, Edge e2){
	return e1.contains(e2.getV1(), e2.getV2());
}

bool operator != (Edge e1, Edge e2){
	return !(e1 == e2);
}