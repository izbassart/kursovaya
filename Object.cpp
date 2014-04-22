#include "StdAfx.h"
#include "Object.h"
#include <assimp\Importer.hpp>
#include <assimp\Scene.h>
#include <assimp\postprocess.h>

Object::Object(void)
{
}


Object::~Object(void)
{
}

bool Object::loadObject(const std::string& filePath)
{
	Assimp::Importer importer;
	const aiScene* scene = importer.ReadFile(filePath, aiProcess_Triangulate);
	if (!scene)
		return false;
	meshesCount = scene->mNumMeshes;
	for (unsigned int i = 0; i < scene->mNumMeshes; i++)
	{
		Mesh currentMesh;
		if (!currentMesh.loadMesh(scene->mMeshes[i]))
			return false;
		meshes.push_back(currentMesh);
	}
	return true;
}

void Object::saveToVRML(const char* fileName)
{
	for (unsigned int i = 0; i < meshes.size(); i++)
	{
		
	}
}
