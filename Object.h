#pragma once
#include <vector>
#include "Mesh.h"
#include <math_lib\math_lib.h>


class Object
{
public:
	Object(void);
	~Object(void);

	bool loadObject(const std::string& filePath);

	void saveToVRML(const char* fileName);

private:
	std::vector<Mesh> meshes;
	int meshesCount;
};

