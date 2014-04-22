// prog.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"

class t
{
public:
	t(){}
	~t(){}

	void f(const int &i, int *y){
		x_ = i;
	}
	void set(int x, int* y){
		x_ = x;
		y_ = new int (*y);
	}
private:
	int *y_;
	int x_;
};

int _tmain(int argc, _TCHAR* argv[])
{
	t X;
	int *z = new int(10);
	X.set(5, z);
	for (int i = 0; i < 1; i++)
	{
		int *y = new int(99);
		X.set(*y, y);
		delete y;
	}
	for (int i = 1; i < 2; i++)
		int y = 100;
	
	X;
	return 0;
}

