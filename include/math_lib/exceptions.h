/********************************************************************************************
Author: Vasiliy Sazonov
date: 13/03/2006
Description: Exceptioins 
********************************************************************************************/

#ifndef __HEXCEPTH__
#define __HEXCEPTH__

#include <stdexcept>
#include <sstream>
#include <stdio.h>

#ifdef _MSC_VER
#define VS_SPRINTF sprintf_s
#define VS_SSCANF  sscanfs_s
#define VS_PRINTF  printf_s
#define VS_FPRINTF  fprintf_s
#else
#define VS_SPRINTF sprintf
#define VS_SSCANF  sscanfs
#define VS_PRINTF  printf
#define VS_FPRINTF  fprintf
#endif

const size_t WrongFile      = 1;
const size_t IOError        = 2;
const size_t BAD_FILE_WRITE = 3;

typedef enum { NoIntersection, NoPlate, ZeroLength, UnknownError } geom_error;
typedef enum { NoIncrement, BadIndexing } grid_error;
enum { ErrorInNameParse, NumberIsExpected, DuplictionDefinition, NoParentDefined, UnexpectedEOF, UnknownType, UnknownParserError } ;

class SBadOutputFile :  public std :: exception
{
  const char* fname;
public:
  SBadOutputFile(const char* name): fname(name) { }
  virtual const char* what() const throw()
  {
    static char buf[200];
    VS_SPRINTF(buf, "Cannot write to file %s. Check permissons.", fname);
    return buf;
  }
};

class SParseCXY : public std :: exception
{
  std :: string fname;
  int type;
  int lineno;
public:
  SParseCXY(const char * fn, int itype, int ln) 
    : fname(fn), type(itype), lineno(ln) { }
  SParseCXY() : type(UnknownParserError) { }
  virtual const char* what() const throw()
  { 
    static char buf[1024];
    switch(type)
    {
      case ErrorInNameParse:
        VS_SPRINTF(buf, "ERROR: %s:%d: Invalid name format.\n", fname.c_str(), lineno);
        break;
      case NumberIsExpected:
        VS_SPRINTF(buf, "ERROR: %s:%d: Number is expected.\n", fname.c_str(), lineno);
        break;
      case DuplictionDefinition:
        VS_SPRINTF(buf, "ERROR: %s:%d: This name is already defined.\n", fname.c_str(), lineno);
        break;
      case NoParentDefined:
        VS_SPRINTF(buf, "ERROR: %s:%d: Parent name is not defined.\n", fname.c_str(), lineno);
        break;    
      case UnexpectedEOF:
        VS_SPRINTF(buf, "ERROR: %s:%d: Unexpected end of file.\n", fname.c_str(), lineno);
        break;
      case UnknownType:    
      case UnknownParserError:
      default: 
        VS_SPRINTF(buf, "ERROR: %s:%d: Uknown error in parser.\n", fname.c_str(), lineno);
    }
    return buf;
  }
};

class SGridException : public std :: exception
{
  grid_error error_type;
public:
  SGridException(grid_error t) : error_type(t) { }
  virtual const char* what() const throw()
  {
    switch(error_type)
    { 
      case NoIncrement: 
        return "No increment in grid";
      case BadIndexing:
        return "Indecies are out of range";
      default:     
        return "Unknown grid exception";
    }
  }
};

class SBadIndeciesException : public SGridException
{
  const char * name;
  int ind_i, ind_xs;
  
public:
  SBadIndeciesException(const char* n, int i, int xs) :
    SGridException(BadIndexing), name(n), ind_i(i), ind_xs(xs) { }
  const char * what() const throw()
  {
    static char buf[200];
    VS_SPRINTF(buf, "Bad indexing %s=%d, where valid index should be [0..%d]", name, ind_i, ind_xs);
    return buf;
  }
};


class SGeomException : public std :: exception
{
  geom_error error_type;
public:
  SGeomException(geom_error type) :   error_type(type) { }
  SGeomException() :   error_type(UnknownError) { }
  virtual const char* what() const throw()
  {
    switch(error_type)
    {
      case NoIntersection: return "Object are not intersected.";
      case NoPlate: return "Error in plate definition.";
      case ZeroLength: return "Vector has zero length. Possible devision by zero.";
      default:
	      return "Unknow error.";
    }
  }

  bool IsZeroLength() const { return error_type == ZeroLength; }

};

class SGeomNoIntersection : public SGeomException
{
public:
  SGeomNoIntersection() : SGeomException(NoIntersection) { }
};

class SGeomNoPlate : public SGeomException
{
public:
  SGeomNoPlate() : SGeomException(NoPlate) { }
};

class SGeomZeroLength : public SGeomException
{
public:
  SGeomZeroLength() : SGeomException(ZeroLength) { }
};


class SIOException : public std :: exception
{
 size_t num;
public:
	SIOException(size_t n) : num(n) { }
 virtual const char* what() const throw()
 {
   switch(num)
   {
    case IOError: return "Error in IO operation.";
    case WrongFile: 
    default:
	   return "Invalid file descriptor. Check io fail.";
   }
 }
};


const size_t DuplicateWalking = 0;
class SEdgeProcessException : public std :: exception
{
 size_t num;
public:
 SEdgeProcessException(size_t n) : num(n) { }
virtual const char* what() const throw()
{
   switch(num)
   {
    case DuplicateWalking: return "Error edge processing: You have already chaged this edge for 3 times. ";
    default:
		return "Error edge processing: Unknown error.";
   }
 }
};

class SInterpolationException : public std :: exception
{
 double val;
 double Vmin;
 double Vmax;

public:
 SInterpolationException(double d, double Imin, double Imax) : val(d), Vmin(Imin), Vmax(Imax) { }
virtual const char* what() const throw()
{
   static char buf[1024];
   VS_PRINTF(buf, "The value %g is out of range [%g, %g]\n", val, Vmin, Vmax);
   return buf;
 }
};

#endif

