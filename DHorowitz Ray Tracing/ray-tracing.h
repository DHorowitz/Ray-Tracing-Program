#include <list>
using namespace std;

class Color;
class Light;
class Figure;
class Plane;
class Sphere;

class Vec
{
	friend Vec operator*(double num, const Vec& v);

  public:
	  double x;
	  double y;
	  double z;
    Vec();
	double dot(const Vec& other) const;
    Vec(ifstream& ifs);
	Vec(double x, double y, double z);
	Vec operator+(const Vec& otherVec) const;
	Vec operator-(const Vec& otherVec) const;
	void normalize();
	double norm() const;

};

class Ray
{
	public:
		Vec point1;
		Vec point2;
		Ray();
		Ray(const Vec& point1, const Vec& point2);
		Vec direction() const;
		Vec position(double t) const;
};

class Color
{
	friend Color operator*(double num, const Color& c);

  public:
	  double red;
	  double green;
	  double blue;
	Color();
    Color(ifstream& ifs);
    Color(double r, double g, double b);
	Color operator+(const Color& c) const;
	Color operator-(const Color& c) const;
	Color operator*(const Color& c) const;
	double norm() const;
	double diff(const Color& c) const;

};

class Light
{
public:
	Light(ifstream& ifs);
	double getAtten(const Vec& intersection) const;
	Vec position;
	Color shading;
	double c0, c1, c2;

};


class Figure
{

public:
	Color ambient;
	Color diffuse;
	Color specular;
	Color reflectivity;
	Color transmissivity;
	double shininess;
	double indexOfRefraction;
	int rFlag, tFlag;
	Figure();
	void initFigure(ifstream& ifs);
	virtual double intersection(const Ray& ray, double minT, double maxT) const = 0;
	virtual pair<Vec, bool> normal(const Vec& vect, const Ray& ray)  const = 0;
	bool transmits() const;
	bool reflects() const;
};



class Plane : public Figure
{

  public:
	  Vec abcVector;
	  double dScalar;
	  Vec direction1;
	  Vec direction2;
    Plane(ifstream& ifs);
	virtual double intersection(const Ray& ray, double minT, double maxT) const;
	virtual pair<Vec, bool> normal(const Vec& i, const Ray& ray) const;
};

class Sphere : public Figure
{

  public:
	  Vec center;
	  double radius;
    Sphere(ifstream& ifs);
	virtual double intersection(const Ray& ray, double minT, double maxT) const;
	virtual pair<Vec, bool> normal(const Vec& i, const Ray& r) const;
};

Vec operator*(double num, const Vec& vector);

Color operator*(double num, const Color& col);

Color RT_trace(const Ray& ray, double depth);
