#include <cstdlib>
#include <cmath>
#include <cfloat>
#include <fstream>
#include <cassert>
#include <list>
#include <iostream>
#include <vector>

using namespace std;

#include "ray-tracing.h"

list<Figure*> shapeList;
list<Light*> lightList;

double epsilon = 0.00000000000001;
double maxT = 1000000;

double cameraX, cameraY, cameraZ;
int horizontalResolution, verticalResolution;
double zCoor;
double minX, maxX;
double minY, maxY;
double pixelWidth;
double pixelHeight;
double pixelStartX;
double pixelStartY;

int maxDepth;
Color backgroundColor;
Color ambient;

Color _DEFAULT_ = Color(0.0, 0.0, 0.0);

vector<vector<Color>> image;

int min(int a, int b) { return (a <= b) ? a : b; }

Vec::Vec() : x(0.0),y(0.0),z(0.0) {}

Vec::Vec(ifstream& ifs)
{
  ifs >> x >> y >> z;
}

Vec::Vec(double x, double y, double z) : x(x), y(y), z(z) {}

double Vec::dot(const Vec& other) const
{
	// sum of product of components
	double sumX = x * (other.x);
	double sumY = y * (other.y);
	double sumZ = z * (other.z);
	return sumX + sumY + sumZ;
}

Vec Vec::operator+(const Vec& other) const
{
	// vector with sum of components of two vectors
	double x2 = x + other.x;
	double y2 = y + other.y;
	double z2 = z + other.z;
	return Vec(x2, y2, z2);
}

Vec Vec::operator-(const Vec& other) const
{
	// vector with difference of components of two vectors
	double x3 = x - other.x;
	double y3 = y - other.y;
	double z3 = z - other.z;
	return Vec(x3, y3, z3);
}

Vec operator*(double num, const Vec& vector)
{
	Vec newVec(num*vector.x, num*vector.y, num*vector.z);
	return newVec;
}

Ray::Ray() : point1(Vec()), point2(Vec()) {};

Ray::Ray(const Vec& p1, const Vec& p2) : point1(p1), point2(p2) {}

Vec Ray::position(double t) const
{
	Vec diff = point2 - point1;
	return (point1 + (t * diff));
}

Vec Ray::direction() const { return point2 - point1; }



Color::Color(): red(0.0), green(0.0), blue(0.0) {};

Color::Color(ifstream& ifs)
{
  ifs >> red >> green >> blue;
}

Color::Color(double r, double g, double b) : red(r), green(g), blue(b)
{}

double Color::norm() const
{
	return fabs(red) + fabs(green) + fabs(blue);
}

double Color::diff(const Color& other) const
{
	return operator-(other).norm();
}

Color Color::operator+ (const Color& c) const
{
	return Color(red + c.red, green + c.green, blue + c.blue);
}

Color Color::operator-(const Color& c) const
{
	return Color(red - c.red, green - c.green, blue - c.blue);
}

Color Color::operator*(const Color& c) const
{
	return Color(red*c.red, green*c.green, blue*c.blue);
}

Color operator*(double num, const Color& col)
{
	Color newColor(num*col.red, num*col.green, num*col.blue);
	return newColor;
}

void Vec::normalize()
{
	// Divide values by length
	double distance = norm();
	x = x / distance;
	y = y / distance;
	z = z / distance;
}

double Vec::norm() const
{
	// sqrt(sum components squared)
	return sqrt(x*x + y * y + z * z);
}

Figure::Figure(){}

void Figure::initFigure(ifstream& ifs)
{
 ambient = Color(ifs);
 diffuse = Color(ifs);
 specular = Color(ifs);
 reflectivity = Color(ifs);
 transmissivity = Color(ifs);
 ifs >> shininess >> indexOfRefraction >> rFlag >> tFlag;
}

bool Figure::transmits() const { return tFlag == 1; }

bool Figure::reflects() const { return rFlag == 1; }

Light::Light(ifstream& ifs) : position(ifs), shading(ifs)
{
  ifs >> c0 >> c1 >> c2;
}


double Light::getAtten(const Vec& intersect) const
{
	if (c1 == 0.0 && c2 == 0.0) return 1.0 / c0;
	Vec separation = intersect - position;
	double distance = separation.norm();
	return 1.0 / (c0 + c1 * distance + c2 * distance*distance);
}

Sphere::Sphere(ifstream& ifs) : center(ifs)
{
  ifs >> radius;
  initFigure(ifs);
}

double Sphere::intersection(const Ray& r, double minT, double maxT) const
{
	Vec vec1 = r.point2 - r.point1;
	Vec vec2 = r.point1 - center;
	double a = vec1.dot(vec1);
	double b = 2.0*(vec1.dot(vec2));
	double c = vec2.dot(vec2) - (radius*radius);
	double twoA = 2.0*a;
	double minusB = -b;
	double disc = (b*b) - (4.0*a*c);
	if (disc >= 0.0)
	{
		double rootDisc = sqrt(disc);
		double t2 = (minusB - rootDisc) / twoA;
		if ((t2 < maxT) && (t2 >= minT)) return t2;
		else
		{
			double t1 = (minusB + rootDisc) / twoA;
			if ((t1 < maxT) && (t1 >= minT)) return t1;
			else return maxT + 1.0;
		}
	}
	else return maxT + 1.0;
}

pair<Vec, bool> Sphere::normal(const Vec& vec, const Ray& ray) const
{
	Vec d = ray.direction();
	Vec n = vec - center;
	if (d.dot(n) < 0.0) return pair<Vec, bool>((1.0 / radius)*n, true);
	else return pair<Vec, bool>((-1.0 / radius)*n, false);
}

Plane::Plane(ifstream& ifs) : abcVector(ifs)
{
  ifs >> dScalar;
  initFigure(ifs);
  direction1 = Vec(ifs);
  direction2 = Vec(ifs);
}

double Plane::intersection(const Ray& ray, double minT, double maxT) const
{
	Vec d = ray.direction();
	Vec n = abcVector;
	double dotProduct = d.dot(n);
	if (dotProduct > -epsilon && dotProduct < epsilon) return maxT + 1.0;
	else {
		double t = (dScalar - ray.point2.dot(n)) / dotProduct;
		if (t >= minT && t <= maxT) return t;
		else return maxT + 1.0;
	}
}

pair<Vec, bool> Plane::normal(const Vec&, const Ray& ray) const
{
	Vec d = ray.direction();
	Vec n = abcVector;
	n.normalize();
	if (d.dot(n) < 0.0) return pair<Vec, bool>(n, true);
	else return pair<Vec, bool>(-1.0*n, false);
}

void parseSceneFile(char* sceneName)
{
  double bgr, bgg, bgb;
  double ar, ag, ab;
  ifstream ifs;
  assert (sceneName != 0);
  ifs.open(sceneName);
  ifs >> cameraX;
  ifs >> cameraY;
  ifs >> cameraZ;
  ifs >> zCoor;
  ifs >> minX >> maxX;
  ifs >> minY >> maxY;
  ifs >> bgr >> bgg >> bgb;
  backgroundColor = Color(bgr,bgg,bgb);
  ifs >> ar >> ag >> ab;
  ambient = Color(ar,ag,ab);
  ifs >> maxDepth;
  ifs >> horizontalResolution >> verticalResolution;
  pixelWidth = (maxX - minX) / (double)horizontalResolution;
  pixelHeight = (maxY - minY) / (double)verticalResolution;
  int numLights, numSpheres, numPlanes;
  ifs >> numLights;
  for (int i=0; i<numLights; ++i) lightList.push_front(new Light(ifs));
  ifs >> numSpheres;
  for (int i=0; i<numSpheres; ++i) shapeList.push_front(new Sphere(ifs));
  ifs >> numPlanes;
  for (int i=0; i<numPlanes; ++i) shapeList.push_front(new Plane(ifs));
  ifs.close();
}

void initializeImage() 
{
	image.reserve(verticalResolution);
	for (int x = 0; x < verticalResolution; x++)
	{
		image.push_back(vector<Color>());
		image[x].reserve(horizontalResolution);
		for (int y = 0; y < horizontalResolution; y++)
		{
			image[x].push_back(_DEFAULT_);
		}
	}
}

Color RT_reflect(Figure* obj, const Ray& ray, const Vec& i,
	const Vec& normal, double depth)
{
	if (depth <= maxDepth && obj->reflects())
	{
		Vec viewDirection = -1.0 * ray.direction();
		viewDirection.normalize();
		double dotProduct = viewDirection.dot(normal);
		Vec reflectDirection = (2.0*dotProduct)*normal - viewDirection;
		Ray r_ray(i, i + reflectDirection);
		Color color = RT_trace(r_ray, depth + 1);
		return obj->reflectivity * color;
	}
	else return _DEFAULT_;
}

Color RT_transmit(Figure* obj, const Ray& ray, const Vec& i, 
	const Vec& normal, bool entering, double depth)
{
	if (depth <= maxDepth && obj->transmits())
	{
		Vec viewDirection = -1.0 * ray.direction();
		viewDirection.normalize();
		double dotProduct = viewDirection.dot(normal);
		double dotProduct2 = dotProduct * dotProduct;
		double indexRatio = entering ? 1.0 / obj->indexOfRefraction : obj->indexOfRefraction;
		double indexRatio2 = indexRatio * indexRatio;
		double temp1 = indexRatio * dotProduct;
		double temp2 = 1.0 - indexRatio2 + indexRatio2 * dotProduct2;
		if (temp2 >= 0.0)
		{
			double temp3 = temp1 - sqrt(temp2);
			Vec t = temp3 * normal - indexRatio * viewDirection;
			Ray t_ray(i, i + t);
			return obj->transmissivity * RT_trace(t_ray, depth + 1);
		}
		else return _DEFAULT_;
	}
	else return _DEFAULT_;
}

Color diffuseShade(Figure* obj, Light* light, double dotProduct)
{
	if (dotProduct > 0.0)
		return dotProduct * (light->shading*obj->diffuse);
	else return _DEFAULT_;
}

Color specularShade(Figure* obj, const Vec& normal,
	Light* light, const Vec& lightDirection, double dotProduct,
	const Ray& ray)
{
	Vec reflectDirection = (2.0*dotProduct)*normal - lightDirection;
	Vec viewDirection = -1.0 * ray.direction();
	viewDirection.normalize();
	double dotProduct2 = viewDirection.dot(reflectDirection);
	if (dotProduct2 > 0.0)
	{
		double shineFactor = pow(dotProduct2, obj->shininess);
		return shineFactor * (light->shading*obj->specular);
	}
	return _DEFAULT_;
}

pair<double, Figure*> nearestIntersection(const Ray& r,
	double minT, double maxT,
	bool mayBeTransparent = true)
{
	double t = maxT + epsilon;
	Figure* f = NULL;
	for (list<Figure*>::iterator iter = shapeList.begin();
		iter != shapeList.end();
		++iter)
	{
		if (mayBeTransparent || !((*iter)->transmits()))
		{
			double newT = (*iter)->intersection(r, minT, maxT);
			if (newT < t && newT >= minT)
			{
				t = newT;
				f = *iter;
			}
		}
	}
	return pair<double, Figure*>(t, f);
}

Color RT_lights(Figure* obj, const Ray& ray, const Vec& i, const Vec& normal)
{

	Color color = Color();
	for (list<Light*>::iterator iter = lightList.begin();
		iter != lightList.end(); iter++)
	{
		Light* light = *iter;
		Ray L_Ray(i, light->position);
		Vec L = L_Ray.direction();
		L.normalize();
		if (L.dot(normal) > 0)
		{
			color = color + (light->getAtten(i)) * 
				diffuseShade(obj, light, L.dot(normal));
			color = color + (light->getAtten(i)) * 
				specularShade(obj, normal, light, L, L.dot(normal), ray);
		}
	}
	return color;
}

Color RT_shade(Figure* obj, const Ray& ray, const Vec& i, 
	const Vec& normal, bool entering, double depth)
{
	Color color = ambient * obj->ambient;
	color = color + RT_lights(obj, ray, i, normal);
	if (depth < maxDepth)
	{
		color = color + RT_reflect(obj, ray, i, normal, depth);
		color = color + RT_transmit(obj, ray, i, normal, entering, depth);
	}
	return color;
}

Color RT_trace(const Ray& r, double depth) 
{
	pair<double, Figure*> intersection = nearestIntersection(r, epsilon, maxT);
	if (intersection.second == NULL)
		return backgroundColor;
	else
	{
		Figure* figure = intersection.second;
		Vec i = (r.point1 + (intersection.first * (r.point2 - r.point1)));
		pair<Vec, bool> normalInfo = figure->normal(i, r);
		Vec normal = normalInfo.first;
		bool entering = normalInfo.second;
		return RT_shade(figure, r, i, normal, entering, depth);
	}
}

pair<double, double> pixelCenter(int i, int j)
{
	double x = minX + (pixelWidth * (0.5 + i));
	double y = minY + (pixelHeight * (0.5 + j));
	return pair<double, double>(x, y);
}

void RT_algorithm()
{
	Vec origin(cameraX, cameraY, cameraZ);
	for (int x = 0; x < verticalResolution; x++)
	{
		for (int y = 0; y < horizontalResolution; y++)
		{
			pair<double, double> pixCent = pixelCenter(x, y);
			const double x1 = pixCent.first;
			const double y1 = pixCent.second;
			Vec p(x1, y1, zCoor);
			Ray R(origin, p);
			image[x][y] = RT_trace(R, 1);
		}
	}

}

void writeImageFile()
{
	ofstream img("sample.ppm");
	img << "P3" << endl;
	img << horizontalResolution << " " << verticalResolution << endl;
	img << "255" << endl;

	for (int x = 0; x < verticalResolution; x++)
	{
		for (int y = 0; y < horizontalResolution; y++)
		{
			img << min(image[x][y].red * 255, 255) << " " 
				<< min(image[x][y].green * 255, 255) << " " 
				<< min(image[x][y].blue * 255, 255) << endl;
		}
	}
	img.close();

}

int main(int, char *argv[])
{
    parseSceneFile(argv[1]);
	initializeImage();
	RT_algorithm();
	writeImageFile();
    return 0;
}
