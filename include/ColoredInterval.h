#ifndef COLORED_INTERVAL_H
#define COLORED_INTERVAL_H

#include <iostream>
#include <vector>
#include <algorithm>
#define _USE_MATH_DEFINES
#include <math.h>
#include <ctgmath>
#include <float.h>

class ColoredRange;
class ColoredInterval
{
protected:
	long double color_;
	long double left_;
	long double right_;
public:
	ColoredInterval(long double left, long double right, long double color);

	friend class ColoredRange;
	friend ColoredRange operator+ (ColoredRange l, const ColoredInterval& r);
	friend ColoredRange operator+ (ColoredInterval l, const ColoredRange& r);
	friend ColoredRange operator+ (ColoredInterval l, const ColoredInterval& r);
};

class ColoredRange
{
protected:
	std::vector<ColoredInterval> arr_;
public:
	ColoredRange(void);
	ColoredRange(const ColoredInterval &inter);

	long int NumOfIndices(void);
	long double Value (long int index);
	void Trim (long double left, long double right);
	void Print(std::ostream & str);

	friend ColoredRange operator+ (ColoredRange l, const ColoredRange& r);
	friend ColoredRange operator+ (ColoredRange l, const ColoredInterval& r);
	friend ColoredRange operator+ (ColoredInterval l, const ColoredRange& r);
	friend ColoredRange operator+ (ColoredInterval l, const ColoredInterval& r);
};

ColoredRange operator+ (ColoredRange l, const ColoredRange& r);
ColoredRange operator+ (ColoredRange l, const ColoredInterval& r);
ColoredRange operator+ (ColoredInterval l, const ColoredRange& r);
ColoredRange operator+ (ColoredInterval l, const ColoredInterval& r);

#endif
