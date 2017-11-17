//#pragma once

#ifndef CORRECT_H
#define CORRECT_H

#include "bmplib.h"

//灰度图像
typedef struct
{
	int width;
	int height;
	unsigned char* grayData;
}GrayImage;

/*//鱼眼图像参数
typedef struct
{
	double radius;	//圆半径
	double centerX;	//圆心坐标
	double centerY;
}ImageParam;*/

//矫正参数
typedef struct
{
	double radius;	//圆半径

	double fx;	//相机焦距横坐标
	double fy;	//相机焦距纵坐标
	double cx;	//图像基准点横坐标	//圆心坐标
	double cy;	//图像基准点纵坐标	//圆心坐标
	double k1;	//径向畸变参数1
	double k2;	//径向畸变参数2
	double p1;	//切向畸变参数1
	double p2;	//切向畸变参数2
}CorrectParam;

GrayImage *Image2Gray(ClImage *img);
int ObtainParam(CorrectParam *param, GrayImage *gray);
int EasyObtainParam(CorrectParam *param, ClImage *img);
void ImageCoord2RectangularCoord(int x0, int y0, int u, int v, int *x, int *y);
void CorrectMap(int R, double f, int x, int y, int *u, int *v);
ClImage *Correct1(ClImage *mImg, CorrectParam *param, double f);
ClImage *ShowValidArea(ClImage *mImg, CorrectParam *param);
ClImage *Correct2(ClImage *bmpIn, CorrectParam *param);
ClImage *Correct3(ClImage *bmpIn, CorrectParam *param);

#endif  