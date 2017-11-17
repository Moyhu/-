#include "correctlib.h"
#include "light_matrix.h"
#include <stdio.h>  
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define PI 3.14159265358979323846
//图像的像素直接提取  
#define        _I(img,x,y,n) ((unsigned char*)((img)->imageData))[(y)*(img)->width*(img)->channels + (x)*(img)->channels + (n)]
//亚像素级灰度值  
#define        _IF(image,x,y,n)    ( ((int)(x+1)-(x))*((int)(y+1)-(y))*_I((image), (int)(x),(int)(y), (int)(n))\
									+ ((int)(x+1)-(x))*((y)-(int)(y))*_I((image), (int)(x),(int)(y+1), (int)(n))\
									+ ((x)-(int)(x))*((int)(y+1)-(y))*_I((image), (int)(x+1),(int)(y), (int)(n))\
									+ ((x)-(int)(x))*((y)-(int)(y))*_I((image),(int)(x+1),(int)(y+1), (int)(n)) )//插值后的像素值(IN表示interpolation),x、y可以为小数  

/*
图像灰度化
采用最大值法
将彩色图像中的三分量亮度的最大值作为灰度图的灰度值。
f(i, j) = max(R(i, j), G(i, j), B(i, j))
*/
GrayImage *Image2Gray(ClImage *img)
{
	GrayImage *gray = (GrayImage*)malloc(sizeof(GrayImage));
	gray->height = img->height;
	gray->width = img->width;
	gray->grayData = (unsigned char*)malloc(gray->height*gray->width);
	int w, h;
	for (h = 0; h < gray->height; h++)
		for (w = 0; w < gray->width; w++)
			gray->grayData[h*gray->width + w] = 
				img->imageData[h*img->width * 3 + w] > img->imageData[h*img->width * 3 + w + 1]
					? img->imageData[h*img->width * 3 + w] > img->imageData[h*img->width * 3 + w + 2]
						? img->imageData[h*img->width * 3 + w] : img->imageData[h*img->width * 3 + w + 2]
					: img->imageData[h*img->width * 3 + w + 1] > img->imageData[h*img->width * 3 + w + 2]
						? img->imageData[h*img->width * 3 + w + 1] : img->imageData[h*img->width * 3 + w + 2];
	return gray;
}

/*
提取鱼眼图像有效轮廓参数
半径、圆心坐标
使用 变角度线扫描法
*/
int ObtainParam(CorrectParam *param, GrayImage *gray)
{
	double k;	//切线斜率
	int tangentPoint[6][2][2];	//12个切点坐标  [6对][正向/反向][x/y]
	bool isOk;	//是否找到切点
	int x, y, b, i;
	char grayPoint;
	for (i = 0; i < 6; i++)
	{
		k = sin((PI / 6) * i);
		y = 0;
		isOk = false;
		for (b = 0; b < gray->height; b++)
		{
			for (x = 0; x < gray->width; x++)
			{
				y = k * x + b;
				if (y >= gray->height || y < 0)
					break;
				grayPoint = gray->grayData[y * gray->width + x];
				if (grayPoint > 50)
				{
					tangentPoint[i][0][0] = x;
					tangentPoint[i][0][1] = y;
					isOk = true;
					break;
				}
			}
			if (isOk)
				break;
		}
		
		isOk = false;
		for (b = gray->height - 1; b >= 0; b--)
		{
			for (x = gray->width - 1; x >= 0; x--)
			{
				y = k * x + b;
				if (y >= gray->height || y < 0)
					break;
				grayPoint = gray->grayData[y * gray->width + x];
				if (grayPoint > 50)
				{
					tangentPoint[i][1][0] = x;
					tangentPoint[i][1][1] = y;
					isOk = true;
					break;
				}
			}
			if (isOk)
				break;
		}
	}

	//剔除无效切点
	double length[6], avglength;
	int n = 0;	//有效坐标对个数
	for (i = 0; i < 6; i++)
		length[i] = sqrt(pow((tangentPoint[i][1][1] - tangentPoint[i][0][1]), 2) + pow((tangentPoint[i][1][0] - tangentPoint[i][0][0]), 2));
	
	avglength = (length[0] + length[1] + length[2] + length[3] + length[4] + length[5]) / 6;

	for (i = 0; i < 6; i++)
	{
		if (length[i] <= avglength)
		{
			tangentPoint[n][0][0] = tangentPoint[i][0][0];
			tangentPoint[n][0][1] = tangentPoint[i][0][1];
			tangentPoint[n][1][0] = tangentPoint[i][1][0];
			tangentPoint[n][1][1] = tangentPoint[i][1][1];
			n++;
		}
	}

	//使用kasa圆拟合法

	Mat mA, mAinv, mB, mP;
	float valA[36], valB[12], valP[3];
	for (i = 0; i <= n; i++)
	{
		valA[i * 6 + 0] = tangentPoint[i][0][0];
		valA[i * 6 + 1] = tangentPoint[i][0][1];
		valA[i * 6 + 2] = 1;
		valA[i * 6 + 3] = tangentPoint[i][1][0];
		valA[i * 6 + 4] = tangentPoint[i][1][1];
		valA[i * 6 + 5] = 1;

		valB[i * 2 + 0] = pow(tangentPoint[i][0][0], 2) + pow(tangentPoint[i][0][1], 2);
		valB[i * 2 + 1] = pow(tangentPoint[i][1][0], 2) + pow(tangentPoint[i][1][1], 2);
	}

	MatCreate(&mA, n, 3);
	MatCopy(&mA, &mAinv);
	MatCreate(&mB, n, 1);
	MatCreate(&mP, 3, 1);

	MatSetVal(&mA, valA);
	MatSetVal(&mB, valB);

	MatInv(&mA, &mAinv);
	MatMul(&mAinv, &mB, &mP);

	param->cx = mP.element[0][0] / 2;
	param->cy = mP.element[1][0] / 2;
	param->radius = sqrt((pow(mP.element[0][0], 2.0) + pow(mP.element[1][0], 2.0)) / 4 + mP.element[2][0]);

	return 1;
}

int EasyObtainParam(CorrectParam *param, ClImage *img)
{
	param->cx = img->width / 2;
	param->cy = img->height / 2;
	param->radius = param->cx < param->cy ? param->cx : param->cy;
	return 1;
}

/*
图像坐标转直角坐标
x0,y0:直角坐标原点在图像坐标系中的位置
u,v:图像坐标
x,y:直角坐标
*/
void ImageCoord2RectangularCoord(int x0, int y0, int u, int v, int *x, int *y)
{
	*x = u - x0;
	*y = -1 * (v - y0);
}

/*
矫正前后坐标映射关系
R:鱼眼图像半径
f:鱼眼图像焦距
x,y:源图像坐标（直角坐标系）
u,v:目标图像
*/
void CorrectMap(int R, double f, int x, int y, int *u, int *v)
{
	double temp1 = R * tan(sqrt(x * x + y * y) / f);
	double temp2 = atan(y / x);
	*u = temp1 * cos(temp2);
	*v = temp1 * sin(temp2);
}

/*
显示有效区域
*/
ClImage *ShowValidArea(ClImage *mImg, CorrectParam *param)
{
	ClImage *cImg = (ClImage*)malloc(sizeof(ClImage));
	cImg->height = 2 * param->radius;
	cImg->width = 2 * param->radius;
	cImg->channels = 3;
	cImg->imageData = (unsigned char*)malloc(sizeof(unsigned char) * cImg->height * cImg->width * cImg->channels);
	int w, h, u, v;
	for (h = param->cy - param->radius, v = 0; h < param->cy + param->radius; h++, v++)
	{
		if (h < 0 || h >= mImg->height)	//判断该点是否在源图像内
			continue;
		for (w = param->cx - param->radius, u = 0; w < param->cx + param->radius; w++, u++)
		{
			if (w < 0 || w >= mImg->width)	//判断该点是否在源图像内
				continue;
			memcpy(&cImg->imageData[v*cImg->width*cImg->channels + u*cImg->channels], &mImg->imageData[h*mImg->width*mImg->channels + w*mImg->channels], cImg->channels);
		}
	}

	return cImg;
}

/*
矫正
*/
ClImage *Correct1(ClImage *mImg, CorrectParam *param, double f)
{
	ClImage *cImg = (ClImage*)malloc(sizeof(ClImage));
	cImg->height = 2 * param->radius;
	cImg->width = 2 * param->radius;
	cImg->channels = 3;
	cImg->imageData = (unsigned char*)malloc(sizeof(unsigned char) * cImg->height * cImg->width * cImg->channels);
	//memset(cImg->imageData, 0, sizeof(sizeof(unsigned char) * cImg->height * cImg->width * cImg->channels));
	int w, h, x, y, u, v;
	for (h = param->cy - param->radius; h < param->cy + param->radius; h++)
	{
		if (h < 0 || h >= mImg->height)	//判断该点是否在源图像内
			continue;
		for (w = param->cx - param->radius; w < param->cx + param->radius; w++)
		{
			if (w < 0 || w >= mImg->width)	//判断该点是否在源图像内
				continue;
			ImageCoord2RectangularCoord(param->cx, param->cy, w, h, &x, &y);
			if (x == 0)
				continue;
			CorrectMap(param->radius, f, x, y, &u, &v);
			//CorrectMap(param->radius, f, w + 1, h + 1, &u, &v);
			if (u < 0 || u >= cImg->width || v < 0 || v >= cImg->height)
				continue;
			memcpy(&cImg->imageData[v*cImg->width*cImg->channels + u*cImg->channels], &mImg->imageData[h*mImg->width*mImg->channels + w*mImg->channels], cImg->channels);
		}
	}

	return cImg;
}

ClImage *Correct2(ClImage *bmpIn, CorrectParam *param)
{
	ClImage *bmpOut = (ClImage*)malloc(sizeof(ClImage));
	bmpOut->height = 4 * param->radius;
	bmpOut->width = 4 * param->radius;
	bmpOut->channels = 3;
	bmpOut->imageData = (unsigned char*)malloc(sizeof(unsigned char) * bmpOut->height * bmpOut->width * bmpOut->channels);

	double lens = param->radius * 2 / PI;
	//assume the fov is 180
	//R = f*theta

	int x, y, src_x, src_y, i, j;
	double r, theta, k;
	for (j = 0; j < bmpOut->height; j++)
	{
		for (i = 0; i < bmpOut->width; i++)
		{

			//offset to center
			ImageCoord2RectangularCoord(param->cx, param->cx, i, j, &x, &y);
			y = -y;
			r = sqrt(x*x + y*y);
			theta = atan(r / param->radius);
			if (theta<0.00001)
				k = 1;
			else
				k = lens*theta / r;

			src_x = (int)(x*k + param->cx);
			src_y = (int)(y*k + param->cy);
			if (src_x < 0 || src_x >= bmpIn->width || src_y < 0 || src_y >= bmpIn->height)
				continue;

			memcpy(&bmpOut->imageData[j*bmpOut->width*bmpOut->channels + i*bmpOut->channels], &bmpIn->imageData[src_y*bmpIn->width*bmpIn->channels + src_x*bmpIn->channels], bmpIn->channels);
			//printf("%d ", i);
		}
	}
	return bmpOut;
}


ClImage *Correct3(ClImage *bmpIn, CorrectParam *param)
{
	ClImage *bmpOut = (ClImage*)malloc(sizeof(ClImage));
	bmpOut->height = 2 * param->radius; //bmpIn->height + 100;
	bmpOut->width = 2 * param->radius; //bmpIn->width + 100;
	bmpOut->channels = 3;
	bmpOut->imageData = (unsigned char*)malloc(sizeof(unsigned char) * bmpOut->height * bmpOut->width * bmpOut->channels);
	int w, h, u, v;
	double x, y, xx, yy, xxx, yyy, r2, r4;
	for (h = 0; h < bmpOut->height; h++)
	{
		for (w = 0; w < bmpOut->width; w++)
		{
			x = w;// -50;
			y = h;// -50;
			xx = (x - param->cx) / param->fx;
			yy = (y - param->cy) / param->fy;
			r2 = pow(xx, 2) + pow(yy, 2);
			r4 = pow(r2, 2);
			xxx = xx*(1 + param->k1*r2 + param->k2*r4) + 2 * param->p1*xx*yy + param->p2*(r2 + 2 * xx*xx);
			yyy = yy*(1 + param->k1*r2 + param->k2*r4) + 2 * param->p2*xx*yy + param->p1*(r2 + 2 * yy*yy);
			u = xxx*param->fx + param->cx;
			v = yyy*param->fy + param->cy;
			if (u > 0 && u < bmpIn->width && v > 0 && v < bmpIn->height)
			{
				//memcpy(&bmpOut->imageData[h*bmpOut->width*bmpOut->channels + w*bmpOut->channels], &bmpIn->imageData[v*bmpIn->width*bmpIn->channels + u*bmpIn->channels], bmpIn->channels);
				//bmpOut->imageData[h*bmpOut->width*bmpOut->channels + w*bmpOut->channels] = bmpIn->imageData[v*bmpIn->width*bmpIn->channels + u*bmpIn->channels];
				//bmpOut->imageData[h*bmpOut->width*bmpOut->channels + w*bmpOut->channels + 1] = bmpIn->imageData[v*bmpIn->width*bmpIn->channels + u*bmpIn->channels + 1];
				//bmpOut->imageData[h*bmpOut->width*bmpOut->channels + w*bmpOut->channels + 2] = bmpIn->imageData[v*bmpIn->width*bmpIn->channels + u*bmpIn->channels + 2];
				_I(bmpOut, w, h, 0) = (int)_IF(bmpIn, u, v, 0);
				_I(bmpOut, w, h, 1) = (int)_IF(bmpIn, u, v, 1);
				_I(bmpOut, w, h, 2) = (int)_IF(bmpIn, u, v, 2);
			}
			else
			{
				_I(bmpOut, w, h, 0) = 0;
				_I(bmpOut, w, h, 1) = 0;
				_I(bmpOut, w, h, 2) = 0;
			}
		}
	}

	return bmpOut;
}