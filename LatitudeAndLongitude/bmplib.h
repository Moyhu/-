#ifndef BMPLIB_H
#define BMPLIB_H

typedef struct  
{  
    //unsigned short    bfType;  
    unsigned long    bfSize;  
    unsigned short    bfReserved1;  
    unsigned short    bfReserved2;  
    unsigned long    bfOffBits;  
} ClBitMapFileHeader;  
  
typedef struct  
{  
    unsigned long  biSize;   
    long   biWidth;   
    long   biHeight;   
    unsigned short   biPlanes;   
    unsigned short   biBitCount;  
    unsigned long  biCompression;   
    unsigned long  biSizeImage;   
    long   biXPelsPerMeter;   
    long   biYPelsPerMeter;   
    unsigned long   biClrUsed;   
    unsigned long   biClrImportant;   
} ClBitMapInfoHeader;  
  
typedef struct   
{  
    unsigned char rgbBlue; //????????  
    unsigned char rgbGreen; //????????  
    unsigned char rgbRed; //????????  
    unsigned char rgbReserved; //???  
} ClRgbQuad;  
  
typedef struct  
{  
    int width;  
    int height;  
    int channels;  
    unsigned char* imageData;  
}ClImage;  
  
ClImage* clLoadImage(char* path);  
bool clSaveImage(char* path, ClImage* bmpImg);  
  
#endif  