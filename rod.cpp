#include <OpenImageIO/imageio.h>
#include <iostream>
#include <stdlib.h>
#include <cstdlib>
#include <GL/glut.h>
#include <fstream>
#include <math.h>
#include <vector>
#include <queue>
#include <ctime>
#include <unordered_set>

#define MAX_SIZE 250
#define paletteSize 8

using namespace std;
OIIO_NAMESPACE_USING

// Strings to store the paths of the images and the kernel
string imageFileName;
string kernelPath;

struct Pixel{ // defines a pixel structure
	unsigned char r,g,b,a;
};

class UniqueIdentifier {
public:

	int operator() (const Pixel &pixel) const {
		return (int) pixel.r*50 + (int) pixel.g*40 + (int) pixel.b*30 + (int) pixel.a*20;
	}
};

int randomization;

int colors[paletteSize][4];
bool isImageSet = false;
//
// Global variables and constants
//
const int DEFAULTWIDTH = 600;	// default window dimensions if no image
const int DEFAULTHEIGHT = 600;

int WinWidth, WinHeight;	// window width and height
int ImWidth, ImHeight;		// image width and height
int ImChannels;           // number of channels per image pixel

int VpWidth, VpHeight;		// viewport width and height
int Xoffset, Yoffset;     // viewport offset from lower left corner of window

Pixel **pixmap = NULL;	// the image pixmap used for OpenGL
Pixel **dupPixmap = NULL; // duplicate pixmap containing read-only values for initial operation
int pixformat; 			// the pixel format used to correctly  draw the image

// Vector to store the fixed-size palette of colors
vector<Pixel> palette;

//
//  Routines to cleanup the memory. Note that GLUT and OpenImageIO v1.6 lead to memory leaks that cannot be suppressed.
//
void destroy(){
 if (pixmap){
     delete pixmap[0];
	 delete pixmap;
  }
}

//
//  Routine to read an image file and store in a pixmap
//  returns the size of the image in pixels if correctly read, or 0 if failure
//
int readimage(string infilename){
  // Create the oiio file handler for the image, and open the file for reading the image.
  // Once open, the file spec will indicate the width, height and number of channels.
  ImageInput *infile = ImageInput::open(infilename);


  if(!infile){
    cerr << "Could not input image file " << infilename << ", error = " << geterror() << endl;
    return 0;
  }

  // Record image width, height and number of channels in global variables
  ImWidth = infile->spec().width;
  ImHeight = infile->spec().height;
  ImChannels = infile->spec().nchannels;


  // allocate temporary structure to read the image
  unsigned char tmp_pixels[ImWidth * ImHeight * ImChannels];

  // read the image into the tmp_pixels from the input file, flipping it upside down using negative y-stride,
  // since OpenGL pixmaps have the bottom scanline first, and
  // oiio expects the top scanline first in the image file.
  int scanlinesize = ImWidth * ImChannels * sizeof(unsigned char);
  if(!infile->read_image(TypeDesc::UINT8, &tmp_pixels[0] + (ImHeight - 1) * scanlinesize, AutoStride, -scanlinesize)){
    cerr << "Could not read image from " << infilename << ", error = " << geterror() << endl;
    ImageInput::destroy(infile);
    return 0;
  }

 // get rid of the old OpenGL pixmap and make a new one of the new size
  destroy();

 // allocate space for the Pixmap (contiguous approach)
  pixmap = new Pixel*[ImHeight];
	dupPixmap = new Pixel*[ImHeight];

	if(pixmap != NULL)
		pixmap[0] = new Pixel[ImWidth * ImHeight];

	// Similar operations performed on dupPixmap to preserve the original pixmap values
	// Operations to be performed on original pixmap, hence it is prone to
	// changes in its values throughout the execution of the program
	if(dupPixmap != NULL)
		dupPixmap[0] = new Pixel[ImWidth * ImHeight];

	for(int i = 1; i < ImHeight; i++)
		pixmap[i] = pixmap[i - 1] + ImWidth;

	for(int i = 1; i < ImHeight; i++)
		dupPixmap[i] = dupPixmap[i - 1] + ImWidth;

 //  assign the read pixels to the Pixmap
 int index;
  for(int row = 0; row < ImHeight; ++row) {
    for(int col = 0; col < ImWidth; ++col) {
		index = (row*ImWidth+col)*ImChannels;

		if (ImChannels==1){
			pixmap[row][col].r = tmp_pixels[index];
			dupPixmap[row][col].r = tmp_pixels[index];

			pixmap[row][col].g = tmp_pixels[index];
			dupPixmap[row][col].g = tmp_pixels[index];

			pixmap[row][col].b = tmp_pixels[index];
			dupPixmap[row][col].b = tmp_pixels[index];

			pixmap[row][col].a = 255;
			dupPixmap[row][col].a = 255;

		}
		else{
			pixmap[row][col].r = tmp_pixels[index];
			dupPixmap[row][col].r = tmp_pixels[index];

			pixmap[row][col].g = tmp_pixels[index+1];
			dupPixmap[row][col].g = tmp_pixels[index + 1];

			pixmap[row][col].b = tmp_pixels[index+2];
			dupPixmap[row][col].b = tmp_pixels[index + 2];

			if (ImChannels <4) {// no alpha value is present so set it to 255
				pixmap[row][col].a = 255;
				dupPixmap[row][col].a = 255;
			}
			else { // read the alpha value
				pixmap[row][col].a = tmp_pixels[index+3];
				dupPixmap[row][col].a = tmp_pixels[index+3];
			}
		}
    }
  }

  // close the image file after reading, and free up space for the oiio file handler
  infile->close();
  ImageInput::destroy(infile);

  // set the pixel format to GL_RGBA and fix the # channels to 4
  pixformat = GL_RGBA;
  ImChannels = 4;

  // return image size in pixels
  return ImWidth * ImHeight;
}


//
// Routine to display a pixmap in the current window
//
void displayimage(){
  // if the window is smaller than the image, scale it down, otherwise do not scale
  if(WinWidth < ImWidth  || WinHeight < ImHeight)
    glPixelZoom(float(VpWidth) / ImWidth, float(VpHeight) / ImHeight);
  else
    glPixelZoom(1.0, 1.0);

  // display starting at the lower lefthand corner of the viewport
  glRasterPos2i(0, 0);

  glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
  glDrawPixels(ImWidth, ImHeight, pixformat, GL_UNSIGNED_BYTE, pixmap[0]);
}

//
// Routine to write the current framebuffer to an image file
//
void writeimage(string outfilename){
  // make a pixmap that is the size of the window and grab OpenGL framebuffer into it
   unsigned char local_pixmap[WinWidth * WinHeight * ImChannels];
   glReadPixels(0, 0, WinWidth, WinHeight, pixformat, GL_UNSIGNED_BYTE, local_pixmap);

	 for(int i=0; i<ImHeight; ++i) {
		 for(int j=0; j<ImWidth; ++j) {
			 local_pixmap[(i*ImWidth + j)*ImChannels + 0] = pixmap[i][j].r;
			 local_pixmap[(i*ImWidth + j)*ImChannels + 1] = pixmap[i][j].g;
			 local_pixmap[(i*ImWidth + j)*ImChannels + 2] = pixmap[i][j].b;
			 local_pixmap[(i*ImWidth + j)*ImChannels + 3] = pixmap[i][j].a;
		 }
	 }

  // create the oiio file handler for the image
  ImageOutput *outfile = ImageOutput::create(outfilename);
  if(!outfile){
    cerr << "Could not create output image for " << outfilename << ", error = " << geterror() << endl;
    return;
  }

  // Open a file for writing the image. The file header will indicate an image of
  // width WinWidth, height WinHeight, and ImChannels channels per pixel.
  // All channels will be of type unsigned char
  ImageSpec spec(WinWidth, WinHeight, ImChannels, TypeDesc::UINT8);
  if(!outfile->open(outfilename, spec)){
    cerr << "Could not open " << outfilename << ", error = " << geterror() << endl;
    ImageOutput::destroy(outfile);
    return;
  }

  // Write the image to the file. All channel values in the pixmap are taken to be
  // unsigned chars. While writing, flip the image upside down by using negative y stride,
  // since OpenGL pixmaps have the bottom scanline first, and oiio writes the top scanline first in the image file.
  int scanlinesize = WinWidth * ImChannels * sizeof(unsigned char);
  if(!outfile->write_image(TypeDesc::UINT8, local_pixmap + (WinHeight - 1) * scanlinesize, AutoStride, -scanlinesize)){
    cerr << "Could not write image to " << outfilename << ", error = " << geterror() << endl;
    ImageOutput::destroy(outfile);
    return;
  }

  // close the image file after the image is written and free up space for the
  // ooio file handler
  outfile->close();
  ImageOutput::destroy(outfile);
}

unsigned char findColorLimit(float colorElement) {
	if(colorElement > 255) {
		return 255;
	}

	else if(colorElement < 0) {
		return 0;
	}

	else {
		return colorElement;
	}
}

// Finds the closest color of a pixel value with respect to the provided palette
void findClosestColor(int colors[][4]) {
	// indices to store the values where the color difference is the least
	int minRed, minGreen, minBlue, minAlpha;

	for(int row=0; row<ImHeight; ++row) {
		for(int col=0; col<ImWidth; ++col) {
			minRed = 0;
			minGreen = 0;
			minBlue = 0;
			minAlpha = 0;

			// mininum difference operations performed for each channel
			for (int r_=0; r_<paletteSize; ++r_) {
				if (abs(dupPixmap[row][col].r - colors[r_][0]) < abs(dupPixmap[row][col].r - colors[minRed][0])) {
					minRed = (r_ + randomization*randomization)%paletteSize;
				}
			}

			// Value of pixmap replaced with the one found in the palette
			pixmap[row][col].r = colors[minRed][0];

			for (int g_=0; g_<paletteSize; ++g_) {
				if (abs(dupPixmap[row][col].g - colors[g_][1]) < abs(dupPixmap[row][col].g - colors[minGreen][1])) {
					minGreen = (g_ + randomization*randomization)%paletteSize;
				}
			}

			pixmap[row][col].g = colors[minGreen][1];
			//pixmap[row][col].g = colors[minRed][1];

			for (int b_=0; b_<paletteSize; ++b_) {
				if (abs(dupPixmap[row][col].b - colors[b_][2]) < abs(dupPixmap[row][col].b - colors[minBlue][2])) {
					minBlue = (b_ + randomization*randomization)%paletteSize;
				}
			}

			pixmap[row][col].b = colors[minBlue][2];
			//pixmap[row][col].b = colors[minRed][2];

			// for (int a_=0; a_<paletteSize; ++a_) {
			// 	if (abs(dupPixmap[row][col].a - colors[a_][3]) < abs(dupPixmap[row][col].a - colors[minBlue][3])) {
			// 		minAlpha = a_;
			// 	}
			// }
      //
			// pixmap[row][col].a = colors[minAlpha][3];

		}
	}
}

// Function that performs Floyd-Steinberg dithering
// Incorporates error-diffusion by including values of colors of its neighbors
void performDithering() {
	// variables for storing the old and new pixel values

	int oldRed, oldGreen, oldBlue; // store values found before operation
	int newRed, newGreen, newBlue; // store values found after operation
	float redDiff, greenDiff, blueDiff; // store color differences by each channel
	float neighborRed, neighborGreen, neighborBlue; // stores channel-wise value of a neighbor

	for(int row=0; row<ImHeight-1; ++row) {
		for(int col=1; col<ImWidth-1; ++col) {
			// old values obtained from original values stored in dupPixmap
			oldRed = dupPixmap[row][col].r;
			oldGreen = dupPixmap[row][col].g;
			oldBlue = dupPixmap[row][col].b;

			// new values obtained from values in pixmap
			newRed = pixmap[row][col].r;
			newGreen = pixmap[row][col].g;
			newBlue = pixmap[row][col].b;
			//cout << "oldRed: " << oldRed << "\tnewRed: " << newRed << "\n";

			// calculating difference of each channel
			redDiff = oldRed - newRed;
			greenDiff = oldGreen - newGreen;
			blueDiff = oldBlue - newBlue;

			// Right neighbor
			// First neighboring element
			neighborRed = pixmap[row][col + 1].r;
			neighborGreen = pixmap[row][col + 1].g;
			neighborBlue = pixmap[row][col + 1].b;

			// error-diffusion of neighboring pixel, channel-wise
			pixmap[row][col + 1].r = findColorLimit(neighborRed + redDiff * 7/16);
			pixmap[row][col + 1].g = findColorLimit(neighborGreen + greenDiff * 7/16);
			pixmap[row][col + 1].b = findColorLimit(neighborBlue + blueDiff * 7/16);

			// Bottom-left neighbor
			// Second neighboring element
			neighborRed = pixmap[row + 1][col - 1].r;
			neighborGreen = pixmap[row + 1][col - 1].g;
			neighborBlue = pixmap[row + 1][col - 1].b;

			pixmap[row + 1][col - 1].r = findColorLimit(neighborRed + redDiff * 3/16);
			pixmap[row + 1][col - 1].g = findColorLimit(neighborGreen + greenDiff * 3/16);
			pixmap[row + 1][col - 1].b = findColorLimit(neighborBlue + blueDiff * 3/16);

			// Bottom neighbor
			// Third neighboring element
			neighborRed = pixmap[row + 1][col].r;
			neighborGreen = pixmap[row + 1][col].g;
			neighborBlue = pixmap[row + 1][col].b;

			pixmap[row + 1][col].r = findColorLimit(neighborRed + redDiff * 5/16);
			pixmap[row + 1][col].g = findColorLimit(neighborGreen + greenDiff * 5/16);
			pixmap[row + 1][col].b = findColorLimit(neighborBlue + blueDiff * 5/16);

			// Bottom-right neighbor
			// Fourth neighboring element
			neighborRed = pixmap[row + 1][col + 1].r;
			neighborGreen = pixmap[row + 1][col + 1].g;
			neighborBlue = pixmap[row + 1][col + 1].b;

			pixmap[row + 1][col + 1].r = findColorLimit(neighborRed + redDiff * 1/16);
			pixmap[row + 1][col + 1].g = findColorLimit(neighborGreen + greenDiff * 1/16);
			pixmap[row + 1][col + 1].b = findColorLimit(neighborBlue + blueDiff * 1/16);
		}
	}
}

// The next few functions try to work on
// the implementation of Median-Cut algorithm
// Alternate function to find the nearest color
// Returns the index of the palette which yields the nearest color
int findNearestColor(Pixel &currentColor, vector<Pixel> palette) {
	float smallestValue = 255;
	int paletteIndex = -1;

	for(int i=0; i < paletteSize; ++i) {
		Pixel paletteColor = palette[i];

		Pixel colorDifference;

		// calculating channel-wise color difference
		colorDifference.r = abs(currentColor.r - paletteColor.r);
		colorDifference.g = abs(currentColor.g - paletteColor.g);
		colorDifference.b = abs(currentColor.b - paletteColor.b);

		float diff = (float) pow(colorDifference.r, 2) + (float) pow(colorDifference.g, 2) + (float) pow(colorDifference.b, 2);

		if (diff < smallestValue) {
			smallestValue = diff;
			paletteIndex = i;
		}
	}

	return paletteIndex;
}

// Function that replaces the current pixel's color value
// with ones found in the fixed-size palette
void reduceCurrentPalette(vector<Pixel> &palette) {
	for(int row=0; row<ImHeight; ++row) {
		for(int col=0; col<ImWidth; ++col) {
			Pixel currentPixel = pixmap[row][col];
			int paletteIndex = findNearestColor(currentPixel, palette);
			pixmap[row][col].r = palette[paletteIndex].r;
			pixmap[row][col].g = palette[paletteIndex].g;
			pixmap[row][col].b = palette[paletteIndex].b;
			pixmap[row][col].a = palette[paletteIndex].a;
		}
	}
}

// Function that calculates the root mean square of two pixel values
void rootMeanSquare(vector<Pixel> &pixels, int startPoint, int endPoint, double &sumRed, double &sumGreen, double &sumBlue) {
	double n = (double) endPoint - startPoint + 1;

	for(int i=startPoint; i<endPoint+1; ++i) {
		double red = (double) pixels[i].r;
		double green = (double) pixels[i].g;
		double blue = (double) pixels[i].b;

		sumRed += red * red;
		sumGreen += green * green;
		sumBlue += blue * blue;
	}

	sumRed = (unsigned char) sqrt(sumRed/n);
	sumGreen = (unsigned char) sqrt(sumGreen/n);
	sumBlue = (unsigned char) sqrt(sumBlue/n);
}

// Function that returns the channel that has the most difference in
// a given range
char getWidestChannel(vector<Pixel> &pixels, int start, int end) {
	unsigned char minRed, maxRed, minGreen, maxGreen, minBlue, maxBlue;

	minRed = 255;
	minGreen = 255;
	minBlue = 255;

	maxRed = 0;
	maxGreen = 0;
	maxBlue = 0;

	for(int i=0; i<end+1; +i) {
		if (pixels[i].r < minRed) {
			minRed = pixels[i].r;
		}

		if (pixels[i].g < minGreen) {
			minGreen = pixels[i].g;
		}

		if (pixels[i].b < minBlue) {
			minBlue = pixels[i].b;
		}

		if (pixels[i].r > maxRed) {
			maxRed = pixels[i].r;
		}

		if (pixels[i].g > maxRed) {
			maxGreen = pixels[i].g;
		}

		if (pixels[i].b > maxBlue) {
			maxBlue = pixels[i].b;
		}
	}

	// range variables to store differences between maximum and minimum values
	unsigned char redRange = maxRed - minRed;
	unsigned char greenRange = maxGreen - minGreen;
	unsigned char blueRange = maxBlue - minBlue;

	// Returns the specific channel with the most difference
	if (redRange > greenRange && redRange > blueRange) {
		return 'r';
	}

	if (blueRange > redRange && blueRange > greenRange) {
		return 'b';
	}

	if (greenRange > redRange && greenRange > blueRange) {
		return 'g';
	}
}

bool compareR (Pixel colorA, Pixel colorB) {
	return colorA.r < colorB.r;
}

bool compareG (Pixel colorA, Pixel colorB) {
	return colorA.g < colorB.g;
}

bool compareB (Pixel colorA, Pixel colorB) {
	return colorA.b < colorB.b;
}

// Function that recursively calls itself around a median point
// Sorts the obtained values
void medianCutAlgorithm(vector<Pixel> &pixels, vector<Pixel> &palette, int startPoint, int endPoint, int paletteIndex, int length) {

	// For empty case
	if (length == 0) {
		double sumRed = 0, sumGreen = 0, sumBlue = 0;
		Pixel &color = palette[paletteIndex++];
		rootMeanSquare(pixels, startPoint, endPoint, sumRed, sumGreen, sumBlue);
		color.r = (unsigned char) sumRed;
		color.g = (unsigned char) sumGreen;
		color.b = (unsigned char) sumBlue;
		color.a = 255;
	}

	else {
		char channelInfo = getWidestChannel(pixels, startPoint, endPoint);
		if (channelInfo = 'r') {
			sort(pixels.begin() + startPoint, pixels.begin() + endPoint + 1, compareR);
		}

		else if (channelInfo = 'g') {
			sort(pixels.begin() + startPoint, pixels.begin() + endPoint + 1, compareG);
		}

		else if (channelInfo = 'b') {
			sort(pixels.begin() + startPoint, pixels.begin() + endPoint + 1, compareB);
		}

		int midPoint = (startPoint + endPoint) / 2;

		// Recursion around the midPoint (median point)
		medianCutAlgorithm(pixels, palette, startPoint, midPoint, paletteIndex, length-1);
		medianCutAlgorithm(pixels, palette, midPoint+1, endPoint, paletteIndex, length-1);
	}
}

// Function to perform the median cut algorithm
void performMedianCut(vector<Pixel> &pixels, vector<Pixel> &palette) {
	int length = palette.size();
	double logLength = log2(length);

	int paletteIndex = 0;

	medianCutAlgorithm(pixels, palette, 0, pixels.size()-1, paletteIndex, logLength);
}

// To be completed (for median cut algorithm implementation)

/*
// Function that returns true if a certain color is present in the given data structure
bool PixelExists(unordered_set<Pixel> &unique, Pixel &pixelValue) {
	for(int i=0; i<unique.size(); ++i) {
		// Compares channel value of given pixel value with each item's corresponding channel values
		if (pixelValue.r == (int) unique[i].r && pixelValue.g == (int) unique[i].g && pixelValue.b == (int) unique[i].b) {
			return true;
		}
	}

	// Returns false if the pixel does not exist
	return false;
}

// Data-type error
void findReducedPalette(vector<Pixel> &palette) {
	unordered_set<Pixel, UniqueIdentifier> unique;

	for(int row=0; row<ImHeight; ++row) {
		for(int col=0; col<ImWidth; ++col) {
			Pixel currentPixel = pixmap[row][col];

			if(!PixelExists(unique, currentPixel)) {
				unique.insert(currentPixel);
			}
		}
	}
}
*/

//
//   Display Callback Routine: clear the screen and draw the current image
//
void handleDisplay(){

  // specify window clear (background) color to be opaque black
  glClearColor(0, 0, 0, 1);
  // clear window to background color
  glClear(GL_COLOR_BUFFER_BIT);

  // only draw the image if it is of a valid size
  if(ImWidth > 0 && ImHeight > 0)
    displayimage();

  // flush the OpenGL pipeline to the viewport
  glFlush();
}

//
//  Keyboard Callback Routine: 'r' - read and display a new image,
//  'w' - write the current window to an image file, 'q' or ESC - quit
//
void handleKey(unsigned char key, int x, int y){
  string infilename, outfilename;
  int ok;

  switch(key){
    case 'r':		// 'r' - read an image from a file
    case 'R':
		if(!isImageSet) {
      cout << "Enter the path to the input image: ";	  // prompt user for input filename
      cin >> infilename;
			imageFileName = infilename;
      ok = readimage(infilename);
      if(ok){
        glutReshapeWindow(ImWidth, ImHeight); // OpenGL window should match new image
				isImageSet = true; // Store image path; prevents repeated prompts to read image
			}
			glutPostRedisplay();
		}
		else {
			ok = readimage(imageFileName);
			if(ok){
				glutReshapeWindow(ImWidth, ImHeight); // OpenGL window should match new image
			}
			glutPostRedisplay();
		}
      break;


		case 'c':
		case 'C':
			findClosestColor(colors);
			glutPostRedisplay();
			break;

		case 'd':
		case 'D':
			performDithering();
			glutPostRedisplay();
			break;

    case 'w':		// 'w' - write the image to a file
    case 'W':
      cout << "Enter the path for the output image: ";  // prompt user for output filename
      cin >> outfilename;
      writeimage(outfilename);
      break;

    case 'q':		// q or ESC - quit
    case 'Q':
    case 27:
      destroy();
      exit(0);

    default:		// not a valid key -- just ignore it
      return;
  }
}

//
//  Reshape Callback Routine: If the window is too small to fit the image,
//  make a viewport of the maximum size that maintains the image proportions.
//  Otherwise, size the viewport to match the image size. In either case, the
//  viewport is centered in the window.
//
void handleReshape(int w, int h){
  float imageaspect = (float)ImWidth / (float)ImHeight;	// aspect ratio of image
  float newaspect = (float)w / (float)h; // new aspect ratio of window

  // record the new window size in global variables for easy access
  WinWidth = w;
  WinHeight = h;

  // if the image fits in the window, viewport is the same size as the image
  if(w >= ImWidth && h >= ImHeight){
    Xoffset = (w - ImWidth) / 2;
    Yoffset = (h - ImHeight) / 2;
    VpWidth = ImWidth;
    VpHeight = ImHeight;
  }
  // if the window is wider than the image, use the full window height
  // and size the width to match the image aspect ratio
  else if(newaspect > imageaspect){
    VpHeight = h;
    VpWidth = int(imageaspect * VpHeight);
    Xoffset = int((w - VpWidth) / 2);
    Yoffset = 0;
  }
  // if the window is narrower than the image, use the full window width
  // and size the height to match the image aspect ratio
  else{
    VpWidth = w;
    VpHeight = int(VpWidth / imageaspect);
    Yoffset = int((h - VpHeight) / 2);
    Xoffset = 0;
  }

  // center the viewport in the window
  glViewport(Xoffset, Yoffset, VpWidth, VpHeight);

  // viewport coordinates are simply pixel coordinates
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluOrtho2D(0, VpWidth, 0, VpHeight);
  glMatrixMode(GL_MODELVIEW);
}

// Function to define the palette with values found
// in fixed-intervals
void definePalette() {
	// Populating the color palette with values
	for (int i = 0; i < paletteSize; ++i) {
		for (int j = 0; j < 4; ++j) {
				colors[i][j] = ((256/paletteSize)*(i+1)) - 1;
		}
	}
}

// Function to display the contents of the color palette
void displayPalette() {
	// Displaying the color palette's values
	for(int i=0; i < paletteSize; ++i) {
		for(int j=0; j<4; ++j) {
			cout << colors[i][j] << "\t";
		}
		cout << "\n";
	}
}

//
// Main program to scan the commandline, set up GLUT and OpenGL, and start Main Loop
//
int main(int argc, char* argv[]){

	// to prevent the same number from being generated everytime
	srand(time(NULL)); // setting the seed for the random number

	randomization = rand() % 20; // setting the randomization factor to be less than 20

	definePalette();

	cout << "The RGBA color values are given as follows: " << "\n\n";
	cout << "R\tG\tB\tA" << "\n";
	displayPalette();
	cout << "\n";

  // scan command line and process
  // only one parameter allowed, an optional image filename and extension

  // set up the default window and empty pixmap if no image or image fails to load
  WinWidth = DEFAULTWIDTH;
  WinHeight = DEFAULTHEIGHT;
  ImWidth = 0;
  ImHeight = 0;

  // load the image if present, and size the window to match
  if(argc == 2){
		imageFileName = argv[1];
		isImageSet = true;
    if(readimage(argv[1])){
      WinWidth = ImWidth;
      WinHeight = ImHeight;
    }
  }

	if (argc > 2) { // Info to be displayed in case the number of command line arguments exceed 5
		cout << "Usage: ./rod <inputImage.ext>" <<"\n";
		cout << "Images are present in images/" <<"\n";
		exit(1);
	}
  // start up GLUT
  glutInit(&argc, argv);

  // create the graphics window, giving width, height, and title text
  glutInitDisplayMode(GLUT_SINGLE | GLUT_RGBA);
  glutInitWindowSize(WinWidth, WinHeight);
  glutCreateWindow("Random Order Dithering");

  // set up the callback routines
  glutDisplayFunc(handleDisplay); // display update callback
  glutKeyboardFunc(handleKey);	  // keyboard key press callback
  glutReshapeFunc(handleReshape); // window resize callback


  // Enter GLUT's event loop
  glutMainLoop();
  return 0;
}
