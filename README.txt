This project is called Random Order Dithering. The source code for the project is located in the same folder, and can be found in 'rod.cpp'. It uses the principle of color quantization, and retrieves colors from a palette with limited entries.

A few sample images have been provided in the images/ folder. Similarly, a few output images have been placed in the outputs/ folder. The palette size has been defined as 8. The user can change it by changing line 14 found in 'rod.cpp' from the following:

#define paletteSize 8

to something like this:

#define paletteSize <yourNumber>

where <yourNumber> can be a user-defined number, although it has to be of the form 2^n, with n ranging from 1 to 8. This means that the paletteSize should be numbers that can be expressed as exponents of 2. The numbers that are valid are:

[2, 4, 8, 16, 32, 64, 128, 256]

However, the higher the number is, the slower the execution of the operations will be (Setting it to 2 might yield a completely gray image! The program would then be similar in concept to thresholding. However, subsequent operations as given in this document would still work as intended).

Before running the program, the following command has to be run into the terminal:

$ make

This will compile the C++ program and generate an object file and an executable by the name of 'rod'.

There are two ways of running the program:

1. $ ./rod
2. $ ./rod <inputImage.ext>

In both the methods, the first argument will run the executable.

For the first method, once the program has been run, an OpenGL window will be opened as well. In order to read an image, the user needs to select the OpenGL window and press 'r' or 'R'. This will prompt the user to provide the path to the image that he wishes to read. Along with the image name, the extension should be provided as well.

For the second method, the image path has to be provided while running the program.

Examples are given as:

$ make
$ ./rod

Alternatively,

$ make
$ ./rod images/image2.jpg

Once the OpenGL window opens, the user can press 'c' or 'C' to perform the color quantization operation. Pressing 'd' or 'D' will perform the Floyd-Steinberg dithering operation. The user can press 'd' or 'D' repeatedly to observe the results of repeated operations. In such cases, more noise will be introduced each time, but the image will start to resemble the original image. At any time in the procedure, pressing 'r' or 'R' will restore the original image again. Pressing 'w' or 'W' will prompt the user to write the image onto disk. Here, the correct extension should be provided as well.

Pressing 'q' or 'Q' or 'Esc' will close the OpenGL window and quit the program.

In order to observe the best results, after the image has been read, the user should press 'c' on the OpenGL window to perform the color quantization operation (The same operation can yield different results every time, based on how long the user waits before pressing 'c'. For example, press 'c' 1 second after the window is active for the first attempt; run the program again but this time wait 3 seconds before pressing 'c' and continuing with the rest of the steps). Once this step has been performed, the user should prompt 'd' or 'D' on the OpenGL window to perform the dithering operation. The user can press 'd' or 'D' multiple times to repeatedly perform the operation, observing more and more noise with each iteration.

The example of the above can be given as:

1. $ make (Enter on the terminal)
2. $ ./rod images/image4.jpg (Enter on the terminal)
3. Press 'c' on the OpenGL window
4. Press 'd' on the OpenGL window (applied once)
5. Press 'd' on the OpenGL window (applied twice)
6. Press 'd' on the OpenGL window (applied thrice)
7. Press 'w' on the OpenGL window to write the output image onto disk
8. Provide a valid file name (for instance, "sample.png")
9. Press 'q' on the OpenGL window (exits the program)

After the executable has been run, the user can delete the executable and the object files by running the following command in the terminal:

$ make clean

Issues with current iteration:
- Data-type errors observed when using std::unordered_set<>
- Implementation of the Median Cut algorithm can be possible if a data structure can be used that can optimally store the unique pixels occurring in the image
