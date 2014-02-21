OCR_formula2latex
=================

Convert a written mathematical formula to latex.

Our code base can be found in the src folder. We broke our code into three folders based on functionality: image_segmentation, feature_extraction, and machine_learning. It is possible to call the individual functions found in these folders on images. However, it will probably be of most intrest to see the complete system work together by calling the function "convert" in the top level of the source folder.

convert(img, V)
  -img: filename to image with mathematical formula
  -V: optional argument to run the program in verbose mode, if set to true (or 1) the program will display images showing each step of the process
  -returns a string representation of the original image 

The "demo" folder contains many interesting testing images.
  Ex.) convert('demo/img_000.jpg', 1) % run in verbose mode