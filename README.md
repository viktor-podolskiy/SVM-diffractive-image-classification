# SVM-diffractive-image-classification
Machine learning-based diffractive image analysis with subwavelength resolution

The enclosed files illustrate Support-Vector-Machine-based diffractive imaging with subwavelength resolution, 
see A.Ghosh et.al., ACS Photonics (2021)

(c) 2019-2021 V. Podolskiy, U Mass Lowell
--- for non-commercial use only ---

Typical workflow: 

(1) Run makeImgArr.m to generate an array of diffractive signatures of phantom objects
(2) Run img2rphiTheory.m to read the array of images, add optional noise, and save data in polar coordinates
(3) Run rphiTheorySVM.m to train a variety of SVMs and analyze the best combinations of parameters l,m,j to use for object classification
