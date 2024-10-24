# BregmanCookbook
 Different image processing algorithms implemented using the split Bregman iterations
 
Current version: 
3.2 - Convergence criteria are updated (fix some bugs) and homogenize for all functions + test script for 1D case

Previous version:
3.0 - New toolbox organization in different subfolders
2.0 - Isotropic ROF is added and both Isotropic and Anisotropic ROF are now solved in the Fourier domain.
1.2 - Deconvolution functions now accept both kernels defined in spatial or Fourier domains.
1.0
 
Date: 06/23/2013

This toolbox provides the source code associated with the Bregman Cookbook available at
http://jegilles.sdsu.edu

This software is free and should be used only for nonprofit purposes. Any
unauthorized use of this software for industrial or profit-oriented
activities is expressively prohibited.

Some functions need some extra toolbox:
- Framelet toolbox, written by Jianfeng Cai and available at https://jegilles.sdsu.edu/code/Framelet.zip
- Curvelab toolbox, freely available for non-profit purposes at http://www.curvelet.org

These toolbox are assumed to be installed properly on your computer (ie the corresponding
Paths are added into Matlab).

This toolbox contains the following folders and files:

Root folder:
-README : this file!

Doc:
- BregmanCookbook.pdf

1D:
- L1_SplitBregmanIteration.m : performs the recovery of a sparse signal affected by a known linear operator

2D:
- ATV_NB_Deconvolution.m : performs the Nonblind Anisotropic Total Variation Deconvolution
- ATV_ROF.m : performs the Anisotropic Total Variation Denoising
- ITV_ROF.m : performs the Isotropic Total Variation Denoising
- ITV_NB_Deconvolution.m : performs the Nonblind Isotropic Total Variation Deconvolution
- Curvelet_NB_Deconvolution.m : performs the Nonblind Deconvolution based on Curvelet sparsity
- Framelet_NB_Deconvolution.m : performs the Nonblind Deconvolution based on Framelet sparsity (Analysis approach)
- Framelet_NB_Deconvolution2.m : performs the Nonblind Deconvolution based on Framelet sparsity (Synthesis approach)

3D:
- ATV_NB_Deconvolution_3D.m: performs the Nonblind Anisotropic Total Variation Deconvolution
- ITV_NB_Deconvolution_3D.m: performs the Nonblind Isotropic Total Variation Deconvolution
- ATV_ROF_3D.m : performs the Anisotropic Total Variation Denoising
- ITV_ROF_3D.m : performs the Isotropic Total Variation Denoising
- Curvelet_NB_Deconvolution_3D.m : performs the Nonblind Deconvolution based on Curvelet sparsity

Examples:
- lena.mat: Lena image in double format normalized between 0 and 1
- Test1D: customizable script performing a test of the 1D function
- Test2D: customizable script performing a test of the 2D functions
- Test3D: customizable script performing a test of the 3D functions
- TVG_CartoonTexture_Decomposition: function performing the TV-G cartoon+textures decomposition

Utils:
- AddCurveletArray.m : sum the curvelet coefficients of two decomposition structures
- AddFrameletArray.m : sum the framelet coefficients of two decomposition structures
- ShrinkComplexCurvelet.m : performs the shrinkage of complex curvelet coefficients
- ShrinkCurvelet.m : performs the shrinkage of curvelet coefficients
- ShrinkFramelet.m : performs the shrinkage of framelet coefficients
- SubCurveletArray.m : substract the curvelet coefficients of two decomposition structures
- SubFrameletArray.m : substract the framelet coefficients of two decomposition structures

For any questions, bugs report, suggestions, ... Feel free to contact me at jgilles@sdsu.edu
