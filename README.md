# Quick Intro:
This folder contains all the functions and scripts used to implement the research paper titled
"A Hessian-Free Gradient Flow (HFGF) method for the optimization of deep learning neural networks".

# How to run the code?
In order to test the code, and obtain the result which is Table 2 in the main research paper, which can be accessed below,
download this folder, add its path to the Matlab, and simply run Test_HFGF script.

# Note:
Due to the limited access to a high computational power machine, the high dimensional functions were not tested.
Feel free to test them if possible.

# What is happening underneath the hood?
Test_HFGF script will call the main code of the paper named main_HFGF and will send 4 different testing functions 
and return the information needed to construct the table. Inside the main_HFGF, there is a call to Armijo function to obtain the step size alpha.

# Reference:
https://www.sciencedirect.com/science/article/abs/pii/S0098135420303562 

# Contact information:
atrashabdulkarim@gmail.com
