import matplotlib.pyplot as plt
import numpy as np
from PDE_Modules import ReadArray_FortranBinary, MakeGIF_2D

inputfile = "/Users/kieranfitzmaurice/Documents/MATLAB/Spinodal.dat"
outputfile = "/Users/kieranfitzmaurice/Desktop/Spinodal.gif"

A = ReadArray_FortranBinary(inputfile,3)
MakeGIF_2D(outputfile,A)
