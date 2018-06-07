import matplotlib.pyplot as plt
import numpy as np
from PDE_Modules import ReadArray_FortranBinary, MakeGIF_2D

inputfile = "/Users/kieranfitzmaurice/Documents/REU_2018/PDE_data/temp_file.dat"
outputfile = "/Users/kieranfitzmaurice/Documents/REU_2018/PDE_data/temp_file.gif"

A = ReadArray_FortranBinary(inputfile,3)
MakeGIF_2D(outputfile,A)
