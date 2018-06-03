def ReadArray_FortranBinary(filename,D):
    """
    This function allows the user to read in a D-dimensional array
    from unformatted Fortran binary.

    Requires scipy.io.FortranFile
    
    Example input file format for 2D array:

      Line | Entry         |    Data Type
     ------------------------------------
        1  | nrows         |     int32
        2  | ncols         |     int32
        3  | A(1,1)        |     double
        4  | A(2,1)        |     double
        .  |   .           |       .
        .  |   .           |       .
        .  |   .           |       .
       end | A(nrows,ncols)|     double
    """

    f = FortranFile(filename,'r')
    sz = np.zeros(D)

    for i in range(0,D):
        sz[i] = f.read_ints(dtype=np.int32)

    A = f.read_reals(dtype=np.float64)
    f.close()

    A = A.reshape(sz.astype(int),order = 'F')
    A = np.transpose(A)

    return(A)
