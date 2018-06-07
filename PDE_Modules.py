def ReadArray_FortranBinary(filename,D):
    """
    This function allows the user to read in a D-dimensional array
    from unformatted Fortran binary.

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
    import numpy as np
    from scipy.io import FortranFile

    f = FortranFile(filename,'r')
    sz = np.zeros(D)

    for i in range(0,D):
        sz[i] = f.read_ints(dtype=np.int32)

    A = f.read_reals(dtype=np.float64)
    f.close()

    A = A.reshape(sz.astype(int),order = 'F')
    A = np.transpose(A)

    return(A)

def MakeGIF_2D(outputfile,A):
    """
    This function allows the user to create a gif of a 2D system as it
    evolves in time.

    A(t,m,n) is a 3D array which holds 'snapshots' of an m x n array in time.
    """
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.animation import FuncAnimation

    numframe,m,n = A.shape

    height = 6
    width = float(n/m)*height
    delay_in_ms = 50

    fig, ax = plt.subplots(nrows = 1,ncols = 1,figsize = (height,width))
    fig.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=None, hspace=None)

    cmap = plt.get_cmap('jet')

    def animate(frame):
        ax.clear()
        ax.pcolormesh(A[frame,:,:],cmap = cmap, vmin = 0,vmax = 1)
        ax.axis('off')

    anim = FuncAnimation(fig, animate, frames = (numframe), interval = delay_in_ms)
    anim.save(outputfile, writer='imagemagick', fps=30)
