function A = ReadArray_FortranBinary(filename,D)

% This function reads in an D-dimensional array from unformatted FORTRAN
% Skips over record length information to properly extract data
%
% filename = 'myfile.dat'
% D = dimension of array (n=1 for 1D, n=2 for 2D ...)
%
% Example input file format for 2D array:
%
% Line | Entry         |    Data Type
%------------------------------------
%   1  | nrows         |     int32 
%   2  | ncols         |     int32 
%   3  | A(1,1)        |     double 
%   4  | A(2,1)        |     double
%   .  |   .           |       .
%   .  |   .           |       .
%   .  |   .           |       .
%  end | A(nrows,ncols)|     double

size = zeros(1,D);
fileID = fopen(filename,'rb');

% Skip header stuff
fseek(fileID, 4, 'cof'); 

for i = 1:D
    % Read in length of array's ith dimension
    size(i) = fread(fileID,1,'int32');
    
    % More header stuff to skip
    fseek(fileID, 8, 'cof');
end

% Read in data
A = fread(fileID,inf,'double');
fclose(fileID);

% Fortran is column-ordered, so must transpose array
A = reshape(A,size);
A = A.';
