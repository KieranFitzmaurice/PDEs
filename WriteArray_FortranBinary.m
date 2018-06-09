function WriteArray_FortranBinary(filename,A)

% This function writes a D-dimensional array in the format of unformatted
% FORTRAN binary.
%
% Includes record length information to make it easy for FORTRAN to read.  
%
% filename = 'myfile.dat'
% D = dimension of array (n=1 for 1D, n=2 for 2D ...)
%
% Example output file format for 2D array:
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

sz = size(A); % Length in each dimension
D = length(sz); % Number of dimensions of array

fileID = fopen(filename,'w');

% Write dimensions of array at top of output file
for i = 1:D
    fwrite(fileID,4,'int32');     % record length (bytes)
    fwrite(fileID,sz(i),'int32'); % dimension length 
    fwrite(fileID,4,'int32');     % record length
end

% Write array to file in column order
fwrite(fileID,prod(sz)*8,'int32'); % record length (bytes)
fwrite(fileID,A,'real*8');
fwrite(fileID,prod(sz)*8,'int32'); % record length 

fclose(fileID);

end