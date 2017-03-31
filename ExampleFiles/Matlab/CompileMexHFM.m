% Copyright Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay

% Specific project headers
sourceDir.Headers = '../..';

% Shared headers between my projects
sourceDir.JMM_CPPLibs = '../../../../JMM_CPPLibs';
sourceDir.LinearAlgebra =       [sourceDir.JMM_CPPLibs '/LinearAlgebra'];
sourceDir.Output =              [sourceDir.JMM_CPPLibs '/Output'];
sourceDir.DataStructures =      [sourceDir.JMM_CPPLibs '/DataStructures'];


if verLessThan('matlab','8.1')
    cxxFlags = ['CXXFLAGS="-std=c++11" ' ...
        'CXXLIBS="\$CXXLIBS -lc++" ' ]; % This flag is required on some platforms, but must be commented on others...
    outputFlag = '-o ';
else
    cxxFlags = ['CXXFLAGS="-std=c++11" ']; 
    outputFlag = '-output ';
end

compileHFM = @(binaryDir,flag) eval(['mex ' ...
    outputFlag 'MatlabHFM_' flag ' ../../MatlabHFM.cpp' ... %Changing the output name does not work ???
    ' -outdir ' binaryDir ...
    ' ' cxxFlags  '-D' flag ...
    ' -I' sourceDir.Headers ...
    ' -I' sourceDir.LinearAlgebra ...
    ' -I' sourceDir.DataStructures ...
    ' -I' sourceDir.Output ...
    ]);

compileHFMAll = @(binaryDir) cellfun(@(flag) compileHFM(binaryDir,flag),{'Custom','State','TimeVar','Isotropic','Curvature2','Curvature3','Riemann'});
%moveFun = @(source,target) ['movefile ' targetDir '/' source '.' mexext ' ' targetDir '/' target '.' mexext];

fprintf(['\nPlease execute the function compileHFMAll(binaryDir) to build \n'...
'the Hamilton Fast Marching executables in directory binaryDir, \n' ...
'or compileHFM(binaryDir,flag) for a specific model.\n'])

%For me : binaryDir = '/Users/mirebeau/Dropbox/Programmes/MATLAB/MexBin';
