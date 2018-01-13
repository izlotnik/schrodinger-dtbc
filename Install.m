%% get current path
sCurLocation = cd;

%% add local folder to Matlab path
addpath(sCurLocation);
addpath(genpath([ sCurLocation '/module'    ]));
addpath(genpath([ sCurLocation '/parameter' ])); 
addpath(genpath([ sCurLocation '/solver'    ]));

%% clear variale
clear sCurLocation;

%% disable warnings
warning('off','MATLAB:dispatcher:InexactCaseMatch');
warning('off','MATLAB:rmpath:DirNotFound');