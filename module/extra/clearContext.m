%% Script clearContext
clc;                % clear command line
clear all;          % clear memory
close all;          % close all figures
close all hidden;	% close waitbar
fclose('all');      % close all files
%%
rmpath([ cd '/parameter/function/1D' ]);
rmpath([ cd '/parameter/function/2D' ]);