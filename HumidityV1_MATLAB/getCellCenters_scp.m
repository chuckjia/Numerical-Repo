% getCellCenters_sct.m - Read cell centers from humidity output files 
% 
% This script reads cell center coordinates from humidity output files and 
% stores them in vectors cellCentersX and cellCentersP. Both vectors are 
% reshaped in 2D and with dimension [numCellsP, numCellX].
% 
% Note: because of the reshape() function, the two vectors are in 
% column-major order.
%
% To run this file, getPar_sct is required to be executed first
%
% Author: Chuck Jia
% Created on: Oct 20, 2017

matShape = [numCellsP, numCellsX];
folder = "Results/";
cellCentersX = reshape(getVecFromFile_fcn(folder, "cellCentersX"), matShape);
cellCentersP = reshape(getVecFromFile_fcn(folder, "cellCentersP"), matShape);

clear matShape folder

