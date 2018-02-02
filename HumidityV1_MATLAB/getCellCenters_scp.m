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

matShape_getCellCentersScp = [numCellsP, numCellsX];
if (exist(folder_getParScp) == 0)
    folder_getParScp = "Results/";
end
cellCentersX = reshape(getVecFromFile_fcn(folder_getParScp, "cellCentersX"), matShape_getCellCentersScp);
cellCentersP = reshape(getVecFromFile_fcn(folder_getParScp, "cellCentersP"), matShape_getCellCentersScp);

clear matShape_getCellCentersScp

