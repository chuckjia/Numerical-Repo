function stepNo = genActualStepNo(stepNo, resultFolderFullPath, solnName)
% GENACTUALSTEPNO Generate step no. Its main function is to convert negative step numbers to the latest result
%                 in the result folder

if stepNo < 0
    % Automatically get last file in folder
    stepNo = getStepNoOfLastFile(resultFolderFullPath, solnName);
    fprintf("Printing the latest time in simulation: step no. %d\n", stepNo);
else
    fprintf("Printing a specific time in simulation: step no. %d\n", stepNo);
end

end