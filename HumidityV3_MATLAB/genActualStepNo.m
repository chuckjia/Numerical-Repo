function stepNo = genActualStepNo(stepNo, resultFolderFullPath, solnName, printMsg)
% GENACTUALSTEPNO Generate step no. Its main function is to convert negative step numbers to the latest result
%                 in the result folder

if nargin < 4
    printMsg = true;
end

if stepNo < 0
    % Automatically get last file in folder
    stepNo = getStepNoOfLastFile(resultFolderFullPath, solnName);
    if printMsg
        fprintf("Printing the latest time in simulation: step no. %d\n", stepNo);
    end
else
    if printMsg
        fprintf("Printing a specific time in simulation: step no. %d\n", stepNo);
    end
end

end
