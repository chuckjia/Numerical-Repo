function filename = getLastFileInDir(dirFullPath, startChar)
%GETLASTFILEINDIR Summary of this function goes here
%   Detailed explanation goes here

allFileListing = dir(dirFullPath);
filename = allFileListing(end).name;

if nargin >= 2
    allFileListing_len = length(allFileListing);
    maxFilenameLen = 0;
    
    for i = 1:allFileListing_len
        currFilename = allFileListing(i).name;
        currFilenameLen = length(currFilename);
        if currFilename(1) == startChar && currFilenameLen >= maxFilenameLen
            filename = currFilename;
            maxFilenameLen = currFilenameLen;
        end
    end
end

end

