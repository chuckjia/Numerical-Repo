function projectPath = genFolderPathName(projectPath)
%FIXFOLDERPATHNAME Summary of this function goes here
%   Detailed explanation goes here

if ischar(projectPath)
    lastChar = projectPath(end);
    projectPath = string(projectPath);
else
    projectPath_temp = char(projectPath);
    lastChar = projectPath_temp(end);
end

if lastChar ~= '/'
    projectPath = projectPath + "/";
end

end

