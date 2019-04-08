function path = defaultProjectPath()
%DEFAULTPROJECTPATH This function returns the default path for two computers I use regularly

[~, accno] = whichcomputer();
if accno == 1  % On Macbook Pro 2015
    path = "~/Documents/Workspace/Git/Numerical-Repo/HumidityV3/";
elseif accno == 2
    path = "~/git/Numerical-Repo/HumidityV3/";
else
    path = "ThisIsNotAKnownAccount/";
end
    
end

