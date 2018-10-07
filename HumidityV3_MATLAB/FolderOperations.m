%% Set work directory

platform = 0;
if ismac  % Macbook Pro
    cd /Users/chuckjia/Documents/Workspace/git/Numerical-Repo/HumidityV3_MATLAB  % Macbook Pro
elseif isunix
    platform = 1;
    cd ~/git/Numerical-Repo/HumidityV3_MATLAB  % Office Ubuntu
else
    platform = -1;
end

if platform == -1
    fprintf("Unknown platform!!\n");
end

