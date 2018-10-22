function mvcres(destination)
%MVCRES Move the result from C++ computation to temporary storage

[~, accno] = whichcomputer();
projectPath = "~/Documents/Workspace/git/Numerical-Repo/HumidityV3/";
if accno == 1
    fprintf("On Macbook Pro 2015 | chuckjia\n");
    if nargin == 0
        destination = "~/Documents/Workspace/DataStorage/Humidity/temp/";
    end
elseif accno == 2
    projectPath = "~/git/Numerical-Repo/HumidityV3/";  % Path on Ubuntu
    fprintf("On office Ubuntu | chuck\n");
    if nargin == 0
        destination = "~/Documents/workspace/DataStorage/Humidity/temp/";
    end
else
    fprintf("Unknown computer account. Abort. Nothing is changed.\n");
    return
end

fprintf("    Project path = %s\n", projectPath);
MovieFramesPath = projectPath + "MovieFrames/";
OutputPath = projectPath + "Output/";
filelist_movieframes = listfiles(MovieFramesPath);
filelist_output = listfiles(OutputPath);
numFilesToList = 5;
fprintf("    Sample files from 'MovieFrames':  %s, %s, %s, %s, %s\n", filelist_movieframes{1:numFilesToList})
fprintf("    Sample files from 'Output'     :  %s, %s, %s, %s, %s\n", filelist_output{1:numFilesToList})


destination = destination + strjoin(strsplit(num2str(fix(clock))), '_');
fprintf("Moving these 2 folders to %s\n", destination);
prompt = "Press Y to proceed: ";
act = lower(input(char(prompt), 's'));

if strcmp(act, 'y') || strcmp(act, 'yes')
    mkdir(destination);
    movefile(MovieFramesPath, destination);
    movefile(OutputPath, destination);
    mkdir(MovieFramesPath);
    mkdir(OutputPath);
    fprintf("Moved all output files to: " + destination + "\n");
else
    fprintf("Aborted. Nothing is changed.\n");
end
end

