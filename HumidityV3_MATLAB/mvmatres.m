function mvmatres(destination)
%MVMATRES Move all matlab results to an another folder for storage

currFileList = dir("Output/");
if isempty(listfiles("Output/"))
    fprintf("No output files exist. No file changed or moved.\n");
    return
end

username = getenv('USER');
if nargin == 0
    if ismac && strcmp(username, "chuckjia")
    destination = "~/Documents/Workspace/DataStorage/Humidity/temp/";
    else
        fprintf("mvmatres: Unknown computer account. Abort. Nothing is changed.\n");
        return
    end
end
destination = destination + "MATLAB_" + strjoin(strsplit(num2str(fix(clock))), '_');
mkdir(destination);
movefile("Output/*", destination);
fprintf("Moved all output files to: " + destination + "\n");

end