function fileList = listfiles(folder)
%LISTFILES Returns list of all files and folders within a folder

if isunix
    fileList = dir(folder);
    fileList = {fileList.name};
    
    toDelete = ismember(fileList, {'.', '..', '.DS_Store'});
    fileList(toDelete) = [];
    return
end

fileList = false;
fprintf("listfiles: Unknown OS platform!!\n");

end

