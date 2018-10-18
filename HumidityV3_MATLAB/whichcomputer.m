function [accname, accno] = whichcomputer()
%WHICHCOMPUTER Check which computer I am using
%   Currently, the following computer accounts are supported
%       'macOS_chuckjia', accno == 1 : Account name "chuckjia" on Macbook Pro 2015
%       'ubuntu_chuck',   accno == 2 : Account name "chuck"    on office Ubuntu
%

if isunix
    username = getenv('USER');
    if ismac && strcmp(username, "chuckjia")
        accname = 'macOS_chuckjia';
        accno = 1;
        return
    elseif strcmp(username, "chuck")
        accname = 'ubuntu_chuck';
        accno = 2;
        return
    end
end

accname = "unkown_account";

end

