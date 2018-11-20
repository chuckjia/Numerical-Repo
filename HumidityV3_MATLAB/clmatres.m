function clmatres()

prompt = "\n  - Confirm to clear all c result files (Y/N): ";
act = lower(input(char(prompt), 's'));

% if strcmp(act, 'y') || strcmp(act, 'yes')
%     rmdir('Output')
%     mkdir('Output')
%     fprintf("MATLAB results all cleared.\n");
% else
%     fprintf("Aborted. Nothing is changed.\n");
% end

end