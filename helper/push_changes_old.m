function push_changes_old(msg)
%PUSH_CHANGES  Stage, commit, and push all repo changes to GitHub.
%
% Usage:
%   push_changes("Update AHinCBP function")

    if nargin < 1 || isempty(msg)
        msg = "Update from MATLAB";
    end

    % Show current status
    fprintf('--- Git Status BEFORE ---\n');
    system('git status');

    % Stage all changes
    system('git add .');

    % Commit
    commitCmd = sprintf('git commit -m "%s"', msg);
    system(commitCmd);

    % Push
    system('git push');

    % Show status after push
    fprintf('--- Git Status AFTER ---\n');
    system('git status');
end
