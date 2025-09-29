function ok = push_changes(msg, varargin)
%PUSH_CHANGES  Stage, commit, rebase (if needed) and push the current branch.
%
% Usage:
%   push_changes("added plot function")
%   ok = push_changes("msg", 'AutoStash', true, 'Remote', 'origin')

    p = inputParser;
    p.addParameter('AutoStash', true, @(x)islogical(x)||ismember(x,[0 1]));
    p.addParameter('Remote', 'origin', @(x)ischar(x)||isstring(x));
    p.parse(varargin{:});
    autostash = p.Results.AutoStash;
    remote    = char(p.Results.Remote);

    if nargin < 1 || isempty(msg)
        msg = "Update from MATLAB";
    end
    msg = char(msg);

    ok = false;

    fprintf('--- Git Status BEFORE ---\n');
    sys('git status');  %#ok<*NASGU>

    % current branch
    [st, branch] = sys('git rev-parse --abbrev-ref HEAD', true);
    if st~=0
        warn('Could not determine current branch.');
        return;
    end
    branch = strtrim(branch);
    fprintf('Current branch: %s\n', branch);

    % ensure upstream tracking exists
    [hasUpstream] = sys(sprintf('git rev-parse --abbrev-ref %s@{upstream}', branch));
    if hasUpstream ~= 0
        sys(sprintf('git branch --set-upstream-to %s/%s %s', remote, branch, branch));
    end

    % Stage all
    sys('git add -A');

    % Avoid empty commit
    if sys('git diff --cached --quiet') == 0
        fprintf('No staged changes; skipping commit.\n');
    else
        if sys(sprintf('git commit -m "%s"', escape_quotes(msg))) ~= 0
            warn('Commit failed.');
            return;
        end
    end

    % First push try
    pushCmd = sprintf('git push %s %s', remote, branch);
    if sys(pushCmd) == 0
        fprintf('Push succeeded.\n');
        done(); ok = true; return;
    end

    % Push rejected -> fetch + rebase
    fprintf('Push rejected. Fetching and rebasing onto %s/%s ...\n', remote, branch);
    if sys(sprintf('git fetch %s', remote)) ~= 0
        warn('git fetch failed.'); return;
    end

    if autostash
        rb = sprintf('git pull --rebase --autostash %s %s', remote, branch);
    else
        rb = sprintf('git pull --rebase %s %s', remote, branch);
    end
    if sys(rb) ~= 0
        warn('Rebase failed. Resolve conflicts, then: git add <files> ; git rebase --continue');
        sys('git status');
        return;
    end

    % Try push again
    if sys(pushCmd) == 0
        fprintf('Push after rebase succeeded.\n');
        done(); ok = true; return;
    else
        warn('Push still failing (branch protection/CI?).');
    end

    % ---- helpers ----
    function [status, output] = sys(cmd, capture)
        if nargin<2, capture = false; end
        if capture
            [status, output] = system(cmd);
        else
            status = system(cmd);
            if nargout > 1, output = ""; end
        end
    end

    function warn(msg)
        fprintf(2, '[WARN] %s\n', msg);
    end

    function txt = escape_quotes(txt)
        txt = strrep(txt, '"', '\"');
    end

    function done()
        fprintf('--- Git Status AFTER ---\n');
        sys('git status');
    end
end
