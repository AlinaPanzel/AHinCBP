function ok = push_ahin(msg)
% Quick push for AHinCBP project
    if nargin < 1 || isempty(msg)
        msg = "Update from MATLAB";
    end
    
    project_dir = '/Users/alinapanzel/Desktop/Dartmouth_project/Revisions_Annals/AHinCBP';
    
    original_dir = pwd;
    cd(project_dir);
    
    try
        ok = push_changes(msg);
    catch ME
        cd(original_dir);
        rethrow(ME);
    end
    
    cd(original_dir);
end