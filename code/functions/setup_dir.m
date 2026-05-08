function path = setup_dir(path, today)
% SETUP_DIR - Helper function to create a directory if it does not exist yet
% Optional functionality: create a subfolder in the specified path with today's date
% 
% Inputs:
%   path (ch): The full path to the folder
%   today (bool): Decision whether to create a subfolder with today's date
%
% Outputs: 
%   path (ch): The full path to the folder, containing today's subfolder if indicated

% Optional: in the indicated folder path, create a date-specific subfolder
% This is useful when analyses are run multiple times, so as to not overwrite temporary outcomes
if today == true
    % Today's date as string
    today_str = datestr(now, 'yyyy-mm-dd');
    path = fullfile(path, today_str); % add that to the path
end

if ~exist(path, 'dir')
    mkdir(path);
    fprintf("\n✅ Created new directory: %s\n", path);
else
    fprintf("\n✅ This directory exists: %s\n", path);
end
end