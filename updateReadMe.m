function updateReadMe()
%% Update .mlx links to .md links in README or other .md files
% This script finds all instances of (./doc/*.mlx) and replaces them with
% (./doc/md/*.md)

% Specify the markdown file to modify
mdFile = 'README.md';  % <-- change this to your target .md file

% Read file content
txt = fileread(mdFile);

% Regular expression: look for (./doc/filename.mlx)
pattern = '\(./doc/([^)]+)\.mlx\)';

% Replacement: change to (./doc/md/filename.md)
replacement = '(./doc/md/$1.md)';

% Apply replacement globally
newTxt = regexprep(txt, pattern, replacement);

% Write updated text back to file
fid = fopen(mdFile, 'w');
fwrite(fid, newTxt);
fclose(fid);

disp('✅ All .mlx links updated to .md links successfully.');
end