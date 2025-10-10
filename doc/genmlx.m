function genmlx()

%% Get all file names with .m:
thisDir = pwd;
fNames = dir(thisDir);
fNames = {fNames.name}';
fNames = fNames(3:end);
[~,fNames,fExt] = fileparts(fNames);
   idx = strcmp(fNames, mfilename) | ~strcmp(fExt,'.m');
fNames(idx) = [];

%% Convert all file names from .m to .mlx:
for tf=1:size(fNames,1)
        %Create an mlx file:
        matlab.internal.liveeditor.openAndSave([fNames{tf},'.m'], ...
                                               [fNames{tf},'.mlx']);
        %Export it to a markdown file:
        export([fNames{tf},'.mlx'], ...
               ['.',filesep,'md',filesep,fNames{tf},'.md']);
end


end