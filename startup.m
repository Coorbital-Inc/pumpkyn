function startup()
    rootDir = fileparts(mfilename('fullpath'));
    addpath(rootDir);
    addpath([rootDir,filesep,'tests']);
    addpath([rootDir,filesep,'doc']);
end
