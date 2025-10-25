function startup()
    rootDir = fileparts(mfilename('fullpath'));
    addpath([rootDir,filesep,'src']);
    addpath([rootDir,filesep,'tests']);
    addpath([rootDir,filesep,'doc']);
end
