function packagePumpkynToolbox()
%% Purpose:
%
%  This routine will re-run the mlx -> md scripts and other
%  miscellanious dependencies to prepare and package up the pumpkyn toolbox
%  for distribution.
%
%% Revision History:
%  Darin C. Koblick                                         (c) 10/25/2025
%  Coorbital Inc.
%% ---------------------- Begin Code Sequence -----------------------------

%% Generate MLX and Markdown files:
fprintf(1,'%s\n','Generating MLX and MD Files');
genmlx();

%% Update Links in MD to other MD:
fprintf(1,'%s\n','Updating hyperlinks in MD Files');
updateReadMe()

%% Package up toolbox:
fprintf(1,'%s\n','Packaging MATLAB Toolbox');
            uuid = "pumpkyn";
   toolboxFolder = "C:\GitHub\pumpkyn\";
            opts = matlab.addons.toolbox.ToolboxOptions(toolboxFolder, uuid);
opts.ToolboxName = "Pumpkyn Toolbox";
opts.Description = "Toolbox for exploring tulip-shaped orbits in the CR3BP";
opts.Summary = append("Designed as both a training resource and a preliminary mission",  ...
               "design toolkit, pumpkyn makes  three-body dynamics accessible ", ...
               "through a curated set of routines for orbit continuation ", ...
               "stability assessment, invariant manifold generation, ", ...
               "station-keeping cost estimation, and Earth-centered/J2000", ...
               "frame conversion. Its tutorial scripts and examples allow ", ...
               "students and researchers to build intuition for multi-body ", ...
               "dynamics, while its engineering-grade functions give mission", ...
               "designers the ability to quickly explore architectures and ", ...
               "trade studies in the cislunar regime.");
                opts.AuthorName = "coorbital inc";
             opts.AuthorCompany = "coorbital inc";
               opts.AuthorEmail = "info@coorbital.com";
          opts.ToolboxImageFile = "logo.png";
         opts.ToolboxMatlabPath = {'src','tests','doc'};
opts.ToolboxGettingStartedGuide ="ReadMe.mlx";
              opts.ToolboxFiles = {'src','tests','doc','ReadMe.mlx', 'startup.m', 'LICENSE'};
      opts.MinimumMatlabRelease = "";
      opts.MaximumMatlabRelease = "";
                opts.OutputFile = 'C:\GitHub\pumpkyn\pumpkyn.mltbx';
      matlab.addons.toolbox.packageToolbox(opts);
end