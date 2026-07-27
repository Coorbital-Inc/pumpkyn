function [h,globe] = showEarth(jd0,lStar,muStar,hIn)
%% Purpose:
%
%  This routine will properly place the Earth in dimensionless coordinates
%  with the correct scaling and position at [-mu,0,0]
%
%% Inputs:
%
%  jd0                  double              Julian Date (days)
%
%  lStar                double              Characteristic Length (km)
%
%  muStar               double              Mass Ratio of Primaries
%                                           muStar = mu2/(mu1+mu2)
%
%  hIn                  handle              Optional figure or globe handle
%
%% Outputs:
%
%  h                    handle              Handle to current figure
%
%  globe                handle              Handle to Globe graphics
%
%% Revision History:
%  Darin C. Koblick                                              06/18/2026
%  Copyright 2026 Coorbital, Inc.
%% --------------------------- Begin Code Sequence ------------------------

if nargin == 0
          jd0 = pumpkyn.util.juliandate('03/21/2000 00:00:00');
        lStar = 389703.264829278;
       muStar = 0.012150585609624;

    [h,globe] = pumpkyn.cr3bp.showEarth(jd0,lStar,muStar);

    jdVec = jd0 + (0:3600:86400)'./86400;
    view(90,0);

    for tt = 1:numel(jdVec)
        pumpkyn.cr3bp.showEarth( ...
            jdVec(tt),lStar,muStar,globe);

        drawnow;
        pause(0.1);
    end

    return;
end

if ~exist('hIn','var')
    hIn = figure('color',[0 0 0]);
    set(gca(hIn),'color','k');
    axis(gca(hIn),'off','equal');
    hold on;
end
             opts.posOffset = [-muStar,0,0];
                 opts.scale = 1./lStar;
                  opts.type = 'day';
               opts.overlay = true;
                 opts.atmos = false;
            opts.AddShading = false;
                 opts.WGS84 = true;
% If hIn is an existing globe surface, update that globe directly.
% Otherwise, create a new globe in the supplied figure.
if isgraphics(hIn,'surface')
    globe = hIn;
    h = ancestor(globe,'figure');
else
    [h,globe] = pumpkyn.util.earth3D(opts,hIn);
end
% Place this globe under its own transform object. This permits the
% orientation of this Earth to change without affecting other objects.
if isa(globe.Parent,'matlab.graphics.primitive.Transform')
    hTransform = globe.Parent;
else
    hAxes = ancestor(globe,'axes');
    hTransform = hgtransform('Parent',hAxes);
    globe.Parent = hTransform;
end

if ~isempty(jd0)
    % Get three Earth-fixed basis directions expressed in CR3BP frame
    xObs_x = pumpkyn.cr3bp.fromLLA(jd0,[0 0 0],muStar,lStar,1);
    xObs_y = pumpkyn.cr3bp.fromLLA(jd0,[0 pi/2 0],muStar,lStar,1);
    xObs_z = pumpkyn.cr3bp.fromLLA(jd0,[pi/2 0 0],muStar,lStar,1);
    % Convert absolute positions into Earth-centered direction vectors
    rE = [-muStar,0,0];
    ex = xObs_x(1:3) - rE;
    ey = xObs_y(1:3) - rE;
    ez = xObs_z(1:3) - rE;
    ex = ex ./ norm(ex);
    ey = ey ./ norm(ey);
    ez = ez ./ norm(ez);
    % Construct the sampled Earth-to-CR3BP orientation matrix.
    DCME2C = [ex(:),ey(:),ez(:)];
    % Project the sampled matrix onto the nearest orthonormal matrix.
    % This removes numerical scale and shear that hgtransform rejects.
    [U,~,V] = svd(DCME2C);
    DCME2C = U*V';
    % Force the result to be a proper right-handed rotation.
    if det(DCME2C) < 0
        U(:,3) = -U(:,3);
        DCME2C = U*V';
    end
    % Construct an absolute rotation about the center of the Earth.
      TtoEarth = makehgtform('translate',rE);
    TfromEarth = makehgtform('translate',-rE);
             R = eye(4);
    R(1:3,1:3) = DCME2C;
    % Replace the existing transform rather than accumulating rotations.
    hTransform.Matrix = TtoEarth*R*TfromEarth;
end
% Modify only the axes containing this globe.
set(ancestor(globe,'axes'),'clipping','off');
end