function R = Rx(theta)
%% Purpose:
%  Generate a 3-D rotation matrix to rotate another DCM or vector about the
%  x-axis.
%
%% Inputs:
%  theta                [N x 1]                         Rotation Angle in
%                                                       degrees.
%
%% Outputs:
%  R                    [3 x 3 x N]                     Rotation Matrix
%
%% Revision History:
%  Darin C. Koblick                             (c) 2013
%% ---------------------- Begin Code Sequence -----------------------------
theta = theta(:);
R = zeros(3,3,length(theta));
R(1,1,:) = 1;
R(2,2,:) = cosd(theta); R(2,3,:) = -sind(theta);
R(3,2,:) = sind(theta); R(3,3,:) = cosd(theta);
end