function R = Rz(theta)
%% Purpose:
%  Generate a 3-D rotation matrix to rotate another DCM or vector about the
%  z-axis.
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
R(1,1,:) = cosd(theta);  R(1,2,:) = -sind(theta);
R(2,1,:) = sind(theta);  R(2,2,:) =  cosd(theta); 
R(3,3,:) = 1;
end