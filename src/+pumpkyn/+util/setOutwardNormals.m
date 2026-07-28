function normals = setOutwardNormals(globe,center)
%% setOutwardNormals
% Set a surface's vertex normals to the outward orientation while
% preserving local mesh slopes such as lunar topography.

if ~isscalar(globe) || ~isgraphics(globe,'surface')
    error( ...
        'pumpkyn:setOutwardNormals:InvalidSurface', ...
        'globe must be a scalar surface graphics handle.');
end

validateattributes( ...
    center,{'numeric'},{'real','finite','vector','numel',3}, ...
    'pumpkyn.util.setOutwardNormals','center');

center = reshape(center,1,3);

x = globe.XData;
y = globe.YData;
z = globe.ZData;

[nx,ny,nz] = surfnorm(x,y,z);

normals = cat(3,nx,ny,nz);
normalMagnitude = sqrt(sum(normals.^2,3));

radial = cat( ...
    3, ...
    x-center(1), ...
    y-center(2), ...
    z-center(3));

radialMagnitude = sqrt(sum(radial.^2,3));

validNormal = ...
    isfinite(normalMagnitude) & normalMagnitude > 0;

for dimension = 1:3
    component = normals(:,:,dimension);
    component(validNormal) = ...
        component(validNormal)./normalMagnitude(validNormal);
    normals(:,:,dimension) = component;
end

validRadial = ...
    isfinite(radialMagnitude) & radialMagnitude > 0;

invalidNormal = ~validNormal & validRadial;

for dimension = 1:3
    component = normals(:,:,dimension);
    radialComponent = radial(:,:,dimension);
    component(invalidNormal) = ...
        radialComponent(invalidNormal)./radialMagnitude(invalidNormal);
    normals(:,:,dimension) = component;
end

alignment = sum(normals.*radial,3);
validAlignment = isfinite(alignment) & validRadial;

if any(validAlignment,'all') && ...
        median(alignment(validAlignment)) < 0

    normals = -normals;
end

globe.VertexNormals = normals;
globe.VertexNormalsMode = 'manual';

end
