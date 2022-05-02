function area = voronoidens(kx,ky)
% function area = voronoidens(kx,ky);
% input: kx, ky are k-space trajectories
% output: area of cells for each point

% combine kx andd ky trajectories
kxy = [kx(:),ky(:)];

% create complex k-space locations and find unique points
kxy_complex = kx + 1i*ky;
[C, ia, ic] = unique(kxy_complex);
nunique_pts = numel(C);

% returns vertices and cells of
% voronoi diagram
[V,C] = voronoin(kxy(ia,:));

% unpack cell array, compute area of each ploygon
areas_unique = zeros(nunique_pts,1);
for jj = 1:nunique_pts
    x = V(C{jj},1);
    y = V(C{jj},2);
    areas_unique(jj) = polyarea(x,y);
end

area = areas_unique(ic);

end