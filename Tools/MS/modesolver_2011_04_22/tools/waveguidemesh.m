function [x,y,xc,yc,nx,ny,eps,varargout] = waveguidemesh(n,h,rh,rw,side,dx,dy);

% This function creates an index mesh for the finite-difference
% mode solver.  The function will accommodate a generalized three
% layer rib waveguide structure.  (Note: channel waveguides can
% also be treated by selecting the parameters appropriately.) 
% 
% USAGE:
% 
% [x,y,xc,yc,nx,ny,eps] = waveguidemesh(n,h,rh,rw,side,dx,dy)
% [x,y,xc,yc,nx,ny,eps,edges] = waveguidemesh(n,h,rh,rw,side,dx,dy)
%
% INPUT
%
% n - indices of refraction for layers in waveguide
% h - height of each layer in waveguide
% rh - height of waveguide feature
% rw - half-width of waveguide
% side - excess space to the right of waveguide
% dx - horizontal grid spacing
% dy - vertical grid spacing
% 
% OUTPUT
% 
% x,y - vectors specifying mesh coordinates
% xc,yc - vectors specifying grid-center coordinates
% nx,ny - size of index mesh
% eps - index mesh (n^2)
% edges - (optional) list of edge coordinates, to be used later
%   with the line() command to plot the waveguide edges
%
% AUTHOR:  Thomas E. Murphy (tem@umd.edu)

% dx = 0.0125, dy = 0.0125 -> resolution = 80 cells/unit length
ih = round(h/dy); % h=[h1,h2,h3]=[2.0,1.3,0.5] -> irh=[ih1,ih2,ih3]=[160,104,40]
irh = round (rh/dy); % rh=1.1 -> irh = 88
irw = round (rw/dx); % rw = 1.0 -> irw = 80
iside = round (side/dx); % side = 1.5 -> iside = 120
nlayers = length(h); % nlayers = 3

% Total span of the grid to be created
nx = irw+iside+1; % nx = irw+iside+1 = 80 + 120 + 1 = 201
ny = sum(ih)+1; % ny = ih1+ih2+ih3+1 = 160 + 104 + 40 + 1 = 305

xc = (1:(nx-1))'*dx - dx/2; % 'center points' for the x direction, in SI units (size=(200,1))
yc = (1:(ny-1))*dy - dy/2; % 'center points' for the y direction, in SI units (size=(1,305))
x = (0:(nx-1))'*dx; % Linearly interpolated segments for X (size=(201,1))
y = (0:(ny-1))*dy; % Linearly interpolated segments for Y (size=(1,305))

% Permittivity matrix -  size = (nx-1,ny-1) = (200, 304)
eps = zeros(nx-1,ny-1);

iy = 1;

for jj = 1:nlayers,
  for i = 1:ih(jj),
	eps(:,iy) = n(jj)^2*ones(nx-1,1);
	iy = iy+1;
  end
end

iy = sum(ih)-ih(nlayers);
for i = 1:irh,
   eps(irw+1:irw+iside,iy) = n(nlayers)^2*ones(iside,1);
   iy = iy-1;
end

nx = length(xc);
ny = length(yc);

if (nargout == 8)
  iyp = cumsum(ih);
  for jj = 1:nlayers-2,
    if (iyp(jj) >= (iyp(nlayers-1)-irh)) % checking if the waveguide dips into another layer
  edges{1,jj} = dx*[0,irw];
    else
      edges{1,jj} = dx*[0,irw+iside];
    end
    edges{2,jj} = dy*[1,1]*iyp(jj);
  end
  jj = nlayers-1;
  edges{1,jj} = dx*[0,irw,irw,irw+iside];
  edges{2,jj} = dy*[iyp(jj),iyp(jj),iyp(jj)-irh,iyp(jj)-irh];
  varargout(1) = {edges};
end