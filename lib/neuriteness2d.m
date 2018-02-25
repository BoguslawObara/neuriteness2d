function [imn,vx,vy] = neuriteness2d(im,sigma,gamma)
%%  neuriteness2d - neuriteness vessel enhancement filtering
%   
%   REFERENCE:
%       E. Meijering et al. 
%       Design and validation of a tool for neurite tracing and analysis 
%       in fluorescence microscopy images
%       Cytometry Part A, 58A, 167-176, 2004
%
%   INPUT:
%       im      - 2D gray image
%       sigma   - Gaussian kernel sigma
%       gamma   - parameter
%
%   OUTPUT:
%       imv     - neuriteness
%
%   AUTHOR:
%       Boguslaw Obara
%

%% default parameters
if isempty(sigma);  sigma = 1;  end
if isempty(gamma);  gamma = 2;  end

%% normalize
im = double(im); im = (im - min(im(:))) / (max(im(:)) - min(im(:)));

%% convert image to grey-scale
im = im2uint8(im); % I assume that Vesselness was defined for grey-scale images

%% convert image to double
im = double(im);

%% second derivatives - Hessian
[hxx,hxy,hyy] = hessian2d(im,sigma);

%% normalized derivative - scale
hxx = power(sigma,gamma)*hxx;
hxy = power(sigma,gamma)*hxy;
hyy = power(sigma,gamma)*hyy;

%% eigen values and vectors
[l1,l2,v1,v2,v3,v4] = eigen2d_m(hxx,hxy,hxy,hyy);
l2(l2==0) = eps;
    
%% modified Hessian
alfa = -1/3;
l1 = l1 + alfa*l2;
l2 = l2 + alfa*l1;

%% sort l1s > l2s
index = abs(l1)>abs(l2); %<
l1s = l1; 
l2s = l2;
l1s(index) = l2(index);
l2s(index) = l1(index);
l1 = l1s;
l2 = l2s;
l2(l2==0) = eps;

v3s = v3;
v4s = v4;
v3s(index) = v1(index);
v4s(index) = v2(index);
vy = v4s; 
vx = v3s;

%% neuriteness
lmin = min(l1s(:));
imn = zeros(size(l1));
imn(l1s<0) = l1s(l1s<0)./lmin;

end