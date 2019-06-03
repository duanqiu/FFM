function [ I, R] = FractionalEnhence2( src,V1,V2, debug ,alpha, beta, lambda, vareps, r, r0, K)
addpath('./utils');
if (~exist('V1','var'))	% alpha -- parameter for shape
    V1 =2.22;
end
if (~exist('V2','var'))	% alpha -- parameter for shape
    V2 = 2.11;
end
if (~exist('alpha','var'))	% alpha -- parameter for shape
    alpha = 0.001;
end
if (~exist('beta','var'))	% beta -- parameter for texture
    beta = 0.005;
end
if (~exist('lambda','var'))	% lambda -- parameter for illumination
    lambda = 0.25;
end
if (~exist('vareps','var')) % vareps -- stopping parameter
    vareps = 0.01;
end
if (~exist('K','var'))      % K -- maximum iterations
    K = 8;
end
if (~exist('r','var'))      % r -- the size of Omega in Eq.(3)
    r = 3;
end
if (~exist('r0','var'))     % r0 -- the size of Omega in Eq.(7)
    r0 = 5;
end
if (~exist('debug','var'))  % debug -- set debug/release
    debug = true;
end
r = (r-1)/2;
r0 = (r0-1)/2;
eps=0.0001;
if size(src,3)==1
    S = src;
    gray = src;
else
    hsv = rgb2hsv(src);
    S = hsv(:,:,3);
    gray = rgb2gray(src);
end

B = convMax(single(S),r0);                              % Eq.(7)
B = guidedfilter(gray, B, 20, eps);                     % use guided filer to refine bright channel
%B = tsmooth(B,0.015,3);                                 % use tsmooth filer to refine bright channel
I=S;                                                    % initialize I_0
R=ones(size(S));                                        % initialize R_0
if debug == true
    fprintf('-- Stop iteration until eplison < %02f or K > %d\n', vareps, K);
end

[h, w] = size(S);
hw = h*w;
[MY,MX] = FractionalMatrix(hw,V1,h);
VV1  =  FractionalProduct(V1);
fprintf('\t========== %02f½×¾ØÕóÒÑ¹¹½¨==========\n', V1);
for iter = 1:K
    preI=I;
    preR=R;

    I=S./R;
    Ix = diff(I,1,2); Ix = padarray(Ix, [0 1], 'post');
    Iy = diff(I,1,1); Iy = padarray(Iy, [1 0], 'post');
    
    ksize = bitor(round(5*r),1);
    g = fspecial('gaussian', [1,ksize], r); 
    fbin = conv2(I,g,'same');
    fbin = conv2(fbin,g','same');     
    gfx = diff(fbin,1,2);
    gfx = padarray(gfx, [0 1], 'post');
    gfy = diff(fbin,1,1);
    gfy = padarray(gfy, [1 0], 'post');
  


    ux = max(abs(gfx.*Ix),eps).^(-1);                     % ux in Eq.(11)
    uy = max(abs(gfy.*Iy),eps).^(-1);                     % uy in Eq.(11)
 
    ux(:,end) = 0;
    uy(end,:) = 0;
    
    I = solveLinearSystem(S, R, ux, uy, alpha, B, lambda);  % Eq.(12)
    eplisonI = norm(I-preI, 'fro')/norm(preI, 'fro');       % iterative error of I
    

    R=S./I;

    fVV1 = flip(VV1);
    DXR1 = conv2(R,VV1,'same');
    DXR2 = conv2(R,fVV1,'same');
    DXR = max(DXR1,DXR2);
    vx = power(max(abs(DXR),eps),2-V2);
    vx = vx .*beta;
    
    DYR1 =  conv2(R',VV1,'same')';
    DYR2 = conv2(R',fVV1,'same')';
    DYR = max(DYR1,DYR2);
    vy = power(max(abs(DYR),eps),2-V2);
    vy = vy .*beta;
    
    IF = power(I,2);
    IFM = spdiags(IF(:),0,hw,hw);
    Matrix = IFM + MX'*spdiags(vx(:),0,hw,hw)*MX +MY'*spdiags(vy(:),0,hw,hw)*MY;
	
	
    RT = (I).*(S);
    [dst,~] = minres(Matrix,RT(:), 0.01, 40);
    R = reshape(abs(dst), h, w);
    eplisonR = norm(R-preR, 'fro')/norm(preR, 'fro');
	
    if debug == true
        fprintf('Iter #%d : eplisonI = %f; eplisonR = %f\n', iter, eplisonI, eplisonR);
    end
    if(eplisonI<vareps||eplisonR<vareps)
        %        I=S./R;
        break;
    end
end
end


