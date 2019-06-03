function [MyV1,MyV2,R] = FFM_Original_V1(im,v1,v2)
addpath('./BM3D');
addpath('./pyramid')
if (~exist('v1','var'))	% alpha -- parameter for shape
    v1 = 1.25;
end
if (~exist('v2','var'))	% alpha -- parameter for shape
    v2 = 1.25;
end

param = {-0.3293    1.1258};
[a, b] = deal(param{:});
f = @(x)exp((1-x.^a).*b);
g = @(I,k)real(I.^(k.^a).*f(k));
% g = @(I,k)real(I.^(k.^a).*f(k));

ratioMax = 7;


[I1, R1] = FractionalEnhence(im,v1,v2);
I1 = tsmooth(I1,0.015,3);
hsv = rgb2hsv(im);
kRatio = min(1./I1,ratioMax);
J = g(hsv(:,:,3), kRatio.^(1/1.5));
hsv(:,:,3) = J;
MyV1 = hsv2rgb(hsv);


YUV = rgb2ycbcr(MyV1);
Y = YUV(:,:,1);
sigma_BM3D = 10;
[~, Y_d] = BM3D(Y,Y,sigma_BM3D,'lc',0);
I_d = ycbcr2rgb(cat(3,Y_d,YUV(:,:,2:3)));
% IB1 = min(power(I1,1/2.2),1);
% IB1 = max(IB1,0);
IB1 = I1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
w1(:,:,1)=IB1;
w1(:,:,2)=1-IB1;

I11(:,:,:,1) = MyV1;
I11(:,:,:,2) = I_d;

r = size(I11,1);
c = size(I11,2);
N = size(I11,4);

pyr1 = gaussian_pyramid(zeros(r,c,3));
nlev = length(pyr1);


for i = 1:N
    % construct pyramid from each input image
    pyrW1 = gaussian_pyramid(w1(:,:,i));
    pyrI = laplacian_pyramid(I11(:,:,:,i));
    
    % blend
    for l = 1:nlev
        W = repmat(pyrW1{l},[1 1 3]);
        pyr1{l} = pyr1{l} + W.*pyrI{l};
    end
end

% reconstruct
I_f = abs(reconstruct_laplacian_pyramid(pyr1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[I2, R2] = FractionalEnhence(I_f,v1,v2);
I2 = tsmooth(I2,0.015,3);
hsv = rgb2hsv(I_f);
kRatio = min(1./I2,ratioMax);
J = g(hsv(:,:,3), kRatio.^(1/1.3));
hsv(:,:,3) = J;
MyV2 = hsv2rgb(hsv);



% I1B = repmat(power(I1,1/2.2).*1.5,[1,1,3]);
% I1B = min(I1B,1);
% MyV3 = I1B.*MyV1+ (1-I1B).*MyV2;
% figure;imshow(MyV1);title('MyV1');
% figure;imshow(MyV2);title('MyV2');
% figure;imshow(MyV3);title('MyV3');

I(:,:,:,1) = im;
I(:,:,:,2) = MyV1;
I(:,:,:,3) = MyV2;

r = size(I,1);
c = size(I,2);
N = size(I,4);

pyr = gaussian_pyramid(zeros(r,c,3));
nlev = length(pyr);

%%%%%%%%%%%%%%%% use gaussian%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% original_illumination = tsmooth( max(im,[],3),0.015,3);
% iteration1_illumination = tsmooth( max(MyV1,[],3),0.015,3);
% iteration2_illumination = tsmooth( max(MyV2,[],3),0.015,3);
% w(:,:,1)= exp(-(original_illumination-0.5).^2/0.125);
% w(:,:,2)= exp(-(iteration1_illumination-0.5).^2/0.125);
% w(:,:,3)= exp(-(iteration2_illumination-0.5).^2/0.125);
% w = w./repmat(sum(w,3),[1 1 3]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

original_illumination = tsmooth( max(im,[],3),0.015,3);
iteration1_illumination = tsmooth( max(MyV1,[],3),0.015,3);
iteration2_illumination = tsmooth( max(MyV2,[],3),0.015,3);
w(:,:,1)= 1./(1+ exp(-4*original_illumination+2));
w(:,:,2)= exp(-(iteration1_illumination-0.5).^2/0.125);
w(:,:,3)=  1./(1+ exp(4*iteration2_illumination-2)).*0.8;
w = w./repmat(sum(w,3),[1 1 3]);


for i = 1:N
    % construct pyramid from each input image
    pyrW = gaussian_pyramid(w(:,:,i));
    pyrI = laplacian_pyramid(I(:,:,:,i));
    
    % blend
    for l = 1:nlev
        W = repmat(pyrW{l},[1 1 3]);
        pyr{l} = pyr{l} + W.*pyrI{l};
    end
end

% reconstruct
R = reconstruct_laplacian_pyramid(pyr);


