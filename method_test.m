
img_path ='./pic/9.png';

im=imread(img_path);



im = im2double(im);



%[iteration1,iteration2,Synthesis] = FFM_Original_V1(im,2.2,3.8);
%[iteration1,iteration2,Synthesis] = FFM_Original_V2(im,2.2,3.8);
[iteration1,iteration2,Synthesis] = FFM_SpeedUp(im,2.2,3.8);


figure;imshow(iteration1);title('iteration1');
figure;imshow(iteration2);title('iteration2');
figure;imshow(Synthesis);title('Synthesis');

