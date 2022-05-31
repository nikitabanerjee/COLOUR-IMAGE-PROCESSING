%%Taking a standard colour image and extract its R, G and B components
clear;
clc;
I=imread('flower.jfif');
 Red = I; 
 Red(:,:,2:3)=0; %Removing Green and Blue components
 Green = I;
 Green(:,:,1)=0; %Removing Red components
 Green(:,:,3)=0; %Removing Blue components
 Blue = I;
 Blue(:,:,1:2)=0;  %Removing Green and Red components
 OpImg=Red+Green+Blue; %Adding the images with individual components
 %Display of the images
subplot(3,3,[1,2,3])
imshow(I)
subplot(3,3,4)
imshow(Red)
colormap([[0:1/255:1]', zeros(256,1), zeros(256,1)])
subplot(3,3,5)
imshow(Green)
subplot(3,3,6)
imshow(Blue)
subplot(3,3,[7,8,9])
imshow(OpImg)
%% Finding its (colour image) histogram
clear;
clc;
I=imread('flower.jfif');
 Red = I(:,:,1);
  Green = I(:,:,2);
  Blue = I(:,:,3);
       [yRed, r] = imhist(Red);
       [yGreen, g] = imhist(Green);
       [yBlue, b] = imhist(Blue);
   subplot(2, 3, 1); 
          plot(r, yRed, 'Red')
          xlabel('Histogram of Input image(Red Component) ');
   subplot(2, 3, 2); 
          plot(g, yGreen, 'Green')
          xlabel('Histogram of Input image(Green Component) ');
   subplot(2, 3, 3); 
          plot(b, yBlue, 'Blue')
          xlabel('Histogram of Input image(Blue Component) ');
   subplot(2, 3, [4,5,6]); 
         bar(yRed, 'Red')
         hold on
         bar(yGreen, 'Green')
         hold on
         bar( yBlue, 'Blue')
         xlabel('Colour Histogram of image ');
         % From the above histogram we can interpret that when taken 
         % individually all the colour components behave like 
         % gray scale images and the sum of all three colour componets
         % and their intensity as weights we get the final colour of the
         % image.

%% the chromaticity diagram
clear;
clc;
I=imread('flower.jfif');
A=rgb2xyz(I);
X=A(:,:,1); % X component
Y=A(:,:,2); % Y component
Z=A(:,:,3); % Z component
[row,col]=size(X);
for i=1:row
    for j=1:col
        qx(i,j)=(X(i,j))/(X(i,j)+Y(i,j)+Z(i,j)); % qx component
        qy(i,j)=(Y(i,j))/(X(i,j)+Y(i,j)+Z(i,j)); % qy component
        qz(i,j)=(Z(i,j))/(X(i,j)+Y(i,j)+Z(i,j)); % qz component
    end
end
figure()
plotChromaticity
hold on
plot(qx,qy,'*',MarkerSize=2)

%% Convert the RGB image to HSI, CMY and CMYK model.
clc
clear
i=imread('flower.jfif');
imahsi=rgb2hsv(i); %Converts RGB to HSV format
c1=makecform('srgb2cmyk'); %Converts RGB to CMYK format
I_CMYK=applycform(i,c1);
c2=makecform('srgb2cmy');  %Converts RGB to CMY format
I_CMY=applycform(i,c2);  
figure(),imshow(imahsi)

%% histogram equalization of colour image without using any built-in MATLAB function
clc;
clear;
rgbImage = imread('flower.jfif');
I=double(rgbImage)/255;
R=I(:,:,1);
G=I(:,:,2);
B=I(:,:,3);
A=R+G+B;
%Hue
numi=1/2*((R-G)+(R-B));
denom=((R-G).^2+((R-B).*(G-B))).^0.5;
%To avoid arithamatic exception
H=acosd(numi./(denom+0.001));
%If B>G then H= 360-Theta
H(B>G)=360-H(B>G);
%Normalize to the range [0 1]
H=H/360;
%Saturation
S=1- (3./(sum(I,3)+0.001)).*min(I,[],3);
%Intensity
I=A./3;
HSI=zeros(size(rgbImage));
HSI(:,:,1)=H;
HSI(:,:,2)=S;
HSI(:,:,3)=I;
I=round(I*255);
[r,c]=size(I);
num_of_pixels=r*c;
freq=zeros(512,1);
prob=zeros(512,1);
cdf=zeros(512,1);
probc=zeros(512,1);
HI=uint8(zeros(r,c));
for i=1:r
    for j=1:c
        pixint=I(i,j);
        freq(pixint+1)=freq(pixint+1)+1;
        prob(pixint+1)=freq(pixint+1)/num_of_pixels;
    end
end
sum=0;
no_bins=511;
%The cumulative distribution probability is calculated.
for i=1:size(prob)
    sum=sum+freq(i);
    cdf(i)=sum;
    probc(i)=cdf(i)/num_of_pixels;
    output(i)=round(probc(i)*no_bins);
end
for i=1:r
    for j=1:c
        HI(i,j)=output(I(i,j)+1);
    end
end
I=HI/255;
figure,imshow(rgbImage),title('Original Image')
A=hsv2rgb(HSI);
figure,imshow(A);title('Image after histogram equalization');
