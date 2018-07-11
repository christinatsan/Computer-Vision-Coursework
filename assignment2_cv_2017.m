InputImage = 'IDPicture.bmp';
A = imread(InputImage);

%  QUESTION 1

%calculate intensity image
I = 0.299 * A(:,:,1) + 0.587 * A(:,:,2) + 0.114 * A(:,:,3);

%Generate histogram  - GO BACK AND FIX INDEXING
hist = zeros([255,1]);
for i=1:size(I,1)
    for j=1:size(I,2)
        hist(I(i,j)) = hist(I(i,j)) + 1;
    end
end

figure,histogram = bar(hist,0.2)
title('Histogram')
xlim([0 255])

%Perform histogram operations

% 1. Linear Stretch
% Use the function g(f) = [(f-s1)(t2-t1)]/[(s2-s1) + t2] to stretch the
% range of the image from [s1,s2] to [t1,t2].

L = double(I) % convert to double;
Output = ((L-20)*255)/210;

for i=1:size(Output,1)
    for j=1:size(Output,2)
        if Output(i,j) < 0
            Output(i,j) = 0;
        elseif Output(i,j) > 255
            Output(i,j) = 255;
        end
    end
end
               

Output = uint8(Output); % Convert to uint8

figure,imshow(Output)
title('Linear Stretch');

%Generate histogram 
hist = zeros([256,1]);
for i=1:250
    for j=1:250 
        hist(Output(i,j)+1) = hist(Output(i,j)+1) + 1;
    end
end

figure,histogram = bar(hist,0.2)
xlim([0 255])
title('Linear Stretch Histogram');

% 2. Thresholding
% Transverse through each pixel and if pixel is greater than threshold
% value of 100 make pixel value equal to 255 and if pixel value is less
% than threshold value of 100 make pixel value equal to 0

binary_image = I 

for i=1:size(binary_image,1)
    for j=1:size(binary_image,2)
        if binary_image(i,j) > 100
            binary_image(i,j) = 255;
        else 
            binary_image(i,j) = 0;
        end
    end  
end

%Generate histogram 
hist = zeros([256,1]);
for i=1:size(binary_image,1)
    for j=1:size(binary_image,2)
        hist(binary_image(i,j)+1) = hist(binary_image(i,j)+1) + 1;
    end
end

figure,imshow(binary_image)
title('Threshold Image');
figure,histogram = bar(hist,0.2)
xlim([0 255])
title('Threshold Histogram');

%%

% 3. Histogram Equalization

Img=imhist(I)

%Divide the result by number of pixels
Output=Img/numel(I);
%Calculate the Cumlative sum
CSum=cumsum(Output);

%Perform the transformation S=T(R) where S and R in the range [ 0 1]
equalized_image=CSum(I+1);
equalized_image=uint8(equalized_image*bin);

for i=1:size(equalized_image,1)
    for j=1:size(equalized_image,2)
        if equalized_image(i,j) < 0
            equalized_image(i,j) = 0;
        elseif equalized_image(i,j) > 255
            equalized_image(i,j) = 255;
        end
    end
end

figure,imshow(equalized_image);

title('Histogram Equalization');

%Generate histogram 
hist = zeros([256,1]);
for i=1:size(equalized_image,1)
    for j=1:size(equalized_image,2)
        hist(equalized_image(i,j)+1) = hist(equalized_image(i,j)+1) + 1;
    end
end

figure,histogram = bar(hist,0.2)
xlim([0 255])
title('Equalized Histogram');

%%
% QUESTION 2

% 1. Combined Sobel Operator

% Kernel: x: [-2,-1,-2;0,0,0;1,2,1] y: [-1,0,1;-2,0,-2;-1,0,1]

C=double(I); % Convert to double


for i=1:size(C,1)-2
    for j=1:size(C,2)-2
        %Find gradient for x direction
        X=((2*C(i+2,j+1)+C(i+2,j)+C(i+2,j+2))-(2*C(i,j+1)+C(i,j)+C(i,j+2)));
        %Find gradient for y direction
        Y=((2*C(i+1,j+2)+C(i,j+2)+C(i+2,j+2))-(2*C(i+1,j)+C(i,j)+C(i+2,j)));
      
        sobel_image(i,j)=sqrt(X.^2+Y.^2); 
      
    end
end
sobel_image = uint8(sobel_image); 

figure,imshow(sobel_image); 
title('Sobel gradient');

%Generate histogram 
hist = zeros([256,1]);
for i=1:size(sobel_image,1)
    for j=1:size(sobel_image,2)
        hist(sobel_image(i,j)+1) = hist(sobel_image(i,j)+1) + 1;
    end
end

figure,histogram = bar(hist,0.2)
xlim([0 255])
title('Sobel Histogram');
%%

% 2. Horizontal Sobel Operator

C=double(I); % Convert image to double


for i=1:size(C,1)-2
    for j=1:size(C,2)-2
        %Find gradient for x direction
        X=((2*C(i+2,j+1)+C(i+2,j)+C(i+2,j+2))-(2*C(i,j+1)+C(i,j)+C(i,j+2)));
       
       hsobel_image(i,j)=sqrt(X.^2);
      
    end
end

hsobel_image = uint8(hsobel_image);

figure,imshow(hsobel_image);
title('Horizontal Sobel gradient');

%%

% 3. Vertical Sobel Gradient

C=double(I); % Convert image to double


for i=1:size(C,1)-2
    for j=1:size(C,2)-2
        
       %Find gradient for y direction
        Y=((2*C(i+1,j+2)+C(i,j+2)+C(i+2,j+2))-(2*C(i+1,j)+C(i,j)+C(i+2,j)));
      
        vsobel_image(i,j)=sqrt(Y.^2);
      
    end
end
vsobel_image = uint8(vsobel_image);

figure,imshow(vsobel_image); 
title('Vertical Sobel gradient');


%%

% 4. 1x2 Combined edge 

% Kernel: x: [-1,1] y: [-1;1]

C = double(I); % Convert image to a double

for i=1:size(C,1)-2
    for j=1:size(C,2)-2
        X = (C(i,j+1) - C(i,j));
        Y = (C(i+1,j) - C(i,j));
        
        onetwo_image(i,j) = sqrt(X.^2+Y.^2);
    end
end

onetwo_image = uint8(onetwo_image);

figure,imshow(onetwo_image); 
title('1x2');

%Generate histogram 
hist = zeros([256,1]);
for i=1:size(onetwo_image,1)
    for j=1:size(onetwo_image,2)
        hist(onetwo_image(i,j)+1) = hist(onetwo_image(i,j)+1) + 1;
    end
end

figure,histogram = bar(hist,0.2)
xlim([0 255])
title('1x2 Histogram');

%%

% 5. 1x2 Horizontal edge 

C = double(I); % Convert to a double

for i=1:size(C,1)-2
    for j=1:size(C,2)-2
        % Find gradient for x direction
        X = (C(i+1,j) - C(i,j));
        
        h_onetwo(i,j) = sqrt(X.^2);
    end
end

h_onetwo = uint8(h_onetwo);

figure,imshow(h_onetwo);
title('1x2 Horizontal gradient');

%%

% 6. 1x2 Vertical edge

C = double(I); % Convert to a double

for i=1:size(C,1)-2
    for j=1:size(C,2)-2
       %Find gradient for y direction
        Y = (C(i,j+1) - C(i,j));
        
        v_onetwo(i,j) = sqrt(Y.^2);
    end
end

v_onetwo = uint8(v_onetwo);

figure,imshow(v_onetwo); 
title('1x2 Vertical gradient');

%%

% QUESTION 4

% Sobel Edge Map

% Use histogram to choose threshold value

thresh=120;
sobel_thresh = sobel_image;

for i=1:size(sobel_image,1)
    for j=1:size(sobel_image,2)
        if sobel_image(i,j) > thresh
            sobel_thresh(i,j) = 255;
        else
            sobel_thresh(i,j) = 0;
        end
    end
end

figure,imshow(sobel_thresh)
title('Sobel Threshold');

%Generate histogram 
hist = zeros([256,1]);
for i=1:size(sobel_thresh,1)
    for j=1:size(sobel_thresh,2)
        hist(sobel_thresh(i,j)+1) = hist(sobel_thresh(i,j)+1) + 1;
    end
end

figure,histogram = bar(hist,0.2)
xlim([0 255])
title('Sobel Threshold Histogram');




%%

% 1x2 Edge Map

% Use histogram to choose threshold value

thresh = 30;
onetwo_thresh = onetwo_image;

for i=1:size(onetwo_image,1)
    for j=1:size(onetwo_image,2)
        if onetwo_image(i,j) > thresh
            onetwo_thresh(i,j) = 255;
        else
            onetwo_thresh(i,j) = 0;
        end
    end
end

figure,imshow(onetwo_thresh)
title('1x2 Threshold');

%Generate histogram 
hist = zeros([256,1]);
for i=1:size(onetwo_thresh,1)
    for j=1:size(onetwo_thresh,2)
        hist(onetwo_thresh(i,j)+1) = hist(onetwo_thresh(i,j)+1) + 1;
    end
end

figure,histogram = bar(hist,0.2)
xlim([0 255])
title('1x2 Threshold Histogram');

%%

% Subtract gradient maps

subtract_gradient = sobel_image - onetwo_image;
figure,imshow(subtract_gradient)
title('Sobel - 1x2');

%% 

% Color edge detector

C = double(A); % Convert to double

% Apply sobel operator to each color band

R = C(:,:,1); % Red color band
for i=1:size(C,1)-2
    for j=1:size(C,2)-2
        %Find gradient in x direction
        X=((2*R(i+2,j+1)+R(i+2,j)+R(i+2,j+2))-(2*R(i,j+1)+R(i,j)+R(i,j+2)));
        %Find gradient in y direction
        Y=((2*R(i+1,j+2)+R(i,j+2)+R(i+2,j+2))-(2*R(i+1,j)+R(i,j)+R(i+2,j)));
      
        sobel_image_red(i,j)=sqrt(X.^2+Y.^2);
        
    end
end

sobel_image_red = uint8(sobel_image_red);

G = C(:,:,2); % Green color band
for i=1:size(C,1)-2
    for j=1:size(C,2)-2
        %Find gradient in x direction
        X=((2*G(i+2,j+1)+G(i+2,j)+G(i+2,j+2))-(2*G(i,j+1)+G(i,j)+G(i,j+2)));
        %Find gradient in y direction
        Y=((2*G(i+1,j+2)+G(i,j+2)+G(i+2,j+2))-(2*G(i+1,j)+G(i,j)+G(i+2,j)));
     
        sobel_image_green(i,j)=sqrt(X.^2+Y.^2);
        
    end
end

sobel_image_green = uint8(sobel_image_green);

B = C(:,:,3); % Blue color band
for i=1:size(C,1)-2
    for j=1:size(C,2)-2
        %Find gradient in x direction
        X=((2*B(i+2,j+1)+B(i+2,j)+B(i+2,j+2))-(2*B(i,j+1)+B(i,j)+B(i,j+2)));
        %Find gradient in y direction
        Y=((2*B(i+1,j+2)+B(i,j+2)+B(i+2,j+2))-(2*B(i+1,j)+B(i,j)+B(i+2,j)));
     
        sobel_image_blue(i,j)=sqrt(X.^2+Y.^2);
        
    end
end

sobel_image_blue = uint8(sobel_image_blue);

%concatenate 3 color band images
rgb_sobel = cat(3,sobel_image_red,sobel_image_green,sobel_image_blue);

figure,imshow(rgb_sobel);
title('Color Sobel');

%Intensity Image of rgb sobel

rgb_sobel_I = 0.299 * rgb_sobel(:,:,1) + 0.587 * rgb_sobel(:,:,2) + 0.114 * rgb_sobel(:,:,3);

figure,imshow(rgb_sobel_I);
title('RGB Sobel Intensity');

%Generate histogram 
hist = zeros([256,1]);
for i=1:size(rgb_sobel_I,1)
    for j=1:size(rgb_sobel_I,2)
        hist(rgb_sobel_I(i,j)+1) = hist(rgb_sobel_I(i,j)+1) + 1;
    end
end

figure,histogram = bar(hist,0.2)
xlim([0 255])
title('RGB Sobel Intensity Histogram');

%Threshold of rgb sobel intensity image

thresh = 120;

rgbsobel_I_thresh = rgb_sobel_I;

for i=1:size(rgb_sobel_I,1)
    for j=1:size(rgb_sobel_I,2)
        if rgb_sobel_I(i,j) > thresh
            rgbsobel_I_thresh(i,j) = 255;
        else
            rgbsobel_I_thresh(i,j) = 0;
        end
    end
end

figure,imshow(rgbsobel_I_thresh);
title('RGB Sobel Intensity Threshold');

%Threshold of rgb sobel

thresh = 120;

rgbsobel_thresh = rgb_sobel;

for i=1:size(rgb_sobel,1)
    for j=1:size(rgb_sobel,2)
        for k=1:size(rgb_sobel,3)

            if rgb_sobel(i,j,k) > thresh
                rgb_sobel(i,j,k) = 255;
            else
                rgb_sobel(i,j,k) = 0;
            end
        end
    end
end

figure,imshow(rgb_sobel);
title('RGB Sobel Threshold');
