%===================================================
% Computer Vision Programming Assignment 1
% @Zhigang Zhu, 2003-2009
% City College of New York
%===================================================

% ---------------- Step 1 ------------------------
% Read in an image, get information
% type help imread for more information

InputImage = 'IDPicture.bmp'; 
%OutputImage1 = 'IDPicture_bw.bmp';

C1 = imread(InputImage);
[ROWS COLS CHANNELS] = size(C1);

% ---------------- Step 2 ------------------------
% If you want to display the three separate bands
% with the color image in one window, here is 
% what you need to do
% Basically you generate three "color" images
% using the three bands respectively
% and then use [] operator to concatenate the four images
% the orignal color, R band, G band and B band

% First, generate a blank image. Using "uinit8" will 
% give you an image of 8 bits for each pixel in each channel
% Since the Matlab will generate everything as double by default
CR1 =uint8(zeros(ROWS, COLS, CHANNELS));

% Note how to put the Red band of the color image C1 into 
% each band of the three-band grayscale image CR1
for band = 1 : CHANNELS,
    CR1(:,:,band) = (C1(:,:,1));
end

% Do the same thing for G
CG1 =uint8(zeros(ROWS, COLS, CHANNELS));
for band = 1 : CHANNELS,
    CG1(:,:,band) = (C1(:,:,2));
end

% and for B
CB1 =uint8(zeros(ROWS, COLS, CHANNELS));
for band = 1 : CHANNELS,
    CB1(:,:,band) = (C1(:,:,3));
end

% Whenever you use figure, you generate a new figure window 
No1 = figure;  % Figure No. 1

%This is what I mean by concatenation
disimg = [C1, CR1;CG1, CB1]; 

% Then "image" will do the display for you!
image(disimg);

% ---------------- Step 3 ------------------------
% Now we can calculate its intensity image from 
% the color image. Don't forget to use "uint8" to 
% covert the double results to unsigned 8-bit integers

I1    = uint8(round(sum(C1,3)/3));

% You can definitely display the black-white (grayscale)
% image directly without turn it into a three-band thing,
% which is a waste of memeory space

No2 = figure;  % Figure No. 2
image(I1);


% If you just stop your program here, you will see a 
% false color image since the system need a colormap to 
% display a 8-bit image  correctly. 
% The above display uses a default color map
% which is not correct. It is beautiful, though

% ---------------- Step 4 ------------------------
% So we need to generate a color map for the grayscale
% I think Matlab should have a function to do this,
% but I am going to do it myself anyway.

% Colormap is a 256 entry table, each index has three entries 
% indicating the three color components of the index

MAP =zeros(256, 3);

% For a gray scale C[i] = (i, i, i)
% But Matlab use color value from 0 to 1 
% so I scale 0-255 into 0-1 (and note 
% that I do not use "unit8" for MAP

for i = 1 : 256,  % a comma means pause 
    for band = 1:CHANNELS,
        MAP(i,band) = (i-1)/255;
    end 
end

%call colormap to enfore the MAP
colormap(MAP);

% I forgot to mention one thing: the index of Matlab starts from
% 1 instead 0.

% Is it correct this time? Remember the color table is 
% enforced for the current one, which is  the one we 
% just displayed.

% You can test if I am right by try to display the 
% intensity image again:

No3 = figure; % Figure No. 3
image(I1);


% See???
% You can actually check the color map using 
% the edit menu of each figure window

% ---------------- Step 5 ------------------------
% Use imwrite save any image
% check out image formats supported by Matlab
% by typing "help imwrite
% imwrite(I1, OutputImage1, 'BMP');


% ---------------- Step 6 and ... ------------------------
% Students need to do the rest of the jobs from c to g.
% Write code and comments - turn it in both in hard copies and 
% soft copies (electronically)


% ------------------ Part C - Display Intensity Image --------------------
intensity_image = 0.299 * C1(:,:,1) + 0.587 * C1(:,:,2) + 0.114 * C1(:,:,3);
figure,imshow(intensity_image)
title('Intensity Image using NTSC')

 
% ------------------ Part D - Quantize Intensity Image -------------------
 
% Convert the image to a double and divide by 255 so the range of values 
% are [0,1]
% Multiply by k-level you want to quantize by 
% Convert the image to uint8
% To view convert to double and divide by k-level
 
level = 4;
quantized_image1 = double(intensity_image) / 255;
quantized_image1 = uint8(quantized_image1 * level);
quantized_image1 = double(quantized_image1) / level;
 
level = 16;
quantized_image2 = double(intensity_image) / 255;
quantized_image2 = uint8(quantized_image2 * level);
quantized_image2 = double(quantized_image2) / level;
 
level = 32;
quantized_image3 = double(intensity_image) / 255;
quantized_image3 = uint8(quantized_image3 * level);
quantized_image3 = double(quantized_image3) / level;
 
level = 64;
quantized_image4 = double(intensity_image) / 255;
quantized_image4 = uint8(quantized_image4 * level);
quantized_image4 = double(quantized_image4) / level;

%concatenate images
quant_img = [quantized_image1, quantized_image2; 
             quantized_image3, quantized_image4]; 
         
%display figure         
figure,imshow(quant_img); 
title('Quantized Intensity Image')

% ------------------- Part E - Quantize Three-band color image -----------------

level = 2;
quantized_image5 = double(C1) / 255;
quantized_image5 = uint8(quantized_image5 * level);
quantized_image5 = double(quantized_image5) / level;
 
level = 4;
quantized_image6 = double(C1) / 255;
quantized_image6 = uint8(quantized_image6 * level);
quantized_image6 = double(quantized_image6) / level;

%concatenate images
quantcolor_img = [quantized_image5, quantized_image6];  
              
%display
figure,imshow(quantcolor_img); 
title('Quantized Color Image')

% ---------- Part F - Quantize original 3-band into color image -----------

% Quantize original three color band image with a logarthmic function
% I' = Cln(I+1)
 
% Convert the image to a double and divide by 255 so the range of values 
% are [0,1]
% Multiply by k-level you want to quantize by 
% Convert the image to uint8
% To view convert to double and divide by k-level
 
level = 85*log(1+1); % I=1
log_quant1 = double(C1) / 255;
log_quant1 = uint8(log_quant1 * level);
log_quant1 = double(log_quant1) / level;
 
level = 85*log(4+1); % I=4
log_quant2 = double(C1) / 255;
log_quant2 = uint8(log_quant2 * level);
log_quant2 = double(log_quant2) / level;
 
level = 85*log(8+1); % I=8
log_quant3 = double(C1) / 255;
log_quant3 = uint8(log_quant3 * level);
log_quant3 = double(log_quant3) / level;
 
level = 85*log(256+1); %I=256
log_quant4 = double(C1) / 255;
log_quant4 = uint8(log_quant4 * level);
log_quant4 = double(log_quant4) / level;
 
%concatenate images 
log_image = [log_quant1 log_quant2; 
             log_quant3 log_quant4]; 
         
%display images
figure,imshow(log_image); 
title('Logarithmically Quantized Image')


