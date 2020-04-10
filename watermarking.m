clc;
clear;

%DWT embedding
I = imread("peppers.tif"); %get the input image "peppers.tif"
W = load("watermark.mat"); %load the watermark
w1 = W.w1; %watermark 1
w2 = W.w2; %watermark 2
w3 = W.w3; %watermark 3
a = 0.4; %watermark strength
c = 0.4; %shape operator value

%display the original image
figure(1), imshow(I), title("original image");

%decomposition
%first decomposition of the original image
[LL1, HL1, LH1, HH1] = dwt2(double(I), 'haar', 'mode', 'per'); 
% second decomposition
[LL2, HL2, LH2, HH2] = dwt2(double(LL1), 'haar', 'mode', 'per'); 

%embed the watermarks into the appropriate sub-bands of the dwt
%transformed image using the embed function 
HL2_mark = embed(HL2, w1, a);
HH2_mark = embed(HH2, w2, a);
LH2_mark = embed(LH2, w3, a);

%rebuild image using inverse DWT
R_LL1 = idwt2(LL2, HL2_mark, LH2_mark, HH2_mark, 'haar', 'mode', 'per');
dwt_watermark = idwt2(R_LL1, HL1, LH1, HH1, 'haar', 'mode', 'per');
figure(2), imshow(uint8(dwt_watermark)), title("DWT domain watermark");
%write the image to a file
imwrite(uint8(round(dwt_watermark)), 'dwt_watermark.tif');

%DWT detection
%first decomposition of the watermarked image
[LL1_D, HL1_D, LH1_D, HH1_D] = dwt2(double(dwt_watermark), 'haar', 'mode', 'per'); 
% second decomposition 
[LL2_D, HL2_D, LH2_D, HH2_D] = dwt2(double(LL1_D), 'haar', 'mode', 'per');

%calculate phi value for each suband and respective watermark, using the
%phi function
dwtphi1 = phi(HL2_D, w1, a, c);
dwtphi2 = phi(HH2_D, w2, a, c);
dwtphi3 = phi(LH2_D, w3, a, c);

%calculate mu value for each watermark, using the mu function
dwtmu1 = mu(w1, a, c);
dwtmu2 = mu(w2, a, c);
dwtmu3 = mu(w3, a, c);

%calculate lambda value for each watermark, using lamba function
dwtlambda1 = lambda(w1, a, c);
dwtlambda2 = lambda(w2, a, c);
dwtlambda3 = lambda(w3, a, c);

%run the detect function on the DWT watermarked image to detect w1 
disp("DWT w1 detection");
disp(detect(dwtphi1, dwtlambda1, dwtmu1));

%run the detect function on the DWT watermarked image to detect w2
disp("DWT w2 detection");
disp(detect(dwtphi2, dwtlambda2, dwtmu2));

%run the detect function on the DWT watermarked image to detect w3
disp("DWT w3 detection");
disp(detect(dwtphi3, dwtlambda3, dwtmu3));

%--------------------------------------------------------------------------

%DCT embedding
I_dct = dct2(I);

%split image into 4 quadrants
Q1 = I_dct(1:256, 1:256); %top left quadrant
Q2 = I_dct(1:256, 257:512); %top right quadrant
Q3 = I_dct(257:512, 1:256); %bottom left quadrant
Q4 = I_dct(257:512, 257:512); %bottom right quadrant

%split Q1 into 4 sub-quadrants
SQ1 = Q1(1:128, 1:128); %top left quadrant of Q1
SQ2 = Q1(1:128, 129:256); %top right quadrant of Q1
SQ3 = Q1(129:256, 1:128); %bottom left quadrant of Q1
SQ4 = Q1(129:256, 129:256); %bottom right quadrant of Q1

%embed the watermarks into the appropriate sub quadrants of the dct
%transformed image using the embed function 
SQ2_mark = embed(SQ2, w1, a);
SQ4_mark = embed(SQ4, w2, a);
SQ3_mark = embed(SQ3, w3, a);

%reconstruct Q1
Q1_mark = [SQ1, SQ2_mark ; SQ3_mark, SQ4_mark];

%reconstruct the DCT image
I_dct_mark = [Q1_mark, Q2 ; Q3, Q4];

%inverse the DCT image to get the watermarked image
dct_watermark = idct2(I_dct_mark);

%display the reconstructed image that is watermarked
figure(3), imshow(uint8(dct_watermark)), title("DTC domain watermark");
%write the image to a file
imwrite(uint8(round(dct_watermark)), 'dct_watermark.tif');

%DCT detection

%DCT transform the watermarked image
dct_I = dct2(dct_watermark);

%select the upper left quadrant, where the watermark is 
Q1_D = dct_I(1:256, 1:256); %top left quadrant

%select the sub-quadrants where the watermark is held
SQ2_w1 = Q1_D(1:128, 129:256); %top right quadrant of Q1_D
SQ4_w2 = Q1_D(129:256, 129:256); %bottom right quadrant of Q1_D
SQ3_w3 = Q1_D(129:256, 1:128); %bottom left quadrant of Q1_D

%determine the phi values for DCT watermarked image
dctphi1 = phi(SQ2_w1, w1, a, c);
dctphi2 = phi(SQ4_w2, w2, a, c);
dctphi3 = phi(SQ3_w3, w3, a, c);

%determine the mu values for DCT
dctmu1 = mu(w1, a, c);
dctmu2 = mu(w2, a, c);
dctmu3 = mu(w3, a, c);

%determine the lambda values for DCT
dctlambda1 = lambda(w1, a, c);
dctlambda2 = lambda(w2, a, c);
dctlambda3 = lambda(w3, a, c);

%run the detect function on the DCT watermarked image to detect w1
disp("DWT w1 detection");
disp(detect(dctphi1, dctlambda1, dctmu1));

%run the detect function on the DCT watermarked image to detect w2
disp("DWT w2 detection");
disp(detect(dctphi2, dctlambda2, dctmu2));

%run the detect function on the DCT watermarked image to detect w3
disp("DWT w3 detection");
disp(detect(dctphi3, dctlambda3, dctmu3));

%--------------------------------------------------------------------------

%Section C

%circleshift of 1,1 of the DWT and DCT watermarked images
dwt_watermark_cs = circshift(dwt_watermark, [1,1]);
dct_watermark_cs = circshift(dct_watermark, [1,1]);

%3x3 average filter
filter = fspecial('average', [3,3]);

%average the DWT and DCT watermarked images with filter
dwt_watermark_avg = imfilter(dwt_watermark, filter);
dct_watermark_avg = imfilter(dct_watermark, filter);

%decomposition
%first decomposition of the original image
[LL1_cs, HL1_cs, LH1_cs, HH1_cs] = dwt2(double(dwt_watermark_cs), 'haar', 'mode', 'per'); 
% second decomposition
[LL2_cs, HL2_cs, LH2_cs, HH2_cs] = dwt2(double(LL1_cs), 'haar', 'mode', 'per');

%determine the phi values for the DWT watermarked image that has been circle
%shifted
dwtphi1_cs = phi(HL2_cs, w1, a, c);
dwtphi2_cs = phi(HH2_cs, w2, a, c);
dwtphi3_cs = phi(LH2_cs, w3, a, c);

%determine the mu values for DWT watermarked image
dwtmu1_cs = mu(w1, a, c);
dwtmu2_cs = mu(w2, a, c);
dwtmu3_cs = mu(w3, a, c);

%determine the lambda values for DWT watermarked image
dwtlambda1_cs = lambda(w1, a, c);
dwtlambda2_cs = lambda(w2, a, c);
dwtlambda3_cs = lambda(w3, a, c);

%run the detect function on the DCT watermarked image that has been
%circleshifted to detect w1
disp("DWT w1 detection, cricleshift");
disp(detect(dwtphi1_cs, dwtlambda1_cs, dwtmu1_cs));

%run the detect function on the DCT watermarked image that has been
%circleshifted to detect w2
disp("DWT w2 detection, cricleshift");
disp(detect(dwtphi2_cs, dwtlambda2_cs, dwtmu2_cs));

%run the detect function on the DCT watermarked image that has been
%circleshifted to detect w3
disp("DWT w3 detection, cricleshift");
disp(detect(dwtphi3_cs, dwtlambda3_cs, dwtmu3_cs));

%--------------------------------------------------------------------------

%decomposition
%first decomposition of the original image
[LL1_avg, HL1_avg, LH1_avg, HH1_avg] = dwt2(double(dwt_watermark_avg), 'haar', 'mode', 'per'); 
% second decomposition
[LL2_avg, HL2_avg, LH2_avg, HH2_avg] = dwt2(double(LL1_avg), 'haar', 'mode', 'per');

%determine the phi values for the DWT watermarked image that has been 3x3
%averaged
dwtphi1_avg = phi(HL2_avg, w1, a, c);
dwtphi2_avg = phi(HH2_avg, w2, a, c);
dwtphi3_avg = phi(LH2_avg, w3, a, c);

%determine the mu values for DWT watermarked image
dwtmu1_avg = mu(w1, a, c);
dwtmu2_avg = mu(w2, a, c);
dwtmu3_avg = mu(w3, a, c);

%determine the lambda values for DWT watermarked image
dwtlambda1_avg = lambda(w1, a, c);
dwtlambda2_avg = lambda(w2, a, c);
dwtlambda3_avg = lambda(w3, a, c);

%run the detect function on the DWT watermarked image that has been 3x3
%averaged to detect w1
disp("DWT w1 detection, 3x3 average");
disp(detect(dwtphi1_avg, dwtlambda1_avg, dwtmu1_avg));

%run the detect function on the DWT watermarked image that has been 3x3
%averaged to detect w2
disp("DWT w2 detection, 3x3 average");
disp(detect(dwtphi2_avg, dwtlambda2_avg, dwtmu2_avg));

%run the detect function on the DWT watermarked image that has been 3x3
%averaged to detect w3
disp("DWT w3 detection, 3x3 average");
disp(detect(dwtphi3_avg, dwtlambda3_avg, dwtmu3_avg));

%--------------------------------------------------------------------------

dct_I_cs = dct2(dct_watermark_cs);

%select the upper left quadrant, where the watermark is 

Q1_D_cs = dct_I_cs(1:256, 1:256); %top left quadrant

%select the sub-quadrants where the watermark is held

SQ2_w1_cs = Q1_D_cs(1:128, 129:256); %top right quadrant of Q1_D
SQ4_w2_cs = Q1_D_cs(129:256, 129:256); %bottom right quadrant of Q1_D
SQ3_w3_cs = Q1_D_cs(129:256, 1:128); %bottom left quadrant of Q1_D

%determine the phi values for the DCT watermarked image that has been
%circle shifted
dctphi1_cs = phi(SQ2_w1_cs, w1, a, c);
dctphi2_cs = phi(SQ4_w2_cs, w2, a, c);
dctphi3_cs = phi(SQ3_w3_cs, w3, a, c);

%determine the mu values for the DCT watermarked image
dctmu1_cs = mu(w1, a, c);
dctmu2_cs = mu(w2, a, c);
dctmu3_cs = mu(w3, a, c);

%determine the lambda values for the DCT watermarked image
dctlambda1_cs = lambda(w1, a, c);
dctlambda2_cs = lambda(w2, a, c);
dctlambda3_cs = lambda(w3, a, c);

%run the detect function on the DCT watermarked image that has been circle
%shifted to detect w1
disp("DWT w1 detection, cricleshift");
disp(detect(dctphi1_cs, dctlambda1_cs, dctmu1_cs));

%run the detect function on the DCT watermarked image that has been circle
%shifted to detect w2
disp("DWT w2 detection, cricleshift");
disp(detect(dctphi2_cs, dctlambda2_cs, dctmu2_cs));

%run the detect function on the DCT watermarked image that has been circle
%shifted to detect w3
disp("DWT w3 detection, cricleshift");
disp(detect(dctphi3_cs, dctlambda3_cs, dctmu3_cs));

%--------------------------------------------------------------------------

dct_I_avg = dct2(dct_watermark_avg);

%select the upper left quadrant, where the watermark is 
Q1_D_avg = dct_I_avg(1:256, 1:256); %top left quadrant

%select the sub-quadrants where the watermark is held
SQ2_w1_avg = Q1_D_avg(1:128, 129:256); %top right quadrant of Q1_D
SQ4_w2_avg = Q1_D_avg(129:256, 129:256); %bottom right quadrant of Q1_D
SQ3_w3_avg = Q1_D_avg(129:256, 1:128); %bottom left quadrant of Q1_D

%determine the phi values for the DCT watermarked image that has been 3x3
%averaged
dctphi1_avg = phi(SQ2_w1_avg, w1, a, c);
dctphi2_avg = phi(SQ4_w2_avg, w2, a, c);
dctphi3_avg = phi(SQ3_w3_avg, w3, a, c);

%determine the mu values for the DCT watermarked image
dctmu1_avg = mu(w1, a, c);
dctmu2_avg = mu(w2, a, c);
dctmu3_avg = mu(w3, a, c);

%determine the lambda values for the DCT watermarked image
dctlambda1_avg = lambda(w1, a, c);
dctlambda2_avg = lambda(w2, a, c);
dctlambda3_avg = lambda(w3, a, c);

%run the detect function on the DCT watermarked image that has been 3x3
%averaged to detect w1
disp("DWT w1 detection, 3x3 average");
disp(detect(dctphi1_avg, dctlambda1_avg, dctmu1_avg));

%run the detect function on the DCT watermarked image that has been 3x3
%averaged to detect w2
disp("DWT w2 detection, 3x3 average");
disp(detect(dctphi2_avg, dctlambda2_avg, dctmu2_avg));

%run the detect function on the DCT watermarked image that has been 3x3
%averaged to detect w3
disp("DWT w3 detection, 3x3 average");
disp(detect(dctphi3_avg, dctlambda3_avg, dctmu3_avg));

%--------------------------------------------------------------------------

%phi function
function [output] = phi(x, y, z, d)
    % x = sub-band, 
    % y = watermark, 
    % z = watermark strength value, 
    % d = shape operator value
    
    %initialize a temporary matrix of zeros
    temporary = zeros(128,128);
    
    %loop through each index and update the pixel at ith and jth location
    for i = 1:size(y)
        for j = 1:size(y)
            
            %store the values in a temporary matrix
            temporary(i, j) = sum((abs(37.5 * x(i, j)) ./ std(x(i, j)))^d * (1 - (1 ./ (1 + (z * y(i, j)))^d)));
        end
    end
    
    %output the result of the loop into the variable that called the
    %function
    output = temporary;
end

%mu function
function [output] = mu(x, y, z)
   % x = watermark,
   % y = watermark strength value,
   % z = shape operator value
   
   %initialize a temporary matrix of zeros
   temporary = zeros(128,128);
   
   %loop through each index and update the ith and jth pixel
   for i = 1:size(x)
       for j = 1:size(x)
           
           %store the values in a temporary matrix
           temporary(i, j) = sum((1 ./ z) * (1 - (1 ./ (1 + (y * x(i, j)))^z)));
       end
   end
   
   %output the result of the loop into the variable that called the
   %function
   output = temporary;
end

%lambda function
function [output] = lambda(x, y, z)
    % x = watermark
    % y = watermark strength value,
    % z = shape operator value
    
    %initialize a temporary matrix of zeros   
    temporary = zeros(128, 128);
    
    %loop through each index and update the pixel at i,j
    for i = 1:size(x)
        for j = 1:size(x)
            
            %store the values in a temporary matrix
            temporary(i, j) = sum((1 ./ z) * (1 - (1 ./ (1 + (y * x(i, j))^z)))^2);
        end
    end
    
    %output the result of the loop into the variable that called the
    %function
    output = temporary;
end

%detection function
function [output] = detect(x, y, z)
    
    b = sqrt(y);
       
    if(x > (4.75 * b + z))
        
        output = 1;
    else
        
        output = 0;
    end
end

%embedding function
function [output] = embed(x, y, z)
    % x = sub-band,
    % y = watermark,
    % z = watermark strength value

    %loop through the ith and jth indexes and update each index with the
    %multiplicative embedding rule
    for i = 1:size(x)
        for j = 1:size(x)
            x(i, j) = x(i, j) + (z * y(i, j) * x(i, j));
        end
    end
            
    %output the result of the loop into the variable that called the
    %function
    output = x;
end
