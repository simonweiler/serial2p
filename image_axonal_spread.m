cd ('C:\Users\simonw\S1-V1 interaction Dropbox\Simon Weiler\Callosal_L6\NWG2023\Figures\anatomy\axonal_dist');
i1 = imread('axonal_distibution_max_zproject_SW_01531_100.tif'); 
i2 = imread('axonal_distibution_max_zproject_SW_01532_100micron.tif'); 
%%
%show image and crop

i1_mod= imrotate(i1,-30);
i1crop = imcrop(i1_mod);
i1crop=imadjust(i1crop);
imshow(i1crop);

%%

imshow(i1crop);
i1_profile = improfile(i1crop)