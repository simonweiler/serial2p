cd ('C:\Users\simonw\S1-V1 interaction Dropbox\Simon Weiler\Callosal_L6\NWG2023\Figures\anatomy\axonal_dist');
i1 = imread('axonal_distibution_max_zproject_SW_01531_100.tif'); 
i2 = imread('axonal_distibution_max_zproject_SW_01532_100micron.tif'); 
%%
%show image and crop
i1crop =[];
i1_mod= imrotate(i1,-30);
i1crop = imcrop(i1_mod);
i1crop=imadjust(i1crop);
imshow(i1crop);
imwrite (i1crop, 'i1crop.png');
%%
i2crop =[];
i2_mod= imrotate(i2,-30);
i2crop = imcrop(i2_mod);
i2crop=imadjust(i2crop);
imshow(i2crop);
imwrite (i2crop, 'i2crop.png');
%%
%read image 
cd ('C:\Users\simonw\S1-V1 interaction Dropbox\Simon Weiler\Callosal_L6\NWG2023\Figures\anatomy\axonal_dist');
i1crop = imread('i1crop.png');
i2crop = imread('i2crop.png');
%%
%choose points from L1 to WM
x_off=[0 20 40];
imshow(i1crop);
[x,y] = ginput(2);
%choose points from L1 to WM
for i=1:3
i1_profile = improfile(i1crop,[x(1)+x_off(i) ;x(1)+x_off(i)],y);
i1b_profile = improfile(i1crop,[x(1)+x_off(i) ;x(1)+x_off(i)],...
    [y(2)-200;y(2)-250]);
i1_profile_s(:,i)=i1_profile-nanmean(i1b_profile);
end

%%
%choose points from L1 to WM
x_off=[0 20 40];
imshow(i2crop);
[x,y] = ginput(2);
%choose points from L1 to WM
for i=1:3
i2_profile = improfile(i2crop,[x(1)+x_off(i) ;x(1)+x_off(i)],y);
i2b_profile = improfile(i2crop,[x(1)+x_off(i) ;x(1)+x_off(i)],...
    [y(2)-200;y(2)-250]);
i2_profile_s(:,i)=i2_profile-nanmean(i2b_profile);
end

%%
cd ('C:\Users\simonw\S1-V1 interaction Dropbox\Simon Weiler\Callosal_L6\NWG2023\Figures\anatomy\axonal_dist');
load('i1profile.mat');load('i2profile.mat');
i1pr=[];i2pr=[];
i1pr=nanmean(i1_profile_s,2);
i2pr=nanmean(i2_profile_s,2);
i1pr=i1pr./max(i1pr);
i2pr=i2pr./max(i2pr);
i1pr=smoothdata(i1pr);
i2pr=smoothdata(i2pr);
i1pr_s=i1pr(1:5:end);
i2pr_s=i2pr(1:5:end);
xnorm=[1:length(i1pr_s)]./max([1:length(i1pr_s)]);
% figure;plot(i1pr(1:1:end),1:1:length(i1pr));hold on;
% plot(i2pr(1:1:end),1:1:length(i2pr));
fig1= figure;set(fig1, 'Name', 'Intensity profiles');set(fig1, 'Position', [400, 500, 350, 400]);set(gcf,'color','w');
shadedErrorBar(xnorm,nanmean([i1pr_s i2pr_s(1:length(i1pr_s))],2)...
    ,nanstd([i1pr_s i2pr_s(1:length(i1pr_s))],[],2)/sqrt(2),'lineProps','k');
set(gca,'XDir','reverse');
view([90 -90]);ylim([-0.3 1]);xlim([0 1]);ylabel('Normalized Pixel intensity');
xlabel('Normalized depth from pia');set(gca,'FontSize',14);
saveas(gcf,'Normalized_Pixel_intensity_axons.emf')
%hold on;line([0.125 0.125],[-0.3 1],'linewidth',1,'Color',[0.5 0.5 0.5],'LineStyle','--');
%%
cd ('C:\Users\simonw\S1-V1 interaction Dropbox\Simon Weiler\Callosal_L6\NWG2023\Figures\anatomy\axonal_dist');
i1crop = imread('i1crop.png');
imshow(i1crop);
% hold on;line([0 587],[75 75],'linewidth',1,'Color','w','LineStyle','--');
% hold on;line([0 587],[175 175],'linewidth',1,'Color','w','LineStyle','--');
% hold on;line([0 587],[225 225],'linewidth',1,'Color','w','LineStyle','--');
% hold on;line([0 587],[340 340],'linewidth',1,'Color','w','LineStyle','--');
i1crop_figure=imcrop(i1crop);
imshow(i1crop_figure);
%%

imwrite (i1crop_figure, 'i1crop_figure.png');