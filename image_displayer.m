cd ('C:\Users\simonw\S1-V1 interaction Dropbox\Simon_Manuel\Manuscripts\L6_dominance\Injections\retro\V1\CT_overlap\211026_SW5070_NTSR1\V1\slice1');
c1 =  tiffreadVolume('ch1.tif',...
    'PixelRegion', {[1 2 inf], [1 2 inf], [10 50]});
c2 =  tiffreadVolume('ch2.tif',...
    'PixelRegion', {[1 2 inf], [1 2 inf], [10 50]});
 c1_a= max(c1,[],3);
 c2_a= max(c2,[],3);
%% 

 ch_merge = imfuse(c1_a,c2_a,'falsecolor','Scaling','joint','ColorChannels',[1 2 0]);