function [p_bin p_bin_c] = anatomy_testarea(hemi_idx,ipsi_val,contra_val)

for i=1:length(hemi_idx);
    h=[];p=[];pi=[];
   
   [h p] = ttest(hemi_idx(i,:),0);
    %[h p] = ttest(contra_val(i,:),ipsi_val(i,:),"Tail","right");
    try
    %[p h] = signrank(hemi_idx(i,:),0,"Tail","right");
     [jj(i) pp] = adtest(hemi_idx(i,:));
    catch
        %p=0.9;
         jj(i)=0;
    end
    p_bin_c(i)=p<0.05/(length(hemi_idx)-1);

    [h pi] = ttest(hemi_idx(i,:),0);
    %[h pi] = ttest(contra_val(i,:),ipsi_val(i,:),"Tail","left");
   try
     %[pi h] = signrank(hemi_idx(i,:),0,"Tail","left");
   catch
        %pi=0.9;
   end
    p_bin(i)=pi<0.05/(length(hemi_idx)-1);
end
end