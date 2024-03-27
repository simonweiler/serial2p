function [data3_m data4_m] = zero_area_layers(data3,data4,thr,an)
data=squeeze(nansum(data3(:,:,:)));
data_m=data;

data3_m=data3;
data4_m=data4;

for i=1:size(data,2)
        data_m(find(data(:,i)<thr),i)=0;
end

for i=1:size(data,2)
        data3_m(:,find(data(:,i)<thr),i)=zeros(6,length(find(data(:,i)<thr)));  
        data4_m(:,find(data(:,i)<thr),i)=zeros(6,length(find(data(:,i)<thr)));  

end



ind_a=find(sum(data_m>0,2)<an);
for i=1:size(data_m,2)
data3_m(:,ind_a,i)=zeros(6,length(ind_a));
data4_m(:,ind_a,i)=zeros(6,length(ind_a));
end

idxnan=find(sum(isnan(data3(3,:,:)),3)==size(data,2));
for i=1:size(data,2)
data3_m(3,idxnan,i)=NaN;
data4_m(3,idxnan,i)=NaN;
end

end