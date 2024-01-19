function [data_m data2_m data3_m] = zero_area(data,data2,data3,thr)
data_m=data;
data2_m=data2;
data3_m=data3;

    for i=1:size(data,2)
        data_m(find(data(:,i)<thr),i)=0;
        data2_m(find(data(:,i)<thr),i)=0;
    end

for i=1:size(data,2)
        data3_m(:,find(data(:,i)<thr),i)=zeros(6,length(find(data(:,i)<thr)));
        
    end
end