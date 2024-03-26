function [data_ilni data_ilnc order_i] = iln_change(rm1,rm2,idx)

for i=1:length(rm1)
temp_l=[];temp_l2=[];
temp_l=squeeze(rm1(:,i,:));
temp_l2=squeeze(rm2(:,i,:));
[i_metric(:,i,:)] = anatomy_indexcalc(temp_l);
[c_metric(:,i,:)] = anatomy_indexcalc(temp_l2);
end 

temp_metric1=i_metric(:,:,idx);temp_metric2=c_metric(:,:,idx);
temp_rm1=[];
temp_rm2=[];
temp_rm1=nanmean(temp_metric1);
temp_rm2=nanmean(temp_metric2);
[data_ilni order_i]=sort(temp_rm1,'ascend');
data_ilnc=temp_rm2;
data_ilni=temp_rm1;
end