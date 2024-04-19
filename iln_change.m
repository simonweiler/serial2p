function [data_ilni data_ilnc rm1_error rm2_error order_i order_c temp_metric1 temp_metric2] = iln_change(rm1,rm2,idx)

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
rm1_error=nanstd(temp_metric1)/sqrt(size(temp_metric1,1));
rm2_error=nanstd(temp_metric2)/sqrt(size(temp_metric2,1));
[data_xxi order_i]=sort(temp_rm1,'ascend');
[data_xxc order_c]=sort(temp_rm2,'ascend');
data_ilnc=temp_rm2;
data_ilni=temp_rm1;
end