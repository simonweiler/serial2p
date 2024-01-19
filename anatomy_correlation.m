%calculates correlation values and average/SEM of data points to be correlted
function [value1_all value2_all value1 value2 xerr yerr r_a p_a r_avg p_avg] = anatomy_correlation(rm1, rm2, idx)
%input to function: 
%rm1=
%rm2=
%idx=indexes to be correlated 
%output
%r_a/p_a= individual animal/modules based correlation/pvalues
%avg p_avg= correlation on average data points
%dependencies: anatomy_indexcalc.m

i_metric=[];c_metric=[];
%calculate indexes across areas
for i=1:size(rm1,2)
temp_l=[];temp_l2=[];
temp_l=squeeze(rm1(:,i,:));
temp_l2=squeeze(rm2(:,i,:));
[i_metric(:,i,:)] = anatomy_indexcalc(temp_l);
[c_metric(:,i,:)] = anatomy_indexcalc(temp_l2);
end 
temp_metric1=i_metric(:,:,idx);temp_metric2=c_metric(:,:,idx);

%Calculate correaltions across Ns
r_a=[];p_a=[];
for i=1:size(temp_metric1,1)
r=[];p=[];
[r p] = corr(temp_metric1(i,:)',(temp_metric2(i,:)-temp_metric1(i,:))','Type','Spearman','Rows','complete');
r_a(i)=r(1);
p_a(i)=p(1);
end
%average values 

value1=nanmean(temp_metric1);
value2=nanmean(temp_metric2)-nanmean(temp_metric1);
value1_all=temp_metric1;
value2_all=temp_metric2-temp_metric1;
%xSEM
xerr=[];
xerr=nanstd(temp_metric1);
xerr=xerr/sqrt(length(xerr(~isnan(xerr))));

%ySEM
yerr=[];
yerr=nanstd(temp_metric2-temp_metric1);
yerr=yerr/sqrt(length(yerr(~isnan(yerr))));
%average correlation 
[r_avg p_avg] = corr(value1',value2','Type','Spearman','Rows','complete');
end