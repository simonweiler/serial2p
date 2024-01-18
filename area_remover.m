function [remove_idx remove_idx_com] = area_remover(frac_nr, tot_nr,prct,cortex_abb,title_panel)

fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 500, 2200, 800]);set(gcf,'color','w');tiledlayout("horizontal")
t.TileSpacing = 'compact';t.Padding = 'compact';
%title_panel={'VISp ipsi','SSp ipsi','MOp ipsi','VISp contra','SSp contra','MOp contra'}
%V1
keep_idx=[];
for i=1:size(frac_nr,2)
temp_s=[];
[temp_s kk]=sort(frac_nr(:,i),'descend');
temp_st=sort(tot_nr(:,i),'descend');
nexttile
plot(temp_s,'o','Color','k');box off;
hold on;title(title_panel{i});
idx_prct=[];
idx_prct=find(temp_s>prctile(temp_s,prct))
hold on;plot(idx_prct(end),prctile(temp_s,prct)+0.01,'v');
hold on;text(idx_prct(end),prctile(temp_s,prct)+0.02,[num2str(temp_st(idx_prct(end))) ' prct=' num2str(prct)]);
remove_idx(:,i)=find(frac_nr(:,i)<prctile(frac_nr(:,i),prct));
xticks([1:1:45]);
xticklabels(cortex_abb(kk));
xtickangle(45);

end
nexttile
temp_s=[];
[temp_s kk]=sort(frac_nr(:),'descend');
temp_st=sort(tot_nr(:),'descend');
plot(temp_s,'o','Color','k');box off;
idx_prct=[];
idx_prct=find(temp_s>prctile(temp_s,prct))
idx_prct=find(temp_s>prctile(temp_s,prct))
hold on;plot(idx_prct(end),prctile(temp_s,prct)+0.01,'v');

title('Combined')
all_thrsh=prctile(frac_nr(:),prct);
for i=1:size(frac_nr,2)
remove_idx_com{:,i}=find(frac_nr(:,i)<all_thrsh);
end
end