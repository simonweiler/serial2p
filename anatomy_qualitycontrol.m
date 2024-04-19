function anatomy_qualitycontrol(ipsi_all,contra_all,thresholds,thresh_ani,prct)

fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 0, 1500, 800]);set(gcf,'color','w');
tiledlayout(2,3);

%color1='m';
% nexttile
% b1=bar([sum(ipsi_all{1}); sum(contra_all{1})]',0.7);b1(1).FaceColor=[0.8 0.8 0.8];b1(2).FaceColor=[0.5 0.5 0.5];set(b1,'ShowBaseLine','off')
% ylabel('Total cell number');xlabel('Animal');set(gca,'FontSize',8);set(gca,'TickDir','out');box off;xtickangle(45);
% offsetAxes;title('VISp ipsi vs contra');
panel_tit={'VISp','SSp-bfd','MOp'};
for i=1:3
nexttile
pp=[];pp1=sort(nanmean(ipsi_all{i},2),'descend');
hold on;line([0 45],[prctile(pp1,prct) prctile(pp1,prct)],'LineStyle','--','Color',[0.5 0.5 0.5]);
%text(46,prctile(pp,prct),num2str(round(prctile(pp,prct)),2))
ylim([0 max(pp1)]);hold on,p1=plot(pp1,'-o','Color',[0.5 0.5 0.5]);p1.MarkerSize=3;

hold on
 pp2=[];pp2=sort(nanmean(contra_all{i},2),'descend');
 hold on;line([0 45],[prctile(pp2,prct) prctile(pp2,prct)],'LineStyle','--','Color',[0.5 0.5 0.5]);
 %text(53,prctile(pp,prct),num2str(round(prctile(pp,prct)),2))
 ylim([0 max(pp2)]);hold on,p1=plot(pp2,'-o','Color',[0.8 0.8 0.8]);p1.MarkerSize=3;

hold on
 pp3=[];pp3=sort(nanmean([ipsi_all{i} contra_all{i}],2),'descend');
 hold on;line([0 45],[prctile(pp3,prct) prctile(pp3,prct)],'LineStyle','--','Color',[0.5 0.5 0.5]);
%text(58,prctile(pp,prct),num2str(round(prctile(pp,prct)),2))

 ylim([0 max(pp3)]);hold on,p1=plot(pp3,'-o','Color',[0.1 0.1 0.1]);p1.MarkerSize=3;
 ylabel('Average cell number');xlabel('Areas');set(gca,'FontSize',8);set(gca,'TickDir','out');box off;xtickangle(45);
 offsetAxes;
 title(panel_tit{i});
 prctile_across(:,i)=[prctile(pp1,prct) prctile(pp2,prct) prctile(pp3,prct)];
end

nexttile
 plot(prctile_across,'-o');legend({'VISp','SSp','MOp'})
 ylabel('Cell numbers');;set(gca,'FontSize',8);set(gca,'TickDir','out');box off;xtickangle(45);
 offsetAxes;title('All animals ipsi vs contra');xticks([1 2 3]);xticklabels({'ipsi','contra','together'})


nexttile
 b1=bar([sum(ipsi_all{1}) sum(ipsi_all{2}) sum(ipsi_all{3}); sum(contra_all{1}) sum(contra_all{2}) sum(contra_all{3})]',1);b1(1).FaceColor=[0.8 0.8 0.8];b1(2).FaceColor=[0.5 0.5 0.5];set(b1,'ShowBaseLine','off')
 ylabel('Total cell number');xlabel('Animal');set(gca,'FontSize',8);set(gca,'TickDir','out');box off;xtickangle(45);
 offsetAxes;title('All animals ipsi vs contra');

nexttile
 %add all average counts across areas and animals 
all_total=[];all_total=[ipsi_all{1}(:) ipsi_all{2}(:) ipsi_all{3}(:) contra_all{1}(:) contra_all{2}(:) contra_all{3}(:)];
[temp_s kk]=sort(all_total(:),'descend');
%Plot 
pz=plot(temp_s,'-o','Color','k');pz.MarkerFaceColor=[1 1 1];pz.MarkerEdgeColor=[0 0 0];pz.MarkerSize=2;idx_prct=[];idx_prct=find(temp_s>prctile(temp_s,prct));
hold on;plot(idx_prct(end),prctile(temp_s,prct)+3000,'v');hold on;text(idx_prct(end)+5,prctile(temp_s,prct)+4000,['cells= ' num2str(temp_s(idx_prct(end)))]);
text(idx_prct(end)+5,prctile(temp_s,prct)+8000,[' prct=' num2str(prct)])
ylabel('Total cell number');xlabel('Sorted areas across all');set(gca,'FontSize',11);set(gca,'TickDir','out');box off;xtickangle(45);ylim([-1500 max(temp_s)]);xlim([-5 length(temp_s)+1]);
offsetAxes;


%% 

fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 0, 1500, 800]);set(gcf,'color','w');
tiledlayout(2,3);
for i=1:3
nexttile
keeper = ([ipsi_all{1} ipsi_all{2} ipsi_all{3}]>thresholds(i));

pcolor(keeper)
title(['t= ' num2str(thresholds(i))])
end
%keeper_all(:,i)=sum(keeper,2);
for i=1:3
nexttile
keeper = ([contra_all{1} contra_all{2} contra_all{3}]>thresholds(i));

pcolor(keeper)
title(['t= ' num2str(thresholds(i))])
end

% colorbar
% 
% nexttile
% pcolor(keeper_all>=thresh_ani);colorbar
end