%This scipt reads out the saved mat data structure with all animals and
%injection inlcuded and plots figures/ performs analysis 
%% Dependendencies 
%uipickfiles: https://uk.mathworks.com/matlabcentral/fileexchange/10867-uipickfiles-uigetfile-on-steroids
%Matlab structure
%% Color scheme for the plots
v1_color=[173 7 227]./256;%magenta V1
s1_color=[5 190 120]./256;%green S1
m1_color=[87 87 244]./256;%blue M1
%% Folder where mat structure is (CHANGE ACCORDINGLY)
str   = 'C:\Users\simonw\S1-V1 interaction Dropbox\Simon_Manuel\Manuscripts\L6_dominance\Injections\output_structure';
% load structure (mat file) using uipickfiles
folder_list = uipickfiles('FilterSpec',str);load(char(folder_list));
%% Load csv with cortex names and assign 6 larger group areas: Frontal, Lateral, Somatomotor, Visual, Medial, Auditory
cortex_info = readtable([str '\names_abbreviations.csv']);
%read out nr, abbreviation and long cortex name of ABI
cortex_names=table2cell(cortex_info);
%change here index of categories
frontal_idx=[1:8];lateral_idx=[9:16 44:45];somamo_idx=[17:26];visual_idx=[27:36];medial_idx=[37:39];aud_idx=[40:43];
%injected brain areas VISp, SSp, MOp, (AUDp)
visp_idx=[31];ssp_idx=[18];mop_idx=[25];audp_idx=[41];
%% STARTS HERE

%% Find the cells in a given area in their given layer RETRO
% select all retro injection regardeless of their injection type
temp1=[];temp1=cell_selecter(data,'virustype',1);
retro=[];retro=find(temp1==1);
%call function anatomy_cellnr
[i_animal c_animal i_areas_animal c_areas_animal i_areas_animalf c_areas_animalf  i_areas_animalcd c_areas_animalcd] = anatomy_cellnr(data,retro,cortex_names);
%% Plot all cell numbers ipsi contra 
dat=[];dat=[squeeze(nansum(nansum(i_areas_animal)))./(squeeze(nansum(nansum(i_areas_animal)))+squeeze(nansum(nansum(c_areas_animal))))...
    squeeze(nansum(nansum(c_areas_animal)))./(squeeze(nansum(nansum(i_areas_animal)))+squeeze(nansum(nansum(c_areas_animal))))]
cl={[0.8 0.8 0.8],[0.3 0.3 0.3]};
paired_plot_box(dat,cl);
xticklabels({'ipsi','contra'});ylabel('Fraction of neurons');hold on;title([]);set(gca,'FontSize',11);set(gca,'TickDir','out'); 
%[p1]=signrank(dat(:,1),dat(:,2))
[u p1]=ttest(dat(:,1),dat(:,2))
hold on;text(1.25,max(dat(:,1)),['***'],'FontSize',18);
%% Plot barplot for all regions all areas RETRO (animals=18, areas=45)
%means
means_layers=[nanmean(i_animal(2,:)) nanmean(c_animal(2,:));nanmean(i_animal(3,:)) nanmean(c_animal(3,:));...
    nanmean(i_animal(4,:)) nanmean(c_animal(4,:));nanmean(i_animal(5,:)) nanmean(c_animal(5,:));nanmean(i_animal(6,:)) nanmean(c_animal(6,:))];
semd=sqrt(length(i_animal));
%error
err_layers=[[nanstd(i_animal(2,:))/semd nanstd(c_animal(2,:))/semd];[nanstd(i_animal(3,:))/semd nanstd(c_animal(3,:))/semd];...
    [nanstd(i_animal(4,:))/semd nanstd(c_animal(4,:))/semd];[nanstd(i_animal(5,:))/semd nanstd(c_animal(5,:))/semd] ;[nanstd(i_animal(6,:))/semd nanstd(c_animal(6,:))/semd]];
%plotting starts here 
fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 500, 500, 300]);set(gcf,'color','w');
%bars
b=bar(means_layers,1);
b(1).FaceColor=[0.8 0.8 0.8];b(2).FaceColor=[0.3 0.3 0.3];b(1).EdgeColor=[1 1 1];b(2).EdgeColor=[1 1 1];hold on;
%indivdual animals
for r=1:5
rng(1);r1 = r-0.2 + (0.2)*rand(length(retro),1);
sc1=scatter(r1,i_animal(r+1,:),5,'k','filled');hold on;sc1.MarkerEdgeColor=[0.5 0.5 0.5];sc1.MarkerEdgeAlpha=0.5;sc1.MarkerFaceColor=[0.2 0.2 0.2];
end
hold on;
for r=1:5
rng(1);r1 = r+0.1 + (0.2)*rand(length(retro),1);
sc1=scatter(r1,c_animal(r+1,:),5,'k','filled');hold on;sc1.MarkerEdgeColor=[0.5 0.5 0.5];sc1.MarkerEdgeAlpha=0.5;sc1.MarkerFaceColor=[0.2 0.2 0.2];
end
hold on;
%indivdual errorbars
for k = 1:size(means_layers,2)
    % get x positions per group
    xpos = b(k).XData + b(k).XOffset;
    % draw errorbar
    errorbar(xpos, means_layers(:,k), err_layers(:,k), 'LineStyle', 'none', ... 
        'Color', 'k', 'LineWidth', 2);
end
xticklabels({'L2/3','L4','L5','L6a','L6b'});ylabel('Fraction of neurons');set(gca,'FontSize',11);set(gca,'TickDir','out');box off;
%[p1]=signrank(i_animal(6,:),c_animal(6,:))
[u p1]=ttest(i_animal(6,:),c_animal(6,:))
hold on;text(0.75,0.6,['***'],'FontSize',18);hold on;text(3.75,0.6,['***'],'FontSize',18);legend({'ipsi','contra'},"Location","northeast");legend boxoff;
%hold on;text(0.75,0.6,['**'],'FontSize',18);hold on;text(3.75,0.6,['*'],'FontSize',18);legend({'ipsi','contra'},"Location","northeast");legend boxoff;title('MOp','Color',m1_color)
%% Extract and Plot delta contra-ipsi fraction per injeciton area M1 / S1 / V1
temp2=[];temp2=cell_selecter(data,'virustype',1,'area',1);
temp3=[];temp3=cell_selecter(data,'virustype',1,'area',2);
temp4=[];temp4=cell_selecter(data,'virustype',1,'area',3);
v1r=[];v1r=find(temp2==1);s1r=[];s1r=find(temp3==1);m1r=[];m1r=find(temp4==1);
%call function anatomy_cellnr
%V1, S1 and M1
[i_v1a c_v1a i_v1aa c_v1aa i_v1aaf c_v1aaf i_v1aacd c_v1aacd] = anatomy_cellnr(data,v1r,cortex_names);
[i_s1a c_s1a i_s1aa c_s1aa i_s1aaf c_s1aaf i_s1aacd c_s1aacd] = anatomy_cellnr(data,s1r,cortex_names);
[i_m1a c_m1a i_m1aa c_m1aa i_m1aaf c_m1aaf i_m1aacd c_m1aacd] = anatomy_cellnr(data,m1r,cortex_names);
%calculate delta
v1_d=c_v1a-i_v1a;
s1_d=c_s1a-i_s1a;
m1_d=c_m1a-i_m1a;
%Plot figure 
fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 500, 700, 300]);set(gcf,'color','w');tiledlayout("horizontal")
t.TileSpacing = 'compact';t.Padding = 'compact';
%V1
nexttile
b=bar(nanmean(v1_d(2:end,:),2));b.FaceColor=v1_color;box off;
hold on;errorbar([1 2 3 4 5], nanmean(v1_d(2:end,:),2), nanstd(v1_d(2:end,:),[],2)/(sqrt(length(v1_d))), 'LineStyle', 'none', ... 
        'Color', 'k', 'LineWidth', 1);
for r=1:5
rng(1);r1 = r-0.01 + (0.2)*rand(length(v1_d(2:end,:)),1);
sc1=scatter(r1,v1_d(r+1,:),5,'k','filled');hold on;sc1.MarkerEdgeColor=[1 1 1];sc1.MarkerEdgeAlpha=0.5;sc1.MarkerFaceColor=[0.2 0.2 0.2];
end
xticklabels({'L2/3','L4','L5','L6a','L6b'});ylabel('Fractional change (contra - ipsi)');
set(gca,'FontSize',11);set(gca,'TickDir','out');box off;ylim([-0.15 0.25]);
title('VISp','FontWeight','normal','Color',v1_color);
hold on;text(3.75,0.235,['*'],'FontSize',18);hold on;text(0.75,0.235,['*'],'FontSize',18);

%S1
nexttile
b=bar(nanmean(s1_d(2:end,:),2));b.FaceColor=s1_color;box off;
hold on;errorbar([1 2 3 4 5], nanmean(s1_d(2:end,:),2), nanstd(s1_d(2:end,:),[],2)/(sqrt(length(s1_d))), 'LineStyle', 'none', ... 
        'Color', 'k', 'LineWidth', 1);
for r=1:5
rng(1);r1 = r-0.01 + (0.2)*rand(length(s1_d(2:end,:)),1);
sc1=scatter(r1,s1_d(r+1,:),5,'k','filled');hold on;sc1.MarkerEdgeColor=[1 1 1];sc1.MarkerEdgeAlpha=0.5;sc1.MarkerFaceColor=[0.2 0.2 0.2];
end
xticklabels({'L2/3','L4','L5','L6a','L6b'});
set(gca,'FontSize',11);set(gca,'TickDir','out');box off;
ylim([-0.15 0.25]);set(gca,'box','off','ycolor','w');line([0 6],[0 0],'LineWidth', 0.5, 'Color','k');
title('SSpbf','FontWeight','normal','Color',s1_color);
hold on;text(3.75,0.235,['*'],'FontSize',18);hold on;text(0.75,0.235,['*'],'FontSize',18);

%M1
nexttile
b=bar(nanmean(m1_d(2:end,:),2));b.FaceColor=m1_color;box off;
hold on;errorbar([1 2 3 4 5], nanmean(m1_d(2:end,:),2), nanstd(m1_d(2:end,:),[],2)/(sqrt(length(m1_d))), 'LineStyle', 'none', ... 
        'Color', 'k', 'LineWidth', 1);
for r=1:5
rng(1);r1 = r-0.01 + (0.2)*rand(length(m1_d(2:end,:)),1);
sc1=scatter(r1,m1_d(r+1,:),5,'k','filled');hold on;sc1.MarkerEdgeColor=[1 1 1];sc1.MarkerEdgeAlpha=0.5;sc1.MarkerFaceColor=[0.2 0.2 0.2];
end
xticklabels({'L2/3','L4','L5','L6a','L6b'});
set(gca,'FontSize',11);set(gca,'TickDir','out');box off;
ylim([-0.15 0.25]);set(gca,'box','off','ycolor','w');line([0 6],[0 0],'LineWidth', 0.5, 'Color','k');
title('MOp','FontWeight','normal','Color',m1_color);
hold on;text(3.75,0.235,['*'],'FontSize',18);hold on;text(0.75,0.235,['*'],'FontSize',18);
[u p1]=ttest(c_v1a(5,:),i_v1a(5,:))
[u p2]=ttest(c_s1a(5,:),i_s1a(5,:))
[u p3]=ttest(c_m1a(5,:),i_m1a(5,:))
%% calculate indexes: ILN, L6ab, L6a, hindex ALL INJECTIONS TOGETHER
%1= ILN
%2= L6ab dominance index
%3= L6a dominacne index
%4= h-index
%5= L6 pure index 
[i_index] = anatomy_indexcalc(i_animal);
[c_index] = anatomy_indexcalc(c_animal);
%% L6a fraction (pure index)
dat=[];dat=[i_index(:,5) c_index(:,5)];cl={[0.8 0.8 0.8],[0.3 0.3 0.3]};
paired_plot_box(dat,cl);
ylabel('L6a / (L23+L5+L6ab)');xticklabels({'ipsi','contra'});xtickangle(45);ylim([0.08 0.6]);
%pvalue of paired signrank test: 
%[p1]=signrank(dat(:,1),dat(:,2))
[u p1]=ttest(dat(:,1),dat(:,2))
hold on;text(1.25,0.6,['***'],'FontSize',18);
title([]);set(gca,'FontSize',11);set(gca,'TickDir','out');
 
%% calculate delta for specific area: V1, S1, M1
[iv1_index] = anatomy_indexcalc(i_v1a);
[cv1_index] = anatomy_indexcalc(c_v1a);
[is1_index] = anatomy_indexcalc(i_s1a);
[cs1_index] = anatomy_indexcalc(c_s1a);
[im1_index] = anatomy_indexcalc(i_m1a);
[cm1_index] = anatomy_indexcalc(c_m1a);
%% L6a fraction (pure index) average colour coded for V1 S1 M1
temp_color=[v1_color ;s1_color; m1_color];
module_names={'VISp','SSp-bfd','MOp'};
all_dati=[];all_dati=[iv1_index(:,5) is1_index(:,5) im1_index(:,5)];
all_datc=[];all_datc=[cv1_index(:,5) cs1_index(:,5) cm1_index(:,5)];
%plot
fig3= figure;set(fig3, 'Name', 'Paired comp');set(fig3, 'Position', [200, 300, 170, 250]);set(gcf,'color','w');
for j=1:3
dat=[];dat=[all_dati(:,j) all_datc(:,j)];
for i=1:length(data)
     pl=plot([1,2],[dat(:,1),dat(:,2)],'color',[0.5 0.5 0.5]);    
end
hold on;pS=plotSpread([dat(:,1),dat(:,2)],'categoryIdx',[ones(1,length(dat(:,1)))' ones(1,length(dat(:,2)))'*2],...
    'categoryMarkers',{'o','o'},'categoryColors',{temp_color(j,:), temp_color(j,:)});hold on;
hold on;er1=errorbar([0.75],nanmean(dat(:,1)),nanstd(dat(:,1),[],1)/sqrt(length(dat(:,1))),'ok','MarkerFaceColor',temp_color(j,:),'Markersize',8);
hold on;er2=errorbar([2.25],nanmean(dat(:,2)),nanstd(dat(:,2),[],1)/sqrt(length(dat(:,1))),'ok','MarkerFaceColor',temp_color(j,:),'Markersize',8);
text(2.2,1.1-0.07*j,module_names{j},'Color',temp_color(j,:))
end
xticklabels({'ipsi','contra'});ylabel('L6a / (L23+L5+L6ab)');hold on;title([]);xtickangle(45);set(gca,'FontSize',11);set(gca,'TickDir','out');box off;
%% Plot figure delta of indexes contra and ipsi 
inr=5;
fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 500, 170, 250]);set(gcf,'color','w');
b1=bar(1,[nanmean(cv1_index(:,inr)-iv1_index(:,inr))]);hold on;
b2=bar(2,[nanmean(cs1_index(:,inr)-is1_index(:,inr))]);hold on;
b3=bar(3,[nanmean(cm1_index(:,inr)-im1_index(:,inr))]);
b1.FaceColor=v1_color;b2.FaceColor=s1_color;b3.FaceColor=m1_color;b1.LineWidth=1;b2.LineWidth=1;b3.LineWidth=1;
xticks([1:3]);hold on;
hold on;errorbar([1 2 3], [nanmean(cv1_index(:,inr)-iv1_index(:,inr)) nanmean(cs1_index(:,inr)-is1_index(:,inr)) nanmean(cm1_index(:,inr)-im1_index(:,inr))]...
    , [(nanstd(cv1_index(:,inr)-iv1_index(:,inr)))/sqrt(length(cv1_index(:,inr))) (nanstd(cs1_index(:,inr)-is1_index(:,inr)))/sqrt(length(cs1_index(:,inr))) (nanstd(cm1_index(:,inr)-im1_index(:,inr)))/sqrt(length(cm1_index(:,inr))) ]...
    , 'LineStyle', 'none', ... 
        'Color', 'k', 'LineWidth', 1);hold on;
r=1;rng(1);r1 = r-0.01 + (0.2)*rand(length(cv1_index(:,inr)-iv1_index(:,inr)),1);
sc1=scatter(r1,cv1_index(:,inr)-iv1_index(:,inr),10,'k','filled');hold on;sc1.MarkerEdgeColor=[1 1 1];sc1.MarkerEdgeAlpha=0.5;sc1.MarkerFaceColor=[0.2 0.2 0.2];
r=2;rng(1);r1 = r-0.01 + (0.2)*rand(length(cs1_index(:,inr)-is1_index(:,inr)),1);
sc1=scatter(r1,cs1_index(:,inr)-is1_index(:,inr),10,'k','filled');hold on;sc1.MarkerEdgeColor=[1 1 1];sc1.MarkerEdgeAlpha=0.5;sc1.MarkerFaceColor=[0.2 0.2 0.2];
r=3;rng(1);r1 = r-0.01 + (0.2)*rand(length(cm1_index(:,inr)-im1_index(:,inr)),1);
sc1=scatter(r1,cm1_index(:,inr)-im1_index(:,inr),10,'k','filled');hold on;sc1.MarkerEdgeColor=[1 1 1];sc1.MarkerEdgeAlpha=0.5;sc1.MarkerFaceColor=[0.2 0.2 0.2];
xticklabels({'VISp','SSp-bfd','MOp'});ylim([-0.01 0.25])
set(gca,'FontSize',11);set(gca,'TickDir','out');box off;xtickangle(45);
[u p1]=ttest(cv1_index(:,inr),iv1_index(:,inr))
[u p2]=ttest(cs1_index(:,inr),is1_index(:,inr))
[u p3]=ttest(cm1_index(:,inr),im1_index(:,inr))
ylabel('L6a / (L23+L5+L6ab)');
hold on;text(1,max(get(gca,'YLim')),['*'],'FontSize',18);hold on;text(2,max(get(gca,'YLim')),['*'],'FontSize',18);hold on;text(3,max(get(gca,'YLim')),['*'],'FontSize',18);
%% Plot all cell numbers ipsi contra per injection type
%all super imposed
temp_color=[v1_color ;s1_color; m1_color];
module_names={'VISp','SSp-bfd','MOp'};
all_dati=[];all_dati=[squeeze(nansum(nansum(i_v1aa)))./(squeeze(nansum(nansum(i_v1aa)))+squeeze(nansum(nansum(c_v1aa))))...
     squeeze(nansum(nansum(i_s1aa)))./(squeeze(nansum(nansum(i_s1aa)))+squeeze(nansum(nansum(c_s1aa))))...
     squeeze(nansum(nansum(i_m1aa)))./(squeeze(nansum(nansum(i_m1aa)))+squeeze(nansum(nansum(c_m1aa))))];
all_datc=[];all_datc=[squeeze(nansum(nansum(c_v1aa)))./(squeeze(nansum(nansum(i_v1aa)))+squeeze(nansum(nansum(c_v1aa))))...
     squeeze(nansum(nansum(c_s1aa)))./(squeeze(nansum(nansum(i_s1aa)))+squeeze(nansum(nansum(c_s1aa))))...
     squeeze(nansum(nansum(c_m1aa)))./(squeeze(nansum(nansum(i_m1aa)))+squeeze(nansum(nansum(c_m1aa))))];
%plot
fig3= figure;set(fig3, 'Name', 'Paired comp');set(fig3, 'Position', [200, 300, 200, 250]);set(gcf,'color','w');
for j=1:3
dat=[];dat=[all_dati(:,j) all_datc(:,j)];
for i=1:length(data)
     pl=plot([1,2],[dat(:,1),dat(:,2)],'color',[0.5 0.5 0.5]);    
end
hold on;pS=plotSpread([dat(:,1),dat(:,2)],'categoryIdx',[ones(1,length(dat(:,1)))' ones(1,length(dat(:,2)))'*2],...
    'categoryMarkers',{'o','o'},'categoryColors',{temp_color(j,:), temp_color(j,:)});hold on;
hold on;er1=errorbar([0.75],nanmean(dat(:,1)),nanstd(dat(:,1),[],1)/sqrt(length(dat(:,1))),'ok','MarkerFaceColor',temp_color(j,:),'Markersize',8);
hold on;er2=errorbar([2.25],nanmean(dat(:,2)),nanstd(dat(:,2),[],1)/sqrt(length(dat(:,1))),'ok','MarkerFaceColor',temp_color(j,:),'Markersize',8);
text(2.2,1.1-0.07*j,module_names{j},'Color',temp_color(j,:))
end
xticklabels({'ipsi','contra'});ylabel('Fraction');hold on;title([]);set(gca,'FontSize',11);set(gca,'TickDir','out');box off;
%% Plot all cells ipis+contra  TOTAL numbers for supplement
temp_color=[v1_color ;s1_color; m1_color];
p1=[];p1=squeeze(nansum(nansum(i_v1aa)))+squeeze(nansum(nansum(c_v1aa)));
p2=[];p2=squeeze(nansum(nansum(i_s1aa)))+squeeze(nansum(nansum(c_s1aa)));
p3=[];p3=squeeze(nansum(nansum(i_m1aa)))+squeeze(nansum(nansum(c_m1aa)));
temp_p= [p1 p2 p3];
fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 500, 170, 270]);set(gcf,'color','w');tiledlayout("horizontal");
%title(['Total cell nr: ' num2str(nansum(nansum(temp_p)))],'FontWeight','normal')
title(['T: 3.17*10^6' ],'FontWeight','normal','FontSize',5)
for i=1:3
    hold on;
    b1=bar(i,nanmean(temp_p(:,i)));hold on;b1.FaceColor=temp_color(i,:);
      r=i;rng(i);r1 = r-0.01 + (0.2)*rand(length(temp_p(:,i)),1);
sc1=scatter(r1,temp_p(:,i),12,'k','filled');hold on;sc1.MarkerEdgeColor=[1 1 1];sc1.MarkerEdgeAlpha=0.5;sc1.MarkerFaceColor=[0.2 0.2 0.2];
    hold on;errorbar([i],[nanmean(temp_p(:,i))],[nanstd(temp_p(:,i))/sqrt(length(temp_p(:,i)))]...
    , 'LineStyle', 'none', ... 
        'Color', 'k', 'LineWidth', 1);hold on;
end
xticks([1:3]);hold on;box off;xticklabels({'VISp','SSp-bfd','MOp'});ylabel('Total nr of cells');
set(gca,'FontSize',11);set(gca,'TickDir','out');box off;xtickangle(45);ylim([0 410000]);
[p,tbl,stats] = anova1([temp_p])
presults = multcompare(stats)

%% Homotopic vs heterotoppic areas 

%% 
%relative contra and relative ipsi / contra mean
v1a_r=[nanmean(nansum(i_v1aa(:,:,:))./nansum(nansum(i_v1aa(:,:,:))),3)' nanmean(nansum(c_v1aa(:,:,:))./nansum(nansum(c_v1aa(:,:,:))),3)'];
s1a_r=[nanmean(nansum(i_s1aa(:,:,:))./nansum(nansum(i_s1aa(:,:,:))),3)' nanmean(nansum(c_s1aa(:,:,:))./nansum(nansum(c_s1aa(:,:,:))),3)'];
m1a_r=[nanmean(nansum(i_m1aa(:,:,:))./nansum(nansum(i_m1aa(:,:,:))),3)' nanmean(nansum(c_m1aa(:,:,:))./nansum(nansum(c_m1aa(:,:,:))),3)'];
%cell density and relative ipsi / contra mean
v1a_rcd=[nanmean(nansum(i_v1aacd(:,:,:))./nansum(nansum(i_v1aacd(:,:,:))),3)' nanmean(nansum(c_v1aacd(:,:,:))./nansum(nansum(c_v1aacd(:,:,:))),3)'];
s1a_rcd=[nanmean(nansum(i_s1aacd(:,:,:))./nansum(nansum(i_s1aacd(:,:,:))),3)' nanmean(nansum(c_s1aacd(:,:,:))./nansum(nansum(c_s1aacd(:,:,:))),3)'];
m1a_rcd=[nanmean(nansum(i_m1aacd(:,:,:))./nansum(nansum(i_m1aacd(:,:,:))),3)' nanmean(nansum(c_m1aacd(:,:,:))./nansum(nansum(c_m1aacd(:,:,:))),3)'];
%absolute numbers mean
v1a_t=[nanmean(nansum(i_v1aa(:,:,:)),3)' nanmean(nansum(c_v1aa(:,:,:)),3)'];
s1a_t=[nanmean(nansum(i_s1aa(:,:,:)),3)' nanmean(nansum(c_s1aa(:,:,:)),3)'];
m1a_t=[nanmean(nansum(i_m1aa(:,:,:)),3)' nanmean(nansum(c_m1aa(:,:,:)),3)'];
%relative contra and relative ipsi SEM
v1a_r_sem=[nanstd(nansum(i_v1aa(:,:,:))./nansum(nansum(i_v1aa(:,:,:))),[],3)./sqrt(length(i_v1a));...
    nanstd(nansum(c_v1aa(:,:,:))./nansum(nansum(c_v1aa(:,:,:))),[],3)./sqrt(length(c_v1a))]';
s1a_r_sem=[nanstd(nansum(i_s1aa(:,:,:))./nansum(nansum(i_s1aa(:,:,:))),[],3)./sqrt(length(i_s1a));...
    nanstd(nansum(c_s1aa(:,:,:))./nansum(nansum(c_s1aa(:,:,:))),[],3)./sqrt(length(c_s1a))]';
m1a_r_sem=[nanstd(nansum(i_m1aa(:,:,:))./nansum(nansum(i_m1aa(:,:,:))),[],3)./sqrt(length(i_m1a));...
    nanstd(nansum(c_m1aa(:,:,:))./nansum(nansum(c_m1aa(:,:,:))),[],3)./sqrt(length(c_m1a))]';
% total numbers per animal 
v1a_tca=squeeze(nansum(c_v1aa(:,:,:)));
v1a_tia=squeeze(nansum(i_v1aa(:,:,:)));
s1a_tca=squeeze(nansum(c_s1aa(:,:,:)));
s1a_tia=squeeze(nansum(i_s1aa(:,:,:)));
m1a_tca=squeeze(nansum(c_m1aa(:,:,:)));
m1a_tia=squeeze(nansum(i_m1aa(:,:,:)));
% relatives numbers per animal
v1a_rca=squeeze(nansum(c_v1aa(:,:,:))./nansum(nansum(c_v1aa(:,:,:))));
v1a_ria=squeeze(nansum(i_v1aa(:,:,:))./nansum(nansum(i_v1aa(:,:,:))));
s1a_rca=squeeze(nansum(c_s1aa(:,:,:))./nansum(nansum(c_s1aa(:,:,:))));
s1a_ria=squeeze(nansum(i_s1aa(:,:,:))./nansum(nansum(i_s1aa(:,:,:))));
m1a_rca=squeeze(nansum(c_m1aa(:,:,:))./nansum(nansum(c_m1aa(:,:,:))));
m1a_ria=squeeze(nansum(i_m1aa(:,:,:))./nansum(nansum(i_m1aa(:,:,:))));
%% add total numbers in a strcuture/csv for MT
ipsi_contra_cortex =table(cortex_names(:,2),v1a_t(:,1)',v1a_t(:,2)',s1a_t(:,1)',s1a_t(:,2)',m1a_t(:,1)',m1a_t(:,2)',...
    'variablenames',{'Area','V1ipsi','V1contra','S1ipsi','S1contra','M1ipsi','M1contra'});
writetable(ipsi_contra_cortex,'ipsi_contra_cortex.csv')
%% Homotopic areas fractions on contralateral side: SUpplement
%visp_idx=[31];ssp_idx=[18];mop_idx=[25]
fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 500, 170, 250]);set(gcf,'color','w');tiledlayout("horizontal")
b1=bar(1,v1a_r(visp_idx,2));hold on;b1.FaceColor=v1_color;
b2=bar(2,s1a_r(ssp_idx,2));hold on;b2.FaceColor=s1_color;
b3=bar(3,m1a_r(mop_idx,2));b3.FaceColor=m1_color;
xticks([1:3]);hold on;
hold on;errorbar([1 2 3],[v1a_r(visp_idx,2) s1a_r(ssp_idx,2) m1a_r(mop_idx,2)],[v1a_r_sem(visp_idx,2) s1a_r_sem(ssp_idx,2) m1a_r_sem(mop_idx,2)]...
    , 'LineStyle', 'none', ... 
        'Color', 'k', 'LineWidth', 1);hold on;
r=1;rng(1);r1 = r-0.01 + (0.2)*rand(length(v1a_rca(visp_idx,:)),1);
sc1=scatter(r1,v1a_rca(31,:),15,'k','filled');hold on;sc1.MarkerEdgeColor=[1 1 1];sc1.MarkerEdgeAlpha=0.5;sc1.MarkerFaceColor=[0.2 0.2 0.2];
 r=2;rng(1);r1 = r-0.01 + (0.2)*rand(length(s1a_rca(ssp_idx,:)),1);
 sc1=scatter(r1,s1a_rca(18,:),15,'k','filled');hold on;sc1.MarkerEdgeColor=[1 1 1];sc1.MarkerEdgeAlpha=0.5;sc1.MarkerFaceColor=[0.2 0.2 0.2];
 r=3;rng(1);r1 = r-0.01 + (0.2)*rand(length(m1a_rca(mop_idx,:)),1);
 sc1=scatter(r1,m1a_rca(25,:),15,'k','filled');hold on;sc1.MarkerEdgeColor=[1 1 1];sc1.MarkerEdgeAlpha=0.5;sc1.MarkerFaceColor=[0.2 0.2 0.2];
xticklabels({'VISp','SSp-bfd','MOp'});%ylim([-0.05 0.75]);
ylabel('Contra Homotopic fraction')
set(gca,'FontSize',11);set(gca,'TickDir','out');box off;xtickangle(45);
 hold on;line([2 3],[max(get(gca,'YLim')) max(get(gca,'YLim'))],'Color','k');
 hold on;text(2,max(get(gca,'YLim'))+0.03,['***'],'FontSize',18);
  hold on;line([1 3],[max(get(gca,'YLim')+0.1) max(get(gca,'YLim'))+0.1],'Color','k');
  hold on;text(1.5,max(get(gca,'YLim'))-0.06,['***'],'FontSize',18);
[p,tbl,stats] = anova1([v1a_rca(visp_idx,:)' s1a_rca(ssp_idx,:)' m1a_rca(mop_idx,:)'])
 presults = multcompare(stats)
%% Percentile criteria for all total counts using averages
%define percentile
prct=12;
%add all average counts across areas and animals 
all_total=[];all_total=[v1a_t s1a_t m1a_t];
%Plot 
[temp_s kk]=sort(all_total(:),'descend');
fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 500, 500, 400]);set(gcf,'color','w');tiledlayout("horizontal");
plot(temp_s,'-o','Color','k');box off;
idx_prct=[];idx_prct=find(temp_s>prctile(temp_s,prct));
hold on;plot(idx_prct(end),prctile(temp_s,prct)+2000,'v');
hold on;text(idx_prct(end),prctile(temp_s,prct)+3000,[num2str(temp_s(idx_prct(end)))  ' prct=' num2str(prct)]);
ylabel('Total cell number');xlabel('sorted areas across all');set(gca,'FontSize',11);set(gca,'TickDir','out');box off;xtickangle(45);ylim([-100 max(temp_s)]);xlim([-5 length(temp_s)+1]);
%% Use this percentil cut off to set given areas to zero and for the injetion area on ipsi to NaN
%total count threshold 
thr_count=prctile(all_total(:),prct);
%use function zero_area.m to set to zero 
[v1a_tcam v1a_rcam c_v1aam]=zero_area(v1a_tca, v1a_rca,c_v1aaf,thr_count);
[v1a_tiam v1a_riam i_v1aam]=zero_area(v1a_tia, v1a_ria,i_v1aaf,thr_count);
[s1a_tcam s1a_rcam c_s1aam]=zero_area(s1a_tca, s1a_rca,c_s1aaf,thr_count);
[s1a_tiam s1a_riam i_s1aam]=zero_area(s1a_tia, s1a_ria,i_s1aaf,thr_count);
[m1a_tcam m1a_rcam c_m1aam]=zero_area(m1a_tca, m1a_rca,c_m1aaf,thr_count);
[m1a_tiam m1a_riam i_m1aam]=zero_area(m1a_tia, m1a_ria,i_m1aaf,thr_count);
%replace the zeros with NaNs for the injection areas IPSI
%ipsi absolute
v1a_tiam(visp_idx,:)=ones(1,size(v1a_tiam,2))*NaN;
s1a_tiam(ssp_idx,:)=ones(1,size(s1a_tiam,2))*NaN;
m1a_tiam(mop_idx,:)=ones(1,size(m1a_tiam,2))*NaN;
%ispi relative
v1a_riam(visp_idx,:)=ones(1,size(v1a_riam,2))*NaN;
s1a_riam(ssp_idx,:)=ones(1,size(s1a_riam,2))*NaN;
m1a_riam(mop_idx,:)=ones(1,size(m1a_riam,2))*NaN;
%layers relative
for i=1:size(i_v1aam,3)
    i_v1aam(:,visp_idx,i)=ones(1,size(i_v1aam,1))*NaN;
end
for i=1:size(i_s1aam,3)
    i_s1aam(:,ssp_idx,i)=ones(1,size(i_s1aam,1))*NaN;
end
for i=1:size(i_m1aam,3)
    i_m1aam(:,mop_idx,i)=ones(1,size(i_m1aam,1))*NaN;
end
%% Plot the numbers as fraction per area where L2/3, L5 or L6a is dominant
%all areas stacked
cat_temp=[];cat_temp=cat(3,i_v1aam,i_s1aam,i_m1aam);
ind=[];
for i=1:size(cat_temp,3)
[X, ind(:,i)]=nanmax(cat_temp(:,:,i),[],1);
end
p1=[];p1=sum(ind==2)./sum(ind>1);
p2=[];p2=sum(ind==4)./sum(ind>1);
p3=[];p3=sum(ind==5)./sum(ind>1);
temp_p=[];temp_p= [p1' p2' p3'];


fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 500, 150, 270]);set(gcf,'color','w');
b1=bar([1],[nanmean(temp_p(:,1)) nanmean(temp_p(:,2)) nanmean(temp_p(:,3))],'stacked');hold on;
b1(1).FaceColor=[1 1 1];b1(2).FaceColor=[0.7 0.7 0.7];b1(3).FaceColor=[0.55 0.55 0.55];

hold on;errorbar([1],[nanmean(temp_p(:,1)) nanmean(temp_p(:,2))+nanmean(temp_p(:,1)) nanmean(temp_p(:,3))+nanmean(temp_p(:,2))+nanmean(temp_p(:,1))],[nanstd(temp_p(:,1))/sqrt(length(temp_p(:,1))) nanstd(temp_p(:,2))/sqrt(length(temp_p(:,2))) nanstd(temp_p(:,3))/sqrt(length(temp_p(:,3)))]...
      , 'LineStyle', 'none', ... 
         'Color', 'k', 'LineWidth', 1);hold on;

text(0.75,0.75,num2str(round(nanmean(temp_p(:,3)),3)*100),'FontSize',8,'Color','w');
text(0.75,0.4,num2str(round(nanmean(temp_p(:,2)),3)*100),'FontSize',8,'Color','w');
text(0.75,0.1,num2str(round(nanmean(temp_p(:,1)),3)*100),'FontSize',8);

cat_temp=[];cat_temp=cat(3,c_v1aam,c_s1aam,c_m1aam);
ind=[];
for i=1:size(cat_temp,3)
[X, ind(:,i)]=nanmax(cat_temp(:,:,i),[],1);
end
p1=[];p1=sum(ind==2)./sum(ind>1);
p2=[];p2=sum(ind==4)./sum(ind>1);
p3=[];p3=sum(ind==5)./sum(ind>1);
temp_p=[];temp_p= [p1' p2' p3'];

hold on;
b1=bar([2],[nanmean(temp_p(:,1)) nanmean(temp_p(:,2)) nanmean(temp_p(:,3))],'stacked');hold on;
b1(1).FaceColor=[1 1 1];b1(2).FaceColor=[0.7 0.7 0.7];b1(3).FaceColor=[0.55 0.55 0.55];

hold on;errorbar([2],[nanmean(temp_p(:,1)) nanmean(temp_p(:,2))+nanmean(temp_p(:,1)) nanmean(temp_p(:,3))+nanmean(temp_p(:,2))+nanmean(temp_p(:,1))],[nanstd(temp_p(:,1))/sqrt(length(temp_p(:,1))) nanstd(temp_p(:,2))/sqrt(length(temp_p(:,2))) nanstd(temp_p(:,3))/sqrt(length(temp_p(:,3)))]...
      , 'LineStyle', 'none', ... 
         'Color', 'k', 'LineWidth', 1);hold on;

hold on;box off;xticklabels({[]});ylabel('Area layer dominance');
set(gca,'FontSize',11);set(gca,'TickDir','out');


text(1.75,0.75,num2str(round(nanmean(temp_p(:,3)),3)*100),'FontSize',8,'Color','w');
text(1.75,0.3,num2str(round(nanmean(temp_p(:,2)),3)*100),'FontSize',8,'Color','w');
text(1.75,0.05,num2str(round(nanmean(temp_p(:,1)),3)*100),'FontSize',8);


t1=text(2.7,0.8,'L6a','FontSize',11);
set(t1,'Rotation',270);
t1=text(2.7,0.35,'L5','FontSize',11);
set(t1,'Rotation',270);
t1=text(2.7,0.1,'L2/3','FontSize',11);
set(t1,'Rotation',270);
set(gca,'xtick',[])
set(gca,'XColor', 'none')
set(gca,'visible','off')
%Overall
hold on;
text(1.2,1.15,'All','FontSize',11)
text(0.6,1.08,'ipsi','FontSize',11)
text(1.6,1.08,'contra','FontSize',11)
%% Based on area VISp, SSpbf, MOp dominace of L2/3 L5 and L6a
ind=[];temp_c=m1_color;temp_area=[];temp_area=i_v1aam;
for i=1:size(temp_area,3)
[X, ind(:,i)]=nanmax(temp_area(:,:,i),[],1);
end
p1=[];p1=(sum(ind==2)./sum(ind>1))*100;
p2=[];p2=(sum(ind==4)./sum(ind>1))*100;
p3=[];p3=(sum(ind==5)./sum(ind>1))*100;
temp_p= [p1' p2' p3'];
fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 500, 150, 250]);set(gcf,'color','w');tiledlayout("horizontal");
title('MOp','FontWeight','normal')
for i=1:3
    hold on;
    b1=bar(i,nanmean(temp_p(:,i)));hold on;b1.FaceColor=temp_c;
      r=i;rng(i);r1 = r-0.01 + (0.2)*rand(length(temp_p(:,i)),1);
sc1=scatter(r1,temp_p(:,i),12,'k','filled');hold on;sc1.MarkerEdgeColor=[1 1 1];sc1.MarkerEdgeAlpha=0.5;sc1.MarkerFaceColor=[0.2 0.2 0.2];
    hold on;errorbar([i],[nanmean(temp_p(:,i))],[nanstd(temp_p(:,i))/sqrt(length(temp_p(:,i)))]...
    , 'LineStyle', 'none', ... 
        'Color', 'k', 'LineWidth', 1);hold on;
end
xticks([1:3]);hold on;box off;xticklabels({'L2/3','L5','L6a'});ylabel('Layer area dominance (%)');
set(gca,'FontSize',11);set(gca,'TickDir','out');box off;xtickangle(45);
% [p,tbl,stats] = anova1([temp_p])
% presults = multcompare(stats)

%% simply show fraction of area stacked for ipsi - contra 
temp_color=[v1_color ;s1_color; m1_color];
p1=[];p1=sum(v1a_tiam>0);
p2=[];p2=sum(s1a_tiam>0);
p3=[];p3=sum(m1a_tiam>0);
p4=[];p4=sum(v1a_tcam>0);
p5=[];p5=sum(s1a_tcam>0);
p6=[];p6=sum(m1a_tcam>0);
temp_p= [p1' p2' p3'];temp_p2= [p4' p5' p6'];
fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 500, 220, 270]);set(gcf,'color','w');tiledlayout("horizontal");

for i=1:3
    hold on;
    b1=bar([i*2],[nanmean(temp_p2(:,i)) nanmean(temp_p(:,i))],'stacked');hold on;b1(2).FaceColor=[0.8 0.8 0.8];b1(1).FaceColor=[0.3 0.3 0.3];b1(2).EdgeColor=temp_color(i,:);b1(1).EdgeColor=temp_color(i,:);
   

     hold on;errorbar([i*2],[nanmean(temp_p2(:,i)) nanmean(temp_p(:,i))+nanmean(temp_p2(:,i))],[nanstd(temp_p2(:,i))/sqrt(length(temp_p2(:,i))) nanstd(temp_p(:,i))/sqrt(length(temp_p(:,i)))]...
     , 'LineStyle', 'none', ... 
         'Color', 'k', 'LineWidth', 1);hold on;
end

xticks([2:2:6]);hold on;box off;xticklabels({[]});ylabel('Nr of areas');

% ispi contra and numbers text
text(1.7,20,num2str(round(nanmean(p4)),2),'FontSize',10,'Color','w');
text(1.7,55,num2str(round(nanmean(p1)),2),'FontSize',10);
t1=text(2.85,23,'contra','FontSize',11);
set(t1,'Rotation',270);
t1=text(2.85,60,'ipsi','FontSize',11);
set(t1,'Rotation',270);

text(3.7,22,num2str(round(nanmean(p5)),2),'FontSize',10,'Color','w');
text(3.7,62,num2str(round(nanmean(p2)),2),'FontSize',10);
t1=text(4.85,25,'contra','FontSize',11);
set(t1,'Rotation',270);
t1=text(4.85,65,'ipsi','FontSize',11);
set(t1,'Rotation',270);

text(5.7,15,num2str(round(nanmean(p6)),2),'FontSize',10,'Color','w');
text(5.7,45,num2str(round(nanmean(p3)),2),'FontSize',10);
t1=text(6.85,23,'contra','FontSize',11);
set(t1,'Rotation',270);
t1=text(6.85,47,'ipsi','FontSize',11);
set(t1,'Rotation',270);
set(gca,'FontSize',11);set(gca,'TickDir','out');box off;xtickangle(45);
h = gca; h.XAxis.Visible = 'off';
ylim([-1 90])
[p,tbl,stats] = anova1([temp_p])
presults = multcompare(stats)


%% Correlation between hemisphere areas R2 square avergae across injection areas
for i=1:size(v1a_riam,2)
idx = find(~isnan(v1a_riam(:,i)));
mdl = fitlm(v1a_riam(idx,i),v1a_rcam(idx,i));
coefs = polyfit(v1a_riam(idx,i)', v1a_rcam(idx,i)', 1);
slope_v1(i)=coefs(1);
r2_sqv1(i)=mdl.Rsquared.Ordinary;
end
for i=1:size(s1a_riam,2)
idx = find(~isnan(s1a_riam(:,i)));
mdl = fitlm(s1a_riam(idx,i),s1a_rcam(idx,i));
coefs = polyfit(s1a_riam(idx,i)', s1a_rcam(idx,i)', 1);
slope_s1(i)=coefs(1);
r2_sqs1(i)=mdl.Rsquared.Ordinary;
end
for i=1:size(m1a_riam,2)
idx = find(~isnan(m1a_riam(:,i)));
mdl = fitlm(m1a_riam(idx,i),m1a_rcam(idx,i));
coefs = polyfit(m1a_riam(idx,i)', m1a_rcam(idx,i)', 1);
slope_m1(i)=coefs(1);
r2_sqm1(i)=mdl.Rsquared.Ordinary;
end
 tm1=r2_sqv1;
 tm2=r2_sqs1;
 tm3=r2_sqm1;
%tm1=slope_v1;
%tm2=slope_s1;
%tm3=slope_m1;

fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 500, 170, 250]);set(gcf,'color','w');tiledlayout("horizontal")
b1=bar(1,nanmean(tm1));hold on;b1.FaceColor=v1_color;
b2=bar(2,nanmean(tm2));hold on;b2.FaceColor=s1_color;
b3=bar(3,nanmean(tm3));b3.FaceColor=m1_color;
r=1;rng(1);r1 = r-0.01 + (0.2)*rand(length(tm1),1);
sc1=scatter(r1,tm1,15,'k','filled');hold on;sc1.MarkerEdgeColor=[1 1 1];sc1.MarkerEdgeAlpha=0.5;sc1.MarkerFaceColor=[0.2 0.2 0.2];
 r=2;rng(1);r1 = r-0.01 + (0.2)*rand(length(tm2),1);
 sc1=scatter(r1,tm2,15,'k','filled');hold on;sc1.MarkerEdgeColor=[1 1 1];sc1.MarkerEdgeAlpha=0.5;sc1.MarkerFaceColor=[0.2 0.2 0.2];
 r=3;rng(1);r1 = r-0.01 + (0.2)*rand(length(tm3),1);
 sc1=scatter(r1,tm3,15,'k','filled');hold on;sc1.MarkerEdgeColor=[1 1 1];sc1.MarkerEdgeAlpha=0.5;sc1.MarkerFaceColor=[0.2 0.2 0.2];
xticks([1:3]);hold on;
hold on;errorbar([1 2 3],[nanmean(tm1) nanmean(tm2) nanmean(tm3)],[nanstd(tm1)/sqrt(size(tm1,2)) nanstd(tm2)/sqrt(size(tm2,2)) nanstd(tm3)/sqrt(size(tm3,2))]...
    , 'LineStyle', 'none', ... 
        'Color', 'k', 'LineWidth', 1);hold on;
xticklabels({'VISp','SSpbf','MOp'});%ylim([-0.05 0.75]);
ylabel('R2 ipsi contra')
set(gca,'FontSize',11);set(gca,'TickDir','out');box off;xtickangle(45);
[p,tbl,stats] = anova1([tm1' tm2' tm3'])
presults = multcompare(stats)

%% Heterotopic/homotopic areas categroreis 
rm1=[];rm2=[];
rm1=m1a_rcam;
p1=[];p1=nansum(rm1(frontal_idx,:));
p2=[];p2=nansum(rm1(lateral_idx,:));
p3=[];p3=nansum(rm1(somamo_idx,:));
p4=[];p4=nansum(rm1(visual_idx,:));
p5=[];p5=nansum(rm1(medial_idx,:));
p6=[];p6=nansum(rm1(aud_idx,:));
temp_p= [p1' p2' p3' p4' p5' p6'];
temp_c=m1_color;
fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 500, 200, 250]);set(gcf,'color','w');tiledlayout("horizontal");
title('MOp contra','FontWeight','normal')
for i=1:6
    hold on;
    b1=bar(i,nanmean(temp_p(:,i)));hold on;b1.FaceColor=temp_c;
      r=i;rng(i);r1 = r-0.01 + (0.2)*rand(length(temp_p(:,i)),1);
sc1=scatter(r1,temp_p(:,i),12,'k','filled');hold on;sc1.MarkerEdgeColor=[1 1 1];sc1.MarkerEdgeAlpha=0.5;sc1.MarkerFaceColor=[0.2 0.2 0.2];
    hold on;errorbar([i],[nanmean(temp_p(:,i))],[nanstd(temp_p(:,i))/sqrt(length(temp_p(:,i)))]...
    , 'LineStyle', 'none', ... 
        'Color', 'k', 'LineWidth', 1);hold on;

end
xticks([1:6]);hold on;box off;xticklabels({'Fron','Lat','SoMo','Vis','Med','Aud'});ylabel('Fraction of cells');ylim([0 1]);
set(gca,'FontSize',11);set(gca,'TickDir','out');box off;xtickangle(45);
[p,tbl,stats] = anova1([temp_p])
presults = multcompare(stats)
%% Heterotopic/homotopic areas categroreis 
 %use function hemisphere_di.m
rm1=[];rm2=[];
rm1=m1a_rcam;rm2=m1a_riam;
temp_c=v1_color;
p1=[];[p1]=hemisphere_di(rm1(frontal_idx,:),rm2(frontal_idx,:));p1=nanmean(p1);
p2=[];[p2]=hemisphere_di(rm1(lateral_idx,:),rm2(lateral_idx,:));p2=nanmean(p2);
p3=[];[p3]=hemisphere_di(rm1(somamo_idx,:),rm2(somamo_idx,:));p3=nanmean(p3);
p4=[];[p4]=hemisphere_di(rm1(visual_idx,:),rm2(visual_idx,:));p4=nanmean(p4);
p5=[];[p5]=hemisphere_di(rm1(medial_idx,:),rm2(medial_idx,:));p5=nanmean(p5);
p6=[];[p6]=hemisphere_di(rm1(aud_idx,:),rm2(aud_idx,:));p6=nanmean(p6);
temp_p= [p1' p2' p3' p4' p5' p6'];
fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 500, 200, 250]);set(gcf,'color','w');tiledlayout("horizontal");
title('VISp','FontWeight','normal')
for i=1:6
    hold on;
    b1=bar(i,nanmean(temp_p(:,i)));hold on;b1.FaceColor=temp_c;
      r=i;rng(i);r1 = r-0.01 + (0.2)*rand(length(temp_p(:,i)),1);
sc1=scatter(r1,temp_p(:,i),12,'k','filled');hold on;sc1.MarkerEdgeColor=[1 1 1];sc1.MarkerEdgeAlpha=0.5;sc1.MarkerFaceColor=[0.2 0.2 0.2];
    hold on;errorbar([i],[nanmean(temp_p(:,i))],[nanstd(temp_p(:,i))/sqrt(length(temp_p(:,i)))]...
    , 'LineStyle', 'none', ... 
        'Color', 'k', 'LineWidth', 1);hold on;
end
xticks([1:6]);hold on;box off;xticklabels({'Fron','Lat','SoMo','Vis','Med','Aud'});ylabel('Hemispheric DI');
set(gca,'FontSize',11);set(gca,'TickDir','out');box off;xtickangle(45);
[p,tbl,stats] = anova1([temp_p])
presults = multcompare(stats)
%% 
%% ILN all
dat=[];dat=[i_index(:,1) c_index(:,1)];cl={[0.8 0.8 0.8],[0.3 0.3 0.3]};
paired_plot_box(dat,cl);
xticklabels({'ipsi','contra'});xtickangle(45);%ylim([0.08 0.6]);
%pvalue of paired signrank test: 
%[p1]=signrank(dat(:,1),dat(:,2))
[u p1]=ttest(dat(:,1),dat(:,2))
hold on;text(1.25,0.9,['***'],'FontSize',18);
%ylabel({'ILN'; '(L5 + L6ab) / (L2/3 + L5 + L6ab)'})
ylabel('ILN');
title([]);set(gca,'FontSize',11);set(gca,'TickDir','out');
%% ILN with color code
temp_color=[v1_color ;s1_color; m1_color];
module_names={'VISp','SSp-bfd','MOp'};
all_dati=[];all_dati=[iv1_index(:,1) is1_index(:,1) im1_index(:,1)];
all_datc=[];all_datc=[cv1_index(:,1) cs1_index(:,1) cm1_index(:,1)];
%plot
fig3= figure;set(fig3, 'Name', 'Paired comp');set(fig3, 'Position', [200, 300, 170, 250]);set(gcf,'color','w');
for j=1:3
dat=[];dat=[all_dati(:,j) all_datc(:,j)];
for i=1:length(data)
     pl=plot([1,2],[dat(:,1),dat(:,2)],'color',[0.5 0.5 0.5]);    
end
hold on;pS=plotSpread([dat(:,1),dat(:,2)],'categoryIdx',[ones(1,length(dat(:,1)))' ones(1,length(dat(:,2)))'*2],...
    'categoryMarkers',{'o','o'},'categoryColors',{temp_color(j,:), temp_color(j,:)});hold on;
hold on;er1=errorbar([0.75],nanmean(dat(:,1)),nanstd(dat(:,1),[],1)/sqrt(length(dat(:,1))),'ok','MarkerFaceColor',temp_color(j,:),'Markersize',8);
hold on;er2=errorbar([2.25],nanmean(dat(:,2)),nanstd(dat(:,2),[],1)/sqrt(length(dat(:,1))),'ok','MarkerFaceColor',temp_color(j,:),'Markersize',8);
%text(2.2,1.1-0.07*j,module_names{j},'Color',temp_color(j,:))
end
xticklabels({'ipsi','contra'});ylabel('ILN');hold on;title([]);xtickangle(45);set(gca,'FontSize',11);set(gca,'TickDir','out');box off;

%% correlation between indexes per area  ILN
%% all correlations superimposed on each other 
idx=[];idx=1;
temp_all={i_v1aam c_v1aam; i_s1aam c_s1aam; i_m1aam c_m1aam};
temp_color=[v1_color ;s1_color; m1_color];
%plot
fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 500, 350, 350]);set(gcf,'color','w');
for i=1:3
  
    rm1=[];rm1=temp_all{i,1};rm2=[];rm2=temp_all{i,2};
[vali_all valdci_all val_i val_dci xerr yerr r_a p_a r_avg p_avg] = anatomy_correlation(rm1, rm2, idx);
 temp_c=temp_color(i,:)
 
hold on;errorbar(val_i,val_dci,yerr,yerr,xerr,xerr, 'LineStyle', 'none', 'Color', 'k', 'LineWidth', 0.5,'CapSize',0);hold on;
sc1=scatter(val_i,val_dci,50,'filled');sc1.MarkerFaceColor=temp_c;sc1.MarkerEdgeColor='w';
hold on;
hold on;text(0.5,0.7-(0.06*i),...
    ['r= ' num2str(round(nanmean(r_a),2)) ' +- ' num2str(round(nanstd(r_a)/sqrt(length(r_a)),2))],'FontSize',11,'Color',temp_c);
%hold on;text(0.5,0.5-(0.05*i),['p= ' num2str(round(nanmean(p_a),2))],'FontSize',11,'Color',temp_c);
end
%ylabel('L6d_c - L6d_i');xlabel('L6d_i');
ylabel('ILN_c - ILN_i');xlabel('ILN_i');
axis square;set(gca,'FontSize',11);set(gca,'TickDir','out');box off;


%% Sorted modules ILNc-ILNi
idx=[];idx=1;
temp_all={i_v1aam c_v1aam; i_s1aam c_s1aam; i_m1aam c_m1aam};
temp_color=[v1_color ;s1_color; m1_color];
idx_modules={frontal_idx lateral_idx somamo_idx visual_idx medial_idx aud_idx};
module_names={'Fron','Lat','SoMo','Vis','Med','Aud'};
%plot
for i=1:3
fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 500, 180, 250]);set(gcf,'color','w');  
rm1=[];rm1=temp_all{i,1};rm2=[];rm2=temp_all{i,2};
[vali_all valdci_all val_i val_dci xerr yerr r_a p_a r_avg p_avg] = anatomy_correlation(rm1, rm2, idx);
 temp_c=temp_color(i,:);
 hold on
 temp_p=[];p=[];tbl=[];stats=[];err_sem=[];order_err=[];
 for k=1:6
     temp_p(:,k)=nanmean(valdci_all(:,idx_modules{k}),2)
     sort_means(k)=nanmean(nanmean(valdci_all(:,idx_modules{k}),2));
     err_sem(k)=nanstd(temp_p(:,k))/(sqrt(length(~isnan(temp_p(:,k)))))
 end
 [order_m kk]=sort(sort_means,'descend');
order_err=err_sem(kk);
 for k=1:6
     
 hold on
 b1=bar(k,order_m(k));b1.FaceColor=temp_c;
%      r=k;rng(k);r1 = r-0.01 + (0.2)*rand(length(nanmean(valdci_all(:,idx_modules{k}),2)),1);
% sc1=scatter(r1,nanmean(valdci_all(:,idx_modules{k}),2),12,'k','filled');hold on;sc1.MarkerEdgeColor=[1 1 1];sc1.MarkerEdgeAlpha=0.5;sc1.MarkerFaceColor=[0.2 0.2 0.2];
     hold on;errorbar([k],[order_m(k)],[order_err(k)]...
    , 'LineStyle', 'none', ... 
         'Color', 'k', 'LineWidth', 1);hold on;
 end
 xticks([1:6]);hold on;box off;xticklabels(module_names(kk));set(gca,'FontSize',11);set(gca,'TickDir','out');box off;
 xtickangle(45);
 ylim([-0.14 0.3]);
 ylabel('ILNc - ILNi')
 [p,tbl,stats] = anova1([temp_p])
presults = multcompare(stats)

end

%% Sorted areas based on index
idx=[];idx=1;
rm1=[];rm2=[];rm1=i_v1aam;rm2=c_v1aam;
temp_c1=[0.8 0.8 0.8];
temp_c2=[0.3 0.3 0.3];
cl_idx=visp_idx;
title_pan='SSp-bfd';
%calculate and plot 
i_metric=[];c_metric=[];
%calculate indexes across areas
for i=1:length(rm1)
temp_l=[];temp_l2=[];
temp_l=squeeze(rm1(:,i,:));
temp_l2=squeeze(rm2(:,i,:));
[i_metric(:,i,:)] = anatomy_indexcalc(temp_l);
[c_metric(:,i,:)] = anatomy_indexcalc(temp_l2);
end 
temp_metric1=i_metric(:,:,idx);temp_metric2=c_metric(:,:,idx);
rm1=[];rm2=[];
rm1=temp_metric1;rm2=temp_metric2;

temp_rm1=[];temp_rm1=nanmean(rm1);
temp_rm2=[];temp_rm2=nanmean(rm2);
temp_rm1(cl_idx)=[];
temp_rm2(cl_idx)=[];
cortex_abb=cortex_names(:,2);
cortex_abb(cl_idx)=[];
yerr=[];yerr=nanstd(rm1);
yerr2=[];yerr2=nanstd(rm2);
yerr(cl_idx)=[];
yerr2(cl_idx)=[];
sort_d=[];kk=[];
[sort_d kk]=sort(temp_rm1,'ascend');

alpha1=1;
fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 500, 850, 250]);set(gcf,'color','w');
title(title_pan,'FontWeight','normal');
ylabel('ILN')
hold on;h=errorbar([1:44],[sort_d],[yerr/sqrt(length(sort_d(~isnan(sort_d))))]...
    , 'LineStyle', 'none', ... 
        'Color', temp_c1, 'LineWidth', 0.5,'CapSize',0);set([h.Bar, h.Line], 'ColorType', 'truecoloralpha', 'ColorData', [h.Line.ColorData(1:3); 255*alpha1]);
hold on;h=errorbar([1:44],[temp_rm2(kk)],[yerr2(kk)/sqrt(length(temp_rm2(~isnan(temp_rm2))))]...
    , 'LineStyle', 'none', ... 
        'Color', temp_c2, 'LineWidth', 0.5,'CapSize',0);
hold on;pp1=scatter(1:44,sort_d,60,'o');pp1.MarkerEdgeColor=[1 1 1];pp1.MarkerFaceAlpha=1;pp1.MarkerFaceColor=temp_c1;
hold on;pp1=scatter(1:44,temp_rm2(kk),60,'o');pp1.MarkerEdgeColor=[1 1 1];pp1.MarkerFaceAlpha=1;pp1.MarkerFaceColor=temp_c2;box off;
hold on;pp1=scatter(45,nanmean(temp_metric2(:,cl_idx)),60,'o');pp1.MarkerEdgeColor=[1 1 1];pp1.MarkerFaceAlpha=1;pp1.MarkerFaceColor=temp_c2;box off;
xticks([1:45]);ylim([min(get(gca,'YLim')) 1.2])
xticklabels([cortex_abb(kk) ;{title_pan}])
set(gca,'FontSize',11);set(gca,'TickDir','out');
%% Modules showing ILNi and ILNc
idx_modules={frontal_idx lateral_idx somamo_idx visual_idx medial_idx aud_idx};
module_names={'Fron','Lat','SoMo','Vis','Med','Aud'};
idx=[];idx=1;
rm1=[];rm2=[];rm1=i_v1aam;rm2=c_v1aam;

temp_c1=[0.8 0.8 0.8];
temp_c2=[0.3 0.3 0.3];


%calculate and plot 
i_metric=[];c_metric=[];
%calculate indexes across areas
for i=1:length(rm1)
temp_l=[];temp_l2=[];
temp_l=squeeze(rm1(:,i,:));
temp_l2=squeeze(rm2(:,i,:));
[i_metric(:,i,:)] = anatomy_indexcalc(temp_l);
[c_metric(:,i,:)] = anatomy_indexcalc(temp_l2);
end 
temp_metric1=i_metric(:,:,idx);temp_metric2=c_metric(:,:,idx);

fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 500, 400, 350]);set(gcf,'color','w');tiledlayout('vertical');
nexttile
for i=1:6
    hold on;
    yerr1=[];yerr1=nanstd(nanmean(temp_metric1(:,idx_modules{i}),2))/sqrt(length(nanstd(temp_metric1(:,idx_modules{i}),[],2)));
    yerr2=[];yerr2=nanstd(nanmean(temp_metric2(:,idx_modules{i}),2))/sqrt(length(nanstd(temp_metric2(:,idx_modules{i}),[],2)));
    ymean1=[];ymean1=nanmean(temp_metric1(:,idx_modules{i}),2);
    ymean2=[];ymean2=nanmean(temp_metric2(:,idx_modules{i}),2);
    hold on;h=errorbar([i],[nanmean(ymean2)],[yerr2], 'LineStyle', 'none','Color', temp_c1, 'LineWidth', 0.5,'CapSize',0);
    hold on;h=errorbar([i],[nanmean(ymean1)],[yerr1], 'LineStyle', 'none','Color', temp_c1, 'LineWidth', 1,'CapSize',0);
    pp1=scatter(i,nanmean(ymean1),50,'filled');pp1.MarkerEdgeColor=[1 1 1];pp1.MarkerFaceAlpha=1;pp1.MarkerFaceColor=temp_c1;
    pp1=scatter(i,nanmean(ymean2),50,'filled');pp1.MarkerEdgeColor=[1 1 1];pp1.MarkerFaceAlpha=1;pp1.MarkerFaceColor=temp_c2;box off;
end
xlim([0 7]);ylim([0.5 1])
xticks([]);;set(gca,'FontSize',11);set(gca,'TickDir','out');title('VISp','FontWeight','normal');
text(6.5,1.1,'ipsi','Color',temp_c1);text(6.5,1.02,'contra','Color',temp_c2);

rm1=[];rm2=[];rm1=i_s1aam;rm2=c_s1aam;
%calculate and plot 
i_metric=[];c_metric=[];
%calculate indexes across areas
for i=1:length(rm1)
temp_l=[];temp_l2=[];
temp_l=squeeze(rm1(:,i,:));
temp_l2=squeeze(rm2(:,i,:));
[i_metric(:,i,:)] = anatomy_indexcalc(temp_l);
[c_metric(:,i,:)] = anatomy_indexcalc(temp_l2);
end 
temp_metric1=i_metric(:,:,idx);temp_metric2=c_metric(:,:,idx);

nexttile
for i=1:6
    hold on;
    yerr1=[];yerr1=nanstd(nanmean(temp_metric1(:,idx_modules{i}),2))/sqrt(length(nanstd(temp_metric1(:,idx_modules{i}),[],2)));
    yerr2=[];yerr2=nanstd(nanmean(temp_metric2(:,idx_modules{i}),2))/sqrt(length(nanstd(temp_metric2(:,idx_modules{i}),[],2)));
    ymean1=[];ymean1=nanmean(temp_metric1(:,idx_modules{i}),2);
    ymean2=[];ymean2=nanmean(temp_metric2(:,idx_modules{i}),2);
    hold on;h=errorbar([i],[nanmean(ymean2)],[yerr2], 'LineStyle', 'none','Color', temp_c1, 'LineWidth', 0.5,'CapSize',0);
    hold on;h=errorbar([i],[nanmean(ymean1)],[yerr1], 'LineStyle', 'none','Color', temp_c1, 'LineWidth', 1,'CapSize',0);
    pp1=scatter(i,nanmean(ymean1),50,'filled');pp1.MarkerEdgeColor=[1 1 1];pp1.MarkerFaceAlpha=1;pp1.MarkerFaceColor=temp_c1;
    pp1=scatter(i,nanmean(ymean2),50,'filled');pp1.MarkerEdgeColor=[1 1 1];pp1.MarkerFaceAlpha=1;pp1.MarkerFaceColor=temp_c2;box off;
end
xlim([0 7]);ylim([0.5 1])
xticks([]);set(gca,'FontSize',11);ylabel('ILN');set(gca,'TickDir','out');title('SSp-bfd','FontWeight','normal');


rm1=[];rm2=[];rm1=i_m1aam;rm2=c_m1aam;
%calculate and plot 
i_metric=[];c_metric=[];
%calculate indexes across areas
for i=1:length(rm1)
temp_l=[];temp_l2=[];
temp_l=squeeze(rm1(:,i,:));
temp_l2=squeeze(rm2(:,i,:));
[i_metric(:,i,:)] = anatomy_indexcalc(temp_l);
[c_metric(:,i,:)] = anatomy_indexcalc(temp_l2);
end 
temp_metric1=i_metric(:,:,idx);temp_metric2=c_metric(:,:,idx);

nexttile
for i=1:6
    hold on;
    yerr1=[];yerr1=nanstd(nanmean(temp_metric1(:,idx_modules{i}),2))/sqrt(length(nanstd(temp_metric1(:,idx_modules{i}),[],2)));
    yerr2=[];yerr2=nanstd(nanmean(temp_metric2(:,idx_modules{i}),2))/sqrt(length(nanstd(temp_metric2(:,idx_modules{i}),[],2)));
    ymean1=[];ymean1=nanmean(temp_metric1(:,idx_modules{i}),2);
    ymean2=[];ymean2=nanmean(temp_metric2(:,idx_modules{i}),2);
    hold on;h=errorbar([i],[nanmean(ymean2)],[yerr2], 'LineStyle', 'none','Color', temp_c1, 'LineWidth', 0.5,'CapSize',0);
    hold on;h=errorbar([i],[nanmean(ymean1)],[yerr1], 'LineStyle', 'none','Color', temp_c1, 'LineWidth', 1,'CapSize',0);
    pp1=scatter(i,nanmean(ymean1),50,'filled');pp1.MarkerEdgeColor=[1 1 1];pp1.MarkerFaceAlpha=1;pp1.MarkerFaceColor=temp_c1;
    pp1=scatter(i,nanmean(ymean2),50,'filled');pp1.MarkerEdgeColor=[1 1 1];pp1.MarkerFaceAlpha=1;pp1.MarkerFaceColor=temp_c2;box off;
end
xlim([0 7]);ylim([0.5 1])
xticks([1:6]);xticklabels(module_names);set(gca,'FontSize',11);set(gca,'TickDir','out');title('MOp','FontWeight','normal');
%% Layer fraction sorted by L6a using heatmap

%VISp i & c
p1=[];p1=nanmean(i_v1aam,3);
idx_ke=find(sum(v1a_riam>0,2)>1);
kk=[];sort_l6=[];
[sort_l6 kk] = sort(p1(5,idx_ke),'ascend');
p2=[];p2=p1(:,idx_ke);
cortex_abb=[];cortex_abb=cortex_names(:,2);
p3=[];p3=cortex_abb(idx_ke);

fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 500, 1400, 160]);set(gcf,'color','w');tiledlayout("horizontal")
nexttile
imagesc(p2(2:end,kk));
cmap(v1_color,100,5,5);colorbar;box off;xticks([1:45]);
yticks([1:5]);yticklabels({'L2/3','L4','L5','L6a','L6b'});set(gca,'FontSize',10);set(gca,'TickDir','out');caxis([0 1]);
xticklabels(p3(kk));title('VISp ipsi','Fontweight','normal')

nexttile
p1=[];p1=nanmean(c_v1aam,3);
idx_ke=find(sum(v1a_rcam>0,2)>1);
kk=[];sort_l6=[];
[sort_l6 kk] = sort(p1(5,idx_ke),'ascend');
p2=[];p2=p1(:,idx_ke);
cortex_abb=[];cortex_abb=cortex_names(:,2);
p3=[];p3=cortex_abb(idx_ke);
imagesc(p2(2:end,kk));
cmap(v1_color,100,5,5);colorbar;box off;xticks([1:45]);
yticks([1:5]);yticklabels({'L2/3','L4','L5','L6a','L6b'});set(gca,'FontSize',10);set(gca,'TickDir','out');caxis([0 1]);
xticklabels(p3(kk));title('VISp contra','Fontweight','normal')
% 
%SSp i & c
p1=[];p1=nanmean(i_s1aam,3);
idx_ke=find(sum(s1a_riam>0,2)>1);
kk=[];sort_l6=[];
[sort_l6 kk] = sort(p1(5,idx_ke),'ascend');
p2=[];p2=p1(:,idx_ke);
cortex_abb=[];cortex_abb=cortex_names(:,2);
p3=[];p3=cortex_abb(idx_ke);

fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 500, 1400, 160]);set(gcf,'color','w');tiledlayout("horizontal")
nexttile
imagesc(p2(2:end,kk));
cmap(s1_color,100,5,5);colorbar;box off;xticks([1:45]);
yticks([1:5]);yticklabels({'L2/3','L4','L5','L6a','L6b'});set(gca,'FontSize',10);set(gca,'TickDir','out');caxis([0 1]);
xticklabels(p3(kk));title('SSp-bfd ipsi','Fontweight','normal')

nexttile
p1=[];p1=nanmean(c_s1aam,3);
idx_ke=find(sum(s1a_rcam>0,2)>1);
kk=[];sort_l6=[];
[sort_l6 kk] = sort(p1(5,idx_ke),'ascend');
p2=[];p2=p1(:,idx_ke);
cortex_abb=[];cortex_abb=cortex_names(:,2);
p3=[];p3=cortex_abb(idx_ke);
imagesc(p2(2:end,kk));
cmap(s1_color,100,5,5);colorbar;box off;xticks([1:45]);
yticks([1:5]);yticklabels({'L2/3','L4','L5','L6a','L6b'});set(gca,'FontSize',10);set(gca,'TickDir','out');caxis([0 1]);
xticklabels(p3(kk));title('SSp-bfd contra','Fontweight','normal')
 
%MOp i & c
p1=[];p1=nanmean(i_m1aam,3);
idx_ke=find(sum(m1a_riam>0,2)>1);
kk=[];sort_l6=[];
[sort_l6 kk] = sort(p1(5,idx_ke),'ascend');
p2=[];p2=p1(:,idx_ke);
cortex_abb=[];cortex_abb=cortex_names(:,2);
p3=[];p3=cortex_abb(idx_ke);

fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 500, 1400, 160]);set(gcf,'color','w');tiledlayout("horizontal")
nexttile
imagesc(p2(2:end,kk));
cmap(m1_color,100,5,5);colorbar;box off;xticks([1:45]);
yticks([1:5]);yticklabels({'L2/3','L4','L5','L6a','L6b'});set(gca,'FontSize',10);set(gca,'TickDir','out');caxis([0 1]);
xticklabels(p3(kk));title('MOp ipsi','Fontweight','normal')

nexttile
p1=[];p1=nanmean(c_m1aam,3);
idx_ke=find(sum(m1a_rcam>0,2)>1);
kk=[];sort_l6=[];
[sort_l6 kk] = sort(p1(5,idx_ke),'ascend');
p2=[];p2=p1(:,idx_ke);
cortex_abb=[];cortex_abb=cortex_names(:,2);
p3=[];p3=cortex_abb(idx_ke);
imagesc(p2(2:end,kk));
cmap(m1_color,100,5,5);colorbar;box off;xticks([1:45]);
yticks([1:5]);yticklabels({'L2/3','L4','L5','L6a','L6b'});set(gca,'FontSize',10);set(gca,'TickDir','out');caxis([0 1]);
xticklabels(p3(kk));title('MOp contra','Fontweight','normal')
%% L5 vs L6 
idx_modules={frontal_idx lateral_idx somamo_idx visual_idx medial_idx aud_idx};
module_names={'Fron','Lat','SoMo','Vis','Med','Aud'};

rm1=[];rm1=c_m1aam;

fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 500, 400, 200]);set(gcf,'color','w');

for i=1:6
yerr2=[];yerr2=[];
yerr3=[];yerr3=[];
%ymean1=[];ymean1=nanmean(squeeze(rm1(2,idx_modules{i},:)));
ymean2=[];ymean2=nanmean(squeeze(rm1(4,idx_modules{i},:)));
ymean3=[];ymean3=nanmean(squeeze(rm1(5,idx_modules{i},:)));
yerr2=nanstd(ymean2)/sqrt(6)
yerr3=nanstd(ymean3)/sqrt(6)
%ymean_all=[];ymean_all=[nanmean(ymean1) nanmean(ymean2) nanmean(ymean3)];
ymean_all=[];ymean_all=[nanmean(ymean2) nanmean(ymean3)];
b1=bar(i,ymean_all);hold on;b1(1).FaceColor=[1 1 1];b1(2).FaceColor=[0.7 0.7 0.7];%b1(3).FaceColor=[0.55 0.55 0.55];
 hold on;h=errorbar([i-0.16],[nanmean(ymean2)],[yerr2], 'LineStyle', 'none','Color', 'k', 'LineWidth', 1,'CapSize',0);
 hold on;h=errorbar([i+0.16],[nanmean(ymean3)],[yerr3], 'LineStyle', 'none','Color', 'k', 'LineWidth', 1,'CapSize',0);
end
box off;
xticks([1:6]);xticklabels(module_names);set(gca,'FontSize',11);set(gca,'TickDir','out')
legend({'L5','L6'})
title('MOp')
%% 
rm1=[];rm1=c_v1aam;
idx_ke=[];idx_ke=find(sum(v1a_rcam>0,2)>1);
fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 500, 300, 300]);set(gcf,'color','w');
p1=[];p1=nanmean(squeeze(rm1(4,:,:)),2);
p2=[];p2=nanmean(squeeze(rm1(5,:,:)),2);
scatter(p1(idx_ke),p2(idx_ke),40,'filled');
%xlim([0 1]);ylim([0 1])
[r p] = corr(p1(idx_ke),p2(idx_ke),'Type','Spearman','Rows','complete')

%% 

p1=[];p1=nanmean(c_v1aam,3);
idx_ke=find(sum(v1a_rcam>0,2)>1);
kk=[];sort_l6=[];
[sort_l6 kk] = sort(p1(5,idx_ke),'ascend');
p2=[];p2=p1(:,idx_ke);
%% 

rm1=[];rm2=[];
rm1=temp_metric1;rm2=temp_metric2;

temp_rm1=[];temp_rm1=nanmean(rm1);
temp_rm2=[];temp_rm2=nanmean(rm2);
%% 
idx=[];idx=1;
rm1=[];rm2=[];rm1=i_v1aam;rm2=c_v1aam;
temp_c1=[0.8 0.8 0.8];
temp_c2=[0.3 0.3 0.3];
cl_idx=mop_idx;
title_pan='';
%calculate and plot 
i_metric=[];c_metric=[];
%calculate indexes across areas
for i=1:length(rm1)
temp_l=[];temp_l2=[];
temp_l=squeeze(rm1(:,i,:));
temp_l2=squeeze(rm2(:,i,:));
[i_metric(:,i,:)] = anatomy_indexcalc(temp_l);
[c_metric(:,i,:)] = anatomy_indexcalc(temp_l2);
end 
temp_metric1=i_metric(:,:,idx);temp_metric2=c_metric(:,:,idx);
rm1=[];rm2=[];
rm1=temp_metric1;rm2=temp_metric2;

p1=[];[p1]=nanmean(rm2(:,frontal_idx),2)-nanmean(rm1(:,frontal_idx),2);
p2=[];[p2]=nanmean(rm2(:,lateral_idx),2)-nanmean(rm1(:,lateral_idx),2);
p3=[];[p3]=nanmean(rm2(:,somamo_idx),2)-nanmean(rm1(:,somamo_idx),2);
p4=[];[p4]=nanmean(rm2(:,visual_idx),2)-nanmean(rm1(:,visual_idx),2);
p5=[];[p5]=nanmean(rm2(:,medial_idx),2)-nanmean(rm1(:,medial_idx),2);
p6=[];[p6]=nanmean(rm2(:,aud_idx),2)-nanmean(rm1(:,aud_idx),2);

temp_p= [p1 p2 p3 p4 p5 p6];
temp_c=v1_color;
fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 500, 200, 250]);set(gcf,'color','w');tiledlayout("horizontal");
title('VISp','FontWeight','normal')
for i=1:6
    hold on;
    b1=bar(i,nanmean(temp_p(:,i)));hold on;b1.FaceColor=temp_c;
      r=i;rng(i);r1 = r-0.01 + (0.2)*rand(length(temp_p(:,i)),1);
sc1=scatter(r1,temp_p(:,i),12,'k','filled');hold on;sc1.MarkerEdgeColor=[1 1 1];sc1.MarkerEdgeAlpha=0.5;sc1.MarkerFaceColor=[0.2 0.2 0.2];
    hold on;errorbar([i],[nanmean(temp_p(:,i))],[nanstd(temp_p(:,i))/sqrt(length(temp_p(:,i)))]...
    , 'LineStyle', 'none', ... 
        'Color', 'k', 'LineWidth', 1);hold on;

end
xticks([1:6]);hold on;box off;xticklabels({'Fron','Lat','SoMo','Vis','Med','Aud'});ylabel('ILNc-ILNi');%ylim([0 1]);
set(gca,'FontSize',11);set(gca,'TickDir','out');box off;xtickangle(45);
[p,tbl,stats] = anova1([temp_p])
presults = multcompare(stats)
%% 


% %% 
% frac_nr=[];frac_nr=[v1a_rcd];
% tot_nr=[];tot_nr=[v1a_t];
% prct=15;
% title_panel={'VISp ipsi','VISp contra','SSp ipsi','SSp contra','MOp ipsi','MOp contra'}
% [xx rv1_idx] = area_remover(frac_nr, tot_nr,prct,cortex_names(:,2),title_panel)
% %% 
% frac_nr=[];frac_nr=[s1a_rcd];
% tot_nr=[];tot_nr=[s1a_t];
% prct=5;
% title_panel={'SSp ipsi','SSp contra','MOp ipsi','MOp contra'}
% [xx rs1_idx] = area_remover(frac_nr, tot_nr,prct,cortex_names(:,2),title_panel)
% 
% %% 
% frac_nr=[];frac_nr=[m1a_rcd];
% tot_nr=[];tot_nr=[m1a_t];
% prct=10;
% title_panel={'MOp ipsi','MOp contra'}
% [xx rm1_idx] = area_remover(frac_nr, tot_nr,prct,cortex_names(:,2),title_panel)
%% 

%% 



%normalized bu all ipsi+contra
 % v1a_rc=nanmean(nansum(c_v1aa(:,:,:))./(nansum(nansum(c_v1aa(:,:,:)))+nansum(nansum(i_v1aa(:,:,:)))),3);
 % s1a_rc=nanmean(nansum(c_s1aa(:,:,:))./(nansum(nansum(c_s1aa(:,:,:)))+nansum(nansum(i_s1aa(:,:,:)))),3);
 % m1a_rc=nanmean(nansum(c_m1aa(:,:,:))./(nansum(nansum(c_m1aa(:,:,:)))+nansum(nansum(i_m1aa(:,:,:)))),3);
 % v1a_ri=nanmean(nansum(i_v1aa(:,:,:))./(nansum(nansum(c_v1aa(:,:,:)))+nansum(nansum(i_v1aa(:,:,:)))),3);
 % s1a_ri=nanmean(nansum(i_s1aa(:,:,:))./(nansum(nansum(c_s1aa(:,:,:)))+nansum(nansum(i_s1aa(:,:,:)))),3);
 % m1a_ri=nanmean(nansum(i_m1aa(:,:,:))./(nansum(nansum(c_m1aa(:,:,:)))+nansum(nansum(i_m1aa(:,:,:)))),3);



%% 
v1a_rch=squeeze(nansum(c_v1aa(:,[1:30 32:end],:))./nansum(nansum(c_v1aa(:,[1:30 32:end],:))));
s1a_rch=squeeze(nansum(c_s1aa(:,[1:17 19:end],:))./nansum(nansum(c_s1aa(:,[1:17 19:end],:))));
m1a_rch=squeeze(nansum(c_m1aa(:,[1:24 26:end],:))./nansum(nansum(c_m1aa(:,[1:24 26:end],:))));
figure;
plot(sort(nanmean(v1a_rch,2),'descend'),'-o');hold on;
plot(sort(nanmean(s1a_rch,2),'descend'),'-o');hold on;
plot(sort(nanmean(m1a_rch,2),'descend'),'-o');hold on;


%% 
%frac_nr=[v1a_ri' s1a_ri' m1a_ri' v1a_rc' s1a_rc' m1a_rc'];
%tot_nr=[v1a_ti' s1a_ti' m1a_ti' v1a_tc' s1a_tc' m1a_tc'];
frac_nr=[];frac_nr=[v1a_ri'  v1a_rc'];
tot_nr=[];tot_nr=[v1a_ti' v1a_tc'];
prct=20;
title_panel={'VISp ipsi','VISp contra'}
[visp_kidx] = cell_keeper(frac_nr, tot_nr,prct,cortex_names(:,2),title_panel);

%% 
frac_nr=[];frac_nr=[s1a_ri'  s1a_rc'];
tot_nr=[];tot_nr=[s1a_ti' s1a_tc'];
prct=10;
title_panel={'SSp ipsi','SSp contra'}
[ssp_kidx] = cell_keeper(frac_nr, tot_nr,prct,cortex_names(:,2),title_panel)
%% 
frac_nr=[];frac_nr=[m1a_ri'  m1a_rc'];
tot_nr=[];tot_nr=[m1a_ti' m1a_tc'];
prct=30;
title_panel={'MOp ipsi','MOp contra'}
[mp_kidx] = cell_keeper(frac_nr, tot_nr,prct,cortex_names(:,2),title_panel)
%% 

%% 

bwi=0.1;
 idxv1=visp_kidx(:,2);
 idxs1=ssp_kidx(:,2);
 idxm1=mp_kidx(:,2);
p1=[];p1=(v1a_tca(idxv1,:)-v1a_tia(idxv1,:))./(v1a_tca(idxv1,:)+v1a_tia(idxv1,:));
p2=[];p2=(s1a_tca(idxs1,:)-s1a_tia(idxs1,:))./(s1a_tca(idxs1,:)+s1a_tia(idxs1,:));
p3=[];p3=(m1a_tca(idxm1,:)-m1a_tia(idxm1,:))./(m1a_tca(idxm1,:)+m1a_tia(idxm1,:));
fig7=figure;set(gcf,'color','w');set(fig7, 'Position', [1200, 400 ,350, 350]);
h1=histogram(p1,20,'Normalization','probability');hold on;h1.EdgeColor=v1_color;h1.FaceColor='w';h1.FaceAlpha=0.2;h1.LineWidth=1.5;h1.BinWidth = bwi;
h3=histogram(p2,20,'Normalization','probability');hold on;h3.EdgeColor=s1_color;h3.FaceColor='w';h3.FaceAlpha=0.2;h3.LineWidth=1.5;box off;h3.BinWidth = bwi;
h2=histogram(p3,20,'Normalization','probability');hold on;h2.EdgeColor=m1_color;h2.FaceColor='w';h2.FaceAlpha=0.2;h2.LineWidth=1.5;h2.BinWidth = bwi;
xlabel('contra-ipsi index');set(gca,'FontSize',12);ylabel('Fraction');title([]);
g1=[];g2=[];g3=[];par=[];par=[p1(:)' p2(:)' p3(:)'];g1=[];g1=ones(1,length(p1(:)));g2=[];g2=ones(1,length(p2(:)))*2;g3=[];g3=ones(1,length(p3(:)))*3;gro=[];gro=[g1 g2 g3]';
[p,tbl,statsout] = kruskalwallis(par,gro);oo = multcompare(statsout);
%% 
nanmean(nansum(p1>0)/length(v1a_tia))

%% 
sum(p1_t>0)/length(p1_t)
squeeze(nansum(c_v1aa(:,:,:))./nansum(nansum(c_v1aa(:,:,:))))



