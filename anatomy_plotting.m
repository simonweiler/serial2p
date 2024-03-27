%This scipt reads out the saved mat data structure with all animals and
%injection inlcuded and plots figures/ performs analysis 
%% Dependendencies 
%uipickfiles: https://uk.mathworks.com/matlabcentral/fileexchange/10867-uipickfiles-uigetfile-on-steroids
%Matlab structure

%% Folder where mat structure is (CHANGE ACCORDINGLY)
str   = 'C:\Users\simonw\S1-V1 interaction Dropbox\Simon_Manuel\Manuscripts\L6_dominance\Injections\output_structure';
% load structure (mat file) using uipickfiles.m 1)
folder_list = uipickfiles('FilterSpec',str);load(char(folder_list));

%folder location where to save figure panels
save_folder='C:\Users\simonw\S1-V1 interaction Dropbox\Simon_Manuel\Manuscripts\L6_dominance\Vectorgraphics';
% folder where to save values for heatmaps for flatmap for bilateral symmetry
heatmap_folder='C:\Users\simonw\S1-V1 interaction Dropbox\Simon_Manuel\Manuscripts\L6_dominance\Heatmaps_bilateral_symmetry';
%% Perform initial calculations and extraction for all mice together and/or specific for injection 
%arial as default
set(0, 'DefaultAxesFontName', 'Arial'); 
% Color scheme for the plots
v1_color=[173 7 227]./256;%magenta V1
s1_color=[5 190 120]./256;%green S1
m1_color=[87 87 244]./256;%blue M1
% Load csv with cortex names and assign 6 larger group areas: Frontal, Lateral, Somatomotor, Visual, Medial, Auditory
cortex_info = readtable([str '\names_abbreviations.csv']);
%read out nr, abbreviation and long cortex name of ABI
cortex_names=table2cell(cortex_info);
%change here index of categories
frontal_idx=[1:8];lateral_idx=[9:16 44:45];somamo_idx=[17:26];visual_idx=[27:36];medial_idx=[37:39];aud_idx=[40:43];
%injected brain areas VISp, SSp, MOp, (AUDp)
visp_idx=[31];ssp_idx=[18];mop_idx=[25];audp_idx=[41];
%L4 missing areas
l4miss_idx=[1:11 15:16 25:26 37 39 44 45];
% Find the cells in a given area in their given layer RETRO
% select all retro injection regardeless of their injection type using
% cell_selecter.m 2)
temp1=[];temp1=cell_selecter(data,'virustype',1);
retro=[];retro=find(temp1==1);
%call function anatomy_cellnr.m 3)
[i_animal c_animal i_areas_animal c_areas_animal i_areas_animalf c_areas_animalf  i_areas_animalcd c_areas_animalcd] = anatomy_cellnr(data,retro,cortex_names);
% Extract fraction per injeciton area M1 / S1 / V1
temp2=[];temp2=cell_selecter(data,'virustype',1,'area',1);
temp3=[];temp3=cell_selecter(data,'virustype',1,'area',2);
temp4=[];temp4=cell_selecter(data,'virustype',1,'area',3);
v1r=[];v1r=find(temp2==1);s1r=[];s1r=find(temp3==1);m1r=[];m1r=find(temp4==1);
%call function anatomy_cellnr.m
%V1, S1 and M1
[i_v1a c_v1a i_v1aa c_v1aa i_v1aaf c_v1aaf i_v1aacd c_v1aacd] = anatomy_cellnr(data,v1r,cortex_names);
idx_nanv1=find(sum(isnan(i_v1aa(3,:,:)),3)>0 & sum(isnan(i_v1aa(3,:,:)),3)<size(i_v1aa,3));
i_v1aaf(3,idx_nanv1,:)=NaN;i_v1aa(3,idx_nanv1,:)=NaN;i_v1aacd(3,idx_nanv1,:)=NaN;
c_v1aaf(3,idx_nanv1,:)=NaN;c_v1aa(3,idx_nanv1,:)=NaN;c_v1aacd(3,idx_nanv1,:)=NaN;
[i_s1a c_s1a i_s1aa c_s1aa i_s1aaf c_s1aaf i_s1aacd c_s1aacd] = anatomy_cellnr(data,s1r,cortex_names);
idx_nans1=find(sum(isnan(i_s1aa(3,:,:)),3)>0 & sum(isnan(i_s1aa(3,:,:)),3)<size(i_s1aa,3));
i_s1aaf(3,idx_nans1,:)=NaN;i_s1aa(3,idx_nans1,:)=NaN;i_s1aacd(3,idx_nans1,:)=NaN;
c_s1aaf(3,idx_nans1,:)=NaN;c_s1aa(3,idx_nans1,:)=NaN;c_v1aacd(3,idx_nans1,:)=NaN;
[i_m1a c_m1a i_m1aa c_m1aa i_m1aaf c_m1aaf i_m1aacd c_m1aacd] = anatomy_cellnr(data,m1r,cortex_names); 
idx_nanm1=find(sum(isnan(i_m1aa(3,:,:)),3)>0 & sum(isnan(i_m1aa(3,:,:)),3)<size(i_m1aa,3));
i_m1aaf(3,idx_nanm1,:)=NaN;i_m1aa(3,idx_nanm1,:)=NaN;i_m1aacd(3,idx_nanm1,:)=NaN;
c_m1aaf(3,idx_nanm1,:)=NaN;c_m1aa(3,idx_nanm1,:)=NaN;c_m1aacd(3,idx_nanm1,:)=NaN;

% calculate indexes using anatomy_indexcalc.m 4) : ILN, L6ab, L6a, hindex ALL INJECTIONS TOGETHER
    %1= ILN
    %2= L6ab dominance index
    %3= L6a dominacne index
    %4= h-index
    %5= L6 pure index 
%overal indexes
[i_index] = anatomy_indexcalc(i_animal);
[c_index] = anatomy_indexcalc(c_animal);
%indexes for specific injection area: V1, S1, M1
[iv1_index] = anatomy_indexcalc(i_v1a);
[cv1_index] = anatomy_indexcalc(c_v1a);
[is1_index] = anatomy_indexcalc(i_s1a);
[cs1_index] = anatomy_indexcalc(c_s1a);
[im1_index] = anatomy_indexcalc(i_m1a);
[cm1_index] = anatomy_indexcalc(c_m1a);
%% criteria when selecting cellcounts per 45 brain areas, and mice 
%Decided with TWM on the 11.3.2024: 
%when dividing into 45 areas we exclude an area with 10 or less cells and
%this need to be present at least in 3 animals
%cell threshold 
thr_count=10;
%Animal threshold 
an_count=3;
%first output: absolute numbers per layer per area per animal
%second output: relative numbers per layer per area per animal
[i_v1aam i_v1aafm]=zero_area_layers(i_v1aa,i_v1aaf,thr_count,an_count);
[c_v1aam c_v1aafm]=zero_area_layers(c_v1aa,c_v1aaf,thr_count,an_count);
[i_s1aam i_s1aafm]=zero_area_layers(i_s1aa,i_s1aaf,thr_count,an_count);
[c_s1aam c_s1aafm]=zero_area_layers(c_s1aa,c_s1aaf,thr_count,an_count);
[i_m1aam i_m1aafm]=zero_area_layers(i_m1aa,i_m1aaf,thr_count,an_count);
[c_m1aam c_m1aafm]=zero_area_layers(c_m1aa,c_m1aaf,thr_count,an_count);

%absolute across animals per injection after criteria 
v1a_tm=[nanmean(nansum(i_v1aam(:,:,:)),3)' nanmean(nansum(c_v1aam(:,:,:)),3)'];
s1a_tm=[nanmean(nansum(i_s1aam(:,:,:)),3)' nanmean(nansum(c_s1aam(:,:,:)),3)'];
m1a_tm=[nanmean(nansum(i_m1aam(:,:,:)),3)' nanmean(nansum(c_m1aam(:,:,:)),3)'];

%absolute numbers per animal for ipsi and contra 
v1a_tcam=squeeze(nansum(c_v1aam(:,:,:)));
v1a_tiam=squeeze(nansum(i_v1aam(:,:,:)));
s1a_tcam=squeeze(nansum(c_s1aam(:,:,:)));
s1a_tiam=squeeze(nansum(i_s1aam(:,:,:)));
m1a_tcam=squeeze(nansum(c_m1aam(:,:,:)));
m1a_tiam=squeeze(nansum(i_m1aam(:,:,:)));

%in case ine usese the homotopi contralteral set to zero 
%set homotopic to NaN
v1a_tcam2=v1a_tcam;
s1a_tcam2=s1a_tcam;
m1a_tcam2=m1a_tcam;
%absolute with NaN
v1a_tcam2(visp_idx,:)=ones(1,size(v1a_tcam2,2))*NaN;
s1a_tcam2(ssp_idx,:)=ones(1,size(s1a_tcam2,2))*NaN;
m1a_tcam2(mop_idx,:)=ones(1,size(m1a_tcam2,2))*NaN;
%relatives setting contra homotopic to NaN
v1a_rcam2=v1a_tcam2./nansum(v1a_tcam2);
s1a_rcam2=s1a_tcam2./nansum(s1a_tcam2);
m1a_rcam2=m1a_tcam2./nansum(m1a_tcam2);

% relatives numbers per animal for ipsi and contra NOT setting homotopic
% contalteral to NAN
v1a_rcam=squeeze(nansum(c_v1aam(:,:,:))./nansum(nansum(c_v1aam(:,:,:))));
v1a_riam=squeeze(nansum(i_v1aam(:,:,:))./nansum(nansum(i_v1aam(:,:,:))));
s1a_rcam=squeeze(nansum(c_s1aam(:,:,:))./nansum(nansum(c_s1aam(:,:,:))));
s1a_riam=squeeze(nansum(i_s1aam(:,:,:))./nansum(nansum(i_s1aam(:,:,:))));
m1a_rcam=squeeze(nansum(c_m1aam(:,:,:))./nansum(nansum(c_m1aam(:,:,:))));
m1a_riam=squeeze(nansum(i_m1aam(:,:,:))./nansum(nansum(i_m1aam(:,:,:))));

%IPSI: 
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
%% Plot what is left after criteria and thrsholding 
ipsi_all={};ipsi_all={v1a_tiam s1a_tiam m1a_tiam};
contra_all={};contra_all={v1a_tcam2 s1a_tcam2 m1a_tcam2};

fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 0, 1500, 800]);set(gcf,'color','w');
tiledlayout(2,1);nexttile
keeper = ([ipsi_all{1} ipsi_all{2} ipsi_all{3}]>0);
imagesc(keeper);colormap gray;
for i=1:size([ipsi_all{1} ipsi_all{2} ipsi_all{3}],2)
line([i+0.5 i+0.5],[1 45],'Color','k','LineStyle','--');
end
hold on;line([6.5 6.5],[1 45],'Color','r');
hold on;line([12.5 12.5],[1 45],'Color','r');
nexttile
keeper = ([contra_all{1} contra_all{2} contra_all{3}]>0);
colormap gray;imagesc(keeper);
for i=1:size([ipsi_all{1} ipsi_all{2} ipsi_all{3}],2)
line([i+0.5 i+0.5],[1 45],'Color','k','LineStyle','--');
end
hold on;line([6.5 6.5],[1 45],'Color','r');
hold on;line([12.5 12.5],[1 45],'Color','r');

%% add total numbers in a strcuture/csv for MT-> for heat maps in flat map
cd(heatmap_folder)
ipsi_contra_cortex_homo_included =table(cortex_names(:,2),v1a_tm(:,1),v1a_tm(:,2),s1a_tm(:,1),s1a_tm(:,2),m1a_tm(:,1),m1a_tm(:,2),...
    'variablenames',{'Area','V1ipsi','V1contra','S1ipsi','S1contra','M1ipsi','M1contra'});
writetable(ipsi_contra_cortex_homo_included,'ipsi_contra_cortex_homo_included.csv')
%% add total numbers in a strcuture/csv for MT-> for heat maps in flat map
v1a_tm2=[];s1a_tm2=[];m1a_tm2=[];
v1a_tm2=v1a_tm;s1a_tm2=s1a_tm;m1a_tm2=m1a_tm;
v1a_tm2(visp_idx,:)=0;s1a_tm2(ssp_idx,:)=0;m1a_tm2(mop_idx,:)=0;
ipsi_contra_cortex_homo_excluded =table(cortex_names(:,2),v1a_tm2(:,1),v1a_tm2(:,2),s1a_tm2(:,1),s1a_tm2(:,2),m1a_tm2(:,1),m1a_tm2(:,2),...
    'variablenames',{'Area','V1ipsi','V1contra','S1ipsi','S1contra','M1ipsi','M1contra'});
cd(heatmap_folder)
writetable(ipsi_contra_cortex_homo_excluded,'ipsi_contra_cortex_homo_excluded.csv')
%% FIGURES START HERE

%% Figure 1
%% Plot all cell numbers ipsi contra colour-coded per injection type
%all super imposed
temp_color=[v1_color ;s1_color; m1_color];module_names={'VISp','SSp-bfd','MOp'};
all_dati={};all_dati={squeeze(nansum(nansum(i_v1aam)))./(squeeze(nansum(nansum(i_v1aam)))+squeeze(nansum(nansum(c_v1aam))))...
     squeeze(nansum(nansum(i_s1aam)))./(squeeze(nansum(nansum(i_s1aam)))+squeeze(nansum(nansum(c_s1aam))))...
     squeeze(nansum(nansum(i_m1aam)))./(squeeze(nansum(nansum(i_m1aam)))+squeeze(nansum(nansum(c_m1aam))))};
all_datc={};all_datc={squeeze(nansum(nansum(c_v1aam)))./(squeeze(nansum(nansum(i_v1aam)))+squeeze(nansum(nansum(c_v1aam))))...
     squeeze(nansum(nansum(c_s1aam)))./(squeeze(nansum(nansum(i_s1aam)))+squeeze(nansum(nansum(c_s1aam))))...
     squeeze(nansum(nansum(c_m1aam)))./(squeeze(nansum(nansum(i_m1aam)))+squeeze(nansum(nansum(c_m1aam))))};
%plot
fig3= figure;set(fig3, 'Name', 'Paired comp');set(fig3, 'Position', [200, 300, 200, 250]);set(gcf,'color','w');
for j=1:3
dat=[];dat=[all_dati{:,j} all_datc{:,j}];
%lines between paired data points
for i=1:length(data)
     pl=plot([1,2],[dat(:,1),dat(:,2)],'color',[0.5 0.5 0.5]);    
end
%single data points
hold on;pS=plotSpread([dat(:,1),dat(:,2)],'categoryIdx',[ones(1,length(dat(:,1)))' ones(1,length(dat(:,2)))'*2],...
    'categoryMarkers',{'o','o'},'categoryColors',{temp_color(j,:), temp_color(j,:)});hold on;
%mean +- SEM
hold on;er1=errorbar([0.2+j*0.14],nanmean(dat(:,1)),nanstd(dat(:,1),[],1)/sqrt(length(dat(:,1))),'ok','MarkerFaceColor',temp_color(j,:),'Markersize',8);
hold on;er2=errorbar([2.2+j*0.14],nanmean(dat(:,2)),nanstd(dat(:,2),[],1)/sqrt(length(dat(:,1))),'ok','MarkerFaceColor',temp_color(j,:),'Markersize',8);
%legend colour
text(2.2,1.1-0.07*j,module_names{j},'Color',temp_color(j,:));
%stats
[u p1]=ttest(dat(:,1),dat(:,2))
end
xticklabels({'ipsi','contra'});ylabel('Fraction');hold on;title([]);set(gca,'FontSize',11);set(gca,'TickDir','out');box off;
hold on;text(1.25,1,['***'],'FontSize',18);offsetAxes;h = gca;h.XAxis.Visible = 'off'
 t1=text(0.5,-0.1,'ipsi','FontSize',11,'Color',[0.7 0.7 0.7]);set(t1,'Rotation',45);t1=text(1.5,-0.1,'contra','FontSize',11,'Color',[0.3 0.3 0.3]);set(t1,'Rotation',45); ax=gca;ax.LineWidth=1;
% save
cd(save_folder);saveas(gcf, 'ipsi_contra_fraction.pdf');
%% Bar plots showing area numbers ipsi/ contra 
%here contra homotopic needs to be removed so use v1a_tcam2 etc
temp_color=[v1_color ;s1_color; m1_color];
%ipsi_total numbers area
p1=[];p1=sum(v1a_tiam>0);p2=[];p2=sum(s1a_tiam>0);p3=[];p3=sum(m1a_tiam>0);
%contra_total numbers area
p4=[];p4=sum(v1a_tcam2>0);p5=[];p5=sum(s1a_tcam2>0);p6=[];p6=sum(m1a_tcam2>0);
temp_p= {p1' p2' p3'};temp_p2= {p4' p5' p6'};
fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 500, 300, 270]);set(gcf,'color','w');
for i=1:3
    hold on;
     b2=bar([i],[nanmean(temp_p{:,i}) nanmean(temp_p2{:,i})],0.7);
     hold on;b2(1).FaceColor=[0.8 0.8 0.8];b2(1).EdgeColor=temp_color(i,:);b2(1).LineWidth=1.4;set(b2,'ShowBaseLine','off')
     hold on;b2(2).FaceColor=[0.3 0.3 0.3];b2(2).EdgeColor=temp_color(i,:);b2(2).LineWidth=1.4;set(b2,'ShowBaseLine','off')
     hold on;errorbar([i-0.15],[nanmean(temp_p{:,i})],[nanstd(temp_p{:,i})/sqrt(length(temp_p{:,i}))]...
      , 'LineStyle', 'none', ... 
          'Color', 'k', 'LineWidth', 1.2);hold on;
      hold on;errorbar([i+0.15],[nanmean(temp_p2{:,i})],[nanstd(temp_p2{:,i})/sqrt(length(temp_p2{:,i}))]...
      , 'LineStyle', 'none', ... 
          'Color', 'k', 'LineWidth', 1.2);hold on;
      disp(['the range ipsi is ' num2str(min(temp_p{:,i})) ' - ' num2str(max(temp_p{:,i}))])
      disp(['the range contra is ' num2str(min(temp_p2{:,i})) ' - ' num2str(max(temp_p2{:,i}))])
end
t1=text(0.85,18,[num2str(min(p1)) ' - ' num2str(max(p1))],'FontSize',10,'Color','k');set(t1,'Rotation',90);
t1=text(1.15,18,[num2str(min(p4)) ' - ' num2str(max(p4))],'FontSize',10,'Color','w');set(t1,'Rotation',90);t1=text(0.82,45,'ipsi','FontSize',11);set(t1,'Rotation',90);
t1=text(1.15,45,'contra','FontSize',11);set(t1,'Rotation',90);t1=text(1.85,18,[num2str(min(p2)) ' - ' num2str(max(p2))],'FontSize',10,'Color','k');set(t1,'Rotation',90);
t1=text(2.15,18,[num2str(min(p5)) ' - ' num2str(max(p5))],'FontSize',10,'Color','w');set(t1,'Rotation',90);t1=text(1.82,45,'ipsi','FontSize',11);set(t1,'Rotation',90);
t1=text(2.15,45,'contra','FontSize',11);set(t1,'Rotation',90);t1=text(2.85,9,[num2str(min(p3)) ' - ' num2str(max(p3))],'FontSize',10,'Color','k');set(t1,'Rotation',90);
t1=text(3.15,9,[num2str(min(p6)) ' - ' num2str(max(p6))],'FontSize',10,'Color','w');set(t1,'Rotation',90);
t1=text(2.82,45,'ipsi','FontSize',11);set(t1,'Rotation',90);t1=text(3.15,45,'contra','FontSize',11);set(t1,'Rotation',90);
set(gca,'FontSize',11);set(gca,'TickDir','out');box off;xtickangle(45);ylim([-1 52]);axis off

%TEST
x1=[];x2=[];x3=[];x1 = temp_p{:,1};x2 = temp_p{:,2};x3 = temp_p{:,3};
dat_test = [x1' x2' x3']; %// Create row vector with your data
group = {'G1','G1','G1','G1','G1','G1','G2','G2','G2','G2','G2','G2','G3','G3','G3','G3','G3','G3'}; %// set the groups according to the data above
[p,tbl,stats]  = anova1(dat_test, group) %// Use the 'off' option to prevent the table/box plot from showing up.
presults = multcompare(stats)

x1=[];x2=[];x3=[];x1 = temp_p2{:,1};x2 = temp_p2{:,2};x3 = temp_p2{:,3};
dat_test = [x1' x2' x3']; %// Create row vector with your data
[p,tbl,stats]  = anova1(dat_test, group)
presults = multcompare(stats)

%ttest ipsi contra
disp('VISp test ipsi vs contra');[u p11]=ttest(temp_p{:,1},temp_p2{:,1})
nanmean(temp_p{:,1})
nanstd([temp_p{:,1}]/sqrt(length([temp_p{:,1}])))
nanmean(temp_p2{:,1})
nanstd([temp_p2{:,1}]/sqrt(length([temp_p2{:,1}])))
disp('SSp test ipsi vs contra');[u p11]=ttest(temp_p{:,2},temp_p2{:,2})
nanmean(temp_p{:,2})
nanstd([temp_p{:,2}]/sqrt(length([temp_p{:,2}])))
nanmean(temp_p2{:,2})
nanstd([temp_p2{:,2}]/sqrt(length([temp_p2{:,2}])))
disp('MOp test ipsi vs contra');[u p11]=ttest(temp_p{:,3},temp_p2{:,3})
nanmean(temp_p{:,3})
nanstd([temp_p{:,3}]/sqrt(length([temp_p{:,3}])))
nanmean(temp_p2{:,3})
nanstd([temp_p2{:,3}]/sqrt(length([temp_p2{:,3}])))
%Sensory vs motor
%ipsi
disp('ipsi sensory vs motor');[u p11]=ttest2([temp_p{:,1} ;temp_p{:,2}],temp_p{:,3})
nanmean([temp_p{:,1} ;temp_p{:,2}])
nanstd([temp_p{:,1} ;temp_p{:,2}]/sqrt(length([temp_p{:,1} ;temp_p{:,2}])))
nanmean(temp_p{:,3})
nanstd([temp_p{:,3}]/sqrt(length([temp_p{:,3}])))
%contra
disp('contra sensory vs motor');[u p11]=ttest2([temp_p2{:,1} ;temp_p2{:,2}],temp_p2{:,3})
nanmean([temp_p2{:,1} ;temp_p2{:,2}])
nanstd([temp_p2{:,1} ;temp_p2{:,2}]/sqrt(length([temp_p2{:,1} ;temp_p2{:,2}])))
nanmean(temp_p2{:,3})
nanstd([temp_p2{:,3}]/sqrt(length([temp_p2{:,3}])))

%% save
cd(save_folder);saveas(gcf, 'ipsi_contra_are_numbercounts.pdf');
%% 
%% Fraction homo vs fraction hetero per injection target 
temp_color=[v1_color ;s1_color; m1_color];
%contra fraction homo
p1=[];p1=v1a_rcam(visp_idx,:);p2=[];p2=s1a_rcam(ssp_idx,:);p3=[];p3=m1a_rcam(mop_idx,:);
%contra fraction hetero
p4=[];p4=sum(v1a_rcam([1:visp_idx-1 visp_idx+1:45] ,:));p5=[];p5=sum(s1a_rcam([1:ssp_idx-1 ssp_idx+1:45] ,:));p6=[];p6=sum(m1a_rcam([1:mop_idx-1 mop_idx+1:45] ,:));
temp_p= [p1' p2' p3'];temp_p2= [p4' p5' p6'];
%plot
fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 500, 300, 270]);set(gcf,'color','w');
for i=1:3
    hold on;
    b1=bar([i*2],[nanmean(temp_p2(:,i)) nanmean(temp_p(:,i))],'stacked');hold on;b1(2).FaceColor=[1 1 1];b1(1).FaceColor=temp_color(i,:);
    b1(1).FaceAlpha=[0.6];b1(2).EdgeColor=temp_color(i,:);b1(1).LineWidth=1;b1(2).LineWidth=1;b1(1).EdgeColor=temp_color(i,:);
     hold on;errorbar([i*2],[nanmean(temp_p2(:,i)) nanmean(temp_p(:,i))+nanmean(temp_p2(:,i))],[nanstd(temp_p2(:,i))/sqrt(length(temp_p2(:,i))) nanstd(temp_p(:,i))/sqrt(length(temp_p(:,i)))]...
     , 'LineStyle', 'none', ... 
         'Color', 'k', 'LineWidth', 1.2);hold on;
     set(b1,'ShowBaseLine','off')
end
xticks([2:2:6]);hold on;box off;xticklabels({[]});%ylabel('Nr of areas');
set(gca,'FontSize',11);set(gca,'TickDir','out');box off;xtickangle(45);axis off
text(1.7,0.9,[num2str(round(nanmean(p1*100)),2) '%'],'FontSize',9,'Color','k');
text(1.7,0.45,[num2str(round(nanmean(p4*100)),2) '%'],'FontSize',9,'Color','k');
text(3.7,0.9,[num2str(round(nanmean(p2*100)),2) '%'],'FontSize',9,'Color','k');
text(3.7,0.45,[num2str(round(nanmean(p5*100)),2) '%'],'FontSize',9,'Color','k');
text(5.7,0.75,[num2str(round(nanmean(p3*100)),2) '%'],'FontSize',9,'Color','k');
text(5.7,0.27,[num2str(round(nanmean(p6*100)),2) '%'],'FontSize',9,'Color','k');
 t1=text(1,0.23,'He-topic','FontSize',11);set(t1,'Rotation',90);
 t1=text(1,0.75,'Ho-topic','FontSize',11);set(t1,'Rotation',90);
  t1=text(3,0.23,'He-topic','FontSize',11);set(t1,'Rotation',90);
 t1=text(3,0.75,'Ho-topic','FontSize',11);set(t1,'Rotation',90);
  t1=text(5,0.13,'He-topic','FontSize',11);set(t1,'Rotation',90);
 t1=text(5,0.55,'Ho-topic','FontSize',11);set(t1,'Rotation',90);
 %% save
cd(save_folder);saveas(gcf, 'homo_hetero_fraction.pdf');
%% Correlation between hemisphere areas of fraction per area avergae across injection areas
%Calculate R2 and slope of correlation 
%V1
fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 500, 700, 250]);set(gcf,'color','w');tiledlayout("horizontal")
nexttile
%single example 
for i=2:2%size(v1a_riam,2)
idx = find(~isnan(v1a_riam(:,i)));
mdl = fitlm(v1a_riam(idx,i),v1a_rcam(idx,i));
coefs = polyfit(v1a_riam(idx,i)', v1a_rcam(idx,i)', 1);
slope_v1(i)=coefs(1);
r2_sqv1(i)=mdl.Rsquared.Ordinary;
hold on;ph=plot(mdl,'Color',[0.8 0.8 0.8]);hold on;
ph(1).Marker='o';
ph(1).MarkerFaceColor=[0.7 0.7 0.7];
ph(1).MarkerEdgeColor=[1 1 1];
ph(2).Color=[0.7 0.7 0.7];
ph(3).Color=[0.7 0.7 0.7];
ph(4).Color=[0.7 0.7 0.7];
end
hold on;title([]);legend('');xlabel('');ylabel('');
set(gca,'FontSize',11);set(gca,'TickDir','out');legend box off;axis square;offsetAxes;

%average linear r2 line 
hold on;
idx=[];tma1=[];tma2=[];
tma1=nanmean(v1a_riam,2);
tma2=nanmean(v1a_rcam,2);
idx = find(~isnan(tma1));
mdl = fitlm(tma1(idx),tma2(idx));
coefs = polyfit(tma1(idx)',tma2(idx)', 1);
% slope_v1=coefs(1);
% r2_sqv1=mdl.Rsquared.Ordinary;
h = plot(mdl,'Color',v1_color);
delete(h(1));
h(2).Color=v1_color;
h(2).LineWidth=2;
h(3).Color=v1_color;
h(4).Color=v1_color;
hold on;title([]);legend('');xlabel('');ylabel('');
set(gca,'FontSize',11);set(gca,'TickDir','out');legend box off;axis square;offsetAxes;

%S1
nexttile
%single example 
for i=2:2%size(s1a_riam,2)
idx = find(~isnan(s1a_riam(:,i)));
mdl = fitlm(s1a_riam(idx,i),s1a_rcam(idx,i));
coefs = polyfit(s1a_riam(idx,i)', s1a_rcam(idx,i)', 1);
slope_s1(i)=coefs(1);
r2_sqs1(i)=mdl.Rsquared.Ordinary;
hold on;ph=plot(mdl,'Color',[0.8 0.8 0.8]);hold on;
ph(1).Marker='o';
ph(1).MarkerFaceColor=[0.7 0.7 0.7];
ph(1).MarkerEdgeColor=[1 1 1];
ph(2).Color=[0.7 0.7 0.7];
ph(3).Color=[0.7 0.7 0.7];
ph(4).Color=[0.7 0.7 0.7];
end
hold on;title([]);legend('');xlabel('');ylabel('');
set(gca,'FontSize',11);set(gca,'TickDir','out');legend box off;axis square;offsetAxes

%average linear r2 line 
hold on;
idx=[];tma1=[];tma2=[];
tma1=nanmean(s1a_riam,2);
tma2=nanmean(s1a_rcam,2);
idx = find(~isnan(tma1));
mdl = fitlm(tma1(idx),tma2(idx));
coefs = polyfit(tma1(idx)',tma2(idx)', 1);
% slope_v1=coefs(1);
% r2_sqv1=mdl.Rsquared.Ordinary;
h = plot(mdl,'Color',s1_color);
delete(h(1));
h(2).Color=s1_color;
h(2).LineWidth=2;
h(3).Color=s1_color;
h(4).Color=s1_color;
hold on;title([]);legend('');xlabel('Ipsi fraction');ylabel('Contra fraction');
set(gca,'FontSize',11);set(gca,'TickDir','out');legend box off;axis square;offsetAxes;


%M1
nexttile
%single example
for i=5:5%size(m1a_riam,2)
idx = find(~isnan(m1a_riam(:,i)));
mdl = fitlm(m1a_riam(idx,i),m1a_rcam(idx,i));
coefs = polyfit(m1a_riam(idx,i)', m1a_rcam(idx,i)', 1);
slope_m1(i)=coefs(1);
r2_sqm1(i)=mdl.Rsquared.Ordinary;
hold on;ph=plot(mdl,'Color',[0.8 0.8 0.8]);hold on;
ph(1).Marker='o';
ph(1).MarkerFaceColor=[0.7 0.7 0.7];
ph(1).MarkerEdgeColor=[1 1 1];
ph(2).Color=[0.7 0.7 0.7];
ph(3).Color=[0.7 0.7 0.7];
ph(4).Color=[0.7 0.7 0.7];
end
hold on;title([]);legend('');xlabel('');ylabel('');
set(gca,'FontSize',11);set(gca,'TickDir','out');legend box off;axis square;offsetAxes

%average linear r2 line 
hold on;
idx=[];tma1=[];tma2=[];
tma1=nanmean(m1a_riam,2);
tma2=nanmean(m1a_rcam,2);
idx = find(~isnan(tma1));
mdl = fitlm(tma1(idx),tma2(idx));
coefs = polyfit(tma1(idx)',tma2(idx)', 1);
% slope_v1=coefs(1);
% r2_sqv1=mdl.Rsquared.Ordinary;
h = plot(mdl,'Color',m1_color);
delete(h(1));
h(2).Color=m1_color;
h(2).LineWidth=2;
h(3).Color=m1_color;
h(4).Color=m1_color;
hold on;title([]);legend('');xlabel('');ylabel('');
set(gca,'FontSize',11);set(gca,'TickDir','out');legend box off;axis square;offsetAxes;
%% save
cd(save_folder);saveas(gcf, 'correlation_ipsi_contra_lines.pdf');
%% 


%% Bilateral symmetry test an logarthmic plot
v1_hemi=(v1a_rcam2-v1a_riam)./(v1a_rcam2+v1a_riam);
s1_hemi=(s1a_rcam2-s1a_riam)./(s1a_rcam2+s1a_riam);
m1_hemi=(m1a_rcam2-m1a_riam)./(m1a_rcam2+m1a_riam);
% Plot 
temp_color=[v1_color ;s1_color; m1_color];
hemi_all={v1_hemi s1_hemi m1_hemi};
temp_data_i={v1a_riam s1a_riam m1a_riam};
temp_data_c={v1a_rcam2 s1a_rcam2 m1a_rcam};
panel_title={'VISp','SSp-bfd','MOp'};
crt_all1={};
fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 500, 1000, 500]);set(gcf,'color','w');tiledlayout("horizontal")
for i=1:3
nexttile
[p_bin p_bin_c] = anatomy_testarea(hemi_all{i});
disp(['bilateral' num2str(i)])
sum(p_bin==0 & p_bin_c==0)
rm1=[];rm2=[];
rm1=temp_data_i{i};rm2=temp_data_c{i};
%1)
errx_l=[];errx_l=[nanstd(rm1(find(p_bin==0 & p_bin_c==0),:),[],2)]/sqrt(length(6));
erry_l=[];erry_l=[nanstd(rm2(find(p_bin==0 & p_bin_c==0),:),[],2)]/sqrt(length(6));
temp_datax=[];temp_datay=[];diferrx=[];diferry=[];
temp_datax=nanmean(rm1(find(p_bin==0 & p_bin_c==0),:),2);
temp_datay=nanmean(rm2(find(p_bin==0 & p_bin_c==0),:),2);
diferrx=temp_datax-errx_l;
diferry=temp_datay-erry_l;
h=ploterr(temp_datax,temp_datay,errx_l,erry_l,'logx','logy','hhxy',0);
h(1).LineStyle='none';h(2).Color=[0.9 0.9 0.9];h(3).Color=[0.9 0.9 0.9];hold on;
crti=find(diferrx<0 | diferry<0);
crt_all1{i}=crti;
if isempty(crti)==0
for k=1:length(crti)
l1=line([temp_datax(crti(k)) temp_datax(crti(k))+errx_l(crti(k))],[temp_datay(crti(k)) temp_datay(crti(k))],'Color',[0.9 0.9 0.9]);
l1=line([temp_datax(crti(k)) temp_datax(crti(k))],[temp_datay(crti(k)) temp_datay(crti(k))+erry_l(crti(k))],'Color',[0.9 0.9 0.9]);
end
end
if i==1
    l1=line([0.0001 temp_datax(crti(1))],[temp_datay(crti(1)) temp_datay(crti(1))],'Color',[0.9 0.9 0.9]);
l1=line([temp_datax(crti(1)) temp_datax(crti(1))],[0.0001 temp_datay(crti(1))],'Color',[0.9 0.9 0.9]);
end
%2)
errx_l=[];errx_l=[nanstd(rm1(find(p_bin>0),:),[],2)]/sqrt(length(6));
erry_l=[];erry_l=[nanstd(rm2(find(p_bin>0),:),[],2)]/sqrt(length(6));
temp_datax=[];temp_datay=[];diferrx=[];diferry=[];
temp_datax=nanmean(rm1(find(p_bin>0),:),2);
temp_datay=nanmean(rm2(find(p_bin>0),:),2);
diferrx=temp_datax-errx_l;
diferry=temp_datay-erry_l;
h=ploterr(temp_datax,temp_datay,errx_l,erry_l,'logx','logy','hhxy',0);
h(1).LineStyle='none';h(2).Color=[0.9 0.9 0.9];h(3).Color=[0.9 0.9 0.9];hold on;
crti=find(diferrx<0 | diferry<0);
if isempty(crti)==0
for k=1:length(crti)
l1=line([temp_datax(crti(k)) temp_datax(crti(k))+errx_l(crti(k))],[temp_datay(crti(k)) temp_datay(crti(k))],'Color',[0.9 0.9 0.9]);
l1=line([temp_datax(crti(k)) temp_datax(crti(k))],[temp_datay(crti(k)) temp_datay(crti(k))+erry_l(crti(k))],'Color',[0.9 0.9 0.9]);
end
end


hold on;rf=refline(1,0);rf.Color='k';hold on;rf.LineStyle=':';


%3) 
errx_l=[];errx_l=[nanstd(rm1(find(p_bin_c>0),:),[],2)]/sqrt(length(6));
erry_l=[];erry_l=[nanstd(rm2(find(p_bin_c>0),:),[],2)]/sqrt(length(6));
temp_datax=[];temp_datay=[];diferrx=[];diferry=[];
temp_datax=nanmean(rm1(find(p_bin_c>0),:),2);
temp_datay=nanmean(rm2(find(p_bin_c>0),:),2);
diferrx=temp_datax-errx_l;
diferry=temp_datay-erry_l;
h=ploterr(temp_datax,temp_datay,errx_l,erry_l,'logx','logy','hhxy',0);
h(1).LineStyle='none';h(2).Color=[0.9 0.9 0.9];h(3).Color=[0.9 0.9 0.9];hold on;
crti=find(diferrx<0 | diferry<0);
if isempty(crti)==0
for k=1:length(crti)
l1=line([temp_datax(crti(k)) temp_datax(crti(k))+errx_l(crti(k))],[temp_datay(crti(k)) temp_datay(crti(k))],'Color',[0.9 0.9 0.9]);
l1=line([temp_datax(crti(k)) temp_datax(crti(k))],[temp_datay(crti(k)) temp_datay(crti(k))+erry_l(crti(k))],'Color',[0.9 0.9 0.9]);
end
end

loglog(nanmean(rm1(find(p_bin>0),:),2),nanmean(rm2(find(p_bin>0),:),2),'o','MarkerFaceColor',[0.8 0.8 0.8],'MarkerEdgeColor','k');hold on;
loglog(nanmean(rm1(find(p_bin_c>0),:),2),nanmean(rm2(find(p_bin_c>0),:),2),'o','MarkerFaceColor',[0.3 0.3 0.3],'MarkerEdgeColor','k');hold on;
loglog(nanmean(rm1(find(p_bin==0 & p_bin_c==0),:),2),nanmean(rm2(find(p_bin==0 & p_bin_c==0),:),2),'o','MarkerFaceColor',[1 1 1],'MarkerEdgeColor','k');hold on;axis square

xlabel('Ipsi fraction (log)');ylabel('Contra fraction (log)');set(gca,'FontSize',11);set(gca,'TickDir','out');box off;title(panel_title{i},'Color',temp_color(i,:),'FontWeight','normal');
if i==1
xlim([0.0001 0.2]);
ylim([0.0001 0.2]);
set(gca,'XLim',[0.0001 0.2],'XTick',10.^(-4:-1), ...
        'YLim',[0.0001 0.2],'YTick',10.^(-4:-1))

elseif i==2
 xlim([0.0001 0.25]);
 ylim([0.0001 0.25]);
 set(gca,'XLim',[0.0001 0.25],'XTick',10.^(-4:-1), ...
        'YLim',[0.0001 0.25],'YTick',10.^(-4:-1))
else
xlim([0.0001 0.6]);
ylim([0.0001 0.6]);
 set(gca,'XLim',[0.0001 0.6],'XTick',10.^(-4:-1), ...
        'YLim',[0.0001 0.6],'YTick',10.^(-4:-1))
end

set(gca, 'XScale','log', 'YScale','log');
axis square;
end

cd(save_folder);saveas(gcf, 'correlation_log_fraction.pdf');

%% Correlation between hemisphere areas R2 square of fraction per area avergae across injection areas
%Calculate R2 and slope of correlation 
%V1
for i=1:size(v1a_riam,2)
idx = find(~isnan(v1a_riam(:,i)));
mdl = fitlm(v1a_riam(idx,i),v1a_rcam2(idx,i));
coefs = polyfit(v1a_riam(idx,i)', v1a_rcam2(idx,i)', 1);
slope_v1(i)=coefs(1);
r2_sqv1(i)=mdl.Rsquared.Ordinary;
end
%S1
for i=1:size(s1a_riam,2)
idx = find(~isnan(s1a_riam(:,i)));
mdl = fitlm(s1a_riam(idx,i),s1a_rcam2(idx,i));
coefs = polyfit(s1a_riam(idx,i)', s1a_rcam2(idx,i)', 1);
slope_s1(i)=coefs(1);
r2_sqs1(i)=mdl.Rsquared.Ordinary;
end
%M1
for i=1:size(m1a_riam,2)
idx = find(~isnan(m1a_riam(:,i)));
mdl = fitlm(m1a_riam(idx,i),m1a_rcam2(idx,i));
coefs = polyfit(m1a_riam(idx,i)', m1a_rcam2(idx,i)', 1);
slope_m1(i)=coefs(1);
r2_sqm1(i)=mdl.Rsquared.Ordinary;
end
%R2
tm1=[];tm2=[];tm3=[];tm1=r2_sqv1;tm2=r2_sqs1;tm3=r2_sqm1;
%plot
fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 500, 170, 250]);set(gcf,'color','w');tiledlayout("horizontal")
b1=bar(1,nanmean(tm1),0.7);hold on;b1.FaceColor=v1_color;b2=bar(2,nanmean(tm2));hold on;b2.FaceColor=s1_color;
b3=bar(3,nanmean(tm3),0.7);b3.FaceColor=m1_color;
set(b1,'ShowBaseLine','off');set(b2,'ShowBaseLine','off');set(b3,'ShowBaseLine','off');
r=1;rng(1);r1 = r-0.01 + (0.2)*rand(length(tm1),1);
sc1=scatter(r1,tm1,15,'k','filled');hold on;sc1.MarkerEdgeColor=[1 1 1];sc1.MarkerEdgeAlpha=0.5;sc1.MarkerFaceColor=[0.2 0.2 0.2];
 r=2;rng(1);r1 = r-0.01 + (0.2)*rand(length(tm2),1);
 sc1=scatter(r1,tm2,15,'k','filled');hold on;sc1.MarkerEdgeColor=[1 1 1];sc1.MarkerEdgeAlpha=0.5;sc1.MarkerFaceColor=[0.2 0.2 0.2];
 r=3;rng(1);r1 = r-0.01 + (0.2)*rand(length(tm3),1);
 sc1=scatter(r1,tm3,15,'k','filled');hold on;sc1.MarkerEdgeColor=[1 1 1];sc1.MarkerEdgeAlpha=0.5;sc1.MarkerFaceColor=[0.2 0.2 0.2];xticks([1:3]);hold on;
hold on;errorbar([1 2 3],[nanmean(tm1) nanmean(tm2) nanmean(tm3)],[nanstd(tm1)/sqrt(size(tm1,2)) nanstd(tm2)/sqrt(size(tm2,2)) nanstd(tm3)/sqrt(size(tm3,2))]...
    , 'LineStyle', 'none', ... 
        'Color', 'k', 'LineWidth', 1.2);hold on;
xticklabels({'VISp','SSpbf','MOp'});%ylim([-0.05 0.75]);
ylabel('Correlation R2 ipsi-contra');set(gca,'FontSize',11);set(gca,'TickDir','out');box off;xtickangle(45);
 h = gca;h.XAxis.Visible = 'off'
offsetAxes;

  t1=text(0.4,-0.2,'VISp','FontSize',11);set(t1,'Rotation',45);
   t1=text(0.8,-0.3,'SSp-bfd','FontSize',11);set(t1,'Rotation',45);
    t1=text(2.4,-0.2,'MOp','FontSize',11);set(t1,'Rotation',45);

 disp('VISp R2')
nanmean(tm1)
nanstd(tm1)/sqrt(size(tm1,2))
 disp('SSp R2')
nanmean(tm2)
nanstd(tm2)/sqrt(size(tm2,2))
 disp('MOp R2')
nanmean(tm3)
nanstd(tm3)/sqrt(size(tm3,2))
% [p,tbl,stats] = anova1([tm1' tm2' tm3'])
% presults = multcompare(stats)
%Slope
tm1=[];tm2=[];tm3=[];tm1=slope_v1;tm2=slope_s1;tm3=slope_m1;
%% save
cd(save_folder);saveas(gcf, 'r2_correlation_hemi.pdf');

%% Figure 2
%%  ipsi vs contra shown in modules, sorted by largest module
%dat_all={};dat_all={v1a_riam v1a_rcam; s1a_riam s1a_rcam; m1a_riam m1a_rcam};
dat_all={};dat_all={v1a_riam v1a_rcam2; s1a_riam s1a_rcam2; m1a_riam m1a_rcam2};

temp_color=[v1_color ;s1_color; m1_color];
panel_title={'VISp','SSp-bfd','MOp'};
name_sub={'modules_VISp.pdf', 'modules_SSp.pdf' ,'modules_MOp.pdf'};
order_modules={visual_idx  medial_idx lateral_idx aud_idx frontal_idx somamo_idx;...
    somamo_idx visual_idx frontal_idx lateral_idx aud_idx medial_idx;...
    somamo_idx frontal_idx lateral_idx medial_idx visual_idx aud_idx};
order_modules_names={'Vis','Med','Lat','Aud','Prefron','SoMo';...
    'SoMo','Vis','Prefron','Lat','Aud','Med';...
    'SoMo','Prefron','Lat','Med','Vis','Aud'};
fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [200, 500, 800, 240]);set(gcf,'color','w');t=tiledlayout("horizontal",'TileSpacing','Compact');
for m=1:3
    figure(fig7)
    nexttile
%fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [200+m*350, 500, 500, 240]);set(gcf,'color','w');t=tiledlayout("horizontal",'TileSpacing','Compact');
t.Padding = 'compact';

rm1=dat_all{m,1};
rm2=dat_all{m,2};
p1=[];p1=nansum(rm1(order_modules{m,1},:));p2=[];p2=nansum(rm1(order_modules{m,2},:));p3=[];p3=nansum(rm1(order_modules{m,3},:));
p4=[];p4=nansum(rm1(order_modules{m,4},:));p5=[];p5=nansum(rm1(order_modules{m,5},:));p6=[];p6=nansum(rm1(order_modules{m,6},:));
p1_1=[];p1_1=nansum(rm2(order_modules{m,1},:));p2_1=[];p2_1=nansum(rm2(order_modules{m,2},:));p3_1=[];p3_1=nansum(rm2(order_modules{m,3},:));
p4_1=[];p4_1=nansum(rm2(order_modules{m,4},:));p5_1=[];p5_1=nansum(rm2(order_modules{m,5},:));p6_1=[];p6_1=nansum(rm2(order_modules{m,6},:));
temp_p= [p1' p2' p3' p4' p5' p6'];
temp_p2= [p1_1' p2_1' p3_1' p4_1' p5_1' p6_1'];
% [uuj uui]=sort(temp_p)
% [uuj uuc]=sort(temp_p2)
% [pw ,ee, stats]=signrank(uui(:),uuc(:))
title(panel_title{m},'FontWeight','normal','Color',temp_color(m,:));

    for i=1:6
        hold on;
        b1=bar(i,[nanmean(temp_p(:,i));nanmean(temp_p2(:,i))],1);hold on;b1(1).FaceColor=[0.8 0.8 0.8];b1(2).FaceColor=[0.3 0.3 0.3];set(b1,'ShowBaseLine','off');
        b1(1).EdgeColor=[1 1 1];b1(2).EdgeColor=[1 1 1]
       % disp(['stats ipsi vs contra' panel_title{m}])
       % [u p11]=ttest(temp_p(:,i),temp_p2(:,i))
       
         r=i;rng(i);r1 = r-0.01 + (0.2)*rand(length(temp_p(:,i)),1);
       r2=i;rng(i);r2 = r2-0.01 + (0.2)*rand(length(temp_p2(:,i)),1);
         sc1=scatter(r1-0.2,temp_p(:,i),6,'k','filled');hold on;sc1.MarkerEdgeColor=[1 1 1];sc1.MarkerEdgeAlpha=0.5;sc1.MarkerFaceColor=[0.2 0.2 0.2];
         sc1=scatter(r2+0.15,temp_p2(:,i),6,'k','filled');hold on;sc1.MarkerEdgeColor=[1 1 1];sc1.MarkerEdgeAlpha=0.5;sc1.MarkerFaceColor=[0.2 0.2 0.2];

         hold on;h=errorbar([i-0.15],[nanmean(temp_p(:,i))],[nanstd(temp_p(:,i))/sqrt(length(temp_p(:,i)))], 'LineStyle', 'none','Color', 'k', 'LineWidth', 1,'CapSize',3);
         hold on;h=errorbar([i+0.15],[nanmean(temp_p2(:,i))],[nanstd(temp_p2(:,i))/sqrt(length(temp_p2(:,i)))], 'LineStyle', 'none','Color', 'k', 'LineWidth', 1,'CapSize',3);
         %disp(['module 1 vs others stats' panel_title{m}]);
         %[u p12]=ttest2([temp_p(:,1) ;temp_p2(:,1)],temp_p2(:,i))
       
        % hold on;errorbar([i],[nanmean(temp_p(:,i))],[nanstd(temp_p(:,i))/sqrt(length(temp_p(:,i)))]...
        % , 'LineStyle', 'none', ... 
        % 'Color', 'k', 'LineWidth', 1.2);hold on;
        
    end

   
 if m==1
% 
 xticks([1:6]);hold on;box off;xticklabels(order_modules_names(m,:));ylim([0 1]);set(gca,'FontSize',11);set(gca,'TickDir','out');box off;xtickangle(45);
 ylabel('Fraction of cells');
 text(3.5,0.9,'ipsi','Color',[0.7 0.7 0.7]);
 text(3.5,0.8,'contra (ex homotopic)','Color',[0.3 0.3 0.3]);
 else
 xticks([1:6]);xticklabels(order_modules_names(m,:));ylim([0 1]);set(gca,'FontSize',11);set(gca,'TickDir','out');box off;xtickangle(45);
 h = gca;h.YAxis.Visible = 'off';  
 end
   xlim([0 6.5]);
offsetAxes;

%save
%cd(save_folder);saveas(gcf, name_sub{m});
%figure;
 %disp(['modules ipsi+contra against each other ' panel_title{m}])
      % [p,tbl,stats] = anova1([[temp_p(:,1) ;temp_p2(:,1)] [temp_p(:,2) ;temp_p2(:,2)]...
      %     [temp_p(:,3) ;temp_p2(:,3)] [temp_p(:,4) ;temp_p2(:,4)] [temp_p(:,5) ;temp_p2(:,5)] [temp_p(:,6) ;temp_p2(:,6)]])
      % presults = multcompare(stats)
  disp(['modules ipsi against next higher ' panel_title{m}])  
  nanmean(temp_p(:,1))
  [nanstd(temp_p(:,1))/sqrt(length(temp_p(:,1)))]
nanmean(temp_p(:,2))
 [nanstd(temp_p(:,2))/sqrt(length(temp_p(:,2)))]
       [p u]=ttest2(temp_p(:,1),temp_p(:,2))
 disp(['modules contra against next higher ' panel_title{m}])  
 if m==1
       nanmean(temp_p2(:,1))
  [nanstd(temp_p2(:,1))/sqrt(length(temp_p2(:,1)))]
nanmean(temp_p2(:,3))
 [nanstd(temp_p2(:,3))/sqrt(length(temp_p2(:,3)))]
       [p u]=ttest2(temp_p2(:,1),temp_p2(:,3))

 elseif m==2
        nanmean(temp_p2(:,1))
  [nanstd(temp_p2(:,1))/sqrt(length(temp_p2(:,1)))]
nanmean(temp_p2(:,4))
 [nanstd(temp_p2(:,4))/sqrt(length(temp_p2(:,4)))]
       [p u]=ttest2(temp_p2(:,1),temp_p2(:,4))
 else m==3
       nanmean(temp_p2(:,1))
  [nanstd(temp_p2(:,1))/sqrt(length(temp_p2(:,1)))]
nanmean(temp_p2(:,3))
 [nanstd(temp_p2(:,3))/sqrt(length(temp_p2(:,3)))]
       [p u]=ttest2(temp_p2(:,1),temp_p2(:,3))

 end
end
%% Log individual areas RANK ORDER
temp_color=linspecer(6);
cortex_abb=cortex_names(:,2);
idx_modules={frontal_idx lateral_idx somamo_idx visual_idx medial_idx aud_idx};
module_names={'Pfron','Lat','SoMo','Vis','Med','Aud'};
frontal_idx=[1:8];lateral_idx=[9:16];somamo_idx=[17:26];visual_idx=[27:36];medial_idx=[37:39];aud_idx=[40:43];lateral_idx2=[44:45]
areas_idx=ones(45,1);
areas_modul=[];
areas_modul=[areas_idx(frontal_idx)*1; areas_idx(lateral_idx)*2 ; areas_idx(somamo_idx)*3; ...
    areas_idx(visual_idx)*4; areas_idx(medial_idx)*5; areas_idx(aud_idx)*6 ; areas_idx(lateral_idx2)*2];
panel_tit={'VISp','SSp-bfd','MOp'};
dat_all={};dat_all={v1a_riam s1a_riam m1a_riam};
%excluding homotopic
dat_all2={};dat_all2={v1a_rcam2 s1a_rcam2 m1a_rcam2};

fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [100, 300, 1200, 700]);set(gcf,'color','w');tiledlayout("horizontal")
for j=1:3    
    temp_cortex1=[];temp_cortex1=nanmean(dat_all{j},2);
    temp_cortex2=[];temp_cortex2=nanmean(dat_all2{j},2);

    nexttile
    for k=1:45
    if areas_modul(k)==1
    sl=semilogy([1 2],[temp_cortex1(k) temp_cortex2(k)],'-o','Color',temp_color(1,:));
    sl.MarkerFaceColor=temp_color(1,:);
    text(0.75,temp_cortex1(k),cortex_abb{k},'Color',temp_color(1,:),'FontSize',7);
     text(2.05,temp_cortex2(k),cortex_abb{k},'Color',temp_color(1,:),'FontSize',7);
    hold on;
    elseif areas_modul(k)==2
    sl=semilogy([1 2],[temp_cortex1(k) temp_cortex2(k)],'-o','Color',temp_color(2,:));
    sl.MarkerFaceColor=temp_color(2,:);
    text(0.75,temp_cortex1(k),cortex_abb{k},'Color',temp_color(2,:),'FontSize',7);
     text(2.05,temp_cortex2(k),cortex_abb{k},'Color',temp_color(2,:),'FontSize',7);
    hold on;
    elseif areas_modul(k)==3
    sl=semilogy([1 2],[temp_cortex1(k) temp_cortex2(k)],'-o','Color',temp_color(3,:));
    sl.MarkerFaceColor=temp_color(3,:);
     text(0.75,temp_cortex1(k),cortex_abb{k},'Color',temp_color(3,:),'FontSize',7);
     text(2.05,temp_cortex2(k),cortex_abb{k},'Color',temp_color(3,:),'FontSize',7);
    hold on;
     elseif areas_modul(k)==4
    sl=semilogy([1 2],[temp_cortex1(k) temp_cortex2(k)],'-o','Color',temp_color(4,:));
    sl.MarkerFaceColor=temp_color(4,:);
     text(0.75,temp_cortex1(k),cortex_abb{k},'Color',temp_color(4,:),'FontSize',7);
     text(2.05,temp_cortex2(k),cortex_abb{k},'Color',temp_color(4,:),'FontSize',7);
    hold on;
     elseif areas_modul(k)==5
    sl=semilogy([1 2],[temp_cortex1(k) temp_cortex2(k)],'-o','Color',temp_color(5,:));
    sl.MarkerFaceColor=temp_color(5,:);
    text(0.75,temp_cortex1(k),cortex_abb{k},'Color',temp_color(5,:),'FontSize',7);
     text(2.05,temp_cortex2(k),cortex_abb{k},'Color',temp_color(5,:),'FontSize',7);
    hold on;
    else areas_modul(k)==6
    sl=semilogy([1 2],[temp_cortex1(k) temp_cortex2(k)],'-o','Color',temp_color(6,:));
    sl.MarkerFaceColor=temp_color(6,:);
    text(0.75,temp_cortex1(k),cortex_abb{k},'Color',temp_color(6,:),'FontSize',7);
     text(2.05,temp_cortex2(k),cortex_abb{k},'Color',temp_color(6,:),'FontSize',7);
    hold on;
    end
    end
    box off;
    axis off;
    %ylim([0.0001 0.6])
    xlim([0.7 2.3]);
    title(panel_tit{j})

% figure(1)
% histogram(rank_data(:,1)-rank_data(:,2));hold on
if j==1 | j==2
text(0.8,max(temp_cortex2)+0.03,'ipsi','Color',[0.7 0.7 0.7],'FontSize',11);
text(1.8,max(temp_cortex2)+0.03,'contra','Color',[0.3 0.3 0.3],'FontSize',11);
else j==3
  text(0.8,max(temp_cortex2)+0.55,'ipsi','Color',[0.7 0.7 0.7],'FontSize',11);
text(1.8,max(temp_cortex2)+0.55,'contra','Color',[0.3 0.3 0.3],'FontSize',11);  
end

end

cd(save_folder);saveas(gcf, 'log_connecting_areas.pdf');
%% Rank order TWM THIS
%temp_color=[[0 0.6 1] ;[0 0.4 0.4]; [0.6 0.2 0.4]; [1 0.2 0.4];	[1 0.4 1];[0.6 1 1]];
temp_color=linspecer(6);
idx_modules={frontal_idx lateral_idx somamo_idx visual_idx medial_idx aud_idx};
module_names={'Pfron','Lat','SoMo','Vis','Med','Aud'};
frontal_idx=[1:8];lateral_idx=[9:16];somamo_idx=[17:26];visual_idx=[27:36];medial_idx=[37:39];aud_idx=[40:43];lateral_idx2=[44:45]
areas_idx=ones(45,1);
areas_modul=[];
areas_modul=[areas_idx(frontal_idx)*1; areas_idx(lateral_idx)*2 ; areas_idx(somamo_idx)*3; ...
    areas_idx(visual_idx)*4; areas_idx(medial_idx)*5; areas_idx(aud_idx)*6 ; areas_idx(lateral_idx2)*2];
panel_tit={'VISp','SSp-bfd','MOp'};
dat_all={};dat_all={v1a_riam s1a_riam m1a_riam};
%excluding homotopic
dat_all2={};dat_all2={v1a_rcam2 s1a_rcam2 m1a_rcam2};

fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [100, 300, 1400, 700]);set(gcf,'color','w');tiledlayout("horizontal")
for j=1:3    
    temp_cortex1=[];temp_cortex1=nanmean(dat_all{j},2);
    temp_cortex2=[];temp_cortex2=nanmean(dat_all2{j},2);
 so_v1a_riam=[];idxi=[];cs1=[];so_v1a_rcam=[];idxc=[];tcon=[];
 %tcon=find(temp_cortex2>0 & temp_cortex1>0);
 tcon=[];
 tcon=find(~isnan(temp_cortex1));
[so_v1a_riam idxi]=sort(temp_cortex1(tcon),'descend');
[so_v1a_rcam idxc]=sort(temp_cortex2(tcon),'descend');
temp1=[];temp2=[];areas_coo1=[];areas_coo2=[];
temp1=areas_modul(tcon);temp2=areas_modul(tcon);temp3=cortex_names(tcon,2);temp4=cortex_names(tcon,2);
areas_coo1=temp1(idxi);areas_coo2=temp2(idxc);names_sel=temp3(idxi);names_sel2=temp4(idxc);
mom1=[length(tcon):-1:1];color_matrix=[];rank_data=[];
nexttile
for k=1:length(idxi)
    i2=[];c1=[];
    i2=find(idxi(k)==idxc);
    p1=plot([2 3],[mom1(k) mom1(i2)]);hold on;
    p3=plot([2],[mom1(k)],'o');hold on;
    p4=plot([3],[mom1(i2)],'o');
    if so_v1a_riam(k)==0 
        p1.LineStyle='none';
        p3.MarkerSize=0.01;
        values_dot(k)=0;
    else
    p3.MarkerSize=10+log(so_v1a_riam(k));
    values_dot(k)=(so_v1a_riam(k));
    end
      if so_v1a_rcam(i2)==0 
          p4.MarkerSize=0.01;
          values_dot(k)=0;
      else
           p4.MarkerSize=10+log(so_v1a_rcam(i2));
           values_dot(k)=(so_v1a_riam(k));
      end
   hold on;t1=text(1.5,mom1(k),names_sel{k});
      hold on;
     t2=text(3.2,mom1(i2),names_sel2{i2});
    rank_data(k,:)=[mom1(k) mom1(i2)];
    if areas_coo1(k)==1
    col1=temp_color(1,:)
    p1.MarkerFaceColor=col1;p1.MarkerEdgeColor=col1;p3.MarkerFaceColor=col1;p4.MarkerFaceColor=col1;p3.MarkerEdgeColor=col1;p4.MarkerEdgeColor=col1;
    p1.Color=col1;t1.Color=col1;t2.Color=col1;color_matrix(k,:)=col1;
   elseif areas_coo1(k)==2
    col1=temp_color(2,:);
    p1.MarkerFaceColor=col1;p1.MarkerEdgeColor=col1;p3.MarkerFaceColor=col1;p4.MarkerFaceColor=col1;p3.MarkerEdgeColor=col1;p4.MarkerEdgeColor=col1;
    p1.Color=col1;t1.Color=col1;t2.Color=col1;
    color_matrix(k,:)=col1;
        
    elseif areas_coo1(k)==3
          col1=temp_color(3,:);
    p1.MarkerFaceColor=col1;
    p1.MarkerEdgeColor=col1;
      p3.MarkerFaceColor=col1;
    p4.MarkerFaceColor=col1;
       p3.MarkerEdgeColor=col1;
    p4.MarkerEdgeColor=col1;
    p1.Color=col1;
    t1.Color=col1;
    t2.Color=col1;
    color_matrix(k,:)=col1;
        
    elseif areas_coo1(k)==4
         col1=temp_color(4,:);
    p1.MarkerFaceColor=col1;
    p1.MarkerEdgeColor=col1;
      p3.MarkerFaceColor=col1;
    p4.MarkerFaceColor=col1;
    p3.MarkerEdgeColor=col1;
    p4.MarkerEdgeColor=col1;
    p1.Color=col1;
    t1.Color=col1;
    t2.Color=col1;
    color_matrix(k,:)=col1;
    elseif areas_coo1(k)==5
          col1=temp_color(5,:);
    p1.MarkerFaceColor=col1;
    p1.MarkerEdgeColor=col1;
      p3.MarkerFaceColor=col1;
    p4.MarkerFaceColor=col1;
    p3.MarkerEdgeColor=col1;
    p4.MarkerEdgeColor=col1;
    p1.Color=col1;
    t1.Color=col1;
    t2.Color=col1;
    color_matrix(k,:)=col1;
          
    else areas_coo1(k)==6
          col1=temp_color(6,:);
    p1.MarkerFaceColor=col1;
    p1.MarkerEdgeColor=col1;
      p3.MarkerFaceColor=col1;
    p4.MarkerFaceColor=col1;
    p3.MarkerEdgeColor=col1;
    p4.MarkerEdgeColor=col1;
    p1.Color=col1;
    t1.Color=col1;
    t2.Color=col1;
    color_matrix(k,:)=col1;
       
    end
hold on;
box off;

xlim([1.5 3.5])
%offsetAxes;
 h = gca;h.XAxis.Visible = 'off';
 axis off
end
title(panel_tit{j})
diff_rank{j}=rank_data(:,1)-rank_data(:,2);
% figure(1)
% histogram(rank_data(:,1)-rank_data(:,2));hold on
text(1.8,47,'ipsi','Color',[0.7 0.7 0.7],'FontSize',11);
text(2.9,47,'contra','Color',[0.3 0.3 0.3],'FontSize',11);
set(gca,'FontSize',11)
if j==3
 
end
end
temp_scale=[32:2:40];
temp_thick=[2:2:10];
temp_thicktext={'0.0004','','','','0.4'};
nexttile
   
    for oo=1:length(temp_scale)
   pt= plot(1,temp_scale(oo),'o','Color','k');
   pt.MarkerSize=[temp_thick(oo)];hold on;
pt.MarkerFaceColor='k';
text(1.1,temp_scale(oo),temp_thicktext{oo},'FontSize',11);
    end
    text(0.9,43,'log fraction','FontSize',11)
    ylim([1 44])
    box off;axis off;
    set(gca,'FontSize',11)
%% 
cd(save_folder);saveas(gcf, 'rank_order_connectome.pdf');

%% 
temp_color3=[v1_color ;s1_color; m1_color];
fig8= figure;set(fig8, 'Name', 'Barplot groups');set(fig8, 'Position', [200, 1000, 400, 400]);set(gcf,'color','w');
for j=1:3
    hold on;h1=histogram(diff_rank{j},10);h1.FaceColor=temp_color3(j,:);h1.BinWidth = 2;
end
xlabel('Delta rank order ipsi - contra');
ylabel('Counts')




%% FIGURE 3

%% Plot barplot for all regions all areas RETRO (animals=18, areas=45)
%means
means_layers=[nanmean(i_animal(2,:)) nanmean(c_animal(2,:));nanmean(i_animal(3,:)) nanmean(c_animal(3,:));...
    nanmean(i_animal(4,:)) nanmean(c_animal(4,:));nanmean(i_animal(5,:)) nanmean(c_animal(5,:));nanmean(i_animal(6,:)) nanmean(c_animal(6,:))];
semd=sqrt(length(i_animal));
%error
err_layers=[[nanstd(i_animal(2,:))/semd nanstd(c_animal(2,:))/semd];[nanstd(i_animal(3,:))/semd nanstd(c_animal(3,:))/semd];...
    [nanstd(i_animal(4,:))/semd nanstd(c_animal(4,:))/semd];[nanstd(i_animal(5,:))/semd nanstd(c_animal(5,:))/semd] ;[nanstd(i_animal(6,:))/semd nanstd(c_animal(6,:))/semd]];
%plotting starts here 
fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 500, 450, 250]);set(gcf,'color','w');
%bars
b=bar(means_layers,0.7);
b(1).FaceColor=[0.8 0.8 0.8];b(2).FaceColor=[0.3 0.3 0.3];b(1).EdgeColor=[1 1 1];b(2).EdgeColor=[1 1 1];hold on;set(b,'ShowBaseLine','off');
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
;ylabel('Fraction of neurons');set(gca,'FontSize',11);set(gca,'TickDir','out');box off;
%[p1]=signrank(i_animal(6,:),c_animal(6,:))
[u p1]=ttest(i_animal(5,:),c_animal(5,:))
hold on;text(0.75,0.6,['***'],'FontSize',18);hold on;text(3.75,0.6,['***'],'FontSize',18);legend({'ipsi','contra'},"Location","northeast");legend boxoff;xlim([0.5 5.5]);
offsetAxes;
h = gca;h.XAxis.Visible = 'off';  
%xticklabels({'L2/3','L4','L5','L6a','L6b'})
text(0.75,-0.04,'L2/3','FontSize',11);text(1.9,-0.04,'L4','FontSize',11);text(2.9,-0.04,'L5','FontSize',11);text(3.9,-0.04,'L6a','FontSize',11);text(4.9,-0.04,'L6b','FontSize',11);
%hold on;text(0.75,0.6,['**'],'FontSize',18);hold on;text(3.75,0.6,['*'],'FontSize',18);legend({'ipsi','contra'},"Location","northeast");legend boxoff;title('MOp','Color',m1_color)
%save
cd(save_folder);saveas(gcf, 'fraction_cells_layers.pdf');
%% same plot differently sorted for layers 
labels_b={'L2/3','L4','L5','L6a' 'L6b'};
anim_incl={1:6 7:12 13:18};
temp_color=[v1_color ;s1_color; m1_color];
panel_tit={'VISp','SSp-bfd','MOp'};

fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 500, 1200, 250]);set(gcf,'color','w');tiledlayout("horizontal");
for j=1:3
incl1=[];incl1=anim_incl{j}


%means
means_layers=[];err_layers=[];
means_layers=[nanmean(i_animal(2,incl1)) nanmean(i_animal(3,incl1)) nanmean(i_animal(4,incl1)) nanmean(i_animal(5,incl1)) nanmean(i_animal(6,incl1));...
    nanmean(c_animal(2,incl1)) nanmean(c_animal(3,incl1)) nanmean(c_animal(4,incl1)) nanmean(c_animal(5,incl1)) nanmean(c_animal(6,incl1)) ];
semd=sqrt(length(i_animal(:,incl1)));
%error
err_layers=[nanstd(i_animal(2,incl1))/semd nanstd(i_animal(3,incl1))/semd nanstd(i_animal(4,incl1))/semd nanstd(i_animal(5,incl1))/semd nanstd(i_animal(6,incl1))/semd; ...
    nanstd(c_animal(2,incl1))/semd nanstd(c_animal(3,incl1))/semd nanstd(c_animal(4,incl1))/semd nanstd(c_animal(5,incl1))/semd nanstd(c_animal(6,incl1))/semd];

nexttile
%bars
b1=bar(1,means_layers(1,:),0.7);
b1(1).FaceColor=[0.8 0.8 0.8];b1(2).FaceColor=[0.8 0.8 0.8];b1(3).FaceColor=[0.8 0.8 0.8];b1(4).FaceColor=[0.8 0.8 0.8];b1(5).FaceColor=[0.8 0.8 0.8];
b1(1).EdgeColor=[1 1 1];b1(2).EdgeColor=[1 1 1];b1(3).EdgeColor=[1 1 1];b1(4).EdgeColor=[1 1 1];b1(5).EdgeColor=[1 1 1];hold on;set(b1,'ShowBaseLine','off');
b2=bar(2,means_layers(2,:),0.7);
b2(1).FaceColor=[0.3 0.3 0.3];b2(2).FaceColor=[0.3 0.3 0.3];b2(3).FaceColor=[0.3 0.3 0.3];b2(4).FaceColor=[0.3 0.3 0.3];b2(5).FaceColor=[0.3 0.3 0.3];
b2(1).EdgeColor=[1 1 1];b2(2).EdgeColor=[1 1 1];b2(3).EdgeColor=[1 1 1];b2(4).EdgeColor=[1 1 1];b2(5).EdgeColor=[1 1 1];hold on;set(b2,'ShowBaseLine','off');
%indivdual animals
r=1;
rng(1);r1 = r-0.35 + (0.1)*rand(length(incl1),1);
rng(1);r2 = r -0.2+ (0.1)*rand(length(incl1),1);
rng(1);r3 = r-0.05 + (0.1)* rand(length(incl1),1);
rng(1);r4 = r+0.12 + (0.1)*rand(length(incl1),1);
rng(1);r5 = r+0.27 + (0.1)*rand(length(incl1),1);
sc1=scatter([r1 r2 r3 r4 r5],i_animal([2 3 4 5 6],incl1)',5,'k','filled');hold on;
for oo=1:5
sc1(oo).MarkerEdgeColor=[0.5 0.5 0.5];sc1(oo).MarkerFaceColor=[0.2 0.2 0.2];sc1(oo).MarkerFaceAlpha=[0.5];sc1(oo).MarkerEdgeAlpha=[0.5];
end
hold on;
r=2;
rng(1);r1 = r-0.35 + (0.1)*rand(length(incl1),1);
rng(1);r2 = r -0.2+ (0.1)*rand(length(incl1),1);
rng(1);r3 = r-0.05 + (0.1)* rand(length(incl1),1);
rng(1);r4 = r+0.12 + (0.1)*rand(length(incl1),1);
rng(1);r5 = r+0.27 + (0.1)*rand(length(incl1),1);
sc1=scatter([r1 r2 r3 r4 r5],c_animal([2 3 4 5 6],incl1)',5,'k','filled');hold on;
for oo=1:5
sc1(oo).MarkerEdgeColor=[0.5 0.5 0.5];sc1(oo).MarkerFaceColor=[0.2 0.2 0.2];sc1(oo).MarkerFaceAlpha=[0.5];sc1(oo).MarkerEdgeAlpha=[0.5];
end
  hold on;h=errorbar([1-0.3],[means_layers(1,1)],[err_layers(1,1)], 'LineStyle', 'none','Color', 'k', 'LineWidth', 2,'CapSize',0);
   hold on;h=errorbar([1-0.15],[means_layers(1,2)],[err_layers(1,2)], 'LineStyle', 'none','Color', 'k', 'LineWidth', 2,'CapSize',0);
 hold on;h=errorbar([1],[means_layers(1,3)],[err_layers(1,3)], 'LineStyle', 'none','Color', 'k', 'LineWidth', 2,'CapSize',0);
 hold on;h=errorbar([1+0.15],[means_layers(1,4)],[err_layers(1,4)], 'LineStyle', 'none','Color', 'k', 'LineWidth', 2,'CapSize',0);
 hold on;h=errorbar([1+0.3],[means_layers(1,5)],[err_layers(1,5)], 'LineStyle', 'none','Color', 'k', 'LineWidth', 2,'CapSize',0);

   hold on;h=errorbar([2-0.3],[means_layers(2,1)],[err_layers(2,1)], 'LineStyle', 'none','Color', 'k', 'LineWidth', 2,'CapSize',0);
   hold on;h=errorbar([2-0.15],[means_layers(2,2)],[err_layers(2,2)], 'LineStyle', 'none','Color', 'k', 'LineWidth', 2,'CapSize',0);
 hold on;h=errorbar([2],[means_layers(2,3)],[err_layers(2,3)], 'LineStyle', 'none','Color', 'k', 'LineWidth', 2,'CapSize',0);
 hold on;h=errorbar([2+0.15],[means_layers(2,4)],[err_layers(2,4)], 'LineStyle', 'none','Color', 'k', 'LineWidth', 2,'CapSize',0);
 hold on;h=errorbar([2+0.3],[means_layers(2,5)],[err_layers(2,5)], 'LineStyle', 'none','Color', 'k', 'LineWidth', 2,'CapSize',0);
  box off;

ylabel('Fraction of neurons');set(gca,'FontSize',11);set(gca,'TickDir','out');box off;
%[p1]=signrank(i_animal(6,:),c_animal(6,:))
%[u p1]=ttest(i_animal(1,:),i_animal(5,:))
hold on;
%text(0.75,0.6,['***'],'FontSize',18);hold on;text(3.75,0.6,['***'],'FontSize',18);
%legend({'ipsi','contra'},"Location","northeast");legend boxoff;
xlim([0.5 2.5]);
%offsetAxes;
h = gca;h.XAxis.Visible = 'off';  
text(0.68,0.62,'ipsi','FontSize',11,'Color',[0.7 0.7 0.7]);
text(0.68,0.58,'contra','FontSize',11,'Color',[0.3 0.3 0.3]);

t1=text(0.68,-0.09,'L2/3','FontSize',11);set(t1,'Rotation',90);
t1=text(0.84,-0.075,'L4','FontSize',11);set(t1,'Rotation',90);
t1=text(1,-0.075,'L5','FontSize',11);set(t1,'Rotation',90);
t1=text(1.14,-0.09,'L6a','FontSize',11);set(t1,'Rotation',90);
t1=text(1.30,-0.09,'L6b','FontSize',11);set(t1,'Rotation',90);

t1=text(1.68,-0.09,'L2/3','FontSize',11);set(t1,'Rotation',90);
t1=text(1.84,-0.075,'L4','FontSize',11);set(t1,'Rotation',90);
t1=text(2,-0.075,'L5','FontSize',11);set(t1,'Rotation',90);
t1=text(2.14,-0.09,'L6a','FontSize',11);set(t1,'Rotation',90);
t1=text(2.30,-0.09,'L6b','FontSize',11);set(t1,'Rotation',90);
title(panel_tit{j},'FontWeight','normal','Color',temp_color(j,:));
%hold on;text(0.75,0.6,['**'],'FontSize',18);hold on;text(3.75,0.6,['*'],'FontSize',18);legend({'ipsi','contra'},"Location","northeast");legend boxoff;title('MOp','Color',m1_color)
%save  
% line([0.8 1.25],[0.6 0.6],'Color','k');line([1.8 2.25],[0.6 0.6],'Color','k');
% text(1.95,0.63,['***'],'FontSize',18);hold on;text(1,0.64,['n.s'],'FontSize',11);
end
cd(save_folder);saveas(gcf, 'layerFraction_injection_areas.pdf');
%% 

 [p,tbl,stats] = anova1([i_animal(2,:)' i_animal(4,:)' i_animal(5,:)'])
  presults = multcompare(stats)

%% Heatmaps L6a sorted for ipsi and contra
name_sub={'VISpi_heatmaps.pdf', 'SSpi_heatmaps.pdf' ,'MOpi_heatmaps.pdf'};
%IPSI
%dat_all={};dat_all={c_v1aam c_s1aam c_m1aam};
%dat_all2={};dat_all2={v1a_rcam s1a_rcam m1a_rcam};
dat_all={};dat_all={i_v1aafm i_s1aafm i_m1aafm};
dat_all2={};dat_all2={v1a_riam s1a_riam m1a_riam};
temp_color=[v1_color ;s1_color; m1_color];
%panel_tit={'contra VISp','contra SSp-bfd','contra MOp'};
panel_tit={'VISp','SSp-bfd','MOp'};
cortex_abb=[];cortex_abb=cortex_names(:,2);

for j=1:3
p1=[];p1=nanmean(dat_all{j},3);
%p2=[];p2=nanmean(dat_all3{j},3)
idx_ke=find(sum(dat_all2{j}>0,2)>1);
% idx_ke2=find(sum(dat_all4{j}>0,2)>1);
p1e=[];p1e=p1(:,idx_ke);
% p2e=[];p2e=p2(:,idx_ke);
cortex_abb_e=[];
cortex_abb_e=cortex_abb(idx_ke);
ind=[];
[X, ind]=nanmax(p1e,[],1);
l6_an=[];l5_an=[];l23_an=[];
l6_an=find(ind==5);
l5_an=find(ind==4);
l23_an=find(ind==2);
kk1=[];kk2=[];kk3=[];
[l23_p1 kk1 ]=sort(p1e(2,l23_an));
[l5_p1 kk2]=sort(p1e(4,l5_an));
[l6_p1 kk3]=sort(p1e(5,l6_an));

tem_cm=[];tem_cm=[p1e(2:end,l23_an(kk1)) ones(5,2)*1 p1e(2:end,l5_an(kk2)) ones(5,2)*1 p1e(2:end,l6_an(kk3))];
%tem_cm2=[];tem_cm2=[p2e(2:end,l23_an(kk1)) zeros(5,2)*NaN p2e(2:end,l5_an(kk2)) zeros(5,2)*NaN p2e(2:end,l6_an(kk3))];
p3=[];p3=[cortex_abb_e(l23_an(kk1)) ;'.';'.';cortex_abb_e(l5_an(kk2));'.';'.' ;cortex_abb_e(l6_an(kk3))];

fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 1000-250*j, 750, 160]);set(gcf,'color','w');tiledlayout("horizontal")
nexttile
h=imagesc(tem_cm);
set(h, 'AlphaData', ~isnan(tem_cm))
cmap(temp_color(j,:),100,0.5,0.5);colorbar;
[hh oo]=find(isnan(tem_cm));
coord_b=[oo hh];
for u=1:length(coord_b)
    hold on
    text(coord_b(u,1)-0.3,coord_b(u,2)-0.15,'x','FontSize',14,'Color','k')
end
offsetAxes;

box off;xticks([1:length(tem_cm)]);
%ylim([1 5.5])
hold on;xticklabels(p3);yticklabels({'L2/3','L4','L5','L6a','L6b'})
title(panel_tit{j})
%text(-1,-0.3,'L2/3 dom');
text(4,-0.3,'L2/3 dom');

if j==1 
%text(10,-0.3,'L5 dom');
text(15,-0.3,'L5 dom');
%text(31,-0.3,'L6 dom'); 
text(35,-0.3,'L6 dom'); 
elseif j==2
  %text(14,-0.3,'L5 dom');  
  %text(36,-0.3,'L6 dom'); 
  text(17,-0.3,'L5 dom');
  text(40,-0.3,'L6 dom');
else
  %text(6,-0.3,'L5 dom');  
  text(12,-0.3,'L5 dom');  
  text(28,-0.3,'L6 dom');  
end
set(gca,'TickDir','out');box off;set(gca,'FontSize',10)
% fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 500, 1400, 160]);set(gcf,'color','w');tiledlayout("horizontal")
% nexttile
% h=imagesc(tem_cm2);
% set(h, 'AlphaData', ~isnan(tem_cm2))
% cmap(temp_color(j,:),100,5,5);colorbar;box off;xticks([1:length(tem_cm)])
% hold on;xticklabels(p3);
% title(panel_tit{j})
cd(save_folder);saveas(gcf, name_sub{j});
end

%Contra
name_sub={'VISpc_heatmaps.pdf', 'SSpc_heatmaps.pdf' ,'MOpc_heatmaps.pdf'};
dat_all={};dat_all={c_v1aafm c_s1aafm c_m1aafm};
dat_all2={};dat_all2={v1a_rcam s1a_rcam m1a_rcam};
temp_color=[v1_color ;s1_color; m1_color];
panel_tit={'VISp','SSp-bfd','MOp'};
cortex_abb=[];cortex_abb=cortex_names(:,2);

for j=1:3
p1=[];p1=nanmean(dat_all{j},3);
%p2=[];p2=nanmean(dat_all3{j},3)
idx_ke=find(sum(dat_all2{j}>0,2)>1);
% idx_ke2=find(sum(dat_all4{j}>0,2)>1);
p1e=[];p1e=p1(:,idx_ke);
% p2e=[];p2e=p2(:,idx_ke);
cortex_abb_e=[];
cortex_abb_e=cortex_abb(idx_ke);
ind=[];
[X, ind]=nanmax(p1e,[],1);
l6_an=[];l5_an=[];l23_an=[];
l6_an=find(ind==5);
l5_an=find(ind==4);
l23_an=find(ind==2);
kk1=[];kk2=[];kk3=[];
[l23_p1 kk1 ]=sort(p1e(2,l23_an));
[l5_p1 kk2]=sort(p1e(4,l5_an));
[l6_p1 kk3]=sort(p1e(5,l6_an));

tem_cm=[];tem_cm=[p1e(2:end,l23_an(kk1)) ones(5,2)*1 p1e(2:end,l5_an(kk2)) ones(5,2)*1 p1e(2:end,l6_an(kk3))];
%tem_cm2=[];tem_cm2=[p2e(2:end,l23_an(kk1)) zeros(5,2)*NaN p2e(2:end,l5_an(kk2)) zeros(5,2)*NaN p2e(2:end,l6_an(kk3))];
p3=[];p3=[cortex_abb_e(l23_an(kk1)) ;'.';'.';cortex_abb_e(l5_an(kk2));'.';'.' ;cortex_abb_e(l6_an(kk3))];

fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [1300, 1000-250*j, 750, 160]);set(gcf,'color','w');tiledlayout("horizontal")
nexttile
h=imagesc(tem_cm);
set(h, 'AlphaData', ~isnan(tem_cm))
cmap(temp_color(j,:),100,0.5,0.5);colorbar;
[hh oo]=find(isnan(tem_cm));
coord_b=[oo hh];
for u=1:length(coord_b)
    hold on
    text(coord_b(u,1)-0.3,coord_b(u,2)-0.15,'x','FontSize',14,'Color','k')
end
offsetAxes;

box off;;xticks([1:length(tem_cm)]);

hold on;xticklabels(p3);xtickangle(90);yticklabels({'L2/3','L4','L5','L6a','L6b'})
title(panel_tit{j})
text(-1,-0.3,'L2/3 dom');

if j==1 
text(10,-0.3,'L5 dom');
text(31,-0.3,'L6 dom'); 

elseif j==2
  text(16,-0.3,'L5 dom');  
  text(36,-0.3,'L6 dom'); 
 
else
  text(5,-0.3,'L5 dom');  
  text(30,-0.3,'L6 dom');  
   
end
set(gca,'TickDir','out');box off;set(gca,'FontSize',10)
cd(save_folder);saveas(gcf, name_sub{j});
end
%% Plot delta ipsi-contra for VISp SSp and MOp
%calculate delta across VISp SSp and MOp
v1_d=c_v1a-i_v1a;s1_d=c_s1a-i_s1a;m1_d=c_m1a-i_m1a;
%Plot figure 
fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 500, 450, 300]);set(gcf,'color','w');t=tiledlayout("horizontal")
t.TileSpacing = 'compact';t.Padding = 'compact';
%V1
nexttile
b=bar(nanmean(v1_d(2:end,:),2),0.7);b.FaceColor=v1_color;box off;
hold on;errorbar([1 2 3 4 5], nanmean(v1_d(2:end,:),2), nanstd(v1_d(2:end,:),[],2)/(sqrt(length(v1_d))), 'LineStyle', 'none', ... 
        'Color', 'k', 'LineWidth', 1);set(b,'ShowBaseLine','off');
for r=1:5
rng(1);r1 = r-0.01 + (0.2)*rand(length(v1_d(2:end,:)),1);
sc1=scatter(r1,v1_d(r+1,:),5,'k','filled');hold on;sc1.MarkerEdgeColor=[1 1 1];sc1.MarkerEdgeAlpha=0.5;sc1.MarkerFaceColor=[0.2 0.2 0.2];
end
xticklabels({'L2/3','L4','L5','L6a','L6b'});ylabel('Fractional change (contra - ipsi)');
set(gca,'FontSize',11);set(gca,'TickDir','out');box off;ylim([-0.15 0.25]);
%title('VISp','FontWeight','normal','Color',v1_color);
hold on;text(3.75,0.235,['*'],'FontSize',18);
%hold on;text(0.75,0.235,['*'],'FontSize',18);
xtickangle(45);xlim([0.5 5.5]);offsetAxes;
h = gca;h.XAxis.Visible = 'off';  offsetAxes;
t1=text(0.9,-0.21,'L2/3','FontSize',11);set(t1,'Rotation',90);
t1=text(1.9,-0.21,'L4','FontSize',11);set(t1,'Rotation',90);t1=text(2.9,-0.21,'L5','FontSize',11);set(t1,'Rotation',90);
t1=text(3.9,-0.21,'L6a','FontSize',11);set(t1,'Rotation',90);t1=text(4.9,-0.21,'L6b','FontSize',11);set(t1,'Rotation',90);
t3=text(2,-0.24,'VISp','FontWeight','normal','Color',v1_color,'FontSize',11)

%S1
nexttile
b=bar(nanmean(s1_d(2:end,:),2),0.7);b.FaceColor=s1_color;box off;
hold on;errorbar([1 2 3 4 5], nanmean(s1_d(2:end,:),2), nanstd(s1_d(2:end,:),[],2)/(sqrt(length(s1_d))), 'LineStyle', 'none', ... 
        'Color', 'k', 'LineWidth', 1);
for r=1:5
rng(1);r1 = r-0.01 + (0.2)*rand(length(s1_d(2:end,:)),1);
sc1=scatter(r1,s1_d(r+1,:),5,'k','filled');hold on;sc1.MarkerEdgeColor=[1 1 1];sc1.MarkerEdgeAlpha=0.5;sc1.MarkerFaceColor=[0.2 0.2 0.2];
end
xticklabels({'L2/3','L4','L5','L6a','L6b'});
set(gca,'FontSize',11);set(gca,'TickDir','out');box off;
ylim([-0.15 0.25]);set(gca,'box','off','ycolor','w');
%title('SSpbf','FontWeight','normal','Color',s1_color);
hold on;text(3.75,0.235,['*'],'FontSize',18);
%hold on;text(0.75,0.235,['*'],'FontSize',18);
xtickangle(45);offsetAxes;
h = gca;h.XAxis.Visible = 'off';  offsetAxes;
t1=text(0.9,-0.21,'L2/3','FontSize',11);set(t1,'Rotation',90);
t1=text(1.9,-0.21,'L4','FontSize',11);set(t1,'Rotation',90);t1=text(2.9,-0.21,'L5','FontSize',11);set(t1,'Rotation',90);
t1=text(3.9,-0.21,'L6a','FontSize',11);set(t1,'Rotation',90);t1=text(4.9,-0.21,'L6b','FontSize',11);set(t1,'Rotation',90);
set(b,'ShowBaseLine','off');
t3=text(1.5,-0.24,'SSp-bfd','FontWeight','normal','Color',s1_color,'FontSize',11)

%M1
nexttile
b=bar(nanmean(m1_d(2:end,:),2),0.7);b.FaceColor=m1_color;box off;
hold on;errorbar([1 2 3 4 5], nanmean(m1_d(2:end,:),2), nanstd(m1_d(2:end,:),[],2)/(sqrt(length(m1_d))), 'LineStyle', 'none', ... 
        'Color', 'k', 'LineWidth', 1);set(b,'ShowBaseLine','off');
for r=1:5
rng(1);r1 = r-0.01 + (0.2)*rand(length(m1_d(2:end,:)),1);
sc1=scatter(r1,m1_d(r+1,:),5,'k','filled');hold on;sc1.MarkerEdgeColor=[1 1 1];sc1.MarkerEdgeAlpha=0.5;sc1.MarkerFaceColor=[0.2 0.2 0.2];
end
%xticklabels({'L2/3','L4','L5','L6a','L6b'});
set(gca,'FontSize',11);set(gca,'TickDir','out');box off;
ylim([-0.15 0.25]);set(gca,'box','off','ycolor','w');
%title('MOp','FontWeight','normal','Color',m1_color);
hold on;text(3.75,0.235,['*'],'FontSize',18);
%hold on;text(0.75,0.235,['*'],'FontSize',18);
xtickangle(45);h = gca;h.XAxis.Visible = 'off';  offsetAxes;
t1=text(0.9,-0.21,'L2/3','FontSize',11);set(t1,'Rotation',90);
t1=text(1.9,-0.21,'L4','FontSize',11);set(t1,'Rotation',90);t1=text(2.9,-0.21,'L5','FontSize',11);set(t1,'Rotation',90);
t1=text(3.9,-0.21,'L6a','FontSize',11);set(t1,'Rotation',90);t1=text(4.9,-0.21,'L6b','FontSize',11);set(t1,'Rotation',90);
t3=text(2,-0.225,'MOp','FontWeight','normal','Color',m1_color,'FontSize',11)
%stats
[u p1]=ttest(c_v1a(5,:),i_v1a(5,:))
[u p2]=ttest(c_s1a(5,:),i_s1a(5,:))
[u p3]=ttest(c_m1a(5,:),i_m1a(5,:))
%% 
cd(save_folder);saveas(gcf, 'region_specific_layerdiff.pdf');
%% L6a fraction index (pure index) average colour coded for V1 S1 M1
temp_color=[v1_color ;s1_color; m1_color];
module_names={'VISp','SSp-bfd','MOp'};
indx=[];indx=5;
all_dati={};all_dati={iv1_index(:,indx) is1_index(:,indx) im1_index(:,indx)};
all_datc={};all_datc={cv1_index(:,indx) cs1_index(:,indx) cm1_index(:,indx)};
%plot
fig3= figure;set(fig3, 'Name', 'Paired comp');set(fig3, 'Position', [200, 300, 170, 250]);set(gcf,'color','w');
for j=1:3
dat=[];dat=[all_dati{:,j} all_datc{:,j}];
for i=1:length(data)
     pl=plot([1,2],[dat(:,1),dat(:,2)],'color',[0.5 0.5 0.5]);    
end
hold on;pS=plotSpread([dat(:,1),dat(:,2)],'categoryIdx',[ones(1,length(dat(:,1)))' ones(1,length(dat(:,2)))'*2],...
    'categoryMarkers',{'o','o'},'categoryColors',{temp_color(j,:), temp_color(j,:)});hold on;
hold on;er1=errorbar([0.75-j*0.24],nanmean(dat(:,1)),nanstd(dat(:,1),[],1)/sqrt(length(dat(:,1))),'ok','MarkerFaceColor',temp_color(j,:),'Markersize',8);
hold on;er2=errorbar([2.25+j*0.24],nanmean(dat(:,2)),nanstd(dat(:,2),[],1)/sqrt(length(dat(:,1))),'ok','MarkerFaceColor',temp_color(j,:),'Markersize',8);
text(2.2,1.1-0.07*j,module_names{j},'Color',temp_color(j,:))
[u p1]=ttest(dat(:,1),dat(:,2))
end
xticklabels({'ipsi','contra'});ylabel('L6a / (L23+L5+L6ab)');hold on;title([]);xtickangle(45);set(gca,'FontSize',11);set(gca,'TickDir','out');box off;
hold on;text(1.5,0.6,['*'],'FontSize',18);
offsetAxes;h = gca;h.XAxis.Visible = 'off'
 t1=text(0.35,-0.1,'ipsi','FontSize',11);set(t1,'Rotation',45);
  t1=text(1.35,-0.1,'contra','FontSize',11);set(t1,'Rotation',45);

cd(save_folder);saveas(gcf, 'L6a_index.pdf');
%% Plot the numbers as fraction per area where L2/3, L5 or L6a is dominant ACROSS ALL
%all areas stacked
cat_temp=[];cat_temp=cat(3,i_v1aafm,i_s1aafm,i_m1aafm);
ind=[];
for i=1:size(cat_temp,3)
[X, ind(:,i)]=nanmax(cat_temp(:,:,i),[],1);
end
p1=[];p1=sum(ind==2)./sum(ind>1);
p2=[];p2=sum(ind==4)./sum(ind>1);
p3=[];p3=sum(ind==5)./sum(ind>1);
temp_p=[];temp_p= [p1' p2' p3'];


fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 500, 200, 270]);set(gcf,'color','w');
b1=bar([1],[nanmean(temp_p(:,1)) nanmean(temp_p(:,2)) nanmean(temp_p(:,3))],'stacked');hold on;
b1(1).FaceColor=[1 1 1];b1(2).FaceColor=[0.7 0.7 0.7];b1(3).FaceColor=[0.55 0.55 0.55];
b1(1).LineWidth=1.2;b1(2).LineWidth=1.2;b1(3).LineWidth=1.2;

hold on;errorbar([1],[nanmean(temp_p(:,1)) nanmean(temp_p(:,2))+nanmean(temp_p(:,1)) nanmean(temp_p(:,3))+nanmean(temp_p(:,2))+nanmean(temp_p(:,1))],[nanstd(temp_p(:,1))/sqrt(length(temp_p(:,1))) nanstd(temp_p(:,2))/sqrt(length(temp_p(:,2))) nanstd(temp_p(:,3))/sqrt(length(temp_p(:,3)))]...
      , 'LineStyle', 'none', ... 
         'Color', 'k', 'LineWidth', 1.2);hold on;

text(0.665,0.75,num2str(round(nanmean(temp_p(:,3)),3)*100),'FontSize',8,'Color','w');
text(0.665,0.4,num2str(round(nanmean(temp_p(:,2)),3)*100),'FontSize',8,'Color','w');
text(0.665,0.1,num2str(round(nanmean(temp_p(:,1)),3)*100),'FontSize',8);
% [p,tbl,stats] = anova1([temp_p])
%  presults = multcompare(stats)

cat_temp=[];cat_temp=cat(3,c_v1aafm,c_s1aafm,c_m1aafm);
ind=[];
for i=1:size(cat_temp,3)
[X, ind(:,i)]=nanmax(cat_temp(:,:,i),[],1);
end
p1=[];p1=sum(ind==2)./sum(ind>1);
p2=[];p2=sum(ind==4)./sum(ind>1);
p3=[];p3=sum(ind==5)./sum(ind>1);
temp_p=[];temp_p= [p1' p2' p3'];

hold on;
b1=bar([3],[nanmean(temp_p(:,1)) nanmean(temp_p(:,2)) nanmean(temp_p(:,3))],'stacked');hold on;
b1(1).FaceColor=[1 1 1];b1(2).FaceColor=[0.7 0.7 0.7];b1(3).FaceColor=[0.55 0.55 0.55];
b1(1).LineWidth=1.2;b1(2).LineWidth=1.2;b1(3).LineWidth=1.2;

hold on;errorbar([3],[nanmean(temp_p(:,1)) nanmean(temp_p(:,2))+nanmean(temp_p(:,1)) nanmean(temp_p(:,3))+nanmean(temp_p(:,2))+nanmean(temp_p(:,1))],[nanstd(temp_p(:,1))/sqrt(length(temp_p(:,1))) nanstd(temp_p(:,2))/sqrt(length(temp_p(:,2))) nanstd(temp_p(:,3))/sqrt(length(temp_p(:,3)))]...
      , 'LineStyle', 'none', ... 
         'Color', 'k', 'LineWidth', 1.2);hold on;

hold on;box off;xticks([1:2:3]);xticklabels({'ipsi','contra'});t1=text(0,0.1,{'Layer dominance ','per cortical area (%)'},'FontSize',11);set(t1,'Rotation',90);
set(gca,'FontSize',11);set(gca,'TickDir','out');
h = gca;h.YAxis.Visible = 'off';  

text(2.665,0.75,num2str(round(nanmean(temp_p(:,3)),3)*100),'FontSize',8,'Color','w');
text(2.665,0.3,num2str(round(nanmean(temp_p(:,2)),3)*100),'FontSize',8,'Color','w');
text(2.665,0.04,num2str(round(nanmean(temp_p(:,1)),3)*100),'FontSize',8);


t1=text(3.7,0.8,'L6a','FontSize',11);set(t1,'Rotation',270);
t1=text(3.7,0.35,'L5','FontSize',11);set(t1,'Rotation',270);
t1=text(3.7,0.1,'L2/3','FontSize',11);set(t1,'Rotation',270);

%Overall
% hold on;text(1.8,1.15,'All','FontSize',12);
% ax = gca; ax.XColor = 'w'; % Red
% hold on;line([0.55 1.45],[0 0],'Color','k','LineWidth', 1);
% hold on;line([2.55 3.45],[0 0],'Color','k','LineWidth', 1);
h = gca;h.XAxis.Visible = 'off';set(b1,'ShowBaseLine','off');
hold on;text(0.65,-0.065,'ipsi','FontSize',11);
hold on;text(2.4,-0.065,'contra','FontSize',11);
hold on;text(1.88,0.4,'*','FontSize',18);
line([2.3 2.6],[0.7 0.7],'Color','k');
line([2.3 2.6],[0.35 0.35],'Color','k');
line([2.3 2.6],[0.05 0.05],'Color','k');
line([2.3 2.3],[0.05 0.7],'Color','k');

line([1.7 1.4],[0.65 0.65],'Color','k');
line([1.7 1.4],[0.3 0.3],'Color','k');
line([1.7 1.4],[0.1 0.1],'Color','k');
line([1.7 1.7],[0.1 0.65],'Color','k');

cd(save_folder);saveas(gcf,'overall_dominance_2.pdf');

% [p,tbl,stats] = anova1([temp_p])
%  presults = multcompare(stats)
%% Based on area VISp, SSpbf, MOp dominace of L2/3 L5 and L6a CONTRA 
dat_all={};dat_all={c_v1aafm c_s1aafm c_m1aafm};
temp_color=[v1_color ;s1_color; m1_color];
panel_tit={'','contralateral',''};
all_ind={};
%plot
fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 350, 350, 250]);set(gcf,'color','w');t=tiledlayout("horizontal");
t.TileSpacing = 'compact';t.Padding = 'compact';
for j=1:3
    figure(fig7)
ind=[];temp_c=temp_color(j,:);temp_area=[];temp_area=dat_all{j};
    for i=1:size(temp_area,3)
    [X, ind(:,i)]=nanmax(temp_area(:,:,i),[],1);
    end
    all_ind{j}=ind;
   
    p1=[];p1=(sum(ind==2)./sum(ind>1))*100;
    p2=[];p2=(sum(ind==4)./sum(ind>1))*100;
    p3=[];p3=(sum(ind==5)./sum(ind>1))*100;
    temp_p= [p1' p2' p3'];
nexttile
title(panel_tit{j},'FontWeight','normal');
    for i=1:3
    hold on;
    b1=bar(i,nanmean(temp_p(:,i)),0.7);hold on;b1.FaceColor=temp_c;b1.EdgeColor='none';
      r=i;rng(i);r1 = r-0.01 + (0.2)*rand(length(temp_p(:,i)),1);
    sc1=scatter(r1,temp_p(:,i),12,'k','filled');hold on;sc1.MarkerEdgeColor=[1 1 1];sc1.MarkerEdgeAlpha=0.5;sc1.MarkerFaceColor=[0.2 0.2 0.2];
    hold on;errorbar([i],[nanmean(temp_p(:,i))],[nanstd(temp_p(:,i))/sqrt(length(temp_p(:,i)))]...
    , 'LineStyle', 'none', ... 
        'Color', 'k', 'LineWidth', 1);hold on;set(b1,'ShowBaseLine','off');
    end
    if j==1
xticks([1:3]);hold on;xticklabels({'L2/3','L5','L6a'});ylabel('Layer dominance per area (%)');
set(gca,'FontSize',11);set(gca,'TickDir','out');box off;
    else
        h = gca;h.YAxis.Visible = 'off';xticks([1:3]);hold on;xticklabels({'L2/3','L5','L6a'});
        set(gca,'FontSize',11);set(gca,'TickDir','out');box off;xtickangle(45);
    end

xlim([0.5 3.5]);
offsetAxes;h = gca;h.XAxis.Visible = 'off'
t1=text(0.9,-17,'L2/3','FontSize',11);set(t1,'Rotation',90);
t1=text(1.9,-17,'L5','FontSize',11);set(t1,'Rotation',90);
t1=text(2.9,-17,'L6a','FontSize',11);set(t1,'Rotation',90);
ylim([0 80])
[p,tbl,stats] = anova1([temp_p])
 presults = multcompare(stats)
end
cd(save_folder);saveas(gcf, 'dominance_per_injection.pdf');
%% IPSI
dat_all={};dat_all={i_v1aafm i_s1aafm i_m1aafm};
temp_color=[v1_color ;s1_color; m1_color];
panel_tit={'','ipsilateral',''};
all_ind={};
%plot
fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 350, 350, 250]);set(gcf,'color','w');t=tiledlayout("horizontal");
t.TileSpacing = 'compact';t.Padding = 'compact';
for j=1:3
    figure(fig7)
ind=[];temp_c=temp_color(j,:);temp_area=[];temp_area=dat_all{j};
    for i=1:size(temp_area,3)
    [X, ind(:,i)]=nanmax(temp_area(:,:,i),[],1);
    end
    all_ind{j}=ind;
   
    p1=[];p1=(sum(ind==2)./sum(ind>1))*100;
    p2=[];p2=(sum(ind==4)./sum(ind>1))*100;
    p3=[];p3=(sum(ind==5)./sum(ind>1))*100;
    temp_p= [p1' p2' p3'];
nexttile
title(panel_tit{j},'FontWeight','normal');
    for i=1:3
    hold on;
    b1=bar(i,nanmean(temp_p(:,i)),0.7);hold on;b1.FaceColor=temp_c;b1.EdgeColor='none';
      r=i;rng(i);r1 = r-0.01 + (0.2)*rand(length(temp_p(:,i)),1);
    sc1=scatter(r1,temp_p(:,i),12,'k','filled');hold on;sc1.MarkerEdgeColor=[1 1 1];sc1.MarkerEdgeAlpha=0.5;sc1.MarkerFaceColor=[0.2 0.2 0.2];
    hold on;errorbar([i],[nanmean(temp_p(:,i))],[nanstd(temp_p(:,i))/sqrt(length(temp_p(:,i)))]...
    , 'LineStyle', 'none', ... 
        'Color', 'k', 'LineWidth', 1);hold on;set(b1,'ShowBaseLine','off');
    end
    if j==1
xticks([1:3]);hold on;xticklabels({'L2/3','L5','L6a'});ylabel('Layer dominance per area (%)');
set(gca,'FontSize',11);set(gca,'TickDir','out');box off;
    else
        h = gca;h.YAxis.Visible = 'off';xticks([1:3]);hold on;xticklabels({'L2/3','L5','L6a'});
        set(gca,'FontSize',11);set(gca,'TickDir','out');box off;xtickangle(45);
    end

xlim([0.5 3.5]);
offsetAxes;h = gca;h.XAxis.Visible = 'off'
t1=text(0.9,-17,'L2/3','FontSize',11);set(t1,'Rotation',90);
t1=text(1.9,-17,'L5','FontSize',11);set(t1,'Rotation',90);
t1=text(2.9,-17,'L6a','FontSize',11);set(t1,'Rotation',90);
ylim([0 80])

[p,tbl,stats] = anova1([temp_p])
 presults = multcompare(stats)
end

%cd(save_folder);saveas(gcf, 'dominance_per_injection.pdf');
%% Figure 4 HIEARCHY
%% Next to each other ILN ipsi vs contra
idx=1;
temp_color=[v1_color ;s1_color; m1_color];
module_names={'VISp','SSp-bfd','MOp'};
all_dati={};all_dati={iv1_index(:,idx) is1_index(:,idx) im1_index(:,idx)};
all_datc={};all_datc={cv1_index(:,idx) cs1_index(:,idx) cm1_index(:,idx)};
%plot
fig3= figure;set(fig3, 'Name', 'Paired comp');set(fig3, 'Position', [200, 300, 400, 250]);set(gcf,'color','w');t=tiledlayout("horizontal");
t.TileSpacing = 'compact';t.Padding = 'compact';
for j=1:3
    nexttile
 
dat=[];dat=[all_dati{:,j} all_datc{:,j}];
for i=1:length(data)
     pl=plot([1,2],[dat(:,1),dat(:,2)],'color',[0.5 0.5 0.5]);    
end

hold on;pS=plotSpread([dat(:,1),dat(:,2)],'categoryIdx',[ones(1,length(dat(:,1)))' ones(1,length(dat(:,2)))'*2],...
    'categoryMarkers',{'o','o'},'categoryColors',{temp_color(j,:), temp_color(j,:)});hold on;
hold on;er1=errorbar([0.75],nanmean(dat(:,1)),nanstd(dat(:,1),[],1)/sqrt(length(dat(:,1))),'ok','MarkerFaceColor',temp_color(j,:),'Markersize',8);
hold on;er2=errorbar([2.25],nanmean(dat(:,2)),nanstd(dat(:,2),[],1)/sqrt(length(dat(:,1))),'ok','MarkerFaceColor',temp_color(j,:),'Markersize',8);
%stats
[u p1]=ttest(dat(:,1),dat(:,2))

offsetAxes
xticklabels({'ipsi','contra'});ylabel('ILN');hold on;title([]);xtickangle(45);set(gca,'FontSize',11);set(gca,'TickDir','out');box off;
hold on;text(1.5,0.9,['*'],'FontSize',18);ylim([0.4 0.9])
 h = gca;h.XAxis.Visible = 'off'
 t1=text(0.75,0.4,'ipsi','FontSize',11);set(t1,'Rotation',45);
  t1=text(1.75,0.4,'contra','FontSize',11);set(t1,'Rotation',45);
  line([-0.1 3],[0.5 0.5],'Color','k','LineStyle','--')
  if j==2 | j==3
     h = gca;h.YAxis.Visible = 'off'  
  end
end
cd(save_folder);saveas(gcf, 'ILN_global.pdf');
%% All animals all injections together 
idx=1;
all_dati1=[];all_dati1=[iv1_index(:,idx) ;is1_index(:,idx);im1_index(:,idx)];
all_datc1=[];all_datc1=[cv1_index(:,idx) ;cs1_index(:,idx);cm1_index(:,idx)];
idx=7;
all_dati2=[];all_dati2=[iv1_index(:,idx) ;is1_index(:,idx);im1_index(:,idx)];
all_datc2=[];all_datc2=[cv1_index(:,idx) ;cs1_index(:,idx);cm1_index(:,idx)];
%% 


fig3= figure;set(fig3, 'Name', 'Paired comp');set(fig3, 'Position', [200, 300, 200, 250]);set(gcf,'color','w');

dat=[];dat=[all_dati all_datc];
for i=1:length(data)
     pl=plot([1,2],[dat(:,1),dat(:,2)],'color',[0.5 0.5 0.5]);    
end

hold on;pS=plotSpread([dat(:,1),dat(:,2)],'categoryIdx',[ones(1,length(dat(:,1)))' ones(1,length(dat(:,2)))'*2],...
    'categoryMarkers',{'o','o'},'categoryColors',{temp_color(j,:), temp_color(j,:)});hold on;
hold on;er1=errorbar([0.75],nanmean(dat(:,1)),nanstd(dat(:,1),[],1)/sqrt(length(dat(:,1))),'ok','MarkerFaceColor',temp_color(j,:),'Markersize',8);
hold on;er2=errorbar([2.25],nanmean(dat(:,2)),nanstd(dat(:,2),[],1)/sqrt(length(dat(:,1))),'ok','MarkerFaceColor',temp_color(j,:),'Markersize',8);
%stats
[u p1]=ttest(dat(:,1),dat(:,2))

offsetAxes
xticklabels({'ipsi','contra'});ylabel('ILN');hold on;title([]);xtickangle(45);set(gca,'FontSize',11);set(gca,'TickDir','out');box off;
hold on;text(1.5,0.9,['*'],'FontSize',18);ylim([0.4 0.9])
 h = gca;h.XAxis.Visible = 'off'
 t1=text(0.75,0.4,'ipsi','FontSize',11);set(t1,'Rotation',45);
  t1=text(1.75,0.4,'contra','FontSize',11);set(t1,'Rotation',45);
  line([-0.1 3],[0.5 0.5],'Color','k','LineStyle','--')
%% Sorted areas based on index SCHEREN PLOT
dat_all1={};dat_all1={i_v1aafm i_s1aafm i_m1aafm};
dat_all2={};dat_all2={c_v1aafm c_s1aafm c_m1aafm};
title_pan={'VISp','SSp-bfd','MOp'};
cl_idx=[visp_idx ssp_idx mop_idx];
temp_c1=[0.8 0.8 0.8];
temp_c2=[0.3 0.3 0.3];
temp_color=[v1_color ;s1_color; m1_color];
idx=[];idx=1;
aILNi=[];
aILNc=[];


for k=1:3
    rm1=[];rm2=[];rm1=dat_all1{k};rm2=dat_all2{k};
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
aILNi=[aILNi; rm1];
aILNc=[aILNc; rm2];

temp_rm1=[];temp_rm1=nanmean(rm1);
temp_rm2=[];temp_rm2=nanmean(rm2);
temp_rm1(cl_idx(k))=[];
temp_rm2(cl_idx(k))=[];
cortex_abb=cortex_names(:,2);
cortex_abb(cl_idx(k))=[];
yerr=[];yerr=nanstd(rm1);
yerr2=[];yerr2=nanstd(rm2);
yerr3=[];yerr3=nanstd(rm2(:,cl_idx(k)))/sqrt(length(rm2(:,cl_idx(k))));
yerr(cl_idx(k))=[];
yerr2(cl_idx(k))=[];
sort_d=[];kk=[];
[sort_d kk]=sort(temp_rm1,'ascend');

alpha1=1;
fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, -200+k*300, 850, 250]);set(gcf,'color','w');

title(title_pan{k},'FontWeight','normal','Color',temp_color(k,:));
ylabel('ILN')
hold on;h=errorbar([1:44],[sort_d],[yerr/sqrt(length(sort_d(~isnan(sort_d))))]...
    , 'LineStyle', 'none', ... 
        'Color', temp_c1, 'LineWidth', 0.5,'CapSize',2);set([h.Bar, h.Line], 'ColorType', 'truecoloralpha', 'ColorData', [h.Line.ColorData(1:3); 255*alpha1]);
hold on;h=errorbar([1:44],[temp_rm2(kk)],[yerr2(kk)/sqrt(length(temp_rm2(~isnan(temp_rm2))))]...
    , 'LineStyle', 'none', ... 
        'Color', temp_c2, 'LineWidth', 0.5,'CapSize',2);
hold on;h=errorbar([45],[nanmean(temp_metric2(:,cl_idx(k)))],[yerr3]...
    , 'LineStyle', 'none', ... 
        'Color', temp_color(k,:), 'LineWidth', 0.5,'CapSize',2);
hold on;pp1=scatter(1:44,sort_d,30,'o');pp1.MarkerEdgeColor=[1 1 1];pp1.MarkerFaceAlpha=1;pp1.MarkerFaceColor=temp_c1;
hold on;pp1=scatter(1:44,temp_rm2(kk),30,'o');pp1.MarkerEdgeColor=[1 1 1];pp1.MarkerFaceAlpha=1;pp1.MarkerFaceColor=temp_c2;box off;
hold on;pp1=scatter(45,nanmean(temp_metric2(:,cl_idx(k))),30,'o');pp1.MarkerEdgeColor=[1 1 1];pp1.MarkerFaceAlpha=1;pp1.MarkerFaceColor=temp_color(k,:);box off;
xticks([1:45]);ylim([min(get(gca,'YLim')) 1.2])
xticklabels([cortex_abb(kk) ;{title_pan{k}}])
set(gca,'FontSize',11);set(gca,'TickDir','out');
if k==1
cd(save_folder);saveas(gcf, 'ILNipsi_sorted_areasVISP.pdf');
end
if k==2
    cd(save_folder);saveas(gcf, 'ILNipsi_sorted_areasSSP.pdf');
end
if k==3
    cd(save_folder);saveas(gcf, 'ILNipsi_sorted_areasMOP.pdf');
end
end
%% 

% SCHEREN PLOT ALL
temp_rm1=[];temp_rm1=nanmean(aILNi);
temp_rm2=[];temp_rm2=nanmean(aILNc);
temp_c1=[0.8 0.8 0.8];
temp_c2=[0.3 0.3 0.3];

yerr=[];yerr=nanstd(aILNi);
yerr2=[];yerr2=nanstd(aILNc);
%yerr3=[];yerr3=nanstd(rm2(:,cl_idx(k)))/sqrt(length(rm2(:,cl_idx(k))));
sort_d=[];kk=[];
[sort_d kk]=sort(temp_rm1,'ascend');

alpha1=1;
fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 300, 600, 300]);set(gcf,'color','w');

title('sorted by ipsi ILN','FontWeight','normal','Color','k');
ylabel('ILN')
li=line([1 45],[0.5 0.5],'Color',[0.5 0.5 0.5],'LineStyle',':');
hold on;plot(sort_d,'Color',temp_c1);
hold on;plot(temp_rm2(kk),'Color',temp_c2);
hold on;h=errorbar([1:45],[sort_d],[yerr/sqrt(length(sort_d(~isnan(sort_d))))]...
    , 'LineStyle', 'none', ... 
        'Color', temp_c1, 'LineWidth', 0.5,'CapSize',2);set([h.Bar, h.Line], 'ColorType', 'truecoloralpha', 'ColorData', [h.Line.ColorData(1:3); 255*alpha1]);
hold on;h=errorbar([1:45],[temp_rm2(kk)],[yerr2(kk)/sqrt(length(temp_rm2(~isnan(temp_rm2))))]...
    , 'LineStyle', 'none', ... 
        'Color', temp_c2, 'LineWidth', 0.5,'CapSize',2);

hold on;pp1=scatter(1:45,sort_d,30,'o');pp1.MarkerEdgeColor=[1 1 1];pp1.MarkerFaceAlpha=1;pp1.MarkerFaceColor=temp_c1;
hold on;pp1=scatter(1:45,temp_rm2(kk),30,'o');pp1.MarkerEdgeColor=[1 1 1];pp1.MarkerFaceAlpha=1;pp1.MarkerFaceColor=temp_c2;box off;

xticks([1:45]);ylim([0 1])
ylabel('ILN');xlabel('Cortical areas')
set(gca,'FontSize',11);set(gca,'TickDir','out');
text(43,0.7,'ipsi','Color',[0.7 0.7 0.7]);
text(43,0.6,'contra','Color',[0.3 0.3 0.3]);

offsetAxes;

     h = gca;h.XAxis.Visible = 'off'  
     text(20,-0.1,'Cortical areas','FontSize',11);
     cd(save_folder);saveas(gcf, 'ILNipsi_sorted_areasALL.pdf');
  
  %% Histogram across all areas all injections 
bwi=0.05;
dat_all1={};dat_all1={i_v1aafm i_s1aafm i_m1aafm};
dat_all2={};dat_all2={c_v1aafm c_s1aafm c_m1aafm};
idx_modules={frontal_idx lateral_idx somamo_idx visual_idx medial_idx aud_idx};
module_names={'Fron','Lat','SoMo','Vis','Med','Aud'};
pan_title={'VISp','SSp-bfd','MOp'};
temp_color=[v1_color ;s1_color; m1_color];
temp_c1=[0.8 0.8 0.8];
temp_c2=[0.3 0.3 0.3];
idx=[];idx=1;

fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 500, 400, 250]);set(gcf,'color','w');
t1_all=[];t2_all=[];
for j=1:3
rm1=[];rm2=[];rm1=dat_all1{j};rm2=dat_all2{j};
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
[t1_all]=[t1_all;temp_metric1(:)];
[t2_all]=[t2_all;temp_metric2(:)];
end
hold on;
h1=histogram(t1_all,20,'Normalization','probability');hold on;h1.EdgeColor=temp_c1;h1.FaceColor='w';h1.FaceAlpha=0.2;h1.LineWidth=1.5;h1.BinWidth = bwi;
hold on;h2=histogram(t2_all,20,'Normalization','probability');hold on;h2.EdgeColor=temp_c2;h2.FaceColor='w';h2.FaceAlpha=0.2;h2.LineWidth=1.5;h2.BinWidth = bwi;
ylabel('Fraction of areas');xlabel('ILN');set(gca,'FontSize',11);set(gca,'TickDir','out');
text(0,0.25,'ipsi','FontSize',11,'Color',temp_c1);text(0,0.225,'contra','FontSize',11,'Color',temp_c2);
hold on;line([0.5 0.5],[0 0.3],'Color','k','LineStyle','--')
offsetAxes
%% %% Histogram ranked
bwi=0.05;
dat_all1={};dat_all1={i_v1aafm i_s1aafm i_m1aafm};
dat_all2={};dat_all2={c_v1aafm c_s1aafm c_m1aafm};
idx_modules={frontal_idx lateral_idx somamo_idx visual_idx medial_idx aud_idx};
module_names={'Fron','Lat','SoMo','Vis','Med','Aud'};
pan_title={'VISp','SSp-bfd','MOp'};
temp_color=[v1_color ;s1_color; m1_color];
temp_c1=[0.65 0.65 0.65];
temp_c2=[0.3 0.3 0.3];
idx=[];idx=1;
kk_all=[];
fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 500, 500, 300]);set(gcf,'color','w');
t1_all=[];t2_all=[];
for j=1:3
rm1=[];rm2=[];rm1=dat_all1{j};rm2=dat_all2{j};
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
kk_a=[];
for m=1:size(temp_metric2,1)
 sm1=[];sm2=[]; mm=[];oo=[];kk=[];ss=[];ss2=[];
 sm2=temp_metric2(m,:);sm1=temp_metric1(m,:);
 [mm oo] = sort(sm2);
%  [ss kk]=discretize(sm2(oo),45);
% [ss2 kk2]=discretize(sm1(oo),kk);
% kk_a(m,:)=kk;
t1a(m,:)=sm1(oo);
t2a(m,:)=sm2(oo);
end
%kk_all(:,:,j)=kk_a;
% h1=histogram(sm1(oo),kk);hold on;h2=histogram(sm2(oo),kk);
% hold on;h1.EdgeColor=temp_c1;h1.FaceColor='w';h1.FaceAlpha=0.2;h1.LineWidth=1.5;h1.BinWidth = bwi;
% hold on;h2.EdgeColor=temp_c2;h2.FaceColor='w';h2.FaceAlpha=0.2;h2.LineWidth=1.5;h2.BinWidth = bwi;
% t1_all(:,:,j)=t1a;
% t2_all(:,:,j)=t2a;

 [t1_all]=[t1_all;t1a];
 [t2_all]=[t2_all;t2a];
end
yerr2=[];yerr2=nanstd(t2_all);
yerr1=[];yerr1=nanstd(t1_all);

 hold on;h=errorbar([1:45],[nanmean(t1_all)],[yerr1./sqrt(sum(~isnan(t1_all)))]...
    , 'LineStyle', 'none', ... 
        'Color', temp_c1, 'LineWidth', 0.5,'CapSize',2);
  hold on;h=errorbar([1:45],[nanmean(t2_all)],[yerr2./sqrt(sum(~isnan(t2_all)))]...
    , 'LineStyle', 'none', ... 
        'Color', temp_c2, 'LineWidth', 0.5,'CapSize',2);

hold on;pp1=scatter(1:45,nanmean(t1_all),20,'o');pp1.MarkerEdgeColor=[1 1 1];pp1.MarkerFaceAlpha=1;pp1.MarkerFaceColor=temp_c1;
hold on;pp1=scatter(1:45,nanmean(t2_all),20,'o');pp1.MarkerEdgeColor=[1 1 1];pp1.MarkerFaceAlpha=1;pp1.MarkerFaceColor=temp_c2;box off;
xticks([1:1:45]);ylim([min(get(gca,'YLim')) 1]);xlim([-1 42]);yticks([0.3:0.2:1])
ylabel('ILN');xlabel('Cortical areas');
set(gca,'FontSize',11);set(gca,'TickDir','out');xtickangle(45);
 h = gca;h.XAxis.Visible = 'off';
%legend({'ipsi','contra'});
text(32,0.55,'contra ranked','Color',temp_c2,'FontSize',11);
text(32,0.48,'ipsi','Color',temp_c1,'FontSize',11);
text(13,0.26,'Cortical brain areas','FontSize',11);
offsetAxes

cd(save_folder);saveas(gcf, 'ILN_contraRANKED.pdf');
%% 

% edges_avg=nanmean(nanmean(kk_all),3);
% 
% bwi=0.05;
% dat_all1={};dat_all1={i_v1aam i_s1aam i_m1aam};
% dat_all2={};dat_all2={c_v1aam c_s1aam c_m1aam};
% idx_modules={frontal_idx lateral_idx somamo_idx visual_idx medial_idx aud_idx};
% module_names={'Fron','Lat','SoMo','Vis','Med','Aud'};
% pan_title={'VISp','SSp-bfd','MOp'};
% temp_color=[v1_color ;s1_color; m1_color];
% temp_c1=[0.8 0.8 0.8];
% temp_c2=[0.3 0.3 0.3];
% idx=[];idx=1;
% 
% fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 500, 400, 250]);set(gcf,'color','w');
% t1_all=[];t2_all=[];
% for j=1:3
% rm1=[];rm2=[];rm1=dat_all1{j};rm2=dat_all2{j};
% %calculate and plot 
% i_metric=[];c_metric=[];
% %calculate indexes across areas
% for i=1:length(rm1)
% temp_l=[];temp_l2=[];
% temp_l=squeeze(rm1(:,i,:));
% temp_l2=squeeze(rm2(:,i,:));
% [i_metric(:,i,:)] = anatomy_indexcalc(temp_l);
% [c_metric(:,i,:)] = anatomy_indexcalc(temp_l2);
% end 
% temp_metric1=i_metric(:,:,idx);temp_metric2=c_metric(:,:,idx);
% for m=1:size(temp_metric2,1)
%  sm1=[];sm2=[]; mm=[];oo=[];kk=[];ss=[];ss2=[];
%  sm2=temp_metric2(m,:);sm1=temp_metric1(m,:);
%  [mm oo] = sort(sm2);
%  [ss kk]=discretize(sm2(oo),45);
% [ss2 kk2]=discretize(sm1(oo),edges_avg);
% 
% % kk_a(m,:)=kk;
%  t1a(m,:)=sm1(oo);
%  t2a(m,:)=sm2(oo);
% end
% kk_all(:,:,j)=kk_a;
% 
% % t1_all(:,:,j)=t1a;
% % t2_all(:,:,j)=t2a;
% 
%  [t1_all]=[t1_all;t1a(:)];
%  [t2_all]=[t2_all;t2a(:)];
% end
% 
%  h1=histogram(t1_all,edges_avg);hold on;h2=histogram(t2_all,edges_avg);
%  hold on;h1.EdgeColor=temp_c1;h1.FaceColor='w';h1.FaceAlpha=0.2;h1.LineWidth=1.5;
%  hold on;h2.EdgeColor=temp_c2;h2.FaceColor='w';h2.FaceAlpha=0.2;h2.LineWidth=1.5;
% offsetAxes
%% Figure 4

%% Modules showing ILNi and ILNc
dat_all1={};dat_all1={i_v1aafm i_s1aafm i_m1aafm};
dat_all2={};dat_all2={c_v1aafm c_s1aafm c_m1aafm};
%idx_modules={frontal_idx lateral_idx somamo_idx visual_idx medial_idx aud_idx};
idx_modules={somamo_idx visual_idx aud_idx frontal_idx lateral_idx medial_idx };
%module_names={'Pfron','Lat','SoMo','Vis','Med','Aud'};
module_names={'SoMo','Vis','Aud','Pfron','Lat','Med',};
pan_title={'VISp','SSp-bfd','MOp'};
temp_color=[v1_color ;s1_color; m1_color];
temp_c1=[0.8 0.8 0.8];
temp_c2=[0.3 0.3 0.3];
idx=[];idx=1;
p_fig3f=[];
%plot
fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 500, 250, 330]);set(gcf,'color','w');t=tiledlayout('vertical');
t.TileSpacing = 'compact';t.Padding = 'compact';
for j=1:3
rm1=[];rm2=[];rm1=dat_all1{j};rm2=dat_all2{j};
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
title(pan_title{j},'FontWeight','normal','FontSize',7,'Color',temp_color(j,:));
for i=1:6
    hold on;
    yerr1=[];yerr1=nanstd(nanmean(temp_metric1(:,idx_modules{i}),2))/sqrt(length(nanstd(temp_metric1(:,idx_modules{i}),[],2)));
    yerr2=[];yerr2=nanstd(nanmean(temp_metric2(:,idx_modules{i}),2))/sqrt(length(nanstd(temp_metric2(:,idx_modules{i}),[],2)));
    ymean1=[];ymean1=nanmean(temp_metric1(:,idx_modules{i}),2);
    ymean2=[];ymean2=nanmean(temp_metric2(:,idx_modules{i}),2);
    % hold on;h=errorbar([i],[nanmean(ymean2)],[yerr2], 'LineStyle', 'none','Color', temp_c1, 'LineWidth', 0.5,'CapSize',0);
    % hold on;h=errorbar([i],[nanmean(ymean1)],[yerr1], 'LineStyle', 'none','Color', temp_c1, 'LineWidth', 1,'CapSize',0);
    ymean_all=[];ymean_all=[nanmean(ymean1) nanmean(ymean2)];
    b1=bar(i,ymean_all,0.7);hold on;b1(2).FaceColor=temp_c2;b1(1).FaceColor=temp_c1;set(b1,'ShowBaseLine','off');b1(1).EdgeColor='none'
    hold on;h=errorbar([i-0.15],[nanmean(ymean1)],[yerr1], 'LineStyle', 'none','Color', 'k', 'LineWidth', 1,'CapSize',0);
    hold on;h=errorbar([i+0.15],[nanmean(ymean2)],[yerr2], 'LineStyle', 'none','Color', 'k', 'LineWidth', 1,'CapSize',0);
    % pp1=scatter(i,nanmean(ymean1),30,'filled');pp1.MarkerEdgeColor=[1 1 1];pp1.MarkerFaceAlpha=1;pp1.MarkerFaceColor=temp_c1;
    % pp1=scatter(i,nanmean(ymean2),30,'filled');pp1.MarkerEdgeColor=[1 1 1];pp1.MarkerFaceAlpha=1;pp1.MarkerFaceColor=temp_c2;box off;
    box off;
    [u p_fig3f(:,i,j)]=ttest(ymean1,ymean2)
    if p_fig3f(:,i,j)<0.05
        text(i,1,'*','FontSize',12);
    end
end
if j==1 
xlim([0.5 6.5]);
 text(5.5,1.25,'ipsi','Color',[0.7 0.7 0.7]);text(5.5,1.15,'contra','Color',temp_c2);
set(gca,'FontSize',11);set(gca,'TickDir','out');
h = gca;h.XAxis.Visible = 'off';
yticks([0.5:0.2:1]);
elseif j==2
    set(gca,'FontSize',11);set(gca,'TickDir','out');
    xlim([0.5 6.5]);;h = gca;h.XAxis.Visible = 'off'; ylabel('ILN');yticks([0.5:0.2:1])
else
   xlim([0.5 6.5]);
    xticks([1:6]);xticklabels(module_names);xtickangle(45);yticks([0.5:0.2:1])
    set(gca,'FontSize',11);set(gca,'TickDir','out');   
end
offsetAxes
ylim([0.5 1])
end
cd(save_folder);saveas(gcf, 'ILN_modules.pdf');
%% ILN ipsi and ILN contra across all animals 
dat_all1={};dat_all1={i_v1aafm i_s1aafm i_m1aafm};
dat_all2={};dat_all2={c_v1aafm c_s1aafm c_m1aafm};
senso_mo=[idx_modules{3} idx_modules{4} idx_modules{6}];
higher_a=[idx_modules{1} idx_modules{2} idx_modules{5}];
idx_modules={frontal_idx lateral_idx somamo_idx visual_idx medial_idx aud_idx};
module_names={'Pfron','Lat','SoMo','Vis','Med','Aud'};
pan_title={'VISp','SSp-bfd','MOp'};
temp_color=[v1_color ;s1_color; m1_color];
temp_c1='r';
temp_c2=[0.3 0.3 0.3];
idx=[];idx=1;

all_metric1=[];all_metric2=[];
fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 500, 250, 350]);set(gcf,'color','w');
%plot
for j=1:3
rm1=[];rm2=[];rm1=dat_all1{j};rm2=dat_all2{j};
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
all_metric1=[all_metric1; temp_metric1];
all_metric2=[all_metric2; temp_metric2];
end

    hold on;
for i=1:6
    ymean1=[];ymean1=nanmean(all_metric1(:,idx_modules{i}));
    ymean2=[];ymean2=nanmean(all_metric2(:,idx_modules{i}));
     pp1=scatter(1,nanmean(ymean1),80,'filled');     pp1.MarkerEdgeColor=[1 1 1]
     if i==1 | i==2 | i==5
pp1.MarkerFaceColor=[0.7 0.7 0.7];

     else
         pp1.MarkerFaceColor='r';
     end
     if i==2
         hold on;text(1.25,nanmean(ymean1)+0.015,module_names{i});
     elseif i==1
          hold on;text(1.25,nanmean(ymean1)+0.006,module_names{i});
     else
hold on;text(1.25,nanmean(ymean1),module_names{i});
     end
     
        pp1=scatter(3,nanmean(ymean2),80,'filled');pp1.MarkerEdgeColor=[1 1 1];
         if i==1 | i==2 | i==5 
        pp1.MarkerFaceColor=temp_c2;
        pp1.MarkerFaceColor=[0.7 0.7 0.7];
     else
         pp1.MarkerFaceColor='r';
         end
         if i==2
    hold on;text(3.25,nanmean(ymean2)+0.015,module_names{i});
      elseif i==1
          hold on;text(3.25,nanmean(ymean2)-0.006,module_names{i});
         else
hold on;text(3.25,nanmean(ymean2),module_names{i});
         end
    xlim([0 6]);
    ylim([0.55 0.95]);
end
    box off;
    ylabel('Cortical hierarchy (ILN)');set(gca,'FontSize',11);set(gca,'TickDir','out');
    yticks([0.55:0.1:0.95]);xticks([1:2:3]);xticklabels({'ipsi','contra'})
    offsetAxes
h = gca;h.XAxis.Visible = 'off';
t1=text(0.7,0.5,'ipsi','FontSize',11,'Color',[0.7 0.7 0.7]);
 set(t1,'Rotation',45);
 t1=text(2.8,0.5,'contra','FontSize',11,'Color',[0.3 0.3 0.3]);
 set(t1,'Rotation',45); 
line([4.9 4.9],[0.89 0.95],'Color',[0.5 0.5 0.5]);
line([4.9 4.9],[0.75 0.82],'Color','r');
line([4.9 5.3],[0.7850 0.7850],'Color','r');
line([4.9 5.3],[0.92 0.92],'Color',[0.5 0.5 0.5]);
text(5.35,0.7850,'v-sm-a','Color','r');
text(5.35,0.92,'pf-m-l','Color',[0.5 0.5 0.5]);
cd(save_folder);saveas(gcf, 'ILN_modules_hierarchy.pdf');
%% correlation between indexes per area  ILN
%% all correlations superimposed on each other 
idx=[];idx=1;
temp_all={i_v1aafm c_v1aafm; i_s1aafm c_s1aafm; i_m1aafm c_m1aafm};
temp_color=[v1_color ;s1_color; m1_color];
%plot
fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 500, 350, 350]);set(gcf,'color','w');
for i=1:3
  
    rm1=[];rm1=temp_all{i,1};rm2=[];rm2=temp_all{i,2};
[vali_all valdci_all val_i val_dci xerr yerr r_a p_a r_avg p_avg] = anatomy_correlation(rm1, rm2, idx);
 temp_c=temp_color(i,:)
 col_ra2(:,i)=r_a;
hold on;errorbar(val_i,val_dci,yerr,yerr,xerr,xerr, 'LineStyle', 'none', 'Color', 'k', 'LineWidth', 0.5,'CapSize',0);hold on;
sc1=scatter(val_i,val_dci,50,'filled');sc1.MarkerFaceColor=temp_c;sc1.MarkerEdgeColor='w';
hold on;
hold on;text(0.65,0.85-(0.08*i),...
    ['r= ' num2str(round(nanmean(r_a),2)) ' +- ' num2str(round(nanstd(r_a)/sqrt(length(r_a)),2))],'FontSize',11,'Color',temp_c);
%hold on;text(0.5,0.5-(0.05*i),['p= ' num2str(round(nanmean(p_a),2))],'FontSize',11,'Color',temp_c);
end
ylim([-0.4 0.8])
hold on;line([0.5 0.5],[-0.4 0.8],'Color','k','LineStyle','--')
xlim([0.15 1])
%ylabel('L6d_c - L6d_i');xlabel('L6d_i');
ylabel('ILN contra - ILN ipsi');xlabel('ILN ipsi');
axis square;set(gca,'FontSize',11);set(gca,'TickDir','out');box off;
offsetAxes
cd(save_folder);saveas(gcf, 'correlation ILNdelta.pdf');
%% Plot r values under both conditions (ILN with L6 or without)
fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 500, 200, 200]);set(gcf,'color','w');
 b1=bar(2,abs(nanmean(col_ra(:))),0.7);b1.FaceColor=[0.5 0.5 0.5];b1.EdgeColor='none';hold on;
  b1=bar(3,abs(nanmean(col_ra1(:))),0.7);b1.FaceColor='c';b1.EdgeColor='none';
  b1=bar(4,abs(nanmean(col_ra2(:))),0.7);b1.FaceColor='m';b1.EdgeColor='none';
%      r=k;rng(k);r1 = r-0.01 + (0.2)*rand(length(nanmean(valdci_all(:,idx_modules{k}),2)),1);
% sc1=scatter(r1,nanmean(valdci_all(:,idx_modules{k}),2),12,'k','filled');hold on;sc1.MarkerEdgeColor=[1 1 1];sc1.MarkerEdgeAlpha=0.5;sc1.MarkerFaceColor=[0.2 0.2 0.2];
     hold on;errorbar([2],[abs(nanmean(col_ra(:)))],[nanstd(col_ra(:))/sqrt(length(col_ra(:)))]...
    , 'LineStyle', 'none', ... 
         'Color', 'k', 'LineWidth', 1);hold on;set(b1,'ShowBaseLine','off');
      hold on;errorbar([3],[abs(nanmean(col_ra1(:)))],[nanstd(col_ra(:))/sqrt(length(col_ra1(:)))]...
    , 'LineStyle', 'none', ... 
         'Color', 'k', 'LineWidth', 1);hold on;set(b1,'ShowBaseLine','off');
         hold on;errorbar([4],[abs(nanmean(col_ra2(:)))],[nanstd(col_ra(:))/sqrt(length(col_ra2(:)))]...
    , 'LineStyle', 'none', ... 
         'Color', 'k', 'LineWidth', 1);hold on;set(b1,'ShowBaseLine','off');
      ylabel('Correlation r');xticks([2:1:4]);xticklabels({'incl L6+5','excl L5','excl L6'});
axis square;set(gca,'FontSize',11);set(gca,'TickDir','out');box off;
xlim([1 4.5]);
%offsetAxes
%% Sorted modules ILNc-ILNi DIfferences of ILN, sensory have higher difference ipsi contra than 'higher brain areas'
idx=[];idx=1;
order_modules=[];
temp_all={i_v1aafm c_v1aafm; i_s1aafm c_s1aafm; i_m1aafm c_m1aafm};
temp_color=[v1_color ;s1_color; m1_color];idx_modules={frontal_idx lateral_idx somamo_idx visual_idx medial_idx aud_idx};module_names={'Pfron','Lat','SoMo','Vis','Med','Aud'};
%plot
fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 500, 700, 225]);set(gcf,'color','w');t=tiledlayout('horizontal');
t.TileSpacing = 'compact';t.Padding = 'compact';  
hold on;
for i=1:3
nexttile
rm1=[];rm1=temp_all{i,1};rm2=[];rm2=temp_all{i,2};
[vali_all valdci_all val_i val_dci xerr yerr r_a p_a r_avg p_avg] = anatomy_correlation(rm1, rm2, idx);
 temp_c=temp_color(i,:);
 hold on
 temp_p=[];p=[];tbl=[];stats=[];err_sem=[];order_err=[];
 for k=1:6
     temp_p(:,k)=nanmean(valdci_all(:,idx_modules{k}),2);
   %   if length(find(~isnan(temp_p(:,k))))>1 
   % 
   %   else
   % temp_p(find(~isnan(temp_p(:,k))),k);
   %       temp_p(find(~isnan(temp_p(:,k))),k)=NaN;
   %   end
     %sort_means(k)=nanmean(nanmean(valdci_all(:,idx_modules{k}),2));
      sort_means(k)=nanmean(temp_p(:,k));
     err_sem(k)=nanstd(temp_p(:,k))/(sqrt(length(~isnan(temp_p(:,k)))))
     

 end
 [order_m kk]=sort(sort_means,'descend');
 order_modules(:,i)=kk;
order_err=err_sem(kk);
 for k=1:6    
 hold on
 b1=bar(k,order_m(k),0.7);b1.FaceColor=temp_c;b1.EdgeColor='none';
%      r=k;rng(k);r1 = r-0.01 + (0.2)*rand(length(nanmean(valdci_all(:,idx_modules{k}),2)),1);
% sc1=scatter(r1,nanmean(valdci_all(:,idx_modules{k}),2),12,'k','filled');hold on;sc1.MarkerEdgeColor=[1 1 1];sc1.MarkerEdgeAlpha=0.5;sc1.MarkerFaceColor=[0.2 0.2 0.2];
     hold on;errorbar([k],[order_m(k)],[order_err(k)]...
    , 'LineStyle', 'none', ... 
         'Color', 'k', 'LineWidth', 1);hold on;set(b1,'ShowBaseLine','off');
 end
 xticks([1:6]);hold on;box off;xticklabels(module_names(kk))
ax=gca
ax.XTickLabel{1} = ['\color{red}' ax.XTickLabel{1}];
ax.XTickLabel{2} = ['\color{red}' ax.XTickLabel{2}];
if i==1 | i==2
ax.XTickLabel{3} = ['\color{red}' ax.XTickLabel{3}];
end

 set(gca,'FontSize',11);set(gca,'TickDir','out');box off;
 xtickangle(45);
 ylim([-0.14 0.4]);
 ylabel('ILN contra - ILN ipsi');
 if i==1
%h = gca;h.XAxis.Visible = 'off';
 elseif i==2
h = gca;h.YAxis.Visible = 'off';
text(0.1,0.35,'sorted by largest delta');
 else i==3
  h = gca;h.YAxis.Visible = 'off';  
  

 end
%  [p,tbl,stats] = anova1([temp_p])
% presults = multcompare(stats)
offsetAxes
end
cd(save_folder);saveas(gcf, 'ILN_delta_modules_barplot.pdf');

%% Average in two major areas
idx=1;
temp_all={i_v1aafm c_v1aafm; i_s1aafm c_s1aafm; i_m1aafm c_m1aafm};
temp_color=[v1_color ;s1_color; m1_color];
senso_mo=[idx_modules{3} idx_modules{4} idx_modules{6}];
higher_a=[idx_modules{1} idx_modules{2} idx_modules{5}];
pan_title={'VISp','SSp-bfd','MOp'};

fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 500, 520, 220]);set(gcf,'color','w');t=tiledlayout('horizontal');
t.TileSpacing = 'compact';t.Padding = 'compact';  
for i=1:3
    nexttile
rm1=[];rm1=temp_all{i,1};rm2=[];rm2=temp_all{i,2};

[vali_all valdci_all val_i val_dci xerr yerr r_a p_a r_avg p_avg valc_all] = anatomy_correlation(rm1, rm2, idx);

[u p1]=ttest(nanmean(valdci_all(:,higher_a),2),nanmean(valdci_all(:,senso_mo),2))

b1=bar(1,nanmean(nanmean(valdci_all(:,senso_mo),2)),0.5);b1.FaceColor='r';hold on;b1.EdgeColor='none';
b2=bar(2,nanmean(nanmean(valdci_all(:,higher_a),2)),0.5);b2.FaceColor=[0.9 0.9 0.9];set(b1,'ShowBaseLine','off');set(b2,'ShowBaseLine','off');b2.EdgeColor='none';
  r=1;rng(1);r1 = r-0.01 + (0.2)*rand(length(nanmean(valdci_all(:,senso_mo),2)),1);
sc1=scatter(r1,nanmean(valdci_all(:,senso_mo),2),15,'k','filled');hold on;sc1.MarkerEdgeColor=[1 1 1];sc1.MarkerEdgeAlpha=0.5;sc1.MarkerFaceColor=[0.2 0.2 0.2];

r=2;rng(1);r1 = r-0.01 + (0.2)*rand(length(nanmean(valdci_all(:,higher_a),2)),1);
sc1=scatter(r1,nanmean(valdci_all(:,higher_a),2),15,'k','filled');hold on;sc1.MarkerEdgeColor=[1 1 1];sc1.MarkerEdgeAlpha=0.5;sc1.MarkerFaceColor=[0.2 0.2 0.2];
hold on;errorbar([1],[nanmean(nanmean(valdci_all(:,senso_mo),2))],[nanstd(nanmean(valdci_all(:,senso_mo),2))/sqrt(length(nanmean(valdci_all(:,senso_mo),2)))]...
     , 'LineStyle', 'none', ... 
          'Color', 'k', 'LineWidth', 1);hold on;set(b1,'ShowBaseLine','off');
    hold on;errorbar([2],[nanmean(nanmean(valdci_all(:,higher_a),2))],[nanstd(nanmean(valdci_all(:,higher_a),2))/sqrt(length(nanmean(valdci_all(:,higher_a),2)))]...
     , 'LineStyle', 'none', ... 
          'Color', 'k', 'LineWidth', 1);hold on;set(b1,'ShowBaseLine','off');

hold on;box off
xticks([1:2])
xlim([0 2.5])
   
     set(gca,'FontSize',11);set(gca,'TickDir','out');box off;
  
xticklabels({'v-sm-a','pf-m-l'})
xtickangle(45)
 ylim([0 0.25]);
 ax=gca
ax.XTickLabel{1} = ['\color{red}' ax.XTickLabel{1}];
ax.XTickLabel{2} = ['\color{gray}' ax.XTickLabel{2}];
 ylabel('ILN contra - ILN ipsi');
 title(pan_title{i},'Color',temp_color(i,:),'FontWeight','normal')
 if i==1
%h = gca;h.XAxis.Visible = 'off';
text(1.3,0.22,'**','FontSize',18);

 elseif i==2
h = gca;h.YAxis.Visible = 'off';
text(1.3,0.22,'***','FontSize',18);

 else i==3
  h = gca;h.YAxis.Visible = 'off';   
text(1.3,0.22,'**','FontSize',18);
 end

offsetAxes
end
cd(save_folder);saveas(gcf, 'ILN_delta_modules_2largegroups.pdf');
%% two groups ipsi vs contra ILN absolute
idx=1;
rm1=[];rm1=cat(3,i_v1aafm,i_s1aafm,i_m1aafm);
rm2=[];rm2=cat(3,c_v1aafm,c_s1aafm,c_m1aafm);
temp_color=[v1_color ;s1_color; m1_color];
senso_mo=[idx_modules{3} idx_modules{4} idx_modules{6}];
higher_a=[idx_modules{1} idx_modules{2} idx_modules{5}];
[vali_all valdci_all val_i val_dci xerr yerr r_a p_a r_avg p_avg valc_all xerr_c] = anatomy_correlation(rm1, rm2, idx);


% 
fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 500, 350, 220]);set(gcf,'color','w');
b1=bar(1,[nanmean(nanmean(vali_all(:,senso_mo),2)) ;nanmean(nanmean(valc_all(:,senso_mo),2))],0.7);
b1(1).FaceColor=[0.8 0.8 0.8];hold on;b1(1).EdgeColor='r';b1(2).FaceColor=[0.3 0.3 0.3];hold on;b1(2).EdgeColor='r';set(b1,'ShowBaseLine','off');
b2=bar(2,[nanmean(nanmean(vali_all(:,higher_a),2)) ;nanmean(nanmean(valc_all(:,higher_a),2))],0.7);set(b1,'ShowBaseLine','off');set(b2,'ShowBaseLine','off');
b2(1).FaceColor=[0.8 0.8 0.8];hold on;b2(1).EdgeColor=[0.5 0.5 0.5];b2(2).FaceColor=[0.3 0.3 0.3];hold on;b2(2).EdgeColor=[0.5 0.5 0.5]

%indivdual animals
r=1;
rng(1);r1 = r-0.2 + (0.1)*rand(length(nanmean(valc_all(:,senso_mo),2)),1);
rng(1);r2 = r +  (0.25)*rand(length(nanmean(valc_all(:,senso_mo),2)),1);
sc1=scatter([r1],[nanmean(vali_all(:,senso_mo),2)],5,'k','filled');hold on;
sc1(1).MarkerEdgeColor=[0.5 0.5 0.5];sc1(1).MarkerFaceColor=[0.2 0.2 0.2];sc1(1).MarkerFaceAlpha=[0.5];sc1(1).MarkerEdgeAlpha=[0.5];
sc1=scatter([r2],[nanmean(valc_all(:,senso_mo),2)],5,'k','filled');hold on;
sc1(1).MarkerEdgeColor=[0.5 0.5 0.5];sc1(1).MarkerFaceColor=[0.2 0.2 0.2];sc1(1).MarkerFaceAlpha=[0.5];sc1(1).MarkerEdgeAlpha=[0.5];
r=2;
rng(1);r1 = r-0.2 + (0.1)*rand(length(nanmean(valc_all(:,senso_mo),2)),1);
rng(1);r2 = r +  (0.25)*rand(length(nanmean(valc_all(:,senso_mo),2)),1);
sc1=scatter([r1],[nanmean(vali_all(:,higher_a),2)],5,'k','filled');hold on;
sc1(1).MarkerEdgeColor=[0.5 0.5 0.5];sc1(1).MarkerFaceColor=[0.2 0.2 0.2];sc1(1).MarkerFaceAlpha=[0.5];sc1(1).MarkerEdgeAlpha=[0.5];
sc1=scatter([r2],[nanmean(valc_all(:,higher_a),2)],5,'k','filled');hold on;
sc1(1).MarkerEdgeColor=[0.5 0.5 0.5];sc1(1).MarkerFaceColor=[0.2 0.2 0.2];sc1(1).MarkerFaceAlpha=[0.5];sc1(1).MarkerEdgeAlpha=[0.5];
hold on;   
xpos = b1(1).XData + b1(1).XOffset;
errorbar(xpos, b1(1).YData, nanstd(nanmean(vali_all(:,senso_mo),2))/sqrt((length(nanmean(vali_all(:,senso_mo),2)))), 'LineStyle', 'none', ... 
        'Color', 'k', 'LineWidth', 1);
 text(0.77,-0.04,'ipsi','FontSize',11)
xpos = b1(2).XData + b1(2).XOffset;
errorbar(xpos, b1(2).YData, nanstd(nanmean(valc_all(:,senso_mo),2))/sqrt((length(nanmean(valc_all(:,senso_mo),2)))), 'LineStyle', 'none', ... 
        'Color', 'k', 'LineWidth', 1);
text(1,-0.04,'contra','FontSize',11)

xpos = b2(1).XData + b2(1).XOffset;
errorbar(1.85, b2(1).YData, nanstd(nanmean(vali_all(:,higher_a),2))/sqrt((length(nanmean(vali_all(:,higher_a),2)))), 'LineStyle', 'none', ... 
        'Color', 'k', 'LineWidth', 1);
 text(1.77,-0.04,'ipsi','FontSize',11)
xpos = b2(2).XData + b2(2).XOffset;
errorbar(2.15, b2(2).YData, nanstd(nanmean(valc_all(:,higher_a),2))/sqrt((length(nanmean(valc_all(:,higher_a),2)))), 'LineStyle', 'none', ... 
        'Color', 'k', 'LineWidth', 1);
    text(2,-0.04,'contra','FontSize',11)
[u p1]=ttest(nanmean(vali_all(:,senso_mo),2),nanmean(valc_all(:,senso_mo),2))
text(0.9,1.05,'***','FontSize',18,'Color','k')
[u p1]=ttest(nanmean(vali_all(:,higher_a),2),nanmean(valc_all(:,higher_a),2))
text(1.95,1.05,'n.s.','FontSize',11,'Color','k')

   text(0.85,1.2,'v-sm-a','FontSize',11,'Color','r')
   text(1.85,1.2,'pf-m-l','FontSize',11,'Color',[0.5 0.5 0.5])
   ylabel('Cortical hierarchy (ILN)')
box off;set(gca,'FontSize',11);set(gca,'TickDir','out')
h = gca;h.XAxis.Visible = 'off';
ylim([0 1.2])
cd(save_folder);saveas(gcf, 'ILN_large_groups_ipsivscontra.pdf');
%% FIGURE 6

%% L23 vs L5 and L6a 
dat_all1={};dat_all1={i_v1aafm i_s1aafm i_m1aafm};
dat_all2={};dat_all2={c_v1aafm c_s1aafm c_m1aafm};
idx_modules={frontal_idx lateral_idx somamo_idx visual_idx medial_idx aud_idx};
temp_color=[v1_color ;s1_color; m1_color];
senso_mo=[idx_modules{3} idx_modules{4} idx_modules{6}];
higher_a=[idx_modules{1} idx_modules{2} idx_modules{5}];
pan_title={'VISp','SSp-bfd','MOp'};
temp_color=[v1_color ;s1_color; m1_color];
temp_c1=[0.8 0.8 0.8];
temp_c2=[0.3 0.3 0.3];


fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 100, 350, 900]);set(gcf,'color','w');t=tiledlayout('vertical');
t.TileSpacing = 'compact';t.Padding = 'compact';
for j=1:3
    
rm1=[];rm2=[];rm1=dat_all1{j};rm2=dat_all2{j};
hdi_l23=(nanmean(rm2(2,:,:),3)-nanmean(rm1(2,:,:),3));
hdi_l56=(nanmean(sum(rm2(4:5,:,:)),3)-nanmean(sum(rm1(4:5,:,:)),3));
hdi_l5=(nanmean(rm2(4,:,:),3)-nanmean(rm1(4,:,:),3));
hdi_l6=(nanmean(rm2(5,:,:),3)-nanmean(rm1(5,:,:),3));

%ymean1=[];ymean1=nanmean(squeeze(rm1(2,idx_modules{i},:)));
nexttile
order_om=[];
order_om=order_modules(:,j);
for i=1:3
yerr2=[];yerr2=[];yerr3=[];yerr3=[];yerr4=[];
ymean2=[];ymean2=nanmean(squeeze(nanmean(rm2(2,idx_modules{order_om(i)},:)-rm1(2,idx_modules{order_om(i)},:))));
ymean3=[];ymean3=nanmean(squeeze(nanmean(rm2(4,idx_modules{order_om(i)},:)-rm1(4,idx_modules{order_om(i)},:))));
ymean4=[];ymean4=nanmean(squeeze(nanmean(rm2(5,idx_modules{order_om(i)},:)-rm1(5,idx_modules{order_om(i)},:))));
ymean_all=[];ymean_all=[ymean2 ymean3 ymean4];
yerr2=nanstd(squeeze(nanmean(rm2(2,idx_modules{order_om(i)},:)-rm1(2,idx_modules{order_om(i)},:))))/sqrt(size(rm1,3));
yerr3=nanstd(squeeze(nanmean(rm2(4,idx_modules{order_om(i)},:)-rm1(4,idx_modules{order_om(i)},:))))/sqrt(size(rm1,3));
yerr4=nanstd(squeeze(nanmean(rm2(5,idx_modules{order_om(i)},:)-rm1(5,idx_modules{order_om(i)},:))))/sqrt(size(rm1,3));
 % if j==3 & i==2 | j==3 & i==1
 %     b1=bar(i,[0 0 0]);hold on;b1(1).FaceColor=[0.9 0.9 0.9];b1(2).FaceColor=[173 216 230]/256;b1(3).FaceColor=[0.55 0.55 0.55];
 %     b1(1).EdgeColor='none';b1(2).EdgeColor='none';b1(3).EdgeColor='none';
 %     set(b1,'ShowBaseLine','off');
 % hold on;h=errorbar([i-0.26],[0],[0], 'LineStyle', 'none','Color', 'k', 'LineWidth', 1,'CapSize',0);
 % hold on;h=errorbar([i+0],[0],[0], 'LineStyle', 'none','Color', 'k', 'LineWidth', 1,'CapSize',0);
 % hold on;h=errorbar([i+0.26],[0],[0], 'LineStyle', 'none','Color', 'k', 'LineWidth', 1,'CapSize',0);
 % else
b1=bar(i,ymean_all);hold on;b1(1).FaceColor=[0.9 0.9 0.9];b1(2).FaceColor=[173 216 230]/256;b1(3).FaceColor=[0.55 0.55 0.55];
set(b1,'ShowBaseLine','off');b1(1).EdgeColor='none';b1(2).EdgeColor='none';b1(3).EdgeColor='none';
 hold on;h=errorbar([i-0.22],[ymean2],[yerr2], 'LineStyle', 'none','Color', 'k', 'LineWidth', 1,'CapSize',0);
 hold on;h=errorbar([i+0],[ymean3],[yerr3], 'LineStyle', 'none','Color', 'k', 'LineWidth', 1,'CapSize',0);
 hold on;h=errorbar([i+0.22],[ymean4],[yerr4], 'LineStyle', 'none','Color', 'k', 'LineWidth', 1,'CapSize',0);
 %end
end
ylim([-0.25 0.2]);yticks([-0.25:0.1:0.25])
box off;
xticks([1:6]);xticklabels(module_names(order_om));set(gca,'FontSize',11);set(gca,'TickDir','out')
ylabel('Fractional change (contra - ipsi)');
text(1.5,0.2,pan_title{j},'Color',temp_color(j,:),'FontWeight','normal')
 if j==1
%h = gca;h.XAxis.Visible = 'off';
legend({'L2/3','L5','L6a'}); legend box off
 elseif j==2
%h = gca;h.YAxis.Visible = 'off';

 else j==3
  %h = gca;h.YAxis.Visible = 'off';   

 end

offsetAxes
end

cd(save_folder);saveas(gcf, 'layer_change_ILN.pdf');
%% Data for flat maps for MT
dat_all1={};dat_all1={i_v1aafm i_s1aafm i_m1aafm};
dat_all2={};dat_all2={c_v1aafm c_s1aafm c_m1aafm};
senso_mo=[idx_modules{3} idx_modules{4} idx_modules{6}];
higher_a=[idx_modules{1} idx_modules{2} idx_modules{5}];
ymean1=[];
ymean2=[];
ymean3=[];
all_l23p=[];all_l5p=[];all_l6p=[];
for j=1:3
rm1=[];rm2=[];rm1=dat_all1{j};rm2=dat_all2{j};
% ymean1(:,j)=nanmean(rm2(2,:,:)-rm1(2,:,:),3);
% 
% ymean2(:,j)=nanmean(rm2(4,:,:)-rm1(4,:,:),3);
% ymean3(:,j)=nanmean(rm2(5,:,:)-rm1(5,:,:),3);
l23c=squeeze(rm2(2,:,:));
l5c=squeeze(rm2(4,:,:));
l6c=squeeze(rm2(5,:,:));
l23i=squeeze(rm1(2,:,:));
l5i=squeeze(rm1(4,:,:));
l6i=squeeze(rm1(5,:,:));
ymean1(:,j)=nanmean(l23c-l23i,2);
ymean2(:,j)=nanmean(l5c-l5i,2);
ymean3(:,j)=nanmean(l6c-l6i,2);
p_sl23=[];p_sl5=[];p_sl6=[];
    for m=1:45
    [u p_sl23(m)]=ttest(l23c(m,:),l23i(m,:));
    [u p_sl5(m)]=ttest(l5c(m,:),l5i(m,:));
    [u p_sl6(m)]=ttest(l6c(m,:),l6i(m,:));
    end
 all_l23p(:,j)=p_sl23;
 all_l5p(:,j)=p_sl5;
  all_l6p(:,j)=p_sl6;
end
%% save for MT flatmaps/table
pthresh1=[];pthresh2=[];pthresh3=[];
pthresh1=all_l23p<0.05;
pthresh2=all_l5p<0.05;
pthresh3=all_l6p<0.05;
%% 

delta_layers_cortex =table(cortex_names(:,2),ymean1(:,1),ymean2(:,1),ymean3(:,1),ymean1(:,2),ymean2(:,2),ymean3(:,2),ymean1(:,3),ymean2(:,3),ymean3(:,3),...
    pthresh1(:,1),pthresh2(:,1),pthresh3(:,1),pthresh1(:,2),pthresh2(:,2),pthresh3(:,2),pthresh1(:,3),pthresh2(:,3),pthresh3(:,3),...
    'variablenames',{'Area','V1dL2/3','V1dL5','V1dL6','S1dL2/3','S1dL5','S1dL6','M1dL2/3','M1dL5','M1dL6','V1L2/3sig','V1L5sig','V1L6sig'...
    'S1L2/3sig','S1L5sig','S1L6sig','M1L2/3sig','M1L5sig','M1L6sig'});
writetable(delta_layers_cortex,'delta_layers_cortex.csv')
%% show fraction of layers L2/3, L5 and L6 only for the senso- motor module
dat_all1={};dat_all1={i_v1aafm i_s1aafm i_m1aafm};
dat_all2={};dat_all2={c_v1aafm c_s1aafm c_m1aafm};
idx_modules={frontal_idx lateral_idx somamo_idx visual_idx medial_idx aud_idx};
temp_color=[v1_color ;s1_color; m1_color];
senso_mo=[idx_modules{3} idx_modules{4} idx_modules{6}];
higher_a=[idx_modules{1} idx_modules{2} idx_modules{5}];
pan_title={'VISp','SSp-bfd','MOp'};
temp_color=[v1_color ;s1_color; m1_color];
temp_c1=[0.8 0.8 0.8];
temp_c2=[0.3 0.3 0.3];
allsm_l23i=[];allsm_l5i=[];allsm_l6i=[];
allsm_l23c=[];allsm_l5c=[];allsm_l6c=[];
fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 100, 250, 600]);set(gcf,'color','w');t=tiledlayout('vertical');
t.TileSpacing = 'compact';t.Padding = 'compact';
for j=1:3
    
rm1=[];rm2=[];rm1=dat_all1{j};rm2=dat_all2{j};

l23smi=squeeze(nanmean(rm1(2,senso_mo,:),2));
l5smi=squeeze(nanmean(rm1(4,senso_mo,:),2));
l6smi=squeeze(nanmean(rm1(5,senso_mo,:),2));
l23smc=squeeze(nanmean(rm2(2,senso_mo,:),2));
l5smc=squeeze(nanmean(rm2(4,senso_mo,:),2));
l6smc=squeeze(nanmean(rm2(5,senso_mo,:),2));
allsm_l23i=[allsm_l23i; l23smi];
allsm_l5i=[allsm_l5i ;l5smi];
allsm_l6i=[allsm_l6i ;l6smi];
allsm_l23c=[allsm_l23c ;l23smc];
allsm_l5c=[allsm_l5c; l5smc];
allsm_l6c=[allsm_l6c; l6smc];
%ymean1=[];ymean1=nanmean(squeeze(rm1(2,idx_modules{i},:)));
[uu kk]=ttest(l23smi,l23smc)
[uu5 kk5]=ttest(l5smi,l5smc)
[uu6 kk6]=ttest(l6smi,l6smc)
nexttile

b1=bar(1,[mean(l23smi); mean(l23smc)]);hold on;b1(1).FaceColor=temp_c1;b1(2).FaceColor=temp_c2;
set(b1,'ShowBaseLine','off');b1(1).EdgeColor='none';b1(2).EdgeColor='none';
b1=bar(2,[mean(l5smi); mean(l5smc)]);hold on;b1(1).FaceColor=temp_c1;b1(2).FaceColor=temp_c2;
set(b1,'ShowBaseLine','off');b1(1).EdgeColor='none';b1(2).EdgeColor='none';
b1=bar(3,[mean(l6smi); mean(l6smc)]);hold on;b1(1).FaceColor=temp_c1;b1(2).FaceColor=temp_c2;
set(b1,'ShowBaseLine','off');b1(1).EdgeColor='none';b1(2).EdgeColor='none';
r=1;
rng(1);r1 = r-0.2 + (0.1)*rand(length(l23smi),1);
rng(1);r2 = r +  (0.35)*rand(length(l23smc),1);
sc1=scatter([r1],[l23smi],5,'k','filled');hold on;
sc1(1).MarkerEdgeColor=[0.5 0.5 0.5];sc1(1).MarkerFaceColor=[0.2 0.2 0.2];sc1(1).MarkerFaceAlpha=[0.5];sc1(1).MarkerEdgeAlpha=[0.5];
sc1=scatter([r2],[l23smc],5,'k','filled');hold on;
sc1(1).MarkerEdgeColor=[0.5 0.5 0.5];sc1(1).MarkerFaceColor=[0.2 0.2 0.2];sc1(1).MarkerFaceAlpha=[0.5];sc1(1).MarkerEdgeAlpha=[0.5];
r=2;
rng(1);r1 = r-0.2 + (0.1)*rand(length(l5smi),1);
rng(1);r2 = r +  (0.35)*rand(length(l5smc),1);
sc1=scatter([r1],[l5smi],5,'k','filled');hold on;
sc1(1).MarkerEdgeColor=[0.5 0.5 0.5];sc1(1).MarkerFaceColor=[0.2 0.2 0.2];sc1(1).MarkerFaceAlpha=[0.5];sc1(1).MarkerEdgeAlpha=[0.5];
sc1=scatter([r2],[l5smc],5,'k','filled');hold on;
r=3;
rng(1);r1 = r-0.2 + (0.1)*rand(length(l6smi),1);
rng(1);r2 = r +  (0.35)*rand(length(l6smc),1);
sc1=scatter([r1],[l6smi],5,'k','filled');hold on;
sc1(1).MarkerEdgeColor=[0.5 0.5 0.5];sc1(1).MarkerFaceColor=[0.2 0.2 0.2];sc1(1).MarkerFaceAlpha=[0.5];sc1(1).MarkerEdgeAlpha=[0.5];
sc1=scatter([r2],[l6smc],5,'k','filled');hold on;
sc1(1).MarkerEdgeColor=[0.5 0.5 0.5];sc1(1).MarkerFaceColor=[0.2 0.2 0.2];sc1(1).MarkerFaceAlpha=[0.5];sc1(1).MarkerEdgeAlpha=[0.5];
  hold on;h=errorbar([1-0.15],[mean(l23smi)],[nanstd(l23smi)/sqrt(length(l23smi))], 'LineStyle', 'none','Color', 'k', 'LineWidth', 1,'CapSize',0);
  hold on;h=errorbar([1+0.15],[mean(l23smc)],[nanstd(l23smc)/sqrt(length(l23smc))], 'LineStyle', 'none','Color', 'k', 'LineWidth', 1,'CapSize',0);
    hold on;h=errorbar([2-0.15],[mean(l5smi)],[nanstd(l5smi)/sqrt(length(l5smi))], 'LineStyle', 'none','Color', 'k', 'LineWidth', 1,'CapSize',0);
  hold on;h=errorbar([2+0.15],[mean(l5smc)],[nanstd(l5smc)/sqrt(length(l5smc))], 'LineStyle', 'none','Color', 'k', 'LineWidth', 1,'CapSize',0);
    hold on;h=errorbar([3-0.15],[mean(l6smi)],[nanstd(l6smi)/sqrt(length(l6smi))], 'LineStyle', 'none','Color', 'k', 'LineWidth', 1,'CapSize',0);
  hold on;h=errorbar([3+0.15],[mean(l6smc)],[nanstd(l6smc)/sqrt(length(l6smc))], 'LineStyle', 'none','Color', 'k', 'LineWidth', 1,'CapSize',0);
  box off;

 ylim([0 0.5]);
 ylabel('Fraction of neurons');
 xlim([0.5 3.5]);
 xticks([1:3]);xticklabels({'L2/3','L5','L6a'});set(gca,'FontSize',11);set(gca,'TickDir','out')
offsetAxes
 text(1.6,0.5,pan_title{j},'Color',temp_color(j,:),'FontWeight','normal');
 
 text(0.8,0.48,'***','Color','k','FontSize',18);
 if j==1 | j==2
 text(2.8,0.48,'***','Color','k','FontSize',18);
 end
%h = gca;h.XAxis.Visible = 'off';

% 
% %h = gca;h.YAxis.Visible = 'off';
% 
% 
%   %h = gca;h.YAxis.Visible = 'off';   
% 
% 
% 
 
end
%% 
cd(save_folder);saveas(gcf, 'fraction_layer_change_semo.pdf');

%% all together 
fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 100, 400, 250]);set(gcf,'color','w');
b1=bar(1,[mean(allsm_l23i); mean(allsm_l23c)]);hold on;b1(1).FaceColor=temp_c1;b1(2).FaceColor=temp_c2;
set(b1,'ShowBaseLine','off');b1(1).EdgeColor='none';b1(2).EdgeColor='none';
b1=bar(2,[mean(allsm_l5i); mean(allsm_l5c)]);hold on;b1(1).FaceColor=temp_c1;b1(2).FaceColor=temp_c2;
set(b1,'ShowBaseLine','off');b1(1).EdgeColor='none';b1(2).EdgeColor='none';
b1=bar(3,[mean(allsm_l6i); mean(allsm_l6c)]);hold on;b1(1).FaceColor=temp_c1;b1(2).FaceColor=temp_c2;
set(b1,'ShowBaseLine','off');b1(1).EdgeColor='none';b1(2).EdgeColor='none';
box off;
[uu kk]=ttest(l23smi,l23smc)
[uu5 kk5]=ttest(l5smi,l5smc)
[uu6 kk6]=ttest(l6smi,l6smc)
%% Decrease increase percentage of sensory motor areas
pan_title={'VISp','SSp-bfd','MOp'};
temp_color=[v1_color ;s1_color; m1_color];
fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 100, 230, 400]);set(gcf,'color','w');t=tiledlayout('vertical');
t.TileSpacing = 'compact';t.Padding = 'compact';
for j=1:3
nexttile
temp_don1=sum(ymean1(senso_mo(pthresh1(senso_mo,j)),j)>0)/24;
temp_don2=sum(ymean2(senso_mo(pthresh2(senso_mo,j)),j)>0)/24
temp_don3=sum(ymean3(senso_mo(pthresh3(senso_mo,j)),j)>0)/24;

temp_don4=sum(ymean1(senso_mo(pthresh1(senso_mo,j)),j)<0)/24;
temp_don5=sum(ymean2(senso_mo(pthresh2(senso_mo,j)),j)<0)/24
temp_don6=sum(ymean3(senso_mo(pthresh3(senso_mo,j)),j)<0)/24;

%donut([temp_don4,temp_don5,temp_don6;temp_don1,temp_don2,temp_don3],{''},{[ .945 .345 .329],[ .376 .741 .408],[ .365 .647 .855 ]})
%xlim([-4 4]);ylim([-4 4]);

b1=barh([temp_don3,temp_don2,temp_don1],0.5);b1.FaceColor=[246 136 68]/256;b1.EdgeColor='none';
hold on
b1=barh([-temp_don6,-temp_don5,-temp_don4],0.5);b1.FaceColor=[48 65 154]/256;b1.EdgeColor='none';
box off;
xlim([-0.6 0.6])
offsetAxes
  h = gca;h.YAxis.Visible = 'off'; 
  if j==1
      text(-0.4,4,'decrease','FontSize',11);
      text(0.05,4,'increase','FontSize',11);
       h = gca;h.XAxis.Visible = 'off'; 
      % text(-0.5,6,'contra - ipsi fraction','FontSize',11)
     %ylim([0 6])
  end
  if j==2
       h = gca;h.XAxis.Visible = 'off'; 
  end
  if j==3
     xlabel('Percentage areas (%)')
  end
   text(-0.78,3,'L2/3','FontSize',11);
      text(-0.75,2,'L5','FontSize',11);
      text(-0.75,1,'L6a','FontSize',11);
  set(gca,'FontSize',11);set(gca,'TickDir','out');
  xticks([-0.6:0.2:0.6]);xticklabels({'60','40','20','0','20','40','60'});
  title(pan_title{j},'Color',temp_color(j,:),'FontWeight','normal')
end
cd(save_folder);saveas(gcf, 'percentage_layer_change_ILN.pdf');


%% 
   %% Remove L6 and then L5 for ILN ipsi and contra
%% 
senso_mo=[idx_modules{3} idx_modules{4} idx_modules{6}];
higher_a=[idx_modules{1} idx_modules{2} idx_modules{5}];

rm2=[];rm2=cat(3,c_v1aafm(:,senso_mo,:), c_s1aafm(:,senso_mo,:),c_m1aafm(:,senso_mo,:));
rm1=[];rm1=cat(3,i_v1aafm(:,senso_mo,:), i_s1aafm(:,senso_mo,:),i_m1aafm(:,senso_mo,:));

idx=1;
[ilni ilnc ilni_e ilnc_e order_r order_rc] = iln_change(rm1,rm2,idx);
%exL6
idx=7;
[ilni_l6 ilnc_l6 ilni_el6 ilnc_el6 order_l6 order_l6c] = iln_change(rm1,rm2,idx);
%ex L5
idx=8;
[ilni_l5 ilnc_l5 ilni_el5 ilnc_el5 order_l5 order_l5c] = iln_change(rm1,rm2,idx);

%% IPsi ILN excl L5 and L6
temp_c1=[0.8 0.8 0.8];
temp_c2=[0.3 0.3 0.3];


fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 300, 400, 300]);set(gcf,'color','w');
alpha1=1;
title('sorted by ipsi ILN','FontWeight','normal','Color','k');
ylabel('ILN');
li=line([1 size(rm1,2)],[0.5 0.5],'Color',[0.5 0.5 0.5],'LineStyle',':');
 hold on;plot(ilni(order_r),'Color',temp_c1);
 hold on;plot(ilni_l5(order_l5),'Color','b');
 hold on;plot(ilni_l6(order_l6),'Color','r');
hold on;h=errorbar([1:size(rm1,2)],[ilni(order_r)],[ilni_e(order_r)]...
    , 'LineStyle', 'none', ... 
        'Color', temp_c1, 'LineWidth', 0.5,'CapSize',2);set([h.Bar, h.Line], 'ColorType', 'truecoloralpha', 'ColorData', [h.Line.ColorData(1:3); 255*alpha1]);
hold on;h=errorbar([1:size(rm1,2)],[ilni_l5(order_l5)],[ilni_el5(order_l5)]...
    , 'LineStyle', 'none', ... 
        'Color', 'b', 'LineWidth', 0.5,'CapSize',2);
hold on;h=errorbar([1:size(rm1,2)],[ilni_l6(order_l6)],[ilni_el6(order_l6)]...
    , 'LineStyle', 'none', ... 
        'Color', 'r', 'LineWidth', 0.5,'CapSize',2);
hold on;pp1=scatter(1:size(rm1,2),ilni(order_r),30,'o');pp1.MarkerEdgeColor=[1 1 1];pp1.MarkerFaceAlpha=1;pp1.MarkerFaceColor=temp_c1;
hold on;pp1=scatter(1:size(rm1,2),ilni_l5(order_l5),30,'o');pp1.MarkerEdgeColor=[1 1 1];pp1.MarkerFaceAlpha=1;pp1.MarkerFaceColor='b';box off;
hold on;pp1=scatter(1:size(rm1,2),ilni_l6(order_l6),30,'o');pp1.MarkerEdgeColor=[1 1 1];pp1.MarkerFaceAlpha=1;pp1.MarkerFaceColor='r';
xticks([1:size(rm1,2)]);ylim([0 1])
ylabel('Ipsi ILN');xlabel('Cortical areas')
set(gca,'FontSize',11);set(gca,'TickDir','out');
text(size(rm1,2)-2,0.35,'reg','Color',[0.7 0.7 0.7]);
text(size(rm1,2)-2,0.28,'excl L5','Color','b');
text(size(rm1,2)-2,0.21,'excl L6','Color','r');
offsetAxes;

     h = gca;h.XAxis.Visible = 'off'  
     text(size(rm1,2)/2-2,-0.1,'Cortical areas','FontSize',11);
     %cd(save_folder);saveas(gcf, 'ILNipsi_sorted_areasALL.pdf');
% [p,tbl,stats]=anova1([ilni' ilni_l5' ilni_l6'])
% presults = multcompare(stats)
%% 
%% Contra ILN excl L5 and L6
temp_c1=[0.8 0.8 0.8];
temp_c2=[0.3 0.3 0.3];
fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 300, 400, 300]);set(gcf,'color','w');
title('sorted by contra ILN','FontWeight','normal','Color','k');
ylabel('ILN')
li=line([1 size(rm2,2)],[0.5 0.5],'Color',[0.5 0.5 0.5],'LineStyle',':');
 hold on;plot(ilnc(order_rc),'Color',temp_c2);
 hold on;plot(ilnc_l5(order_l5c),'Color','b');
 hold on;plot(ilnc_l6(order_l6c),'Color','r');
hold on;h=errorbar([1:size(rm2,2)],[ilnc(order_rc)],[ilnc_e(order_rc)]...
    , 'LineStyle', 'none', ... 
        'Color', temp_c2, 'LineWidth', 0.5,'CapSize',2);set([h.Bar, h.Line], 'ColorType', 'truecoloralpha', 'ColorData', [h.Line.ColorData(1:3); 255*alpha1]);
hold on;h=errorbar([1:size(rm2,2)],[ilnc_l5(order_l5c)],[ilnc_el5(order_l5c)]...
    , 'LineStyle', 'none', ... 
        'Color', 'b', 'LineWidth', 0.5,'CapSize',2);
hold on;h=errorbar([1:size(rm2,2)],[ilnc_l6(order_l6c)],[ilnc_el6(order_l6c)]...
    , 'LineStyle', 'none', ... 
        'Color', 'r', 'LineWidth', 0.5,'CapSize',2);
hold on;pp1=scatter(1:size(rm2,2),ilnc(order_rc),30,'o');pp1.MarkerEdgeColor=[1 1 1];pp1.MarkerFaceAlpha=1;pp1.MarkerFaceColor=temp_c2;
hold on;pp1=scatter(1:size(rm2,2),ilnc_l5(order_l5c),30,'o');pp1.MarkerEdgeColor=[1 1 1];pp1.MarkerFaceAlpha=1;pp1.MarkerFaceColor='b';box off;
hold on;pp1=scatter(1:size(rm2,2),ilnc_l6(order_l6c),30,'o');pp1.MarkerEdgeColor=[1 1 1];pp1.MarkerFaceAlpha=1;pp1.MarkerFaceColor='r';
xticks([1:size(rm2,2)]);ylim([0 1])
ylabel('Contra ILN');xlabel('Cortical areas')
set(gca,'FontSize',11);set(gca,'TickDir','out');
text(size(rm2,2)-2,0.35,'reg','Color',[0.3 0.3 0.3]);
text(size(rm2,2)-2,0.28,'excl L5','Color','b');
text(size(rm2,2)-2,0.21,'excl L6','Color','r');
offsetAxes;

     h = gca;h.XAxis.Visible = 'off'  
     text(size(rm2,2)/2-2,-0.1,'Cortical areas','FontSize',11);
     %cd(save_folder);saveas(gcf, 'ILNipsi_sorted_areasALL.pdf');
[p,tbl,stats]=anova1([ilni' ilni_l5' ilni_l6'])
presults = multcompare(stats)
% [p,tbl,stats]=anova1([ilnc' ilnc_l5' ilnc_l6'])
% presults = multcompare(stats)
%% Violoinplot
%% Violin plots for IPSI
g1=[];g2=[];g3=[];
p1=[];p2=[];p3=[];
p1=ilni';p2=ilni_l5';p3=ilni_l6';
par=[];par=[p1; p2; p3];
g1=[];g1=ones(1,length(p1));
g2=[];g2=ones(1,length(p2))*2;
g3=[];g3=ones(1,length(p3))*3;
gro=[];gro=[g1 g2 g3]';
fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 500, 200, 300]);set(gcf,'color','w');
%hold on;line([0.8 6.5],[1 1],'LineStyle',':','Color','k','LineWidth',1)
violins = violinplot(par, gro,'ViolinColor',[[0.8 0.8 0.8];[0 0 1];[1 0 0]],'ShowMean', false,'ShowMedian', true);box off;
xlim([0.5 3.5]);ylabel('Ipsi ILN');set(gca,'FontSize',12);
% h = gca;h.XAxis.Visible = 'off';
% set(gca,'xtick',[]);
xticklabels({'Reg','L5excl','L6excl'});
%ylim([-1.2 1.2]);
ylim([-0.05 1]);
set(gca,'TickDir','out'); 
hold on;line([1 2],[0.87 0.87],'LineStyle','-','Color','k','LineWidth',0.5);
hold on;line([1 3],[0.94 0.94],'LineStyle','-','Color','k','LineWidth',0.5);
hold on;line([2 3],[0.83 0.83],'LineStyle','-','Color','k','LineWidth',0.5);
hold on;text(1.8,0.98 ,'***','Color','k','FontSize',18);
hold on;text(2.4,0.85 ,'*','Color','k','FontSize',18);
hold on;text(1.15,0.91 ,'n.s.','Color','k','FontSize',11);
offsetAxes;
%% 
%% Violin plots for CONTRA
g1=[];g2=[];g3=[];
p1=[];p2=[];p3=[];
p1=ilnc';p2=ilnc_l5';p3=ilnc_l6';
par=[];par=[p1; p2; p3];
g1=[];g1=ones(1,length(p1));
g2=[];g2=ones(1,length(p2))*2;
g3=[];g3=ones(1,length(p3))*3;
gro=[];gro=[g1 g2 g3]';
fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 500, 200, 300]);set(gcf,'color','w');
%hold on;line([0.8 6.5],[1 1],'LineStyle',':','Color','k','LineWidth',1)
violins = violinplot(par, gro,'ViolinColor',[[0.8 0.8 0.8];[0 0 1];[1 0 0]],'ShowMean', false,'ShowMedian', true);box off;
xlim([0.5 3.5]);ylabel('Contra ILN');set(gca,'FontSize',12);
% h = gca;h.XAxis.Visible = 'off';
% set(gca,'xtick',[]);
xticklabels({'Reg','L5excl','L6excl'});
%ylim([-1.2 1.2]);
ylim([-0.05 1]);
set(gca,'TickDir','out'); 
%hold on;line([5.5 5.5],[-1.2 1.2],'LineStyle',':','Color','k','LineWidth',1)
offsetAxes;
hold on;line([1 2],[0.92 0.92],'LineStyle','-','Color','k','LineWidth',0.5);
hold on;line([1 3],[0.97 0.97],'LineStyle','-','Color','k','LineWidth',0.5);
hold on;line([2 3],[0.89 0.89],'LineStyle','-','Color','k','LineWidth',0.5);
hold on;text(1.8,0.99 ,'***','Color','k','FontSize',18);
hold on;text(2.4,0.91 ,'*','Color','k','FontSize',18);
hold on;text(1.15,0.95 ,'n.s.','Color','k','FontSize',11);
%% 

%% Violin plots for Delta
g1=[];g2=[];g3=[];
p1=[];p2=[];p3=[];
p1=(ilnc-ilni)';p2=(ilnc_l5-ilni)';p3=(ilnc_l6-ilni)';
par=[];par=[p1; p2; p3];
g1=[];g1=ones(1,length(p1));
g2=[];g2=ones(1,length(p2))*2;
g3=[];g3=ones(1,length(p3))*3;
gro=[];gro=[g1 g2 g3]';
fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 500, 200, 300]);set(gcf,'color','w');
%hold on;line([0.8 6.5],[1 1],'LineStyle',':','Color','k','LineWidth',1)
violins = violinplot(par, gro,'ViolinColor',[[0.8 0.8 0.8];[0 0 1];[1 0 0]],'ShowMean', false,'ShowMedian', true);box off;
xlim([0.5 3.5]);ylabel('Delta Contra-Ipsi ILN');set(gca,'FontSize',12);
% h = gca;h.XAxis.Visible = 'off';
% set(gca,'xtick',[]);
xticklabels({'Reg','L5excontra','L6excontra'});
%ylim([-1.2 1.2]);
ylim([-0.5 0.5]);
set(gca,'TickDir','out'); 
%hold on;line([5.5 5.5],[-1.2 1.2],'LineStyle',':','Color','k','LineWidth',1)
offsetAxes;
hold on;line([1 2],[0.4 0.4],'LineStyle','-','Color','k','LineWidth',0.5);
hold on;line([1 3],[0.5 0.5],'LineStyle','-','Color','k','LineWidth',0.5);
hold on;line([2 3],[0.45 0.45],'LineStyle','-','Color','k','LineWidth',0.5);
hold on;text(1.8,0.52 ,'***','Color','k','FontSize',18);
hold on;text(2.2,0.46 ,'***','Color','k','FontSize',18);
hold on;text(1.15,0.43 ,'n.s.','Color','k','FontSize',11);
%% 
[p,tbl,stats]=anova1([p1 p2 p3])
presults = multcompare(stats)

%% Decrease increase percentage of sensory motor areas
pan_title={'VISp','SSp-bfd','MOp'};
temp_color=[v1_color ;s1_color; m1_color];
fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 100, 230, 400]);set(gcf,'color','w');t=tiledlayout('vertical');
t.TileSpacing = 'compact';t.Padding = 'compact';
for j=1:3
nexttile
temp_don1=sum(ymean1(senso_mo(pthresh1(senso_mo,j)),j)>0)/24;
temp_don2=sum(ymean2(senso_mo(pthresh2(senso_mo,j)),j)>0)/24
temp_don3=sum(ymean3(senso_mo(pthresh3(senso_mo,j)),j)>0)/24;

temp_don4=sum(ymean1(senso_mo(pthresh1(senso_mo,j)),j)<0)/24;
temp_don5=sum(ymean2(senso_mo(pthresh2(senso_mo,j)),j)<0)/24
temp_don6=sum(ymean3(senso_mo(pthresh3(senso_mo,j)),j)<0)/24;

%donut([temp_don4,temp_don5,temp_don6;temp_don1,temp_don2,temp_don3],{''},{[ .945 .345 .329],[ .376 .741 .408],[ .365 .647 .855 ]})
%xlim([-4 4]);ylim([-4 4]);

b1=barh([temp_don3-temp_don6,temp_don2-temp_don5,temp_don1-temp_don4],0.5);b1.FaceColor='k';b1.EdgeColor='none';
hold on
box off;
xlim([-0.5 0.5])
offsetAxes
  h = gca;h.YAxis.Visible = 'off'; 
  if j==1
      text(-0.4,4,'decrease','FontSize',11);
      text(0.05,4,'increase','FontSize',11);
       h = gca;h.XAxis.Visible = 'off'; 
      % text(-0.5,6,'contra - ipsi fraction','FontSize',11)
     %ylim([0 6])
  end
  if j==2
       h = gca;h.XAxis.Visible = 'off'; 
  end
  if j==3
     xlabel('Percentage areas (%)')
  end
   text(-0.75,3,'L2/3','FontSize',11);
      text(-0.72,2,'L5','FontSize',11);
      text(-0.72,1,'L6a','FontSize',11);
  set(gca,'FontSize',11);set(gca,'TickDir','out');
  xticks([-0.6:0.2:0.6]);xticklabels({'60','40','20','0','20','40','60'});
  title(pan_title{j},'Color',temp_color(j,:),'FontWeight','normal')
end

%% END


%% Supplementary
%% Supp1

%% Plot all cells ipis+contra  TOTAL numbers for supplement
temp_color=[v1_color ;s1_color; m1_color];
p1=[];p1=squeeze(nansum(nansum(i_v1aam)))+squeeze(nansum(nansum(c_v1aam)));
p2=[];p2=squeeze(nansum(nansum(i_s1aam)))+squeeze(nansum(nansum(c_s1aam)));
p3=[];p3=squeeze(nansum(nansum(i_m1aam)))+squeeze(nansum(nansum(c_m1aam)));
temp_p= [p1 p2 p3];
fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 500, 150, 270]);set(gcf,'color','w');tiledlayout("horizontal");
%title(['Total cell nr: ' num2str(nansum(nansum(temp_p)))],'FontWeight','normal')
%title(['T: 3.17*10^6' ],'FontWeight','normal','FontSize',5)
for i=1:3
    hold on;
    b1=bar(i,nanmean(temp_p(:,i)),0.7);hold on;b1.FaceColor=temp_color(i,:);  set(b1,'ShowBaseLine','off');b1.EdgeColor='none';
      r=i;rng(i);r1 = r-0.01 + (0.2)*rand(length(temp_p(:,i)),1);
sc1=scatter(r1,temp_p(:,i),12,'k','filled');hold on;sc1.MarkerEdgeColor=[1 1 1];sc1.MarkerEdgeAlpha=0.5;sc1.MarkerFaceColor=[0.2 0.2 0.2];
    hold on;errorbar([i],[nanmean(temp_p(:,i))],[nanstd(temp_p(:,i))/sqrt(length(temp_p(:,i)))]...
    , 'LineStyle', 'none', ... 
        'Color', 'k', 'LineWidth', 1);hold on;
end
xlim([0 3.5]);xticks([1:3]);hold on;box off;;ylabel('Total nr of cells (ipsi + contra)');
set(gca,'FontSize',11);set(gca,'TickDir','out');box off;xtickangle(45);ylim([0 410000]);%xticklabels({'VISp','SSp-bfd','MOp'})
h = gca;h.XAxis.Visible = 'off'; 
t1=text(0.4,-60000,'VISp','FontSize',11);set(t1,'Rotation',45);
   t1=text(0.8,-82000,'SSp-bfd','FontSize',11);set(t1,'Rotation',45);
    t1=text(2.4,-60000,'MOp','FontSize',11);set(t1,'Rotation',45);
offsetAxes
cd(save_folder);saveas(gcf, 'Total_nr_cells.pdf');
[p,tbl,stats] = anova1([temp_p])
presults = multcompare(stats)




%% Homotopic areas fractions on contralateral side: SUpplement
%visp_idx=[31];ssp_idx=[18];mop_idx=[25]
fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 500, 170, 250]);set(gcf,'color','w');tiledlayout("horizontal")
b1=bar(1,v1a_r(visp_idx,2),0.7);hold on;b1.FaceColor=v1_color;set(b1,'ShowBaseLine','off');b1.EdgeColor='none';
b2=bar(2,s1a_r(ssp_idx,2),0.7);hold on;b2.FaceColor=s1_color;set(b2,'ShowBaseLine','off');b2.EdgeColor='none';
b3=bar(3,m1a_r(mop_idx,2),0.7);b3.FaceColor=m1_color;set(b3,'ShowBaseLine','off');b3.EdgeColor='none';
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
ylabel('Contra homotopic area fraction')
set(gca,'FontSize',11);set(gca,'TickDir','out');box off;xtickangle(45);
 hold on;line([2 3],[max(get(gca,'YLim')) max(get(gca,'YLim'))],'Color','k');
 hold on;text(2,max(get(gca,'YLim'))+0.03,['***'],'FontSize',18);
  hold on;line([1 3],[max(get(gca,'YLim')+0.1) max(get(gca,'YLim'))+0.1],'Color','k');
  hold on;text(1.5,max(get(gca,'YLim'))-0.06,['***'],'FontSize',18);
  h = gca;h.XAxis.Visible = 'off'; 
t1=text(0.4,-0.15,'VISp','FontSize',11);set(t1,'Rotation',45);
   t1=text(0.72,-0.22,'SSp-bfd','FontSize',11);set(t1,'Rotation',45);
    t1=text(2.4,-0.15,'MOp','FontSize',11);set(t1,'Rotation',45);offsetAxes

cd(save_folder);saveas(gcf, 'homotopicareas_fraction.pdf');

[p,tbl,stats] = anova1([v1a_rca(visp_idx,:)' s1a_rca(ssp_idx,:)' m1a_rca(mop_idx,:)'])
 presults = multcompare(stats)


%% 
%% Correlation between hemisphere areas SLOPe of fraction per area avergae across injection areas
%Calculate R2 and slope of correlation 
%V1
fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 500, 700, 250]);set(gcf,'color','w');tiledlayout("horizontal")
nexttile
for i=6:6%size(v1a_riam,2)
idx = find(~isnan(v1a_riam(:,i)));
mdl = fitlm(v1a_riam(idx,i),v1a_rcam(idx,i));
coefs = polyfit(v1a_riam(idx,i)', v1a_rcam(idx,i)', 1);
slope_v1(i)=coefs(1);
r2_sqv1(i)=mdl.Rsquared.Ordinary;
hold on;ph=plot(mdl,'Color',v1_color);hold on;
ph(2).Color=[0.8 0.8 0.8];
ph(3).Color=[0.8 0.8 0.8];
ph(4).Color=[0.8 0.8 0.8];
end
hold on;title([]);legend('');xlabel('Ipsi Fraction');ylabel('Contra fraction');
set(gca,'FontSize',11);set(gca,'TickDir','out');legend box off;axis square;offsetAxes

%S1
nexttile
for i=6:6%size(s1a_riam,2)
idx = find(~isnan(s1a_riam(:,i)));
mdl = fitlm(s1a_riam(idx,i),s1a_rcam(idx,i));
coefs = polyfit(s1a_riam(idx,i)', s1a_rcam(idx,i)', 1);
slope_s1(i)=coefs(1);
r2_sqs1(i)=mdl.Rsquared.Ordinary;
hold on;ph=plot(mdl,'Color',s1_color);hold on;
ph(2).Color=[0.8 0.8 0.8];
ph(3).Color=[0.8 0.8 0.8];
ph(4).Color=[0.8 0.8 0.8];
end
hold on;title([]);legend('');xlabel('Ipsi Fraction');ylabel('');
set(gca,'FontSize',11);set(gca,'TickDir','out');legend box off;axis square;offsetAxes
%M1
nexttile
for i=6:6%size(m1a_riam,2)
idx = find(~isnan(m1a_riam(:,i)));
mdl = fitlm(m1a_riam(idx,i),m1a_rcam(idx,i));
coefs = polyfit(m1a_riam(idx,i)', m1a_rcam(idx,i)', 1);
slope_m1(i)=coefs(1);
r2_sqm1(i)=mdl.Rsquared.Ordinary;
hold on;ph=plot(mdl,'Color',m1_color);hold on;
ph(2).Color=[0.8 0.8 0.8];
ph(3).Color=[0.8 0.8 0.8];
ph(4).Color=[0.8 0.8 0.8];
end
hold on;title([]);legend('');xlabel('Ipsi Fraction');ylabel('');
set(gca,'FontSize',11);set(gca,'TickDir','out');legend box off;axis square;offsetAxes

cd(save_folder);saveas(gcf, 'linear regression_plots_i_c.pdf');
%% Barplot correlations
%R2
tm1=[];tm2=[];tm3=[];tm1=slope_v1;tm2=slope_s1;tm3=slope_m1;
%plot
fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 500, 170, 250]);set(gcf,'color','w');tiledlayout("horizontal")
b1=bar(1,nanmean(tm1),0.7);hold on;b1.FaceColor=v1_color;b1.EdgeColor='none';b2=bar(2,nanmean(tm2),0.7);hold on;b2.FaceColor=s1_color;b2.EdgeColor='none';
b3=bar(3,nanmean(tm3),0.7);b3.FaceColor=m1_color;b3.EdgeColor='none';
set(b1,'ShowBaseLine','off');set(b2,'ShowBaseLine','off');set(b3,'ShowBaseLine','off');
r=1;rng(1);r1 = r-0.01 + (0.2)*rand(length(tm1),1);
sc1=scatter(r1,tm1,15,'k','filled');hold on;sc1.MarkerEdgeColor=[1 1 1];sc1.MarkerEdgeAlpha=0.5;sc1.MarkerFaceColor=[0.2 0.2 0.2];
 r=2;rng(1);r1 = r-0.01 + (0.2)*rand(length(tm2),1);
 sc1=scatter(r1,tm2,15,'k','filled');hold on;sc1.MarkerEdgeColor=[1 1 1];sc1.MarkerEdgeAlpha=0.5;sc1.MarkerFaceColor=[0.2 0.2 0.2];
 r=3;rng(1);r1 = r-0.01 + (0.2)*rand(length(tm3),1);
 sc1=scatter(r1,tm3,15,'k','filled');hold on;sc1.MarkerEdgeColor=[1 1 1];sc1.MarkerEdgeAlpha=0.5;sc1.MarkerFaceColor=[0.2 0.2 0.2];xticks([1:3]);hold on;
hold on;errorbar([1 2 3],[nanmean(tm1) nanmean(tm2) nanmean(tm3)],[nanstd(tm1)/sqrt(size(tm1,2)) nanstd(tm2)/sqrt(size(tm2,2)) nanstd(tm3)/sqrt(size(tm3,2))]...
    , 'LineStyle', 'none', ... 
        'Color', 'k', 'LineWidth', 1.2);hold on;
xticklabels({'VISp','SSpbf','MOp'});%ylim([-0.05 0.75]);
ylabel('Slope linear fit ipsi - contra');set(gca,'FontSize',11);set(gca,'TickDir','out');box off;xtickangle(45);
 h = gca;h.XAxis.Visible = 'off'
offsetAxes;

  t1=text(0.4,-0.15,'VISp','FontSize',11);set(t1,'Rotation',45);
   t1=text(0.8,-0.28,'SSp-bfd','FontSize',11);set(t1,'Rotation',45);
    t1=text(2.4,-0.15,'MOp','FontSize',11);set(t1,'Rotation',45);

cd(save_folder);saveas(gcf, 'slope_regression.pdf');

 [p,tbl,stats] = anova1([tm1' tm2' tm3'])
 presults = multcompare(stats)
 %% Supp 2
 %% Plot barplot for all regions all areas RETRO (animals=18, areas=45)
 temp1=[];temp1=cell_selecter(data,'virustype',1,'area',1);
 temp2=[];temp2=cell_selecter(data,'virustype',1,'area',2);
 temp3=[];temp3=cell_selecter(data,'virustype',1,'area',3);
retro=[];retro=find(temp1==1);
%call function anatomy_cellnr.m 3)
retro=[];retro=find(temp1==1);
[i_animalv c_animalv i_areas_animalv c_areas_animalv i_areas_animalfv c_areas_animalfv  i_areas_animalcdv c_areas_animalcdv] = anatomy_cellnr(data,retro,cortex_names);
retro=[];retro=find(temp2==1);
[i_animals c_animals i_areas_animals c_areas_animals i_areas_animalfs c_areas_animalfs  i_areas_animalcds c_areas_animalcds] = anatomy_cellnr(data,retro,cortex_names);
retro=[];retro=find(temp3==1);
[i_animalm c_animalm i_areas_animalm c_areas_animalm i_areas_animalfm c_areas_animalfm  i_areas_animalcdm c_areas_animalcdm] = anatomy_cellnr(data,retro,cortex_names);

%% 
temp_color=[v1_color ;s1_color; m1_color];
pan_title={'VISp','SSp-bfd','MOp'};
data_all1={};data_all1={i_animalv i_animals i_animalm};
data_all2={};data_all2={c_animalv c_animals c_animalm};

%means
for j=1:3
 i_animal_ar=[];i_animal_ar=data_all1{j};
 c_animal_ar=[];c_animal_ar=data_all2{j};
means_layers=[nanmean(i_animal_ar(2,:)) nanmean(c_animal_ar(2,:));nanmean(i_animal_ar(3,:)) nanmean(c_animal_ar(3,:));...
    nanmean(i_animal_ar(4,:)) nanmean(c_animal_ar(4,:));nanmean(i_animal_ar(5,:)) nanmean(c_animal_ar(5,:));nanmean(i_animal_ar(6,:)) nanmean(c_animal_ar(6,:))];
semd=sqrt(length(i_animal_ar));
%error
err_layers=[[nanstd(i_animal_ar(2,:))/semd nanstd(c_animal_ar(2,:))/semd];[nanstd(i_animal_ar(3,:))/semd nanstd(c_animal_ar(3,:))/semd];...
    [nanstd(i_animal_ar(4,:))/semd nanstd(c_animal_ar(4,:))/semd];[nanstd(i_animal_ar(5,:))/semd nanstd(c_animal_ar(5,:))/semd] ;[nanstd(i_animal_ar(6,:))/semd nanstd(c_animal_ar(6,:))/semd]];
%plotting starts here 
fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, -100+300*j, 450, 250]);set(gcf,'color','w');
%bars
b=bar(means_layers,0.7);
b(1).FaceColor=[0.8 0.8 0.8];b(2).FaceColor=[0.3 0.3 0.3];b(1).EdgeColor=[1 1 1];b(2).EdgeColor=[1 1 1];hold on;set(b,'ShowBaseLine','off');
%indivdual animals
for r=1:5
rng(1);r1 = r-0.2 + (0.2)*rand(length(retro),1);
sc1=scatter(r1,i_animal_ar(r+1,:),5,'k','filled');hold on;sc1.MarkerEdgeColor=[0.5 0.5 0.5];sc1.MarkerEdgeAlpha=0.5;sc1.MarkerFaceColor=[0.2 0.2 0.2];
end
hold on;
for r=1:5
rng(1);r1 = r+0.1 + (0.2)*rand(length(retro),1);
sc1=scatter(r1,c_animal_ar(r+1,:),5,'k','filled');hold on;sc1.MarkerEdgeColor=[0.5 0.5 0.5];sc1.MarkerEdgeAlpha=0.5;sc1.MarkerFaceColor=[0.2 0.2 0.2];
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
ylabel('Fraction of neurons');set(gca,'FontSize',11);set(gca,'TickDir','out');box off;
%[p1]=signrank(i_animal(6,:),c_animal(6,:))
[u p1]=ttest(i_animal_ar(5,:),c_animal_ar(5,:))
hold on;text(1,0.6,['*'],'FontSize',18);hold on;text(4,0.6,['*'],'FontSize',18);

xlim([0.5 5.5]);
offsetAxes;
h = gca;h.XAxis.Visible = 'off'; 
title(pan_title{j},'Color',temp_color(j,:),'FontWeight','normal');
%xticklabels({'L2/3','L4','L5','L6a','L6b'})
text(0.75,-0.04,'L2/3','FontSize',11);text(1.9,-0.04,'L4','FontSize',11);text(2.9,-0.04,'L5','FontSize',11);text(3.9,-0.04,'L6a','FontSize',11);text(4.9,-0.04,'L6b','FontSize',11);
%hold on;text(0.75,0.6,['**'],'FontSize',18);hold on;text(3.75,0.6,['*'],'FontSize',18);legend({'ipsi','contra'},"Location","northeast");legend boxoff;title('MOp','Color',m1_color)
if j==1 
    legend({'ipsi','contra'},"Location","northeast");legend boxoff;
    cd(save_folder);saveas(gcf, 'layer_global_injection_areaVISP.pdf');
end
if j==2 
    cd(save_folder);saveas(gcf, 'layer_global_injection_areaSSP.pdf');
end
if j==3
    cd(save_folder);saveas(gcf, 'layer_global_injection_areaMOP.pdf');
end
end
%% Supp3 hierarhcy HVA
dat_all1={};dat_all1={i_v1aafm i_s1aafm i_m1aafm};
dat_all2={};dat_all2={c_v1aafm c_s1aafm c_m1aafm};
%idx_modules={frontal_idx lateral_idx somamo_idx visual_idx medial_idx aud_idx};
idx_modules={somamo_idx visual_idx aud_idx frontal_idx lateral_idx medial_idx };
%module_names={'Pfron','Lat','SoMo','Vis','Med','Aud'};
module_names={'SoMo','Vis','Aud','Pfron','Lat','Med',};
pan_title={'VISp','SSp-bfd','MOp'};
temp_color=[v1_color ;s1_color; m1_color];
temp_c1=[0.8 0.8 0.8];
temp_c2=[0.3 0.3 0.3];
ventral_areas=[];dorsal_areas=[];
ventral_areas=[27 28 29 30];
dorsal_areas=[32 33 34 35 36];

p_fig3f=[];
%plot
%fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 500, 250, 330]);set(gcf,'color','w');t=tiledlayout('vertical');
%t.TileSpacing = 'compact';t.Padding = 'compact';

rm1=[];rm2=[];rm1=dat_all1{1};rm2=dat_all2{1};
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
visp_iln_i=i_metric(:,:,1);visp_hindex_i=i_metric(:,:,4);

vareas_names=cortex_names(27:36,2)
kk1=[];kk2=[];
[ventral_iln kk1]=sort(nanmean(visp_iln_i(:,ventral_areas)),'descend');
kk1=[kk1 5];
[dorsal_iln kk2]=sort(nanmean(visp_iln_i(:,dorsal_areas)),'ascend');
yerr1=[];yerr1=nanstd(visp_iln_i(:,ventral_areas))./sqrt(size(visp_iln_i(:,ventral_areas),1));
yerr1=[yerr1 NaN];
yerr2=[];yerr2=nanstd(visp_iln_i(:,dorsal_areas))./sqrt(size(visp_iln_i(:,dorsal_areas),1));
fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 500, 300, 330]);set(gcf,'color','w')
 hold on;h=errorbar([1:10],[ventral_iln NaN dorsal_iln],[yerr1(kk1) yerr2(kk2)], 'LineStyle', 'none','Color', 'k', 'LineWidth', 1,'CapSize',0);
hold on;sc1=scatter([1:10],[ventral_iln NaN dorsal_iln] ,60,'k','filled');hold on;sc1.MarkerEdgeColor=[1 1 1];sc1.MarkerEdgeAlpha=1;sc1.MarkerFaceColor=[0.65 0.65 0.65];
xticks([1:1:10]);xticklabels(vareas_names([kk1 kk2+5]));
ylabel('ILN ipsi');set(gca,'FontSize',11);set(gca,'TickDir','out');box off;%xtickangle(45);
line([5 5],[0.4 0.9],'LineStyle','--','Color','k');
text(2,0.9,'ventral','FontSize',11);
text(7,0.9,'dorsal','FontSize',11);
offsetAxes;
cd(save_folder);saveas(gcf, 'ILN_HVA_sorted.pdf');
%% Correlation for HVAs between ILN and H-index
%Calculate correaltions across Ns
r_a=[];p_a=[];
for i=1:size(visp_iln_i,1)
r=[];p=[];
[r p] = corr(visp_iln_i(i,[ventral_areas dorsal_areas])',visp_hindex_i(i,[ventral_areas dorsal_areas])','Type','Spearman','Rows','complete');
r_a(i)=r(1);
p_a(i)=p(1);
end
xerr=[];xerr=nanstd(visp_hindex_i(:,[ventral_areas dorsal_areas]))/sqrt(size(visp_hindex_i,1));
yerr=[];yerr=nanstd(visp_iln_i(:,[ventral_areas dorsal_areas]))/sqrt(size(visp_iln_i,1));
fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 500, 250, 250]);set(gcf,'color','w')
hold on;errorbar(nanmean(visp_hindex_i(:,[ventral_areas dorsal_areas])),nanmean(visp_iln_i(:,[ventral_areas dorsal_areas])),yerr,yerr,xerr,xerr, 'LineStyle', 'none', 'Color', 'k', 'LineWidth', 0.5,'CapSize',0);hold on;
sc1=scatter(nanmean(visp_hindex_i(:,[ventral_areas dorsal_areas])),nanmean(visp_iln_i(:,[ventral_areas dorsal_areas])) ,60,'k','filled');sc1.MarkerEdgeColor=[1 1 1];sc1.MarkerEdgeAlpha=1;sc1.MarkerFaceColor=[0.65 0.65 0.65];
xlabel('h-index ipsi');ylabel('ILN ipsi');set(gca,'FontSize',11);set(gca,'TickDir','out');box off;%xtickangle(45);
text(0.15,0.85,['r= ' num2str(round(nanmean(r_a),2)) ' +- ' num2str(round(nanstd(r_a)/sqrt(length(r_a)),2))],'FontSize',11);
text(0.22,0.7,'VISrl','FontSize',11)
offsetAxes;
cd(save_folder);saveas(gcf, 'correlation_HVA_ILN_h-index.pdf');
%% heatmap
dat_all1={};dat_all1={i_v1aam i_s1aam i_m1aam};
rm1=[];rm2=[];rm1=dat_all1{1};rm2=dat_all2{1};
tm_hva=nanmean(rm1(2:6,27:36,:),3);

fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 1000-250*j, 350, 250]);set(gcf,'color','w');
nexttile
h=imagesc(tm_hva(:,[kk1 kk2+5]));
set(h, 'AlphaData', ~isnan(tm_hva(:,[kk1 kk2+5])))
cmap(temp_color(1,:),100,5,5);a1=colorbar;
ylabel(a1,'Fraction of neuron','FontSize',11,'Rotation',270);
[hh oo]=find(isnan(tm_hva(:,[kk1 kk2+5])));
coord_b=[oo hh];
for u=1:length(coord_b)
    hold on
    text(coord_b(u,1)-0.25,coord_b(u,2)-0.1,'x','FontSize',14,'Color','k')
end
box off;set(gca,'TickDir','out');
hold on;yticklabels({'L2/3','L4','L5','L6a','L6b'});
xticks([1:1:10]);xticklabels(vareas_names([kk1 kk2+5]));
cd(save_folder);saveas(gcf, 'ILNHVA_heatmap.pdf');
%offsetAxes;
%% Supp 4

%% 
temp_c1='b';
temp_c2='b';
hold on;
temp_rm1=[];temp_rm1=nanmean(aILNi_l6);
temp_rm2=[];temp_rm2=nanmean(aILNc_l6);


yerr=[];yerr=nanstd(aILNi_l6);
yerr2=[];yerr2=nanstd(aILNc_l6);
%yerr3=[];yerr3=nanstd(rm2(:,cl_idx(k)))/sqrt(length(rm2(:,cl_idx(k))));
sort_d=[];kk=[];
[sort_d kk]=sort(temp_rm1,'ascend');
hold on;h=errorbar([1:45],[sort_d],[yerr/sqrt(length(sort_d(~isnan(sort_d))))]...
    , 'LineStyle', 'none', ... 
        'Color', temp_c1, 'LineWidth', 0.5,'CapSize',2);set([h.Bar, h.Line], 'ColorType', 'truecoloralpha', 'ColorData', [h.Line.ColorData(1:3); 255*alpha1]);
hold on;h=errorbar([1:45],[temp_rm2(kk)],[yerr2(kk)/sqrt(length(temp_rm2(~isnan(temp_rm2))))]...
    , 'LineStyle', 'none', ... 
        'Color', temp_c2, 'LineWidth', 0.5,'CapSize',2);

hold on;pp1=scatter(1:45,sort_d,30,'o');pp1.MarkerEdgeColor=[1 1 1];pp1.MarkerFaceAlpha=1;pp1.MarkerFaceColor=temp_c1;
hold on;pp1=scatter(1:45,temp_rm2(kk),30,'o');pp1.MarkerEdgeColor=[1 1 1];pp1.MarkerFaceAlpha=1;pp1.MarkerFaceColor=temp_c2;box off;
xticks([1:45]);ylim([0 1])
%xtickangle(45);
%% Plot fraction across the modules Supp 4
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











