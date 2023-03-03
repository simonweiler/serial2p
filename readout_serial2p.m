%serial 2p data TRACING

%% Set directories and experimentator
data_dir    = 'W:\margrie\SimonWeiler\AnalyzedData\Tracing_Imaging\serial2p\summaries';%data directory of raw data;change accordingly
save_dir        = 'W:\margrie\SimonWeiler\AnalyzedData\Neuropixels\matfiles\';
ExpXls            = 'C:\Users\simonw\S1-V1 interaction Dropbox\Simon Weiler\Callosal_L6\serial2p_structure\serial2p_readout.xlsx'

batchopt          = parseExperimentsXls_serial2p(ExpXls);%calls the nested function parseExperimentsXls_ephys and considers the user flag (1 or 0)
nummice           = length(batchopt.mouse);%length of experiments to be analyzed

%% Set brain areas (layers) to look at
%use the nomenclature of ABI
list_to_find={};
areas_to_find={};
areas_to_find={'Primary visual area','posteromedial visual area', 'Lateral visual area', 'Rostrolateral area',...
    'Primary somatosensory area, barrel field', 'Primary somatosensory area, lower limb','Primary somatosensory area, upper limb','Primary somatosensory area, nose','Supplemental somatosensory area'...,
    'Primary motor area','Secondary motor area', 'Anterior cingulate area, ventral part','Anterior cingulate area, dorsal part'...
    'Primary auditory area','Retrosplenial area','Visceral area','Agranular insular area'...
    ,'Temporal association areas','Ectorhinal area', 'Entorhinal area', 'Perirhinal area', 'Postsubiculum'};

%Abbreviation for display
areas_abb={'VISp','VISpm','VISal','VISrl',...
    'SSbf','SSll','SSul','SSn','SSs',...
    'MOp','MOs','ACv','ACd'...
    'AUDp','RSP','VISC','Alp',...
    'TEa','ECT','ENTI','PERI','POST'};
% 
% list_to_find={'Primary visual area, layer 1','Primary visual area, layer 2/3',...
%         'Primary visual area, layer 4','Primary visual area, layer 5','Primary visual area, layer 6a',...
%         'Primary visual area, layer 6b'};
% lb={'VISp'};

% list_to_find={'Primary auditory area, layer 1','Primary auditory area, layer 2/3',...
%         'Primary auditory area, layer 4','Primary auditory area, layer 5','Primary auditory area, layer 6a',...
%         'Primary auditory area, layer 6b'}; 
% lb={'AUDp'};

% list_to_find={'Lateral visual area, layer 1','Lateral visual area, layer 2/3',...
%         'Lateral visual area, layer 4','Lateral visual area, layer 5','Lateral visual area, layer 6a',...
%         'Lateral visual area, layer 6b'}; 
% lb={'VISal'};

% list_to_find={'Temporal association areas, layer 1','Temporal association areas, layer 2/3',...
%         'Temporal association areas, layer 4','Temporal association areas, layer 5','Temporal association areas, layer 6a',...
%         'Temporal association areas, layer 6b'}; 
% lb={'TEa'};

% list_to_find={'Retrosplenial area, dorsal part, layer 1','Retrosplenial area, dorsal part, layer 2/3',...
%         'Retrosplenial area, dorsal part, layer 5','Retrosplenial area, dorsal part, layer 6a',...
%         'Retrosplenial area, dorsal part, layer 6b'}; 
% lb={'RSPd'};

% list_to_find={'Ectorhinal area/Layer 1','Ectorhinal area/Layer 2/3',...
%         'Ectorhinal area/Layer 5','Ectorhinal area/Layer 6a',...
%         'Ectorhinal area/Layer 6b'}; 
% lb={'ECT'};
% 
% list_to_find={'Primary somatosensory area, barrel field, layer 1','Primary somatosensory area, barrel field, layer 2/3',...
%         'Primary somatosensory area, barrel field, layer 4','Primary somatosensory area, barrel field, layer 5',...
%         'Primary somatosensory area, barrel field, layer 6a','Primary somatosensory area, barrel field, layer 6b'}; 
% lb={'SSp-bf'};

% list_to_find={'Primary motor area, Layer 1','Primary motor area, Layer 2/3',...
%         'Primary motor area, Layer 5','Primary motor area, Layer 6a',...
%         'Primary motor area, Layer 6b'}; 
% lb={'MOp'};


% list_to_find={'Secondary motor area, layer 1','Secondary motor area, layer 2/3',...
%         'Secondary motor area, layer 5','Secondary motor area, layer 6a',...
%         'Secondary motor area, layer 6b'}; 
% lb={'MOs'};

% list_to_find={'Supplemental somatosensory area, layer 1','Supplemental somatosensory area, layer 2/3',...
%         'Supplemental somatosensory area, layer 4','Supplemental somatosensory area, layer 5',...
%         'Supplemental somatosensory area, layer 6a','Supplemental somatosensory area, layer 6b'}; 
% lb={'MOs'};

%homotopic callsoal areas 
 list_to_find={'Primary visual area, layer 1','Primary visual area, layer 2/3','Primary visual area, layer 4','Primary visual area, layer 5','Primary visual area, layer 6a','Primary visual area, layer 6b',...
         'Primary somatosensory area, barrel field, layer 1','Primary somatosensory area, barrel field, layer 2/3','Primary somatosensory area, barrel field, layer 4','Primary somatosensory area, barrel field, layer 5','Primary somatosensory area, barrel field, layer 6a','Primary somatosensory area, barrel field, layer 6b',...
         'Primary motor area, Layer 1','Primary motor area, Layer 2/3','Primary motor area, Layer 5','Primary motor area, Layer 6a','Primary motor area, Layer 6b',...
         'posteromedial visual area, layer 1','posteromedial visual area, layer 2/3','posteromedial visual area, layer 4','posteromedial visual area, layer 5','posteromedial visual area, layer 6a','posteromedial visual area, layer 6b',...
         'Lateral visual area, layer 1','Lateral visual area, layer 2/3','Lateral visual area, layer 4','Lateral visual area, layer 5','Lateral visual area, layer 6a','Lateral visual area, layer 6b',...
         'Rostrolateral area, layer 1','Rostrolateral area, layer 2/3','Rostrolateral area, layer 4','Rostrolateral area, layer 5','Rostrolateral area, layer 6a','Rostrolateral area, layer 6b',...
         'Primary somatosensory area, lower limb, layer 1','Primary somatosensory area, lower limb, layer 2/3','Primary somatosensory area, lower limb, layer 4','Primary somatosensory area, lower limb, layer 5','Primary somatosensory area, lower limb, layer 6a','Primary somatosensory area, lower limb, layer 6b',...
         'Primary somatosensory area, upper limb, layer 1','Primary somatosensory area, upper limb, layer 2/3','Primary somatosensory area, upper limb, layer 4','Primary somatosensory area, upper limb, layer 5','Primary somatosensory area, upper limb, layer 6a','Primary somatosensory area, upper limb, layer 6b',...
         'Primary somatosensory area, nose, layer 1','Primary somatosensory area, nose, layer 2/3','Primary somatosensory area, nose, layer 4','Primary somatosensory area, nose, layer 5','Primary somatosensory area, nose, layer 6a','Primary somatosensory area, nose, layer 6b',...
         'Supplemental somatosensory area, layer 1','Supplemental somatosensory area, layer 2/3','Supplemental somatosensory area, layer 4','Supplemental somatosensory area, layer 5','Supplemental somatosensory area, layer 6a','Supplemental somatosensory area, layer 6b',...
         'Secondary motor area, layer 1','Secondary motor area, layer 2/3','Secondary motor area, layer 5','Secondary motor area, layer 6a', 'Secondary motor area, layer 6b',...
         'Anterior cingulate area, dorsal part, layer 1','Anterior cingulate area, dorsal part, layer 2/3','Anterior cingulate area, dorsal part, layer 5','Anterior cingulate area, dorsal part, layer 6a','Anterior cingulate area, dorsal part, layer 6b',...
         'Primary auditory area, layer 1','Primary auditory area, layer 2/3','Primary auditory area, layer 4','Primary auditory area, layer 5','Primary auditory area, layer 6a','Primary auditory area, layer 6b',...
         'Retrosplenial area, dorsal part, layer 1','Retrosplenial area, dorsal part, layer 2/3', 'Retrosplenial area, dorsal part, layer 5','Retrosplenial area, dorsal part, layer 6a','Retrosplenial area, dorsal part, layer 6b',...
         'Visceral area, layer 1','Visceral area, layer 2/3','Visceral area, layer 4','Visceral area, layer 5','Visceral area, layer 6a','Visceral area, layer 6b',...
         'Temporal association areas, layer 1','Temporal association areas, layer 2/3','Temporal association areas, layer 4','Temporal association areas, layer 5','Temporal association areas, layer 6a','Temporal association areas, layer 6b',...
         'Ectorhinal area/Layer 1','Ectorhinal area/Layer 2/3','Ectorhinal area/Layer 5','Ectorhinal area/Layer 6a','Ectorhinal area/Layer 6b'};
 list_to_find_layers=[];
 list_to_find_layers=[6 6 5 6 6 6 6 6 6 6 5 5 6 5 6 6 5];
 
 
 lb={'VISp'};

inject_area=[];inject_area=cell2mat(batchopt.injectarea);
visp_layers_all=[];
for i=1:nummice%for loop over experiments across days
     sum_dir=fullfile(data_dir,char(batchopt.mouseID{i}),'summary.csv');
     temp=[];
    temp = readtable(sum_dir);
     outs2={};
    outs2=table2cell(temp);
    hemi=[];hemi=cell2mat(batchopt.hemi_ids(i));
   list_all_areas={outs2{:,1}};   
   idx_areas=[];
   for k=1:length(list_to_find);
   idx_areas(:,k) = ismember(list_all_areas,list_to_find{k});
   end
   idx_brainarea={};
   count_brainareas=[];
   for m=1:length(areas_to_find);
    idx_brainarea{m}=find(contains(list_all_areas,areas_to_find{m})==1)';
    idx_temp=[];
    idx_temp=find(contains(list_all_areas,areas_to_find{m})==1)';
    if hemi==0;
    count_brainareas(:,m)=sum([outs2{idx_temp,3}]');
    count_brain_volume(:,m)=sum([outs2{idx_temp,9}]');
    else
    count_brainareas(:,m)=sum([outs2{idx_temp,2}]');
    count_brain_volume(:,m)=sum([outs2{idx_temp,8}]');
    end
   end

   idx_r=[]; idx_c=[];
   [idx_r idx_c] = find(idx_areas==1);
   visp_layers=[];
 
   
    if hemi==0;
    visp_layers=[outs2{idx_r,3}]';
    else hemi==1
    visp_layers=[outs2{idx_r,2}]';  
    end
    visp_layers_all(:,i)=visp_layers;
    areas_all(:,i)=count_brainareas';
   volume_all(:,i)=count_brain_volume';

end
%% Plot barplots fraction 
%all areasraw count and then fraction
data=[];
data=areas_all./sum(areas_all);
subd=find(inject_area>-1);
fig1= figure;set(fig1, 'Name', 'Barplot groups');set(fig1, 'Position', [400, 500,500, 300]);set(gcf,'color','w');
for i=1:size(data,1)
hold on;
b2=bar(i,nanmean(data(i,subd)));b2.FaceColor=[0.7 0.7 0.7];
plot(ones(1,length(data(i,subd)))*i,data(i,subd),'ko','MarkerEdgeColor',[0,0,0],'MarkerSize',3);
hold on;
er=errorbar(i,nanmean(data(i,subd)),nanstd(data(i,subd))/sqrt(size(data(i,subd),2)));
er.Color = [0 0 0];er.LineWidth=1;er.LineStyle = 'none'; hold on;
end
xticks(1:size(data,1));
xticklabels(areas_abb);
ylabel('Fraction of cells');
title('Contra Hemisphere');
ylim([0 0.6]);
set(gca,'FontSize',10);
%% all areasraw cell per volume and then fraction
data=[];
data=volume_all./sum(volume_all);
subd=find(inject_area>-1);
fig1= figure;set(fig1, 'Name', 'Barplot groups');set(fig1, 'Position', [400, 500, 500, 300]);set(gcf,'color','w');
for i=1:size(data,1)
hold on;
b2=bar(i,nanmean(data(i,subd)));b2.FaceColor=[0.7 0.7 0.7];
plot(ones(1,length(data(i,subd)))*i,data(i,subd),'ko','MarkerEdgeColor',[0,0,0],'MarkerSize',3);
hold on;
er=errorbar(i,nanmean(data(i,subd)),nanstd(data(i,subd))/sqrt(size(data(i,subd),2)));
er.Color = [0 0 0];er.LineWidth=1;er.LineStyle = 'none'; hold on;
end
xticks(1:size(data,1));
xticklabels(areas_abb);
ylabel('Fraction of cells (norm. per Volume)');
title('Contra Hemisphere');
ylim([0 0.4]);
set(gca,'FontSize',10);
%% Layer per area 
data=[];
data=visp_layers_all./sum(visp_layers_all);
subd=find(inject_area>-1);
fig1= figure;set(fig1, 'Name', 'Barplot groups');set(fig1, 'Position', [400, 500, 300, 200]);set(gcf,'color','w');
for i=1:size(data,1)
hold on;
b2=bar(i,nanmean(data(i,subd)));b2.FaceColor=[0.7 0.7 0.7];
plot(ones(1,length(data(i,subd)))*i,data(i,subd),'ko','MarkerEdgeColor',[0,0,0],'MarkerSize',3);
hold on;
er=errorbar(i,nanmean(data(i,subd)),nanstd(data(i,subd))/sqrt(size(data(i,subd),2)));
er.Color = [0 0 0];er.LineWidth=1;er.LineStyle = 'none'; hold on;
end
%xticks(1:size(data,1));xticklabels({'L1','L2/3','L4','L5','L6a','L6b'});
xticks(1:size(data,1));xticklabels({'L1','L2/3','L5','L6a','L6b'});
ylabel('Fraction of cells');
set(gca,'FontSize',10);
title(['Contra Hemisphere' lb]);
ylim([0 0.8]);

%% Plot all different V1, S1 and M1 RETROGRADE in one figure: USE THIS for now
subvis=[];subss=[];subm=[];
subvis=find(cell2mat(batchopt.type)==1 & inject_area<2);
subss=find(cell2mat(batchopt.type)==1 & inject_area==2);
subm=find(cell2mat(batchopt.type)==1 & inject_area==3);
data=[];
data=areas_all./sum(areas_all);
maxyline=0.6;
fig1= figure;set(fig1, 'Name', 'Barplot groups');set(fig1, 'Position', [400, 500,700, 400]);set(gcf,'color','w');
dat_bar=[];dat_barerr=[];
for i=1:size(data,1)
dat_bar(i,:)=[nanmean(data(i,subvis)) nanmean(data(i,subss)) nanmean(data(i,subm))];
dat_barerr(i,:)=[nanstd(data(i,subvis))/sqrt(size(data(i,subvis),2))' nanstd(data(i,subss))/sqrt(size(data(i,subss),2))'...
     nanstd(data(i,subm)/sqrt(size(data(i,subm),2))')];
end
b = bar(dat_bar, 'grouped');b(1).FaceColor=[0.5 0.5 0.5];b(2).FaceColor='m';b(3).FaceColor=[0 0.7 1];
hold on;
% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(dat_bar);
% Get the x coordinate of the bars
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end
% Plot the errorbars
er=errorbar(x',dat_bar,dat_barerr,'k','linestyle','none');
hold off
%er.Color = [0 0 0];er.LineWidth=1;er.LineStyle = 'none'; hold on;
xticks(1:size(data,1));
xticklabels(areas_abb);
ylabel('Fraction of cells');
title('Contralateral origin of CPNs targeting VISp, SSbf and MOp');
ylim([0 maxyline]);
set(gca,'FontSize',10);box off;
hold on;line([4.5 4.5],[0 maxyline],'linewidth',0.5,'Color','k','LineStyle','--');
hold on;line([9.5 9.5],[0 maxyline],'linewidth',0.5,'Color','k','LineStyle','--');
hold on;line([11.5 11.5],[0 maxyline],'linewidth',0.5,'Color','k','LineStyle','--');
leg=legend({['VISp ' 'N=' num2str(length(subvis))], ['SSbf ' 'N=' num2str(length(subss))],['MOp ' 'N=' num2str(length(subm))]});
legend boxoff;set(gca,'FontSize',12);   
%% %% Layer per area for homotopic contra areas comparison: VISp, SSbf, MOp RETROGRADE
datavis=[];datass=[];datam=[];
datavis=visp_layers_all(1:6,subvis)./sum(visp_layers_all(1:6,subvis));
datass=visp_layers_all(7:12,subss)./sum(visp_layers_all(7:12,subss));
datam=visp_layers_all(13:17,subm)./sum(visp_layers_all(13:17,subm));
tempo=nanmean(datam,2);datam_mod=[tempo(1) tempo(2) 0 tempo(3) tempo(4) tempo(5)]';
tempo2=nanstd(datam,[],2)./sqrt(size(datam,1));datam_sem=[tempo2(1) tempo2(2) 0 tempo2(3) tempo2(4) tempo2(5)]';

data_all=[];data_all=[nanmean(datavis,2) nanmean(datass,2) datam_mod];
dat_barerr=[]; dat_barerr=[nanstd(datavis,[],2)./sqrt(size(datavis,1)) nanstd(datass,[],2)./sqrt(size(datass,1))...
    datam_sem];

fig1= figure;set(fig1, 'Name', 'Barplot groups');set(fig1, 'Position', [400, 500, 350, 350]);set(gcf,'color','w');
b = bar(data_all, 'grouped');b(1).FaceColor=[0.5 0.5 0.5];b(2).FaceColor='m';b(3).FaceColor=[0 0.7 1];
hold on;
% Calculate the number of groups and number of bars in each group
[ngroups]=[];[nbars]=[];
[ngroups,nbars] = size(data_all);
% Get the x coordinate of the bars
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end
% Plot the errorbars
er=errorbar(x',data_all,dat_barerr,'k','linestyle','none');
xticks(1:size(data,1));xticklabels({'L1','L2/3','L4','L5','L6a','L6b'});
ylabel('Fraction of cells');
set(gca,'FontSize',10);
title({'Homotopic layer origin of CPNs','to VISp, SSbf and MOp'});
ylim([0 0.7]);box off;
leg=legend({['VISp ' 'N=' num2str(length(subvis))], ['SSbf ' 'N=' num2str(length(subss))],['MOp ' 'N=' num2str(length(subm))]});
legend boxoff;set(gca,'FontSize',12);   
%% Plot in heatmap across all areas
%VISp injeciton 
datavis=[];datass=[];datam=[];endpoints=[];startpoints=[];all_vis=[];all_ss=[];all_m=[];
%threshold for retro
visp_layers_all_sub=[];
visp_layers_all_sub=visp_layers_all;
thrsh=30;
[r m] = find(visp_layers_all_sub<thrsh);
for i=1:length(r)
visp_layers_all_sub(r(i),m(i))=0;
end

find(visp_layers_all>30)
endpoints=cumsum(list_to_find_layers);
startpoints=[1 endpoints(1:end-1)+1];
for i=1:length(endpoints)
    visp_layers_all_sub(startpoints(i))=0;
     d_vis=[];  d_ss=[];  d_m=[];
    datavis=visp_layers_all_sub(startpoints(i):endpoints(i),subvis)./sum(visp_layers_all_sub(startpoints(i):endpoints(i),subvis));
    datass=visp_layers_all_sub(startpoints(i):endpoints(i),subss)./sum(visp_layers_all_sub(startpoints(i):endpoints(i),subss));
    datam=visp_layers_all_sub(startpoints(i):endpoints(i),subm)./sum(visp_layers_all_sub(startpoints(i):endpoints(i),subm));
    d_vis=nanmean(datavis,2);
    d_ss=nanmean(datass,2);
    d_m=nanmean(datam,2);
    if length(d_vis)<6
        f_vis=[d_vis(1:2); 0; d_vis(3:5)];
        f_ss=[d_ss(1:2); 0; d_ss(3:5)];
        f_m=[d_m(1:2); 0; d_m(3:5)];
    else
        f_vis=d_vis;
        f_ss=d_ss;
        f_m=d_m;
    end
  
    all_vis(:,i)=f_vis;
    all_ss(:,i)=f_ss;
    all_m(:,i)=f_m;
end

cmap1='plasma';
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 200 ,600, 300]);
imagesc(all_vis);colormap(cmap1);
yticklabels({'L1','2/3','L4','L5','L6a','L6b'})
xticks(1:1:17)
%xticklabels([areas_abb(1:11) areas_abb(13:16) areas_abb(18:19)]);
xticklabels([areas_abb([1 5 10]) areas_abb(2:4) areas_abb(6:9) areas_abb(11) areas_abb(13:16) areas_abb(18:19)]);
cmap1='plasma';
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 200 ,600, 300]);
imagesc(all_ss);colormap(cmap1);
yticklabels({'L1','L2/3','L4','L5','L6a','L6b'})
xticks(1:1:17)
%xticklabels([areas_abb(1:11) areas_abb(13:16) areas_abb(18:19)]);
xticklabels([areas_abb([1 5 10]) areas_abb(2:4) areas_abb(6:9) areas_abb(11) areas_abb(13:16) areas_abb(18:19)]);
cmap1='plasma';
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 200 ,600, 300]);
imagesc(all_m);colormap(cmap1);
yticklabels({'L1','L2/3','L4','L5','L6a','L6b'})
xticks(1:1:17)
%xticklabels([areas_abb(1:11) areas_abb(13:16) areas_abb(18:19)]);
xticklabels([areas_abb([1 5 10]) areas_abb(2:4) areas_abb(6:9) areas_abb(11) areas_abb(13:16) areas_abb(18:19)]);
%% Plot all different V1, S1 and M1 ANTEROGRADE in one figure
subvis=[];subss=[];subm=[];
subvis=find(cell2mat(batchopt.type)==2 & inject_area<2);
subss=find(cell2mat(batchopt.type)==2 & inject_area==2);
subm=find(cell2mat(batchopt.type)==2 & inject_area==3);
data=[];
data=areas_all./sum(areas_all);
maxyline=0.6;
fig1= figure;set(fig1, 'Name', 'Barplot groups');set(fig1, 'Position', [400, 500,700, 400]);set(gcf,'color','w');
dat_bar=[];dat_barerr=[];
for i=1:size(data,1)
dat_bar(i,:)=[nanmean(data(i,subvis)) nanmean(data(i,subss)) nanmean(data(i,subm))];
dat_barerr(i,:)=[nanstd(data(i,subvis))/sqrt(size(data(i,subvis),2))' nanstd(data(i,subss))/sqrt(size(data(i,subss),2))'...
     nanstd(data(i,subm)/sqrt(size(data(i,subm),2))')];
end
b = bar(dat_bar, 'grouped');b(1).FaceColor=[0.5 0.5 0.5];b(2).FaceColor='m';b(3).FaceColor=[0 0.7 1];
hold on;
% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(dat_bar);
% Get the x coordinate of the bars
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end
% Plot the errorbars
er=errorbar(x',dat_bar,dat_barerr,'k','linestyle','none');
hold off
%er.Color = [0 0 0];er.LineWidth=1;er.LineStyle = 'none'; hold on;
xticks(1:size(data,1));
xticklabels(areas_abb);
ylabel('Fraction of cells');
title('Contralateral projections targets of CPNs orginating from VISp, SSbf and MOp');
ylim([0 0.6]);
hold on;line([4.5 4.5],[0 maxyline],'linewidth',0.5,'Color','k','LineStyle','--');
hold on;line([9.5 9.5],[0 maxyline],'linewidth',0.5,'Color','k','LineStyle','--');
hold on;line([11.5 11.5],[0 maxyline],'linewidth',0.5,'Color','k','LineStyle','--');
leg=legend({['VISp ' 'N=' num2str(length(subvis))], ['SSbf ' 'N=' num2str(length(subss))],['MOp ' 'N=' num2str(length(subm))]});
legend boxoff;set(gca,'FontSize',12);  box off;
%% 
%% %% Layer per area for homtopic contra areas comparison ANTEROGRADE
datavis=[];datass=[];datam=[];
datavis=visp_layers_all(1:6,subvis)./sum(visp_layers_all(1:6,subvis));
datass=visp_layers_all(7:12,subss)./sum(visp_layers_all(7:12,subss));
datam=visp_layers_all(13:17,subm)./sum(visp_layers_all(13:17,subm));
tempo=nanmean(datam,2);datam_mod=[tempo(1) tempo(2) 0 tempo(3) tempo(4) tempo(5)]';
tempo2=nanstd(datam,[],2)./sqrt(size(datam,1));datam_sem=[tempo2(1) tempo2(2) 0 tempo2(3) tempo2(4) tempo2(5)]';

data_all=[];data_all=[nanmean(datavis,2) nanmean(datass,2) datam_mod];
dat_barerr=[]; dat_barerr=[nanstd(datavis,[],2)./sqrt(size(datavis,1)) nanstd(datass,[],2)./sqrt(size(datass,1))...
    datam_sem];

fig1= figure;set(fig1, 'Name', 'Barplot groups');set(fig1, 'Position', [400, 500, 350, 350]);set(gcf,'color','w');
b = bar(data_all, 'grouped');b(1).FaceColor=[0.5 0.5 0.5];b(2).FaceColor='m';b(3).FaceColor=[0 0.7 1];
hold on;
% Calculate the number of groups and number of bars in each group
[ngroups]=[];[nbars]=[];
[ngroups,nbars] = size(data_all);
% Get the x coordinate of the bars
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end
% Plot the errorbars
er=errorbar(x',data_all,dat_barerr,'k','linestyle','none');
xticks(1:size(data,1));xticklabels({'L1','L2/3','L4','L5','L6a','L6b'});
ylabel('Fraction of cells');
set(gca,'FontSize',10);
title(['Homotopic Contra Hemisphere']);
ylim([0 0.7]);box off;
title({'Homotopic layer projection targets of CPNs','from VISp, SSbf and MOp'});
leg=legend({['VISp ' 'N=' num2str(length(subvis))], ['SSbf ' 'N=' num2str(length(subss))],['MOp ' 'N=' num2str(length(subm))]});
legend boxoff;set(gca,'FontSize',12);   
%% 



%% COMPARE ANTERO vs RETRO
subretro=[];subantero=[];
subretro=find(cell2mat(batchopt.type)==1 & inject_area<2);
subantero=find(cell2mat(batchopt.type)==2 & inject_area<2);
% subretro=find(cell2mat(batchopt.type)==1 & inject_area==2);
% subantero=find(cell2mat(batchopt.type)==2 & inject_area==2);
data=[];
data=areas_all./sum(areas_all);

fig1= figure;set(fig1, 'Name', 'Barplot groups');set(fig1, 'Position', [400, 500,700, 400]);set(gcf,'color','w');
dat_bar=[];dat_barerr=[];
for i=1:size(data,1)
dat_bar(i,:)=[nanmean(data(i,subretro)) nanmean(data(i,subantero))];
dat_barerr(i,:)=[nanstd(data(i,subretro))/sqrt(size(data(i,subretro),2))' nanstd(data(i,subantero))/sqrt(size(data(i,subantero),2))'];
end
b = bar(dat_bar, 'grouped');b(1).FaceColor=[0.5 0.5 0.5];b(2).FaceColor='m';
hold on;
% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(dat_bar);
% Get the x coordinate of the bars
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end
% Plot the errorbars
er=errorbar(x',dat_bar,dat_barerr,'k','linestyle','none');
hold off
%er.Color = [0 0 0];er.LineWidth=1;er.LineStyle = 'none'; hold on;
xticks(1:size(data,1));
xticklabels(areas_abb);
ylabel('Fraction of cells');
title(['Contra Hemisphere' lb]);
ylim([0 0.6]);
set(gca,'FontSize',10);box off;
legend({'retro', 'antero'});
%% Show antero vs retro per layer
data=[];
data=visp_layers_all./sum(visp_layers_all);
subretro=[];subantero=[];
subretro=find(cell2mat(batchopt.type)==1 & inject_area<2);
subantero=find(cell2mat(batchopt.type)==2 & inject_area<2);
% subretro=find(cell2mat(batchopt.type)==1 & inject_area==2);
% subantero=find(cell2mat(batchopt.type)==2 & inject_area==2);

fig1= figure;set(fig1, 'Name', 'Barplot groups');set(fig1, 'Position', [400, 500,700, 400]);set(gcf,'color','w');
dat_bar=[];dat_barerr=[];
for i=1:size(data,1)
dat_bar(i,:)=[nanmean(data(i,subretro)) nanmean(data(i,subantero))];
dat_barerr(i,:)=[nanstd(data(i,subretro))/sqrt(size(data(i,subretro),2))' nanstd(data(i,subantero))/sqrt(size(data(i,subantero),2))'];
end

b = bar(dat_bar, 'grouped');b(1).FaceColor=[0.5 0.5 0.5];b(2).FaceColor='m';
hold on;
% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(dat_bar);
% Get the x coordinate of the bars
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end
% Plot the errorbars
er=errorbar(x',dat_bar,dat_barerr,'k','linestyle','none');
hold off
%er.Color = [0 0 0];er.LineWidth=1;er.LineStyle = 'none'; hold on;
xticks(1:size(data,1));xticklabels({'L1','L2/3','L4','L5','L6a','L6b'});
ylabel('Fraction of cells');
set(gca,'FontSize',10);
title(['Homotopic Contra Hemisphere']);
ylim([0 0.9]);box off;
legend({'retro', 'antero'});