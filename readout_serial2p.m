%serial 2p data TRACING

%% Set directories and experimentator
data_dir    = 'W:\margrie\SimonWeiler\AnalyzedData\Tracing_Imaging\serial2p\summaries';%data directory of raw data;change accordingly
save_dir        = 'W:\margrie\SimonWeiler\AnalyzedData\Neuropixels\matfiles\';
ExpXls            = 'C:\Users\simonw\S1-V1 interaction Dropbox\Simon Weiler\Callosal_L6\serial2p_structure\serial2p_readout.xlsx'

batchopt          = parseExperimentsXls_serial2p(ExpXls);%calls the nested function parseExperimentsXls_ephys and considers the user flag (1 or 0)
nummice           = length(batchopt.mouse);%length of experiments to be analyzed

%% Set brain areas (layers) to look at
list_to_find={};
areas_to_find={};
areas_to_find={'Primary visual area','posteromedial visual area', 'Lateral visual area', 'Rostrolateral area','Retrosplenial area',...
    'Primary auditory area', 'Primary somatosensory area, barrel field','Primary motor area','Secondary motor area','Temporal association areas',...
    'Anterior cingulate area','Ectorhinal area', 'Entorhinal area', 'Perirhinal area', 'Postsubiculum'};
areas_abb={'VISp','VISpm','VISal','VISrl','RSP','AUDp','SSbf','MOp','MOs','TEa','AC','ECT','ENTI','PERI','POST'};

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

list_to_find={'Retrosplenial area, dorsal part, layer 1','Retrosplenial area, dorsal part, layer 2/3',...
        'Retrosplenial area, dorsal part, layer 5','Retrosplenial area, dorsal part, layer 6a',...
        'Retrosplenial area, dorsal part, layer 6b'}; 
lb={'RSPd'};

% list_to_find={'Ectorhinal area/Layer 1','Ectorhinal area/Layer 2/3',...
%         'Ectorhinal area/Layer 5','Ectorhinal area/Layer 6a',...
%         'Ectorhinal area/Layer 6b'}; 
% lb={'ECT'};


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