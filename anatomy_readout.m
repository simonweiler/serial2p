%serial 2p data TRACING readoutcsv files and creat structure with all
clear all;
savefile=1;
pointsload=1;
%% Set directories and experimentator
%summary csv folder
%data_dir    = 'C:\Users\simonw\S1-V1 interaction Dropbox\Simon_Manuel\Manuscripts\L6_dominance\Injections\summaries_all';
data_dir    = 'C:\Users\simonw\S1-V1 interaction Dropbox\Simon Weiler\VISpL6_callosal\Anatomy\summaries_all';

%save_dir        = 'C:\Users\simonw\S1-V1 interaction Dropbox\Simon_Manuel\Manuscripts\L6_dominance\Injections\output_structure';
save_dir    =  'C:\Users\simonw\S1-V1 interaction Dropbox\Simon Weiler\VISpL6_callosal\Anatomy\output_structure';

%ExpXls            = 'C:\Users\simonw\S1-V1 interaction Dropbox\Simon_Manuel\Manuscripts\L6_dominance\Injections\excel_sheet_all_info\serial2p_readout.xlsx'
ExpXls            = 'C:\Users\simonw\S1-V1 interaction Dropbox\Simon Weiler\VISpL6_callosal\Anatomy\serial2p_readout_visp.xlsx';

batchopt          = parseExperimentsXls_serial2p(ExpXls);%calls the nested function parseExperimentsXls_ephys and considers the user flag (1 or 0)
nummice           = length(batchopt.mouse);%length of experiments to be analyzed
%% SETING CSV PARAMATERS
startRow = 2;
endRow = inf;
delimiter = ',';
%% 

adder=1;%counting variable
for i=1:nummice%for loop over experiments across days
    filename=fullfile(data_dir,char(batchopt.mouseID{i}),'summary.csv');
    fileID=[];fileID = fopen(filename,'r');    
    formatSpec = '%C%f%f%f%f%f%f%f%f%[^\n\r]';    
    dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    %order always ipsi - contra  (so cell numbers columns always start with ipsi and
    %then contra)
    data_num=[]; 
        if batchopt.hemi_ids{i}==0   
            data_num=[dataArray{1, 2} dataArray{1, 3} dataArray{1, 4} dataArray{1, 5} dataArray{1, 6} dataArray{1, 7} dataArray{1, 8} dataArray{1, 9}];
        else batchopt.hemi_ids{i}==1   
            data_num=[dataArray{1, 3} dataArray{1, 2} dataArray{1, 4} dataArray{1, 6} dataArray{1, 5} dataArray{1, 7} dataArray{1, 9} dataArray{1, 8}];
        end
    fclose(fileID);

  data(adder).animal_name=[char(batchopt.mouseID{i})];
  data(adder).recording_date=batchopt.mouse{i};
  data(adder).virustype=batchopt.type{i};
  data(adder).area=batchopt.injectarea{i};
  data(adder).hemisphere=batchopt.inject_hem{i};
  data(adder).selhemi=batchopt.hemi_ids{i};
  data(adder).gender=batchopt.gender{i};
  data(adder).age=batchopt.age{i};
  %data(adder).area_names=string(dataArray{1, 1});
  data(adder).area_names=dataArray{1, 1};
  data(adder).ipsi_cellnr=data_num(:,1);  
  data(adder).contra_cellnr=data_num(:,2);  
  data(adder).ipsi_celldens=data_num(:,7);  
  data(adder).contra_celldens=data_num(:,8);  

  %add all_points csv as well for KIM atlas VISp
 if pointsload==1
     %collect_b={};collect_m={};
     vcortex_bino={'Primary visual cortex, binocular area, layer1';'Primary visual cortex, binocular area, layer2/3';
    'Primary visual cortex, binocular area, layer4';'Primary visual cortex, binocular area, layer5';
    'Primary visual cortex, binocular area, layer6a';'Primary visual cortex, binocular area, layer6b'};
vcortex_mono={'Primary visual cortex, monocular area, layer1';'Primary visual cortex, monocular area, layer2/3';
    'Primary visual cortex, monocular area, layer4';'Primary visual cortex, monocular area, layer5';
    'Primary visual cortex, monocular area, layer6a';'Primary visual cortex, monocular area, layer16b'};
vcortex_all={'Primary visual cortex, layer1';'Primary visual cortex, layer2/3';
    'Primary visual cortex, layer4';'Primary visual cortex, layer5';
    'Primary visual cortex, layer6a';'Primary visual cortex, layer6b'};

      filename_allp=fullfile(data_dir,char(batchopt.mouseID{i}),'all_points.csv');
      opts = detectImportOptions(filename_allp);
      opts.SelectedVariableNames = {'structure_name','hemisphere','coordinate_atlas_axis_0','coordinate_atlas_axis_1','coordinate_atlas_axis_2'};
      allp_info = readtable(filename_allp,opts);
       allpoints=table2cell(allp_info);
       point_names=[];hemi_names=[];
       point_names=categorical(allpoints(:,1));
       hemi_names=categorical(allpoints(:,2));

       if batchopt.hemi_ids{i}==0
       for k=1:6
       cortex_idxb=[];cortex_idxb=find(matches(string(point_names),vcortex_bino{k})==1 & matches(string(hemi_names),'right')==1);
       cortex_idxm=[];cortex_idxm=find(matches(string(point_names),vcortex_mono{k})==1 & matches(string(hemi_names),'right')==1);
       cortex_idxall=[];cortex_idxall=find(matches(string(point_names),vcortex_all{k})==1 & matches(string(hemi_names),'right')==1);
       collect_b{:,k}=allpoints(cortex_idxb,:);
       collect_m{:,k}=allpoints(cortex_idxm,:);
       collect_all{:,k}=allpoints(cortex_idxall,:);
       end
       else batchopt.hemi_ids{i}==1
           for k=1:6
       cortex_idxb=[];cortex_idxb=find(matches(string(point_names),vcortex_bino{k})==1 & matches(string(hemi_names),'left')==1);
       cortex_idxm=[];cortex_idxm=find(matches(string(point_names),vcortex_mono{k})==1 & matches(string(hemi_names),'left')==1);
       cortex_idxall=[];cortex_idxall=find(matches(string(point_names),vcortex_all{k})==1 & matches(string(hemi_names),'left')==1);
       collect_b{:,k}=allpoints(cortex_idxb,:);
       collect_m{:,k}=allpoints(cortex_idxm,:);
       collect_all{:,k}=allpoints(cortex_idxall,:);
       end
       end
    
     data(adder).points_mono=collect_m;  
     data(adder).points_bino=collect_b;  
     data(adder).points_all=collect_all;  
 end

  %if bolusregion is present
 try
      filename_inj=fullfile(data_dir,char(batchopt.mouseID{i}),'bolus.csv');
      filename_injsum=fullfile(data_dir,char(batchopt.mouseID{i}),'summary_bolus.csv');     
     injection_info = readtable(filename_inj);
     injection_infosum = readtable(filename_injsum);
    %read out nr, abbreviation and long cortex name of ABI
    inj_perc=table2cell(injection_info);
    inj_percsum=table2cell(injection_infosum);
%add bolus info 
 data(adder).bolus=inj_perc;  
 data(adder).bolus_sum=inj_percsum; 
      catch
          disp('no bolus rendering present')
      end
    adder=adder+1;  
   

end


    if savefile==1
        cd(save_dir);
        FileName_save=['Data_','SW','_',datestr(now, 'hh-dd-mmm-yyyy')];
        %         save(FileName,'-struct','LGN','-v7.3');
        save(FileName_save,'data','-v7.3');
        %writetable(struct2table(data), [FileName_save '.csv']);
        disp('FILE SAVED');
        
    else
        disp('FILE NOT SAVED');
    end

