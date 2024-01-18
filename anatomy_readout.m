%serial 2p data TRACING readoutcsv files and creat structure with all
savefile=1;
%% Set directories and experimentator
%summary csv folder
data_dir    = 'C:\Users\simonw\S1-V1 interaction Dropbox\Simon_Manuel\Manuscripts\L6_dominance\Injections\summaries_all';

save_dir        = 'C:\Users\simonw\S1-V1 interaction Dropbox\Simon_Manuel\Manuscripts\L6_dominance\Injections\output_structure';
ExpXls            = 'C:\Users\simonw\S1-V1 interaction Dropbox\Simon_Manuel\Manuscripts\L6_dominance\Injections\excel_sheet_all_info\serial2p_readout.xlsx'

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
  data(adder).gender=batchopt.gender{i};
  data(adder).age=batchopt.age{i};
  %data(adder).area_names=string(dataArray{1, 1});
  data(adder).area_names=dataArray{1, 1};
  data(adder).ipsi_cellnr=data_num(:,1);  
  data(adder).contra_cellnr=data_num(:,2);  
  data(adder).ipsi_celldens=data_num(:,7);  
  data(adder).contra_celldens=data_num(:,8);  
    adder=adder+1;  
   

end


    if savefile==1
        cd(save_dir);
        FileName_save=['Data_','SW','_',datestr(now, 'hh-dd-mmm-yyyy')];
        %         save(FileName,'-struct','LGN','-v7.3');
        save(FileName_save,'data','-v7.3');
        writetable(struct2table(data), [FileName_save '.csv']);
        disp('FILE SAVED');
        
    else
        disp('FILE NOT SAVED');
    end

