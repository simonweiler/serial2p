function [batchopt] = parseExperimentsXls_serial2p(path)

[xls_num,xls_txt]=xlsread(path);
%read out colums of interest
animalname    = find(~cellfun(@isempty, strfind(xls_txt(1,:),'AnimalID')));
loadcol        = find(~cellfun(@isempty, strfind(xls_txt(1,:),'BatchAnalyze')));
datecol       = find(~cellfun(@isempty, strfind(xls_txt(1,:),'ExperimentalDay')));
hemicol       = find(~cellfun(@isempty, strfind(xls_txt(1,:),'Contra_hemisphere')));
typecol= find(~cellfun(@isempty, strfind(xls_txt(1,:),'Type')));
injectiontype= find(~cellfun(@isempty, strfind(xls_txt(1,:),'Injection_area')));
injectionhem= find(~cellfun(@isempty, strfind(xls_txt(1,:),'Injection_hemi')));
gender_col= find(~cellfun(@isempty, strfind(xls_txt(1,:),'Gender')));
age_col= find(~cellfun(@isempty, strfind(xls_txt(1,:),'Age')));



k = 1;

batchopt.XLS.txt = xls_txt;
batchopt.XLS.num = xls_txt;

for i = 2:size(xls_txt,1)
    ana{k}= xls_num(i-1,loadcol-1);
    
    if ~ana{k}
        disp(['skipping experiments ' xls_txt{i,animalname} ' (no batchload flag)']);
        continue
    end
    
    batchopt.mouse{k}     = xls_num(i-1,datecol-1);
    batchopt.mouseID{k}   = xls_txt(i,animalname);

    
    
    expcellids{k}                = xls_txt(i,hemicol);
    expcellids2{k}                = xls_txt(i,typecol);
    expcellids3{k}              = xls_txt(i,injectiontype);
    expcellids4{k}              = xls_txt(i,injectionhem);
    expcellids5{k}              = xls_txt(i,gender_col);
    expcellids6{k}              = xls_txt(i,age_col);
    
    batchopt.hemi_ids{k}        = str2num((expcellids{k}{1}));
    batchopt.type{k}           = str2num((expcellids2{k}{1}));
    batchopt.injectarea{k}         = str2num((expcellids3{k}{1}));
    batchopt.inject_hem{k}         = str2num((expcellids4{k}{1}));
    batchopt.gender{k}         = str2num((expcellids5{k}{1}));
    batchopt.age{k}         = str2num((expcellids6{k}{1}));
 
    %batchopt.loaddrive{k}        = xls_txt(i,loaddrivecol);
    k = k+1;
end

