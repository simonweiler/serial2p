%Check numbers of cortex areas and read out names and displays them
%% Dependendencies 
%uipickfiles: https://uk.mathworks.com/matlabcentral/fileexchange/10867-uipickfiles-uigetfile-on-steroids
%Matlab structure 
%% Folder where structures are (CHANGE ACCORDINGLY)
str   = 'C:\Users\simonw\S1-V1 interaction Dropbox\Simon_Manuel\Manuscripts\L6_dominance\Injections\output_structure';
% load structure (mat file) using uipickfiles
folder_list = uipickfiles('FilterSpec',str);load(char(folder_list));
% select all retro injection regardeless of their injection type
temp1=[];temp1=cell_selecter(data,'virustype',1);
%% Read out names of all cortical areas (using layer 5 as this produces the most cortical areas)
%use retro animal indices only (number of areas should be the same for antero!)
retro=find(temp1==1);
%search which strings have layer 5 (sometimes "L" or "l" so hence search for
%"ayer"
search_str={'ayer 5'};
%loop through animals 
for i=1:length(retro)
cortex_idx=[];cortex_idx=find(contains(string(data(retro(i)).area_names),search_str{1})==1);
cortex_names=[];cortex_names=data(retro(i)).area_names(cortex_idx);
nr_cortex_areas_animals(i)=length(cortex_names);
    for k=1:length(cortex_names)
        temp=[];temp=string(data(retro(i)).area_names(cortex_idx(k)));
       
        %read out cortical areas using comma seperation 
        try
        coma_str=strfind(temp,',');
        name_all{k,:}=temp{1}(1:coma_str(end)-1);
        %as ectorhinal area is corrupted with / and capital L
        catch
        name_all{k,:}='Ectorhinal area';             
        end
    end    
end
disp(['Animals have the following cortex area numbers: ' num2str(nr_cortex_areas_animals)])

names_ordered=sort(name_all)