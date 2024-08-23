function [all_3 all_com all_together all_b all_m allnr_b allnr_m  all_db all_dm all_bolus bolus_center] = serial2p_vispkim(data,type,vcortex_bino,vcortex_mono,vcortex_all)
%function read outs cell fraction per layer for the bino and mono viusal area based on the Kim Atlas
%% Inputs and outputs 
%Inputs
%data = data matlab structure 
%type = experiments selection 
%vcortex_bino = names of bino cortical layers based on Kim Atlas nomenclature
%vcortex_mono = names of mono cortical layers based on Kim Atlas nomenclature
%vcortex_all = names of anterior and posterioir cortical layers which ar not in mono and bino area based on Kim Atlas nomenclature

%Outputs 
%all_3= bino, mono and anterior posterior combined
%all_com=anterior posterior combined
%all_b = fraction per layer bino
%all_m = fraction per layer mono
%allnr_b = absolute nr per layer bino
%allnr_m = absolute nr per layer mono
%all_bolus = bolus injection volume 

midl= 5699;
%% Read out cell nr etc 
all_b=[];all_m=[];all_db=[];all_dm=[];
all_bolus=[];allnr_m=[];allnr_b=[];
bolus_center=[];all_together=[];all_com=[];all_3=[];


for i=1:length(type)
contra_b=[];contra_m=[];contra_dm=[];contra_db=[];contra_all=[];
  for k=1:length(vcortex_bino);
  cortex_idxb=[];cortex_idxb=find(matches(string(data(type(i)).area_names),vcortex_bino{k})==1);
  cortex_idxm=[];cortex_idxm=find(matches(string(data(type(i)).area_names),vcortex_mono{k})==1);
  cortex_idxa=[];cortex_idxa=find(matches(string(data(type(i)).area_names),vcortex_all{k})==1);

  contra_b(k)=data(type(i)).contra_cellnr(cortex_idxb);
  contra_m(k)=data(type(i)).contra_cellnr(cortex_idxm);
  contra_all(k)=data(type(i)).contra_cellnr(cortex_idxa);

  contra_db(k)=data(type(i)).contra_celldens(cortex_idxb);
  contra_dm(k)=data(type(i)).contra_celldens(cortex_idxm);

  end
  both_com=[];both_3=[];
  %combine just the mono and contra portion 
  both_com=contra_b+contra_m;
  %also add the most posterior and most anterior portion of VISp which does
  %not have bino mono based on KIM atlas 
  both_3=both_com+contra_all;
  cb_l=[];cm_l=[];c_both=[];c_all=[];c_3=[];
  for l=2:6
        cb_l(l)=nansum(contra_b(l))/(nansum(nansum(contra_b))-nansum(nansum(contra_b(1))));
        cm_l(l)=nansum(contra_m(l))/(nansum(nansum(contra_m))-nansum(nansum(contra_m(1))));
        c_both(l)=both_com(l)/(nansum(both_com)-(nansum(both_com(1))));
        c_all(l)=nansum(contra_all(l))/(nansum(nansum(contra_all))-nansum(nansum(contra_all(1))));
        c_3(l)=both_3(l)/(nansum(both_3)-(nansum(both_3(1))));
  end

  all_b(:,i)=cb_l;all_m(:,i)=cm_l;all_together(:,i)=c_both;all_com(:,i)=c_all;all_3(:,i)=c_3;
  all_db(:,i)=contra_db;all_dm(:,i)=contra_dm;
  all_bolus(i)=sum([data(type(i)).bolus{:,6}]);
  bolus_center(:,i)=[data(type(i)).bolus_sum{9} data(type(i)).bolus_sum{10} data(type(i)).bolus_sum{11}];
  allnr_m(i)=sum(contra_m);allnr_b(i)=sum(contra_b);
end

for i=1:length(type)
if bolus_center(3,i)<midl
bolus_centerx(i)=bolus_center(3,i);
else
bolus_centerx(i)=midl-(bolus_center(3,i)-midl);
end
end

% midline correction x to one hemisphere

bolus_center(3,:)=bolus_centerx;
end
