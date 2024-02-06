function [i_animal c_animal i_areas_animal c_areas_animal ...
    i_areas_animalf c_areas_animalf i_areas_animalcd c_areas_animalcd] = anatomy_cellnr(data,type,cortex_names)
aa=[18 25 31 41];
for i=1:length(type)
  ipsi_all=[];contra_all=[];ipsi_allcd=[];contra_allcd=[];
  %k=areas
  for k=1:length(cortex_names)
  cortex_idx=[];cortex_idx=find(contains(string(data(type(i)).area_names),cortex_names{k,3})==1);
  a_names=[];a_names=data(type(i)).area_names(cortex_idx);
  %sort from L1 to L6b in a given area
  rb=[];sort_v=[];[rb sort_v]=sort(a_names);
  %save cell numbers with the sorted indx for ipsi and contra
  ipsi_nr=[];ipsi_nr=data(type(i)).ipsi_cellnr(cortex_idx(sort_v));
  contra_nr=[];contra_nr=data(type(i)).contra_cellnr(cortex_idx(sort_v));
  ipsi_cd=[];ipsi_cd=data(type(i)).ipsi_celldens(cortex_idx(sort_v));
  contra_cd=[];contra_cd=data(type(i)).contra_celldens(cortex_idx(sort_v));
  %distinguish areas that dont have a L4 (5 layers only) from the rest
  %five layers
  if length(sort_v)==5
      ipsi_nr=[ipsi_nr(1:2); NaN ;ipsi_nr(3:end)];
      contra_nr=[contra_nr(1:2); NaN ;contra_nr(3:end)];
      ipsi_cd=[ipsi_cd(1:2); NaN ;ipsi_cd(3:end)];
      contra_cd=[contra_cd(1:2); NaN ;contra_cd(3:end)];
  %six layers
  else
      ipsi_nr=ipsi_nr;contra_nr=contra_nr;ipsi_cd=ipsi_cd;contra_cd=contra_cd;
  end
  %save across areas
  ipsi_all(:,k)=ipsi_nr;
  contra_all(:,k)=contra_nr;
  ipsi_allcd(:,k)=ipsi_cd;
  contra_allcd(:,k)=contra_cd;
  end
  %remove respective injection area on the ipsi site
  idx_inj=[];idx_inj=find(nansum(ipsi_all)==max(nansum(ipsi_all)));
  disp([num2str(cortex_names{idx_inj,3}) ' was injected'])
    if  any(idx_inj==aa)==1
        idx_inj=idx_inj;
    else
        idx_inj=idx_inj-1;
    end
  ipsi_all(:,idx_inj)=[NaN;NaN;NaN;NaN;NaN;NaN];
  ipsi_allcd(:,idx_inj)=[NaN;NaN;NaN;NaN;NaN;NaN];

  ihemi=[];chemi=[];
  %calculate fraction across layers (minus L1)
        for l=2:6
        ihemi(l)=nansum(ipsi_all(l,:))/(nansum(nansum(ipsi_all))-nansum(nansum(ipsi_all(1,:))));
        chemi(l)=nansum(contra_all(l,:))/(nansum(nansum(contra_all))-nansum(nansum(contra_all(1,:))));
        end
  %save fractions across animlas regardless of area
  i_animal(:,i)=ihemi;
  c_animal(:,i)=chemi;
    area_li=[];area_lc=[];
    for m=1:length(cortex_names)
       temp_li=[];temp_lc=[];
       for l=2:6
       temp_li(l)=sum(ipsi_all(l,m))/(nansum(ipsi_all(:,m))-nansum(ipsi_all(1,m)));
       temp_lc(l)=sum(contra_all(l,m))/(nansum(contra_all(:,m))-nansum(contra_all(1,m)));
       end
    area_li(:,m)=temp_li;
    area_lc(:,m)=temp_lc;
   end
  i_areas_animal(:,:,i)=ipsi_all;
  c_areas_animal(:,:,i)=contra_all;
  i_areas_animalf(:,:,i)=area_li;
  c_areas_animalf(:,:,i)=area_lc;
  i_areas_animalcd(:,:,i)=ipsi_allcd;
  c_areas_animalcd(:,:,i)=contra_allcd;

end

end