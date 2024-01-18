function [indexcalc] = anatomy_indexcalc(layers)
for i=1:size(layers,2)
%iln=(L5+L6a/b) / (L2/3+L5+L6a/b);
iln(i)=(layers(4,i)+layers(5,i)+layers(6,i)) / (layers(2,i)+layers(4,i)+layers(5,i)+layers(6,i));

%L6i=(L6a/b-L2/3) / (L6a/b+L2/3);
L6i(i)=((layers(5,i)+layers(6,i))-layers(2,i)) / ((layers(5,i)+layers(6,i))+layers(2,i));
%L6ai=(L6a-L2/3) / (L6a+L2/3);
L6ai(i)=(layers(5,i)-layers(2,i)) / (layers(5,i)+layers(2,i));
%h_index=1-((L2/3+L4) / (L2/3+L4+L5));
h_index(i)=1-((layers(2,i)+layers(3,i)) / (layers(2,i)+layers(3,i)+layers(4,i)));
%pure L6a
%(L6a)/(L23 + L5 +L6ab)
pL6a(i)=(layers(5,i)) / (layers(2,i)+layers(4,i)+layers(5,i)+layers(6,i));
end
indexcalc=[iln' L6i' L6ai' h_index' pL6a'];
for k=1:size(indexcalc,2)
    if sum(isnan(indexcalc(:,k)))>=5
    indexcalc(find(~isnan(indexcalc(:,k))),k)=NaN;
    else
    end
end
end