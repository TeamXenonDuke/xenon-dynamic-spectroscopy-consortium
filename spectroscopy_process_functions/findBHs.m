function [BHstart, BHend] = findBHs(t, BHs)

BHstart_all = find(round(t,1) == BHs(1));
BHend_all = find(round(t,1) == BHs(2));

BHstart = BHstart_all(ceil(end/2));
   
if isempty(BHend_all)
    BHend = length(t)-50;
else
    BHend = BHend_all(ceil(end/2));
end 
