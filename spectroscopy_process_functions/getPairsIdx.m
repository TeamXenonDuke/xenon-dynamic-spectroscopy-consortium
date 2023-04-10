spectData = readtable('D:\OneDrive\Documents\Duke\CIVM Research\Siemens Data\Processed Data\spectroscopy_stats_master.xlsx','Sheet','Pairs');
pairs = find(contains(all_twix_dyn,strcat(spectData.Subject,'\')));
length(pairs)

% Add in 000-001G/H pair
in1 = find(contains(all_twix_dyn,"000-001G"))
a1 = find(contains(all_twix_dyn,"000-001H"))
pairs = [in1; pairs];
pairs(2) = a1;

% Add in 002-098/A pair
in1 = find(contains(spectData.Subject,"002-098"))
a1 = find(contains(all_twix_dyn,"002-098A"))
pairs = [pairs(1:2*in1-1);a1;pairs(2*in1:end)]

%remove 005-045s2
rm1 = find(contains(all_twix_dyn,"005-045"))
pairs(pairs == rm1(2)) = [];
all_twix_dyn(pairs)
