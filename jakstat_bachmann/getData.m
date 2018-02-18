function D = getData()

load('data/data_Bachmann.mat','D')

D = getOffsetScalingStd_Bachmann(D);
D = loadInitialConditions(D);

for cond = 1:numel(D)
    D(cond).my(:,[1:10,12:end],:) = 10.^D(cond).my(:,[1:10,12:end],:);
end
D(3).my = D(3).my - 1;

end