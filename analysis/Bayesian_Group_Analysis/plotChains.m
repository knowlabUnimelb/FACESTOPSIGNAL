function plotChains(sampMat, varNames)

sizes = size(sampMat);
if numel(sizes) == 4
    sampMat = permute(sampMat, [1 2 4 3]);
end

nchains = sizes(1);
nsamples = sizes(2);
nvars = prod(sizes(3:end));


[nrs, ncs] = nsubplots(nvars);
figure('WindowStyle', 'docked')
cols = {'r', 'b', 'g', 'k'};
for i = 1:nvars
    subplot(nrs, ncs, i)
    h = plot(sampMat(:,:,i)');
    for k = 1:nchains
        set(h(k), 'Color', cols{k})
    end
    title(varNames{i})
end



