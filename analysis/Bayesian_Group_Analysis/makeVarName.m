function varNames = makeVarName(samples, namestr)

nrows = size(samples, 3);
ncols = size(samples, 4);
vnVal = allcomb(1:nrows, 1:ncols);
for i = 1:size(vnVal,1); 
    varNames{vnVal(i,1), vnVal(i,2)} = sprintf('%s_{%d,%d}', namestr, vnVal(i,1), vnVal(i,2));
end