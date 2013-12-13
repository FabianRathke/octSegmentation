function options = setEvalDefaults(options)

if ~isfield(options,'EdgesEval')
    options.EdgesEval = 1:9;
end

if ~isfield(options,'LayersEval')
    options.LayersEval = 1:10;
end

if ~isfield(options,'numRows')
    options.numRows = 496;
end

if ~isfield(options,'interpolation')
    options.interpolation = 1;
end

options.numLayers = 10;

end
