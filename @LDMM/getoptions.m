function v = getoptions(options, name, v_default)

if nargin<3
    error('Not enough arguments.');
end

if ~isfield(options, name) || isempty(eval(['options.' name ';']))
    warning(['option does not contain field "' name '", set to default value ' name '=' num2str(v_default)]);
    v = v_default;
else
    v = eval(['options.' name ';']);
end

end