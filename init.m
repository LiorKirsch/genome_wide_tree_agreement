
if exist('inited', 'var')
  return
end
parms.dummy=0;

% General parameters
[seed, parms] = take_from_struct(parms, 'seed', 1);

addpath('code-visualization');
addpath('code-general use functions');
inited = true;
