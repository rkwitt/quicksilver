function [iString, iVersionStr] = getInterpreter()

iString = 'octave';
iVersionStr = version;

try
  isempty(OCTAVE_VERSION);
catch
  iString = 'matlab';
end

