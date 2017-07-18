function matlabDataType = decideMETADataType(metaDataType)
% convience function to map META data type to matlab data type
matlabDataType = '';
if strcmpi('MET_FLOAT',metaDataType)
  matlabDataType = 'float32';
elseif strcmpi('MET_DOUBLE',metaDataType)
  matlabDataType = 'float64';
elseif strcmpi('MET_USHORT',metaDataType)
  matlabDataType = 'uint16';
elseif strcmpi('MET_SHORT',metaDataType)
  matlabDataType = 'int16';
elseif strcmpi('MET_UCHAR',metaDataType)
  matlabDataType = 'uint8';
end