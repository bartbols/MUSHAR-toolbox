function value = get_ix(file,prop)

% get_ix retrieves the value of a property from a parameter/transformation file
% 
% Usage:
% value = get_ix(file,property_string)
%
% -----------------------------Inputs--------------------------------------
% file: pathname of the text(ASCII) transformation file
% property_string: name of the requested property
% ------------------------------Outputs----------------------------------
% value: value of the queried property (value can be a string or a vector)
%        value can be a string (for files or, e.g. the type of cost function 
%        in a parameter file, say 'AdvancedMattesMutualInformation') or a 
%        vector (e.g. the value of property 'Size' in a transformation file)
%        (empty value means the name of the property provided by "prop" is
%        not present in the file)
%
% version 1.0b
% 5/3/2011
% 
% Pedro Antonio Valdés Hernández
% Neuroimaging Department
% Cuban Neuroscience Center

fid = fopen(file);
flg = true;
str = [];
c = 0;
while flg
    c = c+1;
    tmp = fgetl(fid);
    if (~ischar(tmp) && tmp == -1)
        flg = false;
    else
        str{c} = deblank(tmp); %#ok<AGROW>
    end
end
fclose(fid);
for i = 1:length(str)
    ind = findstr(str{i},['(' prop ' ']);
    if ~isempty(ind)
        value = str{i}(length(['(' prop])+2:end-1);
        if value(1) == '"'
            value([1 end]) = [];
        else
            value = str2num(value); %#ok<ST2NM>
        end
        break
    end
end
if ~exist('value','var')
    value = '';
    warning('Property ''%s'' is not present in file %s',prop,file); %#ok<WNTAG>
end