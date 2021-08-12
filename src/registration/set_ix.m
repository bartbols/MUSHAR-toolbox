function set_ix(file,prop,value)

% set_ix sets a value to an existing property in a parameter/transformation file
% 
% Usage:
% set_ix(file,property_string,value)
% 
% -----------------------------Inputs--------------------------------------
% file: pathname of the text(ASCII) parameter/transformation file
% property_string: name of the new property
% value: value of the new property
%        value can be a string (for files or, e.g. the type of cost function 
%        in a parameter file, say 'AdvancedMattesMutualInformation') or a 
%        vector (e.g. the value of property 'Size' in a transformation file)
%
% 5/3/2011
% version 1.0b
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
flg = false;
for i = 1:length(str)
    ind = findstr(str{i},['(' prop ' ']);
    if ~isempty(ind)
        if ischar(value)
            if value(1)=='"'
                value(1) = [];
            end
            if value(end)=='"'
                value(end) = [];
            end
            str{i} = ['(' prop ' "' value '")']; %#ok<AGROW>
        elseif isnumeric(value)
            str{i} = ['(' prop ' ' num2str(value) ')']; %#ok<AGROW>
        end
        flg = true;
        break
    end
end
if flg
    fid = fopen(file,'wt');
    for i = 1:length(str)
        fprintf(fid,'%s\n',str{i});
    end
    fclose(fid);
else
    warning('Property ''%s'' is not present in file %s\nUse add_ix.m instead',prop,file); %#ok<WNTAG>
end