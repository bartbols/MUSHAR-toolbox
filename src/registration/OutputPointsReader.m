
function out =  OutputPointsReader(filename,dim)
% OUTPUTPOINTSREADER  Reads in data from stored in 'filename', which should
% be a file with transformed points that is genereated with Transformix

if nargin == 1
    dim = 3;
else
    if ~ismember(dim,[2 3])
        error('''dim'' should be 2 or 3')
    end
end

fid = fopen(filename);
% formatstr = 'Point %d; %s = [ %d %d %d ] ; %s = [ %f %f %f ] ; %s = [ %d %d %d] ; %s = [ %f %f %f ] ; %s = [ %f %f %f ] ; %s = [ %d %d %d ]';
if dim == 2
    formatstr = 'Point %d; %s = [ %d %d ] ; %s = [ %f %f ] ; %s = [ %d %d] ; %s = [ %f %f ] ; %s = [ %f %f ]';
elseif dim == 3
    formatstr = 'Point %d; %s = [ %d %d %d ] ; %s = [ %f %f %f ] ; %s = [ %d %d %d] ; %s = [ %f %f %f ] ; %s = [ %f %f %f ]';
end
mydata = textscan(fid,formatstr);
fclose(fid);

if dim==2
    % Index of points that were transformed
    % PointNumber = mydata{1};
    out.InputIndex(:,1) = mydata{3};
    out.InputIndex(:,2) = mydata{4};
%     out.InputIndex(:,3) = mydata{5};
    
    % Physical location of points that were transformed
    out.InputPoint(:,1) = mydata{6};
    out.InputPoint(:,2) = mydata{7};
%     out.InputPoint(:,3) = mydata{9};
    
    % Index of output points in fixed image coordinates
    out.OutputIndexFixed(:,1) = mydata{9};
    out.OutputIndexFixed(:,2) = mydata{10};
%     out.OutputIndexFixed(:,3) = mydata{13};
    
    % Physical location of output points
    out.OutputPoint(:,1) = mydata{12};
    out.OutputPoint(:,2) = mydata{13};
%     out.OutputPoint(:,3) = mydata{17};
    
    % Deformation of points
    out.Deformation(:,1) = mydata{15};
    out.Deformation(:,2) = mydata{16};
%     out.Deformation(:,3) = mydata{21};
    
elseif dim == 3
    % Index of points that were transformed
    % PointNumber = mydata{1};
    out.InputIndex(:,1) = mydata{3};
    out.InputIndex(:,2) = mydata{4};
    out.InputIndex(:,3) = mydata{5};
    
    % Physical location of points that were transformed
    out.InputPoint(:,1) = mydata{7};
    out.InputPoint(:,2) = mydata{8};
    out.InputPoint(:,3) = mydata{9};
    
    % Index of output points in fixed image coordinates
    out.OutputIndexFixed(:,1) = mydata{11};
    out.OutputIndexFixed(:,2) = mydata{12};
    out.OutputIndexFixed(:,3) = mydata{13};
    
    % Physical location of output points
    out.OutputPoint(:,1) = mydata{15};
    out.OutputPoint(:,2) = mydata{16};
    out.OutputPoint(:,3) = mydata{17};
    
    % Deformation of points
    out.Deformation(:,1) = mydata{19};
    out.Deformation(:,2) = mydata{20};
    out.Deformation(:,3) = mydata{21};
    
end
end % of function
