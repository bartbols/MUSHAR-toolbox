function reg_elastix( fixed, moving, parfile,varargin )
%REG_ELASTIX Calls elastix to perform a registration on the fixed and moving
% image using parameter file 'parfile'. 'parfile' can also be a cell
% with multiple parameter files, in which case a multi-step registration
% will be performed.
%
% Bart Bolsterlee, Neuroscience Research Australia (NeuRA)
% October 2017
%
% ----------------- INPUT -----------------
% ----- REQUIRED -----
% fixed:    filename of the fixed image
% moving:   filename of the moving image
% parfile:  filename (char) or a cell string of filenames of the elastix
%           parameter file(s) with registration parameters. For example, to
%           perform a rigid and then a bspline registration: parfile =
%           {'rigid.txt','bspline.txt'}. Or to only perform a bspline, use
%           'bspline.txt'.
%
% ----- OPTIONAL -----
% Optional inputs are provided as pairs of 'ParameterName',<value> (e.g.
% 'mask','mask.nii.gz')
%
% - transform_file       : filename of the final transformation (as created
%                          by Elastix)
% - mask                 : filename of the binary mask file used as the fixed
%                          mask for registration.
%                          Registration will be optimised for the region in
%                          the mask only; regions outside the mask may have
%                          poor alignment.
% - mask_m               : filename of the binary mask file used as the moving
%                          mask for registration.
%                          Registration will be optimised for the region in
%                          the mask only; regions outside the mask may have
%                          poor alignment.
% - label_number         : label number(s) in the mask used for segmentation
%                          (only required if mask contains multiple labels)
% - foreground_threshold : threshold intensity for foreground. A mask will 
%                          be created from the fixed image using all voxels
%                          with a value above this threshold.
% - stack_f              : if fixed image is 4D, choose which stack is
%                          used for registration.
% - stack_m              : if moving is 4D, choose which stack is
%                          used for registration.
% - surface_in           : STL filename of surface model to be transformed
%                          with the resulting transformation
% - surface_out          : STL filename of surface model after
%                          transformation (will be created)
% - result_image         : filename of the transformed image (moving to
%                          fixed)
% - deformation          : filename of the deformation field (will be 
%                          created using transformix)
% - jacobian             : filename of the jacobian of the deformation
%                          field (will be created using transformix)
% - dilate_mask          : number of voxels to grow the mask by prior to
%                          registration.
% - transform_inv        : filename of the inverse transformation file
% - initial              : filename of the initial transformation
% - ref_points_f         : filename of text-file containing reference point
%                          locations in the fixed image.
% - ref_points_m         : filename of text-file containing reference point
%                          locations in the moving image. Points have to be 
%                          paired with with ref_points_f.
% - log_file             : filename to which a copy of the log-file will be
%                          saved.
% - copy_results_dir     : path to which the temporary working directory
%                          will be saved after registration.

% Change log:
% 11/09/2018, BB: Added option to input corresponding points to guide
%                 registration (ref_points_f and ref_points_m)
% 11/09/2018, BB: Added log file as output.
% 

% Read the inputs
p = inputParser;
addRequired(p,'fixed',@(x) contains(x,'.nii'))
addRequired(p,'moving',@(x) contains(x,'.nii'))
addRequired(p,'parfile',@(x) ischar(x) || iscell(x))
addParameter(p,'transform_file',[],@(x) ischar(x) || isempty(x))
addParameter(p,'mask',[],@(x) isempty(x) || contains(x,'.nii') || isnumeric(x) || iscell(x))
addParameter(p,'mask_m',[],@(x) isempty(x) || contains(x,'.nii') || isnumeric(x) || iscell(x))
addParameter(p,'foreground_threshold',[],@(x) isscalar(x) || isempty(x))
addParameter(p,'label_number',[],@(x) assert(isnumeric(x)))
addParameter(p,'stack_f',[],@(x) assert(isnumeric(x)))
addParameter(p,'stack_m',[],@(x) assert(isnumeric(x)))
addParameter(p,'surface_in',[],@(x) contains(x,'.stl'))
addParameter(p,'surface_out',[],@(x) contains(x,'.stl'))
addParameter(p,'result_image',[],@(x) contains(x,'.nii'))
addParameter(p,'jacobian',[],@(x) contains(x,'.nii'))
addParameter(p,'deformation',[],@(x) contains(x,'.nii'))
addParameter(p,'dilate_mask',[],@(x) isscalar(x))
addParameter(p,'transform_inv',[],@(x) ischar(x) || isempty(x))
addParameter(p,'initial',[],@(x) ischar(x) || isempty(x))
addParameter(p,'ref_points_f',[],@(x) exist(x,'file')==2 || isempty(x))
addParameter(p,'ref_points_m',[],@(x) exist(x,'file')==2 || isempty(x))
addParameter(p,'log_file',[],@(x) ischar(x) || isempty(x))
addParameter(p,'copy_results_dir',[])

parse(p,fixed, moving, parfile,varargin{:});

mask           = p.Results.mask;
mask_m         = p.Results.mask_m;
label_number   = p.Results.label_number;
foreground_threshold = p.Results.foreground_threshold;
stack_f        = p.Results.stack_f;
stack_m        = p.Results.stack_m;
surface_out    = p.Results.surface_out;
surface_in     = p.Results.surface_in;
result_image   = p.Results.result_image;
transform_file = p.Results.transform_file;
dilate_mask    = p.Results.dilate_mask;
transform_inv  = p.Results.transform_inv;
initial        = p.Results.initial;
deformation    = p.Results.deformation;
jacobian       = p.Results.jacobian;
ref_points_f   = p.Results.ref_points_f;
ref_points_m   = p.Results.ref_points_m;
log_file       = p.Results.log_file;
copy_results_dir = p.Results.copy_results_dir;

% Check if elastix is available as an external command.
[status,cmdout] = system('elastix');
if contains(cmdout,'is not recognized as an internal or external command')
    error('elastix not found as an external commant. Please install elastix and try again.')
end


if xor(isempty(ref_points_f), isempty(ref_points_m))
    error('Both fixed and moving reference points should be provided.')
end

% Create temporary working directory.
char_list = char(['a':'z' '0':'9']) ;
tmpdir = [];
while exist(tmpdir,'dir') == 7 || isempty(tmpdir)
    tmpdir = fullfile(tempdir,char_list(ceil(length(char_list)*rand(1,8))));
end
mkdir(tmpdir)

try

    if ~isempty(stack_f)
        % Get the selected stack(s) from the 4D fixed image
        fixed_3D = cell(1,length(stack_f));
        for i = 1 : length(stack_f)
            fprintf('In reg_elastix: Extracting stack %d from 4D fixed image... ',stack_f(i))
            fixed_3D{i} = fullfile(tmpdir,sprintf('fixed_image%02d.nii.gz',i));
            extract_3Dfrom4D(fixed,fixed_3D{i},stack_f(i));
            fprintf('completed.\n')
        end
    else
        fixed_3D{1} = fixed;
    end
    if ~isempty(stack_m)
        % Get the selected stack(s) from the 4D moving image
        moving_3D = cell(1,length(stack_f));
        for i = 1:length(stack_m)
            fprintf('In reg_elastix: Extracting stack %d from 4D moving image... ',stack_m(i))
            moving_3D{i} = fullfile(tmpdir,sprintf('moving_image%02d.nii.gz',i));
            extract_3Dfrom4D(moving,moving_3D{i},stack_m(i));
            fprintf('completed.\n')
        end
    else
        moving_3D{1} = moving;
    end
    
    
    if isempty(mask) && ~isempty(foreground_threshold)
        % Create mask based on foreground_threshold. Use the first fixed
        % image to make the mask.
        fprintf('In reg_elastix: Creating foreground mask with threshold %.1f...',foreground_threshold)
        mask_img = load_untouch_nii(fixed_3D{1});
        mask_img.img = cast(mask_img.img * mask_img.hdr.dime.scl_slope>foreground_threshold,'like',mask_img.img);
        mask_img.hdr.dime.scl_slope = 1;
        mask = fullfile(tmpdir,'mask.nii.gz');
        save_untouch_nii(mask_img,mask);
        fprintf('completed.\n')        
    end
    
    if ~isempty(mask)
        % Load the mask
        M = load_untouch_nii(mask);
        if isempty(label_number)
            % Use all non-zero voxels
            M.img = cast(M.img ~= 0,'like',M.img);
        else
            % Extract the selected label from the mask
            M.img = cast(ismember(M.img, label_number),'like',M.img);
        end
        
        if ~isempty(dilate_mask)        
            % Dilate (grow) the mask by a number of voxels
            M.img = cast(imdilate(M.img,ones(dilate_mask*2+1,dilate_mask*2+1)),'like',M.img);
        end
        
        % Save the mask in the temporary working folder.
        mask = fullfile(tmpdir,'mask.nii.gz');
        save_untouch_nii(M,mask)
    end
        
    % Check how many steps there are in the registration
    if ischar(parfile)
        % Only one parameter file is provided.
        parfile = {parfile};
    end
    
    nSteps = length(parfile);
    % Now, iterate through the steps, using the transformation from the
    % previous step as the initial transformation.
    for stepnr = 1 : nSteps
        mkdir(fullfile(tmpdir,sprintf('step%02d',stepnr)))
        
        % Build up the elastix command
        elastix_cmd = 'elastix ';
        if isempty(stack_f) || numel(stack_f) == 1
            % One fixed image is provided.
            elastix_cmd = [elastix_cmd ' -f ' fixed_3D{1}];
        else
            % Multiple fixed images are provided.
            for i = 1: length(stack_f)
                elastix_cmd = [elastix_cmd ' -f' int2str(i-1) ' ' fixed_3D{i}];
            end
        end
        if isempty(stack_m) || numel(stack_m) == 1
            % One fixed image is provided.
            elastix_cmd = [elastix_cmd ' -m ' moving_3D{1}];
        else
            % Multiple fixed images are provided.
            for i = 1: length(stack_m)
                elastix_cmd = [elastix_cmd ' -m' int2str(i-1) ' ' moving_3D{i}];
            end
        end
        
        % Add output directory and parameter file
        elastix_cmd = [elastix_cmd ' -out ' fullfile(tmpdir,sprintf('step%02d',stepnr)),...
                                   ' -p ' parfile{stepnr}];
        
        if ~isempty(mask)
            % Add fixed mask, if provided.
            elastix_cmd = [elastix_cmd sprintf(' -fMask %s',mask)];
            
        end
        if ~isempty(mask_m)
            % Add moving mask, if provided.
            elastix_cmd = [elastix_cmd sprintf(' -mMask %s',mask_m)];
            
        end
        if ~isempty(ref_points_f)
            % Add fixed and moving reference points, if provided.
             elastix_cmd = [elastix_cmd sprintf(' -fp %s -mp %s',ref_points_f,ref_points_m)];
        end
        if stepnr ~= 1
            % After first step, use previous parameter file as initial
            % transform
            elastix_cmd = [elastix_cmd sprintf(' -t0 %s',...
                fullfile(tmpdir,sprintf('step%02d',stepnr-1),'TransformParameters.0.txt'))];
        else
            if ~isempty(initial)
                % An initial transformation is provided. Use this for the
                % first step.
                elastix_cmd = [elastix_cmd sprintf(' -t0 %s',...
                    initial)];
                
            end
        end
        % Run registration with Elastix
        system(elastix_cmd)
        pause(1)
    end
    transform_file_final = fullfile(tmpdir,sprintf('step%02d',nSteps),'TransformParameters.0.txt');
    
    if ~isempty(transform_file)
        % Copy the final transform file.
        
        if exist(fileparts(transform_file),'dir') ~= 7
            mkdir(fileparts(transform_file))
        end
        copyfile(transform_file_final,transform_file);
        
        % If registration was done in multiple steps, the  transformation 
        % of the last step file will refer to the transformation file from
        % previous steps. Also copy those files and update the initial
        % transforms in the final transformation file.
        
        if nSteps > 1
            [a,b,c] = fileparts(transform_file);
            set_ix(transform_file,'InitialTransformParametersFileName',...
                        fullfile(a,sprintf('%s_step%02d%s',b,nSteps-1,c)));
            for stepnr = 1 : nSteps-1
                fname = fullfile(tmpdir,sprintf('step%02d',stepnr),'TransformParameters.0.txt');
                new_fname = fullfile(a,sprintf('%s_step%02d%s',b,stepnr,c));
                copyfile(fname,new_fname);
                if stepnr == 1
                    set_ix(new_fname,'InitialTransformParametersFileName','NoInitialTransform');
                else
                    set_ix(new_fname,'InitialTransformParametersFileName',...
                        fullfile(a,sprintf('%s_step%02d%s',b,stepnr-1,c)));
                end
                     
            end
        end
    end
    
    if ~isempty(transform_inv)
        % Also calculate the inverse transformation
        % This procedure is described in the Elastix v4.8 manual section 6.1.6
        
        par_file_inv = fullfile(tmpdir,'parfile_for_inv.txt');
                
        copyfile(parfile{nSteps},par_file_inv)
        
        set_ix(par_file_inv,'Metric','DisplacementMagnitudePenalty')
        set_ix(par_file_inv,'HowToCombineTransforms','Compose')
        mkdir(fullfile(tmpdir,'inv'))
        elastix_cmd = sprintf('elastix -f %s -m %s -t0 %s -p %s -out %s',...
            fixed_3D{1},fixed_3D{1},...
            transform_file_final,...
            par_file_inv,fullfile(tmpdir,'inv'));
        
        if ~isempty(mask)
            % Add fixed mask, if provided.
            elastix_cmd = [elastix_cmd sprintf(' -fMask %s',mask)];
        end        
        system(elastix_cmd);
        set_ix(fullfile(tmpdir,'inv','TransformParameters.0.txt'),...
            'InitialTransformParametersFileName','NoInitialTransform')
        
        copyfile(fullfile(tmpdir,'inv','TransformParameters.0.txt'),transform_inv)

    end
    
    if ~isempty(surface_in)
        % Read the surface model
        FV = stlread2(surface_in);        
        
        % Transform the vertices with transformix
        FV.vertices = transformix_points(FV.vertices, transform_file_final);
        
        % Write as stl-file.
        stlwrite2(surface_out,FV)
        fprintf('In reg_elastix: Transformed surface saved as %s\n',surface_out)

    end
    
    if ~isempty(result_image)
        % Write the result image
        img_reg = ls(fullfile(tmpdir,sprintf('step%02d',nSteps),'result.*'));
        
        if exist(fileparts(result_image),'dir') ~= 7
            mkdir(fileparts(result_image))
        end
        if isempty(img_reg)
            % Result image was not written. Create with elastix
            transformix_cmd = sprintf('transformix -in %s -out %s -tp %s',...
                moving_3D{1},fullfile(tmpdir,sprintf('step%02d',nSteps)),transform_file);
            system(transformix_cmd)
            img_reg = ls(fullfile(tmpdir,sprintf('step%02d',nSteps),'result.*'));
        end
        % Check extension
        img_reg = fullfile(tmpdir,sprintf('step%02d',nSteps),img_reg);
        if endsWith(result_image,'.nii.gz') && endsWith(img_reg,'.nii')
            % Result image is currently not zipped, while a zipped file is
            % requested as output.
            gzip(img_reg)
            movefile([img_reg '.gz'],result_image)
        elseif endsWith(result_image,'.nii') && endsWith(img_reg,'.nii.gz')
            % Result image is currently zipped, while a unzipped file
            % is requested.
            gunzip(img_reg)
            movefile(img_reg(1:end-3),result_image)
        elseif endsWith(result_image,'.nii.gz') && endsWith(img_reg,'.nii.gz')
            movefile(img_reg,result_image)
        elseif endsWith(result_image,'.nii') && endsWith(img_reg,'.nii')
            movefile(img_reg,result_image)
        end            
    end

    if ~isempty(deformation)
        fprintf('Calculating deformation field...')
        
        % Create deformation field with transformix.
        transformix_cmd = sprintf('transformix -def all -out %s -tp %s',...
            tmpdir,transform_file);
        system(transformix_cmd);
        
        if exist(fileparts(deformation),'dir') ~= 7
            mkdir(fileparts(deformation))
        end
        % Check extension
        tmpname = ls(fullfile(tmpdir,'deformationField.*'));
        if endsWith(deformation,'.nii.gz') && endsWith(tmpname,'.nii')
            % Result image is currently not zipped, while a zipped file is
            % requested as output.
            gzip(fullfile(tmpdir,tmpname))
            movefile(fullfile(tmpdir,'deformationField.nii.gz'),deformation)
        elseif endsWith(deformation,'.nii') && endsWith(tmpname,'.nii.gz')
            % Result image is currently zipped, while a unzipped file
            % is requested.
            gunzip(fullfile(tmpdir,tmpname))
            movefile(fullfile(tmpdir,'deformationField.nii'),deformation)
        elseif endsWith(deformation,'.nii.gz') && endsWith(tmpname,'.nii.gz')
            movefile(fullfile(tmpdir,'deformationField.nii.gz'),deformation)
        elseif endsWith(deformation,'.nii') && endsWith(tmpname,'.nii')
            movefile(fullfile(tmpdir,'deformationField.nii'),deformation)
        end
        fprintf(' completed. Saved as: %s\n',deformation)
    end
    
    
    if ~isempty(jacobian)
        fprintf('Calculating spatial jacobian...')
        
        % Create spatial jacobian with transformix.
        transformix_cmd = sprintf('transformix -jacmat all -out %s -tp %s',...
            tmpdir,transform_file);
        system(transformix_cmd);
        if exist(fileparts(jacobian),'dir') ~= 7
            mkdir(fileparts(jacobian))
        end
        tmpname = ls(fullfile(tmpdir,'fullSpatialJacobian.*'));
        if endsWith(jacobian,'.nii.gz') && endsWith(tmpname,'.nii')
            % Result image is currently not zipped, while a zipped file is
            % requested as output.
            gzip(fullfile(tmpdir,tmpname))
            movefile(fullfile(tmpdir,'fullSpatialJacobian.nii.gz'),jacobian)
        elseif endsWith(jacobian,'.nii') && endsWith(tmpname,'.nii.gz')
            % Result image is currently zipped, while a unzipped file
            % is requested.
            gunzip(fullfile(tmpdir,tmpname))
            movefile(fullfile(tmpdir,'fullSpatialJacobian.nii'),jacobian)
        elseif endsWith(jacobian,'.nii.gz') && endsWith(tmpname,'.nii.gz')
            movefile(fullfile(tmpdir,'fullSpatialJacobian.nii.gz'),jacobian)
        elseif endsWith(jacobian,'.nii') && endsWith(tmpname,'.nii')
            movefile(fullfile(tmpdir,'fullSpatialJacobian.nii'),jacobian)
        end        
%         movefile(fullfile(tmpdir,'fullSpatialJacobian.nii.gz'),jacobian)
        fprintf(' completed. Saved as: %s\n',jacobian)
    end
    
    % Copy the log file
    if ~isempty(log_file)
        copyfile(fullfile(tmpdir,sprintf('step%02d',nSteps),'elastix.log'),log_file)
    end
    % Copy temporary working directory
    if ~isempty(copy_results_dir)        
        copyfile(tmpdir,copy_results_dir)
    end
    % Delete the temporary working directory.
    rmdir(tmpdir,'s')
    
catch ME
    % remove temporary working directory, then throw error message
    rmdir(tmpdir,'s')
    error(getReport(ME))
end


end

