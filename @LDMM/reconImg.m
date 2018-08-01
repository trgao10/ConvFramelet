function LDMMResult = reconImg(obj)
%RECONIMG Reconstruct the input image from downsampled pixels
%         Function wrapper of the code provided by Z. Shi and R. Yin
%
%   LDMMResult: structure containing outputs of the reconstruction
%       PSNR_list: vector of length obj.LDMMOptions.outerLoop, containing
%                  the PSNR values throughout the iterations
%       img_rec:   cell array of reconstructed images
%
%
%   last modified by Tingran Gao (trgao10@math.duke.edu)
%   last modified on June 23, 2016
%

%%%%%% extract parameters from patchOptions and LDMMOptions
rwType = LDMM.getoptions(obj.LDMMOptions, 'rwType', 'NONE');
outerLoop = LDMM.getoptions(obj.LDMMOptions, 'outerLoop', 100);
innerLoop = LDMM.getoptions(obj.LDMMOptions, 'innerLoop', 1);
lambda = LDMM.getoptions(obj.LDMMOptions, 'BregmanLambda', 0.5);
resetMultiplier = LDMM.getoptions(obj.LDMMOptions, 'resetMultiplier', false);
reInit = LDMM.getoptions(obj.LDMMOptions, 'reInit', false);
verbose = LDMM.getoptions(obj.LDMMOptions, 'verbose', true);
clampPixelVal = LDMM.getoptions(obj.LDMMOptions, 'clampPixelVal', false);
trackDims = LDMM.getoptions(obj.LDMMOptions, 'trackDims', false);
MIN_PIXEL_VAL = min(LDMM.getoptions(obj.LDMMOptions, 'pixelValueBound', [0,255]));
MAX_PIXEL_VAL = max(LDMM.getoptions(obj.LDMMOptions, 'pixelValueBound', [0,255]));

%%%%%% randomly fit in the missing pixels to generate an initial guess
if reInit
    fw = zeros(obj.imgSize);
    fw(obj.downsampledPixelID) = obj.downsampledPixel;
    id0 = setdiff(1:prod(obj.imgSize), obj.downsampledPixelID);
    fw(id0) = mean(double(obj.downsampledPixel))...
        + std(double(obj.downsampledPixel)) * randn(size(id0)); % fill in the missing pixels with random numbers
else
    fw = obj.initialImg;
end

[X,Y] = meshgrid(1:obj.patchOptions.patchStride.x:obj.imgSize(1), 1:obj.patchOptions.patchStride.y:obj.imgSize(2));

F = LDMM.image2patch(fw, X(:), Y(:),...
                     obj.patchOptions.patchSize.x,...
                     obj.patchOptions.patchSize.y);
                 
%%%%%% N = total number of patches, L = patch length
[N, l] = size(F);

LDMMResult.PSNR_list = zeros(outerLoop,1);
LDMMResult.img_rec = cell(outerLoop,1);
if trackDims
    LDMMResult.dims_list = zeros(outerLoop,1);
end

uw = F;
v = 0;
% u_image = zeros(obj.imgSize);
u_image = LDMM.patch2image(uw+v, X(:), Y(:),...
                           obj.patchOptions.patchSize.x,...
                           obj.patchOptions.patchSize.y,...
                           obj.imgSize(1), obj.imgSize(2)); % recover the 2D data from patches              
u_image(obj.downsampledPixelID) = obj.downsampledPixel; % assign the value at sample points to be the given value
if verbose
    fprintf('Step %d, PSNR=%f\n', 0, LDMM.psnr_ldmm(u_image, double(obj.img), MAX_PIXEL_VAL)); % display the PSNR at current step
end

%%%%%% if rw-LDMM + DCT, build fixed basis V
if strcmpi(rwType, 'dct')
    %%%%%% if rw-LDMM + SVD, set svd opt
    %     opts = struct('tol', 1e-8, 'maxit', 150);
    V = gen_DCT2dic(obj.patchOptions.patchSize.x, obj.patchOptions.patchSize.y);
elseif strcmpi(rwType, 'svd')
elseif strcmpi(rwType, 'le')
end

uw_old = zeros(size(uw));

for ii = 1:outerLoop
    W = LDMM.getNonlocalLaplacian(uw, obj.GLOptions);
    if trackDims
        for kk=1:size(uw,2)
            LDMMResult.dims_list(ii) = LDMMResult.dims_list(ii) + uw(:,kk)'*W*uw(:,kk);
        end
    end
    
    %%%%% if rw-LDMM, prepare weights and local basis V
    if strcmpi(rwType, 'SVD')
        r0 = ceil(0.2*obj.patchOptions.patchSize.x*obj.patchOptions.patchSize.y);% take top 20% singular values
        [~,S,V] = RandPCA(uw, r0);
%         [~,S,V] = svds(uw, r0, 'L', opts);
        S = diag(S)';
        S = 1-S/max(S);
    elseif strcmpi(rwType, 'DCT')
        r0 = ceil(.2*obj.patchOptions.patchSize.x*obj.patchOptions.patchSize.y);% take top 20% singular values
        S = 1 - gen_weight(uw, V);
        [S, idx] = sort(S, 'ascend');
        S = S(1:r0);
        V = V(:, idx(1:r0));
    elseif strcmpi(rwType, 'LE')
        r0 = ceil(0.2*obj.patchOptions.patchSize.x*obj.patchOptions.patchSize.y);% take top 20% singular values
        Ws = LDMM.getNonlocalLaplacian(uw', struct('symLaplacian', true, 'NN', 10, 'SelfTuningNN', 5, 'MaxNumCompUpperBound', 2^6));
        Ws = (Ws + Ws')/2;
        invSqrtDs = diag(1./sqrt(sum(Ws)));
        [V, ~] = eigs(invSqrtDs*Ws*invSqrtDs, r0, 'lm', struct('isreal',1,'issym',1,'maxit',1000,'disp',0));
        S = 1 - gen_weight(uw, V);
        [S, idx] = sort(S, 'ascend');
        S = S(idx(1:r0));
        V = V(:, idx(1:r0));
    end
    
    if (resetMultiplier>0) && (mod(ii, resetMultiplier) == 0)
        v = 0;
    end

    for kk=1:innerLoop
        b = lambda*W*(uw-v);
        
        if strcmpi(rwType, 'SVD') || strcmpi(rwType, 'DCT') || strcmpi(rwType, 'LE')
            bt = b*V;
            uwt = uw*V;
            for i=1:r0
                coe_matrix = S(i)*(sparse(1:N, 1:N, sum(W,2), N, N) - W) + lambda*W;
                [uw_old(:,i),~,~] = gmresm(coe_matrix, bt(:,i), [], 1e-2, 50, [], [], uwt(:,i));
            end
            if r0 < size(uw,2) %fast
                % compute the completement of the projection
                b = b - bt*V';
                uw = uw - uwt*V';
                coe_matrix = sparse(1:N, 1:N, sum(W,2), N, N) - W + lambda*W;
                [uw, ~, ~] = gmresm(coe_matrix, b,[],1e-2,50,[],[],uw);
                uw = uw + uw_old(:,1:r0)*V';
            else
                uw = uw_old(:,1:r0)*V';
            end
        elseif strcmpi(rwType, 'NONE')
            coe_matrix = sparse(1:N, 1:N, sum(W,2), N, N) - W + lambda*W;
            [uw, ~, ~] = gmresm(coe_matrix,b,[],1e-3,50,[],[],uw); % solving the linear system using (vectorized) GMRES
        else
            error(['undefined rwType: ' rwType]);
        end
        
        uw_old = uw;
        u_image = LDMM.patch2image(uw+v, X(:), Y(:),...
                                   obj.patchOptions.patchSize.x,...
                                   obj.patchOptions.patchSize.y,...
                                   obj.imgSize(1), obj.imgSize(2)); % recover the 2D data from patches      

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%% experiment: clamping pixel values in u_image
        if clampPixelVal
            u_image = min(MAX_PIXEL_VAL, max(MIN_PIXEL_VAL, u_image));
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        u_image(obj.downsampledPixelID) = obj.downsampledPixel; % assign the value at sample points to be the given value
            
        uw = LDMM.image2patch(u_image, X(:), Y(:),...
                              obj.patchOptions.patchSize.x,...
                              obj.patchOptions.patchSize.y);
                                 
        v = v+uw_old-uw; % update the Lagrangian multiplier

        LDMMResult.PSNR_list(ii) = LDMM.psnr_ldmm(min(MAX_PIXEL_VAL, max(MIN_PIXEL_VAL, u_image)), double(obj.img), MAX_PIXEL_VAL);
%         LDMMResult.PSNR_list(ii) = LDMM.psnr_ldmm(u_image, double(obj.img), MAX_PIXEL_VAL);
        if verbose
            fprintf('Step %d, PSNR=%f\n', ii, LDMMResult.PSNR_list(ii)); % display the PSNR at current step
        end
    end
    LDMMResult.img_rec{ii} = min(MAX_PIXEL_VAL, max(MIN_PIXEL_VAL, u_image));
%     LDMMResult.img_rec{ii} = u_image; % save the data
end

LDMMResult.img = obj.img;
LDMMResult.sampledImg = obj.sampledImg;


obj.LDMMResult = LDMMResult;

end

