classdef LDMM < handle
    %LDMM Class implementing image inpainting algorithms based on Low
    %     Dimensional Manifold Model (LDMM) and Reweighted Low Dimensional
    %     Manifold Model (rw-LDMM)
    %
    %   [1] Stanley Osher, Zuoqiang Shi, and Wei Zhu.
    %       Low Dimensional  Manifold Model for Image Processing.
    %       UCLA, Tech. Rep. CAM report 16-04, 2016.
    %   [2] Rujie Yin, Tingran Gao, Yue Lu, and Ingrid Daubechies.
    %       Regularization on Image Patches: a linear reconstruction from manifold embedding.
    %       arXiv preprint arXiv:1606.01377 (2016).
    %   
    %   The LDMM class properties:
    %       img ---------------------------------------- orignal image
    %       imgSize ------------------------------------ size of the original image, a vector of length 2
    %       downsampleRate  ---------------------------- downsample rate
    %       downsampledPixel --------------------------- downsampled image
    %       downsampledPixelID ------------------------- IDs of the downsampled pixels, index by the linear order in matrices
    %       sampledImg --------------------------------- downsampled image (missing pixels are set to gray scale 0)
    %       mask --------------------------------------- 0-1 matrix encoding pixels sampled
    %       patchOptions
    %                patchSize ------------------------- .x, .y
    %                patchStride ----------------------- .x, .y
    %       LDMMOptions
    %                rwType ---------------------------- { ['None'] | 'SVD' | 'DCT'}
    %                innerLoop ------------------------- { [1] | (any positive integer) }
    %                outerLoop ------------------------- { [100] | (any positive integer) }
    %                resetMultiplier ------------------- { [false] | true } if set to true, 
    %                                                       Bregman iteration will reset v=0 in the
    %                                                       outerloop every time before entering the inner loop
    %                BregmanLambda --------------------- { [0.5] | (between 0 and 1) }
    %                pixelValueBound ------------------- { [ [0, 255] ] | (any integer pair) }
    %                clampPixelVal --------------------- { [false] | true }
    %                trackDims ------------------------- { [false] | true }
    %       graphLaplacianOptions
    %                symLaplacian ---------------------- { [false] | true}
    %                NN -------------------------------- { [50] | (any positive integer) }
    %                                                      number of nearest neighbors
    %                SelfTuningNN ---------------------- { [20] | (any positive integer < NN)}
    %                MaxNumCompUpperBound -------------- { [2^10] | (any positive integer > NN) }
    %                                                      controls maximum number of comparisions
    %                                                      carried out in randomized kdtree search
    %       
    %
    %   Dependency: vlfeat (for kdtree)
    %
    %
    %   last modified by Tingran Gao (trgao10@math.duke.edu)
    %   last modified on June 23, 2016

    properties
        img
        imgSize
        downsampleRate = 0.1;
        downsampledPixel
        downsampledPixelID
        sampledImg
        mask
        initialImg
        patchOptions = struct('patchSize', struct('x', 10, 'y', 10),...
                              'patchStride', struct('x', 1, 'y', 1));
        LDMMOptions = struct(...
            'rwType', 'None',...
            'innerLoop', 1,...
            'outerLoop', 100,...
            'reInit', false,...
            'verbose', true,...
            'resetMultiplier', false,...
            'BregmanLambda', 0.5,...
            'pixelValueBound', [0,255],...
            'clampPixelVal', false,...
            'trackDims', false);
        GLOptions = struct(...
            'symLaplacian', false,...
            'NN', 50,...
            'SelfTuningNN', 20,...
            'MaxNumCompUpperBound', 2^8)
        LDMMResult
    end
    
    methods
        function obj = LDMM(rhs)
            if nargin > 0
                if isa(rhs, 'LDMM')
                    fns = properties(rhs);
                    for i=1:length(fns)
                        obj.(fns{i}) = rhs.(fns{i});
                    end
                else
                    obj = obj.readImg(rhs);
                end
            end
        end
        
        function obj = readImg(obj,imgFileName)
            try
                obj.img = imread(imgFileName);
            catch
                error(['Image file ' imgFileName ' does not exist!']);
            end
            if ~ismatrix(obj.img)
                obj.img = rgb2gray(obj.img);
            end
            obj.imgSize = size(obj.img);
        end
        
        function downsampleImg(obj, downsampleRate)
            if isempty(obj.img)
                error('no image in LDMM class');
            end
            assert(downsampleRate >= 0 && downsampleRate <= 1, 'downsampleRate must be set within range [0, 1]');
            if nargin == 1
                downsampleRate = obj.downsampleRate;
            else
                obj.downsampleRate = downsampleRate;
            end
            [obj.downsampledPixel, obj.downsampledPixelID] = datasample(obj.img(:), floor(numel(obj.img)*downsampleRate), 'Replace', false);
            obj.sampledImg = zeros(obj.imgSize);
            obj.sampledImg(obj.downsampledPixelID) = obj.downsampledPixel;
            [I,J] = ind2sub(obj.imgSize, obj.downsampledPixelID);
            obj.mask = sparse(I, J, ones(size(I)), obj.imgSize(1), obj.imgSize(2));
            
            %%%% prepare an initial guess    
            obj.initializeDownsampledImg();
        end
        
        function initializeDownsampledImg(obj)
            obj.initialImg = obj.sampledImg;
            id0 = setdiff(1:prod(obj.imgSize), obj.downsampledPixelID);
            obj.initialImg(id0) = mean(double(obj.downsampledPixel))...
                                + std(double(obj.downsampledPixel)) * randn(size(id0)); % fill in the missing pixels with random numbers
        end
        
        function setPatchOptions(obj, patchOptions)
            obj.patchOptions = patchOptions;
        end
        
        function setLDMMOptions(obj, LDMMOptions)
            obj.LDMMOptions = LDMMOptions;
        end
        
        function setGLOptions(obj, GLOptions)
            obj.GLOptions = GLOptions;
        end
        
        function setOptionPairs(obj, varargin)
            while length(varargin) >= 2
                if isfield(obj.patchOptions, varargin{1})
                    obj.patchOptions = setfield(obj.patchOptions, varargin{1}, varargin{2});
                    disp(['set patchOptions.' varargin{1} ':']);
                    disp(varargin{2});
                elseif isfield(obj.LDMMOptions, varargin{1})
                    obj.LDMMOptions = setfield(obj.LDMMOptions, varargin{1}, varargin{2});
                    disp(['set LDMMOptions.' varargin{1} ' as ' num2str(varargin{2})]);
                elseif isfield(obj.GLOptions, varargin{1})
                    obj.GLOptions = setfield(obj.GLOptions, varargin{1}, varargin{2});
                    disp(['set GLOptions.' varargin{1} ' as ' num2str(varargin{2})]);
                else
                    warning(['option field ' varargin{1} ' does not exist! ignored....']);
                end
                if length(varargin) < 3
                    break
                else
                    varargin = varargin(3:end);
                end
            end
        end
        
        function resetLDMMOptions(obj)
            obj.LDMMOptions = struct(...
                'rwType', 'None',...
                'innerLoop', 1,...
                'outerLoop', 100,...
                'reInit', false,...
                'resetMultiplier', false,...
                'BregmanLambda', 0.5,...
                'pixelValueBound', [0,255],...
                'clampPixelVal', false,...
                'trackDims', false);
        end
        
        function makeMask(obj)
            [I,J] = ind2sub(obj.imgSize, obj.downsampledPixelID);
            obj.mask = sparse(I, J, ones(size(I)), obj.imgSize(1), obj.imgSize(2));
        end
        
    end
    
    methods (Static)
        v = getoptions(options, name, v_default)
        y = image2patch(f, x1, x2, px, py)
        y = patch2image(f, x1, x2, px, py, mi, ni)
        L = getNonlocalLaplacian(data, graphLaplaicanOptions)
        y = psnr_ldmm(f, g, MAX_PIXEL_VAL)
        [x,flag,relres,iter,resvec] = gmresm(A,b,restart,tol,maxit,M1,M2,x,varargin)
    end    
end

