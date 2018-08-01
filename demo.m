clear vars
close all
path(pathdef);
path(path, './utils/')
run vlfeat-0.9.21/toolbox/vl_setup;

newCaseFlag = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% You can try to generate a random seed and save it to the disk; the
%%% funny thing is that sometimes the effect of resetting the random seed
%%% does not "propagate" into the LDMM class...... use with caution!!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if newCaseFlag
    rng('default')
    seed = randi(1000000000);
    disp(['seed = ' num2str(seed)]);
    rng(seed)
else
    load('seed.mat');
end

OSZ2016 = LDMM('data/barbara_256.png');
OSZ2016.downsampleImg(0.1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% the following are the default settings
%%%% use LDMM.setOptionPairs() to modify option fields
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LDMM.setPatchOptions(struct('patchSize', struct('x', 10, 'y', 10),...
%                             'patchStride', struct('x', 1, 'y', 1)));
% LDMM.setLDMMOptions(struct('rwType', 'None',...
%                            'innerLoop', 1,...
%                            'outerLoop', 100,...
%                            'reInit', false,...
%                            'resetMultiplier', false,...
%                            'BregmanLambda', 0.5,...
%                            'verbose', true));
% LDMM.setGLOptions(struct('symLaplacian', false,...
%                          'NN', 50,...
%                          'SelfTuningNN', 25,...
%                          'MaxNumCompUpperBound', 2^8));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
OSZ2016.setOptionPairs('verbose', false);

%%%% use non-symmetric nonlocal Laplacian (default)
tic;
OSZ2016.reconImg();
fprintf('original LDMM took %f seconds\n', toc);
%%%% use symmetric nonlocal Laplaican
OSZ2016_symGL = LDMM(OSZ2016);
OSZ2016_symGL.setOptionPairs('symLaplacian', true);
tic;
OSZ2016_symGL.reconImg();
fprintf('LDMM with symmetric Laplacian took %f seconds\n', toc);

%%%% use rw-LDMM with SVD basis
YGLD2016_SVD = LDMM(OSZ2016);
YGLD2016_SVD.setOptionPairs('rwType', 'SVD', 'MaxNumCompUpperBound', 2^8);
% YGLD2016_SVD.setOptionPairs('rwType', 'SVD', 'MaxNumCompUpperBound', 2^8, 'SelfTuningNN', 25);
tic;
YGLD2016_SVD.reconImg();
fprintf('rwLDMM-SVD took %f seconds\n', toc);
%%%% use symmetric nonloal Laplaican
YGLD2016_SVD_symGL = LDMM(YGLD2016_SVD);
YGLD2016_SVD_symGL.setOptionPairs('symLaplacian', true);
tic;
YGLD2016_SVD_symGL.reconImg();
fprintf('rwLDMM-SVD with symmetric Laplacian took %f seconds\n', toc);

%%%% use rw-LDMM with DCT basis
YGLD2016_DCT= LDMM(OSZ2016);
YGLD2016_DCT.setOptionPairs('rwType', 'DCT', 'MaxNumCompUpperBound', 2^8);
tic;
YGLD2016_DCT.reconImg();
fprintf('rwLDMM-DCT took %f seconds\n', toc);
%%%% use symmetric nonloal Laplaican
YGLD2016_DCT_symGL = LDMM(YGLD2016_DCT);
YGLD2016_DCT_symGL.setOptionPairs('symLaplacian', true);
tic;
YGLD2016_DCT_symGL.reconImg();
fprintf('rwLDMM-DCT with symmetric Laplacian took %f seconds\n', toc);

%%
%%%%%%% make summary plots
figure;

subplot(2,4,1);
imshow(OSZ2016.img, []);
title('Original');
subplot(2,4,5);
imshow(OSZ2016.sampledImg, []);
title([num2str(OSZ2016.downsampleRate*100) '%  Subsampled']);

subplot(2,4,2);
imshow(OSZ2016.LDMMResult.img_rec{end}, []);
title(sprintf('OSZ2016, PSNR = %0.2f', OSZ2016.LDMMResult.PSNR_list(end)), 'Interpreter', 'none');
subplot(2,4,6);
imshow(OSZ2016_symGL.LDMMResult.img_rec{end}, []);
title(sprintf('OSZ2016_symGL, PSNR = %0.2f', OSZ2016_symGL.LDMMResult.PSNR_list(end)), 'Interpreter', 'none');

subplot(2,4,3);
imshow(YGLD2016_SVD.LDMMResult.img_rec{end}, []);
title(sprintf('YGLD2016_SVD, PSNR = %0.2f', YGLD2016_SVD.LDMMResult.PSNR_list(end)), 'Interpreter', 'none');
subplot(2,4,7);
imshow(YGLD2016_SVD_symGL.LDMMResult.img_rec{end}, []);
title(sprintf('YGLD2016_SVD_symGL, PSNR = %0.2f', YGLD2016_SVD_symGL.LDMMResult.PSNR_list(end)), 'Interpreter', 'none');

subplot(2,4,4);
imshow(YGLD2016_DCT.LDMMResult.img_rec{end}, []);
title(sprintf('YGLD2016_DCT, PSNR = %0.2f', YGLD2016_DCT.LDMMResult.PSNR_list(end)), 'Interpreter', 'none');
subplot(2,4,8);
imshow(YGLD2016_DCT_symGL.LDMMResult.img_rec{end}, []);
title(sprintf('YGLD2016_DCT_symGL, PSNR = %0.2f', YGLD2016_DCT_symGL.LDMMResult.PSNR_list(end)), 'Interpreter', 'none');

%%
%%%%%%% save results in .mat format
save('seed.mat','seed');

save('OSZ2016.mat', 'OSZ2016');
save('OSZ2016_symGL.mat', 'OSZ2016_symGL');

save('YGLD2016_SVD.mat', 'YGLD2016_SVD');
save('YGLD2016_SVD_symGL.mat', 'YGLD2016_SVD_symGL');

save('YGLD2016_DCT.mat', 'YGLD2016_DCT');
save('YGLD2016_DCT_symGL.mat', 'YGLD2016_DCT_symGL');
