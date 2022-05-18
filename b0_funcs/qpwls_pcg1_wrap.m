function [x_recon] = qpwls_pcg1_wrap(n_v,A,W,kdata,C,iter)
% wrapper function for qpwls_pcg1 to output 2D image

x_recon = reshape(qpwls_pcg1(zeros(prod(n_v),1), A, W, kdata(:), C, 'niter', iter),n_v);

end

