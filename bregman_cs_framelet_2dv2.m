function [u, errors, res_errors, errors_per_breg, res_errors_per_breg, iters ] = bregman_cs_framelet_2dv2(f, m, n, A, At, nu, exci, mu, lambda, gammah, maxiters, img, video)
%bregman_cs_framelet
%
%We minimize nu*|nabla u| + exci*|Fu| st Au = f
%F is a framelet transform
%
% A -- a function handle to the sampling operator
% At -- A^t
%
% nu -- tune the TV strength
% exci -- tune the framelet strength
%
% mu -- the exterior constraint split
% gamma -- the framelet split
% lambda -- the TV split
% ie mu,lambda, and gamma are parameters for
% mu*At(A(x)) + lambda*laplacian + gamma*x
% which we solve with pcg
%
% maxiter -- how many times we can apply A before stopping the algorithm
%
% img is the exact image, we only use it to update the errors as you go.
%
k=1;
%don't know if this exists already or not...
norm3 = @(x) sqrt(trapz(trapz(trapz(abs(x).^2))));


%1 is for Haar, 2 is for framelet, 3 is for cubic framelet
%Haar sucks. 2 works best, 3 is slower and worse.
[D,Dt]=GenerateFrameletFilter(2);
%more levels is better? I don't know. Experimentally, more levels is slow
%and worse quality. Hmmm.
n_level = 3;

%set initial guesses for framelets
U = FraDecMultiLevel(zeros(m,n),D,n_level);
dw = U;
bw = U;


%check that we chose sane parameters
if nu == 0
    assert(lambda == 0);
end
if exci == 0
    assert(gammah == 0);
end

%create the AtA operator. for some reason lk convolution like this (not conv2) works best.
%lk = zeros(m,n);
%lk(1,1) = 4;lk(1,2)=-1;lk(2,1)=-1;lk(m,1)=-1;lk(1,n)=-1;
%AtA = @(x) mu*At(A(x)) + lambda*reshape(ifft2(fft2(reshape(x,[m n])).*fft2(lk)), [m*n 1]) + gammah*x;

%2D
Lap = laplacian([m n], {'P','P'}); %negative laplacian
AtA = @(x) mu*At(A(x)) + lambda*Lap*x + gammah*x;

%without correction or pcg
%     uker = zeros(m,n,k);
%     uker(1,1,1) = 6;
%     uker(2,1,1) = -1; uker(1,2,1) = -1; 
%     if k>1
%        uker(1,1,2) = -1;
%        uker(1,1,k) = -1;
%     end
%     uker(m,1,1) = -1; uker(1,n,1)  =-1; 
%     uker = mu*(conj(R).*R)+lambda*fftn(uker)+gammah;

%initial guess.
u = reshape(At(f),[m n]);

%initial splitting vectors
dx = zeros(size(u));
dy = zeros(size(u));
bx = zeros(size(u));
by = zeros(size(u));

%constrained, replace f with fl and update after each unconstrained
fl = f;
res_errors = [1];
errors = [];
errors_per_breg = [];
res_errors_per_breg = [];
iters = [];

% start the optimization.
% the outer loop is constraint enforcement
ell = 0;
while(sum(iters) < maxiters)
    ell = ell+1;

    %unconstrained step, 5 is overkill.
    for kay=1:5
        %update u        
        rhsFrame = FraRecMultiLevel(SubFrameletArray(dw,bw),Dt,n_level);
        rhs = mu.*reshape(At(reshape(fl, [numel(fl) 1])),[m n]);
        rhs = rhs + lambda*(Deriv(dx - bx,1,true) + Deriv(dy - by,2,true));
        rhs = rhs + gammah*rhsFrame;

        
        %this is where all the computation happens
        macpcg = 150;
        [u,flagcg,relres,iterc] = pcg(AtA, reshape(rhs,[numel(rhs) 1]), 1e-6, macpcg);
        
        %Tom's way
  %       iterc = 1;
  %       u = ifftn(fftn(rhs)./uker);

        if flagcg==0
            fprintf('pcg successful, iterc: %i, relres: %d \n',[iterc relres]);        
            iters = [iters iterc];

        else
            fprintf('pcg terrible: %i iterations. relres: %d\n',[macpcg relres]);
            iters = [iters macpcg];
        end
        
        u = reshape(u,[m n]);
        
        %update d's
        Fu = FraDecMultiLevel(u,D,n_level);
        U0 = AddFrameletArray(Fu, bw); %Fu + b_k
        if gammah ~= 0 && exci ~= 0
            dw = ShrinkFramelet(U0, 1/(gammah/exci));
        end
        
        if lambda ~=0 && nu ~= 0
            [dx,dy] = shrink2( Deriv(u,1,false)+bx, Deriv(u,2,false)+by, 1/(lambda/nu));
        end
       
        %update b's       
        bw = SubFrameletArray(U0,dw); %Fu-d_k
        bx = bx + (Deriv(u,1,false) - dx);
        by = by + (Deriv(u,2,false) - dy);
        
        res_errors = [res_errors norm(A(reshape(u, [m*n*k 1])) - f,'fro')/norm(f,'fro')];
        errors = [errors norm3(u./max(u(:)) - img./max(img(:)))];
        
    end
    fl = fl + f - A(reshape(u, [m*n 1]));
    res_errors_per_breg = [res_errors_per_breg norm(A(reshape(u, [m*n 1])) - f)/norm(f)];
    errors_per_breg = [errors_per_breg norm(u - img)/norm(img)];
    
    %
    % Make video
    %
    if (video && randi(1) == 1) %don't plot too much.
        if(ndims(u)==3)
        %%
        hh = figure;
        set(hh,'Visible','off');
        
        whitebg('k');       
        h = vol3d('cdata',abs(u),'texture','3D');
        %view(3);
        view(magicview);
        axis tight;  daspect([1 1 .4])
        alphamap('rampup');
        %colormap bone;
        %%
        %alphamap(.05 .* alphamap);
            
        fnameh = sprintf('3drecon_step%i',ell);
        print(hh,'-dpng',fnameh);
        
        close(hh);
        %pause(.03);
        else if (ndims(u)==2)
                subplot(2,2,1);
                imagesc(real(img));
                title('exact image');
        
                subplot(2,2,2);
                imagesc(real(u));
                title('current reconstruction');
        
                subplot(2,2,4);
                errorfig = abs(u/max(u(:)) - img/max(img(:)));
                imagesc(errorfig);
                title('reconstruction error');
        
                subplot(6,2,7);
                semilogy(errors_per_breg);
                title('image domain error per breg');
                subplot(6,2,11);
                semilogy(res_errors_per_breg(2:end),'r.-');
                title('residual error per breg');
        
                colormap hot;
                pause(0.03);
        
            end
        end
    end
    fprintf('breg step ell = %i error (u-img): %f Ax iter numer: %i \n', [ell errors(end) sum(iters)]);
    fprintf('breg step ell = %i res error : %f Ax iter numer: %i \n---\n', [ell res_errors(end) sum(iters)]);
    
end

