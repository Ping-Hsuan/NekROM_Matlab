%clear all; close all;

hdr;

param = struct;

cmap = colormap(lines);

if exist('ctmp','var') == 0
    [ctmp, mb] = read_tensor('./cu');
else
    fprintf("C tensor is already loaded \n");
end

snap = dlmread("uk");
ns = length(snap)/(mb+1);
snap = reshape(snap,mb+1,ns);

param.csize = ["full","core","skew"];
%param.csize = ["core"];
for j=1:size(param.csize,2);
    nb_list = [100,200,300]
%   nb_list = [100]
    for i=1:size(nb_list,2);
        T1 = table;
        nb = nb_list(i);
        st = "3D turbulent channel flow, "+"$N="+nb+"$";
    
        if param.csize(j) == "full"
            c0 = ctmp(1:nb,1:nb+1,1:nb+1);
            cu  = reshape(ctmp(1:nb,1:nb+1,1:nb+1),nb*(nb+1),nb+1);
        elseif param.csize(j) == "core"
            c0 = ctmp(1:nb,2:nb+1,2:nb+1);
            cu  = reshape(ctmp(1:nb,1:nb+1,1:nb+1),nb*(nb+1),nb+1);
        elseif param.csize(j) == "skew"
            c0 = ctmp(1:nb,2:nb+1,2:nb+1);
            for i=1:nb
                c0(:,1:end,i) = 0.5*(c0(:,1:end,i) - c0(:,1:end,i)');
            end
            cu  = reshape(ctmp(1:nb,1:nb+1,1:nb+1),nb*(nb+1),nb+1);
        end
        
        c_norm = norm(cu);
        c3 = mode_3_unfolding(c0);
        %% How to use mode3 to do tensor contraction
        %% tmp = c3'*snap(1:nb+1,1)
        %% tmp = reshape(tmp,nb,nb+1);
        %% tmpp = tmp*snap(1:nb+1,1)
        [U3,S3,V3] = svd(c3,"econ");
        for r = 1:size(U3,1);
%       for r = 100:100
    
    %       no_approx = reshape(cu*snap(1:nb+1,1),nb,nb+1)*snap(1:nb+1,1);
            app_c = recon_svd(V3,S3,U3,r);
            err_c = app_c-c3;
            relres = norm(err_c)/c_norm;
            perr = max(max((abs(err_c)))/max(max(abs(cu))));
            tt = {nb,r,relres,perr}
            T1 = [T1; tt];
    %       approx = reshape(app_c*snap(1:nb+1,1),nb,nb+1)*snap(1:nb+1,1);
    % tmp = ctmp(1:nb,1,1)*1*1
    % tmp= tmp-reshape(ctmp(1:nb,1,1:nb+1),nb,nb+1)*snap(1:nb+1,1);
    % tmp =tmp-reshape(ctmp(1:nb,1:nb+1,1),nb,nb+1)*snap(1:nb+1,1);
        end
        T1.Properties.VariableNames = {'nb','rank','relres','perr'}
        writetable(T1,"N"+nb+"_svd_approx_"+param.csize(j));
    end
end

function [c0, n_floor] = read_tensor(fname)
   t = dlmread(fname);
   n = nthroot(length(t),3);
   n_ceil = ceil(n);
   n_floor = floor(n);
   c0 = reshape(t,n_floor,n_ceil,n_ceil);
end

function X_1 = mode_1_unfolding(X);
    [I,J,K] = size(X);
    X_1 = reshape(X,I,J*K);
end

function X_2 = mode_2_unfolding(X);
    [I,J,K] = size(X);
    for ii=1:K
        X_hat(:,:,ii) = X(:,:,ii)';
    end
    
    X_2 = reshape(X_hat,J,I*K);
end

function X_3 = mode_3_unfolding(X);
    [I,J,K] = size(X);
    for ii=1:K
       X_3(ii,:) = reshape(X(:,:,ii),I*J,1);
    end
end

function approx_tensor = recon_svd(V3,S3,U3,r);
%   approx_tensor = V3(:,1:r)*S3(1:r,1:r)'*U3(:,1:r)';
    approx_tensor = U3(:,1:r)*S3(1:r,1:r)*V3(:,1:r)';
end
