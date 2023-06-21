function [err1, err2, err3, err4] = c_statistics(tmpua,tmpu2a,t1,t2,t3,istep,dt,T,nb,dir)
    % reconstruct mean profile
    uam = zeros(128,5);
    for j=1:nb+1
       tmp = dlmread("qoi/basavg"+(j-1)+".dat");
       uam(:,3) = uam(:,3) + tmpua(j)*tmp(:,3);
       uam(:,4) = uam(:,4) + tmpua(j)*tmp(:,4);
       uam(:,5) = uam(:,5) + tmpua(j)*tmp(:,5);
    end
    uam(:,2) = tmp(:,2);
    
    uu = zeros(128,1);
    vv = zeros(128,1);
    ww = zeros(128,1);
    for n=1:128
       uu(n)=sum(sum(tmpu2a.*t1(1:nb+1,1:nb+1,n)))-tmpua'*t1(1:nb+1,1:nb+1,n)*tmpua;
       vv(n)=sum(sum(tmpu2a.*t2(1:nb+1,1:nb+1,n)))-tmpua'*t2(1:nb+1,1:nb+1,n)*tmpua;
       ww(n)=sum(sum(tmpu2a.*t3(1:nb+1,1:nb+1,n)))-tmpua'*t3(1:nb+1,1:nb+1,n)*tmpua;
    end
    
    fileID = fopen(dir+"/uam_"+istep+"_"+dt,'w');
    fprintf(fileID,"%24.15e\n",uam);
    fclose(fileID);
    fileID = fopen(dir+"/uu_"+istep+"_"+dt,'w');
    fprintf(fileID,"%24.15e\n",uu);
    fclose(fileID);
    fileID = fopen(dir+"/vv_"+istep+"_"+dt,'w');
    fprintf(fileID,"%24.15e\n",vv);
    fclose(fileID);
    fileID = fopen(dir+"/ww_"+istep+"_"+dt,'w');
    fprintf(fileID,"%24.15e\n",ww);
    fclose(fileID);
    
    if istep == 5000
        fomdir = "fom_25CTU/"
    elseif istep == 10000
        fomdir = "fom_50CTU/"
    elseif istep == 20000
        fomdir = "fom_100CTU/"
    elseif istep == 100000
        fomdir = "fom_500CTU/"
    else
        error("statistics might not be correct\n")
    end
    % Load FOM Reynolds stress
    fomavg   = dlmread(fomdir+"avg.list");
    fomavg2  = dlmread(fomdir+"squa.list");
    fomrm2   = dlmread(fomdir+"rm2.list");
    fomrms   = dlmread(fomdir+"rms.list");
    fomufric = dlmread(fomdir+"u_fric");
    npts = 128; mid = (npts/2)-1;
    
    y = fomavg(:,2);
    err1 = norm(fomavg(:,3)-uam(:,3))/norm(fomavg(:,3));
    
    fomuu = fomrms(:, 3) - fomavg2(:,3);
    err2 = norm(fomuu-uu)/norm(fomuu);
    
    fomvv = fomrms(:, 4) - fomavg2(:,4);
    err3 = norm(fomvv-vv)/norm(fomvv);
    
    fomww = fomrms(:, 5) - fomavg2(:,5);
    err4 = norm(fomww-ww)/norm(fomww);
    
%   fprintf('%s, nb %d # %d relax, radius: %f %f %f %f %f %f \n',roms(kk),nb,nn,relax,dfOrder,dfRadius,err1,err2,err3,err4);
%   tt = {roms(kk),nb,relax,dfOrder,dfRadius,err1,err2,err3,err4};
%   T = [T;tt];
end
