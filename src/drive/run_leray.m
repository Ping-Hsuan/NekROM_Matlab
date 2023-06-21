global ifrom nsteps iostep astep dt nb nu ips ifleray ifcopt;
global mb

if isfile('mor.json')
   % File exists.
   fprintf('reading mor.json... \n');
   nekrom_utils.readconf('mor.json')
else
   % File does not exist.
   error('did not find mor.json, exiting ...');
end

mkdir(date);
CTUs =num2str(nsteps*dt)+"CTUs"
outputdir = fullfile(date,CTUs);
mkdir(outputdir);

global relax dfRadius
% Parameter space
% m = {1, 2, 3, 4}
% chi \in [0.01, 20]
% delta \in [h, 1], where h = 0.001
radius_list = [];
radius_list = [linspace(0.001,0.01,10)'; linspace(0.01,0.1,25)'; linspace(0.1,1,10)'];
radius_list = unique(radius_list,'rows');

roms = ["leray"];

global ifleray ifefr iftr filterType
global dfOrder dfRadius

%nb_list = [10, 50, 100];
nb_list = [10, 20, 30, 40, 60, 70, 80, 90, 100];
for ii=1:size(nb_list,2)
    nb = nb_list(ii)
    for order=2:4
    T = table;
    T1 = table;
    T2 = table;
    dfOrder=order
    for kk=1:size(roms,2)
        ifleray = 0;
        ifefr = 0;
        iftr = 0;
        if roms(kk) == "leray"
           ifleray = 1;
        elseif roms(kk) == "EFR"
           ifefr = 1;
        elseif roms(kk) == "TR"
           iftr = 1;
        end
        ifleray
        ifefr
        iftr
        for nn=1:size(radius_list,1)
            relax = 1;
            dfRadius = radius_list(nn);


            global au0 bu0 cu u0;
            global au bu;
            global uk;
            global u ua u2a ru eu hufac;
            
            nekrom_utils.load_ops(ifrom,ifcopt,"./ops/");
            
            nekrom_utils.initvars()
            
            global time
            global alphas betas
            
            global ifdrag
            if ifdrag
               global fd1 fd3 rdgx rdgy;
            end

            if exist('t1','var') == 0 || exist('t2','var')== 0 || exist('t3','var') ==0
               t1 = zeros(mb+1,mb+1,128);
               t2 = zeros(mb+1,mb+1,128);
               t3 = zeros(mb+1,mb+1,128);
               idx = 0;
               for j=1:mb+1
                  for i=1:mb+1
                     tmp = dlmread("qoi/bascor"+idx+".dat");
                     t1(i,j,:) = tmp(:,3);
                     t2(i,j,:) = tmp(:,4);
                     t3(i,j,:) = tmp(:,5);
                     idx = idx+1;
                  end
               end
               fprintf("done loading operators for Reynolds stress \n");
            end
            
            ucoef=zeros((nsteps/iostep),nb+1);
            tcoef=zeros((nsteps/iostep),nb+1);
            
            for istep=1:nsteps
               ito=min(istep,3);
               if istep<= 3
                  htfac = [];
                  hufac = [];
               end
               if ifrom(1)
                  f=zeros(nb,1);
                  [ru, eu] = nekrom_utils.setr(au0,bu,cu,f,u,u,eu,nu,alphas,betas,dt,ito);
                  aup=au;
                  bup=bu;
                  rup=ru;
                  if ifcopt
                  else
                     [u_new,hufac] = nekrom_utils.step(rup,aup,bup,nu,betas,dt,ito,hufac);
                  end
               end
            
               if ifefr
                  global relax;
                  utmp = u_new;
                  if filterType == "cutoff"
                     global filterModes
                     utmp(1,(nb+1)-filterModes+1:end) = 0;
                  elseif filterType == "df" || filterType == "df_repeat" 
                     global dfHfac;
                     utmp = [1,(dfHfac\(dfHfac'\u_new(1,2:end)'))'];
                  end
                  u_new = (1-relax)*u_new + relax*utmp;
               end
            
               time = time+dt;
               if ifrom(1)
                  u = nekrom_utils.shift(u,u_new,3);
                  if istep >= astep
                     ua=ua+u_new';
                     u2a=u2a+u_new'*u_new;
                     if istep==5000 || istep==10000
                        tmpua=ua/(istep-astep+1);
                        tmpu2a=u2a/(istep-astep+1);

                        fileID = fopen(dir+"/ua_"+istep+"_"+dt,'w');
                        fprintf(fileID,"%24.15e\n",tmpua);
                        fclose(fileID);
                        
                        fileID = fopen(dir+"/u2a_"+istep+"_"+dt,'w');
                        fprintf(fileID,"%24.15e\n",tmpu2a);
                        fclose(fileID);
                        if istep==5000
                            [err1,err2,err3,err4] = c_statistics(tmpua,tmpu2a,t1,t2,t3,istep,dt,T1,nb,dir)
                            tt = {roms(kk),nb,relax,dfOrder,dfRadius,err1,err2,err3,err4};
                            T1 = [T1;tt]
                        elseif istep==10000
                            [err1,err2,err3,err4] = c_statistics(tmpua,tmpu2a,t1,t2,t3,istep,dt,T2,nb,dir)
                            tt = {roms(kk),nb,relax,dfOrder,dfRadius,err1,err2,err3,err4};
                            T2 = [T2;tt]
                        end
                    end
                     if istep==nsteps
                        ua=ua/(nsteps-astep+1);
                        u2a=u2a/(nsteps-astep+1);
                        fileID = fopen(dir+"/ua_"+istep+"_"+dt,'w');
                        fprintf(fileID,"%24.15e\n",tmpua);
                        fclose(fileID);
                        
                        fileID = fopen(dir+"/u2a_"+istep+"_"+dt,'w');
                        fprintf(fileID,"%24.15e\n",tmpu2a);
                        fclose(fileID);
                     end
                  end
               end
            
               if ifdrag
                  nekrom_utils.cdrag(fd1,fd3,rdgx,rdgy,u,nu);
               end
            
               if (mod(istep,iostep) == 0)
                  ucoef(istep/iostep,:)=u(:,1);
               end
            end
            
            casedir= sprintf('%s_%d_%d_%.8f_%.8f',roms(kk),nb,dfOrder,dfRadius,relax)
            dir = fullfile(outputdir,casedir);
            mkdir(dir);

            fileID = fopen(dir+"/ua",'w');
            fprintf(fileID,"%24.15e\n",ua);
            fclose(fileID);
            
            fileID = fopen(dir+"/u2a",'w');
            fprintf(fileID,"%24.15e\n",u2a);
            fclose(fileID);
            
            fileID = fopen(dir+"/ucoef",'w');
            fprintf(fileID,"%24.15e\n",ucoef);
            fclose(fileID);
            
            % mesh plot the coefficients
            i = repmat(linspace(1,nb,nb)',1,nsteps/iostep);
            t = repmat(linspace(1,nsteps*dt,nsteps/iostep),nb,1);
            figure; mesh(i,t,ucoef(:,2:end)'); hold on;
            xlabel('j'); ylabel('t'); zlabel('a\_j');
            saveas(gcf,dir+"/coef.png");
            
            % reconstruct mean profile
            uam = zeros(128,5);
            for j=1:nb+1
               tmp = dlmread("qoi/basavg"+(j-1)+".dat");
               uam(:,3) = uam(:,3) + ua(j)*tmp(:,3);
               uam(:,4) = uam(:,4) + ua(j)*tmp(:,4);
               uam(:,5) = uam(:,5) + ua(j)*tmp(:,5);
            end
            uam(:,2) = tmp(:,2);
            
            uu = zeros(128,1);
            vv = zeros(128,1);
            ww = zeros(128,1);
            for n=1:128
               uu(n)=sum(sum(u2a.*t1(1:nb+1,1:nb+1,n)))-ua'*t1(1:nb+1,1:nb+1,n)*ua;
               vv(n)=sum(sum(u2a.*t2(1:nb+1,1:nb+1,n)))-ua'*t2(1:nb+1,1:nb+1,n)*ua;
               ww(n)=sum(sum(u2a.*t3(1:nb+1,1:nb+1,n)))-ua'*t3(1:nb+1,1:nb+1,n)*ua;
            end
            fileID = fopen(dir+"/uam",'w');
            fprintf(fileID,"%24.15e\n",uam);
            fclose(fileID);
            fileID = fopen(dir+"/uu",'w');
            fprintf(fileID,"%24.15e\n",uu);
            fclose(fileID);
            fileID = fopen(dir+"/vv",'w');
            fprintf(fileID,"%24.15e\n",vv);
            fclose(fileID);
            fileID = fopen(dir+"/ww",'w');
            fprintf(fileID,"%24.15e\n",ww);
            fclose(fileID);
            
            % Load FOM Reynolds stress
            fomavg   = dlmread("fom_100CTU/avg.list");
            fomavg2  = dlmread("fom_100CTU/squa.list");
            fomrm2   = dlmread("fom_100CTU/rm2.list");
            fomrms   = dlmread("fom_100CTU/rms.list");
            fomufric = dlmread("fom_100CTU/u_fric");
            npts = 128; mid = (npts/2)-1;
            
            fn = sprintf('%d_%d_%.8f_%.8f',nb,dfOrder,dfRadius,relax);
            lb = sprintf('$N=%d,~m=%d,~\\delta=%.8f,~\\chi=%.8f$',nb,dfOrder,dfRadius,relax);
            
            y = fomavg(:,2);
            figure; plot(y,fomavg(:,3),'k-'); hold on; plot(y,uam(:,3),'b-');
            legend({'FOM',lb},'Interpreter','latex','FontSize',14);
            ylabel("$\langle\overline{u}\rangle_{xz}$",'Interpreter','latex','fontsize',16); xlabel('$y$','Interpreter','latex','fontsize',16);
            saveas(gcf,dir+"/"+fn+"_uavg.png")
            
            err1 = norm(fomavg(:,3)-uam(:,3))/norm(fomavg(:,3));
            
            fomuu = fomrms(:, 3) - fomavg2(:,3);
            figure; plot(y,fomuu,'k-'); hold on; plot(y,uu,'b-');
            legend({'FOM',lb},'Interpreter','latex','FontSize',14);
            ylabel("$\langle\overline{u'u'}\rangle_{xz}$",'Interpreter','latex','fontsize',16); xlabel('$y$','Interpreter','latex','fontsize',16);
            saveas(gcf,dir+"/"+fn+"_uu.png")
            err2 = norm(fomuu-uu)/norm(fomuu);
            
            fomvv = fomrms(:, 4) - fomavg2(:,4);
            figure; plot(y,fomvv,'k-'); hold on; plot(y,vv,'b-');
            legend({'FOM',lb},'Interpreter','latex','FontSize',14);
            ylabel("$\langle\overline{v'v'}\rangle_{xz}$",'Interpreter','latex','fontsize',16); xlabel('$y$','Interpreter','latex','fontsize',16);
            saveas(gcf,dir+"/"+fn+"_vv.png")
            err3 = norm(fomvv-vv)/norm(fomvv);
            
            fomww = fomrms(:, 5) - fomavg2(:,5);
            figure; plot(y,fomww,'k-'); hold on; plot(y,ww,'b-');
            legend({'FOM',lb},'Interpreter','latex','FontSize',14);
            ylabel("$\langle\overline{w'w'}\rangle_{xz}$",'Interpreter','latex','fontsize',16); xlabel('$y$','Interpreter','latex','fontsize',16);
            saveas(gcf,dir+"/"+fn+"_ww.png")
            err4 = norm(fomww-ww)/norm(fomww);
            
            close all;
            fprintf('%s, nb %d # %d relax, oder, radius, error: %f %d %f %f %f %f %f \n',roms(kk),nb,nn,relax,dfOrder,dfRadius,err1,err2,err3,err4);
            tt = {roms(kk),nb,relax,dfOrder,dfRadius,err1,err2,err3,err4};
            T = [T;tt];
            end
        end
        T.Properties.VariableNames = {'model','nb','relax','order','radius','uavg_err','uu_err','vv_err','ww_err'}
        writetable(T, outputdir+"/table_m"+order+"_leray_N"+nb);
        T1.Properties.VariableNames = {'model','nb','relax','order','radius','uavg_err','uu_err','vv_err','ww_err'}
        writetable(T1,outputdir+"/table_m"+order+"_efr_N"+nb+"_25CTU");
        T2.Properties.VariableNames = {'model','nb','relax','order','radius','uavg_err','uu_err','vv_err','ww_err'}
        writetable(T2,outputdir+"/table_m"+order+"_efr_N"+nb+"_50CTU");
    end
end
