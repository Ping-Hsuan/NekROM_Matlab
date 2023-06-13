clear all; close all;

hdr;

%% Load reduced mass matrix for computing L^2 norm
bu = dlmread("bu");
mb = sqrt(length(bu));
bu = reshape(bu,mb,mb);

%% Load projected coefficients of snapshots
snap = dlmread("uk");
ns = length(snap)/mb;
snap = reshape(snap,mb,ns);
msnap = mean(snap,2);

%% Setting up flags
fileList = dir('g-rom*');
if isempty(fileList)
    ifgrom = false;
else
    ifgrom = true;
end

%% Setup ROM parameters
T = 25;

f=figure(1);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'Units', 'Inches', 'Position', [0, 0, fig_width, fig_height],...
    'PaperUnits', 'Inches', 'PaperSize', [fig_width, fig_height]);
f=figure(2);

set(gcf, 'PaperUnits', 'inches');
set(gcf, 'Units', 'Inches', 'Position', [0, 0, fig_width, fig_height],...
    'PaperUnits', 'Inches', 'PaperSize', [fig_width, fig_height])

nb_list = [10, 20, 40, 80, 100, 200, 300];
cmap = colormap(lines);

for ii=1:size(nb_list,2)
    nb=nb_list(ii); ene_snap = table;
    for jj=1:size(snap,2)
        ene_snap = [ene_snap;...
                    {snap(1:nb+1,jj)'*bu(1:nb+1,1:nb+1)*snap(1:nb+1,jj) - ... 
                    2*snap(1:nb+1,jj)'*bu(1:nb+1,1:nb+1)*msnap(1:nb+1) + ... 
                    msnap(1:nb+1)'*bu(1:nb+1,1:nb+1)*msnap(1:nb+1)}];
    end
    figure(1)
    t=linspace(T/size(ene_snap.b,1),T,size(ene_snap.b,1));
    semilogy(t,ene_snap.b,'-',cr,cmap(ii,:),dispname,"Porjection, $N="+nb+"$"); hold on

    nb = nb_list(ii);
    if (ifgrom)
        grom = dlmread("g-rom_"+nb+"/ucoef");
        ndata = length(grom)/(nb+1);
        grom = reshape(grom,ndata,nb+1);
        ua = dlmread("g-rom_"+nb+"/ua");
        ene = table;
        for jj=1:size(grom,1)
            ene = [ene; ...
                   {grom(jj,:)*bu(1:nb+1,1:nb+1)*grom(jj,:)' - ...
                   2*grom(jj,:)*bu(1:nb+1,1:nb+1)*ua + ...
                   ua'*bu(1:nb+1,1:nb+1)*ua}];
        end
        figure(1)
        t=linspace(T/size(ene.b,1),T,size(ene.b,1));
        semilogy(t,ene.b,'-.',cr,cmap(ii,:),dispname,"G-ROM, $N="+nb+"$"); hold on

        if (ndata == ns)
            figure(2)
            semilogy(t,abs(ene_snap.b-ene.b)./ene_snap.b,...
                     'x',cr,cmap(ii,:),dispname,"$N="+nb+"$",ms,2); hold on
        elseif (ndata ~= ns)
            start_idx = ns/ndata;
            figure(2)
            semilogy(t,abs(ene_snap.b(start_idx:start_idx:end)-ene.b)./ene_snap.b(start_idx:start_idx:end),...
                     'x',cr,cmap(ii,:),dispname,"$N="+nb+"$",ms,2); hold on
        else
            exit
        end
    end
end
figure(1);
ax=gca; ax.FontSize=5; xlim([0, T])
xlabel("$t$",intp,ltx,fs,6); ylabel("$TKE$",intp,ltx,fs,6);
%ylim([0.03 0.04]);
leg = legend({}, fs,5,intp,ltx,'location','best','NumColumns',2);
leg.ItemTokenSize = [12,18];
formatfig(ax); print(gcf,"intke_grom","-dpdf","-r300"); close(1)

figure(2);
ax=gca; ax.FontSize=5;
xlabel("$t$",intp,ltx,fs,6); ylabel("$TKE$",intp,ltx,fs,6);
xlim([0, T]); ylim([1e-5 1e2])
leg = legend({}, fs,5,intp,ltx,'location','best','NumColumns',2);
leg.ItemTokenSize = [12,18]
formatfig(ax); print(gcf,"intke_grom_err","-dpdf","-r300"); close(2)
