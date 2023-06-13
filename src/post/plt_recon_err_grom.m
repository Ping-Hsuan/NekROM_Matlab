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

%% Setting up flags
fileList = dir('g-rom*');
if isempty(fileList)
    ifgrom = false;
else
    ifgrom = true;
end
fileList = dir('cpd_skew*');
if isempty(fileList)
    ifskew = false;
else
    ifskew = true;
end
fileList = dir('cpd_quad*');
if isempty(fileList)
    ifquad = false;
else
    ifquad = true;
end

T = 25;

f=figure(1);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'Units', 'Inches', 'Position', [0, 0, fig_width, fig_height],...
    'PaperUnits', 'Inches', 'PaperSize', [fig_width, fig_height])

nb_list = [10, 20, 40, 80, 100, 200, 300];
cmap = colormap(lines);
for ii=1:size(nb_list,2)
    nb = nb_list(ii);

    grom = dlmread("g-rom_"+nb+"/ucoef");
    ndata = length(grom)/(nb+1);
    grom = reshape(grom,ndata,nb+1);
    if (ndata == ns)
        err  = snap(1:nb+1,:)-grom';
        snapnorm_L2 = sqrt(snap(1:nb+1,:)'*bu(1:nb+1,1:nb+1)*snap(1:nb+1,:));
    elseif (ndata ~= ns)
        start_idx = ns/ndata;
        err  = snap(1:nb+1,start_idx:start_idx:end)-grom';
        snapnorm_L2 = sqrt(snap(1:nb+1,start_idx:start_idx:end)'*bu(1:nb+1,1:nb+1)*snap(1:nb+1,start_idx:start_idx:end));
    end
    recon_err_L2 = sqrt(err'*bu(1:nb+1,1:nb+1)*err);

    t = linspace(T/ndata,T,ndata);
    semilogy(t,diag(recon_err_L2)./diag(snapnorm_L2),'-',cr,cmap(ii,:),dispname,"$N="+nb+"$",ms,2); hold on
end
ax=gca; ax.FontSize=5;
xlabel("$t$",intp,ltx,fs,6);
ylabel("$\|\Pi u_{snapshot}(t) - u_{ROM}(t)\|_{L^2}\|/\|\Pi u_{snapshot}(t)\|_{L^2}$",intp,ltx,fs,6);
xlim([0, T])
leg = legend({}, fs,4,intp,ltx,'location','best','NumColumns',3);
leg.ItemTokenSize = [16,18]
formatfig(ax)
print(gcf,"recon_err_grom","-dpdf","-r300")
