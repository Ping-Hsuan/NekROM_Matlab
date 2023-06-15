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

T = 100;

f=figure(1);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'Units', 'Inches', 'Position', [0, 0, fig_width, fig_height],...
    'PaperUnits', 'Inches', 'PaperSize', [fig_width, fig_height])

nb_list = [300];
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
    semilogy(t,diag(recon_err_L2)./diag(snapnorm_L2),'-',cr,cmap(ii,:),dispname,"G-ROM, $N="+nb+"$",ms,2); hold on

    cpr_list = [10 100 200 400 800];
    for kk=1:size(cpr_list,2)
        cpr = cpr_list(kk);
        cpd = dlmread("cpd_quad_"+nb+"_rank"+cpr+"/ucoef");
        cpd = reshape(cpd,ndata,nb+1);
        if (ndata == ns)
            err  = snap(1:nb+1,:)-cpd';
            snapnorm_L2 = sqrt(snap(1:nb+1,:)'*bu(1:nb+1,1:nb+1)*snap(1:nb+1,:));
        elseif (ndata ~= ns)
            start_idx = ns/ndata;
            err  = snap(1:nb+1,start_idx:start_idx:end)-cpd';
            snapnorm_L2 = sqrt(snap(1:nb+1,start_idx:start_idx:end)'*bu(1:nb+1,1:nb+1)*snap(1:nb+1,start_idx:start_idx:end));
        end
        recon_err_L2 = sqrt(err'*bu(1:nb+1,1:nb+1)*err);

        t=linspace(T/size(err,2),T,size(err,2));
        semilogy(t,diag(recon_err_L2)./diag(snapnorm_L2),':',cr,cmap(kk+1,:),dispname,"ALS-quad, R="+cpr,ms,2); hold on
    end
    for kk=1:size(cpr_list,2)
        cpr = cpr_list(kk);
        cpd = dlmread("cpd_skew_"+nb+"_rank"+cpr+"/ucoef");
        cpd = reshape(cpd,ndata,nb+1);
        if (ndata == ns)
            err  = snap(1:nb+1,:)-cpd';
            snapnorm_L2 = sqrt(snap(1:nb+1,:)'*bu(1:nb+1,1:nb+1)*snap(1:nb+1,:));
        elseif (ndata ~= ns)
            start_idx = ns/ndata;
            err  = snap(1:nb+1,start_idx:start_idx:end)-cpd';
            snapnorm_L2 = sqrt(snap(1:nb+1,start_idx:start_idx:end)'*bu(1:nb+1,1:nb+1)*snap(1:nb+1,start_idx:start_idx:end));
        end
        recon_err_L2 = sqrt(err'*bu(1:nb+1,1:nb+1)*err);
        t=linspace(T/size(err,2),T,size(err,2));
        semilogy(t,diag(recon_err_L2)./diag(snapnorm_L2),'-.',cr,cmap(kk+1,:),dispname,"ALS, R="+cpr,ms,2); hold on
    end
end
ax=gca; ax.FontSize=5;
xlabel("$t$",intp,ltx,fs,4);
ylabel("$\|\Pi u(t) - u_{ROM}(t)\|_{L^2}\|/\|\Pi u(t)\|_{L^2}$",intp,ltx,fs,4);
xlim([0, T])
leg = legend({}, fs,4,intp,ltx,'location','best','NumColumns',3);
leg.ItemTokenSize = [16,18]
formatfig(ax)
print(gcf,"recon_err","-dpdf","-r300")
