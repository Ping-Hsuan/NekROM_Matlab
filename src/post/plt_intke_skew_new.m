clear all; close all;

hdr;

bu = dlmread("bu");
bu = reshape(bu,301,301);
snap = dlmread("uk");
snap = reshape(snap,301,2000);
msnap = mean(snap,2);
umax = dlmread("umax");
umin = dlmread("umin");

f=figure(1);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'Units', 'Inches', 'Position', [0, 0, fig_width, fig_height], 'PaperUnits', 'Inches', 'PaperSize', [fig_width, fig_height])
f=figure(2);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'Units', 'Inches', 'Position', [0, 0, fig_width, fig_height], 'PaperUnits', 'Inches', 'PaperSize', [fig_width, fig_height])
nb_list = [200];
cmap = colormap(lines);
for ii=1:size(nb_list,2)
    nb=nb_list(ii);
    ene_snap = table;
    for jj=1:size(snap,2)
        ene_snap = [ene_snap;{snap(1:nb+1,jj)'*bu(1:nb+1,1:nb+1)*snap(1:nb+1,jj)-2*snap(1:nb+1,jj)'*bu(1:nb+1,1:nb+1)*msnap(1:nb+1)+msnap(1:nb+1)'*bu(1:nb+1,1:nb+1)*msnap(1:nb+1)}];
%           ene_snap = [ene_snap;{snap(1:nb+1,jj)'*bu(1:nb+1,1:nb+1)*snap(1:nb+1,jj)}];
    end
    figure(1)
    t=linspace(25/size(ene_snap.b,1),25,size(ene_snap.b,1));
    semilogy(t,ene_snap.b,'k-',dispname,"Porjection, $N="+nb+"$"); hold on

    nb = nb_list(ii);
    grom = dlmread("g-rom_"+nb+"/ucoef");
    grom = reshape(grom,100,nb+1);
    ua = dlmread("g-rom_"+nb+"/ua");
    ene = table;
    for jj=1:size(grom,1)
        ene = [ene;{grom(jj,:)*bu(1:nb+1,1:nb+1)*grom(jj,:)'-2*grom(jj,:)*bu(1:nb+1,1:nb+1)*ua+ua'*bu(1:nb+1,1:nb+1)*ua}];
    end
    figure(1)
    t=linspace(25/size(ene.b,1),25,size(ene.b,1));
    semilogy(t,ene.b,'-.',cr,cmap(ii,:),dispname,"G-ROM, $N="+nb+"$"); hold on
    figure(2)
    semilogy(t,abs(ene_snap.b(20:20:end)-ene.b)./ene_snap.b(20:20:end),'x',cr,cmap(ii,:),dispname,"$N="+nb+"$",ms,2); hold on
    cpr_list = [10 100 200 400];
    for kk=1:size(cpr_list,2)
        cpr = cpr_list(kk);
        cpd = dlmread("cpd_skew_"+nb+"_rank"+cpr+"/ucoef");
        cpd = reshape(cpd,100,nb+1);
        ua = dlmread("cpd_skew_"+nb+"_rank"+cpr+"/ua");
        ene_cpd = table;
        for jj=1:size(grom,1)
            ene_cpd = [ene_cpd;{cpd(jj,:)*bu(1:nb+1,1:nb+1)*cpd(jj,:)'-2*cpd(jj,:)*bu(1:nb+1,1:nb+1)*ua+ua'*bu(1:nb+1,1:nb+1)*ua}];
        end
        t=linspace(25/size(ene_cpd.b,1),25,size(ene_cpd.b,1));

        figure(1)
        semilogy(t,ene_cpd.b,'-.',cr,cmap(kk+1,:),dispname,"CPD-ROM, R="+cpr); hold on
        figure(2)
        semilogy(t,abs(ene_snap.b(20:20:end)-ene_cpd.b)./ene_snap.b(20:20:end),'x',cr,cmap(kk+1,:),dispname,"CPD-ROM, R="+cpr,ms,2); hold on
    end
end
figure(1)
ax=gca; ax.FontSize=5;
xlabel("$t$",intp,ltx,fs,6);
ylabel("$TKE$",intp,ltx,fs,6);
xlim([0, 25])
%ylim([0.03 0.04]);
leg = legend({}, fs,5,intp,ltx,'location','best','NumColumns',2);
leg.ItemTokenSize = [10,18]
formatfig(ax)
print(gcf,"intke_skew","-dpdf","-r300")
close(1)

figure(2)
ax=gca; ax.FontSize=5;
xlabel("$t$",intp,ltx,fs,6);
ylabel("$TKE$",intp,ltx,fs,6);
xlim([0, 25])
ylim([1e-5 1e2])
leg = legend({}, fs,5,intp,ltx,'location','best','NumColumns',2);
leg.ItemTokenSize = [10,18]
formatfig(ax)
print(gcf,"intke_skew_err","-dpdf","-r300")
close(2)
