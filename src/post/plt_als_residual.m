clear all; close all;

hdr;

cmap = colormap(lines);

%% Load ALS residual
fileList_als = dir('als_skew*.txt');
fileList_quad = dir('als_quad*.txt');

f=figure(1);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'Units', 'Inches', 'Position', [0, 0, fig_width, fig_height],...
    'PaperUnits', 'Inches', 'PaperSize', [fig_width, fig_height]);

N_list = [100];
R_lists = {[10,100,200,400]}
for jj=1:size(N_list,2)
     N = N_list(jj)
     R_list = R_lists{jj};
     T1 = table;
     T2 = table;
     for ii=1:size(R_list,2)
         R=R_list(ii)
         als = readtable("./als_skew_R"+R+"_N"+N);
         als_core = readtable("./als_quad_lamb0p7_R"+R+"_N"+N);
     %%------ Done Reading Tables --------
         tmp = find(als.residual);
         idx1 = tmp(end)
         tmp = find(als_core.residual);
         idx2 = tmp(end)
         T1 = [T1;{R,als.residual(idx1)}]
         T2 = [T2;{R,als_core.residual(idx2)}]
     end
     figure(1)
     semilogy(T2.b1,T2.b2,'-o',cr,cmap(jj,:),ms,2,dispname,"ALS-quad"); hold on
     semilogy(T1.b1,T1.b2,'--x',cr,cmap(jj,:),ms,2,dispname,"ALS"); hold on
 end

figure(1)
ax=gca;
ax.FontSize=5;
xlabel('$R$',intp,ltx,fs,6);
ylabel('Terminated relative residual $\|\mathcal{C}_u-\widetilde{\mathcal{C}}_u\|/\    |\mathcal{C}_u\|$',intp,ltx,fs,6);
grid on;
leg = legend({}, fs,5,intp,ltx,'location','best','NumColumns',2);
leg.ItemTokenSize = [10,18]
formatfig(ax)
print(gcf,"final_residual","-dpng","-r300")
close(1)

f=figure(1);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'Units', 'Inches', 'Position', [0, 0, fig_width, fig_height],...
    'PaperUnits', 'Inches', 'PaperSize', [fig_width, fig_height]);

N_list = [100];
R_lists = {[10,100,200,400]}
for jj=1:size(N_list,2)
     N = N_list(jj)
     R_list = R_lists{jj};
     T1 = table;
     T2 = table;
     for ii=1:size(R_list,2)
         R=R_list(ii)
         als = readtable("./als_skew_R"+R+"_N"+N);
         als_core = readtable("./als_quad_lamb0p7_R"+R+"_N"+N);
     %%------ Done Reading Tables --------
         tmp = find(als.perr);
         idx1 = tmp(end)
         tmp = find(als_core.perr);
         idx2 = tmp(end)
         T1 = [T1;{R,als.perr(idx1)}]
         T2 = [T2;{R,als_core.perr(idx2)}]
     end
     figure(1)
     semilogy(T2.b1,T2.b2,'-o',cr,cmap(jj,:),ms,2,dispname,"ALS-quad"); hold on
     semilogy(T1.b1,T1.b2,'--x',cr,cmap(jj,:),ms,2,dispname,"ALS"); hold on
 end

figure(1)
ax=gca;
ax.FontSize=5;
xlabel('$R$',intp,ltx,fs,6);
ylabel('Terminated relative residual $\|\mathcal{C}_u-\widetilde{\mathcal{C}}_u\|/\    |\mathcal{C}_u\|$',intp,ltx,fs,6);
grid on;
leg = legend({}, fs,5,intp,ltx,'location','best','NumColumns',2);
leg.ItemTokenSize = [10,18]
formatfig(ax)
print(gcf,"final_maxcerr","-dpng","-r300")
