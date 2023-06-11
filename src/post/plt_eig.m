close all; clear all;

hdr;

get_eig

%eig = dlmread("eig.txt");

f=figure(1);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'Units', 'Inches', 'Position', [0, 0, fig_width, fig_height], 'PaperUnits', 'Inches', 'PaperSize', [fig_width, fig_height])
loglog(eig,lw,1); hold on
ax=gca;
ax.FontSize=5;
%legend({"2D Cylinder","2D Lid-Driven","3D Lid-Driven","3D Rayleigh","3D Turbulent Channel"},fs,5,intp,ltx,'Location','Best');
ylabel("$\lambda_i$",intp,ltx,fs,6);
xlabel("$i$",intp,ltx,fs,6);
grid on
%title({"Comparison of the behavior of POD eigenvalues with $N$","between different cases"},intp,ltx,fs,10);
print(gcf,"eig","-dpdf","-r300")

f=figure(2);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'Units', 'Inches', 'Position', [0, 0, fig_width, fig_height], 'PaperUnits', 'Inches', 'PaperSize', [fig_width, fig_height])
loglog(eig/sum(eig),'-',lw,1); hold on
ax=gca;
ax.FontSize=5;
%legend({"2D Cylinder","2D Lid-Driven","3D Lid-Driven","3D Rayleigh","3D Turbulent Channel"},fs,5,intp,ltx,'Location','Best');
ylabel("$\lambda_i/\sum^K_j \lambda_j$",intp,ltx,fs,6);
xlabel("$i$",intp,ltx,fs,6);
grid on
%title({"Comparison of the behavior of normalized POD eigenvalues","with $N$ between different cases"},intp,ltx,fs,10);
print(gcf,"norm_eig","-dpdf","-r300")

f=figure(3);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'Units', 'Inches', 'Position', [0, 0, fig_width, fig_height], 'PaperUnits', 'Inches', 'PaperSize', [fig_width, fig_height])
loglog(eig/(eig(1)),'-',lw,1); hold on
ax=gca;
ax.FontSize=5;
%legend({"2D Cyl","2D LDC","3D LDC","3D Ray","3D TCH"},fs,5,intp,ltx,'Location','Best');
ylabel("$\lambda_i/\lambda_1$",intp,ltx,fs,6);
xlabel("$i$",intp,ltx,fs,6);
grid on
%title({"The behavior of normalized POD eigenvalues"},intp,ltx,fs,6);
print(gcf,"normweig1_eig","-dpdf","-r300")

f=figure(4);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'Units', 'Inches', 'Position', [0, 0, fig_width, fig_height], 'PaperUnits', 'Inches', 'PaperSize', [fig_width, fig_height])
semilogx(cumsum(eig)/sum(eig),'-',lw,1); hold on
semilogx(0.9*ones(1e4,1),'k-',lw,1);
ax=gca;
ax.FontSize=6;
%legend({"2D Cyl, $\rm Re=100$","2D Ldc, $\rm Re=15000$","3D Ldc, $\rm Re=3200$","3D Ray, $\rm Re=394$","3D Tch, $\rm Re=2812$"},fs,8,intp,ltx,'Location','Best');
leg.ItemTokenSize = [10,18]
ylabel("$\sum^i_{j=1}\lambda_j/\sum^K_{j=1}\lambda_j$",intp,ltx,fs,8);
xlabel("$i$",intp,ltx,fs,8);
grid on
%formatfig(4)
print(gcf,"acceig","-dpdf","-r300")
%saveas(gcf,"acceig_compare.png");
