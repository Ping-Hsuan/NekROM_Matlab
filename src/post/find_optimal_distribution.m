clear all; close all;


reprods = {
    readtable("./reprod/reprod_leray.csv"),
    readtable("./reprod/reprod_efr.csv"),
    readtable("./reprod/reprod_tr.csv")
};
short_preds = {
    readtable("./short_pred/spred_leray.csv"),
    readtable("./short_pred/spred_efr.csv"),
    readtable("./short_pred/spred_tr.csv")
};

medium_preds = {
    readtable("./medium_pred/mpred_leray.csv"),
    readtable("./medium_pred/mpred_efr.csv"),
    readtable("./medium_pred/mpred_tr.csv")
};
long_preds = {
    readtable("./long_pred/lpred_leray.csv"),
    readtable("./long_pred/lpred_efr.csv"),
    readtable("./long_pred/lpred_tr.csv")
};


markers = {'o', 's', 'd', '^', 'v', '>', '<', 'p', 'h', 'x'};

for i=1:3
        figure(i)
        T = reprods{i};
        legend_labels = cell(size(T.nb));
        for j=1:numel(T.nb)
            scatter(T.radius(j), T.relax(j), [], 'k', markers{j}); hold on 
            legend_labels{j} = sprintf('nb=%d', T.nb(j));
        end
        T = short_preds{i};
        for j=1:numel(T.nb)
            scatter(T.radius(j), T.relax(j), [], 'b', markers{j}); hold on 
        end
        T = medium_preds{i};
        for j=1:numel(T.nb)
            scatter(T.radius(j), T.relax(j), [], 'r', markers{j}); hold on 
        end
        T = long_preds{i};
        for j=1:numel(T.nb)
            scatter(T.radius(j), T.relax(j), [], 'g', markers{j}); hold on 
        end
        xlim([1e-3 1e-1]);
        ylim([0.01, 10.24]);
        legend(legend_labels, 'Location', 'best');
        set(gca,'xscale','log','yscale','log');
end
