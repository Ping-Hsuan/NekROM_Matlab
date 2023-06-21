function TOut = generate_table(directory, nb_list, save_file, mini_input)
    TOut = table;
    nb_list_all = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100];
    for ii=1:size(nb_list,2)
        nb = nb_list(ii);
        idx = find(nb_list_all==nb)
        T = readtable(directory+"/table_m1_"+directory+"_N"+nb);
        if nargin == 4
            mini = mini_input(idx);
        else
            [minv, mini] = min(T.uu_err);
        end
        tt = {nb, mini, T.relax(mini), T.radius(mini), T.uavg_err(mini), T.uu_err(mini), T.vv_err(mini), T.ww_err(mini)};
        TOut = [TOut; tt];
    end
    TOut.Properties.VariableNames = {'nb','index','relax','radius','uavg_err','uu_err','vv_err','ww_err'}
    TOut.Properties.Description = datestr(now, 'yyyy-mm-dd');
    writetable(TOut, save_file);
end
