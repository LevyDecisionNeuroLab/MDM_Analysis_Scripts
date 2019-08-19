function save_mat(Data, subjectNum, domain, fitpar_out_path)

    if strcmp(domain, 'MON') ==1
        Datamon = Data;
        save(fullfile(fitpar_out_path, ['MDM_' domain '_' num2str(subjectNum) '_fitpar.mat']), 'Datamon')
    else
        Datamed = Data;
        save(fullfile(fitpar_out_path, ['MDM_' domain '_' num2str(subjectNum) '_fitpar.mat']), 'Datamed')
    end