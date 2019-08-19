function Data = load_mat(subjectNum, domain)
    
    fname = sprintf('MDM_%s_%d.mat', domain, subjectNum);
    load(fname) % produces variable `Datamon` or 'Datamed' for convenience, change its name into 'Data'
    
    if strcmp(domain, 'MON') ==1
        Data = Datamon;
    else
        Data = Datamed;
    end