function Convert_sdf_to_xyz(filename)

% Need Open Babel to convert sdf file to xyz format

    [p,f]=fileparts(filename);
    filenameout=fullfile(p,f);
    disp(filename);
    command1 = append('cd ', pwd);
    system(command1);
    command2 = append('obabel ',filename,' -O ', filenameout,'.xyz');
    system(command2);


end