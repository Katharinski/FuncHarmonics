function get_medial_wall_indices()

    homedir = getenv('HOME');
    repodatadir = fullfile(homedir,'repos','FuncHarmonics','data');

    wb_dir = fullfile(homedir,'snap','workbench/bin_linux64/wb_command');
    indfile = fullfile(repodatadir,'HCP','Human.MedialWall_Conte69.32k_fs_LR.dlabel.nii');
    indData = ciftiopen(indfile,wb_dir);
    indices = logical(indData.cdata); % indices contains medial wall locations
    save(fullfile(repodatadir,'HCP_derived','Ind_S900'),'indices')
end