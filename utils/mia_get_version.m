function version = mia_get_version()
%MIA_GET_VERSION Return the current version of MIA
version = 'github-master';
major=0;
minor=0;
patches=0;

path = fileparts(which('mia'));

id = fopen(fullfile(fileparts(path), 'VERSION'));
str=fread(id,'*char' )';

tmp = strsplit(str,'_') ; 

version = tmp{2};

end

