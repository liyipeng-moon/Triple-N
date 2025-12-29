function cmap = colormap_matplotlib(cmapname,n,quiet)
%generate colormaps used in matplotlib. For list of available colormaps, call colormap_matplotlib('list')

if(nargin < 1 || isempty(cmapname))
    cmapname='list';
end
        
if(nargin < 2 || isempty(n))
    n=256;
end

if(nargin < 3 || isempty(quiet))
    quiet=false;
end
if(strcmpi(quiet(1),'q') || strcmpi(n(1),'q'))
    quiet=true;
end

%look for saved colormap file
[mdir,~]=fileparts(mfilename('fullpath'));
cmap_matfile=sprintf('%s/matplotlib_colormaps.mat',mdir);
is_newfile=false;
if(~exist(cmap_matfile,'file'))
    fprintf('%s not found. Trying to generate file via python script\n',cmap_matfile);
    system(sprintf('python %s/save_matplotlib_colormaps.py %s',mdir,cmap_matfile));
    if(~exist(cmap_matfile,'file'))
        error('External python call failed (might not have matplotlib available?). Try the following command:\npython %s/save_matplotlib_colormaps.py %s',...
            mdir,cmap_matfile);
    end
    is_newfile=true;
end

C=load(cmap_matfile);
cmapnames=fieldnames(C);
cmapnames=cmapnames(cellfun(@isempty,regexp(cmapnames,'^version_'))); %python script inserts version info but ignore these

if(~quiet && (is_newfile || strcmpi(cmapname,'list')))
    fprintf('matplotlib colormaps: python version=%s, matplotlib version=%s (disable this message with "quiet" argument)\n',C.version_python,C.version_matplotlib);
end

if(strcmpi(cmapname,'list'))
    cmap=cmapnames;
    return;
end


cmidx=find(strcmpi(cmapnames,cmapname));
if(isempty(cmidx))
    error('Colormap not found: %s\n',cmapname);
end

cmap=C.(cmapnames{cmidx});
if size(cmap,1)~=n
    cmap=interp1(linspace(0,1,size(cmap,1)),cmap,linspace(0,1,n));
end
