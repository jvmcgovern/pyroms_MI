
function nccopy_vars(varnamev,oldnc,newnc,opt)

%function nccopy_vars(varnamev,oldnc,newnc,opt)
%
%Copy NetCDF variables varnamev [1*nvar], including all attributes,
%from file oldnc to newnc.
%Variable data is written for variables with opt.write_data = 1 [1*nvar] (def = 1)
%
%%Phil Wallhead 23/03/2019
%%

if ~iscell(varnamev); varnamev = {varnamev}; end
nvarv = length(varnamev);
if nargin<4; opt = []; end
if isfield(opt,'verbose')==1; verbose = opt.verbose; else verbose = 0; end
if isfield(opt,'write_data')==1; write_data = opt.write_data; else write_data = ones(1,nvarv); end
if isfield(opt,'write_atts')==1; write_atts = opt.write_atts; else write_atts = ones(1,nvarv); end
if isfield(opt,'newvarnamev')==1; newvarnamev = opt.newvarnamev; else newvarnamev = []; end
if isfield(opt,'newdims')==1; newdims = opt.newdims; else newdims = []; end
if isfield(opt,'dimrename')==1; dimrename = opt.dimrename; else dimrename = []; end
if isfield(opt,'newdimname')==1; newdimname = opt.newdimname; else newdimname = []; end
if isfield(opt,'dimsubselname')==1; dimsubselname = opt.dimsubselname; else dimsubselname = []; end
if isfield(opt,'subsel')==1; subsel = opt.subsel; else subsel = []; end
if isfield(opt,'newfillvalue')==1; newfillvalue = opt.newfillvalue; else newfillvalue = NaN*ones(1,nvarv); end
if isfield(opt,'newmissing_value')==1; newmissing_value = opt.newmissing_value; else newmissing_value = NaN*ones(1,nvarv); end
if isfield(opt,'newdatatype')==1; newdatatype = opt.newdatatype; else newdatatype = []; end
if isfield(opt,'newformat')==1; newformat = opt.newformat; else newformat = []; end
if isfield(opt,'newunits')==1; newunits = opt.newunits; else newunits = []; end
if isfield(opt,'newcycle_length')==1; newcycle_length = opt.newcycle_length; else newcycle_length = NaN*ones(1,nvarv); end
if isfield(opt,'newlong_name')==1; newlong_name = opt.newlong_name; else newlong_name = []; end
if isfield(opt,'newfield')==1; newfield = opt.newfield; else newfield = []; end
if isfield(opt,'newfile')==1; newfile = opt.newfile; else newfile = 1; end

if isempty(newvarnamev)==1; newvarnamev = varnamev; end
if length(write_data)==1; write_data = write_data*ones(1,nvarv); end
if length(write_atts)==1; write_atts = write_atts*ones(1,nvarv); end
if length(newfillvalue)==1; newfillvalue = newfillvalue*ones(1,nvarv); end
if length(newmissing_value)==1; newmissing_value = newmissing_value*ones(1,nvarv); end
if length(newcycle_length)==1; newcycle_length = newcycle_length*ones(1,nvarv); end
if length(newfield)==1; newfield = newfield*ones(1,nvarv); end
samev = cell(1,nvarv); for i=1:nvarv; samev{i} = 'same'; end
if isempty(newdatatype)==1; newdatatype = samev; end
if isempty(newunits)==1; newunits = samev; end
if isempty(newlong_name)==1; newlong_name = samev; end
if isempty(newfield)==1; newfield = samev; end
if isempty(subsel)==0; nsubsel = length(subsel); end

if (newfile==0)
    finfonew = ncinfo(newnc);
    varNamesnew = {finfonew.Variables.Name};
else
    varNamesnew = [];
end

for j=1:nvarv   
    vinfo0 = ncinfo(oldnc,varnamev{j});
    
    %Gather attributes and fillvalue from old file
    if (isempty(vinfo0.Attributes)==0)
        attNames0 = {vinfo0.Attributes.Name};
        attValues0 = {vinfo0.Attributes.Value};
        
        if (isnan(newfillvalue(j)))
            fillvalue1 = [];
            for i=1:length(attNames0)
                if strcmpi(attNames0{i},'_FillValue')==1; fillvalue1 = attValues0{i}; end
            end
        else
            fillvalue1 = newfillvalue(j);
        end
    end
    
    %Gather dimensions from old file (or rename dimensions using newdims)
    if (isempty(newdims)==1)
        dims0 = vinfo0.Dimensions; ndims0 = length(dims0);
        dims1 = cell(1,2*ndims0);
        for i=1:ndims0
            dims11 = dims0(i).Name; 
            if strcmp(dims11,dimrename)==1; dims11 = newdimname; end
            dims1{(i-1)*2+1} = dims11; 
            if (vinfo0.Dimensions(i).Unlimited==1)
                dims1{i*2} = Inf;
            else
                if strcmp(dims11,dimsubselname)==1; dims1{i*2} = nsubsel; dimsubsel = i; else dims1{i*2} = dims0(i).Length; end
            end
        end
    else
        dims1 = newdims{j};
    end
    
    %Modify datatype of format here if req'd
    if (strcmpi(newdatatype{j},'same')==0)
        newdatatype1 = newdatatype{j};
    else
        newdatatype1 = vinfo0.Datatype;
    end
    if isempty(newformat)==1; newformat1 = vinfo0.Format; else newformat1 = newformat; end
    
    %Create the NetCDF variable
    if (~any(strcmp(varNamesnew,newvarnamev{j})))
        if (~isempty(fillvalue1) && strcmp(newformat1,'64bit')==0)
            nccreate(newnc,newvarnamev{j},'Dimensions',dims1,'DataType',newdatatype1,'FillValue',fillvalue1,'DeflateLevel',vinfo0.DeflateLevel,'Shuffle',vinfo0.Shuffle,'Format',newformat1)
        else
            nccreate(newnc,newvarnamev{j},'Dimensions',dims1,'DataType',newdatatype1,'DeflateLevel',vinfo0.DeflateLevel,'Shuffle',vinfo0.Shuffle,'Format',newformat1)
        end
        if verbose==1; disp(['Created netcdf variable ',newvarnamev{j}]); end
    end
    
    %Write the data if req'd
    if (write_data(j)==1)
        X1 = ncread(oldnc,varnamev{j});
        if (isempty(dimsubselname)==0)
            indstr = '(';
            for i=1:dimsubsel-1; indstr = [indstr,':,']; end %#ok<AGROW>
            indstr = [indstr,'subsel']; %#ok<AGROW>
            for i=1:ndims0-dimsubsel; indstr = [indstr,',:']; end %#ok<AGROW>
            indstr = [indstr,')']; %#ok<AGROW>
            eval(['X1c = X1',indstr,';'])
            ncwrite(newnc,newvarnamev{j},X1c);
        else
            ncwrite(newnc,newvarnamev{j},X1);
        end
        if verbose==1; disp(['Written data for netcdf variable ',newvarnamev{j}]); end
    end
    
    %Write the attributes if req'd
    if (isempty(vinfo0.Attributes)==0 && write_atts(j)==1)
        attNames1 = attNames0; attValues1 = attValues0;
        for i=1:length(attNames1)
            if (strcmpi(attNames1{i},'_FillValue')==0)
                ncwriteatt(newnc,newvarnamev{j},attNames1{i},attValues1{i});
            end
        end
        if verbose==1; disp(['Written attributes for netcdf variable ',newvarnamev{j}]); end
    end
    
    %Modify attributes as req'd
    if (strcmpi(newunits{j},'same')==0)
        fileattrib(newnc,'+w');
        ncwriteatt(newnc,newvarnamev{j},'units',newunits{j});
        if verbose==1; disp(['Modified units attribute for netcdf variable ',newvarnamev{j}]); end
    end
    if (~isnan(newcycle_length(j)))
        if (newcycle_length(j)==Inf) %Delete the attribute
            ncid = netcdf.open(newnc,'NC_WRITE');
            varid = netcdf.inqVarID(ncid,newvarnamev{j});
            netcdf.reDef(ncid);
            netcdf.delAtt(ncid,varid,'cycle_length');
            netcdf.close(ncid);
        else
            fileattrib(newnc,'+w');
            ncwriteatt(newnc,newvarnamev{j},'cycle_length',newcycle_length(j));
        end
        if verbose==1; disp(['Modified cycle_length attribute for netcdf variable ',newvarnamev{j}]); end
    end
    if (strcmpi(newlong_name{j},'same')==0)
        fileattrib(newnc,'+w');
        ncwriteatt(newnc,newvarnamev{j},'long_name',newlong_name{j});
        if verbose==1; disp(['Modified long_name attribute for netcdf variable ',newvarnamev{j}]); end
    end
    if (~isnan(newmissing_value(j)))
        fileattrib(newnc,'+w');
        ncwriteatt(newnc,newvarnamev{j},'missing_value',newmissing_value(j));
        if verbose==1; disp(['Modified missing_value attribute for netcdf variable ',newvarnamev{j}]); end
    end
    if (strcmpi(newfield{j},'same')==0)
        fileattrib(newnc,'+w');
        ncwriteatt(newnc,newvarnamev{j},'field',newfield{j});
        if verbose==1; disp(['Modified field attribute for netcdf variable ',newvarnamev{j}]); end
    end
end