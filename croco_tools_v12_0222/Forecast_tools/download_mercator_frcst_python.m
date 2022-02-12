function mercator_name=download_mercator_frcst_python(pathMotu,user,password,mercator_type, ...
                                                      motu_url, service_id,product_id, ...
                                                      lh,lf, ...
                                                      lonmin,lonmax,latmin,latmax,zmax, ...
                                                      FRCST_dir,FRCST_prefix,raw_mercator_name,Yorig)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Extract a subgrid from mercator to get a CROCO forcing
% Store that into monthly files.
% Take care of the Greenwitch Meridian.
% 
%
%  Further Information:  
%  http://www.croco-ocean.org
%  
%  This file is part of CROCOTOOLS
%
%  CROCOTOOLS is free software; you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published
%  by the Free Software Foundation; either version 2 of the License,
%  or (at your option) any later version.
%
%  CROCOTOOLS is distributed in the hope that it will be useful, but
%  WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program; if not, write to the Free Software
%  Foundation, Inc., 59 Temple Place, Suite 330, Boston,
%  MA  02111-1307  USA
%
%  Copyright (c) 2005-2006 by Pierrick Penven 
%  e-mail:Pierrick.Penven@ird.fr  
%
%  Updated    6-Sep-2006 by Pierrick Penven
%  Updated    20-Aug-2008 by Matthieu Caillaud & P. Marchesiello
%  Update     23-Oct-2020 by Gildas Cambon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pathMotu is a deprecated parameter ! 
download_raw_data=1;
convert_raw2crocotools=1; % convert -> crocotools format data
%
% Set variable names according to mercator type data
%
if mercator_type==1,
  vars = {'zos' ...
	  'uo' ...
	  'vo' ...
	  'thetao' ...
	  'so'};
else
  vars = {'zos' ...
	  'uo' ...
	  'vo' ...
	  'thetao' ...
	  'so'};
end
%
% Get dates
%
rundate_str=date;
rundate=datenum(rundate_str)-datenum(Yorig,1,1);

for i=1:lh+1
    time1(i)=datenum(rundate_str)-(lh+2-i);
end
time2=datenum(rundate_str);
for j=1:lf+1
    time3(j)=datenum(rundate_str)+j;
end
time=cat(2,time1,time2,time3);
tiempo_inicial = time1(1);
tiempo_final = time3(end);
if (lonmin > 180)
    lonmin = lonmin - 360;
end
if (lonmax > 180)
    lonmax = lonmax - 360;
end
disp([' '])
disp(['Get data for ',rundate_str])
disp(['Minimum Longitude: ',num2str(lonmin)])
disp(['Maximum Longitude: ',num2str(lonmax)])
disp(['Minimum Latitude:  ',num2str(latmin)])
disp(['Maximum Latitude:  ',num2str(latmax)])
disp([' '])

if download_raw_data
    %
    % Get data 
    % (temporarily removing Matlab lib path to avoid conflict with Python, mandatory with python 2.7.X)
    %  for example problem like : 
    %   Execution failed: /usr/lib/python2.7/lib-dynload/pyexpat.x86_64-linux-gnu.so: 
    %                     undefined symbol: XML_SetHashSalt" ) 
    
    pathdyld=getenv('DYLD_LIBRARY_PATH');
    setenv('DYLD_LIBRARY_PATH','');
    pathld=getenv('LD_LIBRARY_PATH');
    setenv('LD_LIBRARY_PATH','');
    
    get_file_python_mercator(pathMotu,mercator_type, ...
                             motu_url,service_id,product_id, ...
                             vars, ...
                             [lonmin-1 lonmax+1 latmin-1 latmax+1 0 zmax], ...
                             {datestr(tiempo_inicial,'yyyy-mm-dd') ...
                             datestr(tiempo_final,  'yyyy-mm-dd')}, ...
                             {user password}, ...
                             raw_mercator_name);
    setenv('DYLD_LIBRARY_PATH',pathdyld); % set back Matlab libs path
    setenv('LD_LIBRARY_PATH',pathld);     % set back Matlab libs path
end  %end download_raw_data


if convert_raw2crocotools 
    %
    % Convert data format and write in a more CROCOTOOLS 
    % compatible input file 
    %
    disp(['Making output data directory ',FRCST_dir]) % create directory
    eval(['!mkdir ',FRCST_dir])
    %
    mercator_name=[FRCST_dir,FRCST_prefix,num2str(rundate),'.cdf'];
    if exist(mercator_name)  
        disp('Mercator file already exist => overwrite it')
    end
    write_mercator_frcst(FRCST_dir,FRCST_prefix,raw_mercator_name, ...
                         mercator_type,vars,time,Yorig); % write data
    eval(['! rm -f ',raw_mercator_name]);
end

end

