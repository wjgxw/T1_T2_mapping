function Para=getpar(pathname,varargin)
%%
% paracell=[at sw np ...]

%%

if strcmpi(pathname(end-5:end),'acqfil') && ~strcmpi(pathname(end-3:end),'.fid')
    pathname=[pathname(1:end-6),'recon'];
end

pathname=strcat(pathname,'/');


parfilefp=fopen(strcat(pathname,'procpar'),'rt');
if parfilefp==-1
	parfilefp=fopen(strcat(pathname,'procpar'),'rt');
	if parfilefp==-1
		tempArr=strfind(pathname,'/');
		parfile=strcat(pathname(1,1:tempArr(length(tempArr)-1)),'procpar');
		parfilefp=fopen(parfile,'rt');
		if parfilefp==-1
			disp('Error:Can''t find parameter file "procpar" in data director.');
			return;
		end
	end
end

N=length(varargin);
for k=1:N
    
    %retrieve varargin{k}  parameter
    bfound=0;
    while(feof(parfilefp)==0)
    	tempArr=fgetl(parfilefp);
        lengthk=length(varargin{k});
    	if(strncmp(tempArr,[varargin{k},' '],lengthk+1)==1) 
    		tempArr=fgetl(parfilefp);
            tempArr=regexp(tempArr,'\s+','split');
            if length(tempArr)==2 && tempArr{2}(1)=='"'
                Para.(varargin{k})=tempArr{2}(2:end-1);
            else
                Para.(varargin{k})=str2double(tempArr(2:end-1));
            end
            bfound=1;
    		break;
    	end
    end
    if(bfound==0)
	disp([varargin{k},' can not be found']);
    end
    frewind(parfilefp);
end
fclose(parfilefp);

end