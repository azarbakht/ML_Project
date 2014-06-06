function Data= dsread(filename)

%fid1 = fopen('../fuzzer_data/d.fcov.ds');
fid1 = fopen(filename);

tline = fgetl(fid1);%reading the header of the features
Data = [];

tline = fgetl(fid1);

while ischar(tline)
    
  %  disp(tline);
%    tline = fgetl(fid1);
    oneValue = [];
    dataOneLine = [];
    
    for i=1: size(tline,2)        
        if(tline(i) == ' ')
            dataOneLine = [dataOneLine str2double(oneValue)];
            oneValue=[];
        else
            oneValue=[oneValue tline(i)];
        end        
    end
    if size(oneValue,2)>0
       Data = [Data ; dataOneLine str2double(oneValue) ];
    else
        Data = [Data ; dataOneLine ];
    end
    tline = fgetl(fid1);
end% end while
fclose(fid1);
end % end fucntion