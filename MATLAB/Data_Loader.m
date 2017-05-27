% Mulitpile Data file import script.
%Hint: If your data has the format: my_data0001.bin, then use the data
%number like this: data_number=num2str(1,'%04i')
%Otherwise just use as normal with 'i'

%To Adapt this to differnt data, you should be able to just change the
%first paragraph

%% Load Hot Film Data
%Assume 1d data read in

%Change these:
hf_data_points=40;              %number of different files to read with these prefixes
folder_string=['Data/'];        %folder location relative to current folder
hf_string=['u_hf_ypos'];        %name string before ID number
bin_string=['.bin'];            %file type / name after ID
data_storage_type=['*float']; %how the data will be read in

%read one image in to find it's size
data_number=num2str(1,'%i'); 
data_loc = strcat(folder_string,hf_string,data_number,bin_string);
fid = fopen(data_loc, 'r');
hf_data1 = fread(fid, data_storage_type);
hf_length=length(hf_data1);

%read the rest of the data in
hf_matrix=zeros(hf_length,hf_data_points)-1;
hf_matrix(:,1)=hf_data1;

for k=2:hf_data_points;
    data_number=num2str(k,'%i'); 
data_loc = strcat(folder_string,hf_string,data_number,bin_string);
fid = fopen(data_loc, 'r');
hf_matrix(:,k) = fread(fid, data_storage_type);
end

%% Load Hot Wire Data
%Assume 1d data read in

hw_data_points=40;              %number of different files to read with these prefixes
folder_string=['Data/'];        %folder location relative to current folder
hf_string=['u_hw_ypos'];        %name string before ID number
bin_string=['.bin'];            %file type / name after ID

%read one image in to find it's size
data_number=num2str(1,'%i'); 
data_loc = strcat(folder_string,hf_string,data_number,bin_string);
fid = fopen(data_loc, 'r');
hw_data1 = fread(fid, data_storage_type);
hw_length=length(hw_data1);

%read the rest of the data in
hw_matrix=zeros(hw_length,hw_data_points)-1;
hw_matrix(:,1)=hw_data1;

for k=2:hw_data_points;
    data_number=num2str(k,'%i'); 
data_loc = strcat(folder_string,hf_string,data_number,bin_string);
fid = fopen(data_loc, 'r');
hw_matrix(:,k) = fread(fid, data_storage_type);
end

%% Load Hw Burst Signals In

burst_hw_data_points=10;              %number of different files to read with these prefixes
folder_string=['Data/'];        %folder location relative to current folder
hf_string=['u_hw_ypos20_burst'];        %name string before ID number
bin_string=['.bin'];            %file type / name after ID

%read one image in to find it's size
data_number=num2str(1,'%i'); 
data_loc = strcat(folder_string,hf_string,data_number,bin_string);
fid = fopen(data_loc, 'r');
burst_hw_data1 = fread(fid, data_storage_type);
burst_hw_length=length(burst_hw_data1);

%read the rest of the data in
burst_hw_matrix=zeros(burst_hw_length,burst_hw_data_points)-1;
burst_hw_matrix(:,1)=burst_hw_data1;

for k=2:burst_hw_data_points;
    data_number=num2str(k,'%i'); 
data_loc = strcat(folder_string,hf_string,data_number,bin_string);
fid = fopen(data_loc, 'r');
burst_hw_matrix(:,k) = fread(fid, data_storage_type);
end

%% Load Hf Burst Signals In

burst_hf_data_points=10;              %number of different files to read with these prefixes
folder_string=['Data/'];        %folder location relative to current folder
hf_string=['u_hf_ypos20_burst'];        %name string before ID number
bin_string=['.bin'];            %file type / name after ID

%read one image in to find it's size
data_number=num2str(1,'%i'); 
data_loc = strcat(folder_string,hf_string,data_number,bin_string);
fid = fopen(data_loc, 'r');
burst_hf_data1 = fread(fid, data_storage_type);
burst_hf_length=length(burst_hf_data1);

%read the rest of the data in
burst_hf_matrix=zeros(burst_hf_length,burst_hf_data_points)-1;
burst_hf_matrix(:,1)=burst_hf_data1;

for k=2:burst_hf_data_points;
    data_number=num2str(k,'%i'); 
data_loc = strcat(folder_string,hf_string,data_number,bin_string);
fid = fopen(data_loc, 'r');
burst_hf_matrix(:,k) = fread(fid, data_storage_type);
end


% save(['hw_hf_data.mat'],'burst_hw_matrix','burst_hf_matrix','hf_matrix','hw_matrix')
% too large, can break up or can just use data_loader.m
%git reset --hard origin/master


