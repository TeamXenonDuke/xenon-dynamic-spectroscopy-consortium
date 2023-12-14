function twix_obj = readtwix(twix_path)

if exist('twix_path') 
    twix_obj = mapVBVD(twix_path);
else 
    twix_obj = mapVBVD('D:\Elly\Documents\Duke\CIVM Research\Siemens Data\');
end 

twix_obj.data =  squeeze(double(twix_obj.image()));

end 