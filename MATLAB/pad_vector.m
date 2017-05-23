function [search_region_paddded_vec] = pad_vector(template_vec,search_region_vec)
%pads the edge of the matrix with zeros
%Takes the input of 
%the 'template' is the smaller vector which is passed over the top of the 
%'search region'
te_len=length(template_vec);
sr_len=length(search_region_vec);

padded_size=(sr_len-1)*2+te_len
padded_matrix=zeros(padded_size,1)

% %determine the extent of the image inside the pad
top_extent=te_len;
bottom_extent=padded_size-(te_len-1);
padded_matrix(top_extent:bottom_extent)=search_region_vec

end


