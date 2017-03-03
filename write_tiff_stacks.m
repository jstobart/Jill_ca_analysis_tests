function fName = write_tiff_stacks(CellScan, filename)
    rawdata = CellScan.rawImg.rawdata;
    ch1 = rawdata(:,:,1,:);
    ch2 = rawdata(:,:,2,:);

    fName{1} = writeTiffStack([filename, '_ch1.tif'], squeeze(ch1));
    fName{2} = writeTiffStack([filename, '_ch2.tif'], squeeze(ch2));
end

function filename = writeTiffStack(filename, stack)
    % writes Tiffstack that can be opened in ImageJ
    if isempty(strfind(filename, '.tif'))
        filename = strcat(filename, '.tif');
    end

    if length(size(stack)) ~= 3
        error('Function only works for 3D stacks')
    end

    stack = uint16(stack);
    for iFrame=1:size(stack,3)
        imwrite(stack(:, :, iFrame), filename, 'WriteMode', 'append', 'Compression','none');
    end
end