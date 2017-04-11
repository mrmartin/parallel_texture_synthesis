function export_terragen(file, terrain)
    f=fopen(file,'w');
    fprintf(f,'TERRAGENTERRAIN SIZE');
    fwrite(f,min(size(terrain,1),size(terrain,2))-1,'int16');
    fwrite(f,0,'int16');
    fprintf(f,'XPTS');
    fwrite(f,size(terrain,1),'int16');
    fwrite(f,0,'int16');
    fprintf(f,'YPTS');
    fwrite(f,size(terrain,2),'int16');
    fwrite(f,0,'int16');
    fprintf(f,'SCAL');
    fwrite(f,0,'int16');
    fprintf(f,' A');
    fwrite(f,0,'int16');
    fprintf(f,' A');
    fwrite(f,0,'int16');
    %according to the documentation (http://planetside.co.uk/wiki/index.php?title=Terragen_.TER_Format)
    % The absolute altitude of a particular point equals BaseHeight Elevation * HeightScale / 65536
    fprintf(f,' AALTW');
    fwrite(f,(max(terrain(:))-min(terrain(:)))/5,'int16');%find what these mean, and get the values from terrain
    fwrite(f,min(terrain(:))/10,'int16');
    terrain=round((terrain-min(terrain(:)))./(max(terrain(:))-min(terrain(:)))*32768);
    for x=1:size(terrain,1)
        for y=1:size(terrain,2)
            fwrite(f,terrain(x,y),'uint16');
        end
    end
    fprintf(f,'EOF ');
    fclose(f);
end
