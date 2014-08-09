function out = ThousandSep(in)
    %add commas to the array of in numbers
    
    if strcmp(class(in),'cell')
        in = cell2mat(in);
    end
    
    import java.text.*
    v = DecimalFormat;

    for i = 1:numel(in)
        out{i} = char(v.format(in(i)));
    end
end
