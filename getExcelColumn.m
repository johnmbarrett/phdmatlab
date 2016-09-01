function column = getExcelColumn(columnIndex)
    if columnIndex < 26
        column = '';
    else
        column = char('A'+floor(columnIndex/26)-1);
    end

    column(end+1) = char('A'+mod(columnIndex,26));
end   