function closeAllBut(fig)
    fh = findall(0, 'type', 'figure');
    for id = 1:length(fh),
        if fh(id) ~= fig,
            close(fh(id));
        end
    end
end