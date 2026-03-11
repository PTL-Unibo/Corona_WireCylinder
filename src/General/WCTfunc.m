function wct = WCTfunc(wct_seconds)

    h = floor(wct_seconds / 3600);
    m = floor(mod(wct_seconds, 3600) / 60);
    s = mod(wct_seconds, 60);

    parts = [];

    if h > 0
        parts = [parts, sprintf('%d h ', round(h))];
    end
    if m > 0 || h > 0
        parts = [parts, sprintf('%d m ', round(m))];
    end
    if s > 0 || (m == 0 && h == 0)  
        parts = [parts, sprintf('%d s', round(s))];
    end

    % Rimuovi spazi finali
    wct = strtrim(parts);
end